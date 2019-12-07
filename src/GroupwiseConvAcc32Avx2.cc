/*
 * Copyright (c) Facebook, Inc. and its affiliates.
 * All rights reserved.
 * This source code is licensed under the BSD-style license found in the
 * LICENSE file in the root directory of this source tree.
 */
#include <asmjit/asmjit.h>
#include <cpuinfo.h>
#include <immintrin.h>
#include <array>
#include <iostream>
#include <map>
#include <stdexcept>
#include <tuple>
#include <type_traits>
#include "./CodeGenHelpers.h"
#include "./GroupwiseConv.h"
#include "./RefImplementations.h"
#include "./TransposeUtils.h"
#include "fbgemm/Fbgemm.h"

namespace fbgemm {

using namespace std;

namespace x86 = asmjit::x86;

template <int SPATIAL_DIM>
void calculateRowOffsets(
    const conv_param_t<SPATIAL_DIM>& conv_param,
    const uint8_t* activations,
    int32_t* rowOffsetBuf,
    int32_t a_zero_point,
    int groupNum) {
  int OH = conv_param.OUT_DIM[0];
  int OW = conv_param.OUT_DIM[1];
  int IH = conv_param.IN_DIM[0];
  int IW = conv_param.IN_DIM[1];
  int G = conv_param.G;
  int C_per_G = conv_param.IC / conv_param.G;
  int H_PAD = conv_param.pad[0];
  int W_PAD = conv_param.pad[1];
  // calculate row offset
  for (int h = 0; h < OH; ++h) {
    for (int w = 0; w < OW; ++w) {
      int32_t sum = 0;
      for (int r = 0; r < conv_param.K[0]; ++r) {
        int h_in = -H_PAD + h * conv_param.stride[0] + r;
        for (int s = 0; s < conv_param.K[1]; ++s) {
          int w_in = -W_PAD + w * conv_param.stride[1] + s;
          for (int c = 0; c < C_per_G; ++c) {
            int a_val;
            if (h_in < 0 || h_in >= IH || w_in < 0 || w_in >= IW) {
              a_val = a_zero_point;
            } else {
              a_val = activations
                  [((h_in * IW + w_in) * G + groupNum) * C_per_G + c];
            }
            sum += a_val;
          }
        }
      }
      rowOffsetBuf[h * OW + w] = sum;
    }
  }
}

template <int SPATIAL_DIM = 2>
kernel_sig_t getKernelSig(
    const conv_param_t<SPATIAL_DIM>& conv_param,
    bool isAZeroPointZero,
    bool needRowOffset,
    bool isTopEdgeIncluded,
    bool isBottomEdgeIncluded,
    bool accum) {
  // kernel is specialized on number of input channels per group, number of
  // output channels per group, whether stride is 1 or stride is 2, whether or
  // not zero point for activations is 0 or not, whether or not row offset
  // calculations are needed, whether or not top edge is included and whether or
  // not bottom edge is included.
  // use_padding_: If false, the right padding on the width side and bottom
  // padding on height side are not used for the case of stride = 2
  // accum: accumulate results for output and rowoffset
  int C_per_G = conv_param.IC / conv_param.G;
  int K_per_G = conv_param.OC / conv_param.G;
  auto kernelSig = make_tuple(
      isAZeroPointZero,
      needRowOffset,
      isTopEdgeIncluded,
      isBottomEdgeIncluded,
      conv_param.IN_DIM[1] % 2 == 0 && conv_param.stride[0] == 2,
      accum,
      conv_param.G,
      conv_param.stride[0],
      C_per_G,
      K_per_G);
  return kernelSig;
}

template <int SPATIAL_DIM = 2>
jit_conv_kernel_fp getOrCreateConvKernel(
    const conv_param_t<SPATIAL_DIM>& conv_param,
    int a_zero_point,
    bool needRowOffset,
    bool isTopEdgeIncluded,
    bool isBottomEdgeIncluded,
    bool accum) {
  // Note: Wrong code is generated if it's not one of the supported convolution
  assert(fbgemmOptimizedGConv<SPATIAL_DIM>(conv_param));
  auto kernelSig = getKernelSig(
      conv_param,
      a_zero_point == 0,
      needRowOffset,
      isTopEdgeIncluded,
      isBottomEdgeIncluded,
      accum);

  if (cpuinfo_initialize()) {
    if (fbgemmHasAvx512Support()) {
      // TODO: Currently dispatching to avx2 code
      return GenConvKernel<SPATIAL_DIM, inst_set_t::avx2>::codeCache_
          .getOrCreate(kernelSig, [&]() {
            auto genObj = GenConvKernel<SPATIAL_DIM, inst_set_t::avx2>(
                conv_param,
                a_zero_point,
                needRowOffset,
                isTopEdgeIncluded,
                isBottomEdgeIncluded,
                accum);
            return genObj.getOrCreate();
          });
    } else if (fbgemmHasAvx2Support()) {
      return GenConvKernel<SPATIAL_DIM, inst_set_t::avx2>::codeCache_
          .getOrCreate(kernelSig, [&]() {
            auto genObj = GenConvKernel<SPATIAL_DIM, inst_set_t::avx2>(
                conv_param,
                a_zero_point,
                needRowOffset,
                isTopEdgeIncluded,
                isBottomEdgeIncluded,
                accum);
            return genObj.getOrCreate();
          });
    } else {
      // TODO: Have default slower path
      assert(0 && "unsupported architecture");
    }
  } else {
    throw runtime_error("Failed to initialize cpuinfo!");
  }
  return nullptr;
}

template <int SPATIAL_DIM>
jit_conv_kernel_fp GenConvKernel<SPATIAL_DIM, inst_set_t::avx2>::getOrCreate() {
  asmjit::CodeHolder code;
  code.init(this->runtime().codeInfo());
  x86::Assembler assembler(&code);
  x86::Emitter* a = assembler.as<x86::Emitter>();

#if defined(FBGEMM_LOG_CODE)
  auto kernelSig = make_tuple(
      this->isAZeroPointZero_,
      this->needRowOffset_,
      this->isTopEdgeIncluded_,
      this->isBottomEdgeIncluded_,
      this->use_padding_,
      this->accum_,
      this->G_,
      this->STRIDE_,
      this->C_per_G_,
      this->K_per_G_);
  // log code to a file
  FILE* codeLogfile = fopen(this->getCodeLoggingFile(kernelSig).c_str(), "w");
  asmjit::FileLogger* codeLogger = new asmjit::FileLogger(codeLogfile);
  if (codeLogger) {
    code.setLogger(codeLogger);
  }
#endif

  // arguments to the function created
#ifdef _MSC_VER
  in_acts_R_ = a->zcx();
  wghts_R_ = a->zdx();
  out_acts_R_ = a->gpz(8);
  a_zero_pt_R_ = a->gpz(9);
  H_start_R_ = a->zdi();
  H_end_R_ = a->zsi();
#else
  in_acts_R_ = a->zdi();
  wghts_R_ = a->zsi();
  out_acts_R_ = a->zdx();
  a_zero_pt_R_ = a->zcx();
  H_start_R_ = a->gpz(8);
  H_end_R_ = a->gpz(9);
#endif
  W_R_ = a->gpz(10);
  row_offset_R_ = a->gpz(11);

  // register for temporary use
  scratchReg1_ = a->gpz(12);
  scratchReg2_ = a->gpz(13);

  func_.init(asmjit::FuncSignatureT<
             void,
             uint8_t*,
             int8_t*,
             int32_t*,
             int32_t,
             int32_t,
             int32_t,
             int32_t,
             int32_t*>(asmjit::CallConv::kIdHost));

  frame_.init(func_);

  frame_.setDirtyRegs(
      x86::Reg::kGroupVec,
      asmjit::Support::bitMask(0, 1, 2, 3, 4, 5, 6, 7) |
          asmjit::Support::bitMask(8, 9, 10, 11, 12, 13, 14, 15));
  frame_.setDirtyRegs(
      x86::Reg::kGroupGp,
      asmjit::Support::bitMask(8, 9, 10, 11, 12, 13, 14, 15));

  asmjit::FuncArgsAssignment args(&func_);
  args.assignAll(
      in_acts_R_,
      wghts_R_,
      out_acts_R_,
      a_zero_pt_R_,
      H_start_R_,
      H_end_R_,
      W_R_,
      row_offset_R_);

  args.updateFuncFrame(frame_);
  frame_.finalize();

  a->emitProlog(frame_);
  a->emitArgsAssignment(frame_, args);

  // We have run out of register so can't keep
  // this in a register. It's generated again at
  // each use. Only used for the case of C_per_G == 2 or 4
  // gen8BitVectorOne(a, oneReg8Bit_V_);
  gen16BitVectorOne(a, oneReg16Bit_V_);

  loopR1_ = a->gpz(14);
  loopR2_ = a->gpz(15);

  if (!this->isAZeroPointZero_) {
    broadcast8Bit(a, a_zero_pt_R_, zeroPTReg_V_);
  }

  genConstForPermutations(a);

  genForLoadingWeights(a);

  // W_R_ is an input to the JIT'ed kernel and is output image width.
  // W_R_ is passed in using stack. We reload it inside kernel where we need it.
  // The following logic calculates the input image width in the same register.
  // Only works for stride == 2
  if (this->STRIDE_ > 1) {
    a->imul(W_R_, W_R_, static_cast<asmjit::Imm>(this->STRIDE_));
    if (this->isImageEven_) {
      a->inc(W_R_);
    }
    a->sub(W_R_, static_cast<asmjit::Imm>(this->STRIDE_ - 1));
  }

  if (this->isTopEdgeIncluded_) {
    genForTopOrBottomEdge(a, true /* isTopEdge */);
  }
  genCoreInsts(a);
  if (this->isBottomEdgeIncluded_) {
    genForTopOrBottomEdge(a, false /* isTopEdge */);
  }

  a->emitEpilog(frame_);

  jit_conv_kernel_fp fn;
  asmjit::Error err;
  {
    unique_lock<mutex> lock(this->rtMutex_);
    err = this->runtime().add(&fn, &code);
  }

  if (err) {
    cout << "Error: in fn add" << endl;
    return nullptr;
  }

#if defined(FBGEMM_LOG_CODE)
  fclose(codeLogfile);
  delete codeLogger;
#endif

  return fn;
}

template <int SPATIAL_DIM>
void GenConvKernel<SPATIAL_DIM, inst_set_t::avx2>::initResultRegs(
    x86::Emitter* a) {
  if (kLoopIters_ > 0) {
    for (int k = 0; k < kLoopIters_; ++k) {
      a->vxorps(x86::Ymm(9 - k), x86::Ymm(9 - k), x86::Ymm(9 - k));
    }
  } else {
    a->vxorps(x86::Ymm(9), x86::Ymm(9), x86::Ymm(9));
  }
}

template <int SPATIAL_DIM>
void GenConvKernel<SPATIAL_DIM, inst_set_t::avx2>::genConstForPermutations(
    x86::Emitter* a) {
  if (this->C_per_G_ == 4) {
    x86::Gp permute_const_reg = a->gpz(12);
    x86::Xmm const_reg_xmm = x86::xmm11;
    // We have 1st group in even lanes and 2nd group in odd lanes.
    // Permute to put 1st group to lower 128-bit and 2nd group in upper
    // 128-bit.
    // load 7, 5, 3, 1, 6, 4, 2, 0 in a 64-bit reg
    a->mov(permute_const_reg, static_cast<asmjit::Imm>(0x0705030106040200));
    a->movq(const_reg_xmm, permute_const_reg);
    // Zero extend 8 packed 8-bit integers in the low 8 bytes of const_reg_xmm
    // to 8 packed 32-bit integers in stPermReg_V_
    a->vpmovzxbd(stPermReg_V_, const_reg_xmm);
  } else {
    // this->C_per_G_ == 2
    x86::Gp permute_const_reg = a->gpz(12);
    x86::Xmm const_reg_xmm = x86::xmm11;
    // We have 1st group in position 0 and 4, 2nd group 1 and 5 and so on.
    // Permute to put 1st group to lower 64-bit and 2nd group to next
    // 64-bit and so on.
    // load 7, 3, 6, 2, 5, 1, 4, 0 in a 64-bit reg
    a->mov(permute_const_reg, static_cast<asmjit::Imm>(0x0703060205010400));
    a->movq(const_reg_xmm, permute_const_reg);
    a->vpmovzxbd(stPermReg_V_, const_reg_xmm);
  }
}

template <int SPATIAL_DIM>
void GenConvKernel<SPATIAL_DIM, inst_set_t::avx2>::genForLoadingWeights(
    x86::Emitter* a) {
  using WRegs = x86::Ymm;
  int paddedICPerG = (this->C_per_G_ + 3) / 4 * 4;
  // load weights
  for (int r = 0; r < this->R_; ++r) {
    for (int s = 0; s < this->S_; ++s) {
      // For other cases, weights are too big to be kept in registers
      // and are loaded as they are used.
      if (this->C_per_G_ == 4 || this->C_per_G_ == 2) {
        a->vmovaps(
            WRegs(r * this->S_ + s),
            x86::dword_ptr(
                wghts_R_,
                (r * this->S_ + s) * this->K_per_G_ * GTogether_ *
                    paddedICPerG * sizeof(int8_t)));
      }
    }
  }
}

template <int SPATIAL_DIM>
void GenConvKernel<SPATIAL_DIM, inst_set_t::avx2>::storeResult(
    x86::Emitter* a) {
  if (GTogether_ > 1) {
    // store with permutation
    a->vpermd(x86::Ymm(9), stPermReg_V_, x86::Ymm(9));
    if (this->accum_) {
      a->vpaddd(x86::Ymm(9), x86::Ymm(9), x86::dword_ptr(out_acts_R_));
    }
    a->vmovups(x86::dword_ptr(out_acts_R_), x86::Ymm(9));
  } else {
    // horizontal add and store
    if (this->C_per_G_ == 8) {
      a->vphaddd(x86::Ymm(9), x86::Ymm(9), x86::Ymm(8));
      a->vpermq(x86::Ymm(9), x86::Ymm(9), static_cast<asmjit::Imm>(0xd8));
      if (this->accum_) {
        a->vpaddd(x86::Ymm(9), x86::Ymm(9), x86::dword_ptr(out_acts_R_));
      }
      a->vmovups(x86::dword_ptr(out_acts_R_), x86::Ymm(9));
    } else if (this->K_per_G_ == 16) {
      a->vphaddd(x86::Ymm(9), x86::Ymm(9), x86::Ymm(8));
      a->vpermq(x86::Ymm(9), x86::Ymm(9), static_cast<asmjit::Imm>(0xd8));

      a->vphaddd(x86::Ymm(7), x86::Ymm(7), x86::Ymm(6));
      a->vpermq(x86::Ymm(7), x86::Ymm(7), static_cast<asmjit::Imm>(0xd8));

      a->vphaddd(x86::Ymm(5), x86::Ymm(5), x86::Ymm(4));
      a->vpermq(x86::Ymm(5), x86::Ymm(5), static_cast<asmjit::Imm>(0xd8));

      a->vphaddd(x86::Ymm(3), x86::Ymm(3), x86::Ymm(2));
      a->vpermq(x86::Ymm(3), x86::Ymm(3), static_cast<asmjit::Imm>(0xd8));

      a->vphaddd(x86::Ymm(9), x86::Ymm(9), x86::Ymm(7));
      a->vpermq(x86::Ymm(9), x86::Ymm(9), static_cast<asmjit::Imm>(0xd8));

      a->vphaddd(x86::Ymm(5), x86::Ymm(5), x86::Ymm(3));
      a->vpermq(x86::Ymm(5), x86::Ymm(5), static_cast<asmjit::Imm>(0xd8));

      if (this->accum_) {
        a->vpaddd(x86::Ymm(9), x86::Ymm(9), x86::dword_ptr(out_acts_R_));
        a->vpaddd(x86::Ymm(5), x86::Ymm(5), x86::dword_ptr(out_acts_R_, 32));
      }
      a->vmovups(x86::dword_ptr(out_acts_R_), x86::Ymm(9));
      a->vmovups(x86::dword_ptr(out_acts_R_, 32), x86::Ymm(5));
    }
  }
}

template <int SPATIAL_DIM>
void GenConvKernel<SPATIAL_DIM, inst_set_t::avx2>::storeOffset(
    x86::Emitter* a) {
  switch (this->C_per_G_) {
    case 2:
      // store 128-bits containing rowoffset for four groups
      if (this->accum_) {
        a->paddd(rowOffsetReg_V_.half(), x86::dword_ptr(row_offset_R_));
      }
      a->vmovdqu(x86::dword_ptr(row_offset_R_), rowOffsetReg_V_.half());
      break;
    case 4:
      // store 64-bits containing rowoffset for two groups
      if (this->accum_) {
        a->vmovq(tmpReg1_V_.half(), x86::dword_ptr(row_offset_R_));
        a->paddd(rowOffsetReg_V_.half(), tmpReg1_V_.half());
      }
      a->vmovq(x86::dword_ptr(row_offset_R_), rowOffsetReg_V_.half());
      break;
    case 8:
      if (this->accum_) {
        a->vmovd(tmpReg1_V_.half(), x86::dword_ptr(row_offset_R_));
        a->paddd(rowOffsetReg_V_.half(), tmpReg1_V_.half());
      }
      a->vmovd(x86::dword_ptr(row_offset_R_), rowOffsetReg_V_.half());
      break;
    case 16:
      // rowOffsetReg_V_[0:63] has sum for first 8 and
      // rowOffsetReg_V_[64:127] has sum for second 8
      // execute vphaddd twice to sum the two
      a->vphaddd(rowOffsetReg_V_, rowOffsetReg_V_, rowOffsetReg_V_);
      a->vphaddd(rowOffsetReg_V_, rowOffsetReg_V_, rowOffsetReg_V_);
      if (this->accum_) {
        a->vmovd(tmpReg1_V_.half(), x86::dword_ptr(row_offset_R_));
        a->paddd(rowOffsetReg_V_.half(), tmpReg1_V_.half());
      }
      a->vmovd(x86::dword_ptr(row_offset_R_), rowOffsetReg_V_.half());
      break;
    default:
      assert(0 && "not supported case");
  }
}

template <int SPATIAL_DIM>
void GenConvKernel<SPATIAL_DIM, inst_set_t::avx2>::genForTopOrBottomEdge(
    x86::Emitter* a,
    bool isTopEdge) {
  if (isTopEdge) {
    // top-left corner code
    genForSingleOutput(
        a,
        true, // isLeft
        false, // isRight
        true, // isTop
        false // isBotom
    );
  } else {
    // bottom-left corner code
    genForSingleOutput(
        a,
        true, // isLeft
        false, // isRight
        false, // isTop
        this->use_padding_ // isBottom
    );
  }

  // edge excluding corners
  asmjit::Label LoopEdge = a->newLabel();
  a->mov(loopR2_, static_cast<asmjit::Imm>(this->W_PAD_));
  a->bind(LoopEdge);

  if (isTopEdge) {
    genForSingleOutput(
        a,
        false, // isLeft
        false, // isRight
        true, // isTop
        false // isBottom
    );
  } else {
    genForSingleOutput(
        a,
        false, // isLeft
        false, // isRight
        false, // isTop
        this->use_padding_ // isBottom
    );
  }

  a->inc(loopR2_);
  // Output width was passed in as the 7th argument (i.e., using stack).
  // Reload it from the same location.
  a->movsxd(
      loopR1_,
      x86::dword_ptr(
          x86::rsp, frame_.saOffsetFromSP() + func_.arg(6).stackOffset()));
  a->sub(loopR1_, static_cast<asmjit::Imm>(this->W_PAD_));
  a->cmp(loopR2_, loopR1_);
  a->jl(LoopEdge);

  if (isTopEdge) {
    // top-right corner code
    genForSingleOutput(
        a,
        false, // isLeft
        this->use_padding_, // isRight
        true, // isTop
        false // isBottom
    );
  } else {
    // bottom-right corner code
    genForSingleOutput(
        a,
        false, // isLeft
        this->use_padding_, // isRight
        false, // isTop
        this->use_padding_ // isBottom
    );
  }

  if (this->STRIDE_ > 1) {
    // STRIDE_ == 2 and even widths,
    // We increase it by C_;
    // STRIDE_ == 2 and odd widths, nothing to do
    // input ptr is already at the right position
    if (this->isImageEven_) {
      a->add(in_acts_R_, static_cast<asmjit::Imm>(this->C_ * sizeof(uint8_t)));
    }

  } else {
    // reset input activation pointer
    a->mov(scratchReg2_, W_R_);
    a->imul(scratchReg2_, static_cast<asmjit::Imm>(this->C_ * sizeof(uint8_t)));
    a->sub(
        scratchReg2_,
        static_cast<asmjit::Imm>(this->W_PAD_ * this->C_ * sizeof(uint8_t)));
    a->sub(in_acts_R_, scratchReg2_);
  }
}

template <int SPATIAL_DIM>
void GenConvKernel<SPATIAL_DIM, inst_set_t::avx2>::genForSingleOutput(
    x86::Emitter* a,
    bool isLeft,
    bool isRight,
    bool isTop,
    bool isBottom) {
  // init result regs
  initResultRegs(a);

  // row offset
  if (this->needRowOffset_) {
    a->vxorps(rowOffsetReg_V_, rowOffsetReg_V_, rowOffsetReg_V_);
  }
<<<<<<< HEAD
#endif

  // arguments to the function created
#ifdef _MSC_VER
  in_acts_R_ = a->zcx();
  a_zero_pt_R_ = a->zdx();
  H_R_ = a->gpz(8);
  W_R_ = a->gpz(9);
  row_offset_R_ = a->zdi();
#else
  in_acts_R_ = a->zdi();
  a_zero_pt_R_ = a->zsi();
  H_R_ = a->zdx();
  W_R_ = a->zcx();
  row_offset_R_ = a->gpz(8);
#endif

  // register for temporary use
  scratchReg1_ = a->gpz(12);
  scratchReg2_ = a->gpz(13);
=======
>>>>>>> upstream

  bool isWidthMiddle = !isLeft && !isRight;
  bool isHeightMiddle = !isTop && !isBottom;
  for (int r = 0; r < this->R_; ++r) {
    int h_in = r;
    if (isTop) {
      h_in = -this->H_PAD_ + r;
    }
    bool in_image_H = (isTop && h_in >= 0) ||
        (isBottom && h_in < (this->R_ - this->H_PAD_)) || isHeightMiddle;
    for (int s = 0; s < this->S_; ++s) {
      int w_in = s;
      if (isLeft) {
        w_in = -this->W_PAD_ + s;
      }
      bool in_image_W = (isLeft && w_in >= 0) ||
          (isRight && w_in < (this->S_ - this->W_PAD_)) || isWidthMiddle;
      if (in_image_H && in_image_W) {
        genForSingleFilterPoint(a, r, s, w_in, false);
      } else {
        if (!this->isAZeroPointZero_) {
          genForSingleFilterPoint(a, r, s, w_in, true);
        }
      }
    }
    if (in_image_H) {
      // advance input pointer by one row
      a->imul(
          scratchReg2_,
          W_R_,
          static_cast<asmjit::Imm>(this->C_ * sizeof(uint8_t)));
      a->add(in_acts_R_, scratchReg2_);
    }
  }

  storeResult(a);

  // row offset
  if (this->needRowOffset_) {
    storeOffset(a);
    a->add(
        row_offset_R_, static_cast<asmjit::Imm>(GTogether_ * sizeof(int32_t)));
  }

  // rewind input ptr
  if (isHeightMiddle) {
    a->imul(
        scratchReg2_,
        W_R_,
        static_cast<asmjit::Imm>(this->R_ * this->C_ * sizeof(uint8_t)));
  } else {
    a->imul(
        scratchReg2_,
        W_R_,
        static_cast<asmjit::Imm>(
            (this->R_ - this->H_PAD_) * this->C_ * sizeof(uint8_t)));
  }
  a->sub(in_acts_R_, scratchReg2_);

  // advance output pointer
  a->add(out_acts_R_, static_cast<asmjit::Imm>(this->K_ * sizeof(int32_t)));

  // advance input ptr
  if (!isLeft) {
    a->add(
        in_acts_R_,
        static_cast<asmjit::Imm>(this->STRIDE_ * this->C_ * sizeof(uint8_t)));
  } else {
    a->add(
        in_acts_R_,
        static_cast<asmjit::Imm>(
            (this->STRIDE_ - this->W_PAD_) * this->C_ * sizeof(uint8_t)));
  }
}

template <int SPATIAL_DIM>
void GenConvKernel<SPATIAL_DIM, inst_set_t::avx2>::genForSingleFilterPoint(
    x86::Emitter* a,
    int r,
    int s,
    int act_s,
    bool use_zero_reg) {
  using WRegs = x86::Ymm;

  if (GTogether_ > 1) {
    if (this->C_per_G_ == 2) {
      if (use_zero_reg) {
        a->vmovapd(actReg_V_, zeroPTReg_V_);
      } else {
        a->vbroadcastsd(
            actReg_V_,
            x86::word_ptr(in_acts_R_, (act_s * this->C_) * sizeof(uint8_t)));
      }
      a->vpmovzxwd(actReg_V_, actReg_V_.half());
    } else if (this->C_per_G_ == 4) {
      if (use_zero_reg) {
        a->vmovapd(actReg_V_, zeroPTReg_V_);
      } else {
        a->vbroadcastsd(
            actReg_V_,
            x86::dword_ptr(in_acts_R_, act_s * this->C_ * sizeof(uint8_t)));
      }
    }
    // row offset
    if (this->needRowOffset_) {
      genU8Sum4(a, actReg_V_, rowOffsetReg_V_, oneReg16Bit_V_, tmpReg1_V_);
    }
    genU8I8S32FMA(
        a,
        actReg_V_,
        WRegs(r * this->S_ + s),
        x86::Ymm(9),
        oneReg16Bit_V_,
        tmpReg1_V_);
  } else {
    if (this->C_per_G_ == 8) {
      if (use_zero_reg) {
        a->vmovapd(actReg_V_, zeroPTReg_V_);
      } else {
        a->vbroadcastsd(
            actReg_V_,
            x86::qword_ptr(in_acts_R_, act_s * this->C_ * sizeof(uint8_t)));
      }
    } else {
      // this->C_per_G_ == 16
      if (use_zero_reg) {
        a->vmovapd(actReg_V_, zeroPTReg_V_);
      } else {
        a->vbroadcastf128(
            actReg_V_,
            x86::oword_ptr(in_acts_R_, act_s * this->C_ * sizeof(uint8_t)));
      }
    }
    // row offset
    if (this->needRowOffset_) {
      genU8Sum8(a, actReg_V_, rowOffsetReg_V_, tmpReg1_V_);
    }
    int kLoopMultiplier = 32 / this->C_per_G_;
    for (int k = 0; k < kLoopIters_; ++k) {
      a->vmovups(
          WRegs(0),
          x86::dword_ptr(
              wghts_R_,
              (((r * this->S_ + s) * this->K_per_G_) + k * kLoopMultiplier) *
                  this->C_per_G_ * sizeof(int8_t)));
      genU8I8S32FMA(
          a, actReg_V_, WRegs(0), x86::Ymm(9 - k), oneReg16Bit_V_, tmpReg1_V_);
    }
  }
}

template <int SPATIAL_DIM>
void GenConvKernel<SPATIAL_DIM, inst_set_t::avx2>::genCoreInsts(
    x86::Emitter* a) {
  // Top edge and bottom edge calculations are done separately
  // so start from next and leave out the last
  if (this->isTopEdgeIncluded_) {
    a->inc(H_start_R_);
  }
  if (this->isBottomEdgeIncluded_) {
    a->dec(H_end_R_);
  }
  // main compute
  asmjit::Label LoopHStart = a->newLabel();
  asmjit::Label LoopHEnd = a->newLabel();
  asmjit::Label LoopWLabel = a->newLabel();

  // H loop
  a->mov(loopR1_, H_start_R_);
  a->jmp(LoopHEnd);
  a->bind(LoopHStart);
  a->inc(loopR1_);

  genForSingleOutput(
      a,
      true, // isLeft,
      false, // isRight
      false, // isTop
      false // isBottom
  );

  a->movsxd(
      scratchReg1_,
      x86::dword_ptr(
          x86::rsp, frame_.saOffsetFromSP() + func_.arg(6).stackOffset()));
  a->sub(scratchReg1_, static_cast<asmjit::Imm>(this->W_PAD_));
  // W loop
  a->mov(loopR2_, static_cast<asmjit::Imm>(this->W_PAD_));
  a->bind(LoopWLabel);

  genForSingleOutput(
      a,
      false, // isLeft,
      false, // isRight
      false, // isTop
      false // isBottom
  );

  a->inc(loopR2_);
  a->cmp(loopR2_, scratchReg1_);
  a->jl(LoopWLabel);

  genForSingleOutput(
      a,
      false, // isLeft
      this->use_padding_, // isRight
      false, // isTop
      false // isBottom
  );

  if (this->STRIDE_ > 1) {
    // STRIDE_ == 2 and even widths,
    // We increase it by extra C_;
    // STRIDE_ == 2 and odd widths, no extra C_
    assert(this->STRIDE_ == 2 && "Not supported case");
    a->mov(scratchReg2_, W_R_);
    if (this->isImageEven_) {
      a->add(scratchReg2_, static_cast<asmjit::Imm>(1));
    }
    a->imul(scratchReg2_, static_cast<asmjit::Imm>(this->C_ * sizeof(uint8_t)));
    a->add(in_acts_R_, scratchReg2_);
  } else {
    a->add(in_acts_R_, static_cast<asmjit::Imm>(this->C_ * sizeof(uint8_t)));
  }

  a->bind(LoopHEnd);
  a->cmp(loopR1_, H_end_R_);
  a->jl(LoopHStart);
}

/*
namespace {

template <
    typename packed_W,
    typename outType,
    typename processOutputType,
    int SPATIAL_DIM>
void fbgemmGroupwiseConvBase_(
    const conv_param_t<SPATIAL_DIM>& conv_param,
    const uint8_t* activations,
    int32_t a_zero_point,
    int32_t* rowOffsetBuf,
    packed_W& packed_weights,
    outType* out,
    int32_t* outBuffer,
    const processOutputType& outProcess,
    int thread_id,
    int num_threads) {
  int MB = conv_param.MB;
  int H = conv_param.OUT_DIM[0];
  int W = conv_param.OUT_DIM[1];
  int G = conv_param.G;
  int K_per_G = conv_param.OC / G;
  int C_per_G = conv_param.IC / G;
  int oh_ow = conv_param.OUT_DIM[0] * conv_param.OUT_DIM[1];
  int ih_iw = conv_param.IN_DIM[0] * conv_param.IN_DIM[1];

  assert(SPATIAL_DIM == 2 && "3D groupwise conv not supported yet");

  int32_t* rowOffsetTrDest = rowOffsetBuf ? rowOffsetBuf + 8 * ih_iw : nullptr;
  if (fbgemmOptimizedGConv<SPATIAL_DIM>(conv_param)) {
    assert(G % 8 == 0);
    // generate convolution kernel
    jit_conv_kernel_fp fpConv =
        getOrCreateConvKernel<SPATIAL_DIM>(conv_param, a_zero_point);
    // generate row offset kernel
    jit_rowoffset_kernel_fp fpRowoffset =
        getOrCreateRowOffsetKernel(conv_param, a_zero_point);
    for (int i = 0; i < MB; ++i) {
      const uint8_t* actStartBatch = activations + i * ih_iw * conv_param.IC;
      for (int gOuter = 0; gOuter < G; gOuter += 8) {
        // row offset is calcualted for 8 groups at a time.
        // The result is row offsets in the format IH*IW x G
        if (rowOffsetBuf) {
          fpRowoffset(
              actStartBatch + gOuter * C_per_G,
              a_zero_point,
              H,
              W,
              rowOffsetBuf);
          // Transpose to get row offsets in the format G x IH*IW
          internal::transpose_avx2(
              ih_iw,
              8,
              reinterpret_cast<const float*>(rowOffsetBuf),
              8,
              reinterpret_cast<float*>(rowOffsetTrDest),
              ih_iw);
        }
        int gLimit = gOuter + 8;
        // Work on 8 output channels at a time (8 * sizeof(int32_t) == 32B VLEN
        // of AVX2), and we need multiple groups if a group has not enough
        // number of channels.
        int gDelta = std::max(8 / C_per_G, 1);
        for (int g = gOuter; g < gLimit; g += gDelta) {
          int32_t* currOutBuf =
              outBuffer + i * oh_ow * conv_param.OC + g * K_per_G;
          const uint8_t* actStartGroup = actStartBatch + g * C_per_G;
          for (int k = 0; k < K_per_G; k += 8) {
            // Don't be confused with k above which refers to output channels.
            // k0 and k1 are filter dimensions (commonly 3 and 3)
            int k0 = conv_param.K[0];
            int k1 = conv_param.K[1];
            fpConv(
                actStartGroup,
                // packed weight is in G (C/4) R S K 4 layout for IC_per_G >= 8
                // in (G/2) R S K (2C) for IC_per_G == 4
                packed_weights.getBuf() +
                    (g * (C_per_G / 4) * k0 * k1 * K_per_G + k) * 4,
                currOutBuf + k,
                a_zero_point,
                H,
                W);
          } // k loop

          // Output processing should be called for each group
          for (int j = 0; j < gDelta; ++j) {
            // calculateRowOffsets(
            // conv_param, actStartGroup, rowOffsetBuf, a_zero_point, j);
            int32_t* rowOffsetForCurG = rowOffsetTrDest
                ? rowOffsetTrDest + ((g - gOuter) + j) * ih_iw
                : nullptr;
            // compare_buffers(rowOffsetBuf, rowOffsetForCurG,
            // conv_param.IN_DIM[0]*conv_param.IN_DIM[1], 1, 1, 100);

            // outProcess expects rowOffsetBuf to contain row offsets for the
            // current group
            memcpy(rowOffsetBuf, rowOffsetForCurG, ih_iw * sizeof(int32_t));

            if (fbgemmHasAvx512Support()) {
              // Currently use avx2 code
              outProcess.template f<inst_set_t::avx2>(
                  out,
                  currOutBuf + j * K_per_G,
                  {i * oh_ow, oh_ow, (g + j) * K_per_G, K_per_G},
                  K_per_G * G,
                  K_per_G * G);
            } else if (fbgemmHasAvx2Support()) {
              outProcess.template f<inst_set_t::avx2>(
                  out,
                  currOutBuf + j * K_per_G,
                  {i * oh_ow, oh_ow, (g + j) * K_per_G, K_per_G},
                  K_per_G * G,
                  K_per_G * G);
            } else {
              // TODO: Have default slower path
              assert(0 && "unsupported architecure");
            }
          } // j loop
        } // g loop
      } // gOuter loop
    } // i loop
  } else {
    // for the not supported cases, just execute the naive C implementation
    conv_ref(
        conv_param,
        activations,
        a_zero_point,
        packed_weights.getBuf(),
        outBuffer);
    for (int i = 0; i < conv_param.MB; ++i) {
      for (int g = 0; g < conv_param.G; ++g) {
        if (rowOffsetBuf) {
          calculateRowOffsets(
              conv_param,
              activations +
                  i * conv_param.IN_DIM[0] * conv_param.IN_DIM[1] *
                      conv_param.IC,
              rowOffsetBuf,
              a_zero_point,
              g);
        }
        outProcess.template f<inst_set_t::anyarch>(
            out,
            outBuffer + i * oh_ow * conv_param.OC + g * K_per_G,
            {i * oh_ow, oh_ow, g * K_per_G, K_per_G},
            K_per_G * G,
            K_per_G * G);
      }
    }
  }
}

} // namespace

template <
    typename packed_W,
    typename outType,
    typename processOutputType,
    int SPATIAL_DIM>
void fbgemmGroupwiseConv(
    const conv_param_t<SPATIAL_DIM>& conv_param,
    const uint8_t* activations,
    int32_t a_zero_point,
    int32_t* rowOffsetBuf,
    packed_W& packed_weights,
    outType* out,
    int32_t* outBuffer,
    const processOutputType& outProcess,
    int thread_id,
    int num_threads) {
  // TODO: Remove this when threading is supported.
  if (thread_id > 0) {
    return;
  }

  return fbgemmGroupwiseConvBase_<
      packed_W,
      outType,
      processOutputType,
      SPATIAL_DIM>(
      conv_param,
      activations,
      a_zero_point,
      rowOffsetBuf,
      packed_weights,
      out,
      outBuffer,
      outProcess,
      thread_id,
      num_threads);
}

*/

/*
 *
 * This function does exactly the same compute as the JIT'ed kernel
 */
template <int SPATIAL_DIM>
void kernel_compute(
    const conv_param_t<SPATIAL_DIM>& conv_p,
    const uint8_t* in_acts,
    int8_t* wghts,
    int32_t* out_acts,
    int32_t a_zero_pt,
    int32_t h_start,
    int32_t h_end,
    int32_t width,
    int32_t* rowOffset,
    bool accum) {
  int IW = conv_p.IN_DIM[1];
  int IC = conv_p.IC;
  int OC = conv_p.OC;
  int G = conv_p.G;
  int R = conv_p.K[0];
  int S = conv_p.K[1];
  int IC_per_G = conv_p.IC / G;
  int OC_per_G = conv_p.OC / G;
  int G_together = PackWeightMatrixForGConv<int8_t, int32_t, SPATIAL_DIM>::
      numOfGroupsTogether(conv_p);
  int paddedICPerG = (IC_per_G + 3) / 4 * 4;
  for (int h = h_start; h < h_end; ++h) {
    for (int w = 0; w < width; ++w) {
      for (int g = 0; g < G_together; ++g) {
        for (int k = 0; k < OC_per_G; ++k) {
          int sum = 0;
          int rowSum = 0;
          for (int r = 0; r < R; ++r) {
            int h_in = -conv_p.pad[0] + h * conv_p.stride[0] + r;
            for (int s = 0; s < S; ++s) {
              int w_in = -conv_p.pad[1] + w * conv_p.stride[1] + s;
              for (int c = 0; c < IC_per_G; ++c) {
                bool out_of_image = h_in < 0 || h_in >= conv_p.IN_DIM[0] ||
                    w_in < 0 || w_in >= conv_p.IN_DIM[1];
                int h_index = h_in;
                if (h_start > 0) {
                  h_index = (h - h_start) * conv_p.stride[1] + r;
                }
                int a = out_of_image
                    ? a_zero_pt
                    : in_acts[(h_index * IW + w_in) * IC + g * IC_per_G + c];
                int idx = (((r * S + s) * OC_per_G + k) * G_together + g) *
                        paddedICPerG +
                    c;
                int b = wghts[idx];
                sum += a * b;
                rowSum += a;
              }
            }
          }
          if (accum) {
            out_acts[((h - h_start) * width + w) * OC + g * OC_per_G + k] +=
                sum;
            if (k == 0) {
              // only accumulate for k == 0
              rowOffset[((h - h_start) * width + w) * G_together + g] += rowSum;
            }
          } else {
            out_acts[((h - h_start) * width + w) * OC + g * OC_per_G + k] = sum;
            rowOffset[((h - h_start) * width + w) * G_together + g] = rowSum;
          }
        }
      }
    }
  }
}

template <typename processOutputType, typename outT, typename inT>
void dispatchOutputProcessing(
    const processOutputType& outProcess,
    int32_t* rowOffsetBuf,
    outT* out,
    const inT* inp,
    const block_type_t& block,
    int ld_out,
    int ld_in,
    int groups,
    int C_per_G,
    true_type) {
  constexpr QuantizationGranularity Q_GRAN = processOutputType::QGRANType;
  constexpr int FUSE_RELU = processOutputType::RELU_FUSED;
  bool b_symmetric = (Q_GRAN == QuantizationGranularity::TENSOR &&
                      outProcess.getBZeroPoint()[0] == 0) ||
      rowOffsetBuf == nullptr;
  int32_t a_zero_point = outProcess.getAZeroPoint();

  // Requantization
  requantizationParams_t<typename processOutputType::BIAS_T> r = {
      a_zero_point,
      outProcess.getBZeroPoint(),
      outProcess.getCZeroPoint(),
      outProcess.getCMultiplier(),
      rowOffsetBuf,
      outProcess.getColOffsets(),
      outProcess.getBias(),
      outProcess.getNCols(),
      groups,
      outProcess.getActWScale()};

  if (C_per_G == 2) {
    if (a_zero_point == 0) {
      if (b_symmetric) {
        if (outProcess.getBias() == nullptr) {
          requantizeOutputProcessingGConvAvx2<
              true,
              true,
              Q_GRAN,
              false,
              FUSE_RELU,
              2>(out, inp, block, ld_out, ld_in, r);
        } else {
          requantizeOutputProcessingGConvAvx2<
              true,
              true,
              Q_GRAN,
              true,
              FUSE_RELU,
              2>(out, inp, block, ld_out, ld_in, r);
        }
      } else {
        if (outProcess.getBias() == nullptr) {
          requantizeOutputProcessingGConvAvx2<
              true,
              false,
              Q_GRAN,
              false,
              FUSE_RELU,
              2>(out, inp, block, ld_out, ld_in, r);
        } else {
          requantizeOutputProcessingGConvAvx2<
              true,
              false,
              Q_GRAN,
              true,
              FUSE_RELU,
              2>(out, inp, block, ld_out, ld_in, r);
        }
      }
    } else {
      if (b_symmetric) {
        if (outProcess.getBias() == nullptr) {
          requantizeOutputProcessingGConvAvx2<
              false,
              true,
              Q_GRAN,
              false,
              FUSE_RELU,
              2>(out, inp, block, ld_out, ld_in, r);
        } else {
          requantizeOutputProcessingGConvAvx2<
              false,
              true,
              Q_GRAN,
              true,
              FUSE_RELU,
              2>(out, inp, block, ld_out, ld_in, r);
        }
      } else {
        if (outProcess.getBias() == nullptr) {
          requantizeOutputProcessingGConvAvx2<
              false,
              false,
              Q_GRAN,
              false,
              FUSE_RELU,
              2>(out, inp, block, ld_out, ld_in, r);
        } else {
          requantizeOutputProcessingGConvAvx2<
              false,
              false,
              Q_GRAN,
              true,
              FUSE_RELU,
              2>(out, inp, block, ld_out, ld_in, r);
        }
      }
    }
  } else if (C_per_G == 4) {
    if (a_zero_point == 0) {
      if (b_symmetric) {
        if (outProcess.getBias() == nullptr) {
          requantizeOutputProcessingGConvAvx2<
              true,
              true,
              Q_GRAN,
              false,
              FUSE_RELU,
              4>(out, inp, block, ld_out, ld_in, r);
        } else {
          requantizeOutputProcessingGConvAvx2<
              true,
              true,
              Q_GRAN,
              true,
              FUSE_RELU,
              4>(out, inp, block, ld_out, ld_in, r);
        }
      } else {
        if (outProcess.getBias() == nullptr) {
          requantizeOutputProcessingGConvAvx2<
              true,
              false,
              Q_GRAN,
              false,
              FUSE_RELU,
              4>(out, inp, block, ld_out, ld_in, r);
        } else {
          requantizeOutputProcessingGConvAvx2<
              true,
              false,
              Q_GRAN,
              true,
              FUSE_RELU,
              4>(out, inp, block, ld_out, ld_in, r);
        }
      }
    } else {
      if (b_symmetric) {
        if (outProcess.getBias() == nullptr) {
          requantizeOutputProcessingGConvAvx2<
              false,
              true,
              Q_GRAN,
              false,
              FUSE_RELU,
              4>(out, inp, block, ld_out, ld_in, r);
        } else {
          requantizeOutputProcessingGConvAvx2<
              false,
              true,
              Q_GRAN,
              true,
              FUSE_RELU,
              4>(out, inp, block, ld_out, ld_in, r);
        }
      } else {
        if (outProcess.getBias() == nullptr) {
          requantizeOutputProcessingGConvAvx2<
              false,
              false,
              Q_GRAN,
              false,
              FUSE_RELU,
              4>(out, inp, block, ld_out, ld_in, r);
        } else {
          requantizeOutputProcessingGConvAvx2<
              false,
              false,
              Q_GRAN,
              true,
              FUSE_RELU,
              4>(out, inp, block, ld_out, ld_in, r);
        }
      }
    }
  } else if (C_per_G == 8) {
    if (a_zero_point == 0) {
      if (b_symmetric) {
        if (outProcess.getBias() == nullptr) {
          requantizeOutputProcessingGConvAvx2<
              true,
              true,
              Q_GRAN,
              false,
              FUSE_RELU,
              8>(out, inp, block, ld_out, ld_in, r);
        } else {
          requantizeOutputProcessingGConvAvx2<
              true,
              true,
              Q_GRAN,
              true,
              FUSE_RELU,
              8>(out, inp, block, ld_out, ld_in, r);
        }
      } else {
        if (outProcess.getBias() == nullptr) {
          requantizeOutputProcessingGConvAvx2<
              true,
              false,
              Q_GRAN,
              false,
              FUSE_RELU,
              8>(out, inp, block, ld_out, ld_in, r);
        } else {
          requantizeOutputProcessingGConvAvx2<
              true,
              false,
              Q_GRAN,
              true,
              FUSE_RELU,
              8>(out, inp, block, ld_out, ld_in, r);
        }
      }
    } else {
      if (b_symmetric) {
        if (outProcess.getBias() == nullptr) {
          requantizeOutputProcessingGConvAvx2<
              false,
              true,
              Q_GRAN,
              false,
              FUSE_RELU,
              8>(out, inp, block, ld_out, ld_in, r);
        } else {
          requantizeOutputProcessingGConvAvx2<
              false,
              true,
              Q_GRAN,
              true,
              FUSE_RELU,
              8>(out, inp, block, ld_out, ld_in, r);
        }
      } else {
        if (outProcess.getBias() == nullptr) {
          requantizeOutputProcessingGConvAvx2<
              false,
              false,
              Q_GRAN,
              false,
              FUSE_RELU,
              8>(out, inp, block, ld_out, ld_in, r);
        } else {
          requantizeOutputProcessingGConvAvx2<
              false,
              false,
              Q_GRAN,
              true,
              FUSE_RELU,
              8>(out, inp, block, ld_out, ld_in, r);
        }
      }
    }
  } else {
    if (a_zero_point == 0) {
      if (b_symmetric) {
        if (outProcess.getBias() == nullptr) {
          requantizeOutputProcessingGConvAvx2<
              true,
              true,
              Q_GRAN,
              false,
              FUSE_RELU,
              16>(out, inp, block, ld_out, ld_in, r);
        } else {
          requantizeOutputProcessingGConvAvx2<
              true,
              true,
              Q_GRAN,
              true,
              FUSE_RELU,
              16>(out, inp, block, ld_out, ld_in, r);
        }
      } else {
        if (outProcess.getBias() == nullptr) {
          requantizeOutputProcessingGConvAvx2<
              true,
              false,
              Q_GRAN,
              false,
              FUSE_RELU,
              16>(out, inp, block, ld_out, ld_in, r);
        } else {
          requantizeOutputProcessingGConvAvx2<
              true,
              false,
              Q_GRAN,
              true,
              FUSE_RELU,
              16>(out, inp, block, ld_out, ld_in, r);
        }
      }
    } else {
      if (b_symmetric) {
        if (outProcess.getBias() == nullptr) {
          requantizeOutputProcessingGConvAvx2<
              false,
              true,
              Q_GRAN,
              false,
              FUSE_RELU,
              16>(out, inp, block, ld_out, ld_in, r);
        } else {
          requantizeOutputProcessingGConvAvx2<
              false,
              true,
              Q_GRAN,
              true,
              FUSE_RELU,
              16>(out, inp, block, ld_out, ld_in, r);
        }
      } else {
        if (outProcess.getBias() == nullptr) {
          requantizeOutputProcessingGConvAvx2<
              false,
              false,
              Q_GRAN,
              false,
              FUSE_RELU,
              16>(out, inp, block, ld_out, ld_in, r);
        } else {
          requantizeOutputProcessingGConvAvx2<
              false,
              false,
              Q_GRAN,
              true,
              FUSE_RELU,
              16>(out, inp, block, ld_out, ld_in, r);
        }
      }
    }
  }
}

template <
    typename packed_W,
    typename outType,
    bool FUSE_RELU,
    QuantizationGranularity Q_GRAN,
    int SPATIAL_DIM,
    typename BIAS_TYPE>
void fbgemmGroupwiseConv(
    const conv_param_t<SPATIAL_DIM>& conv_param,
    const uint8_t* activations,
    int32_t a_zero_point,
    int32_t* rowOffsetBuf,
    packed_W& packed_weights,
    outType* out,
    int32_t* outBuffer,
    const ReQuantizeOutput<FUSE_RELU, Q_GRAN, BIAS_TYPE>& outProcess,
    int thread_id,
    int num_threads) {
  using processOutputType = ReQuantizeOutput<FUSE_RELU, Q_GRAN, BIAS_TYPE>;

  if (!cpuinfo_initialize()) {
    throw runtime_error("Failed to initialize cpuinfo!");
  }

  int MB = conv_param.MB;
  int OT = SPATIAL_DIM == 2 ? 1 : conv_param.OUT_DIM[SPATIAL_DIM - 3];
  int OH = conv_param.OUT_DIM[SPATIAL_DIM - 2];
  int OW = conv_param.OUT_DIM[SPATIAL_DIM - 1];
  int T = SPATIAL_DIM == 2 ? 1 : conv_param.K[SPATIAL_DIM - 3];
  int R = conv_param.K[SPATIAL_DIM - 2];
  int S = conv_param.K[SPATIAL_DIM - 1];
  int G = conv_param.G;
  int OC = conv_param.OC;
  int IC = conv_param.IC;
  int K_per_G = conv_param.OC / G;
  int C_per_G = conv_param.IC / G;
  int OH_OW = OH * OW;
  int OT_OH_OW = OT * OH * OW;
  int IT = SPATIAL_DIM == 2 ? 1 : conv_param.IN_DIM[SPATIAL_DIM - 3];
  int IH = conv_param.IN_DIM[SPATIAL_DIM - 2];
  int IW = conv_param.IN_DIM[SPATIAL_DIM - 1];
  int IH_IW = IH * IW;
  int IT_IH_IW = IT * IH * IW;
  int paddedCPerG = (C_per_G + 3) / 4 * 4;

  bool b_symmetric = (Q_GRAN == QuantizationGranularity::TENSOR &&
                      outProcess.getBZeroPoint()[0] == 0) ||
      rowOffsetBuf == nullptr;
  int G_together = PackWeightMatrixForGConv<int8_t, int32_t, SPATIAL_DIM>::
      numOfGroupsTogether(conv_param);

  if (SPATIAL_DIM == 2) {
    // Parallelization:
    int batch_start = 0;
    int batch_end = MB;
    int oh_start = 0;
    int oh_end = OH;
    if (MB >= num_threads) {
      fbgemmPartition1D(thread_id, num_threads, MB, batch_start, batch_end);
    } else {
      fbgemmPartition1D(thread_id, num_threads, OH, oh_start, oh_end);
    }

    if (batch_start >= batch_end || oh_start >= oh_end) {
      // There is no work for this thread
      return;
    }

    // generate convolution  + rowOffset kernel
    bool calculateRowOffset = !b_symmetric;
    bool isTopEdgeIncluded = oh_start == 0;
    bool isBottomEdgeIncluded = oh_end == OH;
    jit_conv_kernel_fp fpConv = getOrCreateConvKernel<SPATIAL_DIM>(
        conv_param,
        a_zero_point,
        calculateRowOffset,
        isTopEdgeIncluded,
        isBottomEdgeIncluded,
        false);

    int ih_start = 0;
    if (oh_start > 0) {
      ih_start = -conv_param.pad[SPATIAL_DIM - 2] +
          oh_start * conv_param.stride[SPATIAL_DIM - 2];
    }
    int32_t* out_start = outBuffer + oh_start * OW * OC;
    const uint8_t* in_start = activations + ih_start * IW * IC;
    int32_t* rowOffsetBuf_start = rowOffsetBuf + oh_start * OW * G_together;
    for (int i = batch_start; i < batch_end; ++i) {
      const uint8_t* in_start_batch = in_start + i * IH_IW * conv_param.IC;
      int32_t* out_start_batch = out_start + i * OH_OW * OC;
      int32_t* rowOffsetBuf_start_batch =
          rowOffsetBuf_start + i * OH_OW * G_together;
      for (int g = 0; g < G; g += G_together) {
        const uint8_t* in_start_group = in_start_batch + g * C_per_G;
        int8_t* weight_start =
            packed_weights.getBuf() + g * R * S * K_per_G * paddedCPerG;
        int32_t* out_start_group = out_start_batch;
        int32_t* rowOffsetBuf_start_group = rowOffsetBuf_start_batch;
        // Uncomment the following two lines to stop
        // reuse of output and rowoffset buffer
        // out_start_group = out_start_batch + g * K_per_G;
        // rowOffsetBuf_start_group = rowOffsetBuf_start_batch + g * MB * OH_OW;

        // exactly the same compute as the JIT'ed below
        // kernel_compute(
        //    conv_param,
        //    in_start_group,
        //    weight_start,
        //    out_start_group,
        //    a_zero_point,
        //    oh_start,
        //    oh_end,
        //    OW,
        //    rowOffsetBuf_start_group);

        fpConv(
            in_start_group,
            weight_start,
            out_start_group,
            a_zero_point,
            oh_start,
            oh_end,
            OW,
            rowOffsetBuf_start_group);

        const int32_t* inp = out_start_group;
        block_type_t block{i * OT_OH_OW + oh_start * OW,
                           (oh_end - oh_start) * OW,
                           g * K_per_G,
                           G_together * K_per_G};
        int ld_out = G * K_per_G;
        int ld_in = G * K_per_G;

        dispatchOutputProcessing(
            outProcess,
            rowOffsetBuf_start_group,
            out,
            inp,
            block,
            ld_out,
            ld_in,
            G,
            C_per_G,
            is_requantization<processOutputType>());
      }
    }
  } else {
    assert(SPATIAL_DIM == 3 && "Unsupported SPATIAL_DIM");

    conv_param_t<> conv_p_2d(
        conv_param.MB,
        conv_param.IC,
        conv_param.OC,
        {conv_param.IN_DIM[SPATIAL_DIM - 2],
         conv_param.IN_DIM[SPATIAL_DIM - 1]},
        conv_param.G,
        {conv_param.K[SPATIAL_DIM - 2], conv_param.K[SPATIAL_DIM - 1]},
        {conv_param.stride[SPATIAL_DIM - 2],
         conv_param.stride[SPATIAL_DIM - 1]},
        {conv_param.pad[1],
         conv_param.pad[2],
         conv_param.pad[4],
         conv_param.pad[5]});

    // Parallelization:
    int batch_start = 0;
    int batch_end = MB;
    int oh_start = 0;
    int oh_end = OH;
    if (MB >= num_threads) {
      fbgemmPartition1D(thread_id, num_threads, MB, batch_start, batch_end);
    } else {
      fbgemmPartition1D(thread_id, num_threads, OH, oh_start, oh_end);
    }

    if (batch_start >= batch_end || oh_start >= oh_end) {
      // There is no work for this thread
      return;
    }

    // generate convolution  + rowOffset kernel
    bool calculateRowOffset = !b_symmetric;
    bool isTopEdgeIncluded = oh_start == 0;
    bool isBottomEdgeIncluded = oh_end == OH;
    jit_conv_kernel_fp fpConvNoAccum = getOrCreateConvKernel<2>(
        conv_p_2d,
        a_zero_point,
        calculateRowOffset,
        isTopEdgeIncluded,
        isBottomEdgeIncluded,
        false);
    jit_conv_kernel_fp fpConvAccum = getOrCreateConvKernel<2>(
        conv_p_2d,
        a_zero_point,
        calculateRowOffset,
        isTopEdgeIncluded,
        isBottomEdgeIncluded,
        true);
    jit_conv_kernel_fp fpConv;

    int ih_start = 0;
    if (oh_start > 0) {
      ih_start = -conv_p_2d.pad[0] + oh_start * conv_p_2d.stride[0];
    }

    vector<uint8_t> zero_points(IH * IW * IC, a_zero_point);
    int32_t* out_start = outBuffer + oh_start * OW * OC;
    const uint8_t* in_start = activations + ih_start * IW * IC;
    int32_t* rowOffsetBuf_start = rowOffsetBuf + oh_start * OW * G_together;
    for (int i = batch_start; i < batch_end; ++i) {
      const uint8_t* in_start_batch = in_start + i * IT_IH_IW * IC;
      int32_t* out_start_batch = out_start + i * OT_OH_OW * OC;
      int32_t* rowOffsetBuf_start_batch =
          rowOffsetBuf_start + i * OT_OH_OW * G_together;
      for (int g = 0; g < G; g += G_together) {
        const uint8_t* in_start_group = in_start_batch + g * C_per_G;
        int8_t* weight_start =
            packed_weights.getBuf() + g * T * R * S * K_per_G * paddedCPerG;
        int32_t* out_start_group = out_start_batch;
        int32_t* rowOffsetBuf_start_group = rowOffsetBuf_start_batch;
        // Uncomment the following two lines to stop
        // reuse of output and rowoffset buffer
        // out_start_group = out_start_batch + g * K_per_G;
        // rowOffsetBuf_start_group = rowOffsetBuf_start_batch + g * MB *
        // OT_OH_OW;

        for (int ot = 0; ot < OT; ++ot) {
          int32_t* out_start_t = out_start_group + ot * OH_OW * OC;
          int32_t* rowOffsetBuf_start_t =
              rowOffsetBuf_start_group + ot * OH_OW * G_together;
          for (int t = 0; t < T; ++t) {
            int t_in = -conv_param.pad[0] + ot * conv_param.stride[0] + t;
            const uint8_t* in_start_t = in_start_group + t_in * IH_IW * IC;
            int8_t* weight_start_t =
                weight_start + t * R * S * K_per_G * G_together * paddedCPerG;
            if (t_in < 0 || t_in >= IT) {
              in_start_t = zero_points.data();
            }
            // exactly the same compute as the JIT'ed below
            // kernel_compute(
            // conv_p_2d,
            // in_start_t,
            // weight_start_t,
            // out_start_t,
            // a_zero_point,
            // oh_start,
            // oh_end,
            // OW,
            // rowOffsetBuf_start_t,
            // t > 0);

            fpConv = t > 0 ? fpConvAccum : fpConvNoAccum;
            fpConv(
                in_start_t,
                weight_start_t,
                out_start_t,
                a_zero_point,
                oh_start,
                oh_end,
                OW,
                rowOffsetBuf_start_t);
          }

          const int32_t* inp = out_start_t;
          block_type_t block{i * OT_OH_OW + oh_start * OW,
                             (oh_end - oh_start) * OW,
                             g * K_per_G,
                             G_together * K_per_G};
          int ld_out = G * K_per_G;
          int ld_in = G * K_per_G;

          dispatchOutputProcessing(
              outProcess,
              rowOffsetBuf_start_t,
              out + ot * OH_OW * OC,
              inp,
              block,
              ld_out,
              ld_in,
              G,
              C_per_G,
              is_requantization<processOutputType>());
        }
      }
    }
  }
}

template <int SPATIAL_DIM>
int rowOffsetBufferSizeGConv(const conv_param_t<SPATIAL_DIM>& conv_param) {
  // row offset buffer should be a able to hold row offsets for however
  // number of groups we process at a time.
  if (cpuinfo_initialize()) {
    int OT = SPATIAL_DIM == 2 ? 1 : conv_param.OUT_DIM[SPATIAL_DIM - 3];
    int bufferSize = OT * conv_param.OUT_DIM[SPATIAL_DIM - 2] *
        conv_param.OUT_DIM[SPATIAL_DIM - 1];
    if (fbgemmHasAvx512Support()) {
      return conv_param.MB * bufferSize * conv_param.G;
    } else if (fbgemmHasAvx2Support()) {
      return conv_param.MB * bufferSize * conv_param.G;
    } else {
      // TODO: Have default slower path
      assert(0 && "unsupported architecture");
      return -1;
    }
  } else {
    throw runtime_error("Failed to initialize cpuinfo!");
  }
}

template int rowOffsetBufferSizeGConv<2>(const conv_param_t<2>& conv_param);
template int rowOffsetBufferSizeGConv<3>(const conv_param_t<3>& conv_param);

#define INSTANTIATE_BASE(RELU, Q_GRAN, SPATIAL_DIM, BIAS_TYPE)                \
  template void fbgemmGroupwiseConv(                                          \
      const conv_param_t<SPATIAL_DIM>& conv_param,                            \
      const uint8_t* activations,                                             \
      int32_t a_zero_point,                                                   \
      int32_t* rowOffsetBuf,                                                  \
      PackWeightMatrixForGConv<int8_t, int32_t, SPATIAL_DIM>& packed_weights, \
      uint8_t* out,                                                           \
      int32_t* outBuffer,                                                     \
      const ReQuantizeOutput<RELU, Q_GRAN, BIAS_TYPE>& outProcess,            \
      int thread_id,                                                          \
      int num_threads);

#define INSTANTIATE_BIAS_T(RELU, Q_GRAN, SPATIAL_DIM) \
  INSTANTIATE_BASE(RELU, Q_GRAN, SPATIAL_DIM, float); \
  INSTANTIATE_BASE(RELU, Q_GRAN, SPATIAL_DIM, int32_t);

#define INSTANTIATE_SPATIAL_DIM(RELU, Q_GRAN) \
  INSTANTIATE_BIAS_T(RELU, Q_GRAN, 2);        \
  INSTANTIATE_BIAS_T(RELU, Q_GRAN, 3);

#define INSTANTIATE_Q_GRANS(RELU)                                 \
  INSTANTIATE_SPATIAL_DIM(RELU, QuantizationGranularity::TENSOR); \
  INSTANTIATE_SPATIAL_DIM(RELU, QuantizationGranularity::GROUP);  \
  INSTANTIATE_SPATIAL_DIM(RELU, QuantizationGranularity::OUT_CHANNEL);

INSTANTIATE_Q_GRANS(false);
INSTANTIATE_Q_GRANS(true);

#undef INSTANTIATE_Q_GRANS
#undef INSTANTIATE_SPATIAL_DIM
#undef INSTANTIATE_BIAS_T
#undef INSTANTIATE_BASE

/*
template void fbgemmGroupwiseConv(
    const conv_param_t<2>& conv_param,
    const uint8_t* activations,
    int32_t a_zero_point,
    int32_t* rowOffsetBuf,
    PackWeightMatrixForGConv<int8_t, int32_t, 2>& packed_weights,
    int32_t* out,
    int32_t* outBuffer,
    const DoNothing<int32_t, int32_t>& outProcess,
    int thread_id,
    int num_threads);
*/

} // namespace fbgemm
