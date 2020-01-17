/*
 * Copyright (c) Facebook, Inc. and its affiliates.
 * All rights reserved.
 * This source code is licensed under the BSD-style license found in the
 * LICENSE file in the root directory of this source tree.
 */
#include "./FbgemmFP16UKernelsAvx512.h"
#include <immintrin.h>
#include <iostream>

namespace fbgemm {

void gemmkernel_Avx512_fp16_fA0fB0fC0(GemmParamsFP16* gp, size_t kernel_nrow) {
  // register buffer
  __m512 zmmSum[28];
  size_t idxA = 0, idxB = 0, idxC = 0;
  // load beta
  __m128 xmmBeta = _mm_broadcast_ss(&gp->beta);
  __m512 zmmBeta = _mm512_broadcastss_ps(xmmBeta);

  // ldc in float size
  size_t ldc_floatsize = gp->ldc / sizeof(float);

  // Initialize register
  for (size_t ii = 0; ii < kernel_nrow; ii++) {
    zmmSum[2 * ii] = _mm512_setzero_ps();
    zmmSum[2 * ii + 1] = _mm512_setzero_ps();
  }

  for(uint64_t ii = 0; ii < gp->b_block_cols; ii++) {
    // reset index
    idxA = 0;

    // one inner loop
    uint64_t kk = 1;
    // load B
    __m256i b_load0 = _mm256_load_si256((__m256i*)(gp->B + idxB));
    __m256i b_load1 = _mm256_load_si256((__m256i*)(gp->B + idxB + 16));
    idxB += 32;
    __m512 zmmB0 = _mm512_cvtph_ps(b_load0);
    __m512 zmmB1 = _mm512_cvtph_ps(b_load1);

    if(gp->beta != 0) {  // accumulate
      for(size_t jj = 0; jj < kernel_nrow; jj++) {
        // load C
        __m512 zmmC0 = _mm512_load_ps(gp->C + idxC + jj * ldc_floatsize);
        __m512 zmmC1 = _mm512_load_ps(gp->C + idxC + 16 + jj * ldc_floatsize);
        zmmSum[2 * jj] = _mm512_mul_ps(zmmBeta, zmmC0);
        zmmSum[2 * jj + 1] = _mm512_mul_ps(zmmBeta, zmmC1);

        // load A
        __m128 a_load = _mm_broadcast_ss(gp->A + idxA + jj);
        __m512 zmmA = _mm512_broadcastss_ps(a_load);
        idxA += 1;
        zmmSum[2 * jj] = _mm512_fmadd_ps(zmmA, zmmB0, zmmSum[2 * jj]);
        zmmSum[2 * jj + 1] = _mm512_fmadd_ps(zmmA, zmmB1, zmmSum[2 * jj + 1]);
      }
    } else {  // set zero
      for(size_t jj = 0; jj < kernel_nrow; jj++) {
        // load A
        __m128 a_load = _mm_broadcast_ss(gp->A + idxA + jj);
        __m512 zmmA = _mm512_broadcastss_ps(a_load);
        idxA += 1;
        zmmSum[2 * jj] = _mm512_mul_ps(zmmA, zmmB0);
        zmmSum[2 * jj + 1] = _mm512_mul_ps(zmmA, zmmB1);
      }
    }

    // inner loop
    for(; kk < gp->k; kk++) {
      // load B
      b_load0 = _mm256_load_si256((__m256i*)(gp->B + idxB));
      b_load1 = _mm256_load_si256((__m256i*)(gp->B + idxB + 16));
      idxB += 32;
      zmmB0 = _mm512_cvtph_ps(b_load0);
      zmmB1 = _mm512_cvtph_ps(b_load1);
      for(size_t jj = 0; jj < kernel_nrow; jj++) {
        // load A
        __m128 a_load = _mm_broadcast_ss(gp->A + idxA + jj);
        __m512 zmmA = _mm512_broadcastss_ps(a_load);
        idxA += 1;
        zmmSum[2 * jj] = _mm512_fmadd_ps(zmmA, zmmB0, zmmSum[2 * jj]);
        zmmSum[2 * jj + 1] = _mm512_fmadd_ps(zmmA, zmmB1, zmmSum[2 * jj + 1]);
      }
    }

    // store C
    for(size_t jj = 0; jj < kernel_nrow; jj++) {
      _mm512_store_ps(gp->C + idxC + jj * ldc_floatsize, zmmSum[2 * jj]);
      _mm512_store_ps(gp->C + idxC + 16 + jj * ldc_floatsize, zmmSum[2 * jj + 1]);
    }

    idxC += 32;
  }
}

void NOINLINE
gemmkernel_1x2_Avx512_fp16_fA0fB0fC0(GemmParamsFP16* gp) {
//   asm volatile(
// #if !defined(__clang__)
//       "mov r14, %[gp]\t\n"
// #else
//       "mov %[gp], %%r14\t\n"
//       ".intel_syntax noprefix\t\n"
// #endif

//       // Copy parameters
//       // k
//       "mov r8, [r14 + 0]\t\n"
//       "dec r8\t\n"
//       // A
//       "mov r9, [r14 + 8]\t\n"
//       // B
//       "mov r10, [r14 + 16]\t\n"
//       // beta
//       "lea r15, [r14 + 24]\t\n"
//       // C
//       "mov r12, [r14 + 32]\t\n"
//       // ldc
//       "mov r13, [r14 + 40]\t\n"
//       // b_block_cols
//       "mov rdi, [r14 + 48]\t\n"
//       // b_block_size
//       "mov rsi, [r14 + 56]\t\n"

//       // Make copies of A and C
//       "mov rax, r9\t\n"
//       "mov rcx, r12\t\n"

//       "mov rbx, 0\t\n"
//       "loop_outter%=:\t\n"
//       "mov r14, r8\t\n"
//       "vbroadcastss zmm31,DWORD PTR [r15]\t\n"
//       "vcvtph2ps zmm3,YMMWORD PTR [r10 + 0]\t\n"
//       "vcvtph2ps zmm4,YMMWORD PTR [r10 + 32]\t\n"
//       "vxorps xmm0, xmm0, xmm0\t\n"
//       "vcomiss xmm31, xmm0\t\n"
//       "jz zero_regs%=\t\n"

//       // Setup values with beta multiplication
//       "vmulps zmm0, zmm31, [r12 + 0]\t\n"
//       "vmulps zmm1, zmm31, [r12 + 64]\t\n"
//       "test r14,r14\t\n"
//       "jz skip_preload%=\t\n"
//       "vcvtph2ps zmm31,YMMWORD PTR [r10 + 64]\t\n"
//       "skip_preload%=:\t\n"
//       "vbroadcastss zmm2,DWORD PTR [r9+0]\t\n"
//       "vfmadd231ps zmm0,zmm3,zmm2\t\n"
//       "vfmadd231ps zmm1,zmm4,zmm2\t\n"
//       "test r14,r14\t\n"
//       "jnz next_inner%=\t\n"
//       "add r10,64\t\n"
//       "jmp dump_C%=\t\n"

//       "zero_regs%=:\t\n"

//       "test r14,r14\t\n"
//       "jz skip_preload_b_zero%=\t\n"
//       "vcvtph2ps zmm31,YMMWORD PTR [r10 + 64]\t\n"
//       "skip_preload_b_zero%=:\t\n"
//       "vbroadcastss zmm2,DWORD PTR [r9+0]\t\n"
//       "vmulps zmm0,zmm3,zmm2\t\n"
//       "vmulps zmm1,zmm4,zmm2\t\n"
//       "test r14,r14\t\n"
//       "jnz next_inner%=\t\n"
//       "add r10,64\t\n"
//       "jmp dump_C%=\t\n"

//       "loop_inner%=:\t\n"

//       "vmovaps zmm3,zmm31\t\n"
//       "vcvtph2ps zmm4,YMMWORD PTR [r10 + 32]\t\n"
//       "vcvtph2ps zmm31,YMMWORD PTR [r10 + 64]\t\n"
//       "vbroadcastss zmm2,DWORD PTR [r9+0]\t\n"
//       "vfmadd231ps zmm0,zmm3,zmm2\t\n"
//       "vfmadd231ps zmm1,zmm4,zmm2\t\n"

//       "next_inner%=:\t\n"
//       "add r9,4\t\n"
//       "add r10,64\t\n"
//       "dec r14\t\n"
//       "jnz loop_inner%=\t\n"

//       "vmovaps zmm3,zmm31\t\n"
//       "vcvtph2ps zmm4,YMMWORD PTR [r10 + 32]\t\n"
//       "vbroadcastss zmm2,DWORD PTR [r9+0]\t\n"
//       "vfmadd231ps zmm0,zmm3,zmm2\t\n"
//       "vfmadd231ps zmm1,zmm4,zmm2\t\n"
//       "add r9,4\t\n"
//       "add r10,64\t\n"
//       // Dump C
//       "dump_C%=:\t\n"
//       "vmovups zmmword PTR [r12 + 0], zmm0\t\n"
//       "vmovups zmmword PTR [r12 + 64], zmm1\t\n"

//       // next outer iteration
//       "add rcx, 128\t\n"
//       "mov r12, rcx\t\n"
//       "mov r9, rax\t\n"
//       "inc rbx\t\n"
//       "cmp rbx, rdi\t\n"
//       "jl loop_outter%=\t\n"
//       :
//       : [gp] "rm"(gp)
//       : "r8",
//         "r9",
//         "r10",
//         "r11",
//         "r13",
//         "r14",
//         "rax",
//         "rcx",
//         "rsi",
//         "rdi",
//         "rbx",
//         "r12",
//         "r15",
//         "memory");

 gemmkernel_Avx512_fp16_fA0fB0fC0(gp, 1);
}
void NOINLINE
gemmkernel_2x2_Avx512_fp16_fA0fB0fC0(GemmParamsFP16* gp) {
//   asm volatile(
// #if !defined(__clang__)
//       "mov r14, %[gp]\t\n"
// #else
//       "mov %[gp], %%r14\t\n"
//       ".intel_syntax noprefix\t\n"
// #endif

//       // Copy parameters
//       // k
//       "mov r8, [r14 + 0]\t\n"
//       "dec r8\t\n"
//       // A
//       "mov r9, [r14 + 8]\t\n"
//       // B
//       "mov r10, [r14 + 16]\t\n"
//       // beta
//       "lea r15, [r14 + 24]\t\n"
//       // C
//       "mov r12, [r14 + 32]\t\n"
//       // ldc
//       "mov r13, [r14 + 40]\t\n"
//       // b_block_cols
//       "mov rdi, [r14 + 48]\t\n"
//       // b_block_size
//       "mov rsi, [r14 + 56]\t\n"

//       // Make copies of A and C
//       "mov rax, r9\t\n"
//       "mov rcx, r12\t\n"

//       "mov rbx, 0\t\n"
//       "loop_outter%=:\t\n"
//       "mov r14, r8\t\n"
//       "vbroadcastss zmm31,DWORD PTR [r15]\t\n"
//       "vcvtph2ps zmm5,YMMWORD PTR [r10 + 0]\t\n"
//       "vcvtph2ps zmm6,YMMWORD PTR [r10 + 32]\t\n"
//       "vxorps xmm0, xmm0, xmm0\t\n"
//       "vcomiss xmm31, xmm0\t\n"
//       "jz zero_regs%=\t\n"

//       // Setup values with beta multiplication
//       "vmulps zmm0, zmm31, [r12 + 0]\t\n"
//       "vmulps zmm1, zmm31, [r12 + 64]\t\n"
//       "add r12, r13\t\n"
//       "vmulps zmm2, zmm31, [r12 + 0]\t\n"
//       "vmulps zmm3, zmm31, [r12 + 64]\t\n"
//       "test r14,r14\t\n"
//       "jz skip_preload%=\t\n"
//       "vcvtph2ps zmm31,YMMWORD PTR [r10 + 64]\t\n"
//       "skip_preload%=:\t\n"
//       "vbroadcastss zmm4,DWORD PTR [r9+0]\t\n"
//       "vfmadd231ps zmm0,zmm5,zmm4\t\n"
//       "vfmadd231ps zmm1,zmm6,zmm4\t\n"
//       "vbroadcastss zmm4,DWORD PTR [r9+4]\t\n"
//       "vfmadd231ps zmm2,zmm5,zmm4\t\n"
//       "vfmadd231ps zmm3,zmm6,zmm4\t\n"
//       "mov r12, rcx\t\n"
//       "test r14,r14\t\n"
//       "jnz next_inner%=\t\n"
//       "add r10,64\t\n"
//       "jmp dump_C%=\t\n"

//       "zero_regs%=:\t\n"

//       "test r14,r14\t\n"
//       "jz skip_preload_b_zero%=\t\n"
//       "vcvtph2ps zmm31,YMMWORD PTR [r10 + 64]\t\n"
//       "skip_preload_b_zero%=:\t\n"
//       "vbroadcastss zmm4,DWORD PTR [r9+0]\t\n"
//       "vmulps zmm0,zmm5,zmm4\t\n"
//       "vmulps zmm1,zmm6,zmm4\t\n"
//       "add r12, r13\t\n"
//       "vbroadcastss zmm4,DWORD PTR [r9+4]\t\n"
//       "vmulps zmm2,zmm5,zmm4\t\n"
//       "vmulps zmm3,zmm6,zmm4\t\n"
//       "mov r12, rcx\t\n"
//       "test r14,r14\t\n"
//       "jnz next_inner%=\t\n"
//       "add r10,64\t\n"
//       "jmp dump_C%=\t\n"

//       "loop_inner%=:\t\n"

//       "vmovaps zmm5,zmm31\t\n"
//       "vcvtph2ps zmm6,YMMWORD PTR [r10 + 32]\t\n"
//       "vcvtph2ps zmm31,YMMWORD PTR [r10 + 64]\t\n"
//       "vbroadcastss zmm4,DWORD PTR [r9+0]\t\n"
//       "vfmadd231ps zmm0,zmm5,zmm4\t\n"
//       "vfmadd231ps zmm1,zmm6,zmm4\t\n"
//       "vbroadcastss zmm4,DWORD PTR [r9+4]\t\n"
//       "vfmadd231ps zmm2,zmm5,zmm4\t\n"
//       "vfmadd231ps zmm3,zmm6,zmm4\t\n"

//       "next_inner%=:\t\n"
//       "add r9,8\t\n"
//       "add r10,64\t\n"
//       "dec r14\t\n"
//       "jnz loop_inner%=\t\n"

//       "vmovaps zmm5,zmm31\t\n"
//       "vcvtph2ps zmm6,YMMWORD PTR [r10 + 32]\t\n"
//       "vbroadcastss zmm4,DWORD PTR [r9+0]\t\n"
//       "vfmadd231ps zmm0,zmm5,zmm4\t\n"
//       "vfmadd231ps zmm1,zmm6,zmm4\t\n"
//       "vbroadcastss zmm4,DWORD PTR [r9+4]\t\n"
//       "vfmadd231ps zmm2,zmm5,zmm4\t\n"
//       "vfmadd231ps zmm3,zmm6,zmm4\t\n"
//       "add r9,8\t\n"
//       "add r10,64\t\n"
//       // Dump C
//       "dump_C%=:\t\n"
//       "vmovups zmmword PTR [r12 + 0], zmm0\t\n"
//       "vmovups zmmword PTR [r12 + 64], zmm1\t\n"
//       "add r12, r13\t\n"
//       "vmovups zmmword PTR [r12 + 0], zmm2\t\n"
//       "vmovups zmmword PTR [r12 + 64], zmm3\t\n"

//       // next outer iteration
//       "add rcx, 128\t\n"
//       "mov r12, rcx\t\n"
//       "mov r9, rax\t\n"
//       "inc rbx\t\n"
//       "cmp rbx, rdi\t\n"
//       "jl loop_outter%=\t\n"
//       :
//       : [gp] "rm"(gp)
//       : "r8",
//         "r9",
//         "r10",
//         "r11",
//         "r13",
//         "r14",
//         "rax",
//         "rcx",
//         "rsi",
//         "rdi",
//         "rbx",
//         "r12",
//         "r15",
//         "memory");

 gemmkernel_Avx512_fp16_fA0fB0fC0(gp, 2);
}
void NOINLINE
gemmkernel_3x2_Avx512_fp16_fA0fB0fC0(GemmParamsFP16* gp) {
  asm volatile(
#if !defined(__clang__)
      "mov r14, %[gp]\t\n"
#else
      "mov %[gp], %%r14\t\n"
      ".intel_syntax noprefix\t\n"
#endif

      // Copy parameters
      // k
      "mov r8, [r14 + 0]\t\n"
      "dec r8\t\n"
      // A
      "mov r9, [r14 + 8]\t\n"
      // B
      "mov r10, [r14 + 16]\t\n"
      // beta
      "lea r15, [r14 + 24]\t\n"
      // C
      "mov r12, [r14 + 32]\t\n"
      // ldc
      "mov r13, [r14 + 40]\t\n"
      // b_block_cols
      "mov rdi, [r14 + 48]\t\n"
      // b_block_size
      "mov rsi, [r14 + 56]\t\n"

      // Make copies of A and C
      "mov rax, r9\t\n"
      "mov rcx, r12\t\n"

      "mov rbx, 0\t\n"
      "loop_outter%=:\t\n"
      "mov r14, r8\t\n"
      "vbroadcastss zmm31,DWORD PTR [r15]\t\n"
      "vcvtph2ps zmm7,YMMWORD PTR [r10 + 0]\t\n"
      "vcvtph2ps zmm8,YMMWORD PTR [r10 + 32]\t\n"
      "vxorps xmm0, xmm0, xmm0\t\n"
      "vcomiss xmm31, xmm0\t\n"
      "jz zero_regs%=\t\n"

      // Setup values with beta multiplication
      "vmulps zmm0, zmm31, [r12 + 0]\t\n"
      "vmulps zmm1, zmm31, [r12 + 64]\t\n"
      "add r12, r13\t\n"
      "vmulps zmm2, zmm31, [r12 + 0]\t\n"
      "vmulps zmm3, zmm31, [r12 + 64]\t\n"
      "add r12, r13\t\n"
      "vmulps zmm4, zmm31, [r12 + 0]\t\n"
      "vmulps zmm5, zmm31, [r12 + 64]\t\n"
      "test r14,r14\t\n"
      "jz skip_preload%=\t\n"
      "vcvtph2ps zmm31,YMMWORD PTR [r10 + 64]\t\n"
      "skip_preload%=:\t\n"
      "vbroadcastss zmm6,DWORD PTR [r9+0]\t\n"
      "vfmadd231ps zmm0,zmm7,zmm6\t\n"
      "vfmadd231ps zmm1,zmm8,zmm6\t\n"
      "vbroadcastss zmm6,DWORD PTR [r9+4]\t\n"
      "vfmadd231ps zmm2,zmm7,zmm6\t\n"
      "vfmadd231ps zmm3,zmm8,zmm6\t\n"
      "vbroadcastss zmm6,DWORD PTR [r9+8]\t\n"
      "vfmadd231ps zmm4,zmm7,zmm6\t\n"
      "vfmadd231ps zmm5,zmm8,zmm6\t\n"
      "mov r12, rcx\t\n"
      "test r14,r14\t\n"
      "jnz next_inner%=\t\n"
      "add r10,64\t\n"
      "jmp dump_C%=\t\n"

      "zero_regs%=:\t\n"

      "test r14,r14\t\n"
      "jz skip_preload_b_zero%=\t\n"
      "vcvtph2ps zmm31,YMMWORD PTR [r10 + 64]\t\n"
      "skip_preload_b_zero%=:\t\n"
      "vbroadcastss zmm6,DWORD PTR [r9+0]\t\n"
      "vmulps zmm0,zmm7,zmm6\t\n"
      "vmulps zmm1,zmm8,zmm6\t\n"
      "add r12, r13\t\n"
      "vbroadcastss zmm6,DWORD PTR [r9+4]\t\n"
      "vmulps zmm2,zmm7,zmm6\t\n"
      "vmulps zmm3,zmm8,zmm6\t\n"
      "add r12, r13\t\n"
      "vbroadcastss zmm6,DWORD PTR [r9+8]\t\n"
      "vmulps zmm4,zmm7,zmm6\t\n"
      "vmulps zmm5,zmm8,zmm6\t\n"
      "mov r12, rcx\t\n"
      "test r14,r14\t\n"
      "jnz next_inner%=\t\n"
      "add r10,64\t\n"
      "jmp dump_C%=\t\n"

      "loop_inner%=:\t\n"

      "vmovaps zmm7,zmm31\t\n"
      "vcvtph2ps zmm8,YMMWORD PTR [r10 + 32]\t\n"
      "vcvtph2ps zmm31,YMMWORD PTR [r10 + 64]\t\n"
      "vbroadcastss zmm6,DWORD PTR [r9+0]\t\n"
      "vfmadd231ps zmm0,zmm7,zmm6\t\n"
      "vfmadd231ps zmm1,zmm8,zmm6\t\n"
      "vbroadcastss zmm6,DWORD PTR [r9+4]\t\n"
      "vfmadd231ps zmm2,zmm7,zmm6\t\n"
      "vfmadd231ps zmm3,zmm8,zmm6\t\n"
      "vbroadcastss zmm6,DWORD PTR [r9+8]\t\n"
      "vfmadd231ps zmm4,zmm7,zmm6\t\n"
      "vfmadd231ps zmm5,zmm8,zmm6\t\n"

      "next_inner%=:\t\n"
      "add r9,12\t\n"
      "add r10,64\t\n"
      "dec r14\t\n"
      "jnz loop_inner%=\t\n"

      "vmovaps zmm7,zmm31\t\n"
      "vcvtph2ps zmm8,YMMWORD PTR [r10 + 32]\t\n"
      "vbroadcastss zmm6,DWORD PTR [r9+0]\t\n"
      "vfmadd231ps zmm0,zmm7,zmm6\t\n"
      "vfmadd231ps zmm1,zmm8,zmm6\t\n"
      "vbroadcastss zmm6,DWORD PTR [r9+4]\t\n"
      "vfmadd231ps zmm2,zmm7,zmm6\t\n"
      "vfmadd231ps zmm3,zmm8,zmm6\t\n"
      "vbroadcastss zmm6,DWORD PTR [r9+8]\t\n"
      "vfmadd231ps zmm4,zmm7,zmm6\t\n"
      "vfmadd231ps zmm5,zmm8,zmm6\t\n"
      "add r9,12\t\n"
      "add r10,64\t\n"
      // Dump C
      "dump_C%=:\t\n"
      "vmovups zmmword PTR [r12 + 0], zmm0\t\n"
      "vmovups zmmword PTR [r12 + 64], zmm1\t\n"
      "add r12, r13\t\n"
      "vmovups zmmword PTR [r12 + 0], zmm2\t\n"
      "vmovups zmmword PTR [r12 + 64], zmm3\t\n"
      "add r12, r13\t\n"
      "vmovups zmmword PTR [r12 + 0], zmm4\t\n"
      "vmovups zmmword PTR [r12 + 64], zmm5\t\n"

      // next outer iteration
      "add rcx, 128\t\n"
      "mov r12, rcx\t\n"
      "mov r9, rax\t\n"
      "inc rbx\t\n"
      "cmp rbx, rdi\t\n"
      "jl loop_outter%=\t\n"
      :
      : [gp] "rm"(gp)
      : "r8",
        "r9",
        "r10",
        "r11",
        "r13",
        "r14",
        "rax",
        "rcx",
        "rsi",
        "rdi",
        "rbx",
        "r12",
        "r15",
        "memory");
}
void NOINLINE
gemmkernel_4x2_Avx512_fp16_fA0fB0fC0(GemmParamsFP16* gp) {
  asm volatile(
#if !defined(__clang__)
      "mov r14, %[gp]\t\n"
#else
      "mov %[gp], %%r14\t\n"
      ".intel_syntax noprefix\t\n"
#endif

      // Copy parameters
      // k
      "mov r8, [r14 + 0]\t\n"
      "dec r8\t\n"
      // A
      "mov r9, [r14 + 8]\t\n"
      // B
      "mov r10, [r14 + 16]\t\n"
      // beta
      "lea r15, [r14 + 24]\t\n"
      // C
      "mov r12, [r14 + 32]\t\n"
      // ldc
      "mov r13, [r14 + 40]\t\n"
      // b_block_cols
      "mov rdi, [r14 + 48]\t\n"
      // b_block_size
      "mov rsi, [r14 + 56]\t\n"

      // Make copies of A and C
      "mov rax, r9\t\n"
      "mov rcx, r12\t\n"

      "mov rbx, 0\t\n"
      "loop_outter%=:\t\n"
      "mov r14, r8\t\n"
      "vbroadcastss zmm31,DWORD PTR [r15]\t\n"
      "vcvtph2ps zmm9,YMMWORD PTR [r10 + 0]\t\n"
      "vcvtph2ps zmm10,YMMWORD PTR [r10 + 32]\t\n"
      "vxorps xmm0, xmm0, xmm0\t\n"
      "vcomiss xmm31, xmm0\t\n"
      "jz zero_regs%=\t\n"

      // Setup values with beta multiplication
      "vmulps zmm0, zmm31, [r12 + 0]\t\n"
      "vmulps zmm1, zmm31, [r12 + 64]\t\n"
      "add r12, r13\t\n"
      "vmulps zmm2, zmm31, [r12 + 0]\t\n"
      "vmulps zmm3, zmm31, [r12 + 64]\t\n"
      "add r12, r13\t\n"
      "vmulps zmm4, zmm31, [r12 + 0]\t\n"
      "vmulps zmm5, zmm31, [r12 + 64]\t\n"
      "add r12, r13\t\n"
      "vmulps zmm6, zmm31, [r12 + 0]\t\n"
      "vmulps zmm7, zmm31, [r12 + 64]\t\n"
      "test r14,r14\t\n"
      "jz skip_preload%=\t\n"
      "vcvtph2ps zmm31,YMMWORD PTR [r10 + 64]\t\n"
      "skip_preload%=:\t\n"
      "vbroadcastss zmm8,DWORD PTR [r9+0]\t\n"
      "vfmadd231ps zmm0,zmm9,zmm8\t\n"
      "vfmadd231ps zmm1,zmm10,zmm8\t\n"
      "vbroadcastss zmm8,DWORD PTR [r9+4]\t\n"
      "vfmadd231ps zmm2,zmm9,zmm8\t\n"
      "vfmadd231ps zmm3,zmm10,zmm8\t\n"
      "vbroadcastss zmm8,DWORD PTR [r9+8]\t\n"
      "vfmadd231ps zmm4,zmm9,zmm8\t\n"
      "vfmadd231ps zmm5,zmm10,zmm8\t\n"
      "vbroadcastss zmm8,DWORD PTR [r9+12]\t\n"
      "vfmadd231ps zmm6,zmm9,zmm8\t\n"
      "vfmadd231ps zmm7,zmm10,zmm8\t\n"
      "mov r12, rcx\t\n"
      "test r14,r14\t\n"
      "jnz next_inner%=\t\n"
      "add r10,64\t\n"
      "jmp dump_C%=\t\n"

      "zero_regs%=:\t\n"

      "test r14,r14\t\n"
      "jz skip_preload_b_zero%=\t\n"
      "vcvtph2ps zmm31,YMMWORD PTR [r10 + 64]\t\n"
      "skip_preload_b_zero%=:\t\n"
      "vbroadcastss zmm8,DWORD PTR [r9+0]\t\n"
      "vmulps zmm0,zmm9,zmm8\t\n"
      "vmulps zmm1,zmm10,zmm8\t\n"
      "add r12, r13\t\n"
      "vbroadcastss zmm8,DWORD PTR [r9+4]\t\n"
      "vmulps zmm2,zmm9,zmm8\t\n"
      "vmulps zmm3,zmm10,zmm8\t\n"
      "add r12, r13\t\n"
      "vbroadcastss zmm8,DWORD PTR [r9+8]\t\n"
      "vmulps zmm4,zmm9,zmm8\t\n"
      "vmulps zmm5,zmm10,zmm8\t\n"
      "add r12, r13\t\n"
      "vbroadcastss zmm8,DWORD PTR [r9+12]\t\n"
      "vmulps zmm6,zmm9,zmm8\t\n"
      "vmulps zmm7,zmm10,zmm8\t\n"
      "mov r12, rcx\t\n"
      "test r14,r14\t\n"
      "jnz next_inner%=\t\n"
      "add r10,64\t\n"
      "jmp dump_C%=\t\n"

      "loop_inner%=:\t\n"

      "vmovaps zmm9,zmm31\t\n"
      "vcvtph2ps zmm10,YMMWORD PTR [r10 + 32]\t\n"
      "vcvtph2ps zmm31,YMMWORD PTR [r10 + 64]\t\n"
      "vbroadcastss zmm8,DWORD PTR [r9+0]\t\n"
      "vfmadd231ps zmm0,zmm9,zmm8\t\n"
      "vfmadd231ps zmm1,zmm10,zmm8\t\n"
      "vbroadcastss zmm8,DWORD PTR [r9+4]\t\n"
      "vfmadd231ps zmm2,zmm9,zmm8\t\n"
      "vfmadd231ps zmm3,zmm10,zmm8\t\n"
      "vbroadcastss zmm8,DWORD PTR [r9+8]\t\n"
      "vfmadd231ps zmm4,zmm9,zmm8\t\n"
      "vfmadd231ps zmm5,zmm10,zmm8\t\n"
      "vbroadcastss zmm8,DWORD PTR [r9+12]\t\n"
      "vfmadd231ps zmm6,zmm9,zmm8\t\n"
      "vfmadd231ps zmm7,zmm10,zmm8\t\n"

      "next_inner%=:\t\n"
      "add r9,16\t\n"
      "add r10,64\t\n"
      "dec r14\t\n"
      "jnz loop_inner%=\t\n"

      "vmovaps zmm9,zmm31\t\n"
      "vcvtph2ps zmm10,YMMWORD PTR [r10 + 32]\t\n"
      "vbroadcastss zmm8,DWORD PTR [r9+0]\t\n"
      "vfmadd231ps zmm0,zmm9,zmm8\t\n"
      "vfmadd231ps zmm1,zmm10,zmm8\t\n"
      "vbroadcastss zmm8,DWORD PTR [r9+4]\t\n"
      "vfmadd231ps zmm2,zmm9,zmm8\t\n"
      "vfmadd231ps zmm3,zmm10,zmm8\t\n"
      "vbroadcastss zmm8,DWORD PTR [r9+8]\t\n"
      "vfmadd231ps zmm4,zmm9,zmm8\t\n"
      "vfmadd231ps zmm5,zmm10,zmm8\t\n"
      "vbroadcastss zmm8,DWORD PTR [r9+12]\t\n"
      "vfmadd231ps zmm6,zmm9,zmm8\t\n"
      "vfmadd231ps zmm7,zmm10,zmm8\t\n"
      "add r9,16\t\n"
      "add r10,64\t\n"
      // Dump C
      "dump_C%=:\t\n"
      "vmovups zmmword PTR [r12 + 0], zmm0\t\n"
      "vmovups zmmword PTR [r12 + 64], zmm1\t\n"
      "add r12, r13\t\n"
      "vmovups zmmword PTR [r12 + 0], zmm2\t\n"
      "vmovups zmmword PTR [r12 + 64], zmm3\t\n"
      "add r12, r13\t\n"
      "vmovups zmmword PTR [r12 + 0], zmm4\t\n"
      "vmovups zmmword PTR [r12 + 64], zmm5\t\n"
      "add r12, r13\t\n"
      "vmovups zmmword PTR [r12 + 0], zmm6\t\n"
      "vmovups zmmword PTR [r12 + 64], zmm7\t\n"

      // next outer iteration
      "add rcx, 128\t\n"
      "mov r12, rcx\t\n"
      "mov r9, rax\t\n"
      "inc rbx\t\n"
      "cmp rbx, rdi\t\n"
      "jl loop_outter%=\t\n"
      :
      : [gp] "rm"(gp)
      : "r8",
        "r9",
        "r10",
        "r11",
        "r13",
        "r14",
        "rax",
        "rcx",
        "rsi",
        "rdi",
        "rbx",
        "r12",
        "r15",
        "memory");
}
void NOINLINE
gemmkernel_5x2_Avx512_fp16_fA0fB0fC0(GemmParamsFP16* gp) {
  asm volatile(
#if !defined(__clang__)
      "mov r14, %[gp]\t\n"
#else
      "mov %[gp], %%r14\t\n"
      ".intel_syntax noprefix\t\n"
#endif

      // Copy parameters
      // k
      "mov r8, [r14 + 0]\t\n"
      "dec r8\t\n"
      // A
      "mov r9, [r14 + 8]\t\n"
      // B
      "mov r10, [r14 + 16]\t\n"
      // beta
      "lea r15, [r14 + 24]\t\n"
      // C
      "mov r12, [r14 + 32]\t\n"
      // ldc
      "mov r13, [r14 + 40]\t\n"
      // b_block_cols
      "mov rdi, [r14 + 48]\t\n"
      // b_block_size
      "mov rsi, [r14 + 56]\t\n"

      // Make copies of A and C
      "mov rax, r9\t\n"
      "mov rcx, r12\t\n"

      "mov rbx, 0\t\n"
      "loop_outter%=:\t\n"
      "mov r14, r8\t\n"
      "vbroadcastss zmm31,DWORD PTR [r15]\t\n"
      "vcvtph2ps zmm11,YMMWORD PTR [r10 + 0]\t\n"
      "vcvtph2ps zmm12,YMMWORD PTR [r10 + 32]\t\n"
      "vxorps xmm0, xmm0, xmm0\t\n"
      "vcomiss xmm31, xmm0\t\n"
      "jz zero_regs%=\t\n"

      // Setup values with beta multiplication
      "vmulps zmm0, zmm31, [r12 + 0]\t\n"
      "vmulps zmm1, zmm31, [r12 + 64]\t\n"
      "add r12, r13\t\n"
      "vmulps zmm2, zmm31, [r12 + 0]\t\n"
      "vmulps zmm3, zmm31, [r12 + 64]\t\n"
      "add r12, r13\t\n"
      "vmulps zmm4, zmm31, [r12 + 0]\t\n"
      "vmulps zmm5, zmm31, [r12 + 64]\t\n"
      "add r12, r13\t\n"
      "vmulps zmm6, zmm31, [r12 + 0]\t\n"
      "vmulps zmm7, zmm31, [r12 + 64]\t\n"
      "add r12, r13\t\n"
      "vmulps zmm8, zmm31, [r12 + 0]\t\n"
      "vmulps zmm9, zmm31, [r12 + 64]\t\n"
      "test r14,r14\t\n"
      "jz skip_preload%=\t\n"
      "vcvtph2ps zmm31,YMMWORD PTR [r10 + 64]\t\n"
      "skip_preload%=:\t\n"
      "vbroadcastss zmm10,DWORD PTR [r9+0]\t\n"
      "vfmadd231ps zmm0,zmm11,zmm10\t\n"
      "vfmadd231ps zmm1,zmm12,zmm10\t\n"
      "vbroadcastss zmm10,DWORD PTR [r9+4]\t\n"
      "vfmadd231ps zmm2,zmm11,zmm10\t\n"
      "vfmadd231ps zmm3,zmm12,zmm10\t\n"
      "vbroadcastss zmm10,DWORD PTR [r9+8]\t\n"
      "vfmadd231ps zmm4,zmm11,zmm10\t\n"
      "vfmadd231ps zmm5,zmm12,zmm10\t\n"
      "vbroadcastss zmm10,DWORD PTR [r9+12]\t\n"
      "vfmadd231ps zmm6,zmm11,zmm10\t\n"
      "vfmadd231ps zmm7,zmm12,zmm10\t\n"
      "vbroadcastss zmm10,DWORD PTR [r9+16]\t\n"
      "vfmadd231ps zmm8,zmm11,zmm10\t\n"
      "vfmadd231ps zmm9,zmm12,zmm10\t\n"
      "mov r12, rcx\t\n"
      "test r14,r14\t\n"
      "jnz next_inner%=\t\n"
      "add r10,64\t\n"
      "jmp dump_C%=\t\n"

      "zero_regs%=:\t\n"

      "test r14,r14\t\n"
      "jz skip_preload_b_zero%=\t\n"
      "vcvtph2ps zmm31,YMMWORD PTR [r10 + 64]\t\n"
      "skip_preload_b_zero%=:\t\n"
      "vbroadcastss zmm10,DWORD PTR [r9+0]\t\n"
      "vmulps zmm0,zmm11,zmm10\t\n"
      "vmulps zmm1,zmm12,zmm10\t\n"
      "add r12, r13\t\n"
      "vbroadcastss zmm10,DWORD PTR [r9+4]\t\n"
      "vmulps zmm2,zmm11,zmm10\t\n"
      "vmulps zmm3,zmm12,zmm10\t\n"
      "add r12, r13\t\n"
      "vbroadcastss zmm10,DWORD PTR [r9+8]\t\n"
      "vmulps zmm4,zmm11,zmm10\t\n"
      "vmulps zmm5,zmm12,zmm10\t\n"
      "add r12, r13\t\n"
      "vbroadcastss zmm10,DWORD PTR [r9+12]\t\n"
      "vmulps zmm6,zmm11,zmm10\t\n"
      "vmulps zmm7,zmm12,zmm10\t\n"
      "add r12, r13\t\n"
      "vbroadcastss zmm10,DWORD PTR [r9+16]\t\n"
      "vmulps zmm8,zmm11,zmm10\t\n"
      "vmulps zmm9,zmm12,zmm10\t\n"
      "mov r12, rcx\t\n"
      "test r14,r14\t\n"
      "jnz next_inner%=\t\n"
      "add r10,64\t\n"
      "jmp dump_C%=\t\n"

      "loop_inner%=:\t\n"

      "vmovaps zmm11,zmm31\t\n"
      "vcvtph2ps zmm12,YMMWORD PTR [r10 + 32]\t\n"
      "vcvtph2ps zmm31,YMMWORD PTR [r10 + 64]\t\n"
      "vbroadcastss zmm10,DWORD PTR [r9+0]\t\n"
      "vfmadd231ps zmm0,zmm11,zmm10\t\n"
      "vfmadd231ps zmm1,zmm12,zmm10\t\n"
      "vbroadcastss zmm10,DWORD PTR [r9+4]\t\n"
      "vfmadd231ps zmm2,zmm11,zmm10\t\n"
      "vfmadd231ps zmm3,zmm12,zmm10\t\n"
      "vbroadcastss zmm10,DWORD PTR [r9+8]\t\n"
      "vfmadd231ps zmm4,zmm11,zmm10\t\n"
      "vfmadd231ps zmm5,zmm12,zmm10\t\n"
      "vbroadcastss zmm10,DWORD PTR [r9+12]\t\n"
      "vfmadd231ps zmm6,zmm11,zmm10\t\n"
      "vfmadd231ps zmm7,zmm12,zmm10\t\n"
      "vbroadcastss zmm10,DWORD PTR [r9+16]\t\n"
      "vfmadd231ps zmm8,zmm11,zmm10\t\n"
      "vfmadd231ps zmm9,zmm12,zmm10\t\n"

      "next_inner%=:\t\n"
      "add r9,20\t\n"
      "add r10,64\t\n"
      "dec r14\t\n"
      "jnz loop_inner%=\t\n"

      "vmovaps zmm11,zmm31\t\n"
      "vcvtph2ps zmm12,YMMWORD PTR [r10 + 32]\t\n"
      "vbroadcastss zmm10,DWORD PTR [r9+0]\t\n"
      "vfmadd231ps zmm0,zmm11,zmm10\t\n"
      "vfmadd231ps zmm1,zmm12,zmm10\t\n"
      "vbroadcastss zmm10,DWORD PTR [r9+4]\t\n"
      "vfmadd231ps zmm2,zmm11,zmm10\t\n"
      "vfmadd231ps zmm3,zmm12,zmm10\t\n"
      "vbroadcastss zmm10,DWORD PTR [r9+8]\t\n"
      "vfmadd231ps zmm4,zmm11,zmm10\t\n"
      "vfmadd231ps zmm5,zmm12,zmm10\t\n"
      "vbroadcastss zmm10,DWORD PTR [r9+12]\t\n"
      "vfmadd231ps zmm6,zmm11,zmm10\t\n"
      "vfmadd231ps zmm7,zmm12,zmm10\t\n"
      "vbroadcastss zmm10,DWORD PTR [r9+16]\t\n"
      "vfmadd231ps zmm8,zmm11,zmm10\t\n"
      "vfmadd231ps zmm9,zmm12,zmm10\t\n"
      "add r9,20\t\n"
      "add r10,64\t\n"
      // Dump C
      "dump_C%=:\t\n"
      "vmovups zmmword PTR [r12 + 0], zmm0\t\n"
      "vmovups zmmword PTR [r12 + 64], zmm1\t\n"
      "add r12, r13\t\n"
      "vmovups zmmword PTR [r12 + 0], zmm2\t\n"
      "vmovups zmmword PTR [r12 + 64], zmm3\t\n"
      "add r12, r13\t\n"
      "vmovups zmmword PTR [r12 + 0], zmm4\t\n"
      "vmovups zmmword PTR [r12 + 64], zmm5\t\n"
      "add r12, r13\t\n"
      "vmovups zmmword PTR [r12 + 0], zmm6\t\n"
      "vmovups zmmword PTR [r12 + 64], zmm7\t\n"
      "add r12, r13\t\n"
      "vmovups zmmword PTR [r12 + 0], zmm8\t\n"
      "vmovups zmmword PTR [r12 + 64], zmm9\t\n"

      // next outer iteration
      "add rcx, 128\t\n"
      "mov r12, rcx\t\n"
      "mov r9, rax\t\n"
      "inc rbx\t\n"
      "cmp rbx, rdi\t\n"
      "jl loop_outter%=\t\n"
      :
      : [gp] "rm"(gp)
      : "r8",
        "r9",
        "r10",
        "r11",
        "r13",
        "r14",
        "rax",
        "rcx",
        "rsi",
        "rdi",
        "rbx",
        "r12",
        "r15",
        "memory");
}
void NOINLINE
gemmkernel_6x2_Avx512_fp16_fA0fB0fC0(GemmParamsFP16* gp) {
  asm volatile(
#if !defined(__clang__)
      "mov r14, %[gp]\t\n"
#else
      "mov %[gp], %%r14\t\n"
      ".intel_syntax noprefix\t\n"
#endif

      // Copy parameters
      // k
      "mov r8, [r14 + 0]\t\n"
      "dec r8\t\n"
      // A
      "mov r9, [r14 + 8]\t\n"
      // B
      "mov r10, [r14 + 16]\t\n"
      // beta
      "lea r15, [r14 + 24]\t\n"
      // C
      "mov r12, [r14 + 32]\t\n"
      // ldc
      "mov r13, [r14 + 40]\t\n"
      // b_block_cols
      "mov rdi, [r14 + 48]\t\n"
      // b_block_size
      "mov rsi, [r14 + 56]\t\n"

      // Make copies of A and C
      "mov rax, r9\t\n"
      "mov rcx, r12\t\n"

      "mov rbx, 0\t\n"
      "loop_outter%=:\t\n"
      "mov r14, r8\t\n"
      "vbroadcastss zmm31,DWORD PTR [r15]\t\n"
      "vcvtph2ps zmm13,YMMWORD PTR [r10 + 0]\t\n"
      "vcvtph2ps zmm14,YMMWORD PTR [r10 + 32]\t\n"
      "vxorps xmm0, xmm0, xmm0\t\n"
      "vcomiss xmm31, xmm0\t\n"
      "jz zero_regs%=\t\n"

      // Setup values with beta multiplication
      "vmulps zmm0, zmm31, [r12 + 0]\t\n"
      "vmulps zmm1, zmm31, [r12 + 64]\t\n"
      "add r12, r13\t\n"
      "vmulps zmm2, zmm31, [r12 + 0]\t\n"
      "vmulps zmm3, zmm31, [r12 + 64]\t\n"
      "add r12, r13\t\n"
      "vmulps zmm4, zmm31, [r12 + 0]\t\n"
      "vmulps zmm5, zmm31, [r12 + 64]\t\n"
      "add r12, r13\t\n"
      "vmulps zmm6, zmm31, [r12 + 0]\t\n"
      "vmulps zmm7, zmm31, [r12 + 64]\t\n"
      "add r12, r13\t\n"
      "vmulps zmm8, zmm31, [r12 + 0]\t\n"
      "vmulps zmm9, zmm31, [r12 + 64]\t\n"
      "add r12, r13\t\n"
      "vmulps zmm10, zmm31, [r12 + 0]\t\n"
      "vmulps zmm11, zmm31, [r12 + 64]\t\n"
      "test r14,r14\t\n"
      "jz skip_preload%=\t\n"
      "vcvtph2ps zmm31,YMMWORD PTR [r10 + 64]\t\n"
      "skip_preload%=:\t\n"
      "vbroadcastss zmm12,DWORD PTR [r9+0]\t\n"
      "vfmadd231ps zmm0,zmm13,zmm12\t\n"
      "vfmadd231ps zmm1,zmm14,zmm12\t\n"
      "vbroadcastss zmm12,DWORD PTR [r9+4]\t\n"
      "vfmadd231ps zmm2,zmm13,zmm12\t\n"
      "vfmadd231ps zmm3,zmm14,zmm12\t\n"
      "vbroadcastss zmm12,DWORD PTR [r9+8]\t\n"
      "vfmadd231ps zmm4,zmm13,zmm12\t\n"
      "vfmadd231ps zmm5,zmm14,zmm12\t\n"
      "vbroadcastss zmm12,DWORD PTR [r9+12]\t\n"
      "vfmadd231ps zmm6,zmm13,zmm12\t\n"
      "vfmadd231ps zmm7,zmm14,zmm12\t\n"
      "vbroadcastss zmm12,DWORD PTR [r9+16]\t\n"
      "vfmadd231ps zmm8,zmm13,zmm12\t\n"
      "vfmadd231ps zmm9,zmm14,zmm12\t\n"
      "vbroadcastss zmm12,DWORD PTR [r9+20]\t\n"
      "vfmadd231ps zmm10,zmm13,zmm12\t\n"
      "vfmadd231ps zmm11,zmm14,zmm12\t\n"
      "mov r12, rcx\t\n"
      "test r14,r14\t\n"
      "jnz next_inner%=\t\n"
      "add r10,64\t\n"
      "jmp dump_C%=\t\n"

      "zero_regs%=:\t\n"

      "test r14,r14\t\n"
      "jz skip_preload_b_zero%=\t\n"
      "vcvtph2ps zmm31,YMMWORD PTR [r10 + 64]\t\n"
      "skip_preload_b_zero%=:\t\n"
      "vbroadcastss zmm12,DWORD PTR [r9+0]\t\n"
      "vmulps zmm0,zmm13,zmm12\t\n"
      "vmulps zmm1,zmm14,zmm12\t\n"
      "add r12, r13\t\n"
      "vbroadcastss zmm12,DWORD PTR [r9+4]\t\n"
      "vmulps zmm2,zmm13,zmm12\t\n"
      "vmulps zmm3,zmm14,zmm12\t\n"
      "add r12, r13\t\n"
      "vbroadcastss zmm12,DWORD PTR [r9+8]\t\n"
      "vmulps zmm4,zmm13,zmm12\t\n"
      "vmulps zmm5,zmm14,zmm12\t\n"
      "add r12, r13\t\n"
      "vbroadcastss zmm12,DWORD PTR [r9+12]\t\n"
      "vmulps zmm6,zmm13,zmm12\t\n"
      "vmulps zmm7,zmm14,zmm12\t\n"
      "add r12, r13\t\n"
      "vbroadcastss zmm12,DWORD PTR [r9+16]\t\n"
      "vmulps zmm8,zmm13,zmm12\t\n"
      "vmulps zmm9,zmm14,zmm12\t\n"
      "add r12, r13\t\n"
      "vbroadcastss zmm12,DWORD PTR [r9+20]\t\n"
      "vmulps zmm10,zmm13,zmm12\t\n"
      "vmulps zmm11,zmm14,zmm12\t\n"
      "mov r12, rcx\t\n"
      "test r14,r14\t\n"
      "jnz next_inner%=\t\n"
      "add r10,64\t\n"
      "jmp dump_C%=\t\n"

      "loop_inner%=:\t\n"

      "vmovaps zmm13,zmm31\t\n"
      "vcvtph2ps zmm14,YMMWORD PTR [r10 + 32]\t\n"
      "vcvtph2ps zmm31,YMMWORD PTR [r10 + 64]\t\n"
      "vbroadcastss zmm12,DWORD PTR [r9+0]\t\n"
      "vfmadd231ps zmm0,zmm13,zmm12\t\n"
      "vfmadd231ps zmm1,zmm14,zmm12\t\n"
      "vbroadcastss zmm12,DWORD PTR [r9+4]\t\n"
      "vfmadd231ps zmm2,zmm13,zmm12\t\n"
      "vfmadd231ps zmm3,zmm14,zmm12\t\n"
      "vbroadcastss zmm12,DWORD PTR [r9+8]\t\n"
      "vfmadd231ps zmm4,zmm13,zmm12\t\n"
      "vfmadd231ps zmm5,zmm14,zmm12\t\n"
      "vbroadcastss zmm12,DWORD PTR [r9+12]\t\n"
      "vfmadd231ps zmm6,zmm13,zmm12\t\n"
      "vfmadd231ps zmm7,zmm14,zmm12\t\n"
      "vbroadcastss zmm12,DWORD PTR [r9+16]\t\n"
      "vfmadd231ps zmm8,zmm13,zmm12\t\n"
      "vfmadd231ps zmm9,zmm14,zmm12\t\n"
      "vbroadcastss zmm12,DWORD PTR [r9+20]\t\n"
      "vfmadd231ps zmm10,zmm13,zmm12\t\n"
      "vfmadd231ps zmm11,zmm14,zmm12\t\n"

      "next_inner%=:\t\n"
      "add r9,24\t\n"
      "add r10,64\t\n"
      "dec r14\t\n"
      "jnz loop_inner%=\t\n"

      "vmovaps zmm13,zmm31\t\n"
      "vcvtph2ps zmm14,YMMWORD PTR [r10 + 32]\t\n"
      "vbroadcastss zmm12,DWORD PTR [r9+0]\t\n"
      "vfmadd231ps zmm0,zmm13,zmm12\t\n"
      "vfmadd231ps zmm1,zmm14,zmm12\t\n"
      "vbroadcastss zmm12,DWORD PTR [r9+4]\t\n"
      "vfmadd231ps zmm2,zmm13,zmm12\t\n"
      "vfmadd231ps zmm3,zmm14,zmm12\t\n"
      "vbroadcastss zmm12,DWORD PTR [r9+8]\t\n"
      "vfmadd231ps zmm4,zmm13,zmm12\t\n"
      "vfmadd231ps zmm5,zmm14,zmm12\t\n"
      "vbroadcastss zmm12,DWORD PTR [r9+12]\t\n"
      "vfmadd231ps zmm6,zmm13,zmm12\t\n"
      "vfmadd231ps zmm7,zmm14,zmm12\t\n"
      "vbroadcastss zmm12,DWORD PTR [r9+16]\t\n"
      "vfmadd231ps zmm8,zmm13,zmm12\t\n"
      "vfmadd231ps zmm9,zmm14,zmm12\t\n"
      "vbroadcastss zmm12,DWORD PTR [r9+20]\t\n"
      "vfmadd231ps zmm10,zmm13,zmm12\t\n"
      "vfmadd231ps zmm11,zmm14,zmm12\t\n"
      "add r9,24\t\n"
      "add r10,64\t\n"
      // Dump C
      "dump_C%=:\t\n"
      "vmovups zmmword PTR [r12 + 0], zmm0\t\n"
      "vmovups zmmword PTR [r12 + 64], zmm1\t\n"
      "add r12, r13\t\n"
      "vmovups zmmword PTR [r12 + 0], zmm2\t\n"
      "vmovups zmmword PTR [r12 + 64], zmm3\t\n"
      "add r12, r13\t\n"
      "vmovups zmmword PTR [r12 + 0], zmm4\t\n"
      "vmovups zmmword PTR [r12 + 64], zmm5\t\n"
      "add r12, r13\t\n"
      "vmovups zmmword PTR [r12 + 0], zmm6\t\n"
      "vmovups zmmword PTR [r12 + 64], zmm7\t\n"
      "add r12, r13\t\n"
      "vmovups zmmword PTR [r12 + 0], zmm8\t\n"
      "vmovups zmmword PTR [r12 + 64], zmm9\t\n"
      "add r12, r13\t\n"
      "vmovups zmmword PTR [r12 + 0], zmm10\t\n"
      "vmovups zmmword PTR [r12 + 64], zmm11\t\n"

      // next outer iteration
      "add rcx, 128\t\n"
      "mov r12, rcx\t\n"
      "mov r9, rax\t\n"
      "inc rbx\t\n"
      "cmp rbx, rdi\t\n"
      "jl loop_outter%=\t\n"
      :
      : [gp] "rm"(gp)
      : "r8",
        "r9",
        "r10",
        "r11",
        "r13",
        "r14",
        "rax",
        "rcx",
        "rsi",
        "rdi",
        "rbx",
        "r12",
        "r15",
        "memory");
}
void NOINLINE
gemmkernel_7x2_Avx512_fp16_fA0fB0fC0(GemmParamsFP16* gp) {
  asm volatile(
#if !defined(__clang__)
      "mov r14, %[gp]\t\n"
#else
      "mov %[gp], %%r14\t\n"
      ".intel_syntax noprefix\t\n"
#endif

      // Copy parameters
      // k
      "mov r8, [r14 + 0]\t\n"
      "dec r8\t\n"
      // A
      "mov r9, [r14 + 8]\t\n"
      // B
      "mov r10, [r14 + 16]\t\n"
      // beta
      "lea r15, [r14 + 24]\t\n"
      // C
      "mov r12, [r14 + 32]\t\n"
      // ldc
      "mov r13, [r14 + 40]\t\n"
      // b_block_cols
      "mov rdi, [r14 + 48]\t\n"
      // b_block_size
      "mov rsi, [r14 + 56]\t\n"

      // Make copies of A and C
      "mov rax, r9\t\n"
      "mov rcx, r12\t\n"

      "mov rbx, 0\t\n"
      "loop_outter%=:\t\n"
      "mov r14, r8\t\n"
      "vbroadcastss zmm31,DWORD PTR [r15]\t\n"
      "vcvtph2ps zmm15,YMMWORD PTR [r10 + 0]\t\n"
      "vcvtph2ps zmm16,YMMWORD PTR [r10 + 32]\t\n"
      "vxorps xmm0, xmm0, xmm0\t\n"
      "vcomiss xmm31, xmm0\t\n"
      "jz zero_regs%=\t\n"

      // Setup values with beta multiplication
      "vmulps zmm0, zmm31, [r12 + 0]\t\n"
      "vmulps zmm1, zmm31, [r12 + 64]\t\n"
      "add r12, r13\t\n"
      "vmulps zmm2, zmm31, [r12 + 0]\t\n"
      "vmulps zmm3, zmm31, [r12 + 64]\t\n"
      "add r12, r13\t\n"
      "vmulps zmm4, zmm31, [r12 + 0]\t\n"
      "vmulps zmm5, zmm31, [r12 + 64]\t\n"
      "add r12, r13\t\n"
      "vmulps zmm6, zmm31, [r12 + 0]\t\n"
      "vmulps zmm7, zmm31, [r12 + 64]\t\n"
      "add r12, r13\t\n"
      "vmulps zmm8, zmm31, [r12 + 0]\t\n"
      "vmulps zmm9, zmm31, [r12 + 64]\t\n"
      "add r12, r13\t\n"
      "vmulps zmm10, zmm31, [r12 + 0]\t\n"
      "vmulps zmm11, zmm31, [r12 + 64]\t\n"
      "add r12, r13\t\n"
      "vmulps zmm12, zmm31, [r12 + 0]\t\n"
      "vmulps zmm13, zmm31, [r12 + 64]\t\n"
      "test r14,r14\t\n"
      "jz skip_preload%=\t\n"
      "vcvtph2ps zmm31,YMMWORD PTR [r10 + 64]\t\n"
      "skip_preload%=:\t\n"
      "vbroadcastss zmm14,DWORD PTR [r9+0]\t\n"
      "vfmadd231ps zmm0,zmm15,zmm14\t\n"
      "vfmadd231ps zmm1,zmm16,zmm14\t\n"
      "vbroadcastss zmm14,DWORD PTR [r9+4]\t\n"
      "vfmadd231ps zmm2,zmm15,zmm14\t\n"
      "vfmadd231ps zmm3,zmm16,zmm14\t\n"
      "vbroadcastss zmm14,DWORD PTR [r9+8]\t\n"
      "vfmadd231ps zmm4,zmm15,zmm14\t\n"
      "vfmadd231ps zmm5,zmm16,zmm14\t\n"
      "vbroadcastss zmm14,DWORD PTR [r9+12]\t\n"
      "vfmadd231ps zmm6,zmm15,zmm14\t\n"
      "vfmadd231ps zmm7,zmm16,zmm14\t\n"
      "vbroadcastss zmm14,DWORD PTR [r9+16]\t\n"
      "vfmadd231ps zmm8,zmm15,zmm14\t\n"
      "vfmadd231ps zmm9,zmm16,zmm14\t\n"
      "vbroadcastss zmm14,DWORD PTR [r9+20]\t\n"
      "vfmadd231ps zmm10,zmm15,zmm14\t\n"
      "vfmadd231ps zmm11,zmm16,zmm14\t\n"
      "vbroadcastss zmm14,DWORD PTR [r9+24]\t\n"
      "vfmadd231ps zmm12,zmm15,zmm14\t\n"
      "vfmadd231ps zmm13,zmm16,zmm14\t\n"
      "mov r12, rcx\t\n"
      "test r14,r14\t\n"
      "jnz next_inner%=\t\n"
      "add r10,64\t\n"
      "jmp dump_C%=\t\n"

      "zero_regs%=:\t\n"

      "test r14,r14\t\n"
      "jz skip_preload_b_zero%=\t\n"
      "vcvtph2ps zmm31,YMMWORD PTR [r10 + 64]\t\n"
      "skip_preload_b_zero%=:\t\n"
      "vbroadcastss zmm14,DWORD PTR [r9+0]\t\n"
      "vmulps zmm0,zmm15,zmm14\t\n"
      "vmulps zmm1,zmm16,zmm14\t\n"
      "add r12, r13\t\n"
      "vbroadcastss zmm14,DWORD PTR [r9+4]\t\n"
      "vmulps zmm2,zmm15,zmm14\t\n"
      "vmulps zmm3,zmm16,zmm14\t\n"
      "add r12, r13\t\n"
      "vbroadcastss zmm14,DWORD PTR [r9+8]\t\n"
      "vmulps zmm4,zmm15,zmm14\t\n"
      "vmulps zmm5,zmm16,zmm14\t\n"
      "add r12, r13\t\n"
      "vbroadcastss zmm14,DWORD PTR [r9+12]\t\n"
      "vmulps zmm6,zmm15,zmm14\t\n"
      "vmulps zmm7,zmm16,zmm14\t\n"
      "add r12, r13\t\n"
      "vbroadcastss zmm14,DWORD PTR [r9+16]\t\n"
      "vmulps zmm8,zmm15,zmm14\t\n"
      "vmulps zmm9,zmm16,zmm14\t\n"
      "add r12, r13\t\n"
      "vbroadcastss zmm14,DWORD PTR [r9+20]\t\n"
      "vmulps zmm10,zmm15,zmm14\t\n"
      "vmulps zmm11,zmm16,zmm14\t\n"
      "add r12, r13\t\n"
      "vbroadcastss zmm14,DWORD PTR [r9+24]\t\n"
      "vmulps zmm12,zmm15,zmm14\t\n"
      "vmulps zmm13,zmm16,zmm14\t\n"
      "mov r12, rcx\t\n"
      "test r14,r14\t\n"
      "jnz next_inner%=\t\n"
      "add r10,64\t\n"
      "jmp dump_C%=\t\n"

      "loop_inner%=:\t\n"

      "vmovaps zmm15,zmm31\t\n"
      "vcvtph2ps zmm16,YMMWORD PTR [r10 + 32]\t\n"
      "vcvtph2ps zmm31,YMMWORD PTR [r10 + 64]\t\n"
      "vbroadcastss zmm14,DWORD PTR [r9+0]\t\n"
      "vfmadd231ps zmm0,zmm15,zmm14\t\n"
      "vfmadd231ps zmm1,zmm16,zmm14\t\n"
      "vbroadcastss zmm14,DWORD PTR [r9+4]\t\n"
      "vfmadd231ps zmm2,zmm15,zmm14\t\n"
      "vfmadd231ps zmm3,zmm16,zmm14\t\n"
      "vbroadcastss zmm14,DWORD PTR [r9+8]\t\n"
      "vfmadd231ps zmm4,zmm15,zmm14\t\n"
      "vfmadd231ps zmm5,zmm16,zmm14\t\n"
      "vbroadcastss zmm14,DWORD PTR [r9+12]\t\n"
      "vfmadd231ps zmm6,zmm15,zmm14\t\n"
      "vfmadd231ps zmm7,zmm16,zmm14\t\n"
      "vbroadcastss zmm14,DWORD PTR [r9+16]\t\n"
      "vfmadd231ps zmm8,zmm15,zmm14\t\n"
      "vfmadd231ps zmm9,zmm16,zmm14\t\n"
      "vbroadcastss zmm14,DWORD PTR [r9+20]\t\n"
      "vfmadd231ps zmm10,zmm15,zmm14\t\n"
      "vfmadd231ps zmm11,zmm16,zmm14\t\n"
      "vbroadcastss zmm14,DWORD PTR [r9+24]\t\n"
      "vfmadd231ps zmm12,zmm15,zmm14\t\n"
      "vfmadd231ps zmm13,zmm16,zmm14\t\n"

      "next_inner%=:\t\n"
      "add r9,28\t\n"
      "add r10,64\t\n"
      "dec r14\t\n"
      "jnz loop_inner%=\t\n"

      "vmovaps zmm15,zmm31\t\n"
      "vcvtph2ps zmm16,YMMWORD PTR [r10 + 32]\t\n"
      "vbroadcastss zmm14,DWORD PTR [r9+0]\t\n"
      "vfmadd231ps zmm0,zmm15,zmm14\t\n"
      "vfmadd231ps zmm1,zmm16,zmm14\t\n"
      "vbroadcastss zmm14,DWORD PTR [r9+4]\t\n"
      "vfmadd231ps zmm2,zmm15,zmm14\t\n"
      "vfmadd231ps zmm3,zmm16,zmm14\t\n"
      "vbroadcastss zmm14,DWORD PTR [r9+8]\t\n"
      "vfmadd231ps zmm4,zmm15,zmm14\t\n"
      "vfmadd231ps zmm5,zmm16,zmm14\t\n"
      "vbroadcastss zmm14,DWORD PTR [r9+12]\t\n"
      "vfmadd231ps zmm6,zmm15,zmm14\t\n"
      "vfmadd231ps zmm7,zmm16,zmm14\t\n"
      "vbroadcastss zmm14,DWORD PTR [r9+16]\t\n"
      "vfmadd231ps zmm8,zmm15,zmm14\t\n"
      "vfmadd231ps zmm9,zmm16,zmm14\t\n"
      "vbroadcastss zmm14,DWORD PTR [r9+20]\t\n"
      "vfmadd231ps zmm10,zmm15,zmm14\t\n"
      "vfmadd231ps zmm11,zmm16,zmm14\t\n"
      "vbroadcastss zmm14,DWORD PTR [r9+24]\t\n"
      "vfmadd231ps zmm12,zmm15,zmm14\t\n"
      "vfmadd231ps zmm13,zmm16,zmm14\t\n"
      "add r9,28\t\n"
      "add r10,64\t\n"
      // Dump C
      "dump_C%=:\t\n"
      "vmovups zmmword PTR [r12 + 0], zmm0\t\n"
      "vmovups zmmword PTR [r12 + 64], zmm1\t\n"
      "add r12, r13\t\n"
      "vmovups zmmword PTR [r12 + 0], zmm2\t\n"
      "vmovups zmmword PTR [r12 + 64], zmm3\t\n"
      "add r12, r13\t\n"
      "vmovups zmmword PTR [r12 + 0], zmm4\t\n"
      "vmovups zmmword PTR [r12 + 64], zmm5\t\n"
      "add r12, r13\t\n"
      "vmovups zmmword PTR [r12 + 0], zmm6\t\n"
      "vmovups zmmword PTR [r12 + 64], zmm7\t\n"
      "add r12, r13\t\n"
      "vmovups zmmword PTR [r12 + 0], zmm8\t\n"
      "vmovups zmmword PTR [r12 + 64], zmm9\t\n"
      "add r12, r13\t\n"
      "vmovups zmmword PTR [r12 + 0], zmm10\t\n"
      "vmovups zmmword PTR [r12 + 64], zmm11\t\n"
      "add r12, r13\t\n"
      "vmovups zmmword PTR [r12 + 0], zmm12\t\n"
      "vmovups zmmword PTR [r12 + 64], zmm13\t\n"

      // next outer iteration
      "add rcx, 128\t\n"
      "mov r12, rcx\t\n"
      "mov r9, rax\t\n"
      "inc rbx\t\n"
      "cmp rbx, rdi\t\n"
      "jl loop_outter%=\t\n"
      :
      : [gp] "rm"(gp)
      : "r8",
        "r9",
        "r10",
        "r11",
        "r13",
        "r14",
        "rax",
        "rcx",
        "rsi",
        "rdi",
        "rbx",
        "r12",
        "r15",
        "memory");
}
void NOINLINE
gemmkernel_8x2_Avx512_fp16_fA0fB0fC0(GemmParamsFP16* gp) {
  asm volatile(
#if !defined(__clang__)
      "mov r14, %[gp]\t\n"
#else
      "mov %[gp], %%r14\t\n"
      ".intel_syntax noprefix\t\n"
#endif

      // Copy parameters
      // k
      "mov r8, [r14 + 0]\t\n"
      "dec r8\t\n"
      // A
      "mov r9, [r14 + 8]\t\n"
      // B
      "mov r10, [r14 + 16]\t\n"
      // beta
      "lea r15, [r14 + 24]\t\n"
      // C
      "mov r12, [r14 + 32]\t\n"
      // ldc
      "mov r13, [r14 + 40]\t\n"
      // b_block_cols
      "mov rdi, [r14 + 48]\t\n"
      // b_block_size
      "mov rsi, [r14 + 56]\t\n"

      // Make copies of A and C
      "mov rax, r9\t\n"
      "mov rcx, r12\t\n"

      "mov rbx, 0\t\n"
      "loop_outter%=:\t\n"
      "mov r14, r8\t\n"
      "vbroadcastss zmm31,DWORD PTR [r15]\t\n"
      "vcvtph2ps zmm17,YMMWORD PTR [r10 + 0]\t\n"
      "vcvtph2ps zmm18,YMMWORD PTR [r10 + 32]\t\n"
      "vxorps xmm0, xmm0, xmm0\t\n"
      "vcomiss xmm31, xmm0\t\n"
      "jz zero_regs%=\t\n"

      // Setup values with beta multiplication
      "vmulps zmm0, zmm31, [r12 + 0]\t\n"
      "vmulps zmm1, zmm31, [r12 + 64]\t\n"
      "add r12, r13\t\n"
      "vmulps zmm2, zmm31, [r12 + 0]\t\n"
      "vmulps zmm3, zmm31, [r12 + 64]\t\n"
      "add r12, r13\t\n"
      "vmulps zmm4, zmm31, [r12 + 0]\t\n"
      "vmulps zmm5, zmm31, [r12 + 64]\t\n"
      "add r12, r13\t\n"
      "vmulps zmm6, zmm31, [r12 + 0]\t\n"
      "vmulps zmm7, zmm31, [r12 + 64]\t\n"
      "add r12, r13\t\n"
      "vmulps zmm8, zmm31, [r12 + 0]\t\n"
      "vmulps zmm9, zmm31, [r12 + 64]\t\n"
      "add r12, r13\t\n"
      "vmulps zmm10, zmm31, [r12 + 0]\t\n"
      "vmulps zmm11, zmm31, [r12 + 64]\t\n"
      "add r12, r13\t\n"
      "vmulps zmm12, zmm31, [r12 + 0]\t\n"
      "vmulps zmm13, zmm31, [r12 + 64]\t\n"
      "add r12, r13\t\n"
      "vmulps zmm14, zmm31, [r12 + 0]\t\n"
      "vmulps zmm15, zmm31, [r12 + 64]\t\n"
      "test r14,r14\t\n"
      "jz skip_preload%=\t\n"
      "vcvtph2ps zmm31,YMMWORD PTR [r10 + 64]\t\n"
      "skip_preload%=:\t\n"
      "vbroadcastss zmm16,DWORD PTR [r9+0]\t\n"
      "vfmadd231ps zmm0,zmm17,zmm16\t\n"
      "vfmadd231ps zmm1,zmm18,zmm16\t\n"
      "vbroadcastss zmm16,DWORD PTR [r9+4]\t\n"
      "vfmadd231ps zmm2,zmm17,zmm16\t\n"
      "vfmadd231ps zmm3,zmm18,zmm16\t\n"
      "vbroadcastss zmm16,DWORD PTR [r9+8]\t\n"
      "vfmadd231ps zmm4,zmm17,zmm16\t\n"
      "vfmadd231ps zmm5,zmm18,zmm16\t\n"
      "vbroadcastss zmm16,DWORD PTR [r9+12]\t\n"
      "vfmadd231ps zmm6,zmm17,zmm16\t\n"
      "vfmadd231ps zmm7,zmm18,zmm16\t\n"
      "vbroadcastss zmm16,DWORD PTR [r9+16]\t\n"
      "vfmadd231ps zmm8,zmm17,zmm16\t\n"
      "vfmadd231ps zmm9,zmm18,zmm16\t\n"
      "vbroadcastss zmm16,DWORD PTR [r9+20]\t\n"
      "vfmadd231ps zmm10,zmm17,zmm16\t\n"
      "vfmadd231ps zmm11,zmm18,zmm16\t\n"
      "vbroadcastss zmm16,DWORD PTR [r9+24]\t\n"
      "vfmadd231ps zmm12,zmm17,zmm16\t\n"
      "vfmadd231ps zmm13,zmm18,zmm16\t\n"
      "vbroadcastss zmm16,DWORD PTR [r9+28]\t\n"
      "vfmadd231ps zmm14,zmm17,zmm16\t\n"
      "vfmadd231ps zmm15,zmm18,zmm16\t\n"
      "mov r12, rcx\t\n"
      "test r14,r14\t\n"
      "jnz next_inner%=\t\n"
      "add r10,64\t\n"
      "jmp dump_C%=\t\n"

      "zero_regs%=:\t\n"

      "test r14,r14\t\n"
      "jz skip_preload_b_zero%=\t\n"
      "vcvtph2ps zmm31,YMMWORD PTR [r10 + 64]\t\n"
      "skip_preload_b_zero%=:\t\n"
      "vbroadcastss zmm16,DWORD PTR [r9+0]\t\n"
      "vmulps zmm0,zmm17,zmm16\t\n"
      "vmulps zmm1,zmm18,zmm16\t\n"
      "add r12, r13\t\n"
      "vbroadcastss zmm16,DWORD PTR [r9+4]\t\n"
      "vmulps zmm2,zmm17,zmm16\t\n"
      "vmulps zmm3,zmm18,zmm16\t\n"
      "add r12, r13\t\n"
      "vbroadcastss zmm16,DWORD PTR [r9+8]\t\n"
      "vmulps zmm4,zmm17,zmm16\t\n"
      "vmulps zmm5,zmm18,zmm16\t\n"
      "add r12, r13\t\n"
      "vbroadcastss zmm16,DWORD PTR [r9+12]\t\n"
      "vmulps zmm6,zmm17,zmm16\t\n"
      "vmulps zmm7,zmm18,zmm16\t\n"
      "add r12, r13\t\n"
      "vbroadcastss zmm16,DWORD PTR [r9+16]\t\n"
      "vmulps zmm8,zmm17,zmm16\t\n"
      "vmulps zmm9,zmm18,zmm16\t\n"
      "add r12, r13\t\n"
      "vbroadcastss zmm16,DWORD PTR [r9+20]\t\n"
      "vmulps zmm10,zmm17,zmm16\t\n"
      "vmulps zmm11,zmm18,zmm16\t\n"
      "add r12, r13\t\n"
      "vbroadcastss zmm16,DWORD PTR [r9+24]\t\n"
      "vmulps zmm12,zmm17,zmm16\t\n"
      "vmulps zmm13,zmm18,zmm16\t\n"
      "add r12, r13\t\n"
      "vbroadcastss zmm16,DWORD PTR [r9+28]\t\n"
      "vmulps zmm14,zmm17,zmm16\t\n"
      "vmulps zmm15,zmm18,zmm16\t\n"
      "mov r12, rcx\t\n"
      "test r14,r14\t\n"
      "jnz next_inner%=\t\n"
      "add r10,64\t\n"
      "jmp dump_C%=\t\n"

      "loop_inner%=:\t\n"

      "vmovaps zmm17,zmm31\t\n"
      "vcvtph2ps zmm18,YMMWORD PTR [r10 + 32]\t\n"
      "vcvtph2ps zmm31,YMMWORD PTR [r10 + 64]\t\n"
      "vbroadcastss zmm16,DWORD PTR [r9+0]\t\n"
      "vfmadd231ps zmm0,zmm17,zmm16\t\n"
      "vfmadd231ps zmm1,zmm18,zmm16\t\n"
      "vbroadcastss zmm16,DWORD PTR [r9+4]\t\n"
      "vfmadd231ps zmm2,zmm17,zmm16\t\n"
      "vfmadd231ps zmm3,zmm18,zmm16\t\n"
      "vbroadcastss zmm16,DWORD PTR [r9+8]\t\n"
      "vfmadd231ps zmm4,zmm17,zmm16\t\n"
      "vfmadd231ps zmm5,zmm18,zmm16\t\n"
      "vbroadcastss zmm16,DWORD PTR [r9+12]\t\n"
      "vfmadd231ps zmm6,zmm17,zmm16\t\n"
      "vfmadd231ps zmm7,zmm18,zmm16\t\n"
      "vbroadcastss zmm16,DWORD PTR [r9+16]\t\n"
      "vfmadd231ps zmm8,zmm17,zmm16\t\n"
      "vfmadd231ps zmm9,zmm18,zmm16\t\n"
      "vbroadcastss zmm16,DWORD PTR [r9+20]\t\n"
      "vfmadd231ps zmm10,zmm17,zmm16\t\n"
      "vfmadd231ps zmm11,zmm18,zmm16\t\n"
      "vbroadcastss zmm16,DWORD PTR [r9+24]\t\n"
      "vfmadd231ps zmm12,zmm17,zmm16\t\n"
      "vfmadd231ps zmm13,zmm18,zmm16\t\n"
      "vbroadcastss zmm16,DWORD PTR [r9+28]\t\n"
      "vfmadd231ps zmm14,zmm17,zmm16\t\n"
      "vfmadd231ps zmm15,zmm18,zmm16\t\n"

      "next_inner%=:\t\n"
      "add r9,32\t\n"
      "add r10,64\t\n"
      "dec r14\t\n"
      "jnz loop_inner%=\t\n"

      "vmovaps zmm17,zmm31\t\n"
      "vcvtph2ps zmm18,YMMWORD PTR [r10 + 32]\t\n"
      "vbroadcastss zmm16,DWORD PTR [r9+0]\t\n"
      "vfmadd231ps zmm0,zmm17,zmm16\t\n"
      "vfmadd231ps zmm1,zmm18,zmm16\t\n"
      "vbroadcastss zmm16,DWORD PTR [r9+4]\t\n"
      "vfmadd231ps zmm2,zmm17,zmm16\t\n"
      "vfmadd231ps zmm3,zmm18,zmm16\t\n"
      "vbroadcastss zmm16,DWORD PTR [r9+8]\t\n"
      "vfmadd231ps zmm4,zmm17,zmm16\t\n"
      "vfmadd231ps zmm5,zmm18,zmm16\t\n"
      "vbroadcastss zmm16,DWORD PTR [r9+12]\t\n"
      "vfmadd231ps zmm6,zmm17,zmm16\t\n"
      "vfmadd231ps zmm7,zmm18,zmm16\t\n"
      "vbroadcastss zmm16,DWORD PTR [r9+16]\t\n"
      "vfmadd231ps zmm8,zmm17,zmm16\t\n"
      "vfmadd231ps zmm9,zmm18,zmm16\t\n"
      "vbroadcastss zmm16,DWORD PTR [r9+20]\t\n"
      "vfmadd231ps zmm10,zmm17,zmm16\t\n"
      "vfmadd231ps zmm11,zmm18,zmm16\t\n"
      "vbroadcastss zmm16,DWORD PTR [r9+24]\t\n"
      "vfmadd231ps zmm12,zmm17,zmm16\t\n"
      "vfmadd231ps zmm13,zmm18,zmm16\t\n"
      "vbroadcastss zmm16,DWORD PTR [r9+28]\t\n"
      "vfmadd231ps zmm14,zmm17,zmm16\t\n"
      "vfmadd231ps zmm15,zmm18,zmm16\t\n"
      "add r9,32\t\n"
      "add r10,64\t\n"
      // Dump C
      "dump_C%=:\t\n"
      "vmovups zmmword PTR [r12 + 0], zmm0\t\n"
      "vmovups zmmword PTR [r12 + 64], zmm1\t\n"
      "add r12, r13\t\n"
      "vmovups zmmword PTR [r12 + 0], zmm2\t\n"
      "vmovups zmmword PTR [r12 + 64], zmm3\t\n"
      "add r12, r13\t\n"
      "vmovups zmmword PTR [r12 + 0], zmm4\t\n"
      "vmovups zmmword PTR [r12 + 64], zmm5\t\n"
      "add r12, r13\t\n"
      "vmovups zmmword PTR [r12 + 0], zmm6\t\n"
      "vmovups zmmword PTR [r12 + 64], zmm7\t\n"
      "add r12, r13\t\n"
      "vmovups zmmword PTR [r12 + 0], zmm8\t\n"
      "vmovups zmmword PTR [r12 + 64], zmm9\t\n"
      "add r12, r13\t\n"
      "vmovups zmmword PTR [r12 + 0], zmm10\t\n"
      "vmovups zmmword PTR [r12 + 64], zmm11\t\n"
      "add r12, r13\t\n"
      "vmovups zmmword PTR [r12 + 0], zmm12\t\n"
      "vmovups zmmword PTR [r12 + 64], zmm13\t\n"
      "add r12, r13\t\n"
      "vmovups zmmword PTR [r12 + 0], zmm14\t\n"
      "vmovups zmmword PTR [r12 + 64], zmm15\t\n"

      // next outer iteration
      "add rcx, 128\t\n"
      "mov r12, rcx\t\n"
      "mov r9, rax\t\n"
      "inc rbx\t\n"
      "cmp rbx, rdi\t\n"
      "jl loop_outter%=\t\n"
      :
      : [gp] "rm"(gp)
      : "r8",
        "r9",
        "r10",
        "r11",
        "r13",
        "r14",
        "rax",
        "rcx",
        "rsi",
        "rdi",
        "rbx",
        "r12",
        "r15",
        "memory");
}
void NOINLINE
gemmkernel_9x2_Avx512_fp16_fA0fB0fC0(GemmParamsFP16* gp) {
  asm volatile(
#if !defined(__clang__)
      "mov r14, %[gp]\t\n"
#else
      "mov %[gp], %%r14\t\n"
      ".intel_syntax noprefix\t\n"
#endif

      // Copy parameters
      // k
      "mov r8, [r14 + 0]\t\n"
      "dec r8\t\n"
      // A
      "mov r9, [r14 + 8]\t\n"
      // B
      "mov r10, [r14 + 16]\t\n"
      // beta
      "lea r15, [r14 + 24]\t\n"
      // C
      "mov r12, [r14 + 32]\t\n"
      // ldc
      "mov r13, [r14 + 40]\t\n"
      // b_block_cols
      "mov rdi, [r14 + 48]\t\n"
      // b_block_size
      "mov rsi, [r14 + 56]\t\n"

      // Make copies of A and C
      "mov rax, r9\t\n"
      "mov rcx, r12\t\n"

      "mov rbx, 0\t\n"
      "loop_outter%=:\t\n"
      "mov r14, r8\t\n"
      "vbroadcastss zmm31,DWORD PTR [r15]\t\n"
      "vcvtph2ps zmm19,YMMWORD PTR [r10 + 0]\t\n"
      "vcvtph2ps zmm20,YMMWORD PTR [r10 + 32]\t\n"
      "vxorps xmm0, xmm0, xmm0\t\n"
      "vcomiss xmm31, xmm0\t\n"
      "jz zero_regs%=\t\n"

      // Setup values with beta multiplication
      "vmulps zmm0, zmm31, [r12 + 0]\t\n"
      "vmulps zmm1, zmm31, [r12 + 64]\t\n"
      "add r12, r13\t\n"
      "vmulps zmm2, zmm31, [r12 + 0]\t\n"
      "vmulps zmm3, zmm31, [r12 + 64]\t\n"
      "add r12, r13\t\n"
      "vmulps zmm4, zmm31, [r12 + 0]\t\n"
      "vmulps zmm5, zmm31, [r12 + 64]\t\n"
      "add r12, r13\t\n"
      "vmulps zmm6, zmm31, [r12 + 0]\t\n"
      "vmulps zmm7, zmm31, [r12 + 64]\t\n"
      "add r12, r13\t\n"
      "vmulps zmm8, zmm31, [r12 + 0]\t\n"
      "vmulps zmm9, zmm31, [r12 + 64]\t\n"
      "add r12, r13\t\n"
      "vmulps zmm10, zmm31, [r12 + 0]\t\n"
      "vmulps zmm11, zmm31, [r12 + 64]\t\n"
      "add r12, r13\t\n"
      "vmulps zmm12, zmm31, [r12 + 0]\t\n"
      "vmulps zmm13, zmm31, [r12 + 64]\t\n"
      "add r12, r13\t\n"
      "vmulps zmm14, zmm31, [r12 + 0]\t\n"
      "vmulps zmm15, zmm31, [r12 + 64]\t\n"
      "add r12, r13\t\n"
      "vmulps zmm16, zmm31, [r12 + 0]\t\n"
      "vmulps zmm17, zmm31, [r12 + 64]\t\n"
      "test r14,r14\t\n"
      "jz skip_preload%=\t\n"
      "vcvtph2ps zmm31,YMMWORD PTR [r10 + 64]\t\n"
      "skip_preload%=:\t\n"
      "vbroadcastss zmm18,DWORD PTR [r9+0]\t\n"
      "vfmadd231ps zmm0,zmm19,zmm18\t\n"
      "vfmadd231ps zmm1,zmm20,zmm18\t\n"
      "vbroadcastss zmm18,DWORD PTR [r9+4]\t\n"
      "vfmadd231ps zmm2,zmm19,zmm18\t\n"
      "vfmadd231ps zmm3,zmm20,zmm18\t\n"
      "vbroadcastss zmm18,DWORD PTR [r9+8]\t\n"
      "vfmadd231ps zmm4,zmm19,zmm18\t\n"
      "vfmadd231ps zmm5,zmm20,zmm18\t\n"
      "vbroadcastss zmm18,DWORD PTR [r9+12]\t\n"
      "vfmadd231ps zmm6,zmm19,zmm18\t\n"
      "vfmadd231ps zmm7,zmm20,zmm18\t\n"
      "vbroadcastss zmm18,DWORD PTR [r9+16]\t\n"
      "vfmadd231ps zmm8,zmm19,zmm18\t\n"
      "vfmadd231ps zmm9,zmm20,zmm18\t\n"
      "vbroadcastss zmm18,DWORD PTR [r9+20]\t\n"
      "vfmadd231ps zmm10,zmm19,zmm18\t\n"
      "vfmadd231ps zmm11,zmm20,zmm18\t\n"
      "vbroadcastss zmm18,DWORD PTR [r9+24]\t\n"
      "vfmadd231ps zmm12,zmm19,zmm18\t\n"
      "vfmadd231ps zmm13,zmm20,zmm18\t\n"
      "vbroadcastss zmm18,DWORD PTR [r9+28]\t\n"
      "vfmadd231ps zmm14,zmm19,zmm18\t\n"
      "vfmadd231ps zmm15,zmm20,zmm18\t\n"
      "vbroadcastss zmm18,DWORD PTR [r9+32]\t\n"
      "vfmadd231ps zmm16,zmm19,zmm18\t\n"
      "vfmadd231ps zmm17,zmm20,zmm18\t\n"
      "mov r12, rcx\t\n"
      "test r14,r14\t\n"
      "jnz next_inner%=\t\n"
      "add r10,64\t\n"
      "jmp dump_C%=\t\n"

      "zero_regs%=:\t\n"

      "test r14,r14\t\n"
      "jz skip_preload_b_zero%=\t\n"
      "vcvtph2ps zmm31,YMMWORD PTR [r10 + 64]\t\n"
      "skip_preload_b_zero%=:\t\n"
      "vbroadcastss zmm18,DWORD PTR [r9+0]\t\n"
      "vmulps zmm0,zmm19,zmm18\t\n"
      "vmulps zmm1,zmm20,zmm18\t\n"
      "add r12, r13\t\n"
      "vbroadcastss zmm18,DWORD PTR [r9+4]\t\n"
      "vmulps zmm2,zmm19,zmm18\t\n"
      "vmulps zmm3,zmm20,zmm18\t\n"
      "add r12, r13\t\n"
      "vbroadcastss zmm18,DWORD PTR [r9+8]\t\n"
      "vmulps zmm4,zmm19,zmm18\t\n"
      "vmulps zmm5,zmm20,zmm18\t\n"
      "add r12, r13\t\n"
      "vbroadcastss zmm18,DWORD PTR [r9+12]\t\n"
      "vmulps zmm6,zmm19,zmm18\t\n"
      "vmulps zmm7,zmm20,zmm18\t\n"
      "add r12, r13\t\n"
      "vbroadcastss zmm18,DWORD PTR [r9+16]\t\n"
      "vmulps zmm8,zmm19,zmm18\t\n"
      "vmulps zmm9,zmm20,zmm18\t\n"
      "add r12, r13\t\n"
      "vbroadcastss zmm18,DWORD PTR [r9+20]\t\n"
      "vmulps zmm10,zmm19,zmm18\t\n"
      "vmulps zmm11,zmm20,zmm18\t\n"
      "add r12, r13\t\n"
      "vbroadcastss zmm18,DWORD PTR [r9+24]\t\n"
      "vmulps zmm12,zmm19,zmm18\t\n"
      "vmulps zmm13,zmm20,zmm18\t\n"
      "add r12, r13\t\n"
      "vbroadcastss zmm18,DWORD PTR [r9+28]\t\n"
      "vmulps zmm14,zmm19,zmm18\t\n"
      "vmulps zmm15,zmm20,zmm18\t\n"
      "add r12, r13\t\n"
      "vbroadcastss zmm18,DWORD PTR [r9+32]\t\n"
      "vmulps zmm16,zmm19,zmm18\t\n"
      "vmulps zmm17,zmm20,zmm18\t\n"
      "mov r12, rcx\t\n"
      "test r14,r14\t\n"
      "jnz next_inner%=\t\n"
      "add r10,64\t\n"
      "jmp dump_C%=\t\n"

      "loop_inner%=:\t\n"

      "vmovaps zmm19,zmm31\t\n"
      "vcvtph2ps zmm20,YMMWORD PTR [r10 + 32]\t\n"
      "vcvtph2ps zmm31,YMMWORD PTR [r10 + 64]\t\n"
      "vbroadcastss zmm18,DWORD PTR [r9+0]\t\n"
      "vfmadd231ps zmm0,zmm19,zmm18\t\n"
      "vfmadd231ps zmm1,zmm20,zmm18\t\n"
      "vbroadcastss zmm18,DWORD PTR [r9+4]\t\n"
      "vfmadd231ps zmm2,zmm19,zmm18\t\n"
      "vfmadd231ps zmm3,zmm20,zmm18\t\n"
      "vbroadcastss zmm18,DWORD PTR [r9+8]\t\n"
      "vfmadd231ps zmm4,zmm19,zmm18\t\n"
      "vfmadd231ps zmm5,zmm20,zmm18\t\n"
      "vbroadcastss zmm18,DWORD PTR [r9+12]\t\n"
      "vfmadd231ps zmm6,zmm19,zmm18\t\n"
      "vfmadd231ps zmm7,zmm20,zmm18\t\n"
      "vbroadcastss zmm18,DWORD PTR [r9+16]\t\n"
      "vfmadd231ps zmm8,zmm19,zmm18\t\n"
      "vfmadd231ps zmm9,zmm20,zmm18\t\n"
      "vbroadcastss zmm18,DWORD PTR [r9+20]\t\n"
      "vfmadd231ps zmm10,zmm19,zmm18\t\n"
      "vfmadd231ps zmm11,zmm20,zmm18\t\n"
      "vbroadcastss zmm18,DWORD PTR [r9+24]\t\n"
      "vfmadd231ps zmm12,zmm19,zmm18\t\n"
      "vfmadd231ps zmm13,zmm20,zmm18\t\n"
      "vbroadcastss zmm18,DWORD PTR [r9+28]\t\n"
      "vfmadd231ps zmm14,zmm19,zmm18\t\n"
      "vfmadd231ps zmm15,zmm20,zmm18\t\n"
      "vbroadcastss zmm18,DWORD PTR [r9+32]\t\n"
      "vfmadd231ps zmm16,zmm19,zmm18\t\n"
      "vfmadd231ps zmm17,zmm20,zmm18\t\n"

      "next_inner%=:\t\n"
      "add r9,36\t\n"
      "add r10,64\t\n"
      "dec r14\t\n"
      "jnz loop_inner%=\t\n"

      "vmovaps zmm19,zmm31\t\n"
      "vcvtph2ps zmm20,YMMWORD PTR [r10 + 32]\t\n"
      "vbroadcastss zmm18,DWORD PTR [r9+0]\t\n"
      "vfmadd231ps zmm0,zmm19,zmm18\t\n"
      "vfmadd231ps zmm1,zmm20,zmm18\t\n"
      "vbroadcastss zmm18,DWORD PTR [r9+4]\t\n"
      "vfmadd231ps zmm2,zmm19,zmm18\t\n"
      "vfmadd231ps zmm3,zmm20,zmm18\t\n"
      "vbroadcastss zmm18,DWORD PTR [r9+8]\t\n"
      "vfmadd231ps zmm4,zmm19,zmm18\t\n"
      "vfmadd231ps zmm5,zmm20,zmm18\t\n"
      "vbroadcastss zmm18,DWORD PTR [r9+12]\t\n"
      "vfmadd231ps zmm6,zmm19,zmm18\t\n"
      "vfmadd231ps zmm7,zmm20,zmm18\t\n"
      "vbroadcastss zmm18,DWORD PTR [r9+16]\t\n"
      "vfmadd231ps zmm8,zmm19,zmm18\t\n"
      "vfmadd231ps zmm9,zmm20,zmm18\t\n"
      "vbroadcastss zmm18,DWORD PTR [r9+20]\t\n"
      "vfmadd231ps zmm10,zmm19,zmm18\t\n"
      "vfmadd231ps zmm11,zmm20,zmm18\t\n"
      "vbroadcastss zmm18,DWORD PTR [r9+24]\t\n"
      "vfmadd231ps zmm12,zmm19,zmm18\t\n"
      "vfmadd231ps zmm13,zmm20,zmm18\t\n"
      "vbroadcastss zmm18,DWORD PTR [r9+28]\t\n"
      "vfmadd231ps zmm14,zmm19,zmm18\t\n"
      "vfmadd231ps zmm15,zmm20,zmm18\t\n"
      "vbroadcastss zmm18,DWORD PTR [r9+32]\t\n"
      "vfmadd231ps zmm16,zmm19,zmm18\t\n"
      "vfmadd231ps zmm17,zmm20,zmm18\t\n"
      "add r9,36\t\n"
      "add r10,64\t\n"
      // Dump C
      "dump_C%=:\t\n"
      "vmovups zmmword PTR [r12 + 0], zmm0\t\n"
      "vmovups zmmword PTR [r12 + 64], zmm1\t\n"
      "add r12, r13\t\n"
      "vmovups zmmword PTR [r12 + 0], zmm2\t\n"
      "vmovups zmmword PTR [r12 + 64], zmm3\t\n"
      "add r12, r13\t\n"
      "vmovups zmmword PTR [r12 + 0], zmm4\t\n"
      "vmovups zmmword PTR [r12 + 64], zmm5\t\n"
      "add r12, r13\t\n"
      "vmovups zmmword PTR [r12 + 0], zmm6\t\n"
      "vmovups zmmword PTR [r12 + 64], zmm7\t\n"
      "add r12, r13\t\n"
      "vmovups zmmword PTR [r12 + 0], zmm8\t\n"
      "vmovups zmmword PTR [r12 + 64], zmm9\t\n"
      "add r12, r13\t\n"
      "vmovups zmmword PTR [r12 + 0], zmm10\t\n"
      "vmovups zmmword PTR [r12 + 64], zmm11\t\n"
      "add r12, r13\t\n"
      "vmovups zmmword PTR [r12 + 0], zmm12\t\n"
      "vmovups zmmword PTR [r12 + 64], zmm13\t\n"
      "add r12, r13\t\n"
      "vmovups zmmword PTR [r12 + 0], zmm14\t\n"
      "vmovups zmmword PTR [r12 + 64], zmm15\t\n"
      "add r12, r13\t\n"
      "vmovups zmmword PTR [r12 + 0], zmm16\t\n"
      "vmovups zmmword PTR [r12 + 64], zmm17\t\n"

      // next outer iteration
      "add rcx, 128\t\n"
      "mov r12, rcx\t\n"
      "mov r9, rax\t\n"
      "inc rbx\t\n"
      "cmp rbx, rdi\t\n"
      "jl loop_outter%=\t\n"
      :
      : [gp] "rm"(gp)
      : "r8",
        "r9",
        "r10",
        "r11",
        "r13",
        "r14",
        "rax",
        "rcx",
        "rsi",
        "rdi",
        "rbx",
        "r12",
        "r15",
        "memory");
}
void NOINLINE
gemmkernel_10x2_Avx512_fp16_fA0fB0fC0(GemmParamsFP16* gp) {
  asm volatile(
#if !defined(__clang__)
      "mov r14, %[gp]\t\n"
#else
      "mov %[gp], %%r14\t\n"
      ".intel_syntax noprefix\t\n"
#endif

      // Copy parameters
      // k
      "mov r8, [r14 + 0]\t\n"
      "dec r8\t\n"
      // A
      "mov r9, [r14 + 8]\t\n"
      // B
      "mov r10, [r14 + 16]\t\n"
      // beta
      "lea r15, [r14 + 24]\t\n"
      // C
      "mov r12, [r14 + 32]\t\n"
      // ldc
      "mov r13, [r14 + 40]\t\n"
      // b_block_cols
      "mov rdi, [r14 + 48]\t\n"
      // b_block_size
      "mov rsi, [r14 + 56]\t\n"

      // Make copies of A and C
      "mov rax, r9\t\n"
      "mov rcx, r12\t\n"

      "mov rbx, 0\t\n"
      "loop_outter%=:\t\n"
      "mov r14, r8\t\n"
      "vbroadcastss zmm31,DWORD PTR [r15]\t\n"
      "vcvtph2ps zmm21,YMMWORD PTR [r10 + 0]\t\n"
      "vcvtph2ps zmm22,YMMWORD PTR [r10 + 32]\t\n"
      "vxorps xmm0, xmm0, xmm0\t\n"
      "vcomiss xmm31, xmm0\t\n"
      "jz zero_regs%=\t\n"

      // Setup values with beta multiplication
      "vmulps zmm0, zmm31, [r12 + 0]\t\n"
      "vmulps zmm1, zmm31, [r12 + 64]\t\n"
      "add r12, r13\t\n"
      "vmulps zmm2, zmm31, [r12 + 0]\t\n"
      "vmulps zmm3, zmm31, [r12 + 64]\t\n"
      "add r12, r13\t\n"
      "vmulps zmm4, zmm31, [r12 + 0]\t\n"
      "vmulps zmm5, zmm31, [r12 + 64]\t\n"
      "add r12, r13\t\n"
      "vmulps zmm6, zmm31, [r12 + 0]\t\n"
      "vmulps zmm7, zmm31, [r12 + 64]\t\n"
      "add r12, r13\t\n"
      "vmulps zmm8, zmm31, [r12 + 0]\t\n"
      "vmulps zmm9, zmm31, [r12 + 64]\t\n"
      "add r12, r13\t\n"
      "vmulps zmm10, zmm31, [r12 + 0]\t\n"
      "vmulps zmm11, zmm31, [r12 + 64]\t\n"
      "add r12, r13\t\n"
      "vmulps zmm12, zmm31, [r12 + 0]\t\n"
      "vmulps zmm13, zmm31, [r12 + 64]\t\n"
      "add r12, r13\t\n"
      "vmulps zmm14, zmm31, [r12 + 0]\t\n"
      "vmulps zmm15, zmm31, [r12 + 64]\t\n"
      "add r12, r13\t\n"
      "vmulps zmm16, zmm31, [r12 + 0]\t\n"
      "vmulps zmm17, zmm31, [r12 + 64]\t\n"
      "add r12, r13\t\n"
      "vmulps zmm18, zmm31, [r12 + 0]\t\n"
      "vmulps zmm19, zmm31, [r12 + 64]\t\n"
      "test r14,r14\t\n"
      "jz skip_preload%=\t\n"
      "vcvtph2ps zmm31,YMMWORD PTR [r10 + 64]\t\n"
      "skip_preload%=:\t\n"
      "vbroadcastss zmm20,DWORD PTR [r9+0]\t\n"
      "vfmadd231ps zmm0,zmm21,zmm20\t\n"
      "vfmadd231ps zmm1,zmm22,zmm20\t\n"
      "vbroadcastss zmm20,DWORD PTR [r9+4]\t\n"
      "vfmadd231ps zmm2,zmm21,zmm20\t\n"
      "vfmadd231ps zmm3,zmm22,zmm20\t\n"
      "vbroadcastss zmm20,DWORD PTR [r9+8]\t\n"
      "vfmadd231ps zmm4,zmm21,zmm20\t\n"
      "vfmadd231ps zmm5,zmm22,zmm20\t\n"
      "vbroadcastss zmm20,DWORD PTR [r9+12]\t\n"
      "vfmadd231ps zmm6,zmm21,zmm20\t\n"
      "vfmadd231ps zmm7,zmm22,zmm20\t\n"
      "vbroadcastss zmm20,DWORD PTR [r9+16]\t\n"
      "vfmadd231ps zmm8,zmm21,zmm20\t\n"
      "vfmadd231ps zmm9,zmm22,zmm20\t\n"
      "vbroadcastss zmm20,DWORD PTR [r9+20]\t\n"
      "vfmadd231ps zmm10,zmm21,zmm20\t\n"
      "vfmadd231ps zmm11,zmm22,zmm20\t\n"
      "vbroadcastss zmm20,DWORD PTR [r9+24]\t\n"
      "vfmadd231ps zmm12,zmm21,zmm20\t\n"
      "vfmadd231ps zmm13,zmm22,zmm20\t\n"
      "vbroadcastss zmm20,DWORD PTR [r9+28]\t\n"
      "vfmadd231ps zmm14,zmm21,zmm20\t\n"
      "vfmadd231ps zmm15,zmm22,zmm20\t\n"
      "vbroadcastss zmm20,DWORD PTR [r9+32]\t\n"
      "vfmadd231ps zmm16,zmm21,zmm20\t\n"
      "vfmadd231ps zmm17,zmm22,zmm20\t\n"
      "vbroadcastss zmm20,DWORD PTR [r9+36]\t\n"
      "vfmadd231ps zmm18,zmm21,zmm20\t\n"
      "vfmadd231ps zmm19,zmm22,zmm20\t\n"
      "mov r12, rcx\t\n"
      "test r14,r14\t\n"
      "jnz next_inner%=\t\n"
      "add r10,64\t\n"
      "jmp dump_C%=\t\n"

      "zero_regs%=:\t\n"

      "test r14,r14\t\n"
      "jz skip_preload_b_zero%=\t\n"
      "vcvtph2ps zmm31,YMMWORD PTR [r10 + 64]\t\n"
      "skip_preload_b_zero%=:\t\n"
      "vbroadcastss zmm20,DWORD PTR [r9+0]\t\n"
      "vmulps zmm0,zmm21,zmm20\t\n"
      "vmulps zmm1,zmm22,zmm20\t\n"
      "add r12, r13\t\n"
      "vbroadcastss zmm20,DWORD PTR [r9+4]\t\n"
      "vmulps zmm2,zmm21,zmm20\t\n"
      "vmulps zmm3,zmm22,zmm20\t\n"
      "add r12, r13\t\n"
      "vbroadcastss zmm20,DWORD PTR [r9+8]\t\n"
      "vmulps zmm4,zmm21,zmm20\t\n"
      "vmulps zmm5,zmm22,zmm20\t\n"
      "add r12, r13\t\n"
      "vbroadcastss zmm20,DWORD PTR [r9+12]\t\n"
      "vmulps zmm6,zmm21,zmm20\t\n"
      "vmulps zmm7,zmm22,zmm20\t\n"
      "add r12, r13\t\n"
      "vbroadcastss zmm20,DWORD PTR [r9+16]\t\n"
      "vmulps zmm8,zmm21,zmm20\t\n"
      "vmulps zmm9,zmm22,zmm20\t\n"
      "add r12, r13\t\n"
      "vbroadcastss zmm20,DWORD PTR [r9+20]\t\n"
      "vmulps zmm10,zmm21,zmm20\t\n"
      "vmulps zmm11,zmm22,zmm20\t\n"
      "add r12, r13\t\n"
      "vbroadcastss zmm20,DWORD PTR [r9+24]\t\n"
      "vmulps zmm12,zmm21,zmm20\t\n"
      "vmulps zmm13,zmm22,zmm20\t\n"
      "add r12, r13\t\n"
      "vbroadcastss zmm20,DWORD PTR [r9+28]\t\n"
      "vmulps zmm14,zmm21,zmm20\t\n"
      "vmulps zmm15,zmm22,zmm20\t\n"
      "add r12, r13\t\n"
      "vbroadcastss zmm20,DWORD PTR [r9+32]\t\n"
      "vmulps zmm16,zmm21,zmm20\t\n"
      "vmulps zmm17,zmm22,zmm20\t\n"
      "add r12, r13\t\n"
      "vbroadcastss zmm20,DWORD PTR [r9+36]\t\n"
      "vmulps zmm18,zmm21,zmm20\t\n"
      "vmulps zmm19,zmm22,zmm20\t\n"
      "mov r12, rcx\t\n"
      "test r14,r14\t\n"
      "jnz next_inner%=\t\n"
      "add r10,64\t\n"
      "jmp dump_C%=\t\n"

      "loop_inner%=:\t\n"

      "vmovaps zmm21,zmm31\t\n"
      "vcvtph2ps zmm22,YMMWORD PTR [r10 + 32]\t\n"
      "vcvtph2ps zmm31,YMMWORD PTR [r10 + 64]\t\n"
      "vbroadcastss zmm20,DWORD PTR [r9+0]\t\n"
      "vfmadd231ps zmm0,zmm21,zmm20\t\n"
      "vfmadd231ps zmm1,zmm22,zmm20\t\n"
      "vbroadcastss zmm20,DWORD PTR [r9+4]\t\n"
      "vfmadd231ps zmm2,zmm21,zmm20\t\n"
      "vfmadd231ps zmm3,zmm22,zmm20\t\n"
      "vbroadcastss zmm20,DWORD PTR [r9+8]\t\n"
      "vfmadd231ps zmm4,zmm21,zmm20\t\n"
      "vfmadd231ps zmm5,zmm22,zmm20\t\n"
      "vbroadcastss zmm20,DWORD PTR [r9+12]\t\n"
      "vfmadd231ps zmm6,zmm21,zmm20\t\n"
      "vfmadd231ps zmm7,zmm22,zmm20\t\n"
      "vbroadcastss zmm20,DWORD PTR [r9+16]\t\n"
      "vfmadd231ps zmm8,zmm21,zmm20\t\n"
      "vfmadd231ps zmm9,zmm22,zmm20\t\n"
      "vbroadcastss zmm20,DWORD PTR [r9+20]\t\n"
      "vfmadd231ps zmm10,zmm21,zmm20\t\n"
      "vfmadd231ps zmm11,zmm22,zmm20\t\n"
      "vbroadcastss zmm20,DWORD PTR [r9+24]\t\n"
      "vfmadd231ps zmm12,zmm21,zmm20\t\n"
      "vfmadd231ps zmm13,zmm22,zmm20\t\n"
      "vbroadcastss zmm20,DWORD PTR [r9+28]\t\n"
      "vfmadd231ps zmm14,zmm21,zmm20\t\n"
      "vfmadd231ps zmm15,zmm22,zmm20\t\n"
      "vbroadcastss zmm20,DWORD PTR [r9+32]\t\n"
      "vfmadd231ps zmm16,zmm21,zmm20\t\n"
      "vfmadd231ps zmm17,zmm22,zmm20\t\n"
      "vbroadcastss zmm20,DWORD PTR [r9+36]\t\n"
      "vfmadd231ps zmm18,zmm21,zmm20\t\n"
      "vfmadd231ps zmm19,zmm22,zmm20\t\n"

      "next_inner%=:\t\n"
      "add r9,40\t\n"
      "add r10,64\t\n"
      "dec r14\t\n"
      "jnz loop_inner%=\t\n"

      "vmovaps zmm21,zmm31\t\n"
      "vcvtph2ps zmm22,YMMWORD PTR [r10 + 32]\t\n"
      "vbroadcastss zmm20,DWORD PTR [r9+0]\t\n"
      "vfmadd231ps zmm0,zmm21,zmm20\t\n"
      "vfmadd231ps zmm1,zmm22,zmm20\t\n"
      "vbroadcastss zmm20,DWORD PTR [r9+4]\t\n"
      "vfmadd231ps zmm2,zmm21,zmm20\t\n"
      "vfmadd231ps zmm3,zmm22,zmm20\t\n"
      "vbroadcastss zmm20,DWORD PTR [r9+8]\t\n"
      "vfmadd231ps zmm4,zmm21,zmm20\t\n"
      "vfmadd231ps zmm5,zmm22,zmm20\t\n"
      "vbroadcastss zmm20,DWORD PTR [r9+12]\t\n"
      "vfmadd231ps zmm6,zmm21,zmm20\t\n"
      "vfmadd231ps zmm7,zmm22,zmm20\t\n"
      "vbroadcastss zmm20,DWORD PTR [r9+16]\t\n"
      "vfmadd231ps zmm8,zmm21,zmm20\t\n"
      "vfmadd231ps zmm9,zmm22,zmm20\t\n"
      "vbroadcastss zmm20,DWORD PTR [r9+20]\t\n"
      "vfmadd231ps zmm10,zmm21,zmm20\t\n"
      "vfmadd231ps zmm11,zmm22,zmm20\t\n"
      "vbroadcastss zmm20,DWORD PTR [r9+24]\t\n"
      "vfmadd231ps zmm12,zmm21,zmm20\t\n"
      "vfmadd231ps zmm13,zmm22,zmm20\t\n"
      "vbroadcastss zmm20,DWORD PTR [r9+28]\t\n"
      "vfmadd231ps zmm14,zmm21,zmm20\t\n"
      "vfmadd231ps zmm15,zmm22,zmm20\t\n"
      "vbroadcastss zmm20,DWORD PTR [r9+32]\t\n"
      "vfmadd231ps zmm16,zmm21,zmm20\t\n"
      "vfmadd231ps zmm17,zmm22,zmm20\t\n"
      "vbroadcastss zmm20,DWORD PTR [r9+36]\t\n"
      "vfmadd231ps zmm18,zmm21,zmm20\t\n"
      "vfmadd231ps zmm19,zmm22,zmm20\t\n"
      "add r9,40\t\n"
      "add r10,64\t\n"
      // Dump C
      "dump_C%=:\t\n"
      "vmovups zmmword PTR [r12 + 0], zmm0\t\n"
      "vmovups zmmword PTR [r12 + 64], zmm1\t\n"
      "add r12, r13\t\n"
      "vmovups zmmword PTR [r12 + 0], zmm2\t\n"
      "vmovups zmmword PTR [r12 + 64], zmm3\t\n"
      "add r12, r13\t\n"
      "vmovups zmmword PTR [r12 + 0], zmm4\t\n"
      "vmovups zmmword PTR [r12 + 64], zmm5\t\n"
      "add r12, r13\t\n"
      "vmovups zmmword PTR [r12 + 0], zmm6\t\n"
      "vmovups zmmword PTR [r12 + 64], zmm7\t\n"
      "add r12, r13\t\n"
      "vmovups zmmword PTR [r12 + 0], zmm8\t\n"
      "vmovups zmmword PTR [r12 + 64], zmm9\t\n"
      "add r12, r13\t\n"
      "vmovups zmmword PTR [r12 + 0], zmm10\t\n"
      "vmovups zmmword PTR [r12 + 64], zmm11\t\n"
      "add r12, r13\t\n"
      "vmovups zmmword PTR [r12 + 0], zmm12\t\n"
      "vmovups zmmword PTR [r12 + 64], zmm13\t\n"
      "add r12, r13\t\n"
      "vmovups zmmword PTR [r12 + 0], zmm14\t\n"
      "vmovups zmmword PTR [r12 + 64], zmm15\t\n"
      "add r12, r13\t\n"
      "vmovups zmmword PTR [r12 + 0], zmm16\t\n"
      "vmovups zmmword PTR [r12 + 64], zmm17\t\n"
      "add r12, r13\t\n"
      "vmovups zmmword PTR [r12 + 0], zmm18\t\n"
      "vmovups zmmword PTR [r12 + 64], zmm19\t\n"

      // next outer iteration
      "add rcx, 128\t\n"
      "mov r12, rcx\t\n"
      "mov r9, rax\t\n"
      "inc rbx\t\n"
      "cmp rbx, rdi\t\n"
      "jl loop_outter%=\t\n"
      :
      : [gp] "rm"(gp)
      : "r8",
        "r9",
        "r10",
        "r11",
        "r13",
        "r14",
        "rax",
        "rcx",
        "rsi",
        "rdi",
        "rbx",
        "r12",
        "r15",
        "memory");
}
void NOINLINE
gemmkernel_11x2_Avx512_fp16_fA0fB0fC0(GemmParamsFP16* gp) {
  asm volatile(
#if !defined(__clang__)
      "mov r14, %[gp]\t\n"
#else
      "mov %[gp], %%r14\t\n"
      ".intel_syntax noprefix\t\n"
#endif

      // Copy parameters
      // k
      "mov r8, [r14 + 0]\t\n"
      "dec r8\t\n"
      // A
      "mov r9, [r14 + 8]\t\n"
      // B
      "mov r10, [r14 + 16]\t\n"
      // beta
      "lea r15, [r14 + 24]\t\n"
      // C
      "mov r12, [r14 + 32]\t\n"
      // ldc
      "mov r13, [r14 + 40]\t\n"
      // b_block_cols
      "mov rdi, [r14 + 48]\t\n"
      // b_block_size
      "mov rsi, [r14 + 56]\t\n"

      // Make copies of A and C
      "mov rax, r9\t\n"
      "mov rcx, r12\t\n"

      "mov rbx, 0\t\n"
      "loop_outter%=:\t\n"
      "mov r14, r8\t\n"
      "vbroadcastss zmm31,DWORD PTR [r15]\t\n"
      "vcvtph2ps zmm23,YMMWORD PTR [r10 + 0]\t\n"
      "vcvtph2ps zmm24,YMMWORD PTR [r10 + 32]\t\n"
      "vxorps xmm0, xmm0, xmm0\t\n"
      "vcomiss xmm31, xmm0\t\n"
      "jz zero_regs%=\t\n"

      // Setup values with beta multiplication
      "vmulps zmm0, zmm31, [r12 + 0]\t\n"
      "vmulps zmm1, zmm31, [r12 + 64]\t\n"
      "add r12, r13\t\n"
      "vmulps zmm2, zmm31, [r12 + 0]\t\n"
      "vmulps zmm3, zmm31, [r12 + 64]\t\n"
      "add r12, r13\t\n"
      "vmulps zmm4, zmm31, [r12 + 0]\t\n"
      "vmulps zmm5, zmm31, [r12 + 64]\t\n"
      "add r12, r13\t\n"
      "vmulps zmm6, zmm31, [r12 + 0]\t\n"
      "vmulps zmm7, zmm31, [r12 + 64]\t\n"
      "add r12, r13\t\n"
      "vmulps zmm8, zmm31, [r12 + 0]\t\n"
      "vmulps zmm9, zmm31, [r12 + 64]\t\n"
      "add r12, r13\t\n"
      "vmulps zmm10, zmm31, [r12 + 0]\t\n"
      "vmulps zmm11, zmm31, [r12 + 64]\t\n"
      "add r12, r13\t\n"
      "vmulps zmm12, zmm31, [r12 + 0]\t\n"
      "vmulps zmm13, zmm31, [r12 + 64]\t\n"
      "add r12, r13\t\n"
      "vmulps zmm14, zmm31, [r12 + 0]\t\n"
      "vmulps zmm15, zmm31, [r12 + 64]\t\n"
      "add r12, r13\t\n"
      "vmulps zmm16, zmm31, [r12 + 0]\t\n"
      "vmulps zmm17, zmm31, [r12 + 64]\t\n"
      "add r12, r13\t\n"
      "vmulps zmm18, zmm31, [r12 + 0]\t\n"
      "vmulps zmm19, zmm31, [r12 + 64]\t\n"
      "add r12, r13\t\n"
      "vmulps zmm20, zmm31, [r12 + 0]\t\n"
      "vmulps zmm21, zmm31, [r12 + 64]\t\n"
      "test r14,r14\t\n"
      "jz skip_preload%=\t\n"
      "vcvtph2ps zmm31,YMMWORD PTR [r10 + 64]\t\n"
      "skip_preload%=:\t\n"
      "vbroadcastss zmm22,DWORD PTR [r9+0]\t\n"
      "vfmadd231ps zmm0,zmm23,zmm22\t\n"
      "vfmadd231ps zmm1,zmm24,zmm22\t\n"
      "vbroadcastss zmm22,DWORD PTR [r9+4]\t\n"
      "vfmadd231ps zmm2,zmm23,zmm22\t\n"
      "vfmadd231ps zmm3,zmm24,zmm22\t\n"
      "vbroadcastss zmm22,DWORD PTR [r9+8]\t\n"
      "vfmadd231ps zmm4,zmm23,zmm22\t\n"
      "vfmadd231ps zmm5,zmm24,zmm22\t\n"
      "vbroadcastss zmm22,DWORD PTR [r9+12]\t\n"
      "vfmadd231ps zmm6,zmm23,zmm22\t\n"
      "vfmadd231ps zmm7,zmm24,zmm22\t\n"
      "vbroadcastss zmm22,DWORD PTR [r9+16]\t\n"
      "vfmadd231ps zmm8,zmm23,zmm22\t\n"
      "vfmadd231ps zmm9,zmm24,zmm22\t\n"
      "vbroadcastss zmm22,DWORD PTR [r9+20]\t\n"
      "vfmadd231ps zmm10,zmm23,zmm22\t\n"
      "vfmadd231ps zmm11,zmm24,zmm22\t\n"
      "vbroadcastss zmm22,DWORD PTR [r9+24]\t\n"
      "vfmadd231ps zmm12,zmm23,zmm22\t\n"
      "vfmadd231ps zmm13,zmm24,zmm22\t\n"
      "vbroadcastss zmm22,DWORD PTR [r9+28]\t\n"
      "vfmadd231ps zmm14,zmm23,zmm22\t\n"
      "vfmadd231ps zmm15,zmm24,zmm22\t\n"
      "vbroadcastss zmm22,DWORD PTR [r9+32]\t\n"
      "vfmadd231ps zmm16,zmm23,zmm22\t\n"
      "vfmadd231ps zmm17,zmm24,zmm22\t\n"
      "vbroadcastss zmm22,DWORD PTR [r9+36]\t\n"
      "vfmadd231ps zmm18,zmm23,zmm22\t\n"
      "vfmadd231ps zmm19,zmm24,zmm22\t\n"
      "vbroadcastss zmm22,DWORD PTR [r9+40]\t\n"
      "vfmadd231ps zmm20,zmm23,zmm22\t\n"
      "vfmadd231ps zmm21,zmm24,zmm22\t\n"
      "mov r12, rcx\t\n"
      "test r14,r14\t\n"
      "jnz next_inner%=\t\n"
      "add r10,64\t\n"
      "jmp dump_C%=\t\n"

      "zero_regs%=:\t\n"

      "test r14,r14\t\n"
      "jz skip_preload_b_zero%=\t\n"
      "vcvtph2ps zmm31,YMMWORD PTR [r10 + 64]\t\n"
      "skip_preload_b_zero%=:\t\n"
      "vbroadcastss zmm22,DWORD PTR [r9+0]\t\n"
      "vmulps zmm0,zmm23,zmm22\t\n"
      "vmulps zmm1,zmm24,zmm22\t\n"
      "add r12, r13\t\n"
      "vbroadcastss zmm22,DWORD PTR [r9+4]\t\n"
      "vmulps zmm2,zmm23,zmm22\t\n"
      "vmulps zmm3,zmm24,zmm22\t\n"
      "add r12, r13\t\n"
      "vbroadcastss zmm22,DWORD PTR [r9+8]\t\n"
      "vmulps zmm4,zmm23,zmm22\t\n"
      "vmulps zmm5,zmm24,zmm22\t\n"
      "add r12, r13\t\n"
      "vbroadcastss zmm22,DWORD PTR [r9+12]\t\n"
      "vmulps zmm6,zmm23,zmm22\t\n"
      "vmulps zmm7,zmm24,zmm22\t\n"
      "add r12, r13\t\n"
      "vbroadcastss zmm22,DWORD PTR [r9+16]\t\n"
      "vmulps zmm8,zmm23,zmm22\t\n"
      "vmulps zmm9,zmm24,zmm22\t\n"
      "add r12, r13\t\n"
      "vbroadcastss zmm22,DWORD PTR [r9+20]\t\n"
      "vmulps zmm10,zmm23,zmm22\t\n"
      "vmulps zmm11,zmm24,zmm22\t\n"
      "add r12, r13\t\n"
      "vbroadcastss zmm22,DWORD PTR [r9+24]\t\n"
      "vmulps zmm12,zmm23,zmm22\t\n"
      "vmulps zmm13,zmm24,zmm22\t\n"
      "add r12, r13\t\n"
      "vbroadcastss zmm22,DWORD PTR [r9+28]\t\n"
      "vmulps zmm14,zmm23,zmm22\t\n"
      "vmulps zmm15,zmm24,zmm22\t\n"
      "add r12, r13\t\n"
      "vbroadcastss zmm22,DWORD PTR [r9+32]\t\n"
      "vmulps zmm16,zmm23,zmm22\t\n"
      "vmulps zmm17,zmm24,zmm22\t\n"
      "add r12, r13\t\n"
      "vbroadcastss zmm22,DWORD PTR [r9+36]\t\n"
      "vmulps zmm18,zmm23,zmm22\t\n"
      "vmulps zmm19,zmm24,zmm22\t\n"
      "add r12, r13\t\n"
      "vbroadcastss zmm22,DWORD PTR [r9+40]\t\n"
      "vmulps zmm20,zmm23,zmm22\t\n"
      "vmulps zmm21,zmm24,zmm22\t\n"
      "mov r12, rcx\t\n"
      "test r14,r14\t\n"
      "jnz next_inner%=\t\n"
      "add r10,64\t\n"
      "jmp dump_C%=\t\n"

      "loop_inner%=:\t\n"

      "vmovaps zmm23,zmm31\t\n"
      "vcvtph2ps zmm24,YMMWORD PTR [r10 + 32]\t\n"
      "vcvtph2ps zmm31,YMMWORD PTR [r10 + 64]\t\n"
      "vbroadcastss zmm22,DWORD PTR [r9+0]\t\n"
      "vfmadd231ps zmm0,zmm23,zmm22\t\n"
      "vfmadd231ps zmm1,zmm24,zmm22\t\n"
      "vbroadcastss zmm22,DWORD PTR [r9+4]\t\n"
      "vfmadd231ps zmm2,zmm23,zmm22\t\n"
      "vfmadd231ps zmm3,zmm24,zmm22\t\n"
      "vbroadcastss zmm22,DWORD PTR [r9+8]\t\n"
      "vfmadd231ps zmm4,zmm23,zmm22\t\n"
      "vfmadd231ps zmm5,zmm24,zmm22\t\n"
      "vbroadcastss zmm22,DWORD PTR [r9+12]\t\n"
      "vfmadd231ps zmm6,zmm23,zmm22\t\n"
      "vfmadd231ps zmm7,zmm24,zmm22\t\n"
      "vbroadcastss zmm22,DWORD PTR [r9+16]\t\n"
      "vfmadd231ps zmm8,zmm23,zmm22\t\n"
      "vfmadd231ps zmm9,zmm24,zmm22\t\n"
      "vbroadcastss zmm22,DWORD PTR [r9+20]\t\n"
      "vfmadd231ps zmm10,zmm23,zmm22\t\n"
      "vfmadd231ps zmm11,zmm24,zmm22\t\n"
      "vbroadcastss zmm22,DWORD PTR [r9+24]\t\n"
      "vfmadd231ps zmm12,zmm23,zmm22\t\n"
      "vfmadd231ps zmm13,zmm24,zmm22\t\n"
      "vbroadcastss zmm22,DWORD PTR [r9+28]\t\n"
      "vfmadd231ps zmm14,zmm23,zmm22\t\n"
      "vfmadd231ps zmm15,zmm24,zmm22\t\n"
      "vbroadcastss zmm22,DWORD PTR [r9+32]\t\n"
      "vfmadd231ps zmm16,zmm23,zmm22\t\n"
      "vfmadd231ps zmm17,zmm24,zmm22\t\n"
      "vbroadcastss zmm22,DWORD PTR [r9+36]\t\n"
      "vfmadd231ps zmm18,zmm23,zmm22\t\n"
      "vfmadd231ps zmm19,zmm24,zmm22\t\n"
      "vbroadcastss zmm22,DWORD PTR [r9+40]\t\n"
      "vfmadd231ps zmm20,zmm23,zmm22\t\n"
      "vfmadd231ps zmm21,zmm24,zmm22\t\n"

      "next_inner%=:\t\n"
      "add r9,44\t\n"
      "add r10,64\t\n"
      "dec r14\t\n"
      "jnz loop_inner%=\t\n"

      "vmovaps zmm23,zmm31\t\n"
      "vcvtph2ps zmm24,YMMWORD PTR [r10 + 32]\t\n"
      "vbroadcastss zmm22,DWORD PTR [r9+0]\t\n"
      "vfmadd231ps zmm0,zmm23,zmm22\t\n"
      "vfmadd231ps zmm1,zmm24,zmm22\t\n"
      "vbroadcastss zmm22,DWORD PTR [r9+4]\t\n"
      "vfmadd231ps zmm2,zmm23,zmm22\t\n"
      "vfmadd231ps zmm3,zmm24,zmm22\t\n"
      "vbroadcastss zmm22,DWORD PTR [r9+8]\t\n"
      "vfmadd231ps zmm4,zmm23,zmm22\t\n"
      "vfmadd231ps zmm5,zmm24,zmm22\t\n"
      "vbroadcastss zmm22,DWORD PTR [r9+12]\t\n"
      "vfmadd231ps zmm6,zmm23,zmm22\t\n"
      "vfmadd231ps zmm7,zmm24,zmm22\t\n"
      "vbroadcastss zmm22,DWORD PTR [r9+16]\t\n"
      "vfmadd231ps zmm8,zmm23,zmm22\t\n"
      "vfmadd231ps zmm9,zmm24,zmm22\t\n"
      "vbroadcastss zmm22,DWORD PTR [r9+20]\t\n"
      "vfmadd231ps zmm10,zmm23,zmm22\t\n"
      "vfmadd231ps zmm11,zmm24,zmm22\t\n"
      "vbroadcastss zmm22,DWORD PTR [r9+24]\t\n"
      "vfmadd231ps zmm12,zmm23,zmm22\t\n"
      "vfmadd231ps zmm13,zmm24,zmm22\t\n"
      "vbroadcastss zmm22,DWORD PTR [r9+28]\t\n"
      "vfmadd231ps zmm14,zmm23,zmm22\t\n"
      "vfmadd231ps zmm15,zmm24,zmm22\t\n"
      "vbroadcastss zmm22,DWORD PTR [r9+32]\t\n"
      "vfmadd231ps zmm16,zmm23,zmm22\t\n"
      "vfmadd231ps zmm17,zmm24,zmm22\t\n"
      "vbroadcastss zmm22,DWORD PTR [r9+36]\t\n"
      "vfmadd231ps zmm18,zmm23,zmm22\t\n"
      "vfmadd231ps zmm19,zmm24,zmm22\t\n"
      "vbroadcastss zmm22,DWORD PTR [r9+40]\t\n"
      "vfmadd231ps zmm20,zmm23,zmm22\t\n"
      "vfmadd231ps zmm21,zmm24,zmm22\t\n"
      "add r9,44\t\n"
      "add r10,64\t\n"
      // Dump C
      "dump_C%=:\t\n"
      "vmovups zmmword PTR [r12 + 0], zmm0\t\n"
      "vmovups zmmword PTR [r12 + 64], zmm1\t\n"
      "add r12, r13\t\n"
      "vmovups zmmword PTR [r12 + 0], zmm2\t\n"
      "vmovups zmmword PTR [r12 + 64], zmm3\t\n"
      "add r12, r13\t\n"
      "vmovups zmmword PTR [r12 + 0], zmm4\t\n"
      "vmovups zmmword PTR [r12 + 64], zmm5\t\n"
      "add r12, r13\t\n"
      "vmovups zmmword PTR [r12 + 0], zmm6\t\n"
      "vmovups zmmword PTR [r12 + 64], zmm7\t\n"
      "add r12, r13\t\n"
      "vmovups zmmword PTR [r12 + 0], zmm8\t\n"
      "vmovups zmmword PTR [r12 + 64], zmm9\t\n"
      "add r12, r13\t\n"
      "vmovups zmmword PTR [r12 + 0], zmm10\t\n"
      "vmovups zmmword PTR [r12 + 64], zmm11\t\n"
      "add r12, r13\t\n"
      "vmovups zmmword PTR [r12 + 0], zmm12\t\n"
      "vmovups zmmword PTR [r12 + 64], zmm13\t\n"
      "add r12, r13\t\n"
      "vmovups zmmword PTR [r12 + 0], zmm14\t\n"
      "vmovups zmmword PTR [r12 + 64], zmm15\t\n"
      "add r12, r13\t\n"
      "vmovups zmmword PTR [r12 + 0], zmm16\t\n"
      "vmovups zmmword PTR [r12 + 64], zmm17\t\n"
      "add r12, r13\t\n"
      "vmovups zmmword PTR [r12 + 0], zmm18\t\n"
      "vmovups zmmword PTR [r12 + 64], zmm19\t\n"
      "add r12, r13\t\n"
      "vmovups zmmword PTR [r12 + 0], zmm20\t\n"
      "vmovups zmmword PTR [r12 + 64], zmm21\t\n"

      // next outer iteration
      "add rcx, 128\t\n"
      "mov r12, rcx\t\n"
      "mov r9, rax\t\n"
      "inc rbx\t\n"
      "cmp rbx, rdi\t\n"
      "jl loop_outter%=\t\n"
      :
      : [gp] "rm"(gp)
      : "r8",
        "r9",
        "r10",
        "r11",
        "r13",
        "r14",
        "rax",
        "rcx",
        "rsi",
        "rdi",
        "rbx",
        "r12",
        "r15",
        "memory");
}
void NOINLINE
gemmkernel_12x2_Avx512_fp16_fA0fB0fC0(GemmParamsFP16* gp) {
  asm volatile(
#if !defined(__clang__)
      "mov r14, %[gp]\t\n"
#else
      "mov %[gp], %%r14\t\n"
      ".intel_syntax noprefix\t\n"
#endif

      // Copy parameters
      // k
      "mov r8, [r14 + 0]\t\n"
      "dec r8\t\n"
      // A
      "mov r9, [r14 + 8]\t\n"
      // B
      "mov r10, [r14 + 16]\t\n"
      // beta
      "lea r15, [r14 + 24]\t\n"
      // C
      "mov r12, [r14 + 32]\t\n"
      // ldc
      "mov r13, [r14 + 40]\t\n"
      // b_block_cols
      "mov rdi, [r14 + 48]\t\n"
      // b_block_size
      "mov rsi, [r14 + 56]\t\n"

      // Make copies of A and C
      "mov rax, r9\t\n"
      "mov rcx, r12\t\n"

      "mov rbx, 0\t\n"
      "loop_outter%=:\t\n"
      "mov r14, r8\t\n"
      "vbroadcastss zmm31,DWORD PTR [r15]\t\n"
      "vcvtph2ps zmm25,YMMWORD PTR [r10 + 0]\t\n"
      "vcvtph2ps zmm26,YMMWORD PTR [r10 + 32]\t\n"
      "vxorps xmm0, xmm0, xmm0\t\n"
      "vcomiss xmm31, xmm0\t\n"
      "jz zero_regs%=\t\n"

      // Setup values with beta multiplication
      "vmulps zmm0, zmm31, [r12 + 0]\t\n"
      "vmulps zmm1, zmm31, [r12 + 64]\t\n"
      "add r12, r13\t\n"
      "vmulps zmm2, zmm31, [r12 + 0]\t\n"
      "vmulps zmm3, zmm31, [r12 + 64]\t\n"
      "add r12, r13\t\n"
      "vmulps zmm4, zmm31, [r12 + 0]\t\n"
      "vmulps zmm5, zmm31, [r12 + 64]\t\n"
      "add r12, r13\t\n"
      "vmulps zmm6, zmm31, [r12 + 0]\t\n"
      "vmulps zmm7, zmm31, [r12 + 64]\t\n"
      "add r12, r13\t\n"
      "vmulps zmm8, zmm31, [r12 + 0]\t\n"
      "vmulps zmm9, zmm31, [r12 + 64]\t\n"
      "add r12, r13\t\n"
      "vmulps zmm10, zmm31, [r12 + 0]\t\n"
      "vmulps zmm11, zmm31, [r12 + 64]\t\n"
      "add r12, r13\t\n"
      "vmulps zmm12, zmm31, [r12 + 0]\t\n"
      "vmulps zmm13, zmm31, [r12 + 64]\t\n"
      "add r12, r13\t\n"
      "vmulps zmm14, zmm31, [r12 + 0]\t\n"
      "vmulps zmm15, zmm31, [r12 + 64]\t\n"
      "add r12, r13\t\n"
      "vmulps zmm16, zmm31, [r12 + 0]\t\n"
      "vmulps zmm17, zmm31, [r12 + 64]\t\n"
      "add r12, r13\t\n"
      "vmulps zmm18, zmm31, [r12 + 0]\t\n"
      "vmulps zmm19, zmm31, [r12 + 64]\t\n"
      "add r12, r13\t\n"
      "vmulps zmm20, zmm31, [r12 + 0]\t\n"
      "vmulps zmm21, zmm31, [r12 + 64]\t\n"
      "add r12, r13\t\n"
      "vmulps zmm22, zmm31, [r12 + 0]\t\n"
      "vmulps zmm23, zmm31, [r12 + 64]\t\n"
      "test r14,r14\t\n"
      "jz skip_preload%=\t\n"
      "vcvtph2ps zmm31,YMMWORD PTR [r10 + 64]\t\n"
      "skip_preload%=:\t\n"
      "vbroadcastss zmm24,DWORD PTR [r9+0]\t\n"
      "vfmadd231ps zmm0,zmm25,zmm24\t\n"
      "vfmadd231ps zmm1,zmm26,zmm24\t\n"
      "vbroadcastss zmm24,DWORD PTR [r9+4]\t\n"
      "vfmadd231ps zmm2,zmm25,zmm24\t\n"
      "vfmadd231ps zmm3,zmm26,zmm24\t\n"
      "vbroadcastss zmm24,DWORD PTR [r9+8]\t\n"
      "vfmadd231ps zmm4,zmm25,zmm24\t\n"
      "vfmadd231ps zmm5,zmm26,zmm24\t\n"
      "vbroadcastss zmm24,DWORD PTR [r9+12]\t\n"
      "vfmadd231ps zmm6,zmm25,zmm24\t\n"
      "vfmadd231ps zmm7,zmm26,zmm24\t\n"
      "vbroadcastss zmm24,DWORD PTR [r9+16]\t\n"
      "vfmadd231ps zmm8,zmm25,zmm24\t\n"
      "vfmadd231ps zmm9,zmm26,zmm24\t\n"
      "vbroadcastss zmm24,DWORD PTR [r9+20]\t\n"
      "vfmadd231ps zmm10,zmm25,zmm24\t\n"
      "vfmadd231ps zmm11,zmm26,zmm24\t\n"
      "vbroadcastss zmm24,DWORD PTR [r9+24]\t\n"
      "vfmadd231ps zmm12,zmm25,zmm24\t\n"
      "vfmadd231ps zmm13,zmm26,zmm24\t\n"
      "vbroadcastss zmm24,DWORD PTR [r9+28]\t\n"
      "vfmadd231ps zmm14,zmm25,zmm24\t\n"
      "vfmadd231ps zmm15,zmm26,zmm24\t\n"
      "vbroadcastss zmm24,DWORD PTR [r9+32]\t\n"
      "vfmadd231ps zmm16,zmm25,zmm24\t\n"
      "vfmadd231ps zmm17,zmm26,zmm24\t\n"
      "vbroadcastss zmm24,DWORD PTR [r9+36]\t\n"
      "vfmadd231ps zmm18,zmm25,zmm24\t\n"
      "vfmadd231ps zmm19,zmm26,zmm24\t\n"
      "vbroadcastss zmm24,DWORD PTR [r9+40]\t\n"
      "vfmadd231ps zmm20,zmm25,zmm24\t\n"
      "vfmadd231ps zmm21,zmm26,zmm24\t\n"
      "vbroadcastss zmm24,DWORD PTR [r9+44]\t\n"
      "vfmadd231ps zmm22,zmm25,zmm24\t\n"
      "vfmadd231ps zmm23,zmm26,zmm24\t\n"
      "mov r12, rcx\t\n"
      "test r14,r14\t\n"
      "jnz next_inner%=\t\n"
      "add r10,64\t\n"
      "jmp dump_C%=\t\n"

      "zero_regs%=:\t\n"

      "test r14,r14\t\n"
      "jz skip_preload_b_zero%=\t\n"
      "vcvtph2ps zmm31,YMMWORD PTR [r10 + 64]\t\n"
      "skip_preload_b_zero%=:\t\n"
      "vbroadcastss zmm24,DWORD PTR [r9+0]\t\n"
      "vmulps zmm0,zmm25,zmm24\t\n"
      "vmulps zmm1,zmm26,zmm24\t\n"
      "add r12, r13\t\n"
      "vbroadcastss zmm24,DWORD PTR [r9+4]\t\n"
      "vmulps zmm2,zmm25,zmm24\t\n"
      "vmulps zmm3,zmm26,zmm24\t\n"
      "add r12, r13\t\n"
      "vbroadcastss zmm24,DWORD PTR [r9+8]\t\n"
      "vmulps zmm4,zmm25,zmm24\t\n"
      "vmulps zmm5,zmm26,zmm24\t\n"
      "add r12, r13\t\n"
      "vbroadcastss zmm24,DWORD PTR [r9+12]\t\n"
      "vmulps zmm6,zmm25,zmm24\t\n"
      "vmulps zmm7,zmm26,zmm24\t\n"
      "add r12, r13\t\n"
      "vbroadcastss zmm24,DWORD PTR [r9+16]\t\n"
      "vmulps zmm8,zmm25,zmm24\t\n"
      "vmulps zmm9,zmm26,zmm24\t\n"
      "add r12, r13\t\n"
      "vbroadcastss zmm24,DWORD PTR [r9+20]\t\n"
      "vmulps zmm10,zmm25,zmm24\t\n"
      "vmulps zmm11,zmm26,zmm24\t\n"
      "add r12, r13\t\n"
      "vbroadcastss zmm24,DWORD PTR [r9+24]\t\n"
      "vmulps zmm12,zmm25,zmm24\t\n"
      "vmulps zmm13,zmm26,zmm24\t\n"
      "add r12, r13\t\n"
      "vbroadcastss zmm24,DWORD PTR [r9+28]\t\n"
      "vmulps zmm14,zmm25,zmm24\t\n"
      "vmulps zmm15,zmm26,zmm24\t\n"
      "add r12, r13\t\n"
      "vbroadcastss zmm24,DWORD PTR [r9+32]\t\n"
      "vmulps zmm16,zmm25,zmm24\t\n"
      "vmulps zmm17,zmm26,zmm24\t\n"
      "add r12, r13\t\n"
      "vbroadcastss zmm24,DWORD PTR [r9+36]\t\n"
      "vmulps zmm18,zmm25,zmm24\t\n"
      "vmulps zmm19,zmm26,zmm24\t\n"
      "add r12, r13\t\n"
      "vbroadcastss zmm24,DWORD PTR [r9+40]\t\n"
      "vmulps zmm20,zmm25,zmm24\t\n"
      "vmulps zmm21,zmm26,zmm24\t\n"
      "add r12, r13\t\n"
      "vbroadcastss zmm24,DWORD PTR [r9+44]\t\n"
      "vmulps zmm22,zmm25,zmm24\t\n"
      "vmulps zmm23,zmm26,zmm24\t\n"
      "mov r12, rcx\t\n"
      "test r14,r14\t\n"
      "jnz next_inner%=\t\n"
      "add r10,64\t\n"
      "jmp dump_C%=\t\n"

      "loop_inner%=:\t\n"

      "vmovaps zmm25,zmm31\t\n"
      "vcvtph2ps zmm26,YMMWORD PTR [r10 + 32]\t\n"
      "vcvtph2ps zmm31,YMMWORD PTR [r10 + 64]\t\n"
      "vbroadcastss zmm24,DWORD PTR [r9+0]\t\n"
      "vfmadd231ps zmm0,zmm25,zmm24\t\n"
      "vfmadd231ps zmm1,zmm26,zmm24\t\n"
      "vbroadcastss zmm24,DWORD PTR [r9+4]\t\n"
      "vfmadd231ps zmm2,zmm25,zmm24\t\n"
      "vfmadd231ps zmm3,zmm26,zmm24\t\n"
      "vbroadcastss zmm24,DWORD PTR [r9+8]\t\n"
      "vfmadd231ps zmm4,zmm25,zmm24\t\n"
      "vfmadd231ps zmm5,zmm26,zmm24\t\n"
      "vbroadcastss zmm24,DWORD PTR [r9+12]\t\n"
      "vfmadd231ps zmm6,zmm25,zmm24\t\n"
      "vfmadd231ps zmm7,zmm26,zmm24\t\n"
      "vbroadcastss zmm24,DWORD PTR [r9+16]\t\n"
      "vfmadd231ps zmm8,zmm25,zmm24\t\n"
      "vfmadd231ps zmm9,zmm26,zmm24\t\n"
      "vbroadcastss zmm24,DWORD PTR [r9+20]\t\n"
      "vfmadd231ps zmm10,zmm25,zmm24\t\n"
      "vfmadd231ps zmm11,zmm26,zmm24\t\n"
      "vbroadcastss zmm24,DWORD PTR [r9+24]\t\n"
      "vfmadd231ps zmm12,zmm25,zmm24\t\n"
      "vfmadd231ps zmm13,zmm26,zmm24\t\n"
      "vbroadcastss zmm24,DWORD PTR [r9+28]\t\n"
      "vfmadd231ps zmm14,zmm25,zmm24\t\n"
      "vfmadd231ps zmm15,zmm26,zmm24\t\n"
      "vbroadcastss zmm24,DWORD PTR [r9+32]\t\n"
      "vfmadd231ps zmm16,zmm25,zmm24\t\n"
      "vfmadd231ps zmm17,zmm26,zmm24\t\n"
      "vbroadcastss zmm24,DWORD PTR [r9+36]\t\n"
      "vfmadd231ps zmm18,zmm25,zmm24\t\n"
      "vfmadd231ps zmm19,zmm26,zmm24\t\n"
      "vbroadcastss zmm24,DWORD PTR [r9+40]\t\n"
      "vfmadd231ps zmm20,zmm25,zmm24\t\n"
      "vfmadd231ps zmm21,zmm26,zmm24\t\n"
      "vbroadcastss zmm24,DWORD PTR [r9+44]\t\n"
      "vfmadd231ps zmm22,zmm25,zmm24\t\n"
      "vfmadd231ps zmm23,zmm26,zmm24\t\n"

      "next_inner%=:\t\n"
      "add r9,48\t\n"
      "add r10,64\t\n"
      "dec r14\t\n"
      "jnz loop_inner%=\t\n"

      "vmovaps zmm25,zmm31\t\n"
      "vcvtph2ps zmm26,YMMWORD PTR [r10 + 32]\t\n"
      "vbroadcastss zmm24,DWORD PTR [r9+0]\t\n"
      "vfmadd231ps zmm0,zmm25,zmm24\t\n"
      "vfmadd231ps zmm1,zmm26,zmm24\t\n"
      "vbroadcastss zmm24,DWORD PTR [r9+4]\t\n"
      "vfmadd231ps zmm2,zmm25,zmm24\t\n"
      "vfmadd231ps zmm3,zmm26,zmm24\t\n"
      "vbroadcastss zmm24,DWORD PTR [r9+8]\t\n"
      "vfmadd231ps zmm4,zmm25,zmm24\t\n"
      "vfmadd231ps zmm5,zmm26,zmm24\t\n"
      "vbroadcastss zmm24,DWORD PTR [r9+12]\t\n"
      "vfmadd231ps zmm6,zmm25,zmm24\t\n"
      "vfmadd231ps zmm7,zmm26,zmm24\t\n"
      "vbroadcastss zmm24,DWORD PTR [r9+16]\t\n"
      "vfmadd231ps zmm8,zmm25,zmm24\t\n"
      "vfmadd231ps zmm9,zmm26,zmm24\t\n"
      "vbroadcastss zmm24,DWORD PTR [r9+20]\t\n"
      "vfmadd231ps zmm10,zmm25,zmm24\t\n"
      "vfmadd231ps zmm11,zmm26,zmm24\t\n"
      "vbroadcastss zmm24,DWORD PTR [r9+24]\t\n"
      "vfmadd231ps zmm12,zmm25,zmm24\t\n"
      "vfmadd231ps zmm13,zmm26,zmm24\t\n"
      "vbroadcastss zmm24,DWORD PTR [r9+28]\t\n"
      "vfmadd231ps zmm14,zmm25,zmm24\t\n"
      "vfmadd231ps zmm15,zmm26,zmm24\t\n"
      "vbroadcastss zmm24,DWORD PTR [r9+32]\t\n"
      "vfmadd231ps zmm16,zmm25,zmm24\t\n"
      "vfmadd231ps zmm17,zmm26,zmm24\t\n"
      "vbroadcastss zmm24,DWORD PTR [r9+36]\t\n"
      "vfmadd231ps zmm18,zmm25,zmm24\t\n"
      "vfmadd231ps zmm19,zmm26,zmm24\t\n"
      "vbroadcastss zmm24,DWORD PTR [r9+40]\t\n"
      "vfmadd231ps zmm20,zmm25,zmm24\t\n"
      "vfmadd231ps zmm21,zmm26,zmm24\t\n"
      "vbroadcastss zmm24,DWORD PTR [r9+44]\t\n"
      "vfmadd231ps zmm22,zmm25,zmm24\t\n"
      "vfmadd231ps zmm23,zmm26,zmm24\t\n"
      "add r9,48\t\n"
      "add r10,64\t\n"
      // Dump C
      "dump_C%=:\t\n"
      "vmovups zmmword PTR [r12 + 0], zmm0\t\n"
      "vmovups zmmword PTR [r12 + 64], zmm1\t\n"
      "add r12, r13\t\n"
      "vmovups zmmword PTR [r12 + 0], zmm2\t\n"
      "vmovups zmmword PTR [r12 + 64], zmm3\t\n"
      "add r12, r13\t\n"
      "vmovups zmmword PTR [r12 + 0], zmm4\t\n"
      "vmovups zmmword PTR [r12 + 64], zmm5\t\n"
      "add r12, r13\t\n"
      "vmovups zmmword PTR [r12 + 0], zmm6\t\n"
      "vmovups zmmword PTR [r12 + 64], zmm7\t\n"
      "add r12, r13\t\n"
      "vmovups zmmword PTR [r12 + 0], zmm8\t\n"
      "vmovups zmmword PTR [r12 + 64], zmm9\t\n"
      "add r12, r13\t\n"
      "vmovups zmmword PTR [r12 + 0], zmm10\t\n"
      "vmovups zmmword PTR [r12 + 64], zmm11\t\n"
      "add r12, r13\t\n"
      "vmovups zmmword PTR [r12 + 0], zmm12\t\n"
      "vmovups zmmword PTR [r12 + 64], zmm13\t\n"
      "add r12, r13\t\n"
      "vmovups zmmword PTR [r12 + 0], zmm14\t\n"
      "vmovups zmmword PTR [r12 + 64], zmm15\t\n"
      "add r12, r13\t\n"
      "vmovups zmmword PTR [r12 + 0], zmm16\t\n"
      "vmovups zmmword PTR [r12 + 64], zmm17\t\n"
      "add r12, r13\t\n"
      "vmovups zmmword PTR [r12 + 0], zmm18\t\n"
      "vmovups zmmword PTR [r12 + 64], zmm19\t\n"
      "add r12, r13\t\n"
      "vmovups zmmword PTR [r12 + 0], zmm20\t\n"
      "vmovups zmmword PTR [r12 + 64], zmm21\t\n"
      "add r12, r13\t\n"
      "vmovups zmmword PTR [r12 + 0], zmm22\t\n"
      "vmovups zmmword PTR [r12 + 64], zmm23\t\n"

      // next outer iteration
      "add rcx, 128\t\n"
      "mov r12, rcx\t\n"
      "mov r9, rax\t\n"
      "inc rbx\t\n"
      "cmp rbx, rdi\t\n"
      "jl loop_outter%=\t\n"
      :
      : [gp] "rm"(gp)
      : "r8",
        "r9",
        "r10",
        "r11",
        "r13",
        "r14",
        "rax",
        "rcx",
        "rsi",
        "rdi",
        "rbx",
        "r12",
        "r15",
        "memory");
}
void NOINLINE
gemmkernel_13x2_Avx512_fp16_fA0fB0fC0(GemmParamsFP16* gp) {
  asm volatile(
#if !defined(__clang__)
      "mov r14, %[gp]\t\n"
#else
      "mov %[gp], %%r14\t\n"
      ".intel_syntax noprefix\t\n"
#endif

      // Copy parameters
      // k
      "mov r8, [r14 + 0]\t\n"
      "dec r8\t\n"
      // A
      "mov r9, [r14 + 8]\t\n"
      // B
      "mov r10, [r14 + 16]\t\n"
      // beta
      "lea r15, [r14 + 24]\t\n"
      // C
      "mov r12, [r14 + 32]\t\n"
      // ldc
      "mov r13, [r14 + 40]\t\n"
      // b_block_cols
      "mov rdi, [r14 + 48]\t\n"
      // b_block_size
      "mov rsi, [r14 + 56]\t\n"

      // Make copies of A and C
      "mov rax, r9\t\n"
      "mov rcx, r12\t\n"

      "mov rbx, 0\t\n"
      "loop_outter%=:\t\n"
      "mov r14, r8\t\n"
      "vbroadcastss zmm31,DWORD PTR [r15]\t\n"
      "vcvtph2ps zmm27,YMMWORD PTR [r10 + 0]\t\n"
      "vcvtph2ps zmm28,YMMWORD PTR [r10 + 32]\t\n"
      "vxorps xmm0, xmm0, xmm0\t\n"
      "vcomiss xmm31, xmm0\t\n"
      "jz zero_regs%=\t\n"

      // Setup values with beta multiplication
      "vmulps zmm0, zmm31, [r12 + 0]\t\n"
      "vmulps zmm1, zmm31, [r12 + 64]\t\n"
      "add r12, r13\t\n"
      "vmulps zmm2, zmm31, [r12 + 0]\t\n"
      "vmulps zmm3, zmm31, [r12 + 64]\t\n"
      "add r12, r13\t\n"
      "vmulps zmm4, zmm31, [r12 + 0]\t\n"
      "vmulps zmm5, zmm31, [r12 + 64]\t\n"
      "add r12, r13\t\n"
      "vmulps zmm6, zmm31, [r12 + 0]\t\n"
      "vmulps zmm7, zmm31, [r12 + 64]\t\n"
      "add r12, r13\t\n"
      "vmulps zmm8, zmm31, [r12 + 0]\t\n"
      "vmulps zmm9, zmm31, [r12 + 64]\t\n"
      "add r12, r13\t\n"
      "vmulps zmm10, zmm31, [r12 + 0]\t\n"
      "vmulps zmm11, zmm31, [r12 + 64]\t\n"
      "add r12, r13\t\n"
      "vmulps zmm12, zmm31, [r12 + 0]\t\n"
      "vmulps zmm13, zmm31, [r12 + 64]\t\n"
      "add r12, r13\t\n"
      "vmulps zmm14, zmm31, [r12 + 0]\t\n"
      "vmulps zmm15, zmm31, [r12 + 64]\t\n"
      "add r12, r13\t\n"
      "vmulps zmm16, zmm31, [r12 + 0]\t\n"
      "vmulps zmm17, zmm31, [r12 + 64]\t\n"
      "add r12, r13\t\n"
      "vmulps zmm18, zmm31, [r12 + 0]\t\n"
      "vmulps zmm19, zmm31, [r12 + 64]\t\n"
      "add r12, r13\t\n"
      "vmulps zmm20, zmm31, [r12 + 0]\t\n"
      "vmulps zmm21, zmm31, [r12 + 64]\t\n"
      "add r12, r13\t\n"
      "vmulps zmm22, zmm31, [r12 + 0]\t\n"
      "vmulps zmm23, zmm31, [r12 + 64]\t\n"
      "add r12, r13\t\n"
      "vmulps zmm24, zmm31, [r12 + 0]\t\n"
      "vmulps zmm25, zmm31, [r12 + 64]\t\n"
      "test r14,r14\t\n"
      "jz skip_preload%=\t\n"
      "vcvtph2ps zmm31,YMMWORD PTR [r10 + 64]\t\n"
      "skip_preload%=:\t\n"
      "vbroadcastss zmm26,DWORD PTR [r9+0]\t\n"
      "vfmadd231ps zmm0,zmm27,zmm26\t\n"
      "vfmadd231ps zmm1,zmm28,zmm26\t\n"
      "vbroadcastss zmm26,DWORD PTR [r9+4]\t\n"
      "vfmadd231ps zmm2,zmm27,zmm26\t\n"
      "vfmadd231ps zmm3,zmm28,zmm26\t\n"
      "vbroadcastss zmm26,DWORD PTR [r9+8]\t\n"
      "vfmadd231ps zmm4,zmm27,zmm26\t\n"
      "vfmadd231ps zmm5,zmm28,zmm26\t\n"
      "vbroadcastss zmm26,DWORD PTR [r9+12]\t\n"
      "vfmadd231ps zmm6,zmm27,zmm26\t\n"
      "vfmadd231ps zmm7,zmm28,zmm26\t\n"
      "vbroadcastss zmm26,DWORD PTR [r9+16]\t\n"
      "vfmadd231ps zmm8,zmm27,zmm26\t\n"
      "vfmadd231ps zmm9,zmm28,zmm26\t\n"
      "vbroadcastss zmm26,DWORD PTR [r9+20]\t\n"
      "vfmadd231ps zmm10,zmm27,zmm26\t\n"
      "vfmadd231ps zmm11,zmm28,zmm26\t\n"
      "vbroadcastss zmm26,DWORD PTR [r9+24]\t\n"
      "vfmadd231ps zmm12,zmm27,zmm26\t\n"
      "vfmadd231ps zmm13,zmm28,zmm26\t\n"
      "vbroadcastss zmm26,DWORD PTR [r9+28]\t\n"
      "vfmadd231ps zmm14,zmm27,zmm26\t\n"
      "vfmadd231ps zmm15,zmm28,zmm26\t\n"
      "vbroadcastss zmm26,DWORD PTR [r9+32]\t\n"
      "vfmadd231ps zmm16,zmm27,zmm26\t\n"
      "vfmadd231ps zmm17,zmm28,zmm26\t\n"
      "vbroadcastss zmm26,DWORD PTR [r9+36]\t\n"
      "vfmadd231ps zmm18,zmm27,zmm26\t\n"
      "vfmadd231ps zmm19,zmm28,zmm26\t\n"
      "vbroadcastss zmm26,DWORD PTR [r9+40]\t\n"
      "vfmadd231ps zmm20,zmm27,zmm26\t\n"
      "vfmadd231ps zmm21,zmm28,zmm26\t\n"
      "vbroadcastss zmm26,DWORD PTR [r9+44]\t\n"
      "vfmadd231ps zmm22,zmm27,zmm26\t\n"
      "vfmadd231ps zmm23,zmm28,zmm26\t\n"
      "vbroadcastss zmm26,DWORD PTR [r9+48]\t\n"
      "vfmadd231ps zmm24,zmm27,zmm26\t\n"
      "vfmadd231ps zmm25,zmm28,zmm26\t\n"
      "mov r12, rcx\t\n"
      "test r14,r14\t\n"
      "jnz next_inner%=\t\n"
      "add r10,64\t\n"
      "jmp dump_C%=\t\n"

      "zero_regs%=:\t\n"

      "test r14,r14\t\n"
      "jz skip_preload_b_zero%=\t\n"
      "vcvtph2ps zmm31,YMMWORD PTR [r10 + 64]\t\n"
      "skip_preload_b_zero%=:\t\n"
      "vbroadcastss zmm26,DWORD PTR [r9+0]\t\n"
      "vmulps zmm0,zmm27,zmm26\t\n"
      "vmulps zmm1,zmm28,zmm26\t\n"
      "add r12, r13\t\n"
      "vbroadcastss zmm26,DWORD PTR [r9+4]\t\n"
      "vmulps zmm2,zmm27,zmm26\t\n"
      "vmulps zmm3,zmm28,zmm26\t\n"
      "add r12, r13\t\n"
      "vbroadcastss zmm26,DWORD PTR [r9+8]\t\n"
      "vmulps zmm4,zmm27,zmm26\t\n"
      "vmulps zmm5,zmm28,zmm26\t\n"
      "add r12, r13\t\n"
      "vbroadcastss zmm26,DWORD PTR [r9+12]\t\n"
      "vmulps zmm6,zmm27,zmm26\t\n"
      "vmulps zmm7,zmm28,zmm26\t\n"
      "add r12, r13\t\n"
      "vbroadcastss zmm26,DWORD PTR [r9+16]\t\n"
      "vmulps zmm8,zmm27,zmm26\t\n"
      "vmulps zmm9,zmm28,zmm26\t\n"
      "add r12, r13\t\n"
      "vbroadcastss zmm26,DWORD PTR [r9+20]\t\n"
      "vmulps zmm10,zmm27,zmm26\t\n"
      "vmulps zmm11,zmm28,zmm26\t\n"
      "add r12, r13\t\n"
      "vbroadcastss zmm26,DWORD PTR [r9+24]\t\n"
      "vmulps zmm12,zmm27,zmm26\t\n"
      "vmulps zmm13,zmm28,zmm26\t\n"
      "add r12, r13\t\n"
      "vbroadcastss zmm26,DWORD PTR [r9+28]\t\n"
      "vmulps zmm14,zmm27,zmm26\t\n"
      "vmulps zmm15,zmm28,zmm26\t\n"
      "add r12, r13\t\n"
      "vbroadcastss zmm26,DWORD PTR [r9+32]\t\n"
      "vmulps zmm16,zmm27,zmm26\t\n"
      "vmulps zmm17,zmm28,zmm26\t\n"
      "add r12, r13\t\n"
      "vbroadcastss zmm26,DWORD PTR [r9+36]\t\n"
      "vmulps zmm18,zmm27,zmm26\t\n"
      "vmulps zmm19,zmm28,zmm26\t\n"
      "add r12, r13\t\n"
      "vbroadcastss zmm26,DWORD PTR [r9+40]\t\n"
      "vmulps zmm20,zmm27,zmm26\t\n"
      "vmulps zmm21,zmm28,zmm26\t\n"
      "add r12, r13\t\n"
      "vbroadcastss zmm26,DWORD PTR [r9+44]\t\n"
      "vmulps zmm22,zmm27,zmm26\t\n"
      "vmulps zmm23,zmm28,zmm26\t\n"
      "add r12, r13\t\n"
      "vbroadcastss zmm26,DWORD PTR [r9+48]\t\n"
      "vmulps zmm24,zmm27,zmm26\t\n"
      "vmulps zmm25,zmm28,zmm26\t\n"
      "mov r12, rcx\t\n"
      "test r14,r14\t\n"
      "jnz next_inner%=\t\n"
      "add r10,64\t\n"
      "jmp dump_C%=\t\n"

      "loop_inner%=:\t\n"

      "vmovaps zmm27,zmm31\t\n"
      "vcvtph2ps zmm28,YMMWORD PTR [r10 + 32]\t\n"
      "vcvtph2ps zmm31,YMMWORD PTR [r10 + 64]\t\n"
      "vbroadcastss zmm26,DWORD PTR [r9+0]\t\n"
      "vfmadd231ps zmm0,zmm27,zmm26\t\n"
      "vfmadd231ps zmm1,zmm28,zmm26\t\n"
      "vbroadcastss zmm26,DWORD PTR [r9+4]\t\n"
      "vfmadd231ps zmm2,zmm27,zmm26\t\n"
      "vfmadd231ps zmm3,zmm28,zmm26\t\n"
      "vbroadcastss zmm26,DWORD PTR [r9+8]\t\n"
      "vfmadd231ps zmm4,zmm27,zmm26\t\n"
      "vfmadd231ps zmm5,zmm28,zmm26\t\n"
      "vbroadcastss zmm26,DWORD PTR [r9+12]\t\n"
      "vfmadd231ps zmm6,zmm27,zmm26\t\n"
      "vfmadd231ps zmm7,zmm28,zmm26\t\n"
      "vbroadcastss zmm26,DWORD PTR [r9+16]\t\n"
      "vfmadd231ps zmm8,zmm27,zmm26\t\n"
      "vfmadd231ps zmm9,zmm28,zmm26\t\n"
      "vbroadcastss zmm26,DWORD PTR [r9+20]\t\n"
      "vfmadd231ps zmm10,zmm27,zmm26\t\n"
      "vfmadd231ps zmm11,zmm28,zmm26\t\n"
      "vbroadcastss zmm26,DWORD PTR [r9+24]\t\n"
      "vfmadd231ps zmm12,zmm27,zmm26\t\n"
      "vfmadd231ps zmm13,zmm28,zmm26\t\n"
      "vbroadcastss zmm26,DWORD PTR [r9+28]\t\n"
      "vfmadd231ps zmm14,zmm27,zmm26\t\n"
      "vfmadd231ps zmm15,zmm28,zmm26\t\n"
      "vbroadcastss zmm26,DWORD PTR [r9+32]\t\n"
      "vfmadd231ps zmm16,zmm27,zmm26\t\n"
      "vfmadd231ps zmm17,zmm28,zmm26\t\n"
      "vbroadcastss zmm26,DWORD PTR [r9+36]\t\n"
      "vfmadd231ps zmm18,zmm27,zmm26\t\n"
      "vfmadd231ps zmm19,zmm28,zmm26\t\n"
      "vbroadcastss zmm26,DWORD PTR [r9+40]\t\n"
      "vfmadd231ps zmm20,zmm27,zmm26\t\n"
      "vfmadd231ps zmm21,zmm28,zmm26\t\n"
      "vbroadcastss zmm26,DWORD PTR [r9+44]\t\n"
      "vfmadd231ps zmm22,zmm27,zmm26\t\n"
      "vfmadd231ps zmm23,zmm28,zmm26\t\n"
      "vbroadcastss zmm26,DWORD PTR [r9+48]\t\n"
      "vfmadd231ps zmm24,zmm27,zmm26\t\n"
      "vfmadd231ps zmm25,zmm28,zmm26\t\n"

      "next_inner%=:\t\n"
      "add r9,52\t\n"
      "add r10,64\t\n"
      "dec r14\t\n"
      "jnz loop_inner%=\t\n"

      "vmovaps zmm27,zmm31\t\n"
      "vcvtph2ps zmm28,YMMWORD PTR [r10 + 32]\t\n"
      "vbroadcastss zmm26,DWORD PTR [r9+0]\t\n"
      "vfmadd231ps zmm0,zmm27,zmm26\t\n"
      "vfmadd231ps zmm1,zmm28,zmm26\t\n"
      "vbroadcastss zmm26,DWORD PTR [r9+4]\t\n"
      "vfmadd231ps zmm2,zmm27,zmm26\t\n"
      "vfmadd231ps zmm3,zmm28,zmm26\t\n"
      "vbroadcastss zmm26,DWORD PTR [r9+8]\t\n"
      "vfmadd231ps zmm4,zmm27,zmm26\t\n"
      "vfmadd231ps zmm5,zmm28,zmm26\t\n"
      "vbroadcastss zmm26,DWORD PTR [r9+12]\t\n"
      "vfmadd231ps zmm6,zmm27,zmm26\t\n"
      "vfmadd231ps zmm7,zmm28,zmm26\t\n"
      "vbroadcastss zmm26,DWORD PTR [r9+16]\t\n"
      "vfmadd231ps zmm8,zmm27,zmm26\t\n"
      "vfmadd231ps zmm9,zmm28,zmm26\t\n"
      "vbroadcastss zmm26,DWORD PTR [r9+20]\t\n"
      "vfmadd231ps zmm10,zmm27,zmm26\t\n"
      "vfmadd231ps zmm11,zmm28,zmm26\t\n"
      "vbroadcastss zmm26,DWORD PTR [r9+24]\t\n"
      "vfmadd231ps zmm12,zmm27,zmm26\t\n"
      "vfmadd231ps zmm13,zmm28,zmm26\t\n"
      "vbroadcastss zmm26,DWORD PTR [r9+28]\t\n"
      "vfmadd231ps zmm14,zmm27,zmm26\t\n"
      "vfmadd231ps zmm15,zmm28,zmm26\t\n"
      "vbroadcastss zmm26,DWORD PTR [r9+32]\t\n"
      "vfmadd231ps zmm16,zmm27,zmm26\t\n"
      "vfmadd231ps zmm17,zmm28,zmm26\t\n"
      "vbroadcastss zmm26,DWORD PTR [r9+36]\t\n"
      "vfmadd231ps zmm18,zmm27,zmm26\t\n"
      "vfmadd231ps zmm19,zmm28,zmm26\t\n"
      "vbroadcastss zmm26,DWORD PTR [r9+40]\t\n"
      "vfmadd231ps zmm20,zmm27,zmm26\t\n"
      "vfmadd231ps zmm21,zmm28,zmm26\t\n"
      "vbroadcastss zmm26,DWORD PTR [r9+44]\t\n"
      "vfmadd231ps zmm22,zmm27,zmm26\t\n"
      "vfmadd231ps zmm23,zmm28,zmm26\t\n"
      "vbroadcastss zmm26,DWORD PTR [r9+48]\t\n"
      "vfmadd231ps zmm24,zmm27,zmm26\t\n"
      "vfmadd231ps zmm25,zmm28,zmm26\t\n"
      "add r9,52\t\n"
      "add r10,64\t\n"
      // Dump C
      "dump_C%=:\t\n"
      "vmovups zmmword PTR [r12 + 0], zmm0\t\n"
      "vmovups zmmword PTR [r12 + 64], zmm1\t\n"
      "add r12, r13\t\n"
      "vmovups zmmword PTR [r12 + 0], zmm2\t\n"
      "vmovups zmmword PTR [r12 + 64], zmm3\t\n"
      "add r12, r13\t\n"
      "vmovups zmmword PTR [r12 + 0], zmm4\t\n"
      "vmovups zmmword PTR [r12 + 64], zmm5\t\n"
      "add r12, r13\t\n"
      "vmovups zmmword PTR [r12 + 0], zmm6\t\n"
      "vmovups zmmword PTR [r12 + 64], zmm7\t\n"
      "add r12, r13\t\n"
      "vmovups zmmword PTR [r12 + 0], zmm8\t\n"
      "vmovups zmmword PTR [r12 + 64], zmm9\t\n"
      "add r12, r13\t\n"
      "vmovups zmmword PTR [r12 + 0], zmm10\t\n"
      "vmovups zmmword PTR [r12 + 64], zmm11\t\n"
      "add r12, r13\t\n"
      "vmovups zmmword PTR [r12 + 0], zmm12\t\n"
      "vmovups zmmword PTR [r12 + 64], zmm13\t\n"
      "add r12, r13\t\n"
      "vmovups zmmword PTR [r12 + 0], zmm14\t\n"
      "vmovups zmmword PTR [r12 + 64], zmm15\t\n"
      "add r12, r13\t\n"
      "vmovups zmmword PTR [r12 + 0], zmm16\t\n"
      "vmovups zmmword PTR [r12 + 64], zmm17\t\n"
      "add r12, r13\t\n"
      "vmovups zmmword PTR [r12 + 0], zmm18\t\n"
      "vmovups zmmword PTR [r12 + 64], zmm19\t\n"
      "add r12, r13\t\n"
      "vmovups zmmword PTR [r12 + 0], zmm20\t\n"
      "vmovups zmmword PTR [r12 + 64], zmm21\t\n"
      "add r12, r13\t\n"
      "vmovups zmmword PTR [r12 + 0], zmm22\t\n"
      "vmovups zmmword PTR [r12 + 64], zmm23\t\n"
      "add r12, r13\t\n"
      "vmovups zmmword PTR [r12 + 0], zmm24\t\n"
      "vmovups zmmword PTR [r12 + 64], zmm25\t\n"

      // next outer iteration
      "add rcx, 128\t\n"
      "mov r12, rcx\t\n"
      "mov r9, rax\t\n"
      "inc rbx\t\n"
      "cmp rbx, rdi\t\n"
      "jl loop_outter%=\t\n"
      :
      : [gp] "rm"(gp)
      : "r8",
        "r9",
        "r10",
        "r11",
        "r13",
        "r14",
        "rax",
        "rcx",
        "rsi",
        "rdi",
        "rbx",
        "r12",
        "r15",
        "memory");
}
void NOINLINE
gemmkernel_14x2_Avx512_fp16_fA0fB0fC0(GemmParamsFP16* gp) {
  asm volatile(
#if !defined(__clang__)
      "mov r14, %[gp]\t\n"
#else
      "mov %[gp], %%r14\t\n"
      ".intel_syntax noprefix\t\n"
#endif

      // Copy parameters
      // k
      "mov r8, [r14 + 0]\t\n"
      "dec r8\t\n"
      // A
      "mov r9, [r14 + 8]\t\n"
      // B
      "mov r10, [r14 + 16]\t\n"
      // beta
      "lea r15, [r14 + 24]\t\n"
      // C
      "mov r12, [r14 + 32]\t\n"
      // ldc
      "mov r13, [r14 + 40]\t\n"
      // b_block_cols
      "mov rdi, [r14 + 48]\t\n"
      // b_block_size
      "mov rsi, [r14 + 56]\t\n"

      // Make copies of A and C
      "mov rax, r9\t\n"
      "mov rcx, r12\t\n"

      "mov rbx, 0\t\n"
      "loop_outter%=:\t\n"
      "mov r14, r8\t\n"
      "vbroadcastss zmm31,DWORD PTR [r15]\t\n"
      "vcvtph2ps zmm29,YMMWORD PTR [r10 + 0]\t\n"
      "vcvtph2ps zmm30,YMMWORD PTR [r10 + 32]\t\n"
      "vxorps xmm0, xmm0, xmm0\t\n"
      "vcomiss xmm31, xmm0\t\n"
      "jz zero_regs%=\t\n"

      // Setup values with beta multiplication
      "vmulps zmm0, zmm31, [r12 + 0]\t\n"
      "vmulps zmm1, zmm31, [r12 + 64]\t\n"
      "add r12, r13\t\n"
      "vmulps zmm2, zmm31, [r12 + 0]\t\n"
      "vmulps zmm3, zmm31, [r12 + 64]\t\n"
      "add r12, r13\t\n"
      "vmulps zmm4, zmm31, [r12 + 0]\t\n"
      "vmulps zmm5, zmm31, [r12 + 64]\t\n"
      "add r12, r13\t\n"
      "vmulps zmm6, zmm31, [r12 + 0]\t\n"
      "vmulps zmm7, zmm31, [r12 + 64]\t\n"
      "add r12, r13\t\n"
      "vmulps zmm8, zmm31, [r12 + 0]\t\n"
      "vmulps zmm9, zmm31, [r12 + 64]\t\n"
      "add r12, r13\t\n"
      "vmulps zmm10, zmm31, [r12 + 0]\t\n"
      "vmulps zmm11, zmm31, [r12 + 64]\t\n"
      "add r12, r13\t\n"
      "vmulps zmm12, zmm31, [r12 + 0]\t\n"
      "vmulps zmm13, zmm31, [r12 + 64]\t\n"
      "add r12, r13\t\n"
      "vmulps zmm14, zmm31, [r12 + 0]\t\n"
      "vmulps zmm15, zmm31, [r12 + 64]\t\n"
      "add r12, r13\t\n"
      "vmulps zmm16, zmm31, [r12 + 0]\t\n"
      "vmulps zmm17, zmm31, [r12 + 64]\t\n"
      "add r12, r13\t\n"
      "vmulps zmm18, zmm31, [r12 + 0]\t\n"
      "vmulps zmm19, zmm31, [r12 + 64]\t\n"
      "add r12, r13\t\n"
      "vmulps zmm20, zmm31, [r12 + 0]\t\n"
      "vmulps zmm21, zmm31, [r12 + 64]\t\n"
      "add r12, r13\t\n"
      "vmulps zmm22, zmm31, [r12 + 0]\t\n"
      "vmulps zmm23, zmm31, [r12 + 64]\t\n"
      "add r12, r13\t\n"
      "vmulps zmm24, zmm31, [r12 + 0]\t\n"
      "vmulps zmm25, zmm31, [r12 + 64]\t\n"
      "add r12, r13\t\n"
      "vmulps zmm26, zmm31, [r12 + 0]\t\n"
      "vmulps zmm27, zmm31, [r12 + 64]\t\n"
      "test r14,r14\t\n"
      "jz skip_preload%=\t\n"
      "vcvtph2ps zmm31,YMMWORD PTR [r10 + 64]\t\n"
      "skip_preload%=:\t\n"
      "vbroadcastss zmm28,DWORD PTR [r9+0]\t\n"
      "vfmadd231ps zmm0,zmm29,zmm28\t\n"
      "vfmadd231ps zmm1,zmm30,zmm28\t\n"
      "vbroadcastss zmm28,DWORD PTR [r9+4]\t\n"
      "vfmadd231ps zmm2,zmm29,zmm28\t\n"
      "vfmadd231ps zmm3,zmm30,zmm28\t\n"
      "vbroadcastss zmm28,DWORD PTR [r9+8]\t\n"
      "vfmadd231ps zmm4,zmm29,zmm28\t\n"
      "vfmadd231ps zmm5,zmm30,zmm28\t\n"
      "vbroadcastss zmm28,DWORD PTR [r9+12]\t\n"
      "vfmadd231ps zmm6,zmm29,zmm28\t\n"
      "vfmadd231ps zmm7,zmm30,zmm28\t\n"
      "vbroadcastss zmm28,DWORD PTR [r9+16]\t\n"
      "vfmadd231ps zmm8,zmm29,zmm28\t\n"
      "vfmadd231ps zmm9,zmm30,zmm28\t\n"
      "vbroadcastss zmm28,DWORD PTR [r9+20]\t\n"
      "vfmadd231ps zmm10,zmm29,zmm28\t\n"
      "vfmadd231ps zmm11,zmm30,zmm28\t\n"
      "vbroadcastss zmm28,DWORD PTR [r9+24]\t\n"
      "vfmadd231ps zmm12,zmm29,zmm28\t\n"
      "vfmadd231ps zmm13,zmm30,zmm28\t\n"
      "vbroadcastss zmm28,DWORD PTR [r9+28]\t\n"
      "vfmadd231ps zmm14,zmm29,zmm28\t\n"
      "vfmadd231ps zmm15,zmm30,zmm28\t\n"
      "vbroadcastss zmm28,DWORD PTR [r9+32]\t\n"
      "vfmadd231ps zmm16,zmm29,zmm28\t\n"
      "vfmadd231ps zmm17,zmm30,zmm28\t\n"
      "vbroadcastss zmm28,DWORD PTR [r9+36]\t\n"
      "vfmadd231ps zmm18,zmm29,zmm28\t\n"
      "vfmadd231ps zmm19,zmm30,zmm28\t\n"
      "vbroadcastss zmm28,DWORD PTR [r9+40]\t\n"
      "vfmadd231ps zmm20,zmm29,zmm28\t\n"
      "vfmadd231ps zmm21,zmm30,zmm28\t\n"
      "vbroadcastss zmm28,DWORD PTR [r9+44]\t\n"
      "vfmadd231ps zmm22,zmm29,zmm28\t\n"
      "vfmadd231ps zmm23,zmm30,zmm28\t\n"
      "vbroadcastss zmm28,DWORD PTR [r9+48]\t\n"
      "vfmadd231ps zmm24,zmm29,zmm28\t\n"
      "vfmadd231ps zmm25,zmm30,zmm28\t\n"
      "vbroadcastss zmm28,DWORD PTR [r9+52]\t\n"
      "vfmadd231ps zmm26,zmm29,zmm28\t\n"
      "vfmadd231ps zmm27,zmm30,zmm28\t\n"
      "mov r12, rcx\t\n"
      "test r14,r14\t\n"
      "jnz next_inner%=\t\n"
      "add r10,64\t\n"
      "jmp dump_C%=\t\n"

      "zero_regs%=:\t\n"

      "test r14,r14\t\n"
      "jz skip_preload_b_zero%=\t\n"
      "vcvtph2ps zmm31,YMMWORD PTR [r10 + 64]\t\n"
      "skip_preload_b_zero%=:\t\n"
      "vbroadcastss zmm28,DWORD PTR [r9+0]\t\n"
      "vmulps zmm0,zmm29,zmm28\t\n"
      "vmulps zmm1,zmm30,zmm28\t\n"
      "add r12, r13\t\n"
      "vbroadcastss zmm28,DWORD PTR [r9+4]\t\n"
      "vmulps zmm2,zmm29,zmm28\t\n"
      "vmulps zmm3,zmm30,zmm28\t\n"
      "add r12, r13\t\n"
      "vbroadcastss zmm28,DWORD PTR [r9+8]\t\n"
      "vmulps zmm4,zmm29,zmm28\t\n"
      "vmulps zmm5,zmm30,zmm28\t\n"
      "add r12, r13\t\n"
      "vbroadcastss zmm28,DWORD PTR [r9+12]\t\n"
      "vmulps zmm6,zmm29,zmm28\t\n"
      "vmulps zmm7,zmm30,zmm28\t\n"
      "add r12, r13\t\n"
      "vbroadcastss zmm28,DWORD PTR [r9+16]\t\n"
      "vmulps zmm8,zmm29,zmm28\t\n"
      "vmulps zmm9,zmm30,zmm28\t\n"
      "add r12, r13\t\n"
      "vbroadcastss zmm28,DWORD PTR [r9+20]\t\n"
      "vmulps zmm10,zmm29,zmm28\t\n"
      "vmulps zmm11,zmm30,zmm28\t\n"
      "add r12, r13\t\n"
      "vbroadcastss zmm28,DWORD PTR [r9+24]\t\n"
      "vmulps zmm12,zmm29,zmm28\t\n"
      "vmulps zmm13,zmm30,zmm28\t\n"
      "add r12, r13\t\n"
      "vbroadcastss zmm28,DWORD PTR [r9+28]\t\n"
      "vmulps zmm14,zmm29,zmm28\t\n"
      "vmulps zmm15,zmm30,zmm28\t\n"
      "add r12, r13\t\n"
      "vbroadcastss zmm28,DWORD PTR [r9+32]\t\n"
      "vmulps zmm16,zmm29,zmm28\t\n"
      "vmulps zmm17,zmm30,zmm28\t\n"
      "add r12, r13\t\n"
      "vbroadcastss zmm28,DWORD PTR [r9+36]\t\n"
      "vmulps zmm18,zmm29,zmm28\t\n"
      "vmulps zmm19,zmm30,zmm28\t\n"
      "add r12, r13\t\n"
      "vbroadcastss zmm28,DWORD PTR [r9+40]\t\n"
      "vmulps zmm20,zmm29,zmm28\t\n"
      "vmulps zmm21,zmm30,zmm28\t\n"
      "add r12, r13\t\n"
      "vbroadcastss zmm28,DWORD PTR [r9+44]\t\n"
      "vmulps zmm22,zmm29,zmm28\t\n"
      "vmulps zmm23,zmm30,zmm28\t\n"
      "add r12, r13\t\n"
      "vbroadcastss zmm28,DWORD PTR [r9+48]\t\n"
      "vmulps zmm24,zmm29,zmm28\t\n"
      "vmulps zmm25,zmm30,zmm28\t\n"
      "add r12, r13\t\n"
      "vbroadcastss zmm28,DWORD PTR [r9+52]\t\n"
      "vmulps zmm26,zmm29,zmm28\t\n"
      "vmulps zmm27,zmm30,zmm28\t\n"
      "mov r12, rcx\t\n"
      "test r14,r14\t\n"
      "jnz next_inner%=\t\n"
      "add r10,64\t\n"
      "jmp dump_C%=\t\n"

      "loop_inner%=:\t\n"

      "vmovaps zmm29,zmm31\t\n"
      "vcvtph2ps zmm30,YMMWORD PTR [r10 + 32]\t\n"
      "vcvtph2ps zmm31,YMMWORD PTR [r10 + 64]\t\n"
      "vbroadcastss zmm28,DWORD PTR [r9+0]\t\n"
      "vfmadd231ps zmm0,zmm29,zmm28\t\n"
      "vfmadd231ps zmm1,zmm30,zmm28\t\n"
      "vbroadcastss zmm28,DWORD PTR [r9+4]\t\n"
      "vfmadd231ps zmm2,zmm29,zmm28\t\n"
      "vfmadd231ps zmm3,zmm30,zmm28\t\n"
      "vbroadcastss zmm28,DWORD PTR [r9+8]\t\n"
      "vfmadd231ps zmm4,zmm29,zmm28\t\n"
      "vfmadd231ps zmm5,zmm30,zmm28\t\n"
      "vbroadcastss zmm28,DWORD PTR [r9+12]\t\n"
      "vfmadd231ps zmm6,zmm29,zmm28\t\n"
      "vfmadd231ps zmm7,zmm30,zmm28\t\n"
      "vbroadcastss zmm28,DWORD PTR [r9+16]\t\n"
      "vfmadd231ps zmm8,zmm29,zmm28\t\n"
      "vfmadd231ps zmm9,zmm30,zmm28\t\n"
      "vbroadcastss zmm28,DWORD PTR [r9+20]\t\n"
      "vfmadd231ps zmm10,zmm29,zmm28\t\n"
      "vfmadd231ps zmm11,zmm30,zmm28\t\n"
      "vbroadcastss zmm28,DWORD PTR [r9+24]\t\n"
      "vfmadd231ps zmm12,zmm29,zmm28\t\n"
      "vfmadd231ps zmm13,zmm30,zmm28\t\n"
      "vbroadcastss zmm28,DWORD PTR [r9+28]\t\n"
      "vfmadd231ps zmm14,zmm29,zmm28\t\n"
      "vfmadd231ps zmm15,zmm30,zmm28\t\n"
      "vbroadcastss zmm28,DWORD PTR [r9+32]\t\n"
      "vfmadd231ps zmm16,zmm29,zmm28\t\n"
      "vfmadd231ps zmm17,zmm30,zmm28\t\n"
      "vbroadcastss zmm28,DWORD PTR [r9+36]\t\n"
      "vfmadd231ps zmm18,zmm29,zmm28\t\n"
      "vfmadd231ps zmm19,zmm30,zmm28\t\n"
      "vbroadcastss zmm28,DWORD PTR [r9+40]\t\n"
      "vfmadd231ps zmm20,zmm29,zmm28\t\n"
      "vfmadd231ps zmm21,zmm30,zmm28\t\n"
      "vbroadcastss zmm28,DWORD PTR [r9+44]\t\n"
      "vfmadd231ps zmm22,zmm29,zmm28\t\n"
      "vfmadd231ps zmm23,zmm30,zmm28\t\n"
      "vbroadcastss zmm28,DWORD PTR [r9+48]\t\n"
      "vfmadd231ps zmm24,zmm29,zmm28\t\n"
      "vfmadd231ps zmm25,zmm30,zmm28\t\n"
      "vbroadcastss zmm28,DWORD PTR [r9+52]\t\n"
      "vfmadd231ps zmm26,zmm29,zmm28\t\n"
      "vfmadd231ps zmm27,zmm30,zmm28\t\n"

      "next_inner%=:\t\n"
      "add r9,56\t\n"
      "add r10,64\t\n"
      "dec r14\t\n"
      "jnz loop_inner%=\t\n"

      "vmovaps zmm29,zmm31\t\n"
      "vcvtph2ps zmm30,YMMWORD PTR [r10 + 32]\t\n"
      "vbroadcastss zmm28,DWORD PTR [r9+0]\t\n"
      "vfmadd231ps zmm0,zmm29,zmm28\t\n"
      "vfmadd231ps zmm1,zmm30,zmm28\t\n"
      "vbroadcastss zmm28,DWORD PTR [r9+4]\t\n"
      "vfmadd231ps zmm2,zmm29,zmm28\t\n"
      "vfmadd231ps zmm3,zmm30,zmm28\t\n"
      "vbroadcastss zmm28,DWORD PTR [r9+8]\t\n"
      "vfmadd231ps zmm4,zmm29,zmm28\t\n"
      "vfmadd231ps zmm5,zmm30,zmm28\t\n"
      "vbroadcastss zmm28,DWORD PTR [r9+12]\t\n"
      "vfmadd231ps zmm6,zmm29,zmm28\t\n"
      "vfmadd231ps zmm7,zmm30,zmm28\t\n"
      "vbroadcastss zmm28,DWORD PTR [r9+16]\t\n"
      "vfmadd231ps zmm8,zmm29,zmm28\t\n"
      "vfmadd231ps zmm9,zmm30,zmm28\t\n"
      "vbroadcastss zmm28,DWORD PTR [r9+20]\t\n"
      "vfmadd231ps zmm10,zmm29,zmm28\t\n"
      "vfmadd231ps zmm11,zmm30,zmm28\t\n"
      "vbroadcastss zmm28,DWORD PTR [r9+24]\t\n"
      "vfmadd231ps zmm12,zmm29,zmm28\t\n"
      "vfmadd231ps zmm13,zmm30,zmm28\t\n"
      "vbroadcastss zmm28,DWORD PTR [r9+28]\t\n"
      "vfmadd231ps zmm14,zmm29,zmm28\t\n"
      "vfmadd231ps zmm15,zmm30,zmm28\t\n"
      "vbroadcastss zmm28,DWORD PTR [r9+32]\t\n"
      "vfmadd231ps zmm16,zmm29,zmm28\t\n"
      "vfmadd231ps zmm17,zmm30,zmm28\t\n"
      "vbroadcastss zmm28,DWORD PTR [r9+36]\t\n"
      "vfmadd231ps zmm18,zmm29,zmm28\t\n"
      "vfmadd231ps zmm19,zmm30,zmm28\t\n"
      "vbroadcastss zmm28,DWORD PTR [r9+40]\t\n"
      "vfmadd231ps zmm20,zmm29,zmm28\t\n"
      "vfmadd231ps zmm21,zmm30,zmm28\t\n"
      "vbroadcastss zmm28,DWORD PTR [r9+44]\t\n"
      "vfmadd231ps zmm22,zmm29,zmm28\t\n"
      "vfmadd231ps zmm23,zmm30,zmm28\t\n"
      "vbroadcastss zmm28,DWORD PTR [r9+48]\t\n"
      "vfmadd231ps zmm24,zmm29,zmm28\t\n"
      "vfmadd231ps zmm25,zmm30,zmm28\t\n"
      "vbroadcastss zmm28,DWORD PTR [r9+52]\t\n"
      "vfmadd231ps zmm26,zmm29,zmm28\t\n"
      "vfmadd231ps zmm27,zmm30,zmm28\t\n"
      "add r9,56\t\n"
      "add r10,64\t\n"
      // Dump C
      "dump_C%=:\t\n"
      "vmovups zmmword PTR [r12 + 0], zmm0\t\n"
      "vmovups zmmword PTR [r12 + 64], zmm1\t\n"
      "add r12, r13\t\n"
      "vmovups zmmword PTR [r12 + 0], zmm2\t\n"
      "vmovups zmmword PTR [r12 + 64], zmm3\t\n"
      "add r12, r13\t\n"
      "vmovups zmmword PTR [r12 + 0], zmm4\t\n"
      "vmovups zmmword PTR [r12 + 64], zmm5\t\n"
      "add r12, r13\t\n"
      "vmovups zmmword PTR [r12 + 0], zmm6\t\n"
      "vmovups zmmword PTR [r12 + 64], zmm7\t\n"
      "add r12, r13\t\n"
      "vmovups zmmword PTR [r12 + 0], zmm8\t\n"
      "vmovups zmmword PTR [r12 + 64], zmm9\t\n"
      "add r12, r13\t\n"
      "vmovups zmmword PTR [r12 + 0], zmm10\t\n"
      "vmovups zmmword PTR [r12 + 64], zmm11\t\n"
      "add r12, r13\t\n"
      "vmovups zmmword PTR [r12 + 0], zmm12\t\n"
      "vmovups zmmword PTR [r12 + 64], zmm13\t\n"
      "add r12, r13\t\n"
      "vmovups zmmword PTR [r12 + 0], zmm14\t\n"
      "vmovups zmmword PTR [r12 + 64], zmm15\t\n"
      "add r12, r13\t\n"
      "vmovups zmmword PTR [r12 + 0], zmm16\t\n"
      "vmovups zmmword PTR [r12 + 64], zmm17\t\n"
      "add r12, r13\t\n"
      "vmovups zmmword PTR [r12 + 0], zmm18\t\n"
      "vmovups zmmword PTR [r12 + 64], zmm19\t\n"
      "add r12, r13\t\n"
      "vmovups zmmword PTR [r12 + 0], zmm20\t\n"
      "vmovups zmmword PTR [r12 + 64], zmm21\t\n"
      "add r12, r13\t\n"
      "vmovups zmmword PTR [r12 + 0], zmm22\t\n"
      "vmovups zmmword PTR [r12 + 64], zmm23\t\n"
      "add r12, r13\t\n"
      "vmovups zmmword PTR [r12 + 0], zmm24\t\n"
      "vmovups zmmword PTR [r12 + 64], zmm25\t\n"
      "add r12, r13\t\n"
      "vmovups zmmword PTR [r12 + 0], zmm26\t\n"
      "vmovups zmmword PTR [r12 + 64], zmm27\t\n"

      // next outer iteration
      "add rcx, 128\t\n"
      "mov r12, rcx\t\n"
      "mov r9, rax\t\n"
      "inc rbx\t\n"
      "cmp rbx, rdi\t\n"
      "jl loop_outter%=\t\n"
      :
      : [gp] "rm"(gp)
      : "r8",
        "r9",
        "r10",
        "r11",
        "r13",
        "r14",
        "rax",
        "rcx",
        "rsi",
        "rdi",
        "rbx",
        "r12",
        "r15",
        "memory");
}

} // namespace fbgemm
