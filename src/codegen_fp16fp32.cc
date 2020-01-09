/*
 * Copyright (c) Facebook, Inc. and its affiliates.
 * All rights reserved.
 * This source code is licensed under the BSD-style license found in the
 * LICENSE file in the root directory of this source tree.
 */

/* Arguments to code gen
 * --prefetch[=n] - use prefetch in code generation, n == prefetch len
 * --fp32         - generate code for FBGEMM32 (SGEMM)
 */

#include <assert.h>
#include <cpuid.h>
#include <stdlib.h>
#include <string.h>
#include <algorithm>
#include <array>
#include <fstream>
#include <functional>
#include <iostream>
#include <map>
#include <string>
#include <unordered_set>

using namespace std;

void addi(ofstream& of, string i, string asmstr = "", bool disable = false) {
  if (disable == false)
    of << "      " + i + "      //\"" + asmstr + "\\t\\n\"" + "\n";
}
#if 0
void addi(ofstream& of, string i, bool disable = false) {
  if (disable == false)
    of << "      \"" + i + "\\t\\n\"" + "\n";
}
#endif

struct ISA {
  enum class isaType { avx, avx2, avx512, avx512_256 };
  isaType iset;
  string name;
  vector<vector<unsigned>> shapes;
};

enum class DataType { Float32, Float16, BFloat16 };

std::array<std::pair<DataType, std::string>, 2> types_to_gen = {
    std::make_pair(DataType::Float32, "FP32"),
    std::make_pair(DataType::Float16, "FP16")};

constexpr int cache_line_size = 64;
constexpr int prefetch_len_default = cache_line_size * 12;

int parseArgumentInt(
    int argc,
    const char* argv[],
    const char* arg,
    int non_exist_val,
    int def_val) {
  int val = non_exist_val;
  int arg_len = strlen(arg);
  for (auto i = 1; i < argc; ++i) {
    const char* ptr = strstr(argv[i], arg);
    if (ptr) {
      val = (*(ptr + arg_len) == '=') ? atoi(ptr + arg_len + 1) : def_val;
      break;
    }
  }
  return val;
}

bool parseArgumentBool(
    int argc,
    const char* argv[],
    const char* arg,
    bool def_val) {
  int arg_len = strlen(arg);
  for (auto i = 1; i < argc; ++i) {
    const char* ptr = strstr(argv[i], arg);
    if (ptr) {
      return true;
    }
  }
  return def_val;
}

int main(int argc, const char* argv[]) {
  bool iaca = false;
  bool disable = false;
  unordered_set<string> enabledDataType;

  // Always generate FP16
  enabledDataType.insert("FP16");
  if (parseArgumentBool(argc, argv, "--fp32", false)) {
    enabledDataType.insert("FP32");
  }

  // Frefetch 8 cache lines ahead
  const int prefetch_a_len =
      parseArgumentInt(argc, argv, "--prefetch-a", 0, 128);

  // Frefetch 8 cache lines ahead
  const int prefetch_b_len =
      parseArgumentInt(argc, argv, "--prefetch-b", 0, 256);

  const int prefetch_c_len =
      parseArgumentInt(argc, argv, "--prefetch-c", 0, 1024);

  bool fixedA = true, fixedB = true, fixedC = true;

  int eax, ebx, ecx, edx;
  __cpuid(1 /* ecx = vendor string */, eax, ebx, ecx, edx);
  printf("FC16 is %s supported\n", ((ecx & bit_F16C) ? " " : "not"));

  static const string license =
      "/*\n"
      " * Copyright (c) Facebook, Inc. and its affiliates.\n"
      " * All rights reserved.\n"
      " * This source code is licensed under the BSD-style license found in the\n"
      " * LICENSE file in the root directory of this source tree.\n"
      " */\n";

  string comma = ",";

  enum class mult_type {
    fma,
    mul
  };

  vector<ISA> isa = {
      // {1, "AVX", {{4, 1, 0}, {4, 2, 0}, {4, 3, 0}, {3, 1, 0}, {3, 2, 0}, {3,
      // 3, 0}}},
      {ISA::isaType::avx2,
       "Avx2",
       {
           // 4x3 register layout
           // {1, 3, 0},
           // {2, 3, 0},
           // {3, 3, 0},
           // {4, 3, 0},

           // 6x2 register layout
           {1, 2, 0},
           {2, 2, 0},
           {3, 2, 0},
           {4, 2, 0},
           {5, 2, 0},
           {6, 2, 0},

           // 14x1 register layout
           // {1, 1, 0},
           // {2, 1, 0},
           // {3, 1, 0},
           // {4, 1, 0},
           // {5, 1, 0},
           // {6, 1, 0},
           // {7, 1, 0},
           // {8, 1, 0},
           // {9, 1, 0},
           // {10, 1, 0},
           // {11, 1, 0},
           // {12, 1, 0},
           // {13, 1, 0},
           // {14, 1, 0},
       }},
      {ISA::isaType::avx512,
       "Avx512",
       {
           // 14x2 register layout
           {1, 2, 0},
           {2, 2, 0},
           {3, 2, 0},
           {4, 2, 0},
           {5, 2, 0},
           {6, 2, 0},
           {7, 2, 0},
           {8, 2, 0},
           {9, 2, 0},
           {10, 2, 0},
           {11, 2, 0},
           {12, 2, 0},
           {13, 2, 0},
           {14, 2, 0},
       }},
      {ISA::isaType::avx512_256,
       "Avx512_256",
       {
           // 14x2 register layout
           // Implemented by AVX2
           //{1, 2, 0},
           //{2, 2, 0},
           //{3, 2, 0},
           //{4, 2, 0},
           //{5, 2, 0},
           //{6, 2, 0},
           {7, 2, 0},
           {8, 2, 0},
           {9, 2, 0},
           {10, 2, 0},
           {11, 2, 0},
           {12, 2, 0},
           {13, 2, 0},
           {14, 2, 0},
       }}};

<<<<<<< HEAD
  // open all files
  ofstream srcfile;
  srcfile.open("FbgemmFP16UKernelsAvx2.cc");
  srcfile
      << "/*\n"
         " * Copyright (c) Facebook, Inc. and its affiliates.\n"
         " * All rights reserved.\n"
         " * This source code is licensed under the BSD-style license found in the\n"
         " * LICENSE file in the root directory of this source tree.\n"
         " */\n";
  srcfile << "#include \"FbgemmFP16UKernelsAvx2.h\"\n";
  srcfile << "#include <immintrin.h>\n\n";
  srcfile << "namespace fbgemm {\n\n";
  if (iaca) {
    srcfile << "#include \"iacaMarks.h\"\n";
  }

  ofstream hdrfile;
  hdrfile.open("FbgemmFP16UKernelsAvx2.h");
  hdrfile
      << "/*\n"
         " * Copyright (c) Facebook, Inc. and its affiliates.\n"
         " * All rights reserved.\n"
         " * This source code is licensed under the BSD-style license found in the\n"
         " * LICENSE file in the root directory of this source tree.\n"
         " */\n";

  hdrfile << "#ifndef FBGEMM_UKERNELS\n";
  hdrfile << "#define FBGEMM_UKERNELS\n";
  hdrfile << "#include <cstdint>\n";
  hdrfile << "#include \"fbgemm/Types.h\"\n\n";
  hdrfile << "namespace fbgemm {\n\n";
  hdrfile << "using fp16 = float16;\n";
  hdrfile << "using fp32 = float;\n";
  hdrfile << "#ifdef _MSC_VER\n";
  hdrfile << " #define NOINLINE_ATTR __declspec(noinline)\n";
  hdrfile << "#else\n";
  hdrfile << " #define NOINLINE_ATTR __attribute__((noinline))\n";
  hdrfile << "#endif\n";
  hdrfile
      << "struct GemmParams {\n  uint64_t k;\n  float* A;\n  const fp16* B;\n"
         "  float* beta;\n  uint64_t accum;\n  float* C;\n  uint64_t ldc;\n"
         "  uint64_t b_block_cols;\n  uint64_t b_block_size;\n};\n";

  std::map<string, string> fptr_typedef;
  fptr_typedef["fp16"] = "";
  fptr_typedef["fp32"] = "";

  unsigned labelId = 0;
#if 1
  for (auto fixedA : {false})
    for (auto fixedB : {false})
      for (auto fixedC : {false})
#else
  for (auto fixedA : {true})
    for (auto fixedB : {true})
      for (auto fixedC : {true})
#endif
        for (auto s : isa) {
          vector<vector<unsigned>>& ukernel_shape = s.shapes;

          vector<string> funcname(ukernel_shape.size()),
              fheader(ukernel_shape.size());
          string fargs;

          for (auto fp16 : {true}) {
            string B_type = ((fp16) ? "fp16" : "fp32");
            string prefix = s.name + /*"_" + B_type */ +"_" + "fA" +
                to_string(fixedA) + "fB" + to_string(fixedB) + "fC" +
                to_string(fixedC);
            cout << "Generating code for " << s.name << " " << B_type << "\n";

            for (unsigned k = 0; k < ukernel_shape.size(); k++) {
              printf(
                  "shape: %d x %d * 32\n",
                  ukernel_shape[k][0],
                  ukernel_shape[k][1]);

              string p1 = "GemmParams* gp";

              funcname[k] = "gemmkernel_" + to_string(ukernel_shape[k][0]) +
                  "x" + to_string(ukernel_shape[k][1]) + "_";
              funcname[k] += prefix;

              fargs = "(" + p1 + ")";

#if 1
              fheader[k] =
                  "void NOINLINE_ATTR " + funcname[k] + fargs;
              srcfile << fheader[k] << " {\n";

              unsigned last_free_ymmreg = 0;
              // produce register block of C
              vector<vector<string>> vCtile(ukernel_shape[k][0]);
              for (auto r = 0; r < ukernel_shape[k][0]; r++)
                for (auto c = 0; c < ukernel_shape[k][1]; c++) {
                  vCtile[r].push_back("ymm" + to_string(last_free_ymmreg));
                  last_free_ymmreg++;
                }
              assert(last_free_ymmreg <= 14);

              string vAtmp = "ymm" + to_string(last_free_ymmreg++);
              // produce register block of B col
              vector<string> vBcol(ukernel_shape[k][1]);

              for (auto c = 0; c < ukernel_shape[k][1]; c++) {
                vBcol[c] = ("ymm" + to_string(last_free_ymmreg));
                last_free_ymmreg++;
              }

              assert(last_free_ymmreg <= 16);

              //srcfile << "  asm volatile(\n";

              //srcfile << "#if !defined(__clang__)"
                      //<< "\n";
              addi(srcfile, "char* r14 = (char*)gp;", "mov r14, %[gp]");
              //srcfile << "#else\n";
              //addi(srcfile, "mov %[gp], %%r14");
              //addi(srcfile, ".intel_syntax noprefix");
              //srcfile << "#endif\n";

              srcfile << "\n      // Copy parameters\n";
              srcfile << "      // k\n";            addi(srcfile, "uint64_t    r8  = *(uint64_t   *)((char*)r14 + 0 );", "mov r8, [r14 + 0]");
              srcfile << "      // A\n";            addi(srcfile, "float*      r9  = *(float*     *)((char*)r14 + 8 );", "mov r9, [r14 + 8]");
              srcfile << "      // B\n";            addi(srcfile, "const fp16* r10 = *(const fp16**)((char*)r14 + 16);", "mov r10, [r14 + 16]");
              srcfile << "      // beta\n";         addi(srcfile, "float*      r15 = *(float*     *)((char*)r14 + 24);", "mov r15, [r14 + 24]");
              srcfile << "      // accum\n";        addi(srcfile, "uint64_t    rdx = *(uint64_t   *)((char*)r14 + 32);", "mov rdx, [r14 + 32]");
              srcfile << "      // C\n";            addi(srcfile, "float*      r12 = *(float*     *)((char*)r14 + 40);", "mov r12, [r14 + 40]");
              srcfile << "      // ldc\n";          addi(srcfile, "uint64_t    r13 = *(uint64_t   *)((char*)r14 + 48);", "mov r13, [r14 + 48]");
              srcfile << "      // b_block_cols\n"; addi(srcfile, "uint64_t    rdi = *(uint64_t   *)((char*)r14 + 56);", "mov rdi, [r14 + 56]");
              srcfile << "      // b_block_size\n"; addi(srcfile, "uint64_t    rsi = *(uint64_t   *)((char*)r14 + 64);", "mov rsi, [r14 + 64]");
              srcfile << "      // Make copies of A and C\n";
              addi(srcfile, "float* rax = r9;", "mov rax, r9");
              addi(srcfile, "float* rcx = r12;", "mov rcx, r12");
              srcfile << "\n";

              addi(srcfile, "uint64_t rbx = 0;", "mov rbx, 0");

              string exitlabel = "L_exit%=";
              string label2 = "loop_outter%=";
              addi(srcfile, "for (; rbx < rdi; ++rbx) {", "inc rbx; cmp rbx, rdi; jl " + label2);
              addi(srcfile, "// ", label2 + ":");
              addi(srcfile, "  uint64_t r14_i = 0;", "mov r14, 0");

              // set all vCtile regs to zeros
              for (auto r = 0; r < vCtile.size(); r++) {
                for (auto c = 0; c < vCtile[r].size(); c++) {
                  addi(
                      srcfile,
                      "  __m256 " + vCtile[r][c] + " = _mm256_setzero_ps();",
                      "vxorps " + vCtile[r][c] + "," + vCtile[r][c] + "," +
                          vCtile[r][c]);
                }
              }

              // start marker
              //if (iaca) {
              //  addi(srcfile, "mov ebx, 111");
              //  addi(srcfile, ".byte 0x64, 0x67, 0x90");
              //}

              //srcfile << "\n";

              srcfile << "\n";
              string label = "loop_inner%=";
              addi(srcfile, "  for (; r14_i < r8; ++r14_i) {", "inc r14; cmp r14, r8; jl " + label);
              addi(srcfile, "  // " + label + ":");
              //srcfile << "\n";

              for (int c = 0; c < vCtile[0].size(); c++) {
                addi(
                  srcfile,
                  "    auto fp16mem" + to_string(16 * c) + " = _mm_load_si128((__m128i*)((char*)r10 + " + to_string(16 * c) + "));",
                  "vcvtph2ps " + vBcol[c] + ",XMMWORD PTR [r10 + " +
                  to_string(16 * c) + "]");
                addi(
                    srcfile,
                    "    auto " + vBcol[c] + " = _mm256_cvtph_ps(fp16mem" + to_string(16 * c) + ");",
                    "vcvtph2ps " + vBcol[c] + ",XMMWORD PTR [r10 + " +
                        to_string(16 * c) + "]");
              }

              for (int r = 0; r < vCtile.size(); r++) {
                //addi(
                //  srcfile,
                //  ((r == 0) ? "    auto " + vAtmp : "" + vAtmp) + " = _mm256_broadcastss_ps(r9 + " + to_string(4 * r) + ");",
                //  "vbroadcastss " + vAtmp + ",DWORD PTR [r9+" +
                //  to_string(4 * r) + "]");
                addi(
                    srcfile,
                    ((r == 0) ? "    auto " + vAtmp : "    " + vAtmp) + " = _mm256_broadcast_ss((float*)((char*)r9 + " + to_string(4 * r) + "));",
                    "vbroadcastss " + vAtmp + ",DWORD PTR [r9+" +
                        to_string(4 * r) + "]");
                for (int c = 0; c < vCtile[0].size(); c++) {
                  addi(
                      srcfile,
                      "    " + vCtile[r][c] + " = _mm256_fmadd_ps(" + vAtmp + ", " + vBcol[c] + ", " + vCtile[r][c] + ");",
                      "vfmadd231ps " + vCtile[r][c] + "," + vBcol[c] + "," +
                          vAtmp);
                }
              }

              addi(
                  srcfile,
                  "    r9 = (float*)((char*)r9 + " + to_string(4 * ukernel_shape[k][0]) + ");",
                  "add r9," + to_string(4 * ukernel_shape[k][0]),
                  fixedA); // move A ptr

              addi(
                  srcfile,
                  "    r10 = (fp16*)((char*)r10 + " + to_string(16 * ukernel_shape[k][1]) + ");",
                  "add r10," + to_string(16 * ukernel_shape[k][1]),
                  fixedA); // move A ptr

              addi(srcfile, "  }", "inc r14; cmp r14, r8; jl " + label2);
              // move to for loop
              //addi(srcfile, "inc r14");
              //addi(srcfile, "cmp r14, r8");
              //addi(srcfile, "jl " + label);

              //srcfile << "\n";

              //addi(srcfile, exitlabel + ":");

              // addi(srcfile, "add r10, rsi");
              srcfile << "\n";

              // end marker
              if (iaca) {
                addi(srcfile, "mov ebx, 222");
                addi(srcfile, ".byte 0x64, 0x67, 0x90");
              }

              //addi(srcfile, "cmp rdx, 1");
              addi(srcfile, "  if(rdx != 1) {", "cmp rdx, 1; je L_accum%=");

              srcfile << "          // Dump C\n";

              for (auto r = 0; r < vCtile.size(); r++) {
                for (auto c = 0; c < vCtile[r].size(); c++) {
                  addi(
                      srcfile,
                      "    _mm256_storeu_ps((float*)((char*)r12 + " + to_string(32 * c) + "), " + vCtile[r][c] + ");",
                      "vmovups YMMWORD PTR [r12 + " + to_string(32 * c) +
                          "], " + vCtile[r][c],
                      fixedC);
                }
                if (r != vCtile.size() - 1)
                  addi(srcfile, "    r12 = (float*)((char*)r12 + r13);", "add r12, r13", fixedC); // move C ptr
              }
              addi(srcfile, "  } else {", "jmp L_done%=");

              //srcfile << "\n";
              //addi(srcfile, "L_accum%=:");
              srcfile << "          // Dump C with accumulate\n";
=======
   // Labels
   const string label_outer = "loop_outter%=";
   const string label_next_inner = "next_inner%=";
   const string label_inner = "loop_inner%=";
   const string label_zero = "zero_regs%=";
   const string label_dump_C = "dump_C%=";

  for (auto& d_type : types_to_gen) {
    if (enabledDataType.count(d_type.second) == 0) {
      continue;
    }
    for (auto s : isa) {
      bool const isFp16 = d_type.first != DataType::Float32;
      string const B_type = [&]() {
        if (d_type.first == DataType::Float32)
          return "fp32";
        if (d_type.first == DataType::Float16)
          return "fp16";
      }();

      string isa_file_name = "Fbgemm" + d_type.second + "UKernels" + s.name;

      // open all files
      ofstream srcfile;
      srcfile.open(isa_file_name + ".cc");
      srcfile << license;
      srcfile << "#include \"./" + isa_file_name + ".h\"\n\n";
      srcfile << "namespace fbgemm {\n\n";
      if (iaca) {
        srcfile << "#include \"iacaMarks.h\"\n";
      }

      ofstream hdrfile;
      hdrfile.open(isa_file_name + ".h");
      hdrfile << license;

      hdrfile << "#pragma once\n";
      hdrfile << "#include <cstdint>\n";
      hdrfile << "#include \"fbgemm/Types.h\"\n\n";
      hdrfile << "#include \"fbgemm/FbgemmBuild.h\"\n\n";
      hdrfile << "#include \"./Fbgemm" + d_type.second + "Common.h\"\n\n";
      hdrfile << "namespace fbgemm {\n\n";

      unsigned labelId = 0;

      bool fixedA = false, fixedB = false, fixedC = false;

      vector<vector<unsigned>>& ukernel_shape = s.shapes;

      vector<string> funcname(ukernel_shape.size()),
          fheader(ukernel_shape.size());
      string fargs;

      string prefix = s.name + "_" + B_type + "_" + "fA" + to_string(fixedA) +
          "fB" + to_string(fixedB) + "fC" + to_string(fixedC);
      cout << "Generating code for " << s.name << " " << B_type << "\n";

      string vec_reg_prefix = s.iset == ISA::isaType::avx512 ? "zmm" : "ymm";
      int num_vec_regs = s.iset == ISA::isaType::avx2 ? 16 : 32;
      int vec_len_in_bytes = s.iset == ISA::isaType::avx512 ? 64 : 32;

      for (unsigned k = 0; k < ukernel_shape.size(); k++) {
        printf(
            "shape: %d x %d * 32\n", ukernel_shape[k][0], ukernel_shape[k][1]);

        const string A_stride = to_string(4 * ukernel_shape[k][0]);
        const string B_stride = to_string(
            (vec_len_in_bytes >> (int)isFp16) * ukernel_shape[k][1]);

        const string p1 = "GemmParams" + d_type.second + "* gp";

        funcname[k] = "gemmkernel_" + to_string(ukernel_shape[k][0]) + "x" +
            to_string(ukernel_shape[k][1]) + "_";
        funcname[k] += prefix;

        fargs = "(" + p1 + ")";

        fheader[k] = "void NOINLINE\n" + funcname[k] + fargs;
        srcfile << fheader[k] << " {\n";

        unsigned last_free_vecreg = 0;
        // produce register block of C
        vector<vector<string>> vCtile(ukernel_shape[k][0]);
        for (auto r = 0; r < ukernel_shape[k][0]; r++)
          for (auto c = 0; c < ukernel_shape[k][1]; c++) {
            vCtile[r].push_back(vec_reg_prefix + to_string(last_free_vecreg));
            last_free_vecreg++;
          }
        assert(last_free_vecreg <= num_vec_regs - 2);

        string vAtmp = vec_reg_prefix + to_string(last_free_vecreg++);
        // produce register block of B col
        vector<string> vBcol(ukernel_shape[k][1]);

        for (auto c = 0; c < ukernel_shape[k][1]; c++) {
          vBcol[c] = (vec_reg_prefix + to_string(last_free_vecreg));
          last_free_vecreg++;
        }
        assert(last_free_vecreg <= num_vec_regs);
        string r_spare = vec_reg_prefix +
            to_string(num_vec_regs - (s.iset == ISA::isaType::avx ? 2 : 1));

        auto const A_load_mult = [&](int r, mult_type m_type) {
          if (prefetch_a_len && ((4 * r) % cache_line_size == 0)) {
            addi(
                srcfile,
                "prefetcht0 [r9 + " + to_string(prefetch_a_len) + "]",
                fixedC);
          }
          string mul = m_type == mult_type::mul ? "vmulps" : "vfmadd231ps";
          addi(
              srcfile,
              "vbroadcastss " + vAtmp + ",DWORD PTR [r9+" + to_string(4 * r) +
                  "]");
          for (int c = 0; c < vCtile[0].size(); c++) {
            addi(
                srcfile,
                mul + " " + vCtile[r][c] + "," + vBcol[c] + "," + vAtmp);
          }
        };

        // Generate Loads from Matrix B
        auto const B_load = [&](int c, const string& vBcol, int prefetch_len) {
          if (d_type.first == DataType::Float32) {
            addi(
                srcfile,
                "vmovups " + vBcol + "," +
                    (s.iset == ISA::isaType::avx512 ? "ZMM" : "YMM") +
                    "WORD PTR [r10 + " + to_string(vec_len_in_bytes * c) + "]");
          } else if (d_type.first == DataType::Float16) {
            addi(
                srcfile,
                "vcvtph2ps " + vBcol + "," +
                    (s.iset == ISA::isaType::avx512 ? "YMM" : "XMM") +
                    "WORD PTR [r10 + " + to_string(vec_len_in_bytes / 2 * c) +
                    "]");
          }
          if (prefetch_len && ((vec_len_in_bytes * c) % cache_line_size == 0)) {
            addi(
                srcfile,
                "prefetcht0 [r10 + " +
                    to_string(vec_len_in_bytes * c + prefetch_len) + "]",
                fixedC);
          }
        };
>>>>>>> upstream

        auto const C_prefetch = [&](int r) {
          for (auto c = 0; prefetch_c_len && (c < vCtile[r].size()); c++) {
            if ((vec_len_in_bytes * c) % cache_line_size == 0) {
              addi(
                  srcfile,
<<<<<<< HEAD
                  "    auto " + r_spare + " = _mm256_broadcast_ss((float*)r15);",
                  "vbroadcastss " + r_spare + string(",DWORD PTR [r15]"),
                  fixedC);
              // store out C
              for (auto r = 0; r < vCtile.size(); r++) {
                for (auto c = 0; c < vCtile[r].size(); c++) {
                  switch (s.avx) {
                    case 1:
                      addi(
                          srcfile,
                          "not supported",
                          string("vmulps ymm15, ") + r_spare + comma +
                              "YMMWORD PTR [r12 + " + to_string(32 * c) + "]",
                          fixedC);
                      addi(
                          srcfile,
                          "not supported",
                          "vaddps " + vCtile[r][c] + "," + vCtile[r][c] + "," +
                              "ymm15",
                          fixedC);
                      break;
                    case 2:
                      //if (r == 0) {
                        addi(
                          srcfile,
                          ((r == 0) ? "    auto r12_" + to_string(32 * c) : "    r12_" + to_string(32 * c)) + " = _mm256_load_ps((float*)((char*)r12 + " + to_string(32 * c) + "));",
                          "vfmadd231ps " + vCtile[r][c] + "," + r_spare + "," +
                          "YMMWORD PTR [r12 + " + to_string(32 * c) + "]",
                          fixedC);
                      //}
                      addi(
                          srcfile,
                          "    " + vCtile[r][c] + " = _mm256_fmadd_ps(r12_" + to_string(32 * c) + ", " + r_spare + ", " + vCtile[r][c] + ");",
                          "vfmadd231ps " + vCtile[r][c] + "," + r_spare + "," +
                              "YMMWORD PTR [r12 + " + to_string(32 * c) + "]",
                          fixedC);
                      break;
                    default:
                      assert(0);
                  }
                  addi(
                      srcfile,
                      "    _mm256_storeu_ps((float*)((char*)r12 + " + to_string(32 * c) + "), " + vCtile[r][c] + ");",
                      "vmovups YMMWORD PTR [r12 + " + to_string(32 * c) +
                          "], " + vCtile[r][c],
                      fixedC);
                }
                if (r != vCtile.size() - 1)
                  addi(srcfile, "    r12 = (float*)((char*)r12 + r13);", "add r12, r13", fixedC); // move C ptr
              }

              //srcfile << "\n";
              addi(srcfile, "  }", "L_done%=:");

              srcfile << "\n        // next outer iteration\n";
              // C
              addi(
                  srcfile,
                  "  rcx = (float*)((char*)rcx + " + to_string(32 * ukernel_shape[k][1]) + ");",
                  "add rcx, " + to_string(32 * ukernel_shape[k][1]),
                  fixedC);
              addi(srcfile, "  r12 = rcx;", "mov r12, rcx", fixedC);
              // A
              addi(srcfile, "  r9 = rax;", "mov r9, rax");

              // move to top for looop
              //addi(srcfile, "inc rbx");
              //addi(srcfile, "cmp rbx, rdi");
              //addi(srcfile, "jl " + label2);
              addi(srcfile, "}", "inc rbx; cmp rbx, rdi; jl " + label2);

              //// output
              //srcfile << "      :\n";
              //// input
              //srcfile << "      : [gp] \"rm\"(gp)\n";

              //// clobbered
              //srcfile
              //    << "      : \"r8\",\n        \"r9\",\n        \"r10\",\n"
              //       "        \"r11\",\n        \"r15\",\n        \"r13\",\n"
              //       "        \"r14\",\n        \"rax\",\n        \"rcx\",\n"
              //       "        \"rdx\",\n        \"rsi\",\n        \"rdi\",\n"
              //       "        \"rbx\",\n        \"r12\",\n"
              //       "        \"memory\");\n";
              srcfile << "}\n\n";
            }

#else
            fheader[k] =
              "void __attribute__((noinline)) " + funcname[k] + fargs;
            srcfile << fheader[k] << " {\n";

            unsigned last_free_ymmreg = 0;
            // produce register block of C
            vector<vector<string>> vCtile(ukernel_shape[k][0]);
            for (auto r = 0; r < ukernel_shape[k][0]; r++)
              for (auto c = 0; c < ukernel_shape[k][1]; c++) {
                vCtile[r].push_back("ymm" + to_string(last_free_ymmreg));
                last_free_ymmreg++;
              }
            assert(last_free_ymmreg <= 14);

            string vAtmp = "ymm" + to_string(last_free_ymmreg++);
            // produce register block of B col
            vector<string> vBcol(ukernel_shape[k][1]);

            for (auto c = 0; c < ukernel_shape[k][1]; c++) {
              vBcol[c] = ("ymm" + to_string(last_free_ymmreg));
              last_free_ymmreg++;
            }

            assert(last_free_ymmreg <= 16);

            srcfile << "  asm volatile(\n";

            srcfile << "#if !defined(__clang__)"
              << "\n";
            addi(srcfile, "mov r14, %[gp]");
            srcfile << "#else\n";
            addi(srcfile, "mov %[gp], %%r14");
            addi(srcfile, ".intel_syntax noprefix");
            srcfile << "#endif\n";

            srcfile << "\n      // Copy parameters\n";
            srcfile << "      // k\n";
            addi(srcfile, "mov r8, [r14 + 0]");
            srcfile << "      // A\n";
            addi(srcfile, "mov r9, [r14 + 8]");
            srcfile << "      // B\n";
            addi(srcfile, "mov r10, [r14 + 16]");
            srcfile << "      // beta\n";
            addi(srcfile, "mov r15, [r14 + 24]");
            srcfile << "      // accum\n";
            addi(srcfile, "mov rdx, [r14 + 32]");
            srcfile << "      // C\n";
            addi(srcfile, "mov r12, [r14 + 40]");
            srcfile << "      // ldc\n";
            addi(srcfile, "mov r13, [r14 + 48]");
            srcfile << "      // b_block_cols\n";
            addi(srcfile, "mov rdi, [r14 + 56]");
            srcfile << "      // b_block_size\n";
            addi(srcfile, "mov rsi, [r14 + 64]");
            srcfile << "      // Make copies of A and C\n";
            addi(srcfile, "mov rax, r9");
            addi(srcfile, "mov rcx, r12");
            srcfile << "\n";

            addi(srcfile, "mov rbx, 0");

            string exitlabel = "L_exit%=";
            string label2 = "loop_outter%=";
            addi(srcfile, label2 + ":");
            addi(srcfile, "mov r14, 0");

            // set all vCtile regs to zeros
            for (auto r = 0; r < vCtile.size(); r++) {
              for (auto c = 0; c < vCtile[r].size(); c++) {
                addi(
                  srcfile,
                  "vxorps " + vCtile[r][c] + "," + vCtile[r][c] + "," +
                  vCtile[r][c]);
              }
            }

            // start marker
            if (iaca) {
              addi(srcfile, "mov ebx, 111");
              addi(srcfile, ".byte 0x64, 0x67, 0x90");
            }

            srcfile << "\n";

            srcfile << "\n";
            string label = "loop_inner%=";
            addi(srcfile, label + ":");
            srcfile << "\n";

            for (int c = 0; c < vCtile[0].size(); c++) {
              addi(
                srcfile,
                "vcvtph2ps " + vBcol[c] + ",XMMWORD PTR [r10 + " +
                to_string(16 * c) + "]");
=======
                  "prefetcht1 [r12 + " +
                      to_string(
                          /*vec_len_in_bytes * ukernel_shape[k][1] +*/
                          c * cache_line_size + prefetch_c_len) +
                      "]",
                  fixedC);
            }
          }
        };

        // Generate Loads from Matrix C
        auto const C_load = [&](int r) {
          for (auto c = 0; c < vCtile[r].size(); ++c) {
            switch (s.iset) {
              case ISA::isaType::avx:
              case ISA::isaType::avx2:
              case ISA::isaType::avx512:
              case ISA::isaType::avx512_256:
                if (prefetch_c_len &&
                    ((vec_len_in_bytes * c) % cache_line_size == 0)) {
                  addi(
                      srcfile,
                      "prefetcht1 [r12 + " +
                          to_string(
                              /*vec_len_in_bytes * ukernel_shape[k][1] +*/
                              c * cache_line_size + prefetch_c_len) +
                          "]",
                      fixedC);
                }
                addi(
                    srcfile,
                    "vmulps " + vCtile[r][c] + ", " + r_spare + ", " +
                        "[r12 + " + to_string(vec_len_in_bytes * c) + "]",
                    fixedC);
                break;
              default:
                assert(0);
>>>>>>> upstream
            }
          }
        };

        srcfile << "  asm volatile(\n";

        srcfile << "#if !defined(__clang__)"
                << "\n";
        addi(srcfile, "mov r14, %[gp]");
        srcfile << "#else\n";
        addi(srcfile, "mov %[gp], %%r14");
        addi(srcfile, ".intel_syntax noprefix");
        srcfile << "#endif\n";

        srcfile << "\n";
        srcfile << "      // Copy parameters\n";
        srcfile << "      // k\n";
        addi(srcfile, "mov r8, [r14 + 0]");
        // Assuming k >= 1
        addi(srcfile, "dec r8");
        srcfile << "      // A\n";
        addi(srcfile, "mov r9, [r14 + 8]");
        srcfile << "      // B\n";
        addi(srcfile, "mov r10, [r14 + 16]");
        srcfile << "      // beta\n";
        addi(srcfile, "lea r15, [r14 + 24]");
        srcfile << "      // C\n";
        addi(srcfile, "mov r12, [r14 + 32]");
        srcfile << "      // ldc\n";
        addi(srcfile, "mov r13, [r14 + 40]");
        srcfile << "      // b_block_cols\n";
        addi(srcfile, "mov rdi, [r14 + 48]");
        srcfile << "      // b_block_size\n";
        addi(srcfile, "mov rsi, [r14 + 56]");
        srcfile << "\n";
        srcfile << "      // Make copies of A and C\n";
        addi(srcfile, "mov rax, r9");
        addi(srcfile, "mov rcx, r12");
        srcfile << "\n";

        addi(srcfile, "mov rbx, 0");
        addi(srcfile, label_outer + ":");
        addi(srcfile, "mov r14, r8");

        string r_spare_cmp = "xmm" +
            to_string(num_vec_regs - (s.iset == ISA::isaType::avx ? 2 : 1));

        addi(
            srcfile,
            "vbroadcastss " + r_spare + string(",DWORD PTR [r15]"),
            fixedC);
        // Generate first iteration which loads values from C  and interleavs
        // With loads from B and multiplication
        for (auto c = 0; c < vCtile[0].size(); ++c) {
          B_load(c, vBcol[c], prefetch_b_len);
        }
        addi(srcfile, "vxorps xmm0, xmm0, xmm0");
        addi(srcfile, "vcomiss " + r_spare_cmp + ", xmm0");
        addi(srcfile, "jz " + label_zero);

        srcfile << "\n";
        srcfile << "      // Setup values with beta multiplication\n";
        string r_last = vec_reg_prefix + to_string(num_vec_regs - 1);
        for (auto r = 0; r < vCtile.size(); r++) {
          if (r > 0) {
            addi(srcfile, "add r12, r13", fixedC); // move C ptr
          }
          C_load(r);
        }
        // Skip matrix B preload if k == 1 (may OutOfBound access)
        addi(srcfile, "test r14,r14");
        addi(srcfile, "jz skip_preload%=");
        // Preload B index and prefetch with the next iteration
        B_load(vCtile[0].size(), r_spare, prefetch_b_len);
        addi(srcfile, "skip_preload%=:");
        for (auto r = 0; r < vCtile.size(); r++) {
          A_load_mult(r, mult_type::fma);
        }
        if (vCtile.size() > 1) {
          addi(srcfile, "mov r12, rcx");
        }
        addi(srcfile, "test r14,r14");                 // Decrease iterations
        addi(srcfile, "jnz " + label_next_inner);
        addi(srcfile, "add r10," + B_stride, fixedA);  // B stride
        addi(srcfile, "jmp " + label_dump_C);

        //
        // Handle non-accumulate case, the values can be directly stored
        //
        srcfile << "\n";
        addi(srcfile, label_zero + ":");
        srcfile << "\n";
        // Skip matrix B preload if k == 1 (may OutOfBound access)
        addi(srcfile, "test r14,r14");
        addi(srcfile, "jz skip_preload_b_zero%=");
        // Preload B index and with the next iteration
        B_load(vCtile[0].size(), r_spare, prefetch_b_len);
        addi(srcfile, "skip_preload_b_zero%=:");
        // Consider all vCtile regs as zeros, do direct MUL into
        for (auto r = 0; r < vCtile.size(); r++) {
          if (r > 0) {
            addi(srcfile, "add r12, r13", fixedC); // move C ptr
          }
          C_prefetch(r);
          A_load_mult(r, mult_type::mul);
        }
        if (vCtile.size() > 1) {
          addi(srcfile, "mov r12, rcx");
        }
        addi(srcfile, "test r14,r14");                   // Decrease iterations
        addi(srcfile, "jnz " + label_next_inner);
        addi(srcfile, "add r10," + B_stride, fixedA);    // B stride
        addi(srcfile, "jmp " + label_dump_C);

        // start marker
        if (iaca) {
          addi(srcfile, "mov ebx, 111");
          addi(srcfile, ".byte 0x64, 0x67, 0x90");
        }

<<<<<<< HEAD
            for (int r = 0; r < vCtile.size(); r++) {
              addi(
                srcfile,
                "vbroadcastss " + vAtmp + ",DWORD PTR [r9+" +
                to_string(4 * r) + "]");
              for (int c = 0; c < vCtile[0].size(); c++) {
                addi(
                  srcfile,
                  "vfmadd231ps " + vCtile[r][c] + "," + vBcol[c] + "," +
                  vAtmp);
              }
            }

            addi(
              srcfile,
              "add r9," + to_string(4 * ukernel_shape[k][0]),
              fixedA); // move A ptr

            addi(
              srcfile,
              "add r10," + to_string(16 * ukernel_shape[k][1]),
              fixedA); // move A ptr

            addi(srcfile, "inc r14");
            addi(srcfile, "cmp r14, r8");
            addi(srcfile, "jl " + label);

            srcfile << "\n";

            addi(srcfile, exitlabel + ":");

            // addi(srcfile, "add r10, rsi");
            srcfile << "\n";

            // end marker
            if (iaca) {
              addi(srcfile, "mov ebx, 222");
              addi(srcfile, ".byte 0x64, 0x67, 0x90");
            }

            addi(srcfile, "cmp rdx, 1");
            addi(srcfile, "je L_accum%=");
            srcfile << "      // Dump C\n";

            for (auto r = 0; r < vCtile.size(); r++) {
              for (auto c = 0; c < vCtile[r].size(); c++) {
                addi(
                  srcfile,
                  "vmovups YMMWORD PTR [r12 + " + to_string(32 * c) +
                  "], " + vCtile[r][c],
                  fixedC);
              }
              addi(srcfile, "add r12, r13", fixedC); // move C ptr
            }
            addi(srcfile, "jmp L_done%=");

            srcfile << "\n";
            addi(srcfile, "L_accum%=:");
            srcfile << "      // Dump C with accumulate\n";

            string r_spare = (s.avx == 1) ? "ymm14" : "ymm15";
            addi(
              srcfile,
              "vbroadcastss " + r_spare + string(",DWORD PTR [r15]"),
              fixedC);
            // store out C
            for (auto r = 0; r < vCtile.size(); r++) {
              for (auto c = 0; c < vCtile[r].size(); c++) {
                switch (s.avx) {
                case 1:
                  addi(
                    srcfile,
                    string("vmulps ymm15, ") + r_spare + comma +
                    "YMMWORD PTR [r12 + " + to_string(32 * c) + "]",
                    fixedC);
                  addi(
                    srcfile,
                    "vaddps " + vCtile[r][c] + "," + vCtile[r][c] + "," +
                    "ymm15",
                    fixedC);
                  break;
                case 2:
                  addi(
                    srcfile,
                    "vfmadd231ps " + vCtile[r][c] + "," + r_spare + "," +
                    "YMMWORD PTR [r12 + " + to_string(32 * c) + "]",
                    fixedC);
                  break;
                default:
                  assert(0);
                }
                addi(
                  srcfile,
                  "vmovups YMMWORD PTR [r12 + " + to_string(32 * c) +
                  "], " + vCtile[r][c],
                  fixedC);
              }
              addi(srcfile, "add r12, r13", fixedC); // move C ptr
            }

            srcfile << "\n";
            addi(srcfile, "L_done%=:");

            srcfile << "\n      // next outer iteration\n";
            // C
            addi(
              srcfile,
              "add rcx, " + to_string(32 * ukernel_shape[k][1]),
              fixedC);
            addi(srcfile, "mov r12, rcx", fixedC);
            // A
            addi(srcfile, "mov r9, rax");

            addi(srcfile, "inc rbx");
            addi(srcfile, "cmp rbx, rdi");
            addi(srcfile, "jl " + label2);

            // output
            srcfile << "      :\n";
            // input
            srcfile << "      : [gp] \"rm\"(gp)\n";

            // clobbered
            srcfile
              << "      : \"r8\",\n        \"r9\",\n        \"r10\",\n"
              "        \"r11\",\n        \"r15\",\n        \"r13\",\n"
              "        \"r14\",\n        \"rax\",\n        \"rcx\",\n"
              "        \"rdx\",\n        \"rsi\",\n        \"rdi\",\n"
              "        \"rbx\",\n        \"r12\",\n"
              "        \"memory\");\n";
            srcfile << "}\n";
            }

#endif

            for (unsigned k = 0; k < ukernel_shape.size(); k++) {
              hdrfile << fheader[k] << ";\n";
            }
=======
        //
        //  Inner iteration begin
        //
        srcfile << "\n";
        addi(srcfile, label_inner + ":");
        srcfile << "\n";

        // Store preloaded value
        addi(srcfile, "vmovaps " + vBcol[0] + "," + r_spare);
        for (int c = 1; c < vCtile[0].size(); c++) {
          B_load(c, vBcol[c], prefetch_b_len);
        }
        // Preload for next iteration
        B_load(vCtile[0].size(), r_spare, prefetch_b_len);
        for (int r = 0; r < vCtile.size(); r++) {
          A_load_mult(r, mult_type::fma);
        }
>>>>>>> upstream

        // Finish inner iteration
        srcfile << "\n";
        addi(srcfile, label_next_inner + ":");
        addi(srcfile, "add r9," + A_stride, fixedA);  // A stride
        addi(srcfile, "add r10," + B_stride, fixedA); // B stride
        addi(srcfile, "dec r14");                     // Decrease iterations
        addi(srcfile, "jnz " + label_inner);
        srcfile << "\n";

        // end marker
        if (iaca) {
          addi(srcfile, "mov ebx, 222");
          addi(srcfile, ".byte 0x64, 0x67, 0x90");
        }

        // Perform last iteration without preloading B values
        // Store preloaded value
        addi(srcfile, "vmovaps " + vBcol[0] + "," + r_spare);
        for (int c = 1; c < vCtile[0].size(); c++) {
          B_load(c, vBcol[c], 0); // no prefetch
        }
        for (int r = 0; r < vCtile.size(); r++) {
          A_load_mult(r, mult_type::fma);
        }
        addi(srcfile, "add r9," + A_stride, fixedA);  // A stride
        addi(srcfile, "add r10," + B_stride, fixedA); // B stride

        srcfile << "      // Dump C\n";
        addi(srcfile, label_dump_C + ":");
        for (auto r = 0; r < vCtile.size(); r++) {
          if (r > 0) {
            addi(srcfile, "add r12, r13", fixedC); // move C ptr
          }
          for (auto c = 0; c < vCtile[r].size(); c++) {
            addi(
                srcfile,
                "vmovups " + vec_reg_prefix + "word PTR [r12 + " +
                    to_string(vec_len_in_bytes * c) + "], " + vCtile[r][c],
                fixedC);
          }
        }

        srcfile << "\n      // next outer iteration\n";
        // C
        addi(
            srcfile,
            "add rcx, " + to_string(vec_len_in_bytes * ukernel_shape[k][1]),
            fixedC);
        addi(srcfile, "mov r12, rcx", fixedC);
        // A
        addi(srcfile, "mov r9, rax");

        addi(srcfile, "inc rbx");
        addi(srcfile, "cmp rbx, rdi");
        addi(srcfile, "jl " + label_outer);

        // output
        srcfile << "      :\n";
        // input
        srcfile << "      : [gp] \"rm\"(gp)\n";

        // clobbered
        srcfile << "      : \"r8\",\n        \"r9\",\n        \"r10\",\n"
                   "        \"r11\",\n        \"r13\",\n"
                   "        \"r14\",\n        \"rax\",\n        \"rcx\",\n"
                   "        \"rsi\",\n        \"rdi\",\n"
                   "        \"rbx\",\n        \"r12\",\n        \"r15\",\n"
                   "        \"memory\");\n";
        srcfile << "}\n";
      }

      for (unsigned k = 0; k < ukernel_shape.size(); k++) {
        hdrfile << fheader[k] << ";\n";
      }

      srcfile << "\n} // namespace fbgemm\n";
      srcfile.close();
      hdrfile << "\n} // namespace fbgemm\n";
      hdrfile.close();
    } // isa
  }
}
