#pragma once

#include <vector>
#include <CL/cl.h>
// #define __CL_ENABLE_EXCEPTIONS
#define DATA_SIZE (131072)
#define limbs_per_elem (12)
#include <typeinfo>
#include <unistd.h>
#include <string.h>
#include <errno.h>
#include <exception>
#include <chrono>
#include "ocl_kernels/kernel.hpp"
using namespace std::chrono;
using namespace std;





class mnt4753_libsnark {
public:
  class groth16_input;

  class groth16_params;

  struct evaluation_domain;

  struct field;

  struct G1;

  struct G2;

  struct vector_Fr;

  struct vector_G1;

  struct vector_G2;

  static void init_public_params();

  static void print_G1(G1 *a);
  static void print_G2(G2 *a);

  static evaluation_domain *get_evaluation_domain(size_t d);

  static G1 *G1_add(G1 *a, G1 *b);
  static G1 *G1_scale(field *a, G1 *b);

  static void vector_Fr_muleq(vector_Fr *a, vector_Fr *b, size_t size);
  static void vector_Fr_subeq(vector_Fr *a, vector_Fr *b, size_t size);
  static vector_Fr *vector_Fr_offset(vector_Fr *a, size_t offset);
  static void vector_Fr_copy_into(vector_Fr *src, vector_Fr *dst,
                                  size_t length);
  static vector_Fr *vector_Fr_zeros(size_t length);

  static void domain_iFFT(evaluation_domain *domain, vector_Fr *a);
  static void domain_iFFT_GPU(evaluation_domain *domain, vector_Fr *a, Kernel k);
  static void domain_cosetFFT(evaluation_domain *domain, vector_Fr *a);
  static void domain_icosetFFT(evaluation_domain *domain, vector_Fr *a);
  static void domain_divide_by_Z_on_coset(evaluation_domain *domain,
                                          vector_Fr *a);
  static size_t domain_get_m(evaluation_domain *domain);

  static G1 *multiexp_G1(vector_Fr *scalar_start, vector_G1 *g_start,
                         size_t length);
  static G1 *multiexp_G1_GPU(vector_Fr *scalar_start, vector_G1 *g_start,
                         size_t length, Kernel k);
  static G2 *multiexp_G2(vector_Fr *scalar_start, vector_G2 *g_start,
                         size_t length);

  static groth16_input *read_input(const char *path, groth16_params *params);

  static vector_Fr *input_w(groth16_input *input);
  static vector_Fr *input_ca(groth16_input *input);
  static vector_Fr *input_cb(groth16_input *input);
  static vector_Fr *input_cc(groth16_input *input);
  static field *input_r(groth16_input *input);

  static groth16_params *read_params(const char *path);

  static size_t params_d(groth16_params *params);
  static size_t params_m(groth16_params *params);
  static vector_G1 *params_A(groth16_params *params);
  static vector_G1 *params_B1(groth16_params *params);
  static vector_G1 *params_L(groth16_params *params);
  static vector_G1 *params_H(groth16_params *params);
  static vector_G2 *params_B2(groth16_params *params);

  static void delete_G1(G1 *a);
  static void delete_G2(G1 *a);
  static void delete_vector_Fr(vector_Fr *a);
  static void delete_vector_G1(vector_G1 *a);
  static void delete_vector_G2(vector_G2 *a);
  static void delete_groth16_input(groth16_input *a);
  static void delete_groth16_params(groth16_params *a);
  static void delete_evaluation_domain(evaluation_domain *a);

  static void groth16_output_write(G1 *A, G2 *B, G1 *C,
                                   const char *output_path);
};
class mnt6753_libsnark {
public:
  class groth16_input;

  class groth16_params;

  struct evaluation_domain;

  struct field;

  struct G1;

  struct G2;

  struct vector_Fr;

  struct vector_G1;

  struct vector_G2;

  static void init_public_params();

  static void print_G1(G1 *a);
  static void print_G2(G2 *a);

  static evaluation_domain *get_evaluation_domain(size_t d);

  static G1 *G1_add(G1 *a, G1 *b);
  static G1 *G1_scale(field *a, G1 *b);

  static void vector_Fr_muleq(vector_Fr *a, vector_Fr *b, size_t size);
  static void vector_Fr_subeq(vector_Fr *a, vector_Fr *b, size_t size);
  static vector_Fr *vector_Fr_offset(vector_Fr *a, size_t offset);
  static void vector_Fr_copy_into(vector_Fr *src, vector_Fr *dst,
                                  size_t length);
  static vector_Fr *vector_Fr_zeros(size_t length);

  static void domain_iFFT(evaluation_domain *domain, vector_Fr *a);
  static void domain_iFFT_GPU(evaluation_domain *domain, vector_Fr *a, Kernel k);
  static void domain_cosetFFT(evaluation_domain *domain, vector_Fr *a);
  static void domain_icosetFFT(evaluation_domain *domain, vector_Fr *a);
  static void domain_divide_by_Z_on_coset(evaluation_domain *domain,
                                          vector_Fr *a);
  static size_t domain_get_m(evaluation_domain *domain);

  static G1 *multiexp_G1(vector_Fr *scalar_start, vector_G1 *g_start,
                         size_t length);
  static G1 *multiexp_G1_GPU(vector_Fr *scalar_start, vector_G1 *g_start,
                         size_t length, Kernel k);
  static G2 *multiexp_G2(vector_Fr *scalar_start, vector_G2 *g_start,
                         size_t length);

  static groth16_input *read_input(const char *path, groth16_params *params);

  static vector_Fr *input_w(groth16_input *input);
  static vector_Fr *input_ca(groth16_input *input);
  static vector_Fr *input_cb(groth16_input *input);
  static vector_Fr *input_cc(groth16_input *input);
  static field *input_r(groth16_input *input);

  static groth16_params *read_params(const char *path);

  static size_t params_d(groth16_params *params);
  static size_t params_m(groth16_params *params);
  static vector_G1 *params_A(groth16_params *params);
  static vector_G1 *params_B1(groth16_params *params);
  static vector_G1 *params_L(groth16_params *params);
  static vector_G1 *params_H(groth16_params *params);
  static vector_G2 *params_B2(groth16_params *params);

  static void delete_G1(G1 *a);
  static void delete_G2(G1 *a);
  static void delete_vector_Fr(vector_Fr *a);
  static void delete_vector_G1(vector_G1 *a);
  static void delete_vector_G2(vector_G2 *a);
  static void delete_groth16_input(groth16_input *a);
  static void delete_groth16_params(groth16_params *a);
  static void delete_evaluation_domain(evaluation_domain *a);

  static void groth16_output_write(G1 *A, G2 *B, G1 *C,
                                   const char *output_path);
};
