#include <cassert>
#include <cstdio>
#include <fstream>
#include <libff/algebra/curves/mnt753/mnt4753/mnt4753_pp.hpp>
#include <libff/algebra/curves/mnt753/mnt6753/mnt6753_pp.hpp>
#include <libff/algebra/scalar_multiplication/multiexp.hpp>
#include <libff/common/profiling.hpp>
#include <libff/common/rng.hpp>
#include <libff/common/utils.hpp>
#include <libsnark/knowledge_commitment/kc_multiexp.hpp>
#include <libsnark/knowledge_commitment/knowledge_commitment.hpp>
#include <libsnark/reductions/r1cs_to_qap/r1cs_to_qap.hpp>
#include <libsnark/serialization.hpp>
#include <omp.h>

#include <libsnark/zk_proof_systems/ppzksnark/r1cs_gg_ppzksnark/r1cs_gg_ppzksnark.hpp>

#include <libfqfft/evaluation_domain/domains/basic_radix2_domain.hpp>

//#include "prover_reference_include/prover_reference_functions.hpp"

using namespace libff;
using namespace libsnark;

const multi_exp_method method = multi_exp_method_BDLO12;

// Here is where all the FFTs happen.
template <typename ppT>
std::vector<Fr<ppT>> compute_H(size_t d, std::vector<Fr<ppT>> &ca,
                               std::vector<Fr<ppT>> &cb,
                               std::vector<Fr<ppT>> &cc) {
  // Begin witness map
  libff::enter_block("Compute the polynomial H");

  const std::shared_ptr<libfqfft::evaluation_domain<Fr<ppT>>> domain =
      libfqfft::get_evaluation_domain<Fr<ppT>>(d + 1);

  domain->iFFT(ca);
  domain->iFFT(cb);

  domain->cosetFFT(ca, Fr<ppT>::multiplicative_generator);
  domain->cosetFFT(cb, Fr<ppT>::multiplicative_generator);

  libff::enter_block("Compute evaluation of polynomial H on set T");
  std::vector<Fr<ppT>> &H_tmp =
      ca; // can overwrite ca because it is not used later
#ifdef MULTICORE
#pragma omp parallel for
#endif
  for (size_t i = 0; i < domain->m; ++i) {
    H_tmp[i] = ca[i] * cb[i];
  }
  std::vector<Fr<ppT>>().swap(cb); // destroy cb

  domain->iFFT(cc);

  domain->cosetFFT(cc, Fr<ppT>::multiplicative_generator);

#ifdef MULTICORE
#pragma omp parallel for
#endif
  for (size_t i = 0; i < domain->m; ++i) {
    H_tmp[i] = (H_tmp[i] - cc[i]);
  }

  domain->divide_by_Z_on_coset(H_tmp);

  libff::leave_block("Compute evaluation of polynomial H on set T");

  domain->icosetFFT(H_tmp, Fr<ppT>::multiplicative_generator);

  std::vector<Fr<ppT>> coefficients_for_H(domain->m + 1, Fr<ppT>::zero());
#ifdef MULTICORE
#pragma omp parallel for
#endif
  for (size_t i = 0; i < domain->m; ++i) {
    coefficients_for_H[i] = H_tmp[i];
  }

  libff::leave_block("Compute the polynomial H");

  return coefficients_for_H;
}

template <typename G, typename Fr>
G multiexp(typename std::vector<Fr>::const_iterator scalar_start,
           typename std::vector<G>::const_iterator g_start, size_t length) {
#ifdef MULTICORE
  const size_t chunks =
      omp_get_max_threads(); // to override, set OMP_NUM_THREADS env var or call
                             // omp_set_num_threads()
#else
  const size_t chunks = 1;
#endif

  return libff::multi_exp_with_mixed_addition<G, Fr, method>(
      g_start, g_start + length, scalar_start, scalar_start + length, chunks);
}

template <typename ppT> class libsnark_bundle {
  class groth16_input {
  public:
    std::shared_ptr<std::vector<Fr<ppT>>> w;
    std::shared_ptr<std::vector<Fr<ppT>>> ca, cb, cc;
    Fr<ppT> r;

    groth16_input(const char *path, size_t d, size_t m) {
      FILE *inputs = fopen(path, "r");

      for (size_t i = 0; i < m + 1; ++i) {
        w.emplace_back(read_fr<ppT>(inputs));
      }

      for (size_t i = 0; i < d + 1; ++i) {
        ca.emplace_back(read_fr<ppT>(inputs));
      }
      for (size_t i = 0; i < d + 1; ++i) {
        cb.emplace_back(read_fr<ppT>(inputs));
      }
      for (size_t i = 0; i < d + 1; ++i) {
        cc.emplace_back(read_fr<ppT>(inputs));
      }

      r = read_fr<ppT>(inputs);

      fclose(inputs);
    }
  };

  class groth16_params {
  public:
    size_t d;
    size_t m;
    std::shared_ptr<std::vector<libff::G1<ppT>>> A, B1, L, H;
    std::shared_ptr<std::vector<libff::G2<ppT>>> B2;

    groth16_params(const char *path) {
      FILE *params = fopen(path, "r");
      d = read_size_t(params);
      m = read_size_t(params);
      for (size_t i = 0; i <= m; ++i) {
        A.emplace_back(read_g1<ppT>(params));
      }
      for (size_t i = 0; i <= m; ++i) {
        B1.emplace_back(read_g1<ppT>(params));
      }
      for (size_t i = 0; i <= m; ++i) {
        B2.emplace_back(read_g2<ppT>(params));
      }
      for (size_t i = 0; i < m - 1; ++i) {
        L.emplace_back(read_g1<ppT>(params));
      }
      for (size_t i = 0; i < d; ++i) {
        H.emplace_back(read_g1<ppT>(params));
      }
      fclose(params);
    }
  };

  class groth16_output {
  public:
      libff::G1<ppT> A, C;
      libff::G2<ppT> B;

    groth16_output(libff::G1<ppT> &&A, libff::G2<ppT> &&B, libff::G1<ppT> &&C)
        : A(std::move(A)), B(std::move(B)), C(std::move(C)) {}

    void write(const char *path) {
      FILE *out = fopen(path, "w");
      write_g1<ppT>(out, A);
      write_g2<ppT>(out, B);
      write_g1<ppT>(out, C);
      fclose(out);
    }
  };

  struct evaluation_domain {
    std::shared_ptr<libfqfft::evaluation_domain<Fr<ppT>>> data;
  };

  struct field {
    Fr<ppT> data;
  };

  struct G1 {
    libff::G1<ppT> data;
  };

  struct G2 {
    libff::G2<ppT> data;
  };

  struct vector_Fr {
    std::shared_ptr<std::vector<Fr<ppT>>> data;
    size_t offset;
  };

  struct vector_G1 {
    std::shared_ptr<std::vector<libff::G1<ppT>>> data;
  };
  struct vector_G2 {
    std::shared_ptr<std::vector<libff::G2<ppT>>> data;
  };

  static void init_public_params() { ppT::init_public_params(); }

  static evaluation_domain *get_evaluation_domain(size_t d) {
    return new evaluation_domain{
        .data = libfqfft::get_evaluation_domain<Fr<ppT>>(d + 1)};
  }

  static G1 *G1_add(G1 *a, G1 *b) { return new G1{.data = a->data + b->data}; }

  static G1 *G1_scale(field *a, G1 *b) {
    return new G1 { .data = a->data * b->data };
  }

  static void vector_Fr_muleq(vector_Fr *a, vector_Fr *b, size_t size) {
    size_t a_off = a->offset, b_off = b->offset;
#ifdef MULTICORE
#pragma omp parallel for
#endif
    for (size_t i = 0; i < size; i++) {
      a->data[i + a_off] = a->data[i + a_off] * b->data[i + b_off];
    }
  }

  static void vector_Fr_subeq(vector_Fr *a, vector_Fr *b, size_t size) {
    size_t a_off = a->offset, b_off = b->offset;
#ifdef MULTICORE
#pragma omp parallel for
#endif
    for (size_t i = 0; i < size; i++) {
      a->data[i + a_off] = a->data[i + a_off] * b->data[i + b_off];
    }
  }

  static vector_Fr *vector_Fr_offset(vector_Fr *a, size_t offset) {
    return new vector_Fr{.data = a->data, .offset = offset};
  }

  static vector_Fr *vector_Fr_copy(vector_Fr *a, size_t length) {
        auto new_data = std::make_shared<std::vector<Fr<ppT>>>(std::vector<Fr<ppT>>(a->data.begin(), a->data.begin()+length));
        return new vector_Fr { .data= new_data, .offset= 0};
  }

  static void domain_iFFT(evaluation_domain *domain, vector_Fr *a) {
    domain->data->iFFT(&a->data);
  }
  static void domain_cosetFFT(evaluation_domain *domain, vector_Fr *a) {
    domain->data->cosetFFT(&a->data, Fr<ppT>::multiplicative_generator);
  }
  static void domain_icosetFFT(evaluation_domain *domain, vector_Fr *a) {
    domain->data->cosetFFT(&a->data);
  }
  static void domain_divide_by_Z_on_coset(evaluation_domain *domain,
                                          vector_Fr *a) {
    domain->data->divide_by_Z_on_coset(&a->data);
  }
  static size_t domain_get_m(evaluation_domain *domain) {
    return domain->data.m;
  }

  static G1 *multiexp_G1(vector_Fr *scalar_start, vector_G1 *g_start,
                         size_t length) {

    return multiexp<libff::G1<ppT>, Fr<ppT>>(scalar_start->data.begin() +
                                                 scalar_start->offset,
                                             g_start->data.begin(), length);
  }
  static G2 *multiexp_G2(vector_Fr *scalar_start, vector_G2 *g_start,
                         size_t length) {
    return multiexp<libff::G1<ppT>, Fr<ppT>>(scalar_start->data.begin() +
                                                 scalar_start->offset,
                                             g_start->data.begin(), length);
  }

  static groth16_input *read_input(const char *path,
                                   groth16_params *params) {
    return new groth16_input(path, params->data.d, params->data.m);
  }

  static vector_Fr *input_w(groth16_input *input) {
    return new vector_Fr{.data = input->data.w, .offset = 0};
  }
  static vector_Fr *input_ca(groth16_input *input) {
    return new vector_Fr{.data = input->data.ca, .offset = 0};
  }
  static vector_Fr *input_cb(groth16_input *input) {
    return new vector_Fr{.data = input->data.cb, .offset = 0};
  }
  static vector_Fr *input_cc(groth16_input *input) {
    return new vector_Fr{.data = input->data.cc, .offset = 0};
  }
  static field *input_r(groth16_input *input) {
    return new field{.data = input->data.r};
  }

  static groth16_params *read_params(const char *path) {
    return new groth16_params(path);
  }

  static size_t params_d(groth16_params *params) { return params->d; }
  static size_t params_m(groth16_params *params) { return params->m; }
  static vector_G1 *params_A(groth16_params *params) {
    return new vector_G1{.data = params->A};
  }
  static vector_G1 *params_B1(groth16_params *params) {
    return new vector_G1{.data = params->B1};
  }
  static vector_G1 *params_L(groth16_params *params) {
    return new vector_G1{.data = params->L};
  }
  static vector_G1 *params_H(groth16_params *params) {
    return new vector_G1{.data = params->H};
  }
  static vector_G2 *params_B2(groth16_params *params) {
    return new vector_G2{.data = params->B2};
  }

  static void delete_G1(G1 *a) { delete a; }
  static void delete_G2(G1 *a) { delete a; }
  static void delete_vector_Fr(vector_Fr *a) { delete a; }
  static void delete_vector_G1(vector_G1 *a) { delete a; }
  static void delete_vector_G2(vector_G2 *a) { delete a; }
  static void delete_groth16_input(groth16_input *a) { delete a; }
  static void delete_groth16_params(groth16_params *a) { delete a; }
  static void delete_groth16_output(groth16_output *a) { delete a; }
  static void delete_evaluation_domain(evaluation_domain *a) { delete a; }

  static groth16_output *groth16_output_create(G1 *At, G2 *Bt2, G1 *C) {
    return new groth16_output(&At->data, &Bt2->data, &C->data);
  }
  static void groth16_output_write(groth16_output *output,
                                   const char *output_path) {
    output->write(output_path);
  }
};

class mnt4753_libsnark : public libsnark_bundle<mnt4753_pp> { };
class mnt6753_libsnark : public libsnark_bundle<mnt6753_pp> { };
