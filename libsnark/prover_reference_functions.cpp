#include <cassert>
#include <cstdio>
#include <fstream>
#include <fstream>
#include <libff/common/rng.hpp>
#include <libff/common/profiling.hpp>
#include <libff/common/utils.hpp>
#include <libsnark/serialization.hpp>
#include <libff/algebra/curves/mnt753/mnt4753/mnt4753_pp.hpp>
#include <libff/algebra/curves/mnt753/mnt6753/mnt6753_pp.hpp>
#include <omp.h>
#include <libff/algebra/scalar_multiplication/multiexp.hpp>
#include <libsnark/knowledge_commitment/kc_multiexp.hpp>
#include <libsnark/knowledge_commitment/knowledge_commitment.hpp>
#include <libsnark/reductions/r1cs_to_qap/r1cs_to_qap.hpp>

#include <libsnark/zk_proof_systems/ppzksnark/r1cs_gg_ppzksnark/r1cs_gg_ppzksnark.hpp>

#include <libfqfft/evaluation_domain/domains/basic_radix2_domain.hpp>

using namespace libff;
using namespace libsnark;

const multi_exp_method method = multi_exp_method_BDLO12;

template<typename ppT>
class groth16_parameters {
  public:
    size_t d;
    size_t m;
    std::vector<G1<ppT>> A, B1, L, H;
    std::vector<G2<ppT>> B2;

  groth16_parameters(const char* path) {
    FILE* params = fopen(path, "r");
    d = read_size_t(params);
    m = read_size_t(params);
    for (size_t i = 0; i <= m; ++i) { A.emplace_back(read_g1<ppT>(params)); }
    for (size_t i = 0; i <= m; ++i) { B1.emplace_back(read_g1<ppT>(params)); }
    for (size_t i = 0; i <= m; ++i) { B2.emplace_back(read_g2<ppT>(params)); }
    for (size_t i = 0; i < m-1; ++i) { L.emplace_back(read_g1<ppT>(params)); }
    for (size_t i = 0; i < d; ++i) { H.emplace_back(read_g1<ppT>(params)); }
    fclose(params);
  }
};

template<typename ppT>
class groth16_input {
  public:
    std::vector<Fr<ppT>> w;
    std::vector<Fr<ppT>> ca, cb, cc;
    Fr<ppT> r;

  groth16_input(const char* path, size_t d, size_t m) {
    FILE* inputs = fopen(path, "r");

    for (size_t i = 0; i < m + 1; ++i) { w.emplace_back(read_fr<ppT>(inputs)); }

    for (size_t i = 0; i < d + 1; ++i) { ca.emplace_back(read_fr<ppT>(inputs)); }
    for (size_t i = 0; i < d + 1; ++i) { cb.emplace_back(read_fr<ppT>(inputs)); }
    for (size_t i = 0; i < d + 1; ++i) { cc.emplace_back(read_fr<ppT>(inputs)); }

    r = read_fr<ppT>(inputs);

    fclose(inputs);
  }
};

template<typename ppT>
class groth16_output {
  public:
    G1<ppT> A, C;
    G2<ppT> B;

  groth16_output(G1<ppT> &&A, G2<ppT> &&B, G1<ppT> &&C) :
    A(std::move(A)), B(std::move(B)), C(std::move(C)) {}

  void write(const char* path) {
    FILE* out = fopen(path, "w");
    write_g1<ppT>(out, A);
    write_g2<ppT>(out, B);
    write_g1<ppT>(out, C);
    fclose(out);
  }
};

// Here is where all the FFTs happen.
template<typename ppT>
std::vector<Fr<ppT>> compute_H(
    size_t d,
    std::vector<Fr<ppT>> &ca,
    std::vector<Fr<ppT>> &cb,
    std::vector<Fr<ppT>> &cc)
{
    // Begin witness map
    libff::enter_block("Compute the polynomial H");

    const std::shared_ptr<libfqfft::evaluation_domain<Fr<ppT>> > domain = libfqfft::get_evaluation_domain<Fr<ppT>>(d + 1);

    domain->iFFT(ca);
    domain->iFFT(cb);

    domain->cosetFFT(ca, Fr<ppT>::multiplicative_generator);
    domain->cosetFFT(cb, Fr<ppT>::multiplicative_generator);

    libff::enter_block("Compute evaluation of polynomial H on set T");
    std::vector<Fr<ppT>> &H_tmp = ca; // can overwrite ca because it is not used later
#ifdef MULTICORE
#pragma omp parallel for
#endif
    for (size_t i = 0; i < domain->m; ++i)
    {
        H_tmp[i] = ca[i]*cb[i];
    }
    std::vector<Fr<ppT>>().swap(cb); // destroy cb

    domain->iFFT(cc);

    domain->cosetFFT(cc, Fr<ppT>::multiplicative_generator);

#ifdef MULTICORE
#pragma omp parallel for
#endif
    for (size_t i = 0; i < domain->m; ++i)
    {
        H_tmp[i] = (H_tmp[i]-cc[i]);
    }

    domain->divide_by_Z_on_coset(H_tmp);

    libff::leave_block("Compute evaluation of polynomial H on set T");

    domain->icosetFFT(H_tmp, Fr<ppT>::multiplicative_generator);

    std::vector<Fr<ppT>> coefficients_for_H(domain->m+1, Fr<ppT>::zero());
#ifdef MULTICORE
#pragma omp parallel for
#endif
    for (size_t i = 0; i < domain->m; ++i)
    {
        coefficients_for_H[i] = H_tmp[i];
    }

    libff::leave_block("Compute the polynomial H");

    return coefficients_for_H;
}

template<typename G, typename Fr>
G multiexp(typename std::vector<Fr>::const_iterator scalar_start,
           typename std::vector<G>::const_iterator g_start,
           size_t length)
{
#ifdef MULTICORE
    const size_t chunks = omp_get_max_threads(); // to override, set OMP_NUM_THREADS env var or call omp_set_num_threads()
#else
    const size_t chunks = 1;
#endif

    return libff::multi_exp_with_mixed_addition<G,
                                                Fr,
                                                method>(
        g_start,
        g_start + length,
        scalar_start,
        scalar_start + length,
        chunks);

}

template<typename ppT>
int run_prover(
    const char* params_path,
    const char* input_path,
    const char* output_path)
{
    ppT::init_public_params();

    const size_t primary_input_size = 1;

    const groth16_parameters<ppT> parameters(params_path);
    const groth16_input<ppT> input(input_path, parameters.d, parameters.m);

    std::vector<Fr<ppT>> w  = std::move(input.w);
    std::vector<Fr<ppT>> ca = std::move(input.ca);
    std::vector<Fr<ppT>> cb = std::move(input.cb);
    std::vector<Fr<ppT>> cc = std::move(input.cc);

    // End reading of parameters and input

    libff::enter_block("Call to r1cs_gg_ppzksnark_prover");

    std::vector<Fr<ppT>> coefficients_for_H = compute_H<ppT>(
        parameters.d,
        ca, cb, cc);

    libff::enter_block("Compute the proof");
    libff::enter_block("Multi-exponentiations");

    // Now the 5 multi-exponentiations
    G1<ppT> evaluation_At = multiexp<G1<ppT>, Fr<ppT>>(
        w.begin(), parameters.A.begin(), parameters.m + 1);

    G1<ppT> evaluation_Bt1 = multiexp<G1<ppT>, Fr<ppT>>(
        w.begin(), parameters.B1.begin(), parameters.m + 1);

    G2<ppT> evaluation_Bt2 = multiexp<G2<ppT>, Fr<ppT>>(
        w.begin(), parameters.B2.begin(), parameters.m + 1);

    G1<ppT> evaluation_Ht = multiexp<G1<ppT>, Fr<ppT>>(
        coefficients_for_H.begin(), parameters.H.begin(), parameters.d);

    G1<ppT> evaluation_Lt = multiexp<G1<ppT>, Fr<ppT>>(
        w.begin() + primary_input_size + 1,
        parameters.L.begin(),
        parameters.m - 1);

    libff::G1<ppT> C = evaluation_Ht + evaluation_Lt + input.r * evaluation_Bt1; /*+ s *  g1_A  - (r * s) * pk.delta_g1; */

    libff::leave_block("Multi-exponentiations");
    libff::leave_block("Compute the proof");
    libff::leave_block("Call to r1cs_gg_ppzksnark_prover");

    groth16_output<ppT> output(
      std::move(evaluation_At),
      std::move(evaluation_Bt2),
      std::move(C));

    output.write(output_path);

    return 0;
}

