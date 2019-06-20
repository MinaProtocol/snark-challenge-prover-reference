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

// declare the wrapper functions from cuda-fixnum/main.cu
int do_fixnum_example(const char *inputs_file, const char *outputs_file);

int main(int argc, const char * argv[])
{
  setbuf(stdout, NULL);
  std::string curve(argv[1]);
  std::string mode(argv[2]);
  std::string device(argv[3]);


  if (device == "CPU") {
    const char* params_path = argv[4];
    const char* input_path = argv[5];
    const char* output_path = argv[6];

    if (curve == "MNT4753") {
        if (mode == "compute") {
        return run_prover<mnt4753_libsnark>(params_path, input_path, output_path);
        }
    } else if (curve == "MNT6753") {
        if (mode == "compute") {
        return run_prover<mnt6753_libsnark>(params_path, input_path, output_path);
        }
    }
  } else if (device == "GPU") {
    do_fixnum_example(argv[4], argv[5]);
  }
}

template<typename ppT>
void debug(
    Fr<ppT>& r,
    groth16_output<ppT>& output,
    std::vector<Fr<ppT>>& w) {

    const size_t primary_input_size = 1;

    std::vector<Fr<ppT>> primary_input(w.begin() + 1, w.begin() + 1 + primary_input_size);
    std::vector<Fr<ppT>> auxiliary_input(w.begin() + 1 + primary_input_size, w.end() );

    const libff::Fr<ppT> s = libff::Fr<ppT>::random_element();

    r1cs_gg_ppzksnark_proving_key<ppT> pk;
    std::ifstream pk_debug;
    pk_debug.open("proving-key.debug");
    pk_debug >> pk;

    /* A = alpha + sum_i(a_i*A_i(t)) + r*delta */
    libff::G1<ppT> g1_A = pk.alpha_g1 + output.A + r * pk.delta_g1;

    /* B = beta + sum_i(a_i*B_i(t)) + s*delta */
    libff::G2<ppT> g2_B = pk.beta_g2 + output.B + s * pk.delta_g2;

    /* C = sum_i(a_i*((beta*A_i(t) + alpha*B_i(t) + C_i(t)) + H(t)*Z(t))/delta) + A*s + r*b - r*s*delta */
    libff::G1<ppT> g1_C = output.C + s * g1_A + r * pk.beta_g1;

    libff::leave_block("Compute the proof");

    libff::leave_block("Call to r1cs_gg_ppzksnark_prover");

    r1cs_gg_ppzksnark_proof<ppT> proof = r1cs_gg_ppzksnark_proof<ppT>(std::move(g1_A), std::move(g2_B), std::move(g1_C));
    proof.print_size();

    r1cs_gg_ppzksnark_verification_key<ppT> vk;
    std::ifstream vk_debug;
    vk_debug.open("verification-key.debug");
    vk_debug >> vk;
    vk_debug.close();

    assert (r1cs_gg_ppzksnark_verifier_strong_IC<ppT>(vk, primary_input, proof) );

    r1cs_gg_ppzksnark_proof<ppT> proof1=
      r1cs_gg_ppzksnark_prover<ppT>(
          pk, 
          primary_input,
          auxiliary_input);
    assert (r1cs_gg_ppzksnark_verifier_strong_IC<ppT>(vk, primary_input, proof1) );
}
