#include <cassert>
#include <cstdio>
#include <fstream>
#include <fstream>
#include <libff/common/rng.hpp>
#include <libff/common/profiling.hpp>
#include <libff/common/utils.hpp>

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
// const multi_exp_method method = multi_exp_method_bos_coster;

void write_mnt4_fq(FILE* output, Fq<mnt4753_pp> x) {
  fwrite((void *) x.mont_repr.data, libff::mnt4753_q_limbs * sizeof(mp_size_t), 1, output);
}

void write_mnt6_fq(FILE* output, Fq<mnt6753_pp> x) {
  fwrite((void *) x.mont_repr.data, libff::mnt6753_q_limbs * sizeof(mp_size_t), 1, output);
}

void write_mnt4_fq2(FILE* output, Fqe<mnt4753_pp> x) {
  write_mnt4_fq(output, x.c0);
  write_mnt4_fq(output, x.c1);
}

Fq<mnt4753_pp> read_mnt4_fq(FILE* input) {
  Fq<mnt4753_pp> x;
  fread((void *) x.mont_repr.data, libff::mnt4753_q_limbs * sizeof(mp_size_t), 1, input);
  return x;
}

Fq<mnt6753_pp> read_mnt6_fq(FILE* input) {
  Fq<mnt6753_pp> x;
  fread((void *) x.mont_repr.data, libff::mnt6753_q_limbs * sizeof(mp_size_t), 1, input);
  return x;
}

Fqe<mnt4753_pp> read_mnt4_fq2(FILE* input) {
  Fq<mnt4753_pp> c0 = read_mnt4_fq(input);
  Fq<mnt4753_pp> c1 = read_mnt4_fq(input);
  return Fqe<mnt4753_pp>(c0, c1);
}

void write_mnt6_fq3(FILE* output, Fqe<mnt6753_pp> x) {
  write_mnt6_fq(output, x.c0);
  write_mnt6_fq(output, x.c1);
  write_mnt6_fq(output, x.c2);
}

Fqe<mnt6753_pp> read_mnt6_fq3(FILE* input) {
  Fq<mnt6753_pp> c0 = read_mnt6_fq(input);
  Fq<mnt6753_pp> c1 = read_mnt6_fq(input);
  Fq<mnt6753_pp> c2 = read_mnt6_fq(input);
  return Fqe<mnt6753_pp>(c0, c1, c2);
}

G1<mnt4753_pp> read_mnt4_g1(FILE* input) {
  Fq<mnt4753_pp> x = read_mnt4_fq(input);
  Fq<mnt4753_pp> y = read_mnt4_fq(input);
  if (y == Fq<mnt4753_pp>::zero()) {
    return G1<mnt4753_pp>::zero();
  }
  return G1<mnt4753_pp>(x, y, Fq<mnt4753_pp>::one());
}

G1<mnt6753_pp> read_mnt6_g1(FILE* input) {
  Fq<mnt6753_pp> x = read_mnt6_fq(input);
  Fq<mnt6753_pp> y = read_mnt6_fq(input);
  if (y == Fq<mnt6753_pp>::zero()) {
    return G1<mnt6753_pp>::zero();
  }
  return G1<mnt6753_pp>(x, y, Fq<mnt6753_pp>::one());
}

G2<mnt4753_pp> read_mnt4_g2(FILE* input) {
  Fqe<mnt4753_pp> x = read_mnt4_fq2(input);
  Fqe<mnt4753_pp> y = read_mnt4_fq2(input);
  if (y == Fqe<mnt4753_pp>::zero()) {
    return G2<mnt4753_pp>::zero();
  }
  return G2<mnt4753_pp>(x, y, Fqe<mnt4753_pp>::one());
}

G2<mnt6753_pp> read_mnt6_g2(FILE* input) {
  Fqe<mnt6753_pp> x = read_mnt6_fq3(input);
  Fqe<mnt6753_pp> y = read_mnt6_fq3(input);
  if (y == Fqe<mnt6753_pp>::zero()) {
    return G2<mnt6753_pp>::zero();
  }
  return G2<mnt6753_pp>(x, y, Fqe<mnt6753_pp>::one());
}

void write_mnt4_g1(FILE* output, G1<mnt4753_pp> g) {
  g.to_affine_coordinates();
  write_mnt4_fq(output, g.X());
  write_mnt4_fq(output, g.Y());
}

void write_mnt6_g1(FILE* output, G1<mnt6753_pp> g) {
  g.to_affine_coordinates();
  write_mnt6_fq(output, g.X());
  write_mnt6_fq(output, g.Y());
}

void write_mnt4_g2(FILE* output, G2<mnt4753_pp> g) {
  g.to_affine_coordinates();
  write_mnt4_fq2(output, g.X());
  write_mnt4_fq2(output, g.Y());
}

void write_mnt6_g2(FILE* output, G2<mnt6753_pp> g) {
  g.to_affine_coordinates();
  write_mnt6_fq3(output, g.X());
  write_mnt6_fq3(output, g.Y());
}

Fr<mnt4753_pp> read_mnt4_fr(FILE* input) {
  Fr<mnt4753_pp> x;
  fread((void *) x.mont_repr.data, libff::mnt4753_r_limbs * sizeof(mp_size_t), 1, input);
  return x;
}

size_t read_size_t(FILE* input) {
  size_t n;
  fread((void *) &n, sizeof(size_t), 1, input);
  return n;
}

typedef mnt4753_pp ppT;
typedef Fr<ppT> F;

int main(int argc, const char * argv[])
{
    srand(time(NULL));
    setbuf(stdout, NULL);
    ppT::init_public_params();

    // B
    r1cs_gg_ppzksnark_verification_key<ppT> vk;
    std::ifstream vks;
    vks.open("vk");
    vks >> vk;
    vks.close();

    r1cs_gg_ppzksnark_proving_key<ppT> pk;

    std::ifstream parameters;
    parameters.open("parameters");
    parameters >> pk;

    auto params = fopen("params", "r");
    const size_t d = read_size_t(params);
    const size_t d_plus_1 = d + 1;
    const size_t m = read_size_t(params);

    const size_t primary_input_size = 1;

    std::vector<G1<ppT>> A;
    for (size_t i = 0; i <= m; ++i) {
      A.emplace_back(read_mnt4_g1(params));
    }

    std::vector<G1<ppT>> B1;
    for (size_t i = 0; i <= m; ++i) {
      B1.emplace_back(read_mnt4_g1(params));
    }

    std::vector<G2<ppT>> B2;
    for (size_t i = 0; i <= m; ++i) {
      B2.emplace_back(read_mnt4_g2(params));
    }

    std::vector<
      knowledge_commitment<libff::G2<ppT>, libff::G1<ppT> >> Bv;
    for (size_t i = 0; i <= m; ++i) {
      Bv.emplace_back(
        knowledge_commitment<libff::G2<ppT>, libff::G1<ppT> >(B2[i], B1[i]));
    }

    std::vector<G1<ppT>> L;
    for (size_t i = 0; i < m-1; ++i) {
      L.emplace_back(read_mnt4_g1(params));
    }

    std::vector<G1<ppT>> H;
    for (size_t i = 0; i < d; ++i) {
      H.emplace_back(read_mnt4_g1(params));
    }

    knowledge_commitment_vector<libff::G2<ppT>, libff::G1<ppT> > B(std::move(Bv));

    fclose(params);
    assert( m == pk.constraint_system.num_variables() );
    auto inputs = fopen("input", "r");
    std::vector<F> w(m+1);
    for (size_t i = 0; i < m+1; ++i) {
      w[i] = read_mnt4_fr(inputs);
    }

    std::vector<F> ca;
    for (size_t i = 0; i < d_plus_1; ++i) {
      ca.emplace_back(read_mnt4_fr(inputs));
    }
    std::vector<F> cb;
    for (size_t i = 0; i < d_plus_1; ++i) {
      cb.emplace_back(read_mnt4_fr(inputs));
    }
    std::vector<F> cc;
    for (size_t i = 0; i < d_plus_1; ++i) {
      cc.emplace_back(read_mnt4_fr(inputs));
    }

    fclose(inputs);
    std::vector<F> primary_input(
        w.begin() + 1,
        w.begin() + 1 + pk.constraint_system.num_inputs());
    std::vector<F> auxiliary_input(
        w.begin() + 1 + pk.constraint_system.num_inputs(),
        w.end() );

    for (size_t i = 0; i < primary_input.size(); ++i) {
      primary_input[i].print();
    }

    libff::enter_block("Call to r1cs_gg_ppzksnark_prover");

    libff::enter_block("Compute the polynomial H");
    // Begin witness map
    auto d1 = Fr<ppT>::zero();
    auto d2 = Fr<ppT>::zero();
    auto d3 = Fr<ppT>::zero();
    libff::enter_block("Call to r1cs_to_qap_witness_map");

    /* sanity check */

    r1cs_constraint_system<F> cs = pk.constraint_system;
    assert(cs.is_satisfied(primary_input, auxiliary_input));

    const std::shared_ptr<libfqfft::evaluation_domain<F> > domain = libfqfft::get_evaluation_domain<F>(cs.num_constraints() + cs.num_inputs() + 1);

    r1cs_variable_assignment<F> full_variable_assignment = primary_input;
    full_variable_assignment.insert(full_variable_assignment.end(), auxiliary_input.begin(), auxiliary_input.end());

    for (size_t i = 0; i < primary_input_size + auxiliary_input.size(); ++i) {
      if (full_variable_assignment[i] != w[1+i]) {
        printf("Bad!! %d\n", i);
        full_variable_assignment[i].print();
        w[i+1].print();
        assert(false);
      }
    }

    libff::enter_block("Compute coefficients of polynomial A");
    domain->iFFT(ca);
    libff::leave_block("Compute coefficients of polynomial A");

    libff::enter_block("Compute coefficients of polynomial B");
    domain->iFFT(cb);
    libff::leave_block("Compute coefficients of polynomial B");

    std::vector<F> coefficients_for_H(domain->m+1, F::zero());

    libff::enter_block("Compute evaluation of polynomial A on set T");
    domain->cosetFFT(ca, F::multiplicative_generator);
    libff::leave_block("Compute evaluation of polynomial A on set T");

    libff::enter_block("Compute evaluation of polynomial B on set T");
    domain->cosetFFT(cb, F::multiplicative_generator);
    libff::leave_block("Compute evaluation of polynomial B on set T");

    libff::enter_block("Compute evaluation of polynomial H on set T");
    std::vector<F> &H_tmp = ca; // can overwrite ca because it is not used later
#ifdef MULTICORE
#pragma omp parallel for
#endif
    for (size_t i = 0; i < domain->m; ++i)
    {
        H_tmp[i] = ca[i]*cb[i];
    }
    std::vector<F>().swap(cb); // destroy cb

    libff::enter_block("Compute coefficients of polynomial C");
    domain->iFFT(cc);
    libff::leave_block("Compute coefficients of polynomial C");

    libff::enter_block("Compute evaluation of polynomial C on set T");
    domain->cosetFFT(cc, F::multiplicative_generator);
    libff::leave_block("Compute evaluation of polynomial C on set T");

#ifdef MULTICORE
#pragma omp parallel for
#endif
    for (size_t i = 0; i < domain->m; ++i)
    {
        H_tmp[i] = (H_tmp[i]-cc[i]);
    }

    libff::enter_block("Divide by Z on set T");
    domain->divide_by_Z_on_coset(H_tmp);
    libff::leave_block("Divide by Z on set T");

    libff::leave_block("Compute evaluation of polynomial H on set T");

    libff::enter_block("Compute coefficients of polynomial H");
    domain->icosetFFT(H_tmp, F::multiplicative_generator);
    libff::leave_block("Compute coefficients of polynomial H");

    libff::enter_block("Compute sum of H and ZK-patch");
#ifdef MULTICORE
#pragma omp parallel for
#endif
    for (size_t i = 0; i < domain->m; ++i)
    {
        coefficients_for_H[i] = H_tmp[i];
    }
    libff::leave_block("Compute sum of H and ZK-patch");

    libff::leave_block("Call to r1cs_to_qap_witness_map");

    const qap_witness<libff::Fr<ppT> > qap_wit = qap_witness<F>(cs.num_variables(),
                               domain->m,
                               cs.num_inputs(),
                               d1,
                               d2,
                               d3,
                               full_variable_assignment,
                               std::move(coefficients_for_H));



    // End witness map

    /* We are dividing degree 2(d-1) polynomial by degree d polynomial
       and not adding a PGHR-style ZK-patch, so our H is degree d-2 */
    assert(!qap_wit.coefficients_for_H[d-1].is_zero());
    assert(qap_wit.coefficients_for_H[d].is_zero());
    assert(qap_wit.coefficients_for_H[d+1].is_zero());
    libff::leave_block("Compute the polynomial H");

    /* Choose two random field elements for prover zero-knowledge. */
    const libff::Fr<ppT> r = libff::Fr<ppT>::random_element();
    const libff::Fr<ppT> s = libff::Fr<ppT>::random_element();

#ifdef MULTICORE
    const size_t chunks = omp_get_max_threads(); // to override, set OMP_NUM_THREADS env var or call omp_set_num_threads()
#else
    const size_t chunks = 1;
#endif

    libff::enter_block("Compute the proof");

    libff::enter_block("Compute evaluation to A-query", false);
    // TODO: sort out indexing
    libff::Fr_vector<ppT> const_padded_assignment(1, libff::Fr<ppT>::one());
    const_padded_assignment.insert(const_padded_assignment.end(), qap_wit.coefficients_for_ABCs.begin(), qap_wit.coefficients_for_ABCs.end());

    libff::G1<ppT> evaluation_At = libff::multi_exp_with_mixed_addition<libff::G1<ppT>,
                                                                        libff::Fr<ppT>,
                                                                        method>(
        A.begin(),
        A.begin() + m + 1,
        const_padded_assignment.begin(),
        const_padded_assignment.begin() + m + 1,
        chunks);
    libff::leave_block("Compute evaluation to A-query", false);

    libff::enter_block("Compute evaluation to B-query", false);
    knowledge_commitment<libff::G2<ppT>, libff::G1<ppT> > evaluation_Bt = kc_multi_exp_with_mixed_addition<libff::G2<ppT>,
                                                                                                           libff::G1<ppT>,
                                                                                                           libff::Fr<ppT>,
                                                                                                           method>(
        B,
        0,
        m + 1,
        const_padded_assignment.begin(),
        const_padded_assignment.begin() + m + 1,
        chunks);
    libff::leave_block("Compute evaluation to B-query", false);

    libff::enter_block("Compute evaluation to H-query", false);
    libff::G1<ppT> evaluation_Ht = libff::multi_exp<libff::G1<ppT>,
                                                    libff::Fr<ppT>,
                                                    method>(
        H.begin(),
        H.begin() + d,
        qap_wit.coefficients_for_H.begin(),
        qap_wit.coefficients_for_H.begin() + d,
        chunks);
    libff::leave_block("Compute evaluation to H-query", false);

    libff::enter_block("Compute evaluation to L-query", false);
    libff::G1<ppT> evaluation_Lt = libff::multi_exp_with_mixed_addition<libff::G1<ppT>,
                                                                        libff::Fr<ppT>,
                                                                        method>(
        L.begin(),
        L.end(),
        const_padded_assignment.begin() + primary_input_size + 1,
        const_padded_assignment.begin() + m + 1,
        chunks);
    libff::leave_block("Compute evaluation to L-query", false);

    /* A = alpha + sum_i(a_i*A_i(t)) + r*delta */
    libff::G1<ppT> g1_A = pk.alpha_g1 + evaluation_At + r * pk.delta_g1;

    /* B = beta + sum_i(a_i*B_i(t)) + s*delta */
    libff::G1<ppT> g1_B = pk.beta_g1 + evaluation_Bt.h + s * pk.delta_g1;
    libff::G2<ppT> g2_B = pk.beta_g2 + evaluation_Bt.g + s * pk.delta_g2;

    /* C = sum_i(a_i*((beta*A_i(t) + alpha*B_i(t) + C_i(t)) + H(t)*Z(t))/delta) + A*s + r*b - r*s*delta */
    libff::G1<ppT> g1_C = evaluation_Ht + evaluation_Lt + s *  g1_A + r * g1_B - (r * s) * pk.delta_g1;

    libff::leave_block("Compute the proof");

    libff::leave_block("Call to r1cs_gg_ppzksnark_prover");

    r1cs_gg_ppzksnark_proof<ppT> proof = r1cs_gg_ppzksnark_proof<ppT>(std::move(g1_A), std::move(g2_B), std::move(g1_C));
    proof.print_size();

    assert (r1cs_gg_ppzksnark_verifier_strong_IC<ppT>(vk, primary_input, proof) );

    r1cs_gg_ppzksnark_proof<ppT> proof1=
      r1cs_gg_ppzksnark_prover<ppT>(
          pk, 
          primary_input,
          auxiliary_input);
    assert (r1cs_gg_ppzksnark_verifier_strong_IC<ppT>(vk, primary_input, proof1) );

    printf("YO\n");
    proof1.print_size();

    return 0;
}
