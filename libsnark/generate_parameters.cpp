#include <cassert>
#include <cstdio>
#include <fstream>

#include <libff/common/rng.hpp>
#include <libff/common/profiling.hpp>
#include <libff/common/utils.hpp>

#include <libff/algebra/curves/mnt753/mnt4753/mnt4753_pp.hpp>
#include <libff/algebra/curves/mnt753/mnt6753/mnt6753_pp.hpp>
#include <omp.h>
#include <libff/algebra/scalar_multiplication/multiexp.hpp>
#include <libsnark/knowledge_commitment/kc_multiexp.hpp>
#include <libsnark/reductions/r1cs_to_qap/r1cs_to_qap.hpp>
#include <libsnark/relations/constraint_satisfaction_problems/r1cs/examples/r1cs_examples.hpp>
#include <libsnark/zk_proof_systems/ppzksnark/r1cs_gg_ppzksnark/r1cs_gg_ppzksnark.hpp>

using namespace libsnark;
using namespace libff;

void write_size_t(FILE* output, size_t n) {
  fwrite((void *) &n, sizeof(size_t), 1, output);
}

void write_mnt4_fr(FILE* output, Fr<mnt4753_pp> x) {
  fwrite((void *) x.mont_repr.data, libff::mnt4753_r_limbs * sizeof(mp_size_t), 1, output);
}

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

void write_mnt6_fq3(FILE* output, Fqe<mnt6753_pp> x) {
  write_mnt6_fq(output, x.c0);
  write_mnt6_fq(output, x.c1);
  write_mnt6_fq(output, x.c2);
}

void write_mnt4_g1(FILE* output, G1<mnt4753_pp> g) {
  if (g.is_zero())  {
    write_mnt4_fq(output, Fq<mnt4753_pp>::zero());
    write_mnt4_fq(output, Fq<mnt4753_pp>::zero());
    return;
  }

  g.to_affine_coordinates();
  write_mnt4_fq(output, g.X());
  write_mnt4_fq(output, g.Y());
}

void write_mnt6_g1(FILE* output, G1<mnt6753_pp> g) {
  if (g.is_zero())  {
    write_mnt6_fq(output, Fq<mnt6753_pp>::zero());
    write_mnt6_fq(output, Fq<mnt6753_pp>::zero());
    return;
  }

  g.to_affine_coordinates();
  write_mnt6_fq(output, g.X());
  write_mnt6_fq(output, g.Y());
}

void write_mnt4_g2(FILE* output, G2<mnt4753_pp> g) {
  if (g.is_zero())  {
    write_mnt4_fq2(output, Fqe<mnt4753_pp>::zero());
    write_mnt4_fq2(output, Fqe<mnt4753_pp>::zero());
    return;
  }

  g.to_affine_coordinates();
  write_mnt4_fq2(output, g.X());
  write_mnt4_fq2(output, g.Y());
}

void write_mnt6_g2(FILE* output, G2<mnt6753_pp> g) {
  if (g.is_zero())  {
    write_mnt6_fq3(output, Fqe<mnt6753_pp>::zero());
    write_mnt6_fq3(output, Fqe<mnt6753_pp>::zero());
    return;
  }
  g.to_affine_coordinates();
  write_mnt6_fq3(output, g.X());
  write_mnt6_fq3(output, g.Y());
}

template<typename ppT>
std::vector<G1<ppT>> random_g1_vector(size_t n) {
  uint64_t offset = rand();

  std::vector<G1<ppT>> res(n, G1<ppT>::one());

#ifdef MULTICORE
#pragma omp parallel for
#endif
  for (size_t i = 0; i < n; ++i) {
    Fq<ppT> x = SHA512_rng<Fq<ppT>>(offset + i);
    Fq<ppT> y; 

    while (true) {
      Fq<ppT> y2 = x * x.squared() + G1<ppT>::coeff_a * x + G1<ppT>::coeff_b ;
      Fq<ppT> y2e = y2 ^ Fq<ppT>::euler;
      bool y2_is_square = y2e == Fq<ppT>::one();

      if (y2_is_square) {
        y = y2.sqrt();
        break;
      } else {
        x += 1;
      }
    }

    res[i] = G1<ppT>(x, y, Fq<ppT>::one());
  }

  return res;
}

template<typename ppT>
std::vector<G2<ppT>> random_g2_vector(size_t n) {
  std::vector<G2<ppT>> res(n, G2<ppT>::one());

  uint64_t offset0 = rand();

#ifdef MULTICORE
#pragma omp parallel for
#endif
  for (size_t i = 0; i < n; ++i) {
    Fqe<ppT> x;
    Fqe<ppT> y;

    uint64_t offset = offset0 + 128 * i;

    while (true) {
      Fq<ppT> x0 = SHA512_rng<Fq<ppT>>(offset);
      Fq<ppT> x1 = SHA512_rng<Fq<ppT>>(offset + 1);
      offset += 2;

      x = Fqe<ppT>(x0, x1);

      Fqe<ppT> y2 = x * x.squared() + G2<ppT>::coeff_a * x + G2<ppT>::coeff_b ;
      Fqe<ppT> y2e = y2 ^ Fqe<ppT>::euler;
      bool y2_is_square = y2e == Fqe<ppT>::one();

      if (y2_is_square) {
        y = y2.sqrt();
        break;
      }
    }

    res[i] = G2<ppT>(x, y, Fqe<ppT>::one());
  }

  return res;
}

typedef mnt4753_pp pp;
typedef mnt4753_pp ppT;
typedef Fr<pp> F;

int main(int argc, const char * argv[])
{
    srand(time(NULL));
    setbuf(stdout, NULL);

    pp::init_public_params();

    const size_t primary_input_size = 1;

    // auto output = fopen("parameters", "w");

    size_t d_plus_1 = 1 << 14;
    size_t d = d_plus_1 - 1;

    r1cs_example<F> example = generate_r1cs_example_with_field_input<F>(d-1, 1);
    r1cs_gg_ppzksnark_keypair<pp> keypair = r1cs_gg_ppzksnark_generator<pp>(example.constraint_system);

    // Debug
    std::ofstream ovk;
    ovk.open("vk");
    ovk << keypair.vk;
    ovk.close();

    for (size_t i = 0; i < example.primary_input.size(); ++i) {
      example.primary_input[i].print();
    }
    assert (example.primary_input.size() == 1);

    auto params = fopen("params", "w");
    size_t m = example.constraint_system.num_variables();
    write_size_t(params, d);
    write_size_t(params, m);

    for (size_t i = 0; i <= m; ++i) {
      write_mnt4_g1(params, keypair.pk.A_query[i]);
    }

    for (size_t i = 0; i <= m; ++i) {
      write_mnt4_g1(params, keypair.pk.B_query[i].h);
    }

    for (size_t i = 0; i <= m; ++i) {
      write_mnt4_g2(params, keypair.pk.B_query[i].g);
    }

    for (size_t i = 0; i < m-1; ++i) {
      write_mnt4_g1(params, keypair.pk.L_query[i]);
    }

    for (size_t i = 0; i < d; ++i) {
      write_mnt4_g1(params, keypair.pk.H_query[i]);
    }
    fclose(params);

    std::ofstream parameters;
    parameters.open("parameters");
    parameters << keypair.pk;
    parameters.close();
    printf("%d == %d + %d == %d + %d ? \n"
        , keypair.pk.constraint_system.num_variables()
        , keypair.pk.constraint_system.primary_input_size
        , keypair.pk.constraint_system.auxiliary_input_size
        );

    auto input = fopen("input", "w");

    write_mnt4_fr(input, Fr<pp>::one());
    for (size_t i = 0; i < example.primary_input.size(); ++i) {
      write_mnt4_fr(input, example.primary_input[i]);
    }
    for (size_t i = 0; i < example.auxiliary_input.size(); ++i) {
      write_mnt4_fr(input, example.auxiliary_input[i]);
    }

    r1cs_variable_assignment<F> full_variable_assignment = example.primary_input;
    full_variable_assignment.insert(full_variable_assignment.end(), example.auxiliary_input.begin(), example.auxiliary_input.end());

    std::vector<F> aA(d_plus_1, F::zero()), aB(d_plus_1, F::zero()), aC(d_plus_1, F::zero());
    for (size_t i = 0; i <= primary_input_size; ++i)
    {
        aA[i+keypair.pk.constraint_system.num_constraints()] = (i > 0 ? full_variable_assignment[i-1] : F::one());
    }
    for (size_t i = 0; i < keypair.pk.constraint_system.num_constraints(); ++i)
    {
        aA[i] += keypair.pk.constraint_system.constraints[i].a.evaluate(full_variable_assignment);
        aB[i] += keypair.pk.constraint_system.constraints[i].b.evaluate(full_variable_assignment);
    }
    for (size_t i = 0; i < keypair.pk.constraint_system.num_constraints(); ++i)
    {
        aC[i] += keypair.pk.constraint_system.constraints[i].c.evaluate(full_variable_assignment);
    }
    printf("hi my friend is\n");
    aB[0].print();
    aB[1].print();
    aB[2].print();

    for (size_t i = 0; i < d_plus_1; ++i) {
      write_mnt4_fr(input, aA[i]);
    }
    for (size_t i = 0; i < d_plus_1; ++i) {
      write_mnt4_fr(input, aB[i]);
    }
    for (size_t i = 0; i < d_plus_1; ++i) {
      write_mnt4_fr(input, aC[i]);
    }

    fclose(input);

    return 0;
}
