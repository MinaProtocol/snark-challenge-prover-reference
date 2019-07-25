#include <cassert>
#include <chrono>
#include <cstdio>
#include <fstream>
#include <string>

#include <libsnark/prover_reference_include/prover_reference_functions.hpp>

using namespace std::chrono; 
using namespace std;

static inline auto now() -> decltype(std::chrono::high_resolution_clock::now()) {
    return std::chrono::high_resolution_clock::now();
}


template<typename T>
void
print_time(T &t1, const char *str) {
    auto t2 = std::chrono::high_resolution_clock::now();
    auto tim = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
    printf("%s: %ld ms\n", str, tim);
    t1 = t2;
}


// Here is where all the FFTs happen.
template<typename B>
typename B::vector_Fr *compute_H(
    size_t d,
    typename B::vector_Fr *ca,
    typename B::vector_Fr *cb,
    typename B::vector_Fr *cc)
{
    // Begin witness map
    printf("Compute the polynomial H\n");

    Kernel k;

    auto domain = B::get_evaluation_domain(d + 1);

    //domain->iFFT(ca);
    //domain->iFFT(cb);
    B::domain_iFFT_GPU(domain, ca, k);
    B::domain_iFFT(domain, cb);

    B::domain_cosetFFT(domain, ca);
    B::domain_cosetFFT(domain, cb);

    printf("Compute evaluation of polynomial H on set T\n");
    // Use ca to store H
    auto H_tmp = ca;

    // for (size_t i = 0; i < domain->m; ++i)
    // {
    //     H_tmp[i] = ca[i]*cb[i];
    // }

    size_t m = B::domain_get_m(domain);
    B::vector_Fr_muleq(H_tmp, cb, m);

    B::domain_iFFT(domain, cc);
    B::domain_cosetFFT(domain, cc);
    
    m = B::domain_get_m(domain);

    // for (size_t i = 0; i < domain->m; ++i)
    // {
    //     H_tmp[i] = (H_tmp[i]-cc[i]);
    // }

    B::vector_Fr_subeq(H_tmp, cc, m);

    B::domain_divide_by_Z_on_coset(domain, H_tmp);

    printf("Compute evaluation of polynomial H on set T\n");

    B::domain_icosetFFT(domain, H_tmp);

    m = B::domain_get_m(domain);

    typename B::vector_Fr *H_res = B::vector_Fr_zeros(m + 1);

    // for (size_t i = 0; i < domain->m; ++i)
    // {
    //     coefficients_for_H[i] = H_tmp[i];
    // }

    B::vector_Fr_copy_into(H_tmp, H_res, m);

    printf("Compute the polynomial H\n");

    return H_res;
}

template<typename B>
int run_prover(
    const char* params_path,
    const char* input_path,
    const char* output_path)
{
    B::init_public_params();

    auto beginning = now();

    auto t = beginning;

    const size_t primary_input_size = 1;

    auto params = B::read_params(params_path);
    print_time(t, "load params");
    auto t_main = t;
    auto input = B::read_input(input_path, params);
    print_time(t, "load inputs");

    auto d = B::params_d(params);

    // End reading of parameters and input

    printf("Call to r1cs_gg_ppzksnark_prover\n");

    auto coefficients_for_H =
        compute_H<B>(B::params_d(params), B::input_ca(input), B::input_cb(input),
                     B::input_cc(input));

    printf("Compute the proof\n");
    printf("Multi-exponentiations\n");

    // Now the 5 multi-exponentiations
    printf("A G1 multiexp\n");
    typename B::G1 *evaluation_At = B::multiexp_G1(
        B::input_w(input), B::params_A(params), B::params_m(params) + 1);

    printf("B G1 multiexp\n");
    typename B::G1 *evaluation_Bt1 = B::multiexp_G1(
        B::input_w(input), B::params_B1(params), B::params_m(params) + 1);

    printf("B G2 multiexp\n");
    typename B::G2 *evaluation_Bt2 = B::multiexp_G2(
        B::input_w(input), B::params_B2(params), B::params_m(params) + 1);

    printf("H G1 multiexp\n");
    typename B::G1 *evaluation_Ht = B::multiexp_G1(
        coefficients_for_H, B::params_H(params), B::params_d(params));

    printf("L G1 multiexp\n");
    typename B::G1 *evaluation_Lt = B::multiexp_G1(
        B::vector_Fr_offset(B::input_w(input), primary_input_size + 1),
        B::params_L(params), B::params_m(params) - 1);

    auto scaled_Bt1 = B::G1_scale(B::input_r(input), evaluation_Bt1);
    auto Lt1_plus_scaled_Bt1 = B::G1_add(evaluation_Lt, scaled_Bt1);
    auto C = B::G1_add(evaluation_Ht, Lt1_plus_scaled_Bt1);

    print_time(t, "gpu");


    print_time(t, "store");

    B::groth16_output_write(evaluation_At, evaluation_Bt2, C, output_path);

    print_time(t_main, "Total time from input to output: ");

    // free everything
    B::delete_G1(evaluation_Bt1);
    B::delete_G1(evaluation_Ht);
    B::delete_G1(evaluation_Lt);
    B::delete_G1(scaled_Bt1);
    B::delete_G1(Lt1_plus_scaled_Bt1);
    B::delete_vector_Fr(coefficients_for_H);
    B::delete_groth16_input(input);
    B::delete_groth16_params(params);

    return 0;
}

int main(int argc, const char * argv[])
{
  setbuf(stdout, NULL);
  std::string curve(argv[1]);
  std::string mode(argv[2]);

  const char* params_path = argv[3];
  const char* input_path = argv[4];
  const char* output_path = argv[5];

  if (curve == "MNT4753") {
    if (mode == "compute") {
      return run_prover<mnt4753_libsnark>(params_path, input_path, output_path);
    }
  } else if (curve == "MNT6753") {
    if (mode == "compute") {
      return run_prover<mnt6753_libsnark>(params_path, input_path, output_path);
    }
  }
}
