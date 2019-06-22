#include <string>

#include <prover_reference_functions.hpp>

// cuda-fixnum includes

#include <fixnum/warp_fixnum.cu>
#include <array/fixnum_array.h>
#include <functions/modexp.cu>
#include <functions/multi_modexp.cu>
#include <modnum/modnum_monty_redc.cu>
#include <modnum/modnum_monty_cios.cu>

// cuda-fixnum/main.cu has some example code. this declaration lets us run it.
int do_fixnum_example(const char *inputs_file, const char *outputs_file);

// This is where all the FFTs happen

// template over the bundle of types and functions.
// Overwrites ca!
template<typename B>
typename B::vector_Fr *compute_H(size_t d, 
        typename B::vector_Fr *ca,
        typename B::vector_Fr *cb,
        typename B::vector_Fr *cc) {
    auto domain = B::get_evaluation_domain(d + 1);

    B::domain_iFFT(domain, ca);
    B::domain_iFFT(domain, cb);

    B::domain_cosetFFT(domain, ca);
    B::domain_cosetFFT(domain, cb);

    // Use ca to store H
    auto H_tmp = ca;

    size_t m = B::domain_get_m(domain);
    // for i in 0 to m: H_tmp[i] *= cb[i]
    B::vector_Fr_muleq(H_tmp, cb, m);

    B::domain_iFFT(domain, cc);
    B::domain_cosetFFT(domain, cc);

    m = B::domain_get_m(domain);

    // for i in 0 to m: H_tmp[i] -= cc[i]
    B::vector_Fr_subeq(H_tmp, cc, m);

    B::domain_divide_by_Z_on_coset(domain, H_tmp);

    B::domain_icosetFFT(domain, H_tmp);
    
    m = B::domain_get_m(domain);
    typename B::vector_Fr *H_res = B::vector_Fr_zeros(m+1);
    B::vector_Fr_copy_into(H_tmp, H_res, m);
    return H_res;
}

template<typename B>
void run_prover(const char *params_path, const char *input_path, const char *output_path) {
    B::init_public_params();

    size_t primary_input_size = 1;

    auto params = B::read_params(params_path);
    auto input = B::read_input(input_path, params);

    auto d = B::params_d(params);
    auto coefficients_for_H = compute_H<B>(B::params_d(params), B::input_ca(input), B::input_cb(input), B::input_cc(input));

    // Now the 5 multi-exponentiations
    typename B::G1 *evaluation_At = B::multiexp_G1(B::input_w(input), B::params_A(params), B::params_m(params) + 1);
    B::print_G1(evaluation_At);
    typename B::G1 *evaluation_Bt1 = B::multiexp_G1(B::input_w(input), B::params_B1(params), B::params_m(params) + 1);
    B::print_G1(evaluation_Bt1);
    typename B::G2 *evaluation_Bt2 = B::multiexp_G2(B::input_w(input), B::params_B2(params), B::params_m(params) + 1);
    B::print_G2(evaluation_Bt2);
    typename B::G1 *evaluation_Ht = B::multiexp_G1(coefficients_for_H, B::params_H(params), B::params_d(params));
    B::print_G1(evaluation_Ht);
    typename B::G1 *evaluation_Lt = B::multiexp_G1(B::vector_Fr_offset(B::input_w(input), primary_input_size + 1), B::params_L(params), B::params_m(params) - 1);
    B::print_G1(evaluation_Lt);


    auto scaled_Bt1 = B::G1_scale(B::input_r(input), evaluation_Bt1);
    auto Lt1_plus_scaled_Bt1 = B::G1_add(evaluation_Lt, scaled_Bt1);
    auto C = B::G1_add(evaluation_Ht, Lt1_plus_scaled_Bt1);
    B::print_G1(C);

    B::groth16_output_write(evaluation_At, evaluation_Bt2, C, output_path);
    
    // free everything
    B::delete_G1(evaluation_Bt1);
    B::delete_G1(evaluation_Ht);
    B::delete_G1(evaluation_Lt);
    B::delete_G1(scaled_Bt1);
    B::delete_G1(Lt1_plus_scaled_Bt1);
    B::delete_vector_Fr(coefficients_for_H);
    B::delete_groth16_input(input);
    B::delete_groth16_params(params);
}

int main(int argc, char **argv) {
  setbuf(stdout, NULL);
  std::string curve(argv[1]);
  std::string mode(argv[2]);
  std::string device(argv[3]); // this "device" arg 

    const char* params_path = argv[4];
    const char* input_path = argv[5];
    const char* output_path = argv[6];
  if (device == "CPU") {
    if (curve == "MNT4753") {
        if (mode == "compute") {
            run_prover<mnt4753_libsnark>(params_path, input_path, output_path);
        }
    } else if (curve == "MNT6753") {
        if (mode == "compute") {
            run_prover<mnt6753_libsnark>(params_path, input_path, output_path);
        }
    }
  } else if (device == "GPU") {
    do_fixnum_example(argv[4], argv[5]);
  }

  return 0;
}

