#pragma once

#include <vector>

class mnt4753_libsnark {
public:
    struct groth16_input;

    struct groth16_params;

    struct groth16_output;

    struct evaluation_domain;

    struct field;

    struct G1;

    struct G2;

    struct vector_Fr;

    struct vector_G1;

    struct vector_G2;

    static void init_public_params();

    static evaluation_domain *get_evaluation_domain(size_t d);

    static G1 *G1_add(G1 *a, G1 *b);
    static G1 *G1_scale(field *a, G1 *b);

    static void vector_Fr_muleq(vector_Fr *a, vector_Fr *b, size_t size);
    static void vector_Fr_subeq(vector_Fr *a, vector_Fr *b, size_t size);
    static vector_Fr *vector_Fr_offset(vector_Fr *a, size_t offset);
    static vector_Fr *vector_Fr_copy(vector_Fr *a, size_t length);

    static void domain_iFFT(evaluation_domain* domain, vector_Fr *a);
    static void domain_cosetFFT(evaluation_domain* domain, vector_Fr *a);
    static void domain_icosetFFT(evaluation_domain* domain, vector_Fr *a);
    static void domain_divide_by_Z_on_coset(evaluation_domain* domain, vector_Fr *a);
    static size_t domain_get_m(evaluation_domain *domain);

    static G1 *multiexp_G1(vector_Fr *scalar_start, vector_G1 *g_start, size_t length);
    static G2 *multiexp_G2(vector_Fr *scalar_start, vector_G2 *g_start, size_t length);


    static groth16_input *read_input(const char *path, groth16_params *params);

    static vector_Fr *input_w(groth16_input* input);
    static vector_Fr *input_ca(groth16_input* input);
    static vector_Fr *input_cb(groth16_input* input);
    static vector_Fr *input_cc(groth16_input* input);
    static field *input_r(groth16_input* input);

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
    static void delete_groth16_output(groth16_output *a);
    static void delete_evaluation_domain(evaluation_domain *a);

    static groth16_output *groth16_output_create(G1 *At, G2 *Bt2, G1 *C);
    static void groth16_output_write(groth16_output *output, const char *output_path);
};
class mnt6753_libsnark {
public:
    struct groth16_input;

    struct groth16_params;

    struct groth16_output;

    struct evaluation_domain;

    struct field;

    struct G1;

    struct G2;

    struct vector_Fr;

    struct vector_G1;

    struct vector_G2;

    static void init_public_params();

    static evaluation_domain *get_evaluation_domain(size_t d);

    static G1 *G1_add(G1 *a, G1 *b);
    static G1 *G1_scale(field *a, G1 *b);

    static void vector_Fr_muleq(vector_Fr *a, vector_Fr *b, size_t size);
    static void vector_Fr_subeq(vector_Fr *a, vector_Fr *b, size_t size);
    static vector_Fr *vector_Fr_offset(vector_Fr *a, size_t offset);
    static vector_Fr *vector_Fr_copy(vector_Fr *a, size_t length);

    static void domain_iFFT(evaluation_domain* domain, vector_Fr *a);
    static void domain_cosetFFT(evaluation_domain* domain, vector_Fr *a);
    static void domain_icosetFFT(evaluation_domain* domain, vector_Fr *a);
    static void domain_divide_by_Z_on_coset(evaluation_domain* domain, vector_Fr *a);
    static size_t domain_get_m(evaluation_domain *domain);

    static G1 *multiexp_G1(vector_Fr *scalar_start, vector_G1 *g_start, size_t length);
    static G2 *multiexp_G2(vector_Fr *scalar_start, vector_G2 *g_start, size_t length);


    static groth16_input *read_input(const char *path, groth16_params *params);

    static vector_Fr *input_w(groth16_input* input);
    static vector_Fr *input_ca(groth16_input* input);
    static vector_Fr *input_cb(groth16_input* input);
    static vector_Fr *input_cc(groth16_input* input);
    static field *input_r(groth16_input* input);

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
    static void delete_groth16_output(groth16_output *a);
    static void delete_evaluation_domain(evaluation_domain *a);

    static groth16_output *groth16_output_create(G1 *At, G2 *Bt2, G1 *C);
    static void groth16_output_write(groth16_output *output, const char *output_path);
};
