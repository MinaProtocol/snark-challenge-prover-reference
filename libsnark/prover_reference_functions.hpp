#pragma once

#include <vector>

struct groth16_input_mnt4753;

struct evaluation_domain_mnt4753;

struct field_mnt4753;

struct G1_mnt4753;

struct G2_mnt4753;

struct vector_Fr_mnt4753;

struct vector_G1_mnt4753;

struct vector_G2_mnt4753;

void domain_iFFT_mnt4753(evaluation_domain_mnt4753* domain, vector_Fr_mnt4753 *a);

void domain_cosetFFT_mnt4753(evaluation_domain_mnt4753* domain, vector_Fr_mnt4753 *a);

void domain_icosetFFT_mnt4753(evaluation_domain_mnt4753* domain, vector_Fr_mnt4753 *a);

void domain_divide_by_Z_on_coset_mnt4753(evaluation_domain_mnt4753* domain, vector_Fr_mnt4753 *a);

void multiexp_G1_mnt4753(vector_Fr_mnt4753 *scalar_start, vector_G1_mnt4753 *g_start, size_t length);

void multiexp_G2_mnt4753(vector_Fr_mnt4753 *scalar_start, vector_G2_mnt4753 *g_start, size_t length);

vector_Fr_mnt4753 *input_w(groth16_input_mnt4753* params);
vector_Fr_mnt4753 *input_ca(groth16_input_mnt4753* params);
vector_Fr_mnt4753 *input_cb(groth16_input_mnt4753* params);
vector_Fr_mnt4753 *input_cc(groth16_input_mnt4753* params);

field_mnt4753 *input_r(groth16_input_mnt4753* params);

groth16_input_mnt4753 *read_input(const char *path);

vector_G1

void delete_group_mnt4753(group_mnt4753 *a);
void delete_vector_Fr_mnt4753(vector_Fr_mnt4753 *a);
void delete_vector_G1_mnt4753(vector_G1_mnt4753 *a);
void delete_vector_G2_mnt4753(vector_G2_mnt4753 *a);
void delete_groth16_inputs_mnt4753(groth16_input_mnt4753 *a);
void delete_T(T *a);

