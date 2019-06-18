#pragma once

#include <vector>

struct groth16_parameters_mnt4753;

struct evaluation_domain_mnt4753;

struct field_mnt4753;

struct group_mnt4753;

void domain_iFFT_mnt4753(std::vector<field_mnt4753> &a);

void domain_cosetFFT_mnt4753(std::vector<field_mnt4753> &a);

void domain_icosetFFT_mnt4753(std::vector<field_mnt4753> &a);

void domain_divide_by_Z_on_coset_mnt4753(std::vector<field_mnt4753> &a);

void multiexp_mnt4753(std::vector<field_mnt4753>::const_iterator scalar_start, std::vector<group_mnt4753>::const_iterator g_start, size_t length);

struct groth16_parameters_mnt6753;

struct evaluation_domain_mnt6753;

struct field_mnt6753;

struct group_mnt6753;

void domain_iFFT_mnt6753(std::vector<field_mnt6753> &a);

void domain_cosetFFT_mnt6753(std::vector<field_mnt6753> &a);

void domain_icosetFFT_mnt6753(std::vector<field_mnt6753> &a);

void domain_divide_by_Z_on_coset_mnt6753(std::vector<field_mnt6753> &a);

std::unique_ptr<group_mnt6753> multiexp_mnt6753(std::vector<field_mnt6753>::const_iterator scalar_start, std::vector<group_mnt6753>::const_iterator g_start, size_t length);
