#include <prover_reference_functions.hpp>

// Overwrites 
vector_Fr_mnt4753 compute_H_mnt4753(size_t d, 
        vector_Fr_mnt4753 *ca,
        vector_Fr_mnt4753 *cb,
        vector_Fr_mnt4753 *cc) {
    auto domain = get_evaluation_domain_mnt4753(d + 1);

    domain_iFFT_mnt4753(domain, ca);
    domain_iFFT_mnt4753(domain, cb);

    domain_cosetFFT_mnt4753(domain, ca);
    domain_cosetFFT_mnt4753(domain, cb);

    // Use ca to store H
    auto H_tmp = ca;

    size_t m = domain_get_m_mnt4753(domain);
    for (size_t i = 0; i < m; i++) {
        // H_tmp[i] *= cb[i]
        vector_Fr_mnt4753_set_mul(H_tmp, cb, i);
    }
    delete_vector_Fr_mnt4753(cb);

    domain_iFFT_mnt4753(cc);
    domain_cosetFFT_mnt4753(cc);

    size_t m = domain_get_m_mnt4753(domain);
    for (size_t i = 0; i < m; i++) {
        // H_tmp[i] -= cc[i]
        vector_Fr_mnt4753_set_sub(H_tmp, cc, i);
    }

    domain_divide_by_Z_on_coset_mnt4753(H_tmp);

    domain_icosetFFT_mnt4753(H_tmp);
    
    size_t m = domain_get_m_mnt4753(domain);
    auto coefficients_for_H = vector_Fr_mnt4753_zeros(m+1);
}
int main(int argc, char **argv) {
    return 0;
}

