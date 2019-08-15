#include <cassert>
#include <cstdio>
#include <fstream>
#include <libff/algebra/curves/mnt753/mnt4753/mnt4753_pp.hpp>
#include <libff/algebra/curves/mnt753/mnt6753/mnt6753_pp.hpp>
#include <libff/algebra/scalar_multiplication/multiexp.hpp>
#include <libff/common/profiling.hpp>
#include <libff/common/rng.hpp>
#include <libff/common/utils.hpp>
#include <libsnark/knowledge_commitment/kc_multiexp.hpp>
#include <libsnark/knowledge_commitment/knowledge_commitment.hpp>
#include <libsnark/reductions/r1cs_to_qap/r1cs_to_qap.hpp>
#include <libsnark/serialization.hpp>
#include <omp.h>

#include <libsnark/zk_proof_systems/ppzksnark/r1cs_gg_ppzksnark/r1cs_gg_ppzksnark.hpp>

#include <libfqfft/evaluation_domain/domains/basic_radix2_domain.hpp>
//#include "libfqfft/evaluation_domain/evaluation_domain.hpp"
#include "prover_reference_include/prover_reference_functions.hpp"

using namespace libff;
using namespace libsnark;

const multi_exp_method method = multi_exp_method_BDLO12;

//#include "ocl_kernels/kernel.cpp"
// #define __CL_ENABLE_EXCEPTIONS
#include <CL/cl.h>
#define DATA_SIZE (131072)
#define limbs_per_elem (12)
#include <chrono> 

using namespace std::chrono; 
using namespace std;
using namespace libff;





template <typename G, typename Fr>
G multiexp(typename std::vector<Fr>::const_iterator scalar_start,
           typename std::vector<G>::const_iterator g_start, size_t length) {
#ifdef MULTICORE
  const size_t chunks =
      omp_get_max_threads(); // to override, set OMP_NUM_THREADS env var or call
                             // omp_set_num_threads()
#else
  const size_t chunks = 1;
#endif

  return libff::multi_exp_with_mixed_addition<G, Fr, method>(
      g_start, g_start + length, scalar_start, scalar_start + length, chunks);
}

class mnt4753_libsnark::groth16_input {
public:
  std::shared_ptr<std::vector<Fr<mnt4753_pp>>> w;
  std::shared_ptr<std::vector<Fr<mnt4753_pp>>> ca, cb, cc;
  Fr<mnt4753_pp> r;

  groth16_input(const char *path, size_t d, size_t m) {
    w = std::make_shared<std::vector<libff::Fr<mnt4753_pp>>>(
        std::vector<libff::Fr<mnt4753_pp>>());
    ca = std::make_shared<std::vector<libff::Fr<mnt4753_pp>>>(
        std::vector<libff::Fr<mnt4753_pp>>());
    cb = std::make_shared<std::vector<libff::Fr<mnt4753_pp>>>(
        std::vector<libff::Fr<mnt4753_pp>>());
    cc = std::make_shared<std::vector<libff::Fr<mnt4753_pp>>>(
        std::vector<libff::Fr<mnt4753_pp>>());
    FILE *inputs = fopen(path, "r");

    for (size_t i = 0; i < m + 1; ++i) {
      w->emplace_back(read_fr<mnt4753_pp>(inputs));
    }

    for (size_t i = 0; i < d + 1; ++i) {
      ca->emplace_back(read_fr<mnt4753_pp>(inputs));
    }
    for (size_t i = 0; i < d + 1; ++i) {
      cb->emplace_back(read_fr<mnt4753_pp>(inputs));
    }
    for (size_t i = 0; i < d + 1; ++i) {
      cc->emplace_back(read_fr<mnt4753_pp>(inputs));
    }

    r = read_fr<mnt4753_pp>(inputs);

    fclose(inputs);
  }
};

class mnt4753_libsnark::groth16_params {
public:
  size_t d;
  size_t m;
  std::shared_ptr<std::vector<libff::G1<mnt4753_pp>>> A, B1, L, H;
  std::shared_ptr<std::vector<libff::G2<mnt4753_pp>>> B2;

  groth16_params(const char *path) {
    FILE *params = fopen(path, "r");
    d = read_size_t(params);
    m = read_size_t(params);
    A = std::make_shared<std::vector<libff::G1<mnt4753_pp>>>(
        std::vector<libff::G1<mnt4753_pp>>());
    B1 = std::make_shared<std::vector<libff::G1<mnt4753_pp>>>(
        std::vector<libff::G1<mnt4753_pp>>());
    L = std::make_shared<std::vector<libff::G1<mnt4753_pp>>>(
        std::vector<libff::G1<mnt4753_pp>>());
    H = std::make_shared<std::vector<libff::G1<mnt4753_pp>>>(
        std::vector<libff::G1<mnt4753_pp>>());
    B2 = std::make_shared<std::vector<libff::G2<mnt4753_pp>>>(
        std::vector<libff::G2<mnt4753_pp>>());
    for (size_t i = 0; i <= m; ++i) {
      A->emplace_back(read_g1<mnt4753_pp>(params));
    }
    for (size_t i = 0; i <= m; ++i) {
      B1->emplace_back(read_g1<mnt4753_pp>(params));
    }
    for (size_t i = 0; i <= m; ++i) {
      B2->emplace_back(read_g2<mnt4753_pp>(params));
    }
    for (size_t i = 0; i < m - 1; ++i) {
      L->emplace_back(read_g1<mnt4753_pp>(params));
    }
    for (size_t i = 0; i < d; ++i) {
      H->emplace_back(read_g1<mnt4753_pp>(params));
    }
    fclose(params);
  }
};

struct mnt4753_libsnark::evaluation_domain {
  std::shared_ptr<libfqfft::evaluation_domain<Fr<mnt4753_pp>>> data;
};

struct mnt4753_libsnark::field {
  Fr<mnt4753_pp> data;
};

struct mnt4753_libsnark::G1 {
  libff::G1<mnt4753_pp> data;
};

struct mnt4753_libsnark::G2 {
  libff::G2<mnt4753_pp> data;
};

struct mnt4753_libsnark::vector_Fr {
  std::shared_ptr<std::vector<Fr<mnt4753_pp>>> data;
  size_t offset;
};

struct mnt4753_libsnark::vector_G1 {
  std::shared_ptr<std::vector<libff::G1<mnt4753_pp>>> data;
};
struct mnt4753_libsnark::vector_G2 {
  std::shared_ptr<std::vector<libff::G2<mnt4753_pp>>> data;
};

void mnt4753_libsnark::init_public_params() {
  mnt4753_pp::init_public_params();
}

void mnt4753_libsnark::print_G1(mnt4753_libsnark::G1 *a) { a->data.print(); }

void mnt4753_libsnark::print_G2(mnt4753_libsnark::G2 *a) { a->data.print(); }

mnt4753_libsnark::evaluation_domain *
mnt4753_libsnark::get_evaluation_domain(size_t d) {
  return new evaluation_domain{
      .data = libfqfft::get_evaluation_domain<Fr<mnt4753_pp>>(d)};
}

mnt4753_libsnark::G1 *mnt4753_libsnark::G1_add(mnt4753_libsnark::G1 *a,
                                               mnt4753_libsnark::G1 *b) {
  return new mnt4753_libsnark::G1{.data = a->data + b->data};
}

mnt4753_libsnark::G1 *mnt4753_libsnark::G1_scale(field *a, G1 *b) {
  return new G1{.data = a->data * b->data};
}

void mnt4753_libsnark::vector_Fr_muleq(mnt4753_libsnark::vector_Fr *a,
                                       mnt4753_libsnark::vector_Fr *b,
                                       size_t size) {
  size_t a_off = a->offset, b_off = b->offset;
#ifdef MULTICORE
#pragma omp parallel for
#endif
  for (size_t i = 0; i < size; i++) {
    a->data->at(i + a_off) = a->data->at(i + a_off) * b->data->at(i + b_off);
  }
}

void mnt4753_libsnark::vector_Fr_subeq(mnt4753_libsnark::vector_Fr *a,
                                       mnt4753_libsnark::vector_Fr *b,
                                       size_t size) {
  size_t a_off = a->offset, b_off = b->offset;
#ifdef MULTICORE
#pragma omp parallel for
#endif
  for (size_t i = 0; i < size; i++) {
    a->data->at(i + a_off) = a->data->at(i + a_off) - b->data->at(i + b_off);
  }
}

mnt4753_libsnark::vector_Fr *
mnt4753_libsnark::vector_Fr_offset(mnt4753_libsnark::vector_Fr *a,
                                   size_t offset) {
  return new vector_Fr{.data = a->data, .offset = offset};
}

void mnt4753_libsnark::vector_Fr_copy_into(mnt4753_libsnark::vector_Fr *src,
                                           mnt4753_libsnark::vector_Fr *dst,
                                           size_t length) {
  std::cerr << "length is " << length << ", offset is " << src->offset
            << ", size of src is " << src->data->size() << ", size of dst is "
            << dst->data->size() << std::endl;
#ifdef MULTICORE
#pragma omp parallel for
#endif
  for (size_t i = 0; i < length; i++) {
    // std::cerr << "doing iteration " << i << std::endl;
    dst->data->at(i) = src->data->at(i);
  }
  // std::copy(src->data->begin(), src->data->end(), dst->data->begin() );
}

mnt4753_libsnark::vector_Fr *mnt4753_libsnark::vector_Fr_zeros(size_t length) {
  std::vector<Fr<mnt4753_pp>> data(length, Fr<mnt4753_pp>::zero());
  return new mnt4753_libsnark::vector_Fr{
      .data = std::make_shared<std::vector<Fr<mnt4753_pp>>>(data)};
}

void mnt4753_libsnark::domain_iFFT(mnt4753_libsnark::evaluation_domain *domain,
                                   mnt4753_libsnark::vector_Fr *a) {
  std::vector<Fr<mnt4753_pp>> &data = *a->data;
  //printf("===== FFT CPU ===== \n");
  //data[0].print();
 // a->data->at(0).print();
  domain->data->iFFT(data);
  //printf("after FFT\n");
  //Fr<mnt4753_pp> elem = domain->data->get_domain_element(1);

  //a->data->at(0).print();
  //data[0].print();
}

void mnt4753_libsnark::domain_iFFT_GPU(mnt4753_libsnark::evaluation_domain *domain,
                                   mnt4753_libsnark::vector_Fr *a, Kernel kern) {
  std::vector<Fr<mnt4753_pp>> &data = *a->data;
  
  unsigned int MAX_RADIX_DEGREE = 8;
  unsigned int MAX_LOCAL_WORK_SIZE_DEGREE = 7;

  printf("===== FFT GPU ===== \n");
  data[0].print();
  printf("a size:%u \n", data.size());

  size_t m = data.size();
  Fr<mnt4753_pp> omega = Fr<mnt4753_pp>::zero();
  bool err = true;
  printf("gpu omega one:\n");
  Fr<mnt4753_pp>::one().mont_repr.print();

  // get omega for basic radix
  if (m <= 1) {
    err = true;
    omega = Fr<mnt4753_pp>::one();
  } else if (!std::is_same<Fr<mnt4753_pp>, libff::Double>::value) {
      if (Fr<mnt4753_pp>::small_subgroup_defined)
      {
          const size_t q = Fr<mnt4753_pp>::small_subgroup_base;

          const size_t q_adicity = libff::k_adicity(q, m);
          const size_t q_part = libff::pow_int(q, q_adicity);

          const size_t two_adicity = libff::k_adicity(2, m);
          const size_t two_part = 1u << two_adicity;

          if (m != q_part * two_part) {
            err = true;
            omega = Fr<mnt4753_pp>(1);
          } else {
            omega = libff::get_root_of_unity<Fr<mnt4753_pp>>(m, err);
          }
      }
      else
      {
        const size_t logm = libff::log2(m);
        if (logm > (Fr<mnt4753_pp>::s)) {
          err = true;
          omega = Fr<mnt4753_pp>(1);
        } else {
          omega = libff::get_root_of_unity<Fr<mnt4753_pp>>(m, err);
        }
      }
  }

  printf("omega after:\n");
  omega.mont_repr.print();


  cl_kernel kernel;                   // compute kernel
  cl_event event;                     // timing
  cl_ulong time_start;
  cl_ulong time_end;
  unsigned int n = data.size();
  unsigned int lgn = log2(n);
  unsigned int max_deg = std::min(MAX_RADIX_DEGREE, lgn);
  Fr<mnt4753_pp>* src_buff = new Fr<mnt4753_pp>[n];
  Fr<mnt4753_pp>* dst_buff = new Fr<mnt4753_pp>[n];
  Fr<mnt4753_pp>* pq_buff = new Fr<mnt4753_pp>[1 << MAX_RADIX_DEGREE >> 1];
  Fr<mnt4753_pp>* om_buff = new Fr<mnt4753_pp>[32]; 
  Fr<mnt4753_pp>* results = new Fr<mnt4753_pp>[n];
  unsigned int correct;               // number of correct results returned

  cl_mem src;
  cl_mem dst;
  cl_mem pq;
  cl_mem om;
  // Fill our data set with field inputs from param gen
  //
  unsigned int count = n;
  mp_size_t num = 1;

  for(int i = 0; i < count; i++) {
    memcpy(&src_buff[i], &data[i], sizeof(Fr<mnt4753_pp>));
  }
  
  // setup_pq
  Fr<mnt4753_pp> tw = omega^(n >> max_deg);
  pq_buff[0] = Fr<mnt4753_pp>::one();
  if(max_deg > 1) {
      pq_buff[1] = tw;
      for(int i=2; i < (1 << max_deg >> 1); i++) {
          pq_buff[i] = pq_buff[i-1];
          pq_buff[i] = pq_buff[i]*tw;
      }
  }
  om_buff[0] = omega;
  for(int i=1; i<32; i++) {
    om_buff[i] = om_buff[i-1]^2;
  }
  
  // Create the compute kernel in the program we wish to run
  //
  kernel = clCreateKernel(kern.program, "mnt4753_fft", &kern.err);
  if (!kernel || kern.err != CL_SUCCESS)
  {
      printf("Error: Failed to create compute kernel!\n");
      exit(1);
  }

  // Create the input and output arrays in device memory for our calculation
  //
  printf("creating buffer\n");
  src = clCreateBuffer(kern.context,  CL_MEM_READ_WRITE,  sizeof(Fr<mnt4753_pp>) * count, NULL, NULL);
  dst = clCreateBuffer(kern.context,  CL_MEM_READ_WRITE,  sizeof(Fr<mnt4753_pp>) * count, NULL, NULL);
  pq = clCreateBuffer(kern.context, CL_MEM_WRITE_ONLY, sizeof(Fr<mnt4753_pp>) * (1 << MAX_RADIX_DEGREE >> 1), NULL, NULL);
  om = clCreateBuffer(kern.context,  CL_MEM_READ_ONLY,  sizeof(Fr<mnt4753_pp>) * 32, NULL, NULL);

  if (!src || !dst)
  {
      printf("Error: Failed to allocate device memory!\n");
      exit(1);
  }

  // Write our data set into the input array in device memory 
  //
  auto start = high_resolution_clock::now();
  // write PQ
  kern.err = clEnqueueWriteBuffer(kern.commands, pq, CL_TRUE, 0, sizeof(Fr<mnt4753_pp>) *  (1 << MAX_RADIX_DEGREE >> 1), pq_buff, 0, NULL, NULL);
  if (kern.err != CL_SUCCESS)
  {
      printf("Error: Failed to write to pq source array!\n");
      exit(1);
  }

  kern.err = clEnqueueWriteBuffer(kern.commands, src, CL_TRUE, 0, sizeof(Fr<mnt4753_pp>) * count, src_buff, 0, NULL, NULL);
  if (kern.err != CL_SUCCESS)
  {
      printf("Error: Failed to write to source array!\n");
      exit(1);
  }
  kern.err = clEnqueueWriteBuffer(kern.commands, om, CL_TRUE, 0, sizeof(Fr<mnt4753_pp>) * 32, om_buff, 0, NULL, NULL);
  if (kern.err != CL_SUCCESS)
  {
      printf("Error: Failed to write to omega source array!\n");
      exit(1);
  }
  auto stop = high_resolution_clock::now();
  auto duration = duration_cast<microseconds>(stop - start); 
  cout << "Time taken by GPU write function: "
    << duration.count() << " microseconds" << endl;

  unsigned int lgp = 0;
  bool in_src = true;
  while(lgp < lgn) {
      unsigned int deg = std::min(max_deg, lgn-lgp);
      unsigned int lwsd = std::min(deg-1, MAX_LOCAL_WORK_SIZE_DEGREE);

      // Set the arguments to our compute kernel
      //
      kern.err = 0;
      kern.err |= clSetKernelArg(kernel, 0, sizeof(cl_mem), in_src ? &src : &dst);
      kern.err |= clSetKernelArg(kernel, 1, sizeof(cl_mem), in_src ? &dst : &src);
      kern.err |= clSetKernelArg(kernel, 2, sizeof(cl_mem), &pq);
      kern.err |= clSetKernelArg(kernel, 3, sizeof(cl_mem), &om);
      kern.err |= clSetKernelArg(kernel, 4, sizeof(Fr<mnt4753_pp>) * (1 << deg), NULL);
      kern.err |= clSetKernelArg(kernel, 5, sizeof(unsigned int), &count);
      kern.err |= clSetKernelArg(kernel, 6, sizeof(unsigned int), &lgp);
      kern.err |= clSetKernelArg(kernel, 7, sizeof(unsigned int), &deg);
      kern.err |= clSetKernelArg(kernel, 8, sizeof(unsigned int), &max_deg);
      if (kern.err != CL_SUCCESS)
      {
          printf("Error: Failed to set kernel arguments! %d\n", kern.err);
          exit(1);
      }

      // Get the maximum work group size for executing the kernel on the device
      //
      kern.err = clGetKernelWorkGroupInfo(kernel, kern.devices[0], CL_KERNEL_WORK_GROUP_SIZE, sizeof(kern.local), &kern.local, NULL);
      if (kern.err != CL_SUCCESS)
      {
          printf("Error: Failed to retrieve kernel work group info! %d\n", kern.err);
          exit(1);
      }

      printf("Max work size: %u\n", kern.local);

      // Execute the kernel over the entire range of our 1d input data set
      // using the maximum number of work group items for this device
      //
      kern.global = (n >> deg << lwsd);
      kern.local = (1 << lwsd);
      printf("queueing kernel\n");
      kern.err = clEnqueueNDRangeKernel(kern.commands, kernel, 1, NULL, &kern.global, &kern.local, 0, NULL, &event);
      if (kern.err)
      {
          printf("Error: Failed to execute kernel!\n");
          exit(1);
      }

      clWaitForEvents(1, &event);
      clFinish(kern.commands);

      // Time kernel execution time without read/write
      //
      clGetEventProfilingInfo(event, CL_PROFILING_COMMAND_START, sizeof(time_start), &time_start, NULL);
      clGetEventProfilingInfo(event, CL_PROFILING_COMMAND_END, sizeof(time_end), &time_end, NULL);

      double nanoSeconds = time_end-time_start;
      printf("OpenCl Execution time is: %0.3f milliseconds \n",nanoSeconds / 1000000.0);

      lgp+=deg;
      in_src = !in_src;
  }

    // Read back the results from the device to verify the output
    //
    start = high_resolution_clock::now();
    if(in_src) {
        err = clEnqueueReadBuffer(kern.commands, src, CL_TRUE, 0, sizeof(Fr<mnt4753_pp>) * count, results, 0, NULL, NULL );  
        if (err != CL_SUCCESS)
        {
            printf("Error: Failed to read output array! %d\n", err);
            exit(1);
        }
    } else {
        err = clEnqueueReadBuffer(kern.commands, dst, CL_TRUE, 0, sizeof(Fr<mnt4753_pp>) * count, results, 0, NULL, NULL );  
        if (err != CL_SUCCESS)
        {
            printf("Error: Failed to read output array! %d\n", err);
            exit(1);
        }
    }

    stop = high_resolution_clock::now();
    duration = duration_cast<microseconds>(stop - start); 
    cout << "Time taken by GPU read function: "
      << duration.count() << " microseconds" << endl;
    // Validate our results
    //
    printf("Kernel Result \n");
    results[0].print();

    printf("CPU Result\n");

    correct = 0;

    // if(results[0] == _h4_1) {
    //   correct++;
    // }

    
    // Print a brief summary detailing the results
    //
    //printf("Computed '%d/%d' correct fq3 values!\n", correct, count);
    // Shutdown and cleanup
    //

    clReleaseMemObject(src);
    clReleaseMemObject(dst);
    clReleaseMemObject(pq);
    clReleaseMemObject(om);
    //clReleaseProgram(kern.program);
    clReleaseKernel(kernel);
    //clReleaseCommandQueue(kern.commands);
    //clReleaseContext(kern.context);
    // OPENCL END
    //break;   
    //domain->data->iFFT_GPU(data, k);
}

void mnt4753_libsnark::domain_cosetFFT(
    mnt4753_libsnark::evaluation_domain *domain,
    mnt4753_libsnark::vector_Fr *a) {
  domain->data->cosetFFT(*a->data, Fr<mnt4753_pp>::multiplicative_generator);
}


void mnt4753_libsnark::domain_icosetFFT(
    mnt4753_libsnark::evaluation_domain *domain,
    mnt4753_libsnark::vector_Fr *a) {
  domain->data->icosetFFT(*a->data, Fr<mnt4753_pp>::multiplicative_generator);
}
void mnt4753_libsnark::domain_divide_by_Z_on_coset(
    mnt4753_libsnark::evaluation_domain *domain,
    mnt4753_libsnark::vector_Fr *a) {
  domain->data->divide_by_Z_on_coset(*a->data);
}
size_t
mnt4753_libsnark::domain_get_m(mnt4753_libsnark::evaluation_domain *domain) {
  return domain->data->m;
}

mnt4753_libsnark::G1 *
mnt4753_libsnark::multiexp_G1(mnt4753_libsnark::vector_Fr *scalar_start,
                              mnt4753_libsnark::vector_G1 *g_start,
                              size_t length) {

  mnt4753_libsnark::G1 *t = new mnt4753_libsnark::G1{
      multiexp<libff::G1<mnt4753_pp>, Fr<mnt4753_pp>>(
          scalar_start->data->begin() + scalar_start->offset,
          g_start->data->begin(), length)};
  //printf("CPU multiexp result\n");
  //t->data.print();
  return t;
}

mnt4753_libsnark::G1 *
mnt4753_libsnark::multiexp_G1_GPU(mnt4753_libsnark::vector_Fr *scalar_start,
                              mnt4753_libsnark::vector_G1 *g_start,
                              size_t length, Kernel kern) {
  printf("GPU MULTI_EXP STARTING\n");
  size_t NUM_GROUPS  = 224;
  size_t WINDOW_SIZE = 8;
  size_t NUM_WINDOWS = 96;
  size_t BUCKET_LEN = 1 << WINDOW_SIZE;

  std::vector<Fr<mnt4753_pp>> &scalar_data = *scalar_start->data;
  std::vector<libff::G1<mnt4753_pp>> &g_data = *g_start->data;
  
  // compute table with kernel
  cl_kernel kernel;                   // compute kernel
  cl_event event;                     // timing
  cl_ulong time_start;
  cl_ulong time_end;
  unsigned int n = length;
  
  bool skip = 0;
  bool *dm = new bool[length];
  libff::G1<mnt4753_pp> result;
  unsigned int correct;

  cl_mem g1_base_buffer;
  cl_mem g1_result_buffer;
  cl_mem g1_bucket_buffer;
  cl_mem exp_buffer;
  // not sure how to handle the skip buffer, i dont see it in libsnark yet
  // /depends/libff/algebra/scalar_mul/multiexp.tcc appears to have the impl
  //cl_mem dm_buffer;
  //cl_mem res;

  // Fill our data set with inputs from param gen
  //
  libff::G1<mnt4753_pp> *data_bases = new libff::G1<mnt4753_pp>[n];
  bigint<12> *data_scalars = new bigint<12>[n];

  for(int i = 1; i < n; i++) {
    memcpy(&data_bases[i-1], &g_data[i], sizeof(libff::G1<mnt4753_pp>));
  }

  // convert out of montgomery form!!!
  for(int i=1; i<n; i++) {
    data_scalars[i-1] = scalar_data[i].as_bigint();
  }
  
  //printf("copied base: \n");
  //data_scalars[0].print();

  //exit(1);
  unsigned int count = n - 1;


  // Create the compute kernel in the program we wish to run
  //
  kernel = clCreateKernel(kern.program, "G1_batched_lookup_multiexp", &kern.err);
  if (!kernel || kern.err != CL_SUCCESS)
  {
      printf("Error: Failed to create compute kernel!\n");
      exit(1);
  }

  
  // Create the input and output arrays in device memory for our calculation
  //
  printf("creating buffer\n");
  g1_base_buffer = clCreateBuffer(kern.context,  CL_MEM_READ_WRITE,  sizeof(libff::G1<mnt4753_pp>) * count, NULL, NULL);
  g1_result_buffer = clCreateBuffer(kern.context,  CL_MEM_READ_WRITE,  sizeof(libff::G1<mnt4753_pp>) * NUM_WINDOWS * NUM_GROUPS, NULL, NULL);
  g1_bucket_buffer = clCreateBuffer(kern.context,  CL_MEM_READ_WRITE,  sizeof(libff::G1<mnt4753_pp>) * BUCKET_LEN * NUM_WINDOWS * NUM_GROUPS, NULL, NULL);
  exp_buffer = clCreateBuffer(kern.context,  CL_MEM_READ_WRITE,  sizeof(bigint<12>) * count, NULL, NULL);
  //dm_buffer = clCreateBuffer(kern.context, CL_MEM_READ_ONLY, sizeof(bool) * count, NULL, NULL);
  //res = clCreateBuffer(kern.context,  CL_MEM_READ_ONLY,  sizeof(libff::G1<mnt4753_pp>), NULL, NULL);

  // if (!g1_base_buffer || !g1_result_buffer)
  // {
  //     printf("Error: Failed to allocate device memory!\n");
  //     exit(1);
  // }
  // Write our data set into the input array in device memory 
  //
  auto start = high_resolution_clock::now();
  for(int i=0; i<length; i++) {

  }
  kern.err = clEnqueueWriteBuffer(kern.commands, g1_base_buffer, CL_TRUE, 0, sizeof(libff::G1<mnt4753_pp>) * count, data_bases, 0, NULL, NULL);
  if (kern.err != CL_SUCCESS)
  {
      printf("Error: Failed to write to base source array!\n");
      exit(1);
  }

  kern.err = clEnqueueWriteBuffer(kern.commands, exp_buffer, CL_TRUE, 0, sizeof(bigint<12>) * count, data_scalars, 0, NULL, NULL);
  if (kern.err != CL_SUCCESS)
  {
      printf("Error: Failed to write to source array!\n");
      exit(1);
  }
  // kern.err = clEnqueueWriteBuffer(kern.commands, dm_buffer, CL_TRUE, 0, sizeof(bool) * count, dm, 0, NULL, NULL);
  // if (kern.err != CL_SUCCESS)
  // {
  //     printf("Error: Failed to write to omega source array!\n");
  //     exit(1);
  // }
  auto stop = high_resolution_clock::now();
  auto duration = duration_cast<microseconds>(stop - start); 
  cout << "Time taken by GPU write function: "
    << duration.count() << " microseconds" << endl;

  // Set the arguments to our compute kernel
  //
  size_t lnn = length - 1;
  kern.err = 0;
  kern.err |= clSetKernelArg(kernel, 0, sizeof(cl_mem), &g1_base_buffer);
  kern.err |= clSetKernelArg(kernel, 1, sizeof(cl_mem), &g1_bucket_buffer);
  kern.err |= clSetKernelArg(kernel, 2, sizeof(cl_mem), &g1_result_buffer);
  kern.err |= clSetKernelArg(kernel, 3, sizeof(cl_mem), &exp_buffer);
  //kern.err |= clSetKernelArg(kernel, 4, sizeof(cl_mem), &dm_buffer);
  //kern.err |= clSetKernelArg(kernel, 4, sizeof(unsigned int), &skip);
  kern.err |= clSetKernelArg(kernel, 4, sizeof(unsigned int), &lnn);
  if (kern.err != CL_SUCCESS)
  {
      printf("Error: Failed to set kernel arguments! %d\n", kern.err);
      exit(1);
  }

  // Get the maximum work group size for executing the kernel on the device
  //
  kern.err = clGetKernelWorkGroupInfo(kernel, kern.devices[0], CL_KERNEL_WORK_GROUP_SIZE, sizeof(kern.local), &kern.local, NULL);
  if (kern.err != CL_SUCCESS)
  {
      printf("Error: Failed to retrieve kernel work group info! %d\n", kern.err);
      exit(1);
  }

  printf("Max work size: %u\n", kern.local);

  printf("queueing multi exp kernel\n");
  kern.global = NUM_WINDOWS * NUM_GROUPS;
  kern.local = 256;
  kern.err = clEnqueueNDRangeKernel(kern.commands, kernel, 1, NULL, &kern.global, &kern.local, 0, NULL, &event);
  if (kern.err)
  {
      printf("Error: Failed to execute kernel!\n");
      exit(1);
  }


  clWaitForEvents(1, &event);
  clFinish(kern.commands);

  //libff::G1<mnt4753_pp> acc = g_data[0];
  libff::G1<mnt4753_pp> *res = new libff::G1<mnt4753_pp>[NUM_WINDOWS * NUM_GROUPS];

  // Time kernel execution time without read/write
  //
  clGetEventProfilingInfo(event, CL_PROFILING_COMMAND_START, sizeof(time_start), &time_start, NULL);
  clGetEventProfilingInfo(event, CL_PROFILING_COMMAND_END, sizeof(time_end), &time_end, NULL);

  double nanoSeconds = time_end-time_start;
  printf("OpenCl Execution time is: %0.3f milliseconds \n",nanoSeconds / 1000000.0);

  
  // Read back the results from the device to verify the output
  //
  start = high_resolution_clock::now();

  kern.err = clEnqueueReadBuffer(kern.commands, g1_result_buffer, CL_TRUE, 0, sizeof(libff::G1<mnt4753_pp>) * (NUM_WINDOWS * NUM_GROUPS), res, 0, NULL, NULL );  
  if (kern.err != CL_SUCCESS)
  {
      printf("Error: Failed to read output array! %d\n", kern.err);
      exit(1);
  }
  
  stop = high_resolution_clock::now();
  duration = duration_cast<microseconds>(stop - start); 
  cout << "Time taken by GPU read function: "
    << duration.count() << " microseconds" << endl;
  // Validate our results
  //
  //printf("Kernel Result \n");
  
  //res[0].print();
  
  libff::G1<mnt4753_pp> acc = libff::G1<mnt4753_pp>::zero();
  unsigned int bits = 0;
  for(int i=0; i<NUM_WINDOWS; i++) {
    unsigned int w = std::min(static_cast<unsigned int>(WINDOW_SIZE), (768 - bits));
    for(int j=0; j<w; j++) {
      acc = acc.dbl();
    }
    for(int g=0; g<NUM_GROUPS; g++) {
      acc = acc + res[g * NUM_WINDOWS + i];
    }
    bits += w;
  }
  acc = acc + g_data[0];
  //acc.print();

  // if(results[0] == _h4_1) {
  //   correct++;
  // }

  
  // Print a brief summary detailing the results
  //
  //printf("Computed '%d/%d' correct fq3 values!\n", correct, count);
  // Shutdown and cleanup
  //

  clReleaseMemObject(g1_base_buffer);
  clReleaseMemObject(g1_result_buffer);
  clReleaseMemObject(g1_bucket_buffer);
  //clReleaseMemObject(dm_buffer);
  clReleaseMemObject(exp_buffer);
  //clReleaseProgram(kern.program);
  clReleaseKernel(kernel);
  //clReleaseCommandQueue(kern.commands);
  //clReleaseContext(kern.context);
  // OPENCL END
  //break;   
  //domain->data->iFFT_GPU(data, k);


  //libff::G1<mnt4753_pp> h4_1 = libff::G1<mnt4753_pp>::zero();

  //mnt4753_libsnark::G1 *res = mnt4753_libsnark::G1::zero();
  //return res;
  // return new mnt4753_libsnark::G1{
  //     multiexp<libff::G1<mnt4753_pp>, Fr<mnt4753_pp>>(
  //         scalar_start->data->begin() + scalar_start->offset,
  //         g_start->data->begin(), length)};
}

mnt4753_libsnark::G2 *
mnt4753_libsnark::multiexp_G2(mnt4753_libsnark::vector_Fr *scalar_start,
                              mnt4753_libsnark::vector_G2 *g_start,
                              size_t length) {
  return new mnt4753_libsnark::G2{
      multiexp<libff::G2<mnt4753_pp>, Fr<mnt4753_pp>>(
          scalar_start->data->begin() + scalar_start->offset,
          g_start->data->begin(), length)};
}

mnt4753_libsnark::groth16_input *
mnt4753_libsnark::read_input(const char *path,
                             mnt4753_libsnark::groth16_params *params) {
  return new mnt4753_libsnark::groth16_input(path, params->d, params->m);
}

mnt4753_libsnark::vector_Fr *
mnt4753_libsnark::input_w(mnt4753_libsnark::groth16_input *input) {
  return new mnt4753_libsnark::vector_Fr{.data = input->w, .offset = 0};
}
mnt4753_libsnark::vector_Fr *
mnt4753_libsnark::input_ca(mnt4753_libsnark::groth16_input *input) {
  return new mnt4753_libsnark::vector_Fr{.data = input->ca, .offset = 0};
}
mnt4753_libsnark::vector_Fr *mnt4753_libsnark::input_cb(groth16_input *input) {
  return new mnt4753_libsnark::vector_Fr{.data = input->cb, .offset = 0};
}
mnt4753_libsnark::vector_Fr *mnt4753_libsnark::input_cc(groth16_input *input) {
  return new vector_Fr{.data = input->cc, .offset = 0};
}
mnt4753_libsnark::field *mnt4753_libsnark::input_r(groth16_input *input) {
  return new mnt4753_libsnark::field{.data = input->r};
}

mnt4753_libsnark::groth16_params *
mnt4753_libsnark::read_params(const char *path) {
  return new mnt4753_libsnark::groth16_params(path);
}

size_t mnt4753_libsnark::params_d(mnt4753_libsnark::groth16_params *params) {
  return params->d;
}
size_t mnt4753_libsnark::params_m(mnt4753_libsnark::groth16_params *params) {
  return params->m;
}
mnt4753_libsnark::vector_G1 *
mnt4753_libsnark::params_A(mnt4753_libsnark::groth16_params *params) {
  return new mnt4753_libsnark::vector_G1{.data = params->A};
}
mnt4753_libsnark::vector_G1 *
mnt4753_libsnark::params_B1(mnt4753_libsnark::groth16_params *params) {
  return new mnt4753_libsnark::vector_G1{.data = params->B1};
}
mnt4753_libsnark::vector_G1 *
mnt4753_libsnark::params_L(mnt4753_libsnark::groth16_params *params) {
  return new mnt4753_libsnark::vector_G1{.data = params->L};
}
mnt4753_libsnark::vector_G1 *
mnt4753_libsnark::params_H(mnt4753_libsnark::groth16_params *params) {
  return new mnt4753_libsnark::vector_G1{.data = params->H};
}
mnt4753_libsnark::vector_G2 *
mnt4753_libsnark::params_B2(mnt4753_libsnark::groth16_params *params) {
  return new mnt4753_libsnark::vector_G2{.data = params->B2};
}

void mnt4753_libsnark::delete_G1(mnt4753_libsnark::G1 *a) { delete a; }
void mnt4753_libsnark::delete_G2(mnt4753_libsnark::G1 *a) { delete a; }
void mnt4753_libsnark::delete_vector_Fr(mnt4753_libsnark::vector_Fr *a) {
  delete a;
}
void mnt4753_libsnark::delete_vector_G1(mnt4753_libsnark::vector_G1 *a) {
  delete a;
}
void mnt4753_libsnark::delete_vector_G2(mnt4753_libsnark::vector_G2 *a) {
  delete a;
}
void mnt4753_libsnark::delete_groth16_input(
    mnt4753_libsnark::groth16_input *a) {
  delete a;
}
void mnt4753_libsnark::delete_groth16_params(
    mnt4753_libsnark::groth16_params *a) {
  delete a;
}
void mnt4753_libsnark::delete_evaluation_domain(
    mnt4753_libsnark::evaluation_domain *a) {
  delete a;
}

void mnt4753_libsnark::groth16_output_write(mnt4753_libsnark::G1 *A,
                                            mnt4753_libsnark::G2 *B,
                                            mnt4753_libsnark::G1 *C,
                                            const char *output_path) {
  FILE *out = fopen(output_path, "w");
  write_g1<mnt4753_pp>(out, A->data);
  write_g2<mnt4753_pp>(out, B->data);
  write_g1<mnt4753_pp>(out, C->data);
  fclose(out);
}
class mnt6753_libsnark::groth16_input {
public:
  std::shared_ptr<std::vector<Fr<mnt6753_pp>>> w;
  std::shared_ptr<std::vector<Fr<mnt6753_pp>>> ca, cb, cc;
  Fr<mnt6753_pp> r;

  groth16_input(const char *path, size_t d, size_t m) {
    w = std::make_shared<std::vector<libff::Fr<mnt6753_pp>>>(
        std::vector<libff::Fr<mnt6753_pp>>());
    ca = std::make_shared<std::vector<libff::Fr<mnt6753_pp>>>(
        std::vector<libff::Fr<mnt6753_pp>>());
    cb = std::make_shared<std::vector<libff::Fr<mnt6753_pp>>>(
        std::vector<libff::Fr<mnt6753_pp>>());
    cc = std::make_shared<std::vector<libff::Fr<mnt6753_pp>>>(
        std::vector<libff::Fr<mnt6753_pp>>());
    FILE *inputs = fopen(path, "r");

    for (size_t i = 0; i < m + 1; ++i) {
      w->emplace_back(read_fr<mnt6753_pp>(inputs));
    }

    for (size_t i = 0; i < d + 1; ++i) {
      ca->emplace_back(read_fr<mnt6753_pp>(inputs));
    }
    for (size_t i = 0; i < d + 1; ++i) {
      cb->emplace_back(read_fr<mnt6753_pp>(inputs));
    }
    for (size_t i = 0; i < d + 1; ++i) {
      cc->emplace_back(read_fr<mnt6753_pp>(inputs));
    }

    r = read_fr<mnt6753_pp>(inputs);

    fclose(inputs);
  }
};

class mnt6753_libsnark::groth16_params {
public:
  size_t d;
  size_t m;
  std::shared_ptr<std::vector<libff::G1<mnt6753_pp>>> A, B1, L, H;
  std::shared_ptr<std::vector<libff::G2<mnt6753_pp>>> B2;

  groth16_params(const char *path) {
    FILE *params = fopen(path, "r");
    d = read_size_t(params);
    m = read_size_t(params);
    A = std::make_shared<std::vector<libff::G1<mnt6753_pp>>>(
        std::vector<libff::G1<mnt6753_pp>>());
    B1 = std::make_shared<std::vector<libff::G1<mnt6753_pp>>>(
        std::vector<libff::G1<mnt6753_pp>>());
    L = std::make_shared<std::vector<libff::G1<mnt6753_pp>>>(
        std::vector<libff::G1<mnt6753_pp>>());
    H = std::make_shared<std::vector<libff::G1<mnt6753_pp>>>(
        std::vector<libff::G1<mnt6753_pp>>());
    B2 = std::make_shared<std::vector<libff::G2<mnt6753_pp>>>(
        std::vector<libff::G2<mnt6753_pp>>());
    for (size_t i = 0; i <= m; ++i) {
      A->emplace_back(read_g1<mnt6753_pp>(params));
    }
    for (size_t i = 0; i <= m; ++i) {
      B1->emplace_back(read_g1<mnt6753_pp>(params));
    }
    for (size_t i = 0; i <= m; ++i) {
      B2->emplace_back(read_g2<mnt6753_pp>(params));
    }
    for (size_t i = 0; i < m - 1; ++i) {
      L->emplace_back(read_g1<mnt6753_pp>(params));
    }
    for (size_t i = 0; i < d; ++i) {
      H->emplace_back(read_g1<mnt6753_pp>(params));
    }
    fclose(params);
  }
};

struct mnt6753_libsnark::evaluation_domain {
  std::shared_ptr<libfqfft::evaluation_domain<Fr<mnt6753_pp>>> data;
};

struct mnt6753_libsnark::field {
  Fr<mnt6753_pp> data;
};

struct mnt6753_libsnark::G1 {
  libff::G1<mnt6753_pp> data;
};

struct mnt6753_libsnark::G2 {
  libff::G2<mnt6753_pp> data;
};

struct mnt6753_libsnark::vector_Fr {
  std::shared_ptr<std::vector<Fr<mnt6753_pp>>> data;
  size_t offset;
};

struct mnt6753_libsnark::vector_G1 {
  std::shared_ptr<std::vector<libff::G1<mnt6753_pp>>> data;
};
struct mnt6753_libsnark::vector_G2 {
  std::shared_ptr<std::vector<libff::G2<mnt6753_pp>>> data;
};

void mnt6753_libsnark::init_public_params() {
  mnt6753_pp::init_public_params();
}

void mnt6753_libsnark::print_G1(mnt6753_libsnark::G1 *a) { a->data.print(); }

void mnt6753_libsnark::print_G2(mnt6753_libsnark::G2 *a) { a->data.print(); }

mnt6753_libsnark::evaluation_domain *
mnt6753_libsnark::get_evaluation_domain(size_t d) {
  return new evaluation_domain{
      .data = libfqfft::get_evaluation_domain<Fr<mnt6753_pp>>(d)};
}

mnt6753_libsnark::G1 *mnt6753_libsnark::G1_add(mnt6753_libsnark::G1 *a,
                                               mnt6753_libsnark::G1 *b) {
  return new mnt6753_libsnark::G1{.data = a->data + b->data};
}

mnt6753_libsnark::G1 *mnt6753_libsnark::G1_scale(field *a, G1 *b) {
  return new G1{.data = a->data * b->data};
}

void mnt6753_libsnark::vector_Fr_muleq(mnt6753_libsnark::vector_Fr *a,
                                       mnt6753_libsnark::vector_Fr *b,
                                       size_t size) {
  size_t a_off = a->offset, b_off = b->offset;
#ifdef MULTICORE
#pragma omp parallel for
#endif
  for (size_t i = 0; i < size; i++) {
    a->data->at(i + a_off) = a->data->at(i + a_off) * b->data->at(i + b_off);
  }
}

void mnt6753_libsnark::vector_Fr_subeq(mnt6753_libsnark::vector_Fr *a,
                                       mnt6753_libsnark::vector_Fr *b,
                                       size_t size) {
  size_t a_off = a->offset, b_off = b->offset;
#ifdef MULTICORE
#pragma omp parallel for
#endif
  for (size_t i = 0; i < size; i++) {
    a->data->at(i + a_off) = a->data->at(i + a_off) - b->data->at(i + b_off);
  }
}

mnt6753_libsnark::vector_Fr *
mnt6753_libsnark::vector_Fr_offset(mnt6753_libsnark::vector_Fr *a,
                                   size_t offset) {
  return new vector_Fr{.data = a->data, .offset = offset};
}

void mnt6753_libsnark::vector_Fr_copy_into(mnt6753_libsnark::vector_Fr *src,
                                           mnt6753_libsnark::vector_Fr *dst,
                                           size_t length) {
  std::copy(src->data->begin() + src->offset,
            src->data->begin() + src->offset + length, dst->data->begin());
}

mnt6753_libsnark::vector_Fr *mnt6753_libsnark::vector_Fr_zeros(size_t length) {
  return new mnt6753_libsnark::vector_Fr{
      .data = std::make_shared<std::vector<Fr<mnt6753_pp>>>(
          length, Fr<mnt6753_pp>::zero())};
}

void mnt6753_libsnark::domain_iFFT(mnt6753_libsnark::evaluation_domain *domain,
                                   mnt6753_libsnark::vector_Fr *a) {
  std::vector<Fr<mnt6753_pp>> &data = *a->data;
  domain->data->iFFT(data);
}

void mnt6753_libsnark::domain_iFFT_GPU(mnt6753_libsnark::evaluation_domain *domain,
                                   mnt6753_libsnark::vector_Fr *a, Kernel kern) {
  std::vector<Fr<mnt6753_pp>> &data = *a->data;

  //domain->data->iFFT(data);
}
void mnt6753_libsnark::domain_cosetFFT(
    mnt6753_libsnark::evaluation_domain *domain,
    mnt6753_libsnark::vector_Fr *a) {
  domain->data->cosetFFT(*a->data, Fr<mnt6753_pp>::multiplicative_generator);
}
void mnt6753_libsnark::domain_icosetFFT(
    mnt6753_libsnark::evaluation_domain *domain,
    mnt6753_libsnark::vector_Fr *a) {
  domain->data->icosetFFT(*a->data, Fr<mnt6753_pp>::multiplicative_generator);
}
void mnt6753_libsnark::domain_divide_by_Z_on_coset(
    mnt6753_libsnark::evaluation_domain *domain,
    mnt6753_libsnark::vector_Fr *a) {
  domain->data->divide_by_Z_on_coset(*a->data);
}
size_t
mnt6753_libsnark::domain_get_m(mnt6753_libsnark::evaluation_domain *domain) {
  return domain->data->m;
}

mnt6753_libsnark::G1 *
mnt6753_libsnark::multiexp_G1(mnt6753_libsnark::vector_Fr *scalar_start,
                              mnt6753_libsnark::vector_G1 *g_start,
                              size_t length) {

  return new mnt6753_libsnark::G1{
      multiexp<libff::G1<mnt6753_pp>, Fr<mnt6753_pp>>(
          scalar_start->data->begin() + scalar_start->offset,
          g_start->data->begin(), length)};
}
mnt6753_libsnark::G1 *
mnt6753_libsnark::multiexp_G1_GPU(mnt6753_libsnark::vector_Fr *scalar_start,
                              mnt6753_libsnark::vector_G1 *g_start,
                              size_t length, Kernel k) {

  return new mnt6753_libsnark::G1{
      multiexp<libff::G1<mnt6753_pp>, Fr<mnt6753_pp>>(
          scalar_start->data->begin() + scalar_start->offset,
          g_start->data->begin(), length)};
}

mnt6753_libsnark::G2 *
mnt6753_libsnark::multiexp_G2(mnt6753_libsnark::vector_Fr *scalar_start,
                              mnt6753_libsnark::vector_G2 *g_start,
                              size_t length) {
  return new mnt6753_libsnark::G2{
      multiexp<libff::G2<mnt6753_pp>, Fr<mnt6753_pp>>(
          scalar_start->data->begin() + scalar_start->offset,
          g_start->data->begin(), length)};
}

mnt6753_libsnark::groth16_input *
mnt6753_libsnark::read_input(const char *path,
                             mnt6753_libsnark::groth16_params *params) {
  return new mnt6753_libsnark::groth16_input(path, params->d, params->m);
}

mnt6753_libsnark::vector_Fr *
mnt6753_libsnark::input_w(mnt6753_libsnark::groth16_input *input) {
  return new mnt6753_libsnark::vector_Fr{.data = input->w, .offset = 0};
}
mnt6753_libsnark::vector_Fr *
mnt6753_libsnark::input_ca(mnt6753_libsnark::groth16_input *input) {
  return new mnt6753_libsnark::vector_Fr{.data = input->ca, .offset = 0};
}
mnt6753_libsnark::vector_Fr *mnt6753_libsnark::input_cb(groth16_input *input) {
  return new mnt6753_libsnark::vector_Fr{.data = input->cb, .offset = 0};
}
mnt6753_libsnark::vector_Fr *mnt6753_libsnark::input_cc(groth16_input *input) {
  return new vector_Fr{.data = input->cc, .offset = 0};
}
mnt6753_libsnark::field *mnt6753_libsnark::input_r(groth16_input *input) {
  return new mnt6753_libsnark::field{.data = input->r};
}

mnt6753_libsnark::groth16_params *
mnt6753_libsnark::read_params(const char *path) {
  return new mnt6753_libsnark::groth16_params(path);
}

size_t mnt6753_libsnark::params_d(mnt6753_libsnark::groth16_params *params) {
  return params->d;
}
size_t mnt6753_libsnark::params_m(mnt6753_libsnark::groth16_params *params) {
  return params->m;
}
mnt6753_libsnark::vector_G1 *
mnt6753_libsnark::params_A(mnt6753_libsnark::groth16_params *params) {
  return new mnt6753_libsnark::vector_G1{.data = params->A};
}
mnt6753_libsnark::vector_G1 *
mnt6753_libsnark::params_B1(mnt6753_libsnark::groth16_params *params) {
  return new mnt6753_libsnark::vector_G1{.data = params->B1};
}
mnt6753_libsnark::vector_G1 *
mnt6753_libsnark::params_L(mnt6753_libsnark::groth16_params *params) {
  return new mnt6753_libsnark::vector_G1{.data = params->L};
}
mnt6753_libsnark::vector_G1 *
mnt6753_libsnark::params_H(mnt6753_libsnark::groth16_params *params) {
  return new mnt6753_libsnark::vector_G1{.data = params->H};
}
mnt6753_libsnark::vector_G2 *
mnt6753_libsnark::params_B2(mnt6753_libsnark::groth16_params *params) {
  return new mnt6753_libsnark::vector_G2{.data = params->B2};
}

void mnt6753_libsnark::delete_G1(mnt6753_libsnark::G1 *a) { delete a; }
void mnt6753_libsnark::delete_G2(mnt6753_libsnark::G1 *a) { delete a; }
void mnt6753_libsnark::delete_vector_Fr(mnt6753_libsnark::vector_Fr *a) {
  delete a;
}
void mnt6753_libsnark::delete_vector_G1(mnt6753_libsnark::vector_G1 *a) {
  delete a;
}
void mnt6753_libsnark::delete_vector_G2(mnt6753_libsnark::vector_G2 *a) {
  delete a;
}
void mnt6753_libsnark::delete_groth16_input(
    mnt6753_libsnark::groth16_input *a) {
  delete a;
}
void mnt6753_libsnark::delete_groth16_params(
    mnt6753_libsnark::groth16_params *a) {
  delete a;
}
void mnt6753_libsnark::delete_evaluation_domain(
    mnt6753_libsnark::evaluation_domain *a) {
  delete a;
}

void mnt6753_libsnark::groth16_output_write(mnt6753_libsnark::G1 *A,
                                            mnt6753_libsnark::G2 *B,
                                            mnt6753_libsnark::G1 *C,
                                            const char *output_path) {
  FILE *out = fopen(output_path, "w");
  write_g1<mnt6753_pp>(out, A->data);
  write_g2<mnt6753_pp>(out, B->data);
  write_g1<mnt6753_pp>(out, C->data);
  fclose(out);
}
