#include <cassert>
#include <chrono>
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

using namespace libff;
using namespace libsnark;

// const multi_exp_method method = multi_exp_method_BDLO12;
const multi_exp_method method = multi_exp_method_bos_coster;

// #define __CL_ENABLE_EXCEPTIONS
#include <CL/cl.h>
#define DATA_SIZE (131072)
#define limbs_per_elem (12)
#include <typeinfo>
#include <unistd.h>
#include <string.h>
#include <errno.h>
#include <exception>

using namespace std::chrono; 
using namespace std;

char *getcwd(char *buf, size_t size);

static inline auto now() -> decltype(std::chrono::high_resolution_clock::now()) {
    return std::chrono::high_resolution_clock::now();
}

class Kernel: public exception {
  public:
    char* program_source_code;
    size_t program_source_code_size;
    int err;                            // error code returned from api calls
    char name[128];
    size_t global;                      // global domain size for our calculation
    size_t local;                       // local domain size for our calculation
    cl_device_id device_id;             // compute device id 
    cl_context context;                 // compute context
    cl_command_queue commands;          // compute command queue
    cl_program program;                 // compute program
    cl_device_id *devices;

  void init(int n) throw() {
    // OPENCL START
    printf("initializing GPU prover...");

    char cwd[1024];
    if (getcwd(cwd, sizeof(cwd)) != NULL) {
       // printf("Current working dir: %s\n", cwd);
    } else {
       perror("getcwd() error");
       exit(1);
    }

    FILE *fp;
    char *source_str;
    size_t source_size, program_size;

    fp = fopen("ocl_kernels/main.cl", "r");
    if (!fp) {
        fprintf(stderr, "could not open program file\n");
        exit(1);
    }

    program_source_code = (char*)malloc(400000);
    program_source_code_size = fread(program_source_code, 1, 400000, fp);
    fclose(fp);

    // Connect to a compute device
    //

    /* get platform number of OpenCL */
    cl_uint  num_platforms = 0;
    clGetPlatformIDs (0, NULL, &num_platforms);
    printf("num_platforms: %d\n", (int)num_platforms);

    /* allocate a segment of mem space, so as to store supported platform info */
    cl_platform_id *platforms = (cl_platform_id *) malloc (num_platforms * sizeof (cl_platform_id));

    /* get platform info */
    clGetPlatformIDs (num_platforms, platforms, NULL);

    /* get device number on platform */
    cl_uint num_devices = 0;
    clGetDeviceIDs (platforms[0], CL_DEVICE_TYPE_GPU, 0, NULL, &num_devices);
    printf("num_devices: %d\n", (int)num_devices);

    /* allocate a segment of mem space, to store device info, supported by platform */
    devices = (cl_device_id *) malloc (num_devices * sizeof (cl_device_id));

    /* get device info */
    clGetDeviceIDs (platforms[0], CL_DEVICE_TYPE_GPU, num_devices, devices, NULL);

    // int gpu = 1;
    // err = clGetDeviceIDs(NULL, gpu ? CL_DEVICE_TYPE_GPU : CL_DEVICE_TYPE_CPU, 1, &device_id, NULL);
    // if (err != CL_SUCCESS)
    // {
    //     printf("Error: Failed to create a device group!\n");
    //     return EXIT_FAILURE;
    // }

    printf("Device id: %u\n", devices[0]);

    clGetDeviceInfo(devices[0], CL_DEVICE_NAME, 128, name, NULL);
    fprintf(stdout, "Created a dispatch queue using the %s\n", name);

    // Create a compute context 
    //
    printf("creating context\n");
    context = clCreateContext(0, num_devices, devices, NULL, NULL, &err);
    if (!context)
    {
        printf("Error: Failed to create a compute context!\n");
        //return EXIT_FAILURE;
        exit(1);
    }

    // Create a command commands
    //
    commands = clCreateCommandQueue(context, devices[0], CL_QUEUE_PROFILING_ENABLE, &err);
    if (!commands)
    {
        printf("Error: Failed to create a command commands!\n");
        //return EXIT_FAILURE;
        exit(1);
    }

    // Create the compute program from the source buffer
    //
    program = clCreateProgramWithSource(context, 1, (const char **) &program_source_code, &program_source_code_size, &err);
    if (!program)
    {
        printf("Error: Failed to create compute program!\n");
        //return EXIT_FAILURE;
        exit(1);
    }

    // Build the program executable
    //
    printf("building program\n");
    char options[] = "-cl-opt-disable";
    err = clBuildProgram(program, num_devices, devices, options, NULL, NULL);
    if (err != CL_SUCCESS)
    {
        size_t len;
        char buffer[2048];
        //std::cerr << getErrorString(err) << std::endl;
        printf("Error: Failed to build program executable!\n");
        printf ("Message: %s\n",strerror(err));
        clGetProgramBuildInfo(program, device_id, CL_PROGRAM_BUILD_LOG, sizeof(buffer), buffer, &len);
        exit(1);
    }
    // END OCL INIT
  }
};

template<typename T>
void
print_time(T &t1, const char *str) {
    auto t2 = std::chrono::high_resolution_clock::now();
    auto tim = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
    printf("%s: %ld ms\n", str, tim);
    t1 = t2;
}

template<typename ppT>
class groth16_parameters {
  public:
    size_t d;
    size_t m;
    std::vector<G1<ppT>> A, B1, L, H;
    std::vector<G2<ppT>> B2;

  groth16_parameters(const char* path) {
    FILE* params = fopen(path, "r");
    d = read_size_t(params);
    m = read_size_t(params);
    for (size_t i = 0; i <= m; ++i) { A.emplace_back(read_g1<ppT>(params)); }
    for (size_t i = 0; i <= m; ++i) { B1.emplace_back(read_g1<ppT>(params)); }
    for (size_t i = 0; i <= m; ++i) { B2.emplace_back(read_g2<ppT>(params)); }
    for (size_t i = 0; i < m-1; ++i) { L.emplace_back(read_g1<ppT>(params)); }
    for (size_t i = 0; i < d; ++i) { H.emplace_back(read_g1<ppT>(params)); }
    fclose(params);
  }
};

template<typename ppT>
class groth16_input {
  public:
    std::vector<Fr<ppT>> w;
    std::vector<Fr<ppT>> ca, cb, cc;
    Fr<ppT> r;

  groth16_input(const char* path, size_t d, size_t m) {
    FILE* inputs = fopen(path, "r");

    for (size_t i = 0; i < m + 1; ++i) { w.emplace_back(read_fr<ppT>(inputs)); }

    for (size_t i = 0; i < d + 1; ++i) { ca.emplace_back(read_fr<ppT>(inputs)); }
    for (size_t i = 0; i < d + 1; ++i) { cb.emplace_back(read_fr<ppT>(inputs)); }
    for (size_t i = 0; i < d + 1; ++i) { cc.emplace_back(read_fr<ppT>(inputs)); }

    r = read_fr<ppT>(inputs);

    fclose(inputs);
  }
};

template<typename ppT>
class groth16_output {
  public:
    G1<ppT> A, C;
    G2<ppT> B;

  groth16_output(G1<ppT> &&A, G2<ppT> &&B, G1<ppT> &&C) :
    A(std::move(A)), B(std::move(B)), C(std::move(C)) {}

  void write(const char* path) {
    FILE* out = fopen(path, "w");
    write_g1<ppT>(out, A);
    write_g2<ppT>(out, B);
    write_g1<ppT>(out, C);
    fclose(out);
  }
};

// Here is where all the FFTs happen.
template<typename ppT>
std::vector<Fr<ppT>> compute_H(
    size_t d,
    std::vector<Fr<ppT>> &ca,
    std::vector<Fr<ppT>> &cb,
    std::vector<Fr<ppT>> &cc)
{
    // Begin witness map
    libff::enter_block("Compute the polynomial H");

    const std::shared_ptr<libfqfft::evaluation_domain<Fr<ppT>> > domain = libfqfft::get_evaluation_domain<Fr<ppT>>(d + 1);

    domain->iFFT(ca);
    domain->iFFT(cb);

    domain->cosetFFT(ca, Fr<ppT>::multiplicative_generator);
    domain->cosetFFT(cb, Fr<ppT>::multiplicative_generator);

    libff::enter_block("Compute evaluation of polynomial H on set T");
    std::vector<Fr<ppT>> &H_tmp = ca; // can overwrite ca because it is not used later
#ifdef MULTICORE
#pragma omp parallel for
#endif
    for (size_t i = 0; i < domain->m; ++i)
    {
        H_tmp[i] = ca[i]*cb[i];
    }
    std::vector<Fr<ppT>>().swap(cb); // destroy cb

    domain->iFFT(cc);

    domain->cosetFFT(cc, Fr<ppT>::multiplicative_generator);

#ifdef MULTICORE
#pragma omp parallel for
#endif
    for (size_t i = 0; i < domain->m; ++i)
    {
        H_tmp[i] = (H_tmp[i]-cc[i]);
    }

    domain->divide_by_Z_on_coset(H_tmp);

    libff::leave_block("Compute evaluation of polynomial H on set T");

    domain->icosetFFT(H_tmp, Fr<ppT>::multiplicative_generator);

    std::vector<Fr<ppT>> coefficients_for_H(domain->m+1, Fr<ppT>::zero());
#ifdef MULTICORE
#pragma omp parallel for
#endif
    for (size_t i = 0; i < domain->m; ++i)
    {
        coefficients_for_H[i] = H_tmp[i];
    }

    libff::leave_block("Compute the polynomial H");

    return coefficients_for_H;
}

template<typename G, typename Fr>
G multiexp(typename std::vector<Fr>::const_iterator scalar_start,
           typename std::vector<G>::const_iterator g_start,
           size_t length)
{
#ifdef MULTICORE
    const size_t chunks = omp_get_max_threads(); // to override, set OMP_NUM_THREADS env var or call omp_set_num_threads()
#else
    const size_t chunks = 1;
#endif

    return libff::multi_exp_with_mixed_addition<G,
                                                Fr,
                                                method>(
        g_start,
        g_start + length,
        scalar_start,
        scalar_start + length,
        chunks);

}

template<typename ppT>
int run_prover(
    const char* params_path,
    const char* input_path,
    const char* output_path)
{
    ppT::init_public_params();

    auto beginning = now();

    auto t = beginning;

    const size_t primary_input_size = 1;

    const groth16_parameters<ppT> parameters(params_path);
    print_time(t, "load params");

    auto t_main = t;

    const groth16_input<ppT> input(input_path, parameters.d, parameters.m);

    print_time(t, "load inputs");

    std::vector<Fr<ppT>> w  = std::move(input.w);
    std::vector<Fr<ppT>> ca = std::move(input.ca);
    std::vector<Fr<ppT>> cb = std::move(input.cb);
    std::vector<Fr<ppT>> cc = std::move(input.cc);

    // End reading of parameters and input

    libff::enter_block("Call to r1cs_gg_ppzksnark_prover");

    std::vector<Fr<ppT>> coefficients_for_H = compute_H<ppT>(
        parameters.d,
        ca, cb, cc);

    libff::enter_block("Compute the proof");
    libff::enter_block("Multi-exponentiations");

    // Now the 5 multi-exponentiations
    libff::enter_block("A G1 multiexp");
    G1<ppT> evaluation_At = multiexp<G1<ppT>, Fr<ppT>>(
        w.begin(), parameters.A.begin(), parameters.m + 1);
    libff::leave_block("A G1 multiexp");

    libff::enter_block("B G1 multiexp");
    G1<ppT> evaluation_Bt1 = multiexp<G1<ppT>, Fr<ppT>>(
        w.begin(), parameters.B1.begin(), parameters.m + 1);
    libff::leave_block("B G1 multiexp");

    libff::enter_block("B G2 multiexp");
    G2<ppT> evaluation_Bt2 = multiexp<G2<ppT>, Fr<ppT>>(
        w.begin(), parameters.B2.begin(), parameters.m + 1);
    libff::leave_block("B G2 multiexp");

    libff::enter_block("H G1 multiexp");
    G1<ppT> evaluation_Ht = multiexp<G1<ppT>, Fr<ppT>>(
        coefficients_for_H.begin(), parameters.H.begin(), parameters.d);
    libff::leave_block("H G1 multiexp");

    libff::enter_block("L G1 multiexp");
    G1<ppT> evaluation_Lt = multiexp<G1<ppT>, Fr<ppT>>(
        w.begin() + primary_input_size + 1,
        parameters.L.begin(),
        parameters.m - 1);
    libff::leave_block("L G1 multiexp");

    libff::G1<ppT> C = evaluation_Ht + evaluation_Lt + input.r * evaluation_Bt1; /*+ s *  g1_A  - (r * s) * pk.delta_g1; */

    libff::leave_block("Multi-exponentiations");
    libff::leave_block("Compute the proof");
    libff::leave_block("Call to r1cs_gg_ppzksnark_prover");

    print_time(t, "cpu");

    groth16_output<ppT> output(
      std::move(evaluation_At),
      std::move(evaluation_Bt2),
      std::move(C));

    print_time(t, "store");

    output.write(output_path);

    print_time(t_main, "Total time from input to output: ");




    cl_kernel kernel;                   // compute kernel
    cl_event event;                     // timing
    cl_ulong time_start;
    cl_ulong time_end;
    int n = 512;

    G1<mnt4753_pp>* data_x = new G1<mnt4753_pp>[1];              // original data set given to device
    G1<mnt4753_pp>* data_y = new G1<mnt4753_pp>[n];              // original data set given to device
    G1<mnt4753_pp>* results = new G1<mnt4753_pp>[1];           // results returned from device
    
    unsigned int correct;               // number of correct results returned

    cl_mem input_x;                       // device memory used for the input array
    cl_mem input_y;                       // device memory used for the input array
    cl_mem ocl_output;                       // device memory used for the input array
    // Fill our data set with field inputs from param gen
    //
    unsigned int count = n;
    mp_size_t num = 1;
    Kernel kern;
    kern.init(n);

    //memcpy(&data_x[0], &h4_1, sizeof(G1<mnt4753_pp>));
    data_x[0] = G1<mnt4753_pp>::zero();
    printf("count %u\n", n);
    data_x[0].print_coordinates();

    for(int i = 0; i < count; i++) {
      memcpy(&data_y[i], &data_x[0], sizeof(G1<mnt4753_pp>));
    }
    
    
    // Create the compute kernel in the program we wish to run
    //
    kernel = clCreateKernel(kern.program, "multiexp_G1", &kern.err);
    if (!kernel || kern.err != CL_SUCCESS)
    {
        printf("Error: Failed to create compute kernel!\n");
        exit(1);
    }

    // Create the input and output arrays in device memory for our calculation
    //
    printf("creating buffer\n");
    input_x = clCreateBuffer(kern.context,  CL_MEM_READ_ONLY,  sizeof(G1<mnt4753_pp>), NULL, NULL);
    input_y = clCreateBuffer(kern.context,  CL_MEM_READ_ONLY,  sizeof(G1<mnt4753_pp>) * count, NULL, NULL);
    ocl_output = clCreateBuffer(kern.context, CL_MEM_WRITE_ONLY, sizeof(G1<mnt4753_pp>), NULL, NULL);

    if (!input_x || !ocl_output)
    {
        printf("Error: Failed to allocate device memory!\n");
        exit(1);
    }

    // Write our data set into the input array in device memory 
    //
    auto start = high_resolution_clock::now();
    kern.err = clEnqueueWriteBuffer(kern.commands, input_x, CL_TRUE, 0, sizeof(G1<mnt4753_pp>), data_x, 0, NULL, NULL);
    if (kern.err != CL_SUCCESS)
    {
        printf("Error: Failed to write to source array!\n");
        exit(1);
    }
    kern.err = clEnqueueWriteBuffer(kern.commands, input_y, CL_TRUE, 0, sizeof(G1<mnt4753_pp>) * count, data_y, 0, NULL, NULL);
    if (kern.err != CL_SUCCESS)
    {
        printf("Error: Failed to write to source array!\n");
        exit(1);
    }
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(stop - start); 
    cout << "Time taken by GPU write function: "
      << duration.count() << " microseconds" << endl;

    // Set the arguments to our compute kernel
    //
    kern.err = 0;
    kern.err  = clSetKernelArg(kernel, 0, sizeof(cl_mem), &input_x);
    kern.err  = clSetKernelArg(kernel, 1, sizeof(cl_mem), &input_y);
    kern.err |= clSetKernelArg(kernel, 2, sizeof(cl_mem), &ocl_output);
    kern.err |= clSetKernelArg(kernel, 3, sizeof(unsigned int), &count);
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
    kern.global = count;
    //global = 1;
    printf("queueing kernel\n");
    kern.err = clEnqueueNDRangeKernel(kern.commands, kernel, 1, NULL, &kern.global, &kern.local, 0, NULL, &event);
    if (kern.err)
    {
        printf("Error: Failed to execute kernel!\n");
        return EXIT_FAILURE;
    }

    clWaitForEvents(1, &event);
    clFinish(kern.commands);

    // Time kernel execution time without read/write
    //
    clGetEventProfilingInfo(event, CL_PROFILING_COMMAND_START, sizeof(time_start), &time_start, NULL);
    clGetEventProfilingInfo(event, CL_PROFILING_COMMAND_END, sizeof(time_end), &time_end, NULL);

    double nanoSeconds = time_end-time_start;
    printf("OpenCl Execution time is: %0.3f milliseconds \n",nanoSeconds / 1000000.0);

    // Read back the results from the device to verify the output
    //
    start = high_resolution_clock::now();
    kern.err = clEnqueueReadBuffer(kern.commands, ocl_output, CL_TRUE, 0, sizeof(G1<mnt4753_pp>), results, 0, NULL, NULL );  
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
    printf("Kernel Result \n");
    results[0].print();

    // results[0].coeff_a.mont_repr.print_hex();
    // for(int i=0; i<12; i++) {
    //   //std::cout << "Length of array = " << (sizeof(results[1013].non_residue.mont_repr.data)/sizeof(*results[1013].non_residue.mont_repr.data)) << std::endl;
    //   cl_uint x;
    //   cl_uint y;
    //   x = (cl_uint)((results[0].coeff_a.mont_repr.data[i] & 0xFFFFFFFF00000000LL) >> 32);
    //   y = (cl_uint)(results[0].coeff_a.mont_repr.data[i] & 0xFFFFFFFFLL);
    //   gmp_printf("%Mx\n", results[0].coeff_a.mont_repr.data[i]);
    //   printf("%x\n", x);
    //   printf("%x\n", y);
    // }

    results[0].zero().print_coordinates();

    // for(int i=0; i<12; i++) {
    //   //printf("%x\n", results[1013].c0.mod.data[i]);
    //   //std::cout << "Length of array = " << (sizeof(results[1013].non_residue.mont_repr.data)/sizeof(*results[1013].non_residue.mont_repr.data)) << std::endl;
    //   cl_uint x;
    //   cl_uint y;
    //   x = (cl_uint)((results[1013].c0.one().mont_repr.data[i] & 0xFFFFFFFF00000000LL) >> 32);
    //   y = (cl_uint)(results[1013].c0.one().mont_repr.data[i] & 0xFFFFFFFFLL);
    //   gmp_printf("%Mx\n", results[1013].c0.one().mont_repr.data[i]);
    //   printf("%x\n", x);
    //   printf("%x\n", y);
    // }

    printf("CPU Result\n");
    G1<mnt4753_pp> _h4_1 = G1<mnt4753_pp>::zero();

    //for (size_t i = 0; i < n; ++i) { _h4_1 = _h4_1 + g4_1[i]; }
    // _h4_1 = _h4_1 + g4_1[0];
    // _h4_1 = _h4_1 + g4_1[1];
    // _h4_1 = _h4_1 + g4_1[2];
    // _h4_1 = _h4_1 + g4_1[3];
    // _h4_1 = _h4_1 + g4_1[4];
    // _h4_1.print();
    //  g4_1[1].X().print();
    correct = 0;

    // there is some fuckery on the results fqe struct, cant equality check mont_repr
    // if(results[0] == _h4_1) {
    //   correct++;
    // }

    
    // Print a brief summary detailing the results
    //
    //printf("Computed '%d/%d' correct fq3 values!\n", correct, count);
    // Shutdown and cleanup
    //
    clReleaseMemObject(input_x);
    clReleaseMemObject(input_y);
    clReleaseMemObject(ocl_output);
    clReleaseProgram(kern.program);
    clReleaseKernel(kernel);
    clReleaseCommandQueue(kern.commands);
    clReleaseContext(kern.context);

    // OPENCL END
    //break;    
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
      return run_prover<mnt4753_pp>(params_path, input_path, output_path);
    }
  } else if (curve == "MNT6753") {
    if (mode == "compute") {
      return run_prover<mnt6753_pp>(params_path, input_path, output_path);
    }
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
