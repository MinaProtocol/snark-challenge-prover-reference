#pragma once

#include <vector>
#include <CL/cl.h>
// #define __CL_ENABLE_EXCEPTIONS
#define DATA_SIZE (131072)
#define limbs_per_elem (12)
#include <typeinfo>
#include <unistd.h>
#include <string.h>
#include <errno.h>
#include <exception>

char *getcwd(char *buf, size_t size);

class Kernel {
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


class mnt4753_libsnark {
public:
  class groth16_input;

  class groth16_params;

  struct evaluation_domain;

  struct field;

  struct G1;

  struct G2;

  struct vector_Fr;

  struct vector_G1;

  struct vector_G2;

  static void init_public_params();

  static void print_G1(G1 *a);
  static void print_G2(G2 *a);

  static evaluation_domain *get_evaluation_domain(size_t d);

  static G1 *G1_add(G1 *a, G1 *b);
  static G1 *G1_scale(field *a, G1 *b);

  static void vector_Fr_muleq(vector_Fr *a, vector_Fr *b, size_t size);
  static void vector_Fr_subeq(vector_Fr *a, vector_Fr *b, size_t size);
  static vector_Fr *vector_Fr_offset(vector_Fr *a, size_t offset);
  static void vector_Fr_copy_into(vector_Fr *src, vector_Fr *dst,
                                  size_t length);
  static vector_Fr *vector_Fr_zeros(size_t length);

  static void domain_iFFT(evaluation_domain *domain, vector_Fr *a);
  static void domain_iFFT_GPU(evaluation_domain *domain, vector_Fr *a, Kernel k);
  static void domain_cosetFFT(evaluation_domain *domain, vector_Fr *a);
  static void domain_icosetFFT(evaluation_domain *domain, vector_Fr *a);
  static void domain_divide_by_Z_on_coset(evaluation_domain *domain,
                                          vector_Fr *a);
  static size_t domain_get_m(evaluation_domain *domain);

  static G1 *multiexp_G1(vector_Fr *scalar_start, vector_G1 *g_start,
                         size_t length);
  static G2 *multiexp_G2(vector_Fr *scalar_start, vector_G2 *g_start,
                         size_t length);

  static groth16_input *read_input(const char *path, groth16_params *params);

  static vector_Fr *input_w(groth16_input *input);
  static vector_Fr *input_ca(groth16_input *input);
  static vector_Fr *input_cb(groth16_input *input);
  static vector_Fr *input_cc(groth16_input *input);
  static field *input_r(groth16_input *input);

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
  static void delete_evaluation_domain(evaluation_domain *a);

  static void groth16_output_write(G1 *A, G2 *B, G1 *C,
                                   const char *output_path);
};
class mnt6753_libsnark {
public:
  class groth16_input;

  class groth16_params;

  struct evaluation_domain;

  struct field;

  struct G1;

  struct G2;

  struct vector_Fr;

  struct vector_G1;

  struct vector_G2;

  static void init_public_params();

  static void print_G1(G1 *a);
  static void print_G2(G2 *a);

  static evaluation_domain *get_evaluation_domain(size_t d);

  static G1 *G1_add(G1 *a, G1 *b);
  static G1 *G1_scale(field *a, G1 *b);

  static void vector_Fr_muleq(vector_Fr *a, vector_Fr *b, size_t size);
  static void vector_Fr_subeq(vector_Fr *a, vector_Fr *b, size_t size);
  static vector_Fr *vector_Fr_offset(vector_Fr *a, size_t offset);
  static void vector_Fr_copy_into(vector_Fr *src, vector_Fr *dst,
                                  size_t length);
  static vector_Fr *vector_Fr_zeros(size_t length);

  static void domain_iFFT(evaluation_domain *domain, vector_Fr *a);
  static void domain_iFFT_GPU(evaluation_domain *domain, vector_Fr *a, Kernel k);
  static void domain_cosetFFT(evaluation_domain *domain, vector_Fr *a);
  static void domain_icosetFFT(evaluation_domain *domain, vector_Fr *a);
  static void domain_divide_by_Z_on_coset(evaluation_domain *domain,
                                          vector_Fr *a);
  static size_t domain_get_m(evaluation_domain *domain);

  static G1 *multiexp_G1(vector_Fr *scalar_start, vector_G1 *g_start,
                         size_t length);
  static G2 *multiexp_G2(vector_Fr *scalar_start, vector_G2 *g_start,
                         size_t length);

  static groth16_input *read_input(const char *path, groth16_params *params);

  static vector_Fr *input_w(groth16_input *input);
  static vector_Fr *input_ca(groth16_input *input);
  static vector_Fr *input_cb(groth16_input *input);
  static vector_Fr *input_cc(groth16_input *input);
  static field *input_r(groth16_input *input);

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
  static void delete_evaluation_domain(evaluation_domain *a);

  static void groth16_output_write(G1 *A, G2 *B, G1 *C,
                                   const char *output_path);
};
