# cuda-fixnum for snark challenge

Use this code to get started on the snark challenge. It implements the logic needed to do field arithmetic. In particular, for the fields used by mnt4-753 and mnt6-753, this takes the pairwise product of two arrays of field elements. That is, it maps over two arrays.

See `main.cu` for the implementation

## Suggested steps to make a submission for the snark challenge

Here are some suggested steps to solve the tutorial stage and work towards a faster, GPU powered snark prover:

1. [Add Quadratic extension arithmetic (Tutorial, $150)](https://coinlist.co/build/coda/pages/problem-02-quadratic-extension-arithmetic)
2. [Add Cubic extension arithmetic (Tutorial, $150)](https://coinlist.co/build/coda/pages/problem-03-cubic-extension-arithmetic)
3. [Add Curve operations (Tutorial, $200)](https://coinlist.co/build/coda/pages/problem-04-curve-operations)

For each of these, the cuda kernel code will need to be changed (see `main.cu:21-35`)

And here is our best guess on how to effectively make a submission for the [full prover](https://coinlist.co/build/coda/pages/problem-07-groth16-prover-challenges) (up to $70,000, and $7,000 immediately for the first submission to 2x the speed)

The SNARK prover is composed of several FFTs and multiexponentiations. In the C++ reference implementation, the [FFTs are here](https://github.com/CodaProtocol/snark-challenge-prover-reference/blob/master/libsnark/main.cpp#L90) and the [multiexponentiations are here](https://github.com/CodaProtocol/snark-challenge-prover-reference/blob/master/libsnark/main.cpp#L201).

1. Once you've finished the tutorial, try improving the multi-exponentiations with on-GPU versions using the curve operations from the tutorial. Each multi-exponentiation can be seen as a [map-reduce, as explained here](https://youtu.be/81uR9W5PZ5M?t=772). The reduce part may be complicated to implement for GPU, so it may be a good idea to start by implementing the "map" part on GPU and the "reduce" part on CPU.
2. Do the multi-exponentiations entirely on-GPU using an on-GPU [reduce](https://github.com/NVIDIA/cuda-samples/tree/master/Samples/reduction)
3. Use an on-GPU FFT (see for example [cuFFT](https://developer.nvidia.com/cufft)), adapted to finite fields.
   You can find a C++ implementation of a finite-field FFT [here](https://github.com/CodaProtocol/snark-challenge-prover-reference/blob/master/depends/libfqfft/libfqfft/evaluation_domain/domains/basic_radix2_domain_aux.tcc#L46).

## To build and run

To build and run:

1. `./build.sh`
2. `./main compute inputs outputs`
3. `shasum outputs` should be `b0f4a59a4be1c878dd9698fae7f1be86d8261025`

you will need to edit /Makefile:GENCODES to match your GPU [see here](https://arnon.dk/matching-sm-architectures-arch-and-gencode-for-various-nvidia-cards/)
