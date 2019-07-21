## Require
- G++ or Clang++ (support C++14)
- cmake 3.12 
- Boost 1.65.1
- anaconda3(python3.7) with numpy and matplotlib installed.
- GSL(GNU Scientific Library), you can use command 
'sudo apt-get install libgsl-dev' on Ubuntu to install gsl directly.

## Build
```bash
$ cd lsh-cpp
$ mkdir build
$ cd build
$ cmake ..
$ make
```

## Run
```bash
$ cd build
$ ./lsh_cpp
```