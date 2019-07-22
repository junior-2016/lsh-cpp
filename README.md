## Require
- G++ or Clang++ (support C++17)
- cmake 
- Boost
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