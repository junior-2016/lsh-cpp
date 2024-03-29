## Require
- g++ 9.* (support C++17 and C++ parallelism Extension) or clang++ 8.0.1
- Eigen3 3.3.7 with intel mkl and intel tbb
<br> Eigen official download [http://eigen.tuxfamily.org/index.php?title=Main_Page#Download]
<br> To use intel mkl as eigen's backend, refer to [https://eigen.tuxfamily.org/dox/TopicUsingIntelMKL.html]
<br> To use intel tbb as libstdc++ parallelism backend, refer to [https://github.com/intel/tbb]
- cmake 3.12
- Boost
- python2.7 development headers and some packages 
<br> On Ubuntu, install python requirement using the following command.
```bash
$ sudo apt-get install python-matplotlib python-numpy python-tk python2.7-dev
```
- GSL(GNU Scientific Library) 
<br> On Ubuntu, install gsl using the following command.
```bash
$ sudo apt-get install libgsl-dev
```

## Build

```bash
$ git clone 
$ cd lsh-cpp
$ git submodule update --init --recursive
$ mkdir build
$ cd build
$ cmake -DCMAKE_BUILD_TYPE=Release ..
(use cmake -DCMAKE_BUILD_TYPE=Debug .. to build in debug mode)
$ make
```

## Run
```bash
$ cd build
$ ./lsh_cpp_test # run lsh test case
$ ./lsh_benchmark # run lsh benchmark
```

## TODO
- [ ] LSH_CPP implementation
    - [ ] LSH benchmark test
    - [ ] LSH Forest impl
    - [ ] LSH ensemble impl
    - [x] Weight MinHash impl
    - [ ] Lean MinHash impl (减少内存使用/压缩数据及参数/内存缓冲池)
    - [ ] HyperLog / HyperLog++ impl

- [ ] LSH_CPP Large-scale data support
    - [ ] Map reduce compute
    - [ ] data persistence (redis storage) 
    - [ ] data save/read/remove/split/merge/... feature
    
- [ ] Other
    - [ ] 使用C++17 parallelism TS 加速部分for-loop/sort/reduce/find算法.
    - [ ] 加入对 xxhash , parallel_hash_map/set 的可选控制,
    当用户编译不指定这些三方库时改用std:: hash和std:: unordered_map.
    - [ ] Python interface export using pybind11
    - [ ] WeightMinHash 实现中需要 concept 的部分待重构(比如WeightMinHash
    的update()参数需要SetElementType存在weight()和value()接口)