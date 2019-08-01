## Require
- g++ 8.3.0 / clang++ 8.0.1 (support C++17)
- xsimd library \[option\]
<br>Clone xsimd library [https://github.com/QuantStack/xsimd.git] and then
install.
<br>If you want to check whether your cpu support 
the sse/avx instruction,clone [https://github.com/Mysticial/FeatureDetector.git] 
and test on your machine.
- cmake 3.12
- Boost
- anaconda3(python3.7) with numpy and matplotlib installed.
- GSL(GNU Scientific Library) <br>
 You can use command 'sudo apt-get install libgsl-dev' on Ubuntu to install gsl directly.

## Build

```bash
$ git clone 
$ cd lsh-cpp
$ git submodule update --init --recursive
$ mkdir build
$ cd build
$ cmake -DCMAKE_BUILD_TYPE=Release ..
(use cmake -DUSE_SIMD=ON to enable simd-feature)
(use cmake -DCMAKE_BUILD_TYPE=Debug .. to build in debug mode)
$ make
```

## Run
```bash
$ cd build
$ ./lsh_cpp
```
If you occur error : ./lsh_cpp: ~/anaconda3/lib/libstdc++.so.6: 
version `GLIBCXX_3.4.26' not found (required by ./lsh_cpp), please replace
your libstdc++.so.6 in anaconda3 to latest libstdc++ dynamic library. 
```bash
$ find / -name "libstdc++.so*" 
after this command you can find the latest libstdc++ dynamic library such as libstdc++.so.6.0.26
$ cp path_of_libstdc++.so.6.0.26 ~/anaconda3/lib/
$ cd ~/anaconda3/lib/
$ rm -rf libstdc++.so.6
$ ln -s libstdc++.so.6.0.26 libstdc++.so.6
```

## TODO
- [ ] LSH_CPP implementation
    - [ ] LSH benchmark test
    - [ ] LSH Forest impl
    - [ ] LSH ensemble impl
    - [ ] Weight MinHash impl
    - [ ] Lean MinHash impl (减少内存使用/压缩数据及参数/内存缓冲池)
    - [ ] HyperLog / HyperLog++ impl

- [ ] LSH_CPP Large-scale data support
    - [ ] Map reduce compute
    - [ ] data persistence (redis storage) 
    - [ ] data save/read/remove/split/merge/... feature
    
- [ ] Other
    - [ ] Python interface export using pybind11