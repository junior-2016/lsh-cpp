## Require
- G++ or Clang++ (support C++17)
- cmake 
- Boost
- anaconda3(python3.7) with numpy and matplotlib installed.
- GSL(GNU Scientific Library), you can use command 
'sudo apt-get install libgsl-dev' on Ubuntu to install gsl directly.

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