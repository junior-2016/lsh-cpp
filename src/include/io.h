//
// Created by junior on 19-7-22.
//

#ifndef LSH_CPP_IO_H
#define LSH_CPP_IO_H

#include "lsh_cpp.h"

namespace LSH_CPP {
    // 注意C++顺序容器的模板一般是 Sequence_Container<Element,Allocator>
    template<template<typename /*Element*/, typename /*Alloc*/> typename Container,
            typename Element, typename _Alloc = std::allocator<Element>>
    inline void print_sequence_container(const Container<Element, _Alloc> &container) {
        std::cout << "[ ";
        std::for_each(std::begin(container), std::end(container),
                      [](const Element &element) { std::cout << element << " "; });
        std::cout << "]" << std::endl;
    }

    std::vector<std::string> get_document_from_file(const char *path) {
        std::ifstream file(path);
        if (!file.is_open()) {
            fprintf(stderr, "Your file %s don't exist.\n", path);
            std::exit(-1);
        }
        std::string line;
        std::vector<std::string> data;
        while (std::getline(file, line)) {
            data.push_back(line);
        }
        file.close();
        return data;
    }

    inline std::vector<std::string> get_document_from_file(const std::string &path) {
        return get_document_from_file(path.c_str());
    }

    /**
     *  get dna document from fastq file.
     *  Format:
     *  First Line: Head.
     *  Second Line: dna data.
     *  Third Line: +Head.
     *  Forth Line: Quality Score.
     */
    std::vector<std::string> get_document_from_fastq_file(const char *path) {
        std::ifstream file(path);
        if (!file.is_open()) {
            fprintf(stderr, "Your file %s don't exist.\n", path);
            std::exit(-1);
        }
        std::string line;
        std::vector<std::string> data;
        int order = 0;
        while (std::getline(file, line)) {
            if ((order++) % 4 == 1) data.push_back(line);
        }
        file.close();
        return data;
    }
}
#endif //LSH_CPP_IO_H
