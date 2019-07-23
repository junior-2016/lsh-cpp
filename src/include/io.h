//
// Created by junior on 19-7-22.
//

#ifndef LSH_CPP_IO_H
#define LSH_CPP_IO_H

#include "lsh_cpp.h"

namespace LSH_CPP {
    template<typename Element, typename Container>
    inline void print_container(const Container &container) {
        std::cout << "[ ";
        std::for_each(begin(container), end(container), [](const Element &element) { std::cout << element << " "; });
        std::cout << "]" << std::endl;
    }

}
#endif //LSH_CPP_IO_H
