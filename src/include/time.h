//
// Created by junior on 2019/8/2.
//

#ifndef LSH_CPP_TIME_H
#define LSH_CPP_TIME_H

#include "lsh_cpp.h"

namespace LSH_CPP {
    using second = std::chrono::duration<double>;
    using millisecond =  std::chrono::duration<double, std::milli>;

    using TimeVar = std::chrono::high_resolution_clock::time_point;

#define second_duration(a) second(a).count()
#define millisecond_duration(a) millisecond(a).count()
#define duration(a) second_duration(a)
#define timeNow() std::chrono::high_resolution_clock::now()

    template<typename ReturnType, typename F, typename... Args>
    std::pair<ReturnType, double> compute_function_time(F func, Args &&... args) {
        TimeVar t1 = timeNow();
        ReturnType ret = func(std::forward<Args>(args)...);
        return {ret, duration(timeNow() - t1)};
    }

    template<typename F, typename ... Args>
    double compute_function_time(F func, Args &&... args) {
        TimeVar t1 = timeNow();
        func(std::forward<Args>(args)...);
        return duration(timeNow() - t1);
    }
}
#endif //LSH_CPP_TIME_H
