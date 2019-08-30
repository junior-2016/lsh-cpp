//
// Created by junior on 2019/8/30.
//

#ifndef LSH_CPP_PORTABILITY_H
#define LSH_CPP_PORTABILITY_H

/**
 * process platform portability issue
 */

// Generalize warning push/pop to disable some compiler warning or error.
// reference: https://github.com/facebook/folly/blob/master/folly/Portability.h
#if defined(_MSC_VER)
#define LSH_CPP_PUSH_WARNING __pragma(warning(push))
#define LSH_CPP_POP_WARNING __pragma(warning(pop))
// Disable the GCC warnings.
#define LSH_CPP_GNU_DISABLE_WARNING(warningName)
#define LSH_CPP_GCC_DISABLE_WARNING(warningName)
#define LSH_CPP_CLANG_DISABLE_WARNING(warningName)
#define LSH_CPP_MSVC_DISABLE_WARNING(warningNumber) \
  __pragma(warning(disable : warningNumber))
#elif defined(__GNUC__)
// Clang & GCC
#define LSH_CPP_PUSH_WARNING _Pragma("GCC diagnostic push")
#define LSH_CPP_POP_WARNING _Pragma("GCC diagnostic pop")
#define LSH_CPP_GNU_DISABLE_WARNING_INTERNAL2(warningName) #warningName
#define LSH_CPP_GNU_DISABLE_WARNING(warningName) \
  _Pragma(                                     \
      LSH_CPP_GNU_DISABLE_WARNING_INTERNAL2(GCC diagnostic ignored warningName))
#ifdef __clang__
#define LSH_CPP_CLANG_DISABLE_WARNING(warningName) \
  LSH_CPP_GNU_DISABLE_WARNING(warningName)
#define LSH_CPP_GCC_DISABLE_WARNING(warningName)
#else
#define LSH_CPP_CLANG_DISABLE_WARNING(warningName)
#define LSH_CPP_GCC_DISABLE_WARNING(warningName) \
  LSH_CPP_GNU_DISABLE_WARNING(warningName)
#endif
#define LSH_CPP_MSVC_DISABLE_WARNING(warningNumber)
#else
#define LSH_CPP_PUSH_WARNING
#define LSH_CPP_POP_WARNING
#define LSH_CPP_GNU_DISABLE_WARNING(warningName)
#define LSH_CPP_GCC_DISABLE_WARNING(warningName)
#define LSH_CPP_CLANG_DISABLE_WARNING(warningName)
#define LSH_CPP_MSVC_DISABLE_WARNING(warningNumber)
#endif

#endif //LSH_CPP_PORTABILITY_H
