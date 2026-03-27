#ifndef RCPPEIGEN_QUIET_H
#define RCPPEIGEN_QUIET_H

#include <Rcpp.h>
#if defined(__GNUC__) || defined(__clang__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wignored-attributes"
#endif
#include <RcppEigen.h>
#if defined(__GNUC__) || defined(__clang__)
#pragma GCC diagnostic pop
#endif

#endif /* RCPPEIGEN_QUIET_H */
