/*
 * util.h
 *
 * author: tto@andrew.cmu.edu and cdaruwal@andrew.cmu.edu
 *
 */
#pragma once

#include <chrono>

#ifdef DEBUG
#define D(x) x
#else
#define D(x)
#endif

__inline__ int get_chunksize(int num_tasks, int num_threads) {
    return (num_tasks+num_threads-1)/num_threads;
}

using namespace std::chrono;
typedef std::chrono::high_resolution_clock Clock;
typedef std::chrono::duration<double> dsec;
