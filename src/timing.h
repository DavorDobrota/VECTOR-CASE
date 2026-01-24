#ifndef TIMING_H
#define TIMING_H

#ifdef _WIN32
    #include <windows.h>
#else
    #include <time.h>
#endif

typedef struct {
#ifdef _WIN32
    LARGE_INTEGER start;
    LARGE_INTEGER frequency;
#else
    struct timespec start;
#endif
} Timer;

/**
 * Initialize and start the timer
 */
static void timer_start(Timer* timer) {
#ifdef _WIN32
    QueryPerformanceFrequency(&timer->frequency);
    QueryPerformanceCounter(&timer->start);
#else
    clock_gettime(CLOCK_MONOTONIC, &timer->start);
#endif
}

/**
 * Get elapsed time in seconds since timer_start was called
 */
static double timer_elapsed(const Timer* timer) {
#ifdef _WIN32
    LARGE_INTEGER end;
    QueryPerformanceCounter(&end);
    return (double)(end.QuadPart - timer->start.QuadPart) / (double)timer->frequency.QuadPart;
#else
    struct timespec end;
    clock_gettime(CLOCK_MONOTONIC, &end);
    return (double)(end.tv_sec - timer->start.tv_sec) + (double)(end.tv_nsec - timer->start.tv_nsec) * 1e-9;
#endif
}

#endif // TIMING_H
