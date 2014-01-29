/*
*
* cpubench.c - A simple CPU Benchmarking tool
* Author: Suyash Srijan
* Email: suyashsrijan@outlook.com
*
* This program calculates how much time your CPU takes to compute n digits of PI using Chudnovsky Algorithm
* (http://en.wikipedia.org/wiki/Chudnovsky_algorithm) and n prime numbers (http://en.wikipedia.org/wiki/Prime_number)
* and uses the GNU Multiple Precision Arithmetic Library for most of the computations.
*
* Compile using gcc : gcc -O3 -Wall -o cpubench cpubench.c -lgmp -lssl -lcrypto -fopenmp
*
*/

#include <gmp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <sys/resource.h>
#include <sys/utsname.h>
#include <openssl/md5.h>
#include <omp.h>

/* You can't compile this on Windows */
#ifdef _WIN32
#error Sorry, you cannot compile this program on Windows because of some *nix-specific code
#endif

/* Build timestamp */
#define build_time __TIME__
#define build_date __DATE__

/* Define color codes that we will be using */
#define TXTNORMAL  "\x1B[0m"
#define TXTRED     "\x1B[31m"
#define TXTYELLOW  "\x1B[33m"
#define TXTGREEN   "\x1B[32m"

/* Variables we require */
struct timespec start, end;
struct timespec pstart, pend;
unsigned long int i, ti, constant1, constant2, constant3;
unsigned long precision;
unsigned long long x = 0;
unsigned long long y = 0;
unsigned char digest[16];
char *oput;
int pnum = 0;
int tpnums = 0;
int u;
double bits;
mpz_t v1, v2, v3, v4, v5;
mpf_t V1, V2, V3, total, tmp, res;
mp_exp_t exponent;
MD5_CTX context;

/* Calculate log to the base 2 using GCC's bit scan reverse intrinsic */
static __inline__ unsigned int clc_log2(const unsigned int num)
{
    return ((num <= 1) ? 0 : 32 - (__builtin_clz(num - 1)));
}

/* Calculate MD5 checksum for verification */
static __inline__ char *clc_md5(const char *string)
{

    /* Allocate memory to store checksum */
    char *checksum = (char*)malloc(33);
    /* Initialize MD5 context */
    MD5_Init(&context);
    /* Compute MD5 hash */
    MD5_Update(&context, string, strlen(string));
    /* Store MD5 hash */
    MD5_Final(digest, &context);

    /* Format digest */
    for (u = 0; u < 16; ++u)
    {
        snprintf(&(checksum[u*2]), 3, "%02x", (unsigned int)digest[u]);
    }

    /* Return checksum */
    return checksum;
}

/* Calculate prime numbers */
static __inline__ int clc_prime(unsigned long long max)
{
    /* Get high-res time */
    clock_gettime(CLOCK_MONOTONIC_RAW, &pstart);

    #pragma omp parallel for shared (max) private (x, y, pnum) reduction (+:tpnums)

    /* Start computing primes */
    for (x = 2; x <= max; x++)
    {
        pnum = 1;

        for (y = 2; y < x; y++)
        {
            if (x % y == 0)
            {
                pnum = 0;
                break;
            }
        }
        tpnums = tpnums + pnum;
    }

    /* Get high-res time */
    clock_gettime(CLOCK_MONOTONIC_RAW, &pend);

    /* Calculate and print time taken */
    double ptime_taken = (double)(pend.tv_sec - pstart.tv_sec) + (double)(pend.tv_nsec - pstart.tv_nsec) / 1E9;
    printf("Done!\n\nTime taken (seconds): %lf\n", ptime_taken);

    /* Return total primes */
    return tpnums;
}

/* Calculate pi digits main function */
static __inline__ char *clc_pi(unsigned long dgts)
{
    /* Compute required iterations */
    unsigned long iters = (dgts / 15) + 1;

    /* Initialize variables */
    constant1 = 545140134;
    constant2 = 13591409;
    constant3 = 640320;
    bits = clc_log2(10);
    precision = (dgts * bits) + 1;
    mpf_set_default_prec(precision);
    mpz_inits(v1, v2, v3, v4, v5, NULL);
    mpf_inits(res, tmp, V1, V2, V3, total, NULL);
    mpf_set_ui(total, 0);
    mpf_sqrt_ui(tmp, 10005);
    mpf_mul_ui(tmp, tmp, 426880);

    /* Get high-res time */
    clock_gettime(CLOCK_MONOTONIC_RAW, &start);

    /* Print total iterations and start computation of digits */
    printf("Total iterations: %lu\n\n", iters - 1);

    /* Iterate and compute value using Chudnovsky Algorithm */
    for (i = 0x0; i < iters; i++)
    {
        ti = i * 3;
        mpz_fac_ui(v1, 6 * i);
        mpz_set_ui(v2, constant1);
        mpz_mul_ui(v2, v2, i);
        mpz_add_ui(v2, v2, constant2);
        mpz_fac_ui(v3, ti);
        mpz_fac_ui(v4, i);
        mpz_pow_ui(v4, v4, 3);
        mpz_ui_pow_ui(v5, constant3, ti);
        if ((1 & ti) == 1)
        {
            mpz_neg(v5, v5);
        }
        mpz_mul(v1, v1, v2);
        mpf_set_z(V1, v1);
        mpz_mul(v3, v3, v4);
        mpz_mul(v3, v3, v5);
        mpf_set_z(V2, v3);
        mpf_div(V3, V1, V2);
        mpf_add(total, total, V3);

        /* Print interations executed if debugging (I don't like spamming stdout unnecesarily) */
#ifdef DEBUG
        printf("Iteration %lu of %lu successfully executed\n", i, iters - 1);
#endif
    }

    /* Some final computations */
    mpf_ui_div(total, 1, total);
    mpf_mul(total, total, tmp);

    /* Get high-res time */
    clock_gettime(CLOCK_MONOTONIC_RAW, &end);

    /* Calculate and print time taken */
    double time_taken = (double)(end.tv_sec - start.tv_sec) + (double)(end.tv_nsec - start.tv_nsec) / 1E9;
    printf("Done!\n\nTime taken (seconds): %lf\n", time_taken);

    /* Store output */
    oput = mpf_get_str(NULL, &exponent, 10, dgts, total);

    /* Free up space consumed by variables */
    mpz_clears(v1, v2, v3, v4, v5, NULL);
    mpf_clears(res, tmp, V1, V2, V3, total, NULL);

    /* Return value */
    return oput;
}

/* Entry point of program */
int main(int argc, char *argv[])
{

    /* Set number of threads to split the computation between when doing multithreaded bench */
    int numthreads = omp_get_max_threads();
    omp_set_num_threads(numthreads);

    /* Variable declaration and initialization */
    unsigned long cpvalue = 10000;
    unsigned int base = 10;
    char *tmp_ptr;
    int pd = 0;
    int dd = 0;
    int threading = 0;

    /* Try setting process priority to highest */
    int returnvalue = setpriority(PRIO_PROCESS, (id_t)0, -20);
    if (returnvalue == -1)
    {
        printf("%sWARN: Unable to max out priority. Did you not run this app as root?%s\n", TXTYELLOW, TXTNORMAL);
    }

    /* Parse command line */
    if (argc == 4 && ((strcmp(argv[3], "--printdigits") == 0) || (strcmp(argv[3], "--nodigits") == 0) || (strcmp(argv[3], "--dumpdigits") == 0)))
    {
        cpvalue = strtol(argv[1], &tmp_ptr, base);
        threading = (strcmp(argv[2], "--singlethreaded") == 0) ? 1 : 0;
        threading = (strcmp(argv[2], "--multithreaded") == 0) ? 0 : 1;
        pd = (strcmp(argv[3], "--printdigits") == 0) ? 1 : 0;
        dd = (strcmp(argv[3], "--dumpdigits") == 0) ? 1 : 0;
    }

    /* Invalid command line parameters */
    else
    {
        fprintf(stderr, "%sError: Invalid command-line arguments!%s\nUsage: cpubench [value] [threading] [parameter]\nValue: Any number from 1 to 2^32-1\n(in case of single threaded bench, it will be used to compute primes from 1 to n (where n = value between 1 and 2^32-1) or n digits of PI (where n = value between 1 and 2^32-1)\nParameter:\n--printdigits : Prints all digits of PI on console\n--nodigits : Suppresses printing of digits of PI on console (Use this when doing multithreaded bench)\n--dumpdigits : Saves all the digits of PI to a text file\nThreading:\n--singlethreaded : Stresses only one core (PI)\n--multithreaded : Stresses all the cores (PRIMES)\n\nUsage example: cpubench 50000 --singlethreaded --printdigits\n", TXTRED, TXTNORMAL);
        exit(1);
    }

    /* Print introductory text */
    struct utsname uname_ptr;
    uname(&uname_ptr);
    printf("%s\n---------------------------------------------------------------", TXTGREEN);
    printf("\nCPU Bench v1.0 beta (%s)\nBuild date: %s %s\n", uname_ptr.machine, build_date, build_time);
    printf("---------------------------------------------------------------%s\n\n", TXTNORMAL);

    /* Check if digits isnt zero or below */
    if (cpvalue < 1)
    {
        fprintf(stderr, "%sError: Digit cannot be lower than 1%s\n", TXTRED, TXTNORMAL);
        exit(1);
    }

    /* Perform single threaded benchmark */
    if (threading == 1)
    {

        /* Calculate digits of pi */
        printf("Performing single-threaded benchmarking [PI]\nComputing %lu digits of PI...\n", cpvalue);
        char *digits_of_pi = clc_pi(cpvalue);

        /* Print the digits if user specified the --printdigits flag */
        if (pd == 1)
        {
            printf("Here are the digits:\n\n%.1s.%s\n", digits_of_pi, digits_of_pi + 1);
        }

        /* Save digits to text file if user specified the --dumpdigits flag */
        if (dd == 1)
        {
            FILE *file;
            if ((file = fopen("pidigits.txt", "w")) == NULL)
            {
                fprintf(stderr, "%sError while opening file%s\n", TXTRED, TXTNORMAL);
                exit(-1);
            }
            else
            {
                fprintf(file, "%.1s.%s\n", digits_of_pi, digits_of_pi + 1);
                fclose(file);
            }
        }

        /* Print MD5 checksum */
        char *md5 = clc_md5(digits_of_pi);
        printf("MD5 checksum (for verification): %s\n", md5);

        /* Free the memory */
        free(digits_of_pi);
    }

    /* Perform multi-threaded benchmark */
    else if (threading == 0)
    {

        printf("Performing multi-threaded benchmarking [Primes]\nComputing primes under %lu...\n", cpvalue);
        long int tot = clc_prime(cpvalue);
        printf("Total primes found are %lu\n", tot);

        /* Print MD5 checksum */
        char *buffer = malloc(sizeof(tot) + 1);
        sprintf(buffer, "%lu", tot);
        char *md5 = clc_md5(buffer);
        printf("MD5 checksum (for verification): %s\n", md5);

    }

    /* Time to go! */
    printf("Goodbye!\n");
    return 0;
}
