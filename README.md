cpubench.c - A simple CPU Benchmarking tool<br />
Author: Suyash Srijan<br />
Email: suyashsrijan@outlook.com<br />

This program calculates how much time your CPU takes to compute n digits of PI using Chudnovsky Algorithm
(http://en.wikipedia.org/wiki/Chudnovsky_algorithm) and n prime numbers (http://en.wikipedia.org/wiki/Prime_number)
and uses the GNU Multiple Precision Arithmetic Library for most of the computations.</br>

Compile using gcc : gcc -O3 -Wall -o cpubench cpubench.c -lgmp -lssl -lcrypto -fopenmp

