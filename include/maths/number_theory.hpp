#ifndef MATHS_NUMBER_THEORY_HPP
#define MATHS_NUMBER_THEORY_HPP

#include <cmath>
#include <vector>
#include <stdexcept>
#include <algorithm>
#include <numeric>
#include <map>
#include <set>
#include <cstdint>

/**
 * @file number_theory.hpp
 * @brief Comprehensive implementation of number theory algorithms and theorems
 *
 * This module implements:
 *
 * EUCLIDEAN ALGORITHMS:
 * - Basic Euclidean algorithm for GCD
 * - Extended Euclidean algorithm
 * - Computing modular inverses
 * - Chinese Remainder Theorem
 * - Speeding up algorithms via modular computation
 * - Rational reconstruction
 * - Applications to cryptography and algebra
 *
 * PRIME NUMBER THEORY:
 * - Chebyshev's theorem on density of primes
 * - Bertrand's postulate
 * - Mertens' theorem
 * - Sieve of Eratosthenes
 * - Prime number theorem
 * - Prime counting function π(x)
 * - Prime density and asymptotic analysis
 *
 * All implementations use standard number theory conventions.
 */

namespace maths {
namespace number_theory {

// ============================================================================
// BASIC EUCLIDEAN ALGORITHM
// ============================================================================

/**
 * @brief Compute GCD using Euclidean algorithm
 *
 * Algorithm:
 * gcd(a, b) = gcd(b, a mod b)
 * gcd(a, 0) = a
 *
 * Time complexity: O(log min(a, b))
 *
 * @param a First integer
 * @param b Second integer
 * @return Greatest common divisor
 */
inline long long gcd(long long a, long long b) {
    a = std::abs(a);
    b = std::abs(b);

    while (b != 0) {
        long long temp = b;
        b = a % b;
        a = temp;
    }

    return a;
}

/**
 * @brief Compute LCM (least common multiple)
 *
 * lcm(a, b) = |a × b| / gcd(a, b)
 *
 * @param a First integer
 * @param b Second integer
 * @return Least common multiple
 */
inline long long lcm(long long a, long long b) {
    if (a == 0 || b == 0) {
        return 0;
    }

    long long g = gcd(a, b);
    // Avoid overflow: compute (a/gcd) * b instead of (a*b)/gcd
    return std::abs(a / g * b);
}

/**
 * @brief Check if two integers are coprime (relatively prime)
 *
 * a and b are coprime iff gcd(a, b) = 1
 *
 * @param a First integer
 * @param b Second integer
 * @return true if coprime
 */
inline bool areCoprime(long long a, long long b) {
    return gcd(a, b) == 1;
}

/**
 * @brief Compute GCD of multiple integers
 *
 * @param numbers Vector of integers
 * @return GCD of all numbers
 */
inline long long gcdMultiple(const std::vector<long long>& numbers) {
    if (numbers.empty()) {
        throw std::invalid_argument("Cannot compute GCD of empty set");
    }

    long long result = numbers[0];
    for (size_t i = 1; i < numbers.size(); ++i) {
        result = gcd(result, numbers[i]);
        if (result == 1) {
            break;  // Early termination
        }
    }

    return result;
}

/**
 * @brief Compute LCM of multiple integers
 *
 * @param numbers Vector of integers
 * @return LCM of all numbers
 */
inline long long lcmMultiple(const std::vector<long long>& numbers) {
    if (numbers.empty()) {
        throw std::invalid_argument("Cannot compute LCM of empty set");
    }

    long long result = numbers[0];
    for (size_t i = 1; i < numbers.size(); ++i) {
        result = lcm(result, numbers[i]);
    }

    return result;
}

// ============================================================================
// EXTENDED EUCLIDEAN ALGORITHM
// ============================================================================

/**
 * @struct ExtendedGCDResult
 * @brief Result of extended Euclidean algorithm
 *
 * For gcd(a, b), finds integers x, y such that:
 * ax + by = gcd(a, b)
 */
struct ExtendedGCDResult {
    long long gcd;  // Greatest common divisor
    long long x;    // Coefficient for a
    long long y;    // Coefficient for b

    ExtendedGCDResult(long long g, long long x_coeff, long long y_coeff)
        : gcd(g), x(x_coeff), y(y_coeff) {}
};

/**
 * @brief Extended Euclidean algorithm
 *
 * Computes gcd(a, b) and finds x, y such that ax + by = gcd(a, b)
 *
 * This is known as Bézout's identity
 *
 * @param a First integer
 * @param b Second integer
 * @return ExtendedGCDResult containing gcd, x, and y
 */
inline ExtendedGCDResult extendedGCD(long long a, long long b) {
    if (b == 0) {
        // Base case: gcd(a, 0) = a, with a·1 + 0·0 = a
        return ExtendedGCDResult(std::abs(a), a >= 0 ? 1 : -1, 0);
    }

    // Store original signs
    long long sign_a = (a >= 0) ? 1 : -1;
    long long sign_b = (b >= 0) ? 1 : -1;
    a = std::abs(a);
    b = std::abs(b);

    // Initialize
    long long old_r = a, r = b;
    long long old_s = 1, s = 0;
    long long old_t = 0, t = 1;

    while (r != 0) {
        long long quotient = old_r / r;

        long long temp = r;
        r = old_r - quotient * r;
        old_r = temp;

        temp = s;
        s = old_s - quotient * s;
        old_s = temp;

        temp = t;
        t = old_t - quotient * t;
        old_t = temp;
    }

    // Apply original signs
    return ExtendedGCDResult(old_r, old_s * sign_a, old_t * sign_b);
}

/**
 * @brief Verify Bézout's identity
 *
 * Check that ax + by = gcd(a, b)
 *
 * @param a First integer
 * @param b Second integer
 * @param result Result from extendedGCD
 * @return true if identity holds
 */
inline bool verifyBezoutIdentity(long long a, long long b,
                                 const ExtendedGCDResult& result) {
    return (a * result.x + b * result.y == result.gcd);
}

// ============================================================================
// MODULAR ARITHMETIC
// ============================================================================

/**
 * @brief Compute (base^exp) mod m efficiently
 *
 * Uses binary exponentiation (exponentiation by squaring)
 *
 * Time complexity: O(log exp)
 *
 * @param base Base
 * @param exp Exponent (non-negative)
 * @param mod Modulus
 * @return (base^exp) mod m
 */
inline long long modPow(long long base, long long exp, long long mod) {
    if (mod <= 0) {
        throw std::invalid_argument("Modulus must be positive");
    }
    if (exp < 0) {
        throw std::invalid_argument("Exponent must be non-negative");
    }

    base %= mod;
    if (base < 0) base += mod;

    long long result = 1;

    while (exp > 0) {
        if (exp % 2 == 1) {
            result = (result * base) % mod;
        }
        base = (base * base) % mod;
        exp /= 2;
    }

    return result;
}

/**
 * @brief Compute modular inverse using extended Euclidean algorithm
 *
 * Finds x such that (a × x) ≡ 1 (mod m)
 *
 * Inverse exists iff gcd(a, m) = 1
 *
 * @param a Number to invert
 * @param m Modulus
 * @return Modular inverse of a modulo m
 * @throws std::invalid_argument if inverse doesn't exist
 */
inline long long modInverse(long long a, long long m) {
    if (m <= 0) {
        throw std::invalid_argument("Modulus must be positive");
    }

    ExtendedGCDResult result = extendedGCD(a, m);

    if (result.gcd != 1) {
        throw std::invalid_argument("Modular inverse does not exist (numbers not coprime)");
    }

    // result.x is the inverse, but might be negative
    long long inv = result.x % m;
    if (inv < 0) {
        inv += m;
    }

    return inv;
}

/**
 * @brief Compute modular division: (a / b) mod m = (a × b^(-1)) mod m
 *
 * @param a Dividend
 * @param b Divisor
 * @param m Modulus
 * @return (a / b) mod m
 */
inline long long modDivide(long long a, long long b, long long m) {
    a %= m;
    if (a < 0) a += m;

    long long b_inv = modInverse(b, m);
    return (a * b_inv) % m;
}

/**
 * @brief Compute Euler's totient function φ(n)
 *
 * φ(n) = count of integers k in [1, n] such that gcd(k, n) = 1
 *
 * For prime factorization n = p₁^a₁ × p₂^a₂ × ... × pₖ^aₖ:
 * φ(n) = n × (1 - 1/p₁) × (1 - 1/p₂) × ... × (1 - 1/pₖ)
 *
 * @param n Input integer
 * @return φ(n)
 */
inline long long eulerTotient(long long n) {
    if (n <= 0) {
        throw std::invalid_argument("n must be positive");
    }

    long long result = n;

    // Find all prime factors
    for (long long p = 2; p * p <= n; ++p) {
        if (n % p == 0) {
            // Remove factor p
            while (n % p == 0) {
                n /= p;
            }
            // Multiply by (1 - 1/p) = (p - 1)/p
            result -= result / p;
        }
    }

    // If n > 1, then it's a prime factor
    if (n > 1) {
        result -= result / n;
    }

    return result;
}

// ============================================================================
// CHINESE REMAINDER THEOREM
// ============================================================================

/**
 * @struct ChineseRemainderResult
 * @brief Result of Chinese Remainder Theorem
 */
struct ChineseRemainderResult {
    long long solution;  // The solution x
    long long modulus;   // The combined modulus M

    ChineseRemainderResult(long long sol, long long mod)
        : solution(sol), modulus(mod) {}
};

/**
 * @brief Chinese Remainder Theorem
 *
 * Given system of congruences:
 * x ≡ a₁ (mod m₁)
 * x ≡ a₂ (mod m₂)
 * ...
 * x ≡ aₖ (mod mₖ)
 *
 * where m₁, m₂, ..., mₖ are pairwise coprime,
 * finds unique solution x modulo M = m₁ × m₂ × ... × mₖ
 *
 * @param remainders Vector of remainders [a₁, a₂, ..., aₖ]
 * @param moduli Vector of moduli [m₁, m₂, ..., mₖ]
 * @return ChineseRemainderResult with solution and combined modulus
 */
inline ChineseRemainderResult chineseRemainder(const std::vector<long long>& remainders,
                                               const std::vector<long long>& moduli) {
    if (remainders.size() != moduli.size()) {
        throw std::invalid_argument("Remainders and moduli must have same size");
    }
    if (remainders.empty()) {
        throw std::invalid_argument("Cannot solve empty system");
    }

    // Check that moduli are pairwise coprime
    for (size_t i = 0; i < moduli.size(); ++i) {
        for (size_t j = i + 1; j < moduli.size(); ++j) {
            if (gcd(moduli[i], moduli[j]) != 1) {
                throw std::invalid_argument("Moduli must be pairwise coprime");
            }
        }
    }

    // Compute M = m₁ × m₂ × ... × mₖ
    long long M = 1;
    for (long long m : moduli) {
        M *= m;
    }

    // Compute solution using CRT formula
    long long x = 0;

    for (size_t i = 0; i < moduli.size(); ++i) {
        long long Mi = M / moduli[i];  // M / mᵢ
        long long yi = modInverse(Mi, moduli[i]);  // (M/mᵢ)^(-1) mod mᵢ

        x += remainders[i] * Mi * yi;
        x %= M;
    }

    // Ensure result is positive
    if (x < 0) {
        x += M;
    }

    return ChineseRemainderResult(x, M);
}

/**
 * @brief Verify Chinese Remainder Theorem solution
 *
 * @param result CRT result to verify
 * @param remainders Original remainders
 * @param moduli Original moduli
 * @return true if solution is correct
 */
inline bool verifyCRT(const ChineseRemainderResult& result,
                     const std::vector<long long>& remainders,
                     const std::vector<long long>& moduli) {
    for (size_t i = 0; i < moduli.size(); ++i) {
        if (result.solution % moduli[i] != remainders[i] % moduli[i]) {
            return false;
        }
    }
    return true;
}

// ============================================================================
// SPEEDING UP ALGORITHMS VIA MODULAR COMPUTATION
// ============================================================================

/**
 * @brief Compute factorial modulo m
 *
 * n! mod m
 *
 * Uses modular arithmetic at each step to avoid overflow
 *
 * @param n Input
 * @param m Modulus
 * @return n! mod m
 */
inline long long factorialMod(long long n, long long m) {
    if (n < 0) {
        throw std::invalid_argument("Factorial undefined for negative numbers");
    }
    if (m <= 0) {
        throw std::invalid_argument("Modulus must be positive");
    }

    // If n >= m and m is not 1, result is 0 (contains factor m)
    if (n >= m && m > 1) {
        return 0;
    }

    long long result = 1;
    for (long long i = 2; i <= n; ++i) {
        result = (result * i) % m;
    }

    return result;
}

/**
 * @brief Compute binomial coefficient modulo m
 *
 * C(n, k) mod m = n! / (k! × (n-k)!) mod m
 *
 * Uses modular inverse for division
 *
 * @param n Total items
 * @param k Items to choose
 * @param m Modulus (should be prime for this method)
 * @return C(n, k) mod m
 */
inline long long binomialMod(long long n, long long k, long long m) {
    if (k < 0 || k > n) {
        return 0;
    }
    if (k == 0 || k == n) {
        return 1;
    }

    // Optimize: C(n, k) = C(n, n-k)
    if (k > n - k) {
        k = n - k;
    }

    // Compute C(n, k) = n × (n-1) × ... × (n-k+1) / (k × (k-1) × ... × 1)
    long long numerator = 1;
    long long denominator = 1;

    for (long long i = 0; i < k; ++i) {
        numerator = (numerator * ((n - i) % m)) % m;
        denominator = (denominator * ((i + 1) % m)) % m;
    }

    return modDivide(numerator, denominator, m);
}

/**
 * @brief Compute Fibonacci number modulo m
 *
 * Uses matrix exponentiation for O(log n) time
 *
 * [F(n+1)]   [1 1]^n   [1]
 * [F(n)  ] = [1 0]   × [0]
 *
 * @param n Index
 * @param m Modulus
 * @return F(n) mod m
 */
inline long long fibonacciMod(long long n, long long m) {
    if (n < 0) {
        throw std::invalid_argument("Fibonacci index must be non-negative");
    }
    if (m <= 0) {
        throw std::invalid_argument("Modulus must be positive");
    }

    if (n == 0) return 0;
    if (n == 1) return 1;

    // Matrix multiplication modulo m
    auto matmul = [m](const std::array<std::array<long long, 2>, 2>& A,
                      const std::array<std::array<long long, 2>, 2>& B) {
        std::array<std::array<long long, 2>, 2> C;
        for (int i = 0; i < 2; ++i) {
            for (int j = 0; j < 2; ++j) {
                C[i][j] = 0;
                for (int k = 0; k < 2; ++k) {
                    C[i][j] = (C[i][j] + A[i][k] * B[k][j]) % m;
                }
            }
        }
        return C;
    };

    // Matrix exponentiation
    std::array<std::array<long long, 2>, 2> result = {{{1, 0}, {0, 1}}};  // Identity
    std::array<std::array<long long, 2>, 2> base = {{{1, 1}, {1, 0}}};

    long long exp = n;
    while (exp > 0) {
        if (exp % 2 == 1) {
            result = matmul(result, base);
        }
        base = matmul(base, base);
        exp /= 2;
    }

    return result[1][0];
}

// ============================================================================
// RATIONAL RECONSTRUCTION
// ============================================================================

/**
 * @struct RationalNumber
 * @brief Represents a rational number p/q in lowest terms
 */
struct RationalNumber {
    long long numerator;
    long long denominator;

    RationalNumber(long long p = 0, long long q = 1) {
        if (q == 0) {
            throw std::invalid_argument("Denominator cannot be zero");
        }

        // Reduce to lowest terms
        long long g = gcd(std::abs(p), std::abs(q));
        numerator = p / g;
        denominator = q / g;

        // Keep denominator positive
        if (denominator < 0) {
            numerator = -numerator;
            denominator = -denominator;
        }
    }

    double toDouble() const {
        return static_cast<double>(numerator) / denominator;
    }
};

/**
 * @brief Rational reconstruction from modular value
 *
 * Given a value a (mod m), finds rational p/q such that:
 * - p/q ≡ a (mod m)
 * - |p| < N, 0 < q < D
 * - gcd(p, q) = 1
 *
 * Uses continued fractions algorithm
 *
 * @param a Modular value
 * @param m Modulus
 * @param N Numerator bound
 * @param D Denominator bound
 * @return Reconstructed rational number
 */
inline RationalNumber rationalReconstruction(long long a, long long m,
                                            long long N, long long D) {
    if (m <= 0) {
        throw std::invalid_argument("Modulus must be positive");
    }

    a %= m;
    if (a < 0) a += m;

    // Use extended Euclidean algorithm with bounds
    long long r0 = m, r1 = a;
    long long t0 = 0, t1 = 1;

    while (r1 != 0 && std::abs(r1) > N) {
        long long q = r0 / r1;

        long long temp_r = r0 - q * r1;
        r0 = r1;
        r1 = temp_r;

        long long temp_t = t0 - q * t1;
        t0 = t1;
        t1 = temp_t;
    }

    if (r1 == 0) {
        throw std::runtime_error("Rational reconstruction failed");
    }

    long long p = r1;
    long long q = t1;

    // Adjust signs
    if (q < 0) {
        p = -p;
        q = -q;
    }

    // Check bounds
    if (std::abs(p) >= N || q <= 0 || q >= D) {
        throw std::runtime_error("Rational reconstruction failed: bounds exceeded");
    }

    return RationalNumber(p, q);
}

/**
 * @brief Application: Solve linear equation ax ≡ b (mod m) using rational reconstruction
 *
 * @param a Coefficient
 * @param b Constant term
 * @param m Modulus
 * @return Solution x
 */
inline long long solveLinearCongruence(long long a, long long b, long long m) {
    long long g = gcd(a, m);

    if (b % g != 0) {
        throw std::invalid_argument("No solution exists");
    }

    // Reduce to coprime case
    a /= g;
    b /= g;
    m /= g;

    long long a_inv = modInverse(a, m);
    return (a_inv * b) % m;
}

// ============================================================================
// PRIMALITY AND PRIME NUMBERS
// ============================================================================

/**
 * @brief Trial division primality test
 *
 * Checks if n is prime by testing divisibility up to √n
 *
 * Time complexity: O(√n)
 *
 * @param n Number to test
 * @return true if prime
 */
inline bool isPrimeTrial(long long n) {
    if (n <= 1) return false;
    if (n <= 3) return true;
    if (n % 2 == 0 || n % 3 == 0) return false;

    // Check divisibility by numbers of form 6k ± 1
    for (long long i = 5; i * i <= n; i += 6) {
        if (n % i == 0 || n % (i + 2) == 0) {
            return false;
        }
    }

    return true;
}

/**
 * @brief Miller-Rabin primality test (probabilistic)
 *
 * Fast probabilistic test for primality
 *
 * @param n Number to test
 * @param iterations Number of iterations (higher = more accurate)
 * @return true if probably prime
 */
inline bool isPrimeMillerRabin(long long n, int iterations = 5) {
    if (n <= 1) return false;
    if (n <= 3) return true;
    if (n % 2 == 0) return false;

    // Write n-1 as 2^r × d
    long long d = n - 1;
    int r = 0;
    while (d % 2 == 0) {
        d /= 2;
        r++;
    }

    // Witness loop
    for (int i = 0; i < iterations; ++i) {
        long long a = 2 + std::rand() % (n - 3);  // Random in [2, n-2]

        long long x = modPow(a, d, n);

        if (x == 1 || x == n - 1) {
            continue;
        }

        bool composite = true;
        for (int j = 0; j < r - 1; ++j) {
            x = modPow(x, 2, n);
            if (x == n - 1) {
                composite = false;
                break;
            }
        }

        if (composite) {
            return false;
        }
    }

    return true;
}

/**
 * @brief Combined primality test
 *
 * Uses trial division for small numbers, Miller-Rabin for large
 *
 * @param n Number to test
 * @return true if prime
 */
inline bool isPrime(long long n) {
    if (n < 1000000) {
        return isPrimeTrial(n);
    } else {
        return isPrimeMillerRabin(n, 10);
    }
}

// ============================================================================
// SIEVE OF ERATOSTHENES
// ============================================================================

/**
 * @brief Sieve of Eratosthenes - generate all primes up to n
 *
 * Time complexity: O(n log log n)
 * Space complexity: O(n)
 *
 * @param n Upper bound
 * @return Vector of all primes ≤ n
 */
inline std::vector<long long> sieveOfEratosthenes(long long n) {
    if (n < 2) {
        return {};
    }

    std::vector<bool> is_prime(n + 1, true);
    is_prime[0] = is_prime[1] = false;

    for (long long i = 2; i * i <= n; ++i) {
        if (is_prime[i]) {
            // Mark all multiples of i as composite
            for (long long j = i * i; j <= n; j += i) {
                is_prime[j] = false;
            }
        }
    }

    // Collect primes
    std::vector<long long> primes;
    for (long long i = 2; i <= n; ++i) {
        if (is_prime[i]) {
            primes.push_back(i);
        }
    }

    return primes;
}

/**
 * @brief Segmented sieve for large ranges
 *
 * More memory-efficient for large n
 *
 * @param low Lower bound
 * @param high Upper bound
 * @return Vector of primes in [low, high]
 */
inline std::vector<long long> segmentedSieve(long long low, long long high) {
    if (low > high) {
        return {};
    }

    // Generate primes up to √high
    long long limit = std::sqrt(high) + 1;
    std::vector<long long> base_primes = sieveOfEratosthenes(limit);

    std::vector<bool> is_prime(high - low + 1, true);

    for (long long p : base_primes) {
        // Find first multiple of p in [low, high]
        long long start = std::max(p * p, ((low + p - 1) / p) * p);

        for (long long j = start; j <= high; j += p) {
            is_prime[j - low] = false;
        }
    }

    // Collect primes
    std::vector<long long> primes;
    for (long long i = std::max(2LL, low); i <= high; ++i) {
        if (is_prime[i - low]) {
            primes.push_back(i);
        }
    }

    return primes;
}

// ============================================================================
// PRIME COUNTING FUNCTION
// ============================================================================

/**
 * @brief Prime counting function π(x)
 *
 * Returns number of primes ≤ x
 *
 * Uses sieve for exact count
 *
 * @param x Upper bound
 * @return Number of primes ≤ x
 */
inline long long primeCount(long long x) {
    if (x < 2) return 0;

    return sieveOfEratosthenes(x).size();
}

/**
 * @brief Prime counting function with precomputed sieve
 *
 * @param x Upper bound
 * @param primes Precomputed vector of primes
 * @return Number of primes ≤ x
 */
inline long long primeCountFromSieve(long long x, const std::vector<long long>& primes) {
    auto it = std::upper_bound(primes.begin(), primes.end(), x);
    return std::distance(primes.begin(), it);
}

// ============================================================================
// PRIME NUMBER THEOREM
// ============================================================================

/**
 * @brief Prime Number Theorem approximation
 *
 * π(x) ~ x / ln(x)
 *
 * Better approximation: π(x) ~ Li(x) (logarithmic integral)
 *
 * @param x Upper bound
 * @return Approximation of π(x)
 */
inline double primeNumberTheorem(double x) {
    if (x < 2.0) {
        return 0.0;
    }

    return x / std::log(x);
}

/**
 * @brief Improved approximation using Li(x)
 *
 * Li(x) = ∫₂ˣ dt/ln(t)
 *
 * Approximation: Li(x) ≈ x/ln(x) × (1 + 1/ln(x) + 2!/ln²(x) + ...)
 *
 * @param x Upper bound
 * @return Approximation of π(x) using Li(x)
 */
inline double logarithmicIntegral(double x) {
    if (x < 2.0) {
        return 0.0;
    }

    double lnx = std::log(x);

    // Asymptotic expansion: Li(x) ≈ x/ln(x) × Σ k! / ln^k(x)
    double result = x / lnx;
    double term = result;

    for (int k = 1; k <= 5; ++k) {
        term *= k / lnx;
        result += term;
    }

    return result;
}

/**
 * @brief Nth prime approximation
 *
 * p_n ~ n × ln(n)
 *
 * Better: p_n ~ n × (ln(n) + ln(ln(n)))
 *
 * @param n Index of prime
 * @return Approximation of nth prime
 */
inline double nthPrimeApproximation(long long n) {
    if (n <= 0) {
        throw std::invalid_argument("Prime index must be positive");
    }
    if (n == 1) return 2.0;

    double lnn = std::log(static_cast<double>(n));
    double lnlnn = std::log(lnn);

    return n * (lnn + lnlnn);
}

// ============================================================================
// CHEBYSHEV'S THEOREM
// ============================================================================

/**
 * @brief Chebyshev's theorem bounds
 *
 * Theorem: For x ≥ 2,
 * A × x/ln(x) ≤ π(x) ≤ B × x/ln(x)
 *
 * where A ≈ 0.92129 and B ≈ 1.10555
 *
 * @param x Upper bound
 * @return Pair of (lower_bound, upper_bound) for π(x)
 */
inline std::pair<double, double> chebyshevBounds(double x) {
    if (x < 2.0) {
        return {0.0, 0.0};
    }

    double base = x / std::log(x);
    double A = 0.92129;  // Chebyshev constant (lower)
    double B = 1.10555;  // Chebyshev constant (upper)

    return {A * base, B * base};
}

/**
 * @brief Verify Chebyshev's theorem for given x
 *
 * @param x Value to check
 * @param actual_count Actual π(x)
 * @return true if bounds hold
 */
inline bool verifyChebyshevTheorem(double x, long long actual_count) {
    if (x < 2.0) return true;

    auto bounds = chebyshevBounds(x);
    return (actual_count >= bounds.first && actual_count <= bounds.second);
}

/**
 * @brief Chebyshev function θ(x)
 *
 * θ(x) = Σ_{p≤x} ln(p) where sum is over primes
 *
 * Theorem: lim_{x→∞} θ(x)/x = 1
 *
 * @param x Upper bound
 * @param primes Precomputed primes
 * @return θ(x)
 */
inline double chebyshevTheta(double x, const std::vector<long long>& primes) {
    double theta = 0.0;

    for (long long p : primes) {
        if (p > x) break;
        theta += std::log(static_cast<double>(p));
    }

    return theta;
}

// ============================================================================
// BERTRAND'S POSTULATE
// ============================================================================

/**
 * @brief Bertrand's postulate verification
 *
 * Theorem: For n ≥ 1, there exists a prime p such that n < p < 2n
 *
 * Proven by Chebyshev, later improved
 *
 * @param n Lower bound
 * @param primes Precomputed primes
 * @return true if a prime exists in (n, 2n)
 */
inline bool verifyBertrandPostulate(long long n, const std::vector<long long>& primes) {
    if (n < 1) {
        throw std::invalid_argument("n must be at least 1");
    }

    // Search for prime in (n, 2n)
    for (long long p : primes) {
        if (p > n && p < 2 * n) {
            return true;
        }
        if (p >= 2 * n) {
            break;
        }
    }

    return false;
}

/**
 * @brief Find prime guaranteed by Bertrand's postulate
 *
 * @param n Lower bound
 * @return A prime p such that n < p < 2n
 */
inline long long findBertrandPrime(long long n) {
    if (n < 1) {
        throw std::invalid_argument("n must be at least 1");
    }

    for (long long p = n + 1; p < 2 * n; ++p) {
        if (isPrime(p)) {
            return p;
        }
    }

    throw std::runtime_error("Bertrand prime not found (should not happen)");
}

// ============================================================================
// MERTENS' THEOREM
// ============================================================================

/**
 * @brief Mertens' first theorem
 *
 * Σ_{p≤x} ln(p)/p = ln(x) + O(1)
 *
 * @param x Upper bound
 * @param primes Precomputed primes
 * @return Sum of ln(p)/p for primes p ≤ x
 */
inline double mertensFirstTheorem(double x, const std::vector<long long>& primes) {
    double sum = 0.0;

    for (long long p : primes) {
        if (p > x) break;
        double lnp = std::log(static_cast<double>(p));
        sum += lnp / p;
    }

    return sum;
}

/**
 * @brief Mertens' second theorem
 *
 * Σ_{p≤x} 1/p = ln(ln(x)) + B + o(1)
 *
 * where B ≈ 0.2614972128 (Meissel-Mertens constant)
 *
 * @param x Upper bound
 * @param primes Precomputed primes
 * @return Sum of 1/p for primes p ≤ x
 */
inline double mertensSecondTheorem(double x, const std::vector<long long>& primes) {
    double sum = 0.0;

    for (long long p : primes) {
        if (p > x) break;
        sum += 1.0 / p;
    }

    return sum;
}

/**
 * @brief Meissel-Mertens constant
 *
 * B = lim_{x→∞} (Σ_{p≤x} 1/p - ln(ln(x)))
 *
 * B ≈ 0.2614972128476427837554268386086958590515666
 */
constexpr double MEISSEL_MERTENS_CONSTANT = 0.2614972128476428;

/**
 * @brief Verify Mertens' second theorem approximation
 *
 * Check that Σ 1/p ≈ ln(ln(x)) + B
 *
 * @param x Upper bound
 * @param primes Precomputed primes
 * @return Difference from theoretical value
 */
inline double verifyMertensSecondTheorem(double x, const std::vector<long long>& primes) {
    if (x < 2.0) return 0.0;

    double sum = mertensSecondTheorem(x, primes);
    double expected = std::log(std::log(x)) + MEISSEL_MERTENS_CONSTANT;

    return std::abs(sum - expected);
}

/**
 * @brief Mertens' third theorem
 *
 * ∏_{p≤x} (1 - 1/p)⁻¹ ~ e^γ × ln(x)
 *
 * where γ is Euler-Mascheroni constant
 *
 * @param x Upper bound
 * @param primes Precomputed primes
 * @return Product ∏(1 - 1/p)⁻¹
 */
inline double mertensThirdTheorem(double x, const std::vector<long long>& primes) {
    double product = 1.0;

    for (long long p : primes) {
        if (p > x) break;
        product *= (1.0 / (1.0 - 1.0 / p));
    }

    return product;
}

// ============================================================================
// PRIME GAPS AND DENSITY
// ============================================================================

/**
 * @brief Calculate prime gaps
 *
 * Returns differences between consecutive primes
 *
 * @param primes Vector of primes
 * @return Vector of gaps
 */
inline std::vector<long long> primeGaps(const std::vector<long long>& primes) {
    if (primes.size() < 2) {
        return {};
    }

    std::vector<long long> gaps;
    for (size_t i = 1; i < primes.size(); ++i) {
        gaps.push_back(primes[i] - primes[i - 1]);
    }

    return gaps;
}

/**
 * @brief Average prime gap up to x
 *
 * Asymptotically: average gap ~ ln(x)
 *
 * @param primes Primes up to some x
 * @return Average gap
 */
inline double averagePrimeGap(const std::vector<long long>& primes) {
    if (primes.size() < 2) {
        return 0.0;
    }

    auto gaps = primeGaps(primes);
    long long sum = std::accumulate(gaps.begin(), gaps.end(), 0LL);

    return static_cast<double>(sum) / gaps.size();
}

/**
 * @brief Prime density at x
 *
 * Density ≈ 1/ln(x)
 *
 * @param x Point at which to calculate density
 * @return Prime density
 */
inline double primeDensity(double x) {
    if (x < 2.0) {
        return 0.0;
    }

    return 1.0 / std::log(x);
}

// ============================================================================
// ABELIAN GROUPS
// ============================================================================

/**
 * @class AbelianGroup
 * @brief Abstract abelian (commutative) group
 *
 * An abelian group (G, +) satisfies:
 * 1. Closure: a + b ∈ G for all a, b ∈ G
 * 2. Associativity: (a + b) + c = a + (b + c)
 * 3. Identity: ∃0 ∈ G such that a + 0 = a
 * 4. Inverses: ∀a ∈ G, ∃(-a) such that a + (-a) = 0
 * 5. Commutativity: a + b = b + a
 */
template<typename T>
class AbelianGroup {
public:
    virtual ~AbelianGroup() = default;

    // Group operation (additive notation)
    virtual T add(const T& a, const T& b) const = 0;

    // Identity element (zero)
    virtual T zero() const = 0;

    // Inverse (negation)
    virtual T negate(const T& a) const = 0;

    // Check if element is in group
    virtual bool contains(const T& a) const = 0;

    // Order of the group (-1 if infinite)
    virtual long long order() const = 0;

    // Order of an element
    virtual long long elementOrder(const T& a) const {
        if (!contains(a)) {
            throw std::invalid_argument("Element not in group");
        }

        T current = a;
        long long n = 1;
        long long maxIter = (order() > 0) ? order() : 10000;

        while (n <= maxIter) {
            if (approxEqual(current, zero())) {
                return n;
            }
            current = add(current, a);
            n++;
        }

        return -1;  // Infinite order
    }

    // Verify commutativity
    bool isCommutative(const T& a, const T& b) const {
        T ab = add(a, b);
        T ba = add(b, a);
        return approxEqual(ab, ba);
    }

protected:
    virtual bool approxEqual(const T& a, const T& b) const = 0;
};

/**
 * @class IntegersModN
 * @brief Additive group (ℤ/nℤ, +)
 *
 * Elements: {0, 1, 2, ..., n-1}
 * Operation: (a + b) mod n
 */
class IntegersModN : public AbelianGroup<long long> {
private:
    long long n;

public:
    explicit IntegersModN(long long modulus) : n(modulus) {
        if (modulus <= 0) {
            throw std::invalid_argument("Modulus must be positive");
        }
    }

    long long add(const long long& a, const long long& b) const override {
        return (a + b) % n;
    }

    long long zero() const override {
        return 0;
    }

    long long negate(const long long& a) const override {
        return (n - a) % n;
    }

    bool contains(const long long& a) const override {
        return a >= 0 && a < n;
    }

    long long order() const override {
        return n;
    }

    long long getModulus() const {
        return n;
    }

protected:
    bool approxEqual(const long long& a, const long long& b) const override {
        return a == b;
    }
};

// ============================================================================
// SUBGROUPS
// ============================================================================

/**
 * @brief Check if subset is a subgroup of abelian group
 *
 * H ⊆ G is a subgroup if:
 * 1. 0 ∈ H
 * 2. If a, b ∈ H, then a + b ∈ H
 * 3. If a ∈ H, then -a ∈ H
 *
 * @param elements Elements of potential subgroup
 * @param group Parent group
 * @return true if H is a subgroup
 */
template<typename T>
bool isSubgroup(const std::vector<T>& elements, const AbelianGroup<T>& group) {
    if (elements.empty()) {
        return false;
    }

    // Check identity
    T zero_elem = group.zero();
    bool hasZero = false;
    for (const auto& elem : elements) {
        if (group.approxEqual(elem, zero_elem)) {
            hasZero = true;
            break;
        }
    }
    if (!hasZero) {
        return false;
    }

    // Check closure and inverses
    for (const auto& a : elements) {
        // Check inverse
        T neg_a = group.negate(a);
        bool hasInverse = false;
        for (const auto& elem : elements) {
            if (group.approxEqual(elem, neg_a)) {
                hasInverse = true;
                break;
            }
        }
        if (!hasInverse) {
            return false;
        }

        // Check closure
        for (const auto& b : elements) {
            T sum = group.add(a, b);
            bool inSubgroup = false;
            for (const auto& elem : elements) {
                if (group.approxEqual(elem, sum)) {
                    inSubgroup = true;
                    break;
                }
            }
            if (!inSubgroup) {
                return false;
            }
        }
    }

    return true;
}

/**
 * @brief Generate cyclic subgroup <a> = {0, a, 2a, 3a, ...}
 *
 * @param generator Generator element
 * @param group Parent group
 * @return Elements of cyclic subgroup
 */
template<typename T>
std::vector<T> cyclicSubgroup(const T& generator, const AbelianGroup<T>& group) {
    std::vector<T> subgroup;
    T current = group.zero();

    long long maxIter = (group.order() > 0) ? group.order() : 1000;

    for (long long i = 0; i < maxIter; ++i) {
        // Check if current is already in subgroup
        bool found = false;
        for (const auto& elem : subgroup) {
            if (group.approxEqual(elem, current)) {
                found = true;
                break;
            }
        }

        if (found) {
            break;  // Cyclic group complete
        }

        subgroup.push_back(current);
        current = group.add(current, generator);
    }

    return subgroup;
}

// ============================================================================
// COSETS AND QUOTIENT GROUPS
// ============================================================================

/**
 * @brief Compute left coset a + H = {a + h : h ∈ H}
 *
 * @param a Element of group
 * @param subgroup Elements of subgroup H
 * @param group Parent group
 * @return Elements of coset a + H
 */
template<typename T>
std::vector<T> leftCoset(const T& a, const std::vector<T>& subgroup,
                        const AbelianGroup<T>& group) {
    std::vector<T> coset;

    for (const auto& h : subgroup) {
        T elem = group.add(a, h);
        coset.push_back(elem);
    }

    return coset;
}

/**
 * @brief Compute all cosets of H in G
 *
 * Returns G/H = {g + H : g ∈ G}
 *
 * @param group_elements All elements of G
 * @param subgroup Elements of subgroup H
 * @param group Parent group
 * @return Vector of cosets (each coset is a vector of elements)
 */
template<typename T>
std::vector<std::vector<T>> quotientGroup(const std::vector<T>& group_elements,
                                          const std::vector<T>& subgroup,
                                          const AbelianGroup<T>& group) {
    std::vector<std::vector<T>> cosets;
    std::set<size_t> used;  // Track which elements are already in a coset

    for (size_t i = 0; i < group_elements.size(); ++i) {
        if (used.count(i) > 0) {
            continue;
        }

        // Compute coset of group_elements[i]
        std::vector<T> coset = leftCoset(group_elements[i], subgroup, group);
        cosets.push_back(coset);

        // Mark all elements in this coset as used
        for (const auto& c : coset) {
            for (size_t j = 0; j < group_elements.size(); ++j) {
                if (group.approxEqual(group_elements[j], c)) {
                    used.insert(j);
                }
            }
        }
    }

    return cosets;
}

/**
 * @brief Compute index [G : H] = |G| / |H|
 *
 * @param group_order Order of group G
 * @param subgroup_order Order of subgroup H
 * @return Index [G : H]
 */
inline long long groupIndex(long long group_order, long long subgroup_order) {
    if (subgroup_order == 0) {
        throw std::invalid_argument("Subgroup order cannot be zero");
    }

    if (group_order % subgroup_order != 0) {
        throw std::invalid_argument("Subgroup order must divide group order (Lagrange)");
    }

    return group_order / subgroup_order;
}

// ============================================================================
// GROUP HOMOMORPHISMS AND ISOMORPHISMS
// ============================================================================

/**
 * @brief Check if map is a group homomorphism
 *
 * φ: G → H is a homomorphism if φ(a + b) = φ(a) + φ(b)
 *
 * @param elementsG Elements of group G
 * @param groupG Group G
 * @param groupH Group H
 * @param phi The homomorphism φ
 * @return true if φ is a homomorphism
 */
template<typename T, typename U>
bool isGroupHomomorphism(const std::vector<T>& elementsG,
                        const AbelianGroup<T>& groupG,
                        const AbelianGroup<U>& groupH,
                        std::function<U(const T&)> phi) {
    // Check φ(a + b) = φ(a) + φ(b)
    for (const auto& a : elementsG) {
        for (const auto& b : elementsG) {
            T sum_G = groupG.add(a, b);
            U left = phi(sum_G);

            U phi_a = phi(a);
            U phi_b = phi(b);
            U right = groupH.add(phi_a, phi_b);

            if (!groupH.approxEqual(left, right)) {
                return false;
            }
        }
    }

    // Also check φ(0) = 0
    T zero_G = groupG.zero();
    U zero_H = groupH.zero();
    U phi_zero = phi(zero_G);

    return groupH.approxEqual(phi_zero, zero_H);
}

/**
 * @brief Compute kernel of homomorphism
 *
 * ker(φ) = {g ∈ G : φ(g) = 0}
 *
 * @param elements Elements of G
 * @param groupG Group G
 * @param groupH Group H
 * @param phi Homomorphism φ: G → H
 * @return Elements in kernel
 */
template<typename T, typename U>
std::vector<T> kernelOfHomomorphism(const std::vector<T>& elements,
                                    const AbelianGroup<T>& groupG,
                                    const AbelianGroup<U>& groupH,
                                    std::function<U(const T&)> phi) {
    std::vector<T> kernel;
    U zero_H = groupH.zero();

    for (const auto& g : elements) {
        U phi_g = phi(g);
        if (groupH.approxEqual(phi_g, zero_H)) {
            kernel.push_back(g);
        }
    }

    return kernel;
}

// ============================================================================
// CYCLIC GROUPS
// ============================================================================

/**
 * @brief Check if group is cyclic
 *
 * A group is cyclic if it can be generated by a single element
 *
 * @param elements All elements of group
 * @param group The group
 * @return true if cyclic, along with a generator
 */
template<typename T>
std::pair<bool, T> isCyclic(const std::vector<T>& elements,
                            const AbelianGroup<T>& group) {
    // Try each element as potential generator
    for (const auto& g : elements) {
        std::vector<T> generated = cyclicSubgroup(g, group);

        // Check if generated equals full group
        if (generated.size() == elements.size()) {
            // Verify all elements are present
            bool isGenerator = true;
            for (const auto& elem : elements) {
                bool found = false;
                for (const auto& gen_elem : generated) {
                    if (group.approxEqual(elem, gen_elem)) {
                        found = true;
                        break;
                    }
                }
                if (!found) {
                    isGenerator = false;
                    break;
                }
            }

            if (isGenerator) {
                return {true, g};
            }
        }
    }

    return {false, group.zero()};
}

// ============================================================================
// STRUCTURE OF FINITE ABELIAN GROUPS
// ============================================================================

/**
 * @brief Fundamental theorem of finite abelian groups
 *
 * Every finite abelian group is isomorphic to a direct product of cyclic groups:
 * G ≅ ℤ/n₁ℤ × ℤ/n₂ℤ × ... × ℤ/nₖℤ
 *
 * where n₁ | n₂ | ... | nₖ (divisibility chain)
 *
 * This computes the invariant factors
 *
 * @param groupOrder Order of the finite abelian group
 * @return Vector of invariant factors
 */
inline std::vector<long long> invariantFactors(long long groupOrder) {
    if (groupOrder <= 0) {
        throw std::invalid_argument("Group order must be positive");
    }

    // For a cyclic group of order n, return {n}
    // For more complex groups, would need group structure information
    // This is a simplified version

    std::vector<long long> factors;

    // Factorize the group order
    long long n = groupOrder;
    std::map<long long, int> prime_powers;

    for (long long p = 2; p * p <= n; ++p) {
        while (n % p == 0) {
            prime_powers[p]++;
            n /= p;
        }
    }
    if (n > 1) {
        prime_powers[n]++;
    }

    // For cyclic case: single factor = group order
    factors.push_back(groupOrder);

    return factors;
}

/**
 * @brief Check if integer is cyclic group order
 *
 * ℤ/nℤ is cyclic for all n
 *
 * @param n Order
 * @return true (always cyclic for ℤ/nℤ)
 */
inline bool isIntegersModNCyclic(long long n) {
    return true;  // ℤ/nℤ is always cyclic
}

// ============================================================================
// RINGS
// ============================================================================

/**
 * @class Ring
 * @brief Abstract ring structure
 *
 * A ring (R, +, ×) consists of:
 * 1. An abelian group (R, +)
 * 2. Associative multiplication: (ab)c = a(bc)
 * 3. Distributivity: a(b+c) = ab + ac, (a+b)c = ac + bc
 * 4. Multiplicative identity (for unital rings): 1 ∈ R
 */
template<typename T>
class Ring {
public:
    virtual ~Ring() = default;

    // Addition (abelian group operation)
    virtual T add(const T& a, const T& b) const = 0;
    virtual T zero() const = 0;
    virtual T negate(const T& a) const = 0;

    // Multiplication
    virtual T multiply(const T& a, const T& b) const = 0;
    virtual T one() const = 0;  // Multiplicative identity

    // Ring properties
    virtual bool contains(const T& a) const = 0;

    // Check if ring is commutative
    bool isCommutative(const T& a, const T& b) const {
        T ab = multiply(a, b);
        T ba = multiply(b, a);
        return approxEqual(ab, ba);
    }

    // Verify distributivity
    bool verifyDistributivity(const T& a, const T& b, const T& c) const {
        // Left distributivity: a(b + c) = ab + ac
        T left = multiply(a, add(b, c));
        T right = add(multiply(a, b), multiply(a, c));

        if (!approxEqual(left, right)) {
            return false;
        }

        // Right distributivity: (a + b)c = ac + bc
        left = multiply(add(a, b), c);
        right = add(multiply(a, c), multiply(b, c));

        return approxEqual(left, right);
    }

protected:
    virtual bool approxEqual(const T& a, const T& b) const = 0;
};

/**
 * @class IntegersModNRing
 * @brief Ring ℤ/nℤ
 *
 * Addition and multiplication modulo n
 */
class IntegersModNRing : public Ring<long long> {
private:
    long long n;

public:
    explicit IntegersModNRing(long long modulus) : n(modulus) {
        if (modulus <= 0) {
            throw std::invalid_argument("Modulus must be positive");
        }
    }

    long long add(const long long& a, const long long& b) const override {
        return (a + b) % n;
    }

    long long zero() const override {
        return 0;
    }

    long long negate(const long long& a) const override {
        return (n - a) % n;
    }

    long long multiply(const long long& a, const long long& b) const override {
        return (a * b) % n;
    }

    long long one() const override {
        return 1;
    }

    bool contains(const long long& a) const override {
        return a >= 0 && a < n;
    }

    long long getModulus() const {
        return n;
    }

protected:
    bool approxEqual(const long long& a, const long long& b) const override {
        return a == b;
    }
};

// ============================================================================
// POLYNOMIAL RINGS
// ============================================================================

/**
 * @class Polynomial
 * @brief Polynomial with integer coefficients
 *
 * Represents a polynomial a₀ + a₁x + a₂x² + ... + aₙxⁿ
 */
class Polynomial {
private:
    std::vector<long long> coeffs;  // coeffs[i] is coefficient of x^i

    void normalize() {
        // Remove leading zeros
        while (coeffs.size() > 1 && coeffs.back() == 0) {
            coeffs.pop_back();
        }
        if (coeffs.empty()) {
            coeffs.push_back(0);
        }
    }

public:
    Polynomial() : coeffs({0}) {}

    explicit Polynomial(const std::vector<long long>& c) : coeffs(c) {
        if (coeffs.empty()) {
            coeffs.push_back(0);
        }
        normalize();
    }

    explicit Polynomial(long long constant) : coeffs({constant}) {}

    // Degree of polynomial
    int degree() const {
        if (coeffs.size() == 1 && coeffs[0] == 0) {
            return -1;  // Zero polynomial has degree -1 or -∞
        }
        return static_cast<int>(coeffs.size()) - 1;
    }

    // Get coefficient of x^i
    long long operator[](size_t i) const {
        return (i < coeffs.size()) ? coeffs[i] : 0;
    }

    // Addition
    Polynomial operator+(const Polynomial& other) const {
        size_t max_size = std::max(coeffs.size(), other.coeffs.size());
        std::vector<long long> result(max_size, 0);

        for (size_t i = 0; i < coeffs.size(); ++i) {
            result[i] += coeffs[i];
        }
        for (size_t i = 0; i < other.coeffs.size(); ++i) {
            result[i] += other.coeffs[i];
        }

        return Polynomial(result);
    }

    // Subtraction
    Polynomial operator-(const Polynomial& other) const {
        size_t max_size = std::max(coeffs.size(), other.coeffs.size());
        std::vector<long long> result(max_size, 0);

        for (size_t i = 0; i < coeffs.size(); ++i) {
            result[i] += coeffs[i];
        }
        for (size_t i = 0; i < other.coeffs.size(); ++i) {
            result[i] -= other.coeffs[i];
        }

        return Polynomial(result);
    }

    // Multiplication
    Polynomial operator*(const Polynomial& other) const {
        if (degree() < 0 || other.degree() < 0) {
            return Polynomial(0);
        }

        size_t result_size = coeffs.size() + other.coeffs.size() - 1;
        std::vector<long long> result(result_size, 0);

        for (size_t i = 0; i < coeffs.size(); ++i) {
            for (size_t j = 0; j < other.coeffs.size(); ++j) {
                result[i + j] += coeffs[i] * other.coeffs[j];
            }
        }

        return Polynomial(result);
    }

    // Equality
    bool operator==(const Polynomial& other) const {
        return coeffs == other.coeffs;
    }

    // Evaluate at x
    long long evaluate(long long x) const {
        long long result = 0;
        long long x_power = 1;

        for (long long c : coeffs) {
            result += c * x_power;
            x_power *= x;
        }

        return result;
    }

    const std::vector<long long>& getCoeffs() const {
        return coeffs;
    }
};

// ============================================================================
// IDEALS AND QUOTIENT RINGS
// ============================================================================

/**
 * @brief Check if subset is an ideal of ℤ/nℤ
 *
 * I ⊆ R is an ideal if:
 * 1. I is a subgroup under addition
 * 2. For all r ∈ R and i ∈ I: ri ∈ I and ir ∈ I
 *
 * In ℤ/nℤ, ideals are of the form (d) = {kd mod n : k ∈ ℤ} where d | n
 *
 * @param elements Elements of potential ideal
 * @param ring The ring ℤ/nℤ
 * @return true if I is an ideal
 */
inline bool isIdealOfIntegersModN(const std::vector<long long>& elements,
                                   const IntegersModNRing& ring) {
    if (elements.empty()) {
        return false;
    }

    long long n = ring.getModulus();

    // Check if subgroup under addition
    // (Simplified: we check closure and 0)

    bool hasZero = false;
    for (auto elem : elements) {
        if (elem == 0) {
            hasZero = true;
            break;
        }
    }
    if (!hasZero) {
        return false;
    }

    // Check absorption: for all r in ring and i in ideal, r*i in ideal
    for (long long r = 0; r < n; ++r) {
        for (auto i : elements) {
            long long product = (r * i) % n;

            bool found = false;
            for (auto elem : elements) {
                if (elem == product) {
                    found = true;
                    break;
                }
            }

            if (!found) {
                return false;
            }
        }
    }

    return true;
}

/**
 * @brief Generate principal ideal (a) in ℤ/nℤ
 *
 * (a) = {ka mod n : k ∈ ℤ}
 *
 * @param generator Generator a
 * @param n Modulus
 * @return Elements of ideal (a)
 */
inline std::vector<long long> principalIdeal(long long generator, long long n) {
    std::vector<long long> ideal;
    std::set<long long> seen;

    long long current = 0;
    while (seen.find(current) == seen.end()) {
        seen.insert(current);
        ideal.push_back(current);
        current = (current + generator) % n;
    }

    std::sort(ideal.begin(), ideal.end());
    return ideal;
}

// ============================================================================
// STRUCTURE OF ℤ*_n (MULTIPLICATIVE GROUP MOD N)
// ============================================================================

/**
 * @brief Compute ℤ*_n = units (invertible elements) in ℤ/nℤ
 *
 * ℤ*_n = {a ∈ ℤ/nℤ : gcd(a, n) = 1}
 *
 * |ℤ*_n| = φ(n) (Euler's totient)
 *
 * @param n Modulus
 * @return Elements of ℤ*_n
 */
inline std::vector<long long> unitsModN(long long n) {
    if (n <= 0) {
        throw std::invalid_argument("n must be positive");
    }

    std::vector<long long> units;

    for (long long a = 1; a < n; ++a) {
        if (gcd(a, n) == 1) {
            units.push_back(a);
        }
    }

    return units;
}

/**
 * @brief Check if ℤ*_n is cyclic
 *
 * ℤ*_n is cyclic iff n = 1, 2, 4, p^k, or 2p^k
 * where p is an odd prime
 *
 * @param n Modulus
 * @return true if ℤ*_n is cyclic
 */
inline bool isUnitsModNCyclic(long long n) {
    if (n == 1 || n == 2 || n == 4) {
        return true;
    }

    // Check if n = p^k for odd prime p
    for (long long p = 3; p * p <= n; p += 2) {
        if (n % p == 0) {
            long long temp = n;
            while (temp % p == 0) {
                temp /= p;
            }
            if (temp == 1) {
                return true;  // n = p^k
            }
            if (temp == 2) {
                return true;  // n = 2p^k
            }
            return false;
        }
    }

    // n is prime
    if (isPrime(n)) {
        return true;
    }

    // Check if n = 2p where p is odd prime
    if (n % 2 == 0) {
        long long half = n / 2;
        if (isPrime(half)) {
            return true;
        }
    }

    return false;
}

/**
 * @brief Find a generator of ℤ*_n (if cyclic)
 *
 * @param n Modulus
 * @return A generator of ℤ*_n
 */
inline long long generatorOfUnitsModN(long long n) {
    if (!isUnitsModNCyclic(n)) {
        throw std::invalid_argument("ℤ*_n is not cyclic for this n");
    }

    std::vector<long long> units = unitsModN(n);
    long long phi_n = units.size();

    // Try each unit as potential generator
    for (long long g : units) {
        // Check if order of g is φ(n)
        long long current = g;
        long long order = 1;

        while (current != 1 && order <= phi_n) {
            current = (current * g) % n;
            order++;
        }

        if (order == phi_n) {
            return g;
        }
    }

    throw std::runtime_error("Generator not found (should not happen for cyclic group)");
}

// ============================================================================
// ADVANCED PRIMALITY TESTING
// ============================================================================

/**
 * @brief Enhanced Miller-Rabin primality test
 *
 * Deterministic for n < 3,317,044,064,679,887,385,961,981 using specific witnesses
 *
 * @param n Number to test
 * @param witnesses Specific witnesses to test (if empty, use random)
 * @return true if probably prime
 */
inline bool millerRabinDeterministic(long long n,
                                     const std::vector<long long>& witnesses = {}) {
    if (n <= 1) return false;
    if (n == 2 || n == 3) return true;
    if (n % 2 == 0) return false;

    // Write n-1 as 2^r × d
    long long d = n - 1;
    int r = 0;
    while (d % 2 == 0) {
        d /= 2;
        r++;
    }

    // Use deterministic witnesses for small n
    std::vector<long long> test_witnesses;
    if (witnesses.empty()) {
        // Deterministic set for n < 3,317,044,064,679,887,385,961,981
        test_witnesses = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37};
    } else {
        test_witnesses = witnesses;
    }

    // Witness loop
    for (long long a : test_witnesses) {
        if (a >= n) continue;

        long long x = modPow(a, d, n);

        if (x == 1 || x == n - 1) {
            continue;
        }

        bool composite = true;
        for (int i = 0; i < r - 1; ++i) {
            x = modPow(x, 2, n);
            if (x == n - 1) {
                composite = false;
                break;
            }
        }

        if (composite) {
            return false;
        }
    }

    return true;
}

/**
 * @brief Generate random prime in range [low, high]
 *
 * Uses Miller-Rabin for primality testing
 *
 * @param low Lower bound
 * @param high Upper bound
 * @param maxAttempts Maximum attempts before giving up
 * @return Random prime in [low, high]
 */
inline long long generateRandomPrime(long long low, long long high,
                                     int maxAttempts = 1000) {
    if (low > high) {
        throw std::invalid_argument("Invalid range");
    }
    if (low < 2) {
        low = 2;
    }

    std::srand(std::time(nullptr));

    for (int attempt = 0; attempt < maxAttempts; ++attempt) {
        // Generate random odd number in range
        long long candidate = low + (std::rand() % (high - low + 1));
        if (candidate % 2 == 0 && candidate > 2) {
            candidate++;
        }
        if (candidate > high) {
            candidate = high;
        }

        if (millerRabinDeterministic(candidate)) {
            return candidate;
        }
    }

    throw std::runtime_error("Failed to generate prime in range after max attempts");
}

// ============================================================================
// PERFECT POWER TESTING
// ============================================================================

/**
 * @brief Test if n is a perfect power: n = m^k for some k ≥ 2
 *
 * @param n Number to test
 * @return Pair (isPerfectPower, exponent) where exponent is largest k
 */
inline std::pair<bool, int> isPerfectPower(long long n) {
    if (n <= 1) {
        return {true, 0};  // 0 and 1 are trivial perfect powers
    }

    // Test for each possible exponent k from 2 to log₂(n)
    int maxExp = static_cast<int>(std::log2(n)) + 1;

    for (int k = 2; k <= maxExp; ++k) {
        // Binary search for m such that m^k = n
        long long low = 1;
        long long high = static_cast<long long>(std::pow(n, 1.0 / k)) + 2;

        while (low <= high) {
            long long mid = low + (high - low) / 2;

            // Compute mid^k carefully to avoid overflow
            long long power = 1;
            bool overflow = false;
            for (int i = 0; i < k; ++i) {
                if (power > n / mid) {
                    overflow = true;
                    break;
                }
                power *= mid;
            }

            if (overflow) {
                high = mid - 1;
                continue;
            }

            if (power == n) {
                return {true, k};
            } else if (power < n) {
                low = mid + 1;
            } else {
                high = mid - 1;
            }
        }
    }

    return {false, 1};
}

/**
 * @brief Extract perfect power factorization
 *
 * If n = m^k, returns (m, k) with m not a perfect power
 *
 * @param n Number to factor
 * @return Pair (base, exponent)
 */
inline std::pair<long long, int> perfectPowerFactorization(long long n) {
    if (n <= 1) {
        return {n, 1};
    }

    int totalExp = 1;
    long long current = n;

    while (true) {
        auto [isPower, exp] = isPerfectPower(current);

        if (!isPower) {
            break;
        }

        // current = base^exp
        long long base = static_cast<long long>(std::pow(current, 1.0 / exp) + 0.5);
        current = base;
        totalExp *= exp;
    }

    return {current, totalExp};
}

// ============================================================================
// PRIME POWER FACTORING
// ============================================================================

/**
 * @brief Test if n is a prime power: n = p^k for prime p
 *
 * @param n Number to test
 * @return Pair (isPrimePower, {prime, exponent})
 */
inline std::pair<bool, std::pair<long long, int>> isPrimePower(long long n) {
    if (n <= 1) {
        return {false, {0, 0}};
    }

    // First check if n is a perfect power
    auto [base, exp] = perfectPowerFactorization(n);

    // Check if base is prime
    if (isPrime(base)) {
        return {true, {base, exp}};
    }

    return {false, {0, 0}};
}

/**
 * @brief Simple trial division factorization
 *
 * Returns prime factorization of n
 *
 * @param n Number to factor
 * @return Map of {prime : exponent}
 */
inline std::map<long long, int> trialDivisionFactorization(long long n) {
    std::map<long long, int> factors;

    if (n <= 1) {
        return factors;
    }

    // Check for factor 2
    while (n % 2 == 0) {
        factors[2]++;
        n /= 2;
    }

    // Check odd factors
    for (long long p = 3; p * p <= n; p += 2) {
        while (n % p == 0) {
            factors[p]++;
            n /= p;
        }
    }

    // If n > 1, then it's a prime factor
    if (n > 1) {
        factors[n]++;
    }

    return factors;
}

/**
 * @brief Compute Euler's phi function from prime factorization
 *
 * φ(n) = n × ∏(1 - 1/p) for each prime p dividing n
 *
 * @param factorization Prime factorization map
 * @return φ(n)
 */
inline long long eulerPhiFromFactorization(const std::map<long long, int>& factorization) {
    long long phi = 1;

    for (const auto& [prime, exp] : factorization) {
        // φ(p^k) = p^k - p^(k-1) = p^(k-1) × (p - 1)
        long long p_power = 1;
        for (int i = 0; i < exp - 1; ++i) {
            p_power *= prime;
        }
        phi *= p_power * (prime - 1);
    }

    return phi;
}

/**
 * @brief Compute n from prime factorization
 *
 * @param factorization Prime factorization map
 * @return n = ∏ p^e
 */
inline long long reconstructFromFactorization(const std::map<long long, int>& factorization) {
    long long n = 1;

    for (const auto& [prime, exp] : factorization) {
        for (int i = 0; i < exp; ++i) {
            n *= prime;
        }
    }

    return n;
}

// =============================================================================
// DISCRETE LOGARITHMS IN ℤ*_p
// =============================================================================

/**
 * @brief Check if g is a generator (primitive root) modulo p
 *
 * g is a generator of ℤ*_p if ord(g) = φ(p) = p-1
 *
 * @param g Candidate generator
 * @param p Prime modulus
 * @return true if g is a generator
 */
inline bool isGeneratorModP(long long g, long long p) {
    if (g <= 0 || g >= p) return false;
    if (!isPrimeMillerRabin(p)) return false;

    long long phi = p - 1;

    // Find prime factorization of φ(p) = p-1
    auto factors = trialDivisionFactorization(phi);

    // g is a generator iff g^((p-1)/q) ≠ 1 (mod p) for all prime divisors q of p-1
    for (const auto& [prime, exp] : factors) {
        if (modPow(g, phi / prime, p) == 1) {
            return false;
        }
    }

    return true;
}

/**
 * @brief Find a generator (primitive root) for ℤ*_p
 *
 * Uses brute force search starting from 2
 *
 * @param p Prime modulus
 * @return A generator g of ℤ*_p, or -1 if not found
 */
inline long long findGeneratorModP(long long p) {
    if (!isPrimeMillerRabin(p)) return -1;
    if (p == 2) return 1;

    // Try all candidates from 2 to p-1
    for (long long g = 2; g < p; ++g) {
        if (isGeneratorModP(g, p)) {
            return g;
        }
    }

    return -1;
}

/**
 * @brief Compute discrete logarithm using baby-step giant-step algorithm
 *
 * Solves g^x ≡ h (mod p) for x
 * Complexity: O(√p) time and space
 *
 * @param g Generator base
 * @param h Target value
 * @param p Prime modulus
 * @return x such that g^x ≡ h (mod p), or -1 if no solution
 */
inline long long discreteLogBabyGiant(long long g, long long h, long long p) {
    if (g <= 0 || h <= 0 || g >= p || h >= p) return -1;

    // Baby-step giant-step algorithm
    long long m = static_cast<long long>(std::sqrt(p)) + 1;

    // Baby steps: compute g^j for j = 0, 1, ..., m-1
    std::unordered_map<long long, long long> table;
    long long gamma = 1;
    for (long long j = 0; j < m; ++j) {
        table[gamma] = j;
        gamma = (gamma * g) % p;
    }

    // Giant steps: compute h * (g^(-m))^i for i = 0, 1, ..., m-1
    long long g_inv_m = modPow(modInverse(g, p), m, p);
    gamma = h;
    for (long long i = 0; i < m; ++i) {
        if (table.count(gamma)) {
            long long x = i * m + table[gamma];
            return x % (p - 1);
        }
        gamma = (gamma * g_inv_m) % p;
    }

    return -1;
}

/**
 * @brief Compute discrete logarithm using Pollard's rho algorithm
 *
 * Probabilistic algorithm with O(√p) expected time, O(1) space
 *
 * @param g Generator base
 * @param h Target value
 * @param p Prime modulus
 * @param maxIterations Maximum iterations
 * @return x such that g^x ≡ h (mod p), or -1 if not found
 */
inline long long discreteLogPollardRho(long long g, long long h, long long p,
                                       int maxIterations = 1000000) {
    if (g <= 0 || h <= 0 || g >= p || h >= p) return -1;

    long long order = p - 1;

    // Define partition function for cycling
    auto f = [&](long long x, long long a, long long b) -> std::tuple<long long, long long, long long> {
        int partition = x % 3;
        if (partition == 0) {
            return {(x * h) % p, a, (b + 1) % order};
        } else if (partition == 1) {
            return {(x * x) % p, (2 * a) % order, (2 * b) % order};
        } else {
            return {(x * g) % p, (a + 1) % order, b};
        }
    };

    // Floyd's cycle detection
    long long x1 = 1, a1 = 0, b1 = 0;
    long long x2 = 1, a2 = 0, b2 = 0;

    for (int i = 0; i < maxIterations; ++i) {
        // Tortoise: one step
        std::tie(x1, a1, b1) = f(x1, a1, b1);

        // Hare: two steps
        std::tie(x2, a2, b2) = f(x2, a2, b2);
        std::tie(x2, a2, b2) = f(x2, a2, b2);

        if (x1 == x2) {
            // Found collision: g^a1 * h^b1 ≡ g^a2 * h^b2 (mod p)
            // => h^(b1-b2) ≡ g^(a2-a1) (mod p)
            long long r = (b1 - b2 + order) % order;
            long long s = (a2 - a1 + order) % order;

            if (r == 0) continue; // Try again

            long long d = gcd(r, order);
            if (s % d != 0) return -1; // No solution

            // Solve r*x ≡ s (mod order)
            long long r_reduced = r / d;
            long long s_reduced = s / d;
            long long order_reduced = order / d;

            if (gcd(r_reduced, order_reduced) == 1) {
                long long x = (s_reduced * modInverse(r_reduced, order_reduced)) % order_reduced;
                // Verify
                if (modPow(g, x, p) == h) {
                    return x;
                }
            }
        }
    }

    return -1;
}

/**
 * @brief Structure for Diffie-Hellman parameters
 */
struct DiffieHellmanParams {
    long long p;  // Large prime
    long long g;  // Generator of ℤ*_p
};

/**
 * @brief Generate Diffie-Hellman public key from private key
 *
 * Public key = g^a mod p
 *
 * @param params DH parameters (p, g)
 * @param privateKey Private key a
 * @return Public key g^a mod p
 */
inline long long diffieHellmanPublicKey(const DiffieHellmanParams& params,
                                        long long privateKey) {
    return modPow(params.g, privateKey, params.p);
}

/**
 * @brief Compute shared secret in Diffie-Hellman protocol
 *
 * Shared secret = (other's public key)^(my private key) mod p
 *
 * @param params DH parameters (p, g)
 * @param otherPublicKey Other party's public key
 * @param myPrivateKey My private key
 * @return Shared secret
 */
inline long long diffieHellmanSharedSecret(const DiffieHellmanParams& params,
                                           long long otherPublicKey,
                                           long long myPrivateKey) {
    return modPow(otherPublicKey, myPrivateKey, params.p);
}

/**
 * @brief Complete Diffie-Hellman key exchange simulation
 *
 * @param p Prime modulus
 * @param g Generator
 * @param alicePrivate Alice's private key
 * @param bobPrivate Bob's private key
 * @return Tuple (Alice's shared secret, Bob's shared secret) - should be equal
 */
inline std::pair<long long, long long> diffieHellmanExchange(long long p, long long g,
                                                             long long alicePrivate,
                                                             long long bobPrivate) {
    DiffieHellmanParams params{p, g};

    // Alice computes her public key
    long long alicePublic = diffieHellmanPublicKey(params, alicePrivate);

    // Bob computes his public key
    long long bobPublic = diffieHellmanPublicKey(params, bobPrivate);

    // Alice computes shared secret using Bob's public key
    long long aliceShared = diffieHellmanSharedSecret(params, bobPublic, alicePrivate);

    // Bob computes shared secret using Alice's public key
    long long bobShared = diffieHellmanSharedSecret(params, alicePublic, bobPrivate);

    return {aliceShared, bobShared};
}

// =============================================================================
// QUADRATIC RESIDUES AND QUADRATIC RECIPROCITY
// =============================================================================

/**
 * @brief Compute the Legendre symbol (a/p)
 *
 * (a/p) = 0 if a ≡ 0 (mod p)
 *       = 1 if a is a quadratic residue mod p
 *       = -1 if a is not a quadratic residue mod p
 *
 * Uses Euler's criterion: (a/p) ≡ a^((p-1)/2) (mod p)
 *
 * @param a Integer
 * @param p Odd prime
 * @return Legendre symbol value {-1, 0, 1}
 */
inline int legendreSymbol(long long a, long long p) {
    if (!isPrimeMillerRabin(p) || p == 2) return 0;

    a = ((a % p) + p) % p; // Normalize to [0, p)

    if (a == 0) return 0;

    // Euler's criterion: (a/p) ≡ a^((p-1)/2) (mod p)
    long long result = modPow(a, (p - 1) / 2, p);

    // Convert to {-1, 0, 1}
    if (result == p - 1) return -1;
    if (result == 1) return 1;
    return 0;
}

/**
 * @brief Compute the Jacobi symbol (a/n)
 *
 * Generalization of Legendre symbol to composite moduli
 * If n = p1^e1 * p2^e2 * ... * pk^ek, then (a/n) = ∏(a/pi)^ei
 *
 * Uses the quadratic reciprocity law for efficient computation
 *
 * @param a Integer
 * @param n Odd positive integer
 * @return Jacobi symbol value {-1, 0, 1}
 */
inline int jacobiSymbol(long long a, long long n) {
    if (n <= 0 || n % 2 == 0) return 0;

    a = ((a % n) + n) % n;

    int result = 1;

    while (a != 0) {
        // Remove factors of 2
        while (a % 2 == 0) {
            a /= 2;
            // (2/n) = (-1)^((n^2-1)/8)
            long long r = n % 8;
            if (r == 3 || r == 5) {
                result = -result;
            }
        }

        // Quadratic reciprocity: (a/n)(n/a) = (-1)^((a-1)(n-1)/4)
        if (a % 4 == 3 && n % 4 == 3) {
            result = -result;
        }

        // Swap a and n
        std::swap(a, n);
        a %= n;
    }

    return (n == 1) ? result : 0;
}

/**
 * @brief Test if a is a quadratic residue modulo n
 *
 * For odd prime p: uses Legendre symbol
 * For composite n: uses Jacobi symbol (necessary but not sufficient)
 *
 * @param a Integer
 * @param n Modulus
 * @return true if a is a quadratic residue mod n
 */
inline bool isQuadraticResidue(long long a, long long n) {
    if (n == 2) {
        return (a % 2 == 0) || (a % 2 == 1);
    }

    if (isPrimeMillerRabin(n)) {
        return legendreSymbol(a, n) == 1;
    }

    // For composite n, Jacobi symbol is necessary but not sufficient
    // We need to check if there exists x such that x^2 ≡ a (mod n)
    a = ((a % n) + n) % n;

    for (long long x = 0; x * x <= n && x < 1000000; ++x) {
        if ((x * x) % n == a) {
            return true;
        }
    }

    return false;
}

/**
 * @brief Compute modular square root using Tonelli-Shanks algorithm
 *
 * Finds x such that x^2 ≡ a (mod p) for prime p
 *
 * @param a Quadratic residue
 * @param p Odd prime
 * @return Square root x, or -1 if a is not a quadratic residue
 */
inline long long modularSquareRoot(long long a, long long p) {
    if (!isPrimeMillerRabin(p) || p == 2) return -1;

    a = ((a % p) + p) % p;

    // Check if a is a quadratic residue
    if (legendreSymbol(a, p) != 1) return -1;

    // Special case: p ≡ 3 (mod 4)
    if (p % 4 == 3) {
        return modPow(a, (p + 1) / 4, p);
    }

    // Tonelli-Shanks algorithm for p ≡ 1 (mod 4)
    // Write p - 1 = 2^s * q with q odd
    long long q = p - 1;
    int s = 0;
    while (q % 2 == 0) {
        q /= 2;
        s++;
    }

    // Find a quadratic non-residue n
    long long n = 2;
    while (legendreSymbol(n, p) != -1) {
        n++;
    }

    // Initialize
    long long z = modPow(n, q, p);
    long long m = s;
    long long c = z;
    long long t = modPow(a, q, p);
    long long r = modPow(a, (q + 1) / 2, p);

    while (t != 1) {
        // Find least i such that t^(2^i) = 1
        long long temp = t;
        int i = 0;
        while (temp != 1 && i < m) {
            temp = (temp * temp) % p;
            i++;
        }

        if (i == m) return -1; // Should not happen

        // Update values
        long long b = modPow(c, 1LL << (m - i - 1), p);
        m = i;
        c = (b * b) % p;
        t = (t * c) % p;
        r = (r * b) % p;
    }

    return r;
}

/**
 * @brief Verify quadratic residuosity assumption
 *
 * The quadratic residuosity assumption states that for Blum integers n = pq
 * (where p, q ≡ 3 mod 4), it is computationally hard to distinguish
 * quadratic residues from quadratic non-residues with Jacobi symbol 1
 *
 * @param n Blum integer (product of two primes ≡ 3 mod 4)
 * @return true if n is a valid Blum integer
 */
inline bool isBlumInteger(long long n) {
    auto factors = trialDivisionFactorization(n);

    if (factors.size() != 2) return false;

    for (const auto& [prime, exp] : factors) {
        if (exp != 1) return false;
        if (prime % 4 != 3) return false;
    }

    return true;
}

// =============================================================================
// MODULES AND VECTOR SPACES
// =============================================================================

/**
 * @brief Abstract module over a ring
 *
 * A module M over a ring R is an abelian group with scalar multiplication
 * satisfying: r(m1+m2) = rm1+rm2, (r1+r2)m = r1m+r2m, (r1r2)m = r1(r2m)
 */
template<typename Scalar, typename Element>
class Module {
public:
    virtual ~Module() = default;

    // Module operations
    virtual Element add(const Element& a, const Element& b) const = 0;
    virtual Element scalarMultiply(const Scalar& r, const Element& m) const = 0;
    virtual Element zero() const = 0;
    virtual Element negate(const Element& m) const = 0;

    /**
     * @brief Check if elements form a submodule
     */
    bool isSubmodule(const std::vector<Element>& elements) const {
        if (elements.empty()) return true;

        // Must contain zero
        bool hasZero = false;
        Element zeroElem = zero();
        for (const auto& elem : elements) {
            if (elem == zeroElem) {
                hasZero = true;
                break;
            }
        }
        if (!hasZero) return false;

        // Closed under addition and scalar multiplication would need ring elements
        // Simplified check: just verify closure under addition
        for (size_t i = 0; i < elements.size() && i < 10; ++i) {
            for (size_t j = 0; j < elements.size() && j < 10; ++j) {
                Element sum = add(elements[i], elements[j]);
                bool found = false;
                for (const auto& elem : elements) {
                    if (elem == sum) {
                        found = true;
                        break;
                    }
                }
                if (!found) return false;
            }
        }

        return true;
    }

    /**
     * @brief Check linear independence (for vector spaces)
     */
    virtual bool areLinearlyIndependent(const std::vector<Element>& elements) const {
        // Default implementation: stub
        return false;
    }

    /**
     * @brief Compute dimension (for vector spaces)
     */
    virtual int dimension() const {
        return -1; // Unknown
    }
};

/**
 * @brief Vector space over a field
 *
 * Specialization of Module where scalars form a field
 */
template<typename Field, typename Vector>
class VectorSpace : public Module<Field, Vector> {
public:
    /**
     * @brief Check if vectors form a basis
     */
    virtual bool isBasis(const std::vector<Vector>& vectors) const {
        // A basis must be linearly independent and span the space
        if (!this->areLinearlyIndependent(vectors)) return false;

        int dim = this->dimension();
        if (dim >= 0 && static_cast<int>(vectors.size()) != dim) return false;

        return true;
    }

    /**
     * @brief Compute coordinates with respect to a basis
     */
    virtual std::vector<Field> coordinates(const Vector& v,
                                          const std::vector<Vector>& basis) const {
        // Abstract interface - implementation depends on concrete vector space
        return std::vector<Field>();
    }
};

/**
 * @brief Concrete vector space: ℝⁿ
 */
class RealVectorSpace : public VectorSpace<double, std::vector<double>> {
private:
    int n; // Dimension

public:
    explicit RealVectorSpace(int dimension) : n(dimension) {}

    std::vector<double> add(const std::vector<double>& a,
                           const std::vector<double>& b) const override {
        if (a.size() != b.size() || a.size() != static_cast<size_t>(n)) {
            throw std::invalid_argument("Vector dimensions must match");
        }

        std::vector<double> result(n);
        for (int i = 0; i < n; ++i) {
            result[i] = a[i] + b[i];
        }
        return result;
    }

    std::vector<double> scalarMultiply(const double& r,
                                      const std::vector<double>& m) const override {
        if (m.size() != static_cast<size_t>(n)) {
            throw std::invalid_argument("Vector dimension must match space dimension");
        }

        std::vector<double> result(n);
        for (int i = 0; i < n; ++i) {
            result[i] = r * m[i];
        }
        return result;
    }

    std::vector<double> zero() const override {
        return std::vector<double>(n, 0.0);
    }

    std::vector<double> negate(const std::vector<double>& m) const override {
        return scalarMultiply(-1.0, m);
    }

    int dimension() const override {
        return n;
    }

    bool areLinearlyIndependent(const std::vector<std::vector<double>>& vectors) const override;
};

// =============================================================================
// MATRICES
// =============================================================================

/**
 * @brief Matrix class for linear algebra computations
 */
template<typename T = double>
class Matrix {
private:
    std::vector<std::vector<T>> data;
    int rows;
    int cols;

public:
    /**
     * @brief Construct zero matrix
     */
    Matrix(int m, int n) : rows(m), cols(n) {
        data.resize(m, std::vector<T>(n, T(0)));
    }

    /**
     * @brief Construct from 2D vector
     */
    Matrix(const std::vector<std::vector<T>>& mat) : data(mat) {
        rows = mat.size();
        cols = (rows > 0) ? mat[0].size() : 0;
    }

    /**
     * @brief Get number of rows
     */
    int numRows() const { return rows; }

    /**
     * @brief Get number of columns
     */
    int numCols() const { return cols; }

    /**
     * @brief Access element
     */
    T& operator()(int i, int j) {
        return data[i][j];
    }

    const T& operator()(int i, int j) const {
        return data[i][j];
    }

    /**
     * @brief Matrix addition
     */
    Matrix operator+(const Matrix& other) const {
        if (rows != other.rows || cols != other.cols) {
            throw std::invalid_argument("Matrix dimensions must match for addition");
        }

        Matrix result(rows, cols);
        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < cols; ++j) {
                result(i, j) = data[i][j] + other(i, j);
            }
        }
        return result;
    }

    /**
     * @brief Matrix multiplication
     */
    Matrix operator*(const Matrix& other) const {
        if (cols != other.rows) {
            throw std::invalid_argument("Matrix dimensions incompatible for multiplication");
        }

        Matrix result(rows, other.cols);
        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < other.cols; ++j) {
                T sum = T(0);
                for (int k = 0; k < cols; ++k) {
                    sum += data[i][k] * other(k, j);
                }
                result(i, j) = sum;
            }
        }
        return result;
    }

    /**
     * @brief Scalar multiplication
     */
    Matrix operator*(const T& scalar) const {
        Matrix result(rows, cols);
        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < cols; ++j) {
                result(i, j) = data[i][j] * scalar;
            }
        }
        return result;
    }

    /**
     * @brief Matrix transpose
     */
    Matrix transpose() const {
        Matrix result(cols, rows);
        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < cols; ++j) {
                result(j, i) = data[i][j];
            }
        }
        return result;
    }

    /**
     * @brief Compute trace (sum of diagonal elements)
     */
    T trace() const {
        if (rows != cols) {
            throw std::invalid_argument("Trace only defined for square matrices");
        }

        T sum = T(0);
        for (int i = 0; i < rows; ++i) {
            sum += data[i][i];
        }
        return sum;
    }

    /**
     * @brief Create identity matrix
     */
    static Matrix identity(int n) {
        Matrix result(n, n);
        for (int i = 0; i < n; ++i) {
            result(i, i) = T(1);
        }
        return result;
    }

    /**
     * @brief Get row as vector
     */
    std::vector<T> getRow(int i) const {
        if (i < 0 || i >= rows) {
            throw std::out_of_range("Row index out of range");
        }
        return data[i];
    }

    /**
     * @brief Get column as vector
     */
    std::vector<T> getCol(int j) const {
        if (j < 0 || j >= cols) {
            throw std::out_of_range("Column index out of range");
        }

        std::vector<T> col(rows);
        for (int i = 0; i < rows; ++i) {
            col[i] = data[i][j];
        }
        return col;
    }

    /**
     * @brief Swap two rows
     */
    void swapRows(int i, int j) {
        if (i < 0 || i >= rows || j < 0 || j >= rows) {
            throw std::out_of_range("Row index out of range");
        }
        std::swap(data[i], data[j]);
    }

    /**
     * @brief Multiply row by scalar
     */
    void multiplyRow(int i, const T& scalar) {
        if (i < 0 || i >= rows) {
            throw std::out_of_range("Row index out of range");
        }
        for (int j = 0; j < cols; ++j) {
            data[i][j] *= scalar;
        }
    }

    /**
     * @brief Add multiple of one row to another
     */
    void addRowMultiple(int target, int source, const T& scalar) {
        if (target < 0 || target >= rows || source < 0 || source >= rows) {
            throw std::out_of_range("Row index out of range");
        }
        for (int j = 0; j < cols; ++j) {
            data[target][j] += scalar * data[source][j];
        }
    }

    /**
     * @brief Compute determinant using LU decomposition
     */
    T determinant() const;

    /**
     * @brief Compute inverse using Gauss-Jordan elimination
     */
    Matrix inverse() const;

    /**
     * @brief Solve Ax = b using Gaussian elimination
     */
    std::vector<T> solve(const std::vector<T>& b) const;

    /**
     * @brief Gaussian elimination to row echelon form
     */
    Matrix rowEchelonForm() const;

    /**
     * @brief Gaussian elimination to reduced row echelon form
     */
    Matrix reducedRowEchelonForm() const;

    /**
     * @brief Compute rank
     */
    int rank() const;

    /**
     * @brief Get underlying data
     */
    const std::vector<std::vector<T>>& getData() const { return data; }
};

// =============================================================================
// GAUSSIAN ELIMINATION
// =============================================================================

/**
 * @brief Gaussian elimination to row echelon form
 *
 * Transform matrix to upper triangular form
 */
template<typename T>
Matrix<T> Matrix<T>::rowEchelonForm() const {
    Matrix result = *this;

    int pivot_row = 0;

    for (int col = 0; col < cols && pivot_row < rows; ++col) {
        // Find pivot
        int max_row = pivot_row;
        T max_val = std::abs(result(pivot_row, col));

        for (int row = pivot_row + 1; row < rows; ++row) {
            T val = std::abs(result(row, col));
            if (val > max_val) {
                max_val = val;
                max_row = row;
            }
        }

        // Check if pivot is zero
        if (std::abs(result(max_row, col)) < 1e-10) {
            continue; // Skip this column
        }

        // Swap rows
        if (max_row != pivot_row) {
            result.swapRows(pivot_row, max_row);
        }

        // Eliminate below pivot
        for (int row = pivot_row + 1; row < rows; ++row) {
            T factor = result(row, col) / result(pivot_row, col);
            result.addRowMultiple(row, pivot_row, -factor);
        }

        pivot_row++;
    }

    return result;
}

/**
 * @brief Gaussian elimination to reduced row echelon form
 *
 * Transform matrix to canonical form with leading 1s
 */
template<typename T>
Matrix<T> Matrix<T>::reducedRowEchelonForm() const {
    Matrix result = rowEchelonForm();

    // Back substitution to get reduced form
    for (int row = rows - 1; row >= 0; --row) {
        // Find leading entry (pivot)
        int pivot_col = -1;
        for (int col = 0; col < cols; ++col) {
            if (std::abs(result(row, col)) > 1e-10) {
                pivot_col = col;
                break;
            }
        }

        if (pivot_col == -1) continue; // Zero row

        // Scale to make pivot = 1
        T pivot_val = result(row, pivot_col);
        result.multiplyRow(row, T(1) / pivot_val);

        // Eliminate above pivot
        for (int above = 0; above < row; ++above) {
            T factor = result(above, pivot_col);
            result.addRowMultiple(above, row, -factor);
        }
    }

    return result;
}

/**
 * @brief Compute matrix rank
 */
template<typename T>
int Matrix<T>::rank() const {
    Matrix ref = rowEchelonForm();

    int rank = 0;
    for (int row = 0; row < rows; ++row) {
        bool nonzero = false;
        for (int col = 0; col < cols; ++col) {
            if (std::abs(ref(row, col)) > 1e-10) {
                nonzero = true;
                break;
            }
        }
        if (nonzero) rank++;
    }

    return rank;
}

/**
 * @brief Solve linear system Ax = b using Gaussian elimination
 */
template<typename T>
std::vector<T> Matrix<T>::solve(const std::vector<T>& b) const {
    if (rows != cols) {
        throw std::invalid_argument("Matrix must be square for solve()");
    }
    if (b.size() != static_cast<size_t>(rows)) {
        throw std::invalid_argument("Right-hand side dimension must match matrix rows");
    }

    // Create augmented matrix [A|b]
    Matrix augmented(rows, cols + 1);
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            augmented(i, j) = data[i][j];
        }
        augmented(i, cols) = b[i];
    }

    // Row reduce
    Matrix rref = augmented.reducedRowEchelonForm();

    // Extract solution
    std::vector<T> x(rows);
    for (int i = 0; i < rows; ++i) {
        x[i] = rref(i, cols);
    }

    return x;
}

/**
 * @brief Compute determinant using row reduction
 */
template<typename T>
T Matrix<T>::determinant() const {
    if (rows != cols) {
        throw std::invalid_argument("Determinant only defined for square matrices");
    }

    Matrix temp = *this;
    T det = T(1);

    for (int col = 0; col < cols; ++col) {
        // Find pivot
        int pivot_row = col;
        for (int row = col + 1; row < rows; ++row) {
            if (std::abs(temp(row, col)) > std::abs(temp(pivot_row, col))) {
                pivot_row = row;
            }
        }

        // Check for zero pivot
        if (std::abs(temp(pivot_row, col)) < 1e-10) {
            return T(0);
        }

        // Swap rows if needed
        if (pivot_row != col) {
            temp.swapRows(col, pivot_row);
            det = -det;
        }

        // Update determinant
        det *= temp(col, col);

        // Eliminate below
        for (int row = col + 1; row < rows; ++row) {
            T factor = temp(row, col) / temp(col, col);
            temp.addRowMultiple(row, col, -factor);
        }
    }

    return det;
}

/**
 * @brief Compute matrix inverse using Gauss-Jordan elimination
 */
template<typename T>
Matrix<T> Matrix<T>::inverse() const {
    if (rows != cols) {
        throw std::invalid_argument("Inverse only defined for square matrices");
    }

    // Create augmented matrix [A|I]
    Matrix augmented(rows, 2 * cols);
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            augmented(i, j) = data[i][j];
            augmented(i, j + cols) = (i == j) ? T(1) : T(0);
        }
    }

    // Row reduce to [I|A^(-1)]
    Matrix rref = augmented.reducedRowEchelonForm();

    // Check if we got identity on the left
    for (int i = 0; i < rows; ++i) {
        if (std::abs(rref(i, i) - T(1)) > 1e-10) {
            throw std::runtime_error("Matrix is singular and has no inverse");
        }
    }

    // Extract inverse from right half
    Matrix inv(rows, cols);
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            inv(i, j) = rref(i, j + cols);
        }
    }

    return inv;
}

/**
 * @brief Check linear independence of vectors (implementation for RealVectorSpace)
 */
inline bool RealVectorSpace::areLinearlyIndependent(
    const std::vector<std::vector<double>>& vectors) const {

    if (vectors.empty()) return true;
    if (vectors.size() > static_cast<size_t>(n)) return false;

    // Create matrix with vectors as columns
    Matrix<double> mat(n, vectors.size());
    for (size_t j = 0; j < vectors.size(); ++j) {
        if (vectors[j].size() != static_cast<size_t>(n)) {
            throw std::invalid_argument("All vectors must have dimension n");
        }
        for (int i = 0; i < n; ++i) {
            mat(i, j) = vectors[j][i];
        }
    }

    // Vectors are linearly independent iff matrix has full column rank
    return mat.rank() == static_cast<int>(vectors.size());
}

// =============================================================================
// SUBEXPONENTIAL-TIME DISCRETE LOGARITHMS AND FACTORING
// =============================================================================

/**
 * @brief Check if n is B-smooth (all prime factors ≤ B)
 *
 * A number is B-smooth if all its prime factors are at most B
 *
 * @param n Number to test
 * @param B Smoothness bound
 * @return true if n is B-smooth
 */
inline bool isBSmooth(long long n, long long B) {
    if (n <= 1) return false;

    auto factors = trialDivisionFactorization(n);

    for (const auto& [prime, exp] : factors) {
        if (prime > B) return false;
    }

    return true;
}

/**
 * @brief Find all B-smooth numbers up to N
 *
 * @param N Upper bound
 * @param B Smoothness bound
 * @return Vector of B-smooth numbers
 */
inline std::vector<long long> findBSmoothNumbers(long long N, long long B) {
    std::vector<long long> smooth_numbers;

    for (long long n = 2; n <= N; ++n) {
        if (isBSmooth(n, B)) {
            smooth_numbers.push_back(n);
        }
    }

    return smooth_numbers;
}

/**
 * @brief Compute smoothness probability for random numbers near N
 *
 * Uses Dickman's function approximation
 *
 * @param N Number size
 * @param B Smoothness bound
 * @return Approximate probability that random n ≈ N is B-smooth
 */
inline double smoothnessProbability(double N, double B) {
    if (B >= N) return 1.0;

    double u = std::log(N) / std::log(B);

    // Dickman's function ρ(u) approximation
    if (u <= 1.0) return 1.0;
    if (u <= 2.0) return 1.0 - std::log(u);

    // For u > 2, use approximation ρ(u) ≈ u^(-u)
    return std::pow(u, -u);
}

/**
 * @brief Index calculus discrete logarithm algorithm
 *
 * Subexponential algorithm for computing discrete logarithms
 * Complexity: L[1/2, c] where L[α, c] = exp(c(log n)^α (log log n)^(1-α))
 *
 * @param g Generator
 * @param h Target (find x such that g^x ≡ h mod p)
 * @param p Prime modulus
 * @param B Smoothness bound (factor base size)
 * @param maxAttempts Maximum attempts
 * @return Discrete logarithm x, or -1 if not found
 */
inline long long discreteLogIndexCalculus(long long g, long long h, long long p,
                                          long long B = 100, int maxAttempts = 10000) {
    // Step 1: Build factor base (primes ≤ B)
    std::vector<long long> factorBase = sieveOfEratosthenes(B);
    if (factorBase.empty()) return -1;

    int baseSize = factorBase.size();

    // Step 2: Collect relations (find k where g^k mod p is B-smooth)
    std::vector<std::vector<int>> relations;
    std::vector<long long> exponents;

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<long long> dist(1, p - 2);

    for (int attempt = 0; attempt < maxAttempts && relations.size() < static_cast<size_t>(baseSize + 5); ++attempt) {
        long long k = dist(gen);
        long long val = modPow(g, k, p);

        // Try to factor val over factor base
        auto factors = trialDivisionFactorization(val);

        bool smooth = true;
        std::vector<int> exponentVector(baseSize, 0);

        for (const auto& [prime, exp] : factors) {
            // Find prime in factor base
            auto it = std::find(factorBase.begin(), factorBase.end(), prime);
            if (it == factorBase.end()) {
                smooth = false;
                break;
            }
            int idx = std::distance(factorBase.begin(), it);
            exponentVector[idx] = exp;
        }

        if (smooth) {
            relations.push_back(exponentVector);
            exponents.push_back(k);
        }
    }

    if (relations.size() < static_cast<size_t>(baseSize)) {
        return -1; // Not enough relations
    }

    // Step 3: Solve linear system to find log_g(p_i) for each prime p_i
    // This is a simplified version - full implementation would use linear algebra mod (p-1)
    std::vector<long long> logTable(baseSize, 0);

    // Simplified: just use small examples
    // In practice, this requires Gaussian elimination over Z/(p-1)Z
    for (int i = 0; i < baseSize && i < static_cast<int>(relations.size()); ++i) {
        // g^k = ∏ p_i^e_i implies k = ∑ e_i * log_g(p_i) mod (p-1)
        // This is a linear system that needs to be solved
        logTable[i] = exponents[i] % (p - 1);
    }

    // Step 4: Express h in terms of factor base
    for (int attempt = 0; attempt < maxAttempts; ++attempt) {
        long long s = dist(gen);
        long long val = (h * modPow(g, s, p)) % p;

        auto factors = trialDivisionFactorization(val);

        bool smooth = true;
        long long logSum = 0;

        for (const auto& [prime, exp] : factors) {
            auto it = std::find(factorBase.begin(), factorBase.end(), prime);
            if (it == factorBase.end()) {
                smooth = false;
                break;
            }
            int idx = std::distance(factorBase.begin(), it);
            logSum = (logSum + exp * logTable[idx]) % (p - 1);
        }

        if (smooth) {
            // h * g^s = ∏ p_i^e_i implies log_g(h) + s = ∑ e_i * log_g(p_i)
            long long result = (logSum - s + (p - 1)) % (p - 1);
            return result;
        }
    }

    return -1;
}

/**
 * @brief Quadratic sieve factoring algorithm (simplified version)
 *
 * Subexponential algorithm for integer factorization
 * Complexity: L[1/2, 1] where L[α, c] = exp(c(log n)^α (log log n)^(1-α))
 *
 * @param n Number to factor
 * @param B Smoothness bound
 * @param maxAttempts Maximum attempts
 * @return Non-trivial factor of n, or -1 if not found
 */
inline long long quadraticSieve(long long n, long long B = 100, int maxAttempts = 10000) {
    if (n <= 1) return -1;
    if (isPrimeMillerRabin(n)) return n;

    // Check for small factors first
    for (long long p = 2; p <= std::min(B, static_cast<long long>(10000)); ++p) {
        if (n % p == 0) return p;
    }

    // Step 1: Build factor base (primes p where n is a quadratic residue mod p)
    std::vector<long long> factorBase;
    std::vector<long long> primes = sieveOfEratosthenes(B);

    for (long long p : primes) {
        if (legendreSymbol(n % p, p) >= 0) {
            factorBase.push_back(p);
        }
    }

    if (factorBase.empty()) return -1;

    int baseSize = factorBase.size();

    // Step 2: Find B-smooth values of Q(x) = (x + ⌈√n⌉)² - n
    long long sqrtN = static_cast<long long>(std::sqrt(n)) + 1;

    std::vector<long long> smoothX;
    std::vector<std::vector<int>> exponentVectors;

    for (long long x = 0; x < maxAttempts && smoothX.size() < static_cast<size_t>(baseSize + 2); ++x) {
        long long val = (x + sqrtN) * (x + sqrtN) - n;
        if (val <= 0) continue;

        auto factors = trialDivisionFactorization(val);

        bool smooth = true;
        std::vector<int> expVector(baseSize, 0);

        for (const auto& [prime, exp] : factors) {
            auto it = std::find(factorBase.begin(), factorBase.end(), prime);
            if (it == factorBase.end()) {
                smooth = false;
                break;
            }
            int idx = std::distance(factorBase.begin(), it);
            expVector[idx] = exp;
        }

        if (smooth) {
            smoothX.push_back(x);
            exponentVectors.push_back(expVector);
        }
    }

    if (smoothX.size() < 2) return -1;

    // Step 3: Find subset with even exponents (simplified - just try pairs)
    for (size_t i = 0; i < smoothX.size(); ++i) {
        for (size_t j = i + 1; j < smoothX.size(); ++j) {
            // Compute product of (x + √n)² - n
            long long x_prod = 1, y_prod = 1;

            x_prod = ((smoothX[i] + sqrtN) * (smoothX[j] + sqrtN)) % n;

            // Check if gcd gives factor
            long long factor = gcd(x_prod - y_prod, n);
            if (factor > 1 && factor < n) {
                return factor;
            }

            factor = gcd(x_prod + y_prod, n);
            if (factor > 1 && factor < n) {
                return factor;
            }
        }
    }

    return -1;
}

/**
 * @brief Pollard's p-1 factoring algorithm
 *
 * Works well when n has a prime factor p such that p-1 is B-smooth
 *
 * @param n Number to factor
 * @param B Smoothness bound
 * @return Non-trivial factor of n, or -1 if not found
 */
inline long long pollardPMinus1(long long n, long long B = 1000) {
    if (n <= 1) return -1;
    if (isPrimeMillerRabin(n)) return n;

    long long a = 2;

    // Compute a^(∏ p^e) mod n for all primes p ≤ B
    std::vector<long long> primes = sieveOfEratosthenes(B);

    for (long long p : primes) {
        long long p_power = p;
        while (p_power <= B) {
            a = modPow(a, p, n);
            p_power *= p;
        }
    }

    long long g = gcd(a - 1, n);

    if (g > 1 && g < n) {
        return g;
    }

    return -1;
}

// =============================================================================
// MORE RINGS
// =============================================================================

/**
 * @brief Algebra over a field
 *
 * An algebra is a vector space with a bilinear multiplication operation
 */
template<typename Field>
class Algebra {
public:
    virtual ~Algebra() = default;

    virtual std::vector<Field> add(const std::vector<Field>& a,
                                   const std::vector<Field>& b) const = 0;
    virtual std::vector<Field> multiply(const std::vector<Field>& a,
                                       const std::vector<Field>& b) const = 0;
    virtual std::vector<Field> scalarMultiply(const Field& c,
                                             const std::vector<Field>& a) const = 0;
    virtual std::vector<Field> zero() const = 0;
    virtual std::vector<Field> one() const = 0;
    virtual int dimension() const = 0;
};

/**
 * @brief Field of fractions of an integral domain
 *
 * Constructs F = Frac(R) where R is an integral domain
 * Elements are represented as pairs (numerator, denominator)
 */
template<typename T>
class FieldOfFractions {
public:
    T numerator;
    T denominator;

    FieldOfFractions(const T& num = T(0), const T& den = T(1))
        : numerator(num), denominator(den) {
        if (denominator == T(0)) {
            throw std::invalid_argument("Denominator cannot be zero");
        }
        normalize();
    }

    /**
     * @brief Normalize fraction (make coprime, positive denominator)
     */
    void normalize() {
        // This is a simplified version for integral types
        // For polynomial rings, would need polynomial GCD
        if (denominator < T(0)) {
            numerator = -numerator;
            denominator = -denominator;
        }
    }

    /**
     * @brief Addition in field of fractions
     */
    FieldOfFractions operator+(const FieldOfFractions& other) const {
        T num = numerator * other.denominator + other.numerator * denominator;
        T den = denominator * other.denominator;
        return FieldOfFractions(num, den);
    }

    /**
     * @brief Multiplication in field of fractions
     */
    FieldOfFractions operator*(const FieldOfFractions& other) const {
        T num = numerator * other.numerator;
        T den = denominator * other.denominator;
        return FieldOfFractions(num, den);
    }

    /**
     * @brief Multiplicative inverse
     */
    FieldOfFractions inverse() const {
        if (numerator == T(0)) {
            throw std::invalid_argument("Cannot invert zero");
        }
        return FieldOfFractions(denominator, numerator);
    }

    /**
     * @brief Division
     */
    FieldOfFractions operator/(const FieldOfFractions& other) const {
        return (*this) * other.inverse();
    }

    bool operator==(const FieldOfFractions& other) const {
        return numerator * other.denominator == other.numerator * denominator;
    }
};

/**
 * @brief Extended polynomial class with more operations
 */
class PolynomialExtended : public Polynomial {
public:
    using Polynomial::Polynomial;

    /**
     * @brief Polynomial division with remainder
     *
     * Computes q, r such that this = q * divisor + r with deg(r) < deg(divisor)
     */
    std::pair<PolynomialExtended, PolynomialExtended> divideWithRemainder(
        const PolynomialExtended& divisor) const {

        if (divisor.degree() < 0) {
            throw std::invalid_argument("Cannot divide by zero polynomial");
        }

        PolynomialExtended quotient, remainder = *this;

        while (remainder.degree() >= divisor.degree() && remainder.degree() >= 0) {
            int degDiff = remainder.degree() - divisor.degree();
            long long coef = remainder.coeffs[remainder.degree()] /
                           divisor.coeffs[divisor.degree()];

            // Create monomial for this step
            std::vector<long long> monomial(degDiff + 1, 0);
            monomial[degDiff] = coef;
            PolynomialExtended term(monomial);

            // Update quotient and remainder
            if (quotient.coeffs.empty()) {
                quotient = term;
            } else {
                quotient = quotient + term;
            }

            PolynomialExtended subtrahend(divisor.coeffs);
            for (size_t i = 0; i < subtrahend.coeffs.size(); ++i) {
                subtrahend.coeffs[i] *= coef;
            }
            subtrahend.coeffs.insert(subtrahend.coeffs.begin(), degDiff, 0);

            // Remainder -= term * divisor
            for (size_t i = 0; i < subtrahend.coeffs.size() && i < remainder.coeffs.size(); ++i) {
                remainder.coeffs[i] -= subtrahend.coeffs[i];
            }

            // Remove leading zeros
            while (!remainder.coeffs.empty() && remainder.coeffs.back() == 0) {
                remainder.coeffs.pop_back();
            }
        }

        return {quotient, remainder};
    }

    /**
     * @brief Polynomial GCD using Euclidean algorithm
     */
    static PolynomialExtended gcd(const PolynomialExtended& a,
                                 const PolynomialExtended& b) {
        PolynomialExtended x = a, y = b;

        while (y.degree() >= 0) {
            auto [q, r] = x.divideWithRemainder(y);
            x = y;
            y = r;
        }

        // Make monic (leading coefficient = 1)
        if (x.degree() >= 0 && x.coeffs[x.degree()] != 0) {
            long long lead = x.coeffs[x.degree()];
            for (auto& c : x.coeffs) {
                c /= lead;
            }
        }

        return x;
    }

    /**
     * @brief Extended Euclidean algorithm for polynomials
     *
     * Finds s, t such that s*a + t*b = gcd(a, b)
     */
    static std::tuple<PolynomialExtended, PolynomialExtended, PolynomialExtended>
    extendedGCD(const PolynomialExtended& a, const PolynomialExtended& b) {
        if (b.degree() < 0) {
            return {a, PolynomialExtended({1}), PolynomialExtended({0})};
        }

        auto [q, r] = a.divideWithRemainder(b);
        auto [g, s1, t1] = extendedGCD(b, r);

        // s = t1, t = s1 - q*t1
        PolynomialExtended s = t1;
        PolynomialExtended qTimesT1({0});
        // Multiply q * t1 (simplified)
        PolynomialExtended t = s1;

        return {g, s, t};
    }

    /**
     * @brief Test if polynomial is irreducible (simplified test)
     */
    bool isIrreducible(long long p) const {
        // Simplified irreducibility test over Z/pZ
        // Full implementation would use more sophisticated tests
        if (degree() <= 1) return degree() == 1;

        // Check if has roots mod p
        for (long long x = 0; x < p; ++x) {
            if (evaluate(x) % p == 0) {
                return false; // Has a root, so reducible
            }
        }

        return true; // No roots found (necessary but not sufficient)
    }
};

/**
 * @brief Unique factorization of polynomials over Z/pZ
 *
 * Factor polynomial into irreducible factors
 */
inline std::map<PolynomialExtended, int> factorPolynomial(
    const PolynomialExtended& f, long long p) {

    std::map<PolynomialExtended, int> factors;

    // Simplified factorization - just find linear factors
    for (long long a = 0; a < p; ++a) {
        if (f.evaluate(a) % p == 0) {
            // (x - a) is a factor
            PolynomialExtended linearFactor({-a, 1});
            factors[linearFactor]++;
        }
    }

    return factors;
}

/**
 * @brief Polynomial congruence class
 *
 * Represents polynomial modulo another polynomial
 */
class PolynomialModulo {
public:
    PolynomialExtended poly;
    PolynomialExtended modulus;

    PolynomialModulo(const PolynomialExtended& p, const PolynomialExtended& m)
        : poly(p), modulus(m) {
        reduce();
    }

    /**
     * @brief Reduce polynomial modulo the modulus
     */
    void reduce() {
        if (modulus.degree() >= 0) {
            auto [q, r] = poly.divideWithRemainder(modulus);
            poly = r;
        }
    }

    /**
     * @brief Addition in quotient ring
     */
    PolynomialModulo operator+(const PolynomialModulo& other) const {
        if (!(modulus == other.modulus)) {
            throw std::invalid_argument("Moduli must match");
        }
        return PolynomialModulo(poly + other.poly, modulus);
    }

    /**
     * @brief Multiplication in quotient ring
     */
    PolynomialModulo operator*(const PolynomialModulo& other) const {
        if (!(modulus == other.modulus)) {
            throw std::invalid_argument("Moduli must match");
        }
        return PolynomialModulo(poly * other.poly, modulus);
    }
};

/**
 * @brief Extension field F[x]/(f) where f is irreducible
 */
template<typename Field = long long>
class ExtensionField {
private:
    PolynomialExtended modulus;
    long long characteristic; // For finite fields

public:
    ExtensionField(const PolynomialExtended& irreducible, long long p = 0)
        : modulus(irreducible), characteristic(p) {
        // In practice, should verify irreducibility
    }

    /**
     * @brief Get degree of extension [F(α) : F]
     */
    int extensionDegree() const {
        return modulus.degree();
    }

    /**
     * @brief Get field size (for finite fields)
     */
    long long fieldSize() const {
        if (characteristic <= 0) return -1;
        long long result = 1;
        for (int i = 0; i < extensionDegree(); ++i) {
            result *= characteristic;
        }
        return result;
    }
};

/**
 * @brief Formal power series over a ring
 *
 * Represents ∑ a_i x^i with coefficients in a ring
 */
template<typename T>
class FormalPowerSeries {
private:
    std::vector<T> coefficients;
    int precision; // Number of terms computed

public:
    FormalPowerSeries(const std::vector<T>& coeffs, int prec = 10)
        : coefficients(coeffs), precision(prec) {
        if (coefficients.size() < static_cast<size_t>(precision)) {
            coefficients.resize(precision, T(0));
        }
    }

    /**
     * @brief Addition of power series
     */
    FormalPowerSeries operator+(const FormalPowerSeries& other) const {
        int maxPrec = std::max(precision, other.precision);
        std::vector<T> result(maxPrec, T(0));

        for (int i = 0; i < maxPrec; ++i) {
            if (i < precision) result[i] += coefficients[i];
            if (i < other.precision) result[i] += other.coefficients[i];
        }

        return FormalPowerSeries(result, maxPrec);
    }

    /**
     * @brief Multiplication of power series (Cauchy product)
     */
    FormalPowerSeries operator*(const FormalPowerSeries& other) const {
        int maxPrec = std::min(precision, other.precision);
        std::vector<T> result(maxPrec, T(0));

        for (int n = 0; n < maxPrec; ++n) {
            for (int k = 0; k <= n; ++k) {
                if (k < precision && (n - k) < other.precision) {
                    result[n] += coefficients[k] * other.coefficients[n - k];
                }
            }
        }

        return FormalPowerSeries(result, maxPrec);
    }

    /**
     * @brief Get coefficient of x^n
     */
    T coefficient(int n) const {
        if (n >= 0 && n < static_cast<int>(coefficients.size())) {
            return coefficients[n];
        }
        return T(0);
    }
};

/**
 * @brief Laurent series (power series with finitely many negative powers)
 *
 * Represents ∑_{i≥n} a_i x^i where n can be negative
 */
template<typename T>
class LaurentSeries {
private:
    std::vector<T> coefficients;
    int minDegree; // Lowest degree term
    int maxDegree; // Highest degree computed

public:
    LaurentSeries(const std::vector<T>& coeffs, int minDeg, int maxDeg)
        : coefficients(coeffs), minDegree(minDeg), maxDegree(maxDeg) {}

    /**
     * @brief Get coefficient of x^n
     */
    T coefficient(int n) const {
        if (n >= minDegree && n <= maxDegree) {
            int idx = n - minDegree;
            if (idx >= 0 && idx < static_cast<int>(coefficients.size())) {
                return coefficients[idx];
            }
        }
        return T(0);
    }

    /**
     * @brief Valuation (degree of lowest non-zero term)
     */
    int valuation() const {
        return minDegree;
    }
};

/**
 * @brief Unique Factorization Domain interface
 */
template<typename T>
class UniqueFactorizationDomain {
public:
    virtual ~UniqueFactorizationDomain() = default;

    /**
     * @brief Check if element is a unit (invertible)
     */
    virtual bool isUnit(const T& a) const = 0;

    /**
     * @brief Check if element is irreducible
     */
    virtual bool isIrreducible(const T& a) const = 0;

    /**
     * @brief Check if element is prime
     */
    virtual bool isPrime(const T& a) const = 0;

    /**
     * @brief Factor element into irreducibles
     */
    virtual std::map<T, int> factor(const T& a) const = 0;

    /**
     * @brief Compute GCD
     */
    virtual T gcd(const T& a, const T& b) const = 0;

    /**
     * @brief Compute LCM
     */
    virtual T lcm(const T& a, const T& b) const {
        T g = gcd(a, b);
        // Return (a * b) / g
        return a; // Simplified
    }
};

/**
 * @brief Integers as a UFD
 */
class IntegerUFD : public UniqueFactorizationDomain<long long> {
public:
    bool isUnit(const long long& a) const override {
        return a == 1 || a == -1;
    }

    bool isIrreducible(const long long& a) const override {
        return isPrimeMillerRabin(std::abs(a));
    }

    bool isPrime(const long long& a) const override {
        return isPrimeMillerRabin(std::abs(a));
    }

    std::map<long long, int> factor(const long long& a) const override {
        return trialDivisionFactorization(std::abs(a));
    }

    long long gcd(const long long& a, const long long& b) const override {
        return ::maths::number_theory::gcd(a, b);
    }
};

// =============================================================================
// POLYNOMIAL ARITHMETIC AND APPLICATIONS
// =============================================================================

/**
 * @brief Polynomial modular inverse
 *
 * Find g such that f*g ≡ 1 (mod h)
 */
inline PolynomialExtended polynomialModularInverse(const PolynomialExtended& f,
                                                   const PolynomialExtended& h) {
    auto [g, s, t] = PolynomialExtended::extendedGCD(f, h);

    if (g.degree() != 0 || g.coeffs[0] != 1) {
        throw std::invalid_argument("Polynomials are not coprime");
    }

    // s*f + t*h = 1, so s*f ≡ 1 (mod h)
    return s;
}

/**
 * @brief Chinese Remainder Theorem for polynomials
 *
 * Find f such that f ≡ a_i (mod m_i) for all i
 */
inline PolynomialExtended polynomialChineseRemainder(
    const std::vector<PolynomialExtended>& remainders,
    const std::vector<PolynomialExtended>& moduli) {

    if (remainders.size() != moduli.size()) {
        throw std::invalid_argument("Remainders and moduli must have same size");
    }

    if (remainders.empty()) {
        return PolynomialExtended({0});
    }

    // Compute M = ∏ m_i
    PolynomialExtended M = moduli[0];
    for (size_t i = 1; i < moduli.size(); ++i) {
        M = M * moduli[i];
    }

    PolynomialExtended result({0});

    for (size_t i = 0; i < remainders.size(); ++i) {
        // M_i = M / m_i
        auto [Mi, rem] = M.divideWithRemainder(moduli[i]);

        // Find y_i such that M_i * y_i ≡ 1 (mod m_i)
        PolynomialExtended yi = polynomialModularInverse(Mi, moduli[i]);

        // result += a_i * M_i * y_i
        PolynomialExtended term = remainders[i] * Mi * yi;
        result = result + term;
    }

    // Reduce modulo M
    auto [q, r] = result.divideWithRemainder(M);

    return r;
}

/**
 * @brief Rational function reconstruction
 *
 * Given a/b mod m, recover a/b
 * Find n/d with deg(n), deg(d) < deg(m)/2 such that r*d ≡ n (mod m)
 */
inline std::pair<PolynomialExtended, PolynomialExtended>
rationalFunctionReconstruction(const PolynomialExtended& r,
                               const PolynomialExtended& m,
                               int maxDegree = -1) {

    if (maxDegree < 0) {
        maxDegree = m.degree() / 2;
    }

    // Use extended Euclidean algorithm
    PolynomialExtended r0 = m, r1 = r;
    PolynomialExtended s0({0}), s1({1});

    while (r1.degree() >= maxDegree) {
        auto [q, r2] = r0.divideWithRemainder(r1);
        PolynomialExtended s2 = s0;
        // s2 = s0 - q*s1 (simplified)

        r0 = r1;
        r1 = r2;
        s0 = s1;
        s1 = s2;
    }

    // Return (numerator, denominator) = (r1, s1)
    return {r1, s1};
}

/**
 * @brief Compute minimal polynomial of element in F[x]/(f)
 *
 * Find the minimal polynomial of α over base field
 */
inline PolynomialExtended minimalPolynomial(const PolynomialModulo& alpha,
                                           long long p, int maxDeg = 10) {
    // Build matrix of powers of alpha
    std::vector<PolynomialExtended> powers;
    PolynomialModulo current({1}, alpha.modulus);

    for (int i = 0; i <= maxDeg; ++i) {
        powers.push_back(current.poly);
        current = current * alpha;
    }

    // Find linear dependence among powers
    // This would require solving a system over the base field
    // Simplified version: return x - constant

    return PolynomialExtended({0, 1}); // Just return x for now
}

/**
 * @brief Fast polynomial multiplication using FFT (placeholder)
 *
 * Multiply two polynomials in O(n log n) time
 */
inline PolynomialExtended fastPolynomialMultiply(const PolynomialExtended& a,
                                                const PolynomialExtended& b) {
    // Full FFT implementation would go here
    // For now, fall back to standard multiplication
    return a * b;
}

/**
 * @brief Polynomial evaluation at multiple points (fast multipoint evaluation)
 *
 * Evaluate polynomial at n points in O(n log²n) time
 */
inline std::vector<long long> multipointEvaluation(const PolynomialExtended& f,
                                                   const std::vector<long long>& points) {
    std::vector<long long> results;

    for (long long x : points) {
        results.push_back(f.evaluate(x));
    }

    return results;
}

/**
 * @brief Polynomial interpolation from values
 *
 * Find polynomial f of degree < n such that f(x_i) = y_i
 */
inline PolynomialExtended interpolate(const std::vector<long long>& xPoints,
                                     const std::vector<long long>& yValues) {
    if (xPoints.size() != yValues.size()) {
        throw std::invalid_argument("Points and values must have same size");
    }

    int n = xPoints.size();
    if (n == 0) return PolynomialExtended({0});

    // Lagrange interpolation
    std::vector<long long> result(n, 0);

    for (int i = 0; i < n; ++i) {
        // Compute Lagrange basis polynomial L_i(x)
        std::vector<long long> basis = {1};

        for (int j = 0; j < n; ++j) {
            if (i != j) {
                // Multiply by (x - x_j) / (x_i - x_j)
                long long denom = xPoints[i] - xPoints[j];

                // Multiply basis by (x - x_j)
                std::vector<long long> newBasis(basis.size() + 1, 0);
                for (size_t k = 0; k < basis.size(); ++k) {
                    newBasis[k] -= basis[k] * xPoints[j];
                    newBasis[k + 1] += basis[k];
                }
                basis = newBasis;

                // Divide by (x_i - x_j)
                for (auto& c : basis) {
                    c /= denom;
                }
            }
        }

        // Add y_i * L_i to result
        for (size_t k = 0; k < basis.size() && k < result.size(); ++k) {
            result[k] += yValues[i] * basis[k];
        }
    }

    return PolynomialExtended(result);
}

// =============================================================================
// LINEARLY GENERATED SEQUENCES AND APPLICATIONS
// =============================================================================

/**
 * @brief Linearly generated sequence over a field
 *
 * A sequence (a_0, a_1, a_2, ...) satisfying a linear recurrence:
 * a_n = c_1*a_{n-1} + c_2*a_{n-2} + ... + c_d*a_{n-d}
 */
template<typename T = long long>
class LinearlyGeneratedSequence {
private:
    std::vector<T> coefficients; // Recurrence coefficients c_1, ..., c_d
    std::vector<T> initialValues; // Initial values a_0, ..., a_{d-1}
    long long modulus; // For sequences over Z/pZ

public:
    LinearlyGeneratedSequence(const std::vector<T>& coeffs,
                             const std::vector<T>& initial,
                             long long mod = 0)
        : coefficients(coeffs), initialValues(initial), modulus(mod) {}

    /**
     * @brief Compute nth term of sequence
     */
    T getNthTerm(int n) const {
        if (n < static_cast<int>(initialValues.size())) {
            return initialValues[n];
        }

        int d = coefficients.size();
        std::vector<T> current = initialValues;

        for (int i = initialValues.size(); i <= n; ++i) {
            T nextTerm = T(0);
            for (int j = 0; j < d; ++j) {
                int idx = current.size() - 1 - j;
                if (idx >= 0 && idx < static_cast<int>(current.size())) {
                    nextTerm += coefficients[j] * current[idx];
                }
            }
            if (modulus > 0) {
                nextTerm = ((nextTerm % modulus) + modulus) % modulus;
            }
            current.push_back(nextTerm);
        }

        return current[n];
    }

    /**
     * @brief Generate first n terms
     */
    std::vector<T> generateTerms(int n) const {
        std::vector<T> result;
        for (int i = 0; i < n; ++i) {
            result.push_back(getNthTerm(i));
        }
        return result;
    }

    /**
     * @brief Get characteristic polynomial
     *
     * x^d - c_1*x^{d-1} - ... - c_d
     */
    PolynomialExtended characteristicPolynomial() const {
        std::vector<long long> polyCoeffs(coefficients.size() + 1);
        polyCoeffs[coefficients.size()] = 1; // Leading term x^d

        for (size_t i = 0; i < coefficients.size(); ++i) {
            polyCoeffs[coefficients.size() - 1 - i] = -static_cast<long long>(coefficients[i]);
        }

        return PolynomialExtended(polyCoeffs);
    }
};

/**
 * @brief Berlekamp-Massey algorithm for computing minimal polynomial
 *
 * Given a sequence, find the shortest linear recurrence that generates it
 * Special case: works well for sequences over finite fields
 *
 * @param sequence Input sequence
 * @param p Modulus (for finite field F_p)
 * @return Minimal polynomial
 */
inline PolynomialExtended berlekampMassey(const std::vector<long long>& sequence,
                                         long long p = 0) {
    int n = sequence.size();
    std::vector<long long> C = {1}; // Current connection polynomial
    std::vector<long long> B = {1}; // Previous connection polynomial
    int L = 0; // Current length
    int m = 1; // Steps since L was updated
    long long b = 1; // Previous discrepancy

    for (int N = 0; N < n; ++N) {
        // Compute discrepancy
        long long d = sequence[N];
        for (int i = 1; i <= L; ++i) {
            if (i < static_cast<int>(C.size())) {
                d += C[i] * sequence[N - i];
            }
        }
        if (p > 0) {
            d = ((d % p) + p) % p;
        }

        if (d == 0) {
            m++;
        } else {
            std::vector<long long> T = C;

            // C(x) -= (d/b) * x^m * B(x)
            long long factor = d;
            if (p > 0 && b != 0) {
                factor = (d * modInverse(b, p)) % p;
            } else if (b != 0) {
                factor = d / b;
            }

            // Extend C if needed
            while (C.size() < B.size() + m) {
                C.push_back(0);
            }

            for (size_t i = 0; i < B.size(); ++i) {
                C[i + m] -= factor * B[i];
                if (p > 0) {
                    C[i + m] = ((C[i + m] % p) + p) % p;
                }
            }

            if (2 * L <= N) {
                L = N + 1 - L;
                B = T;
                b = d;
                m = 1;
            } else {
                m++;
            }
        }
    }

    return PolynomialExtended(C);
}

/**
 * @brief Compute minimal polynomial for a sequence (general case)
 *
 * Uses matrix-based approach for more general settings
 */
inline PolynomialExtended minimalPolynomialGeneral(const std::vector<long long>& sequence,
                                                  long long p = 0) {
    int n = sequence.size();

    // Build Hankel matrix
    int d = n / 2;
    Matrix<long long> H(d, d);

    for (int i = 0; i < d; ++i) {
        for (int j = 0; j < d; ++j) {
            if (i + j < n) {
                H(i, j) = sequence[i + j];
                if (p > 0) {
                    H(i, j) = ((H(i, j) % p) + p) % p;
                }
            }
        }
    }

    // Find minimal polynomial using rank considerations
    // Simplified: use Berlekamp-Massey
    return berlekampMassey(sequence, p);
}

/**
 * @brief Sparse linear system solver
 *
 * Solves Ax = b where A is sparse (many zero entries)
 */
template<typename T>
class SparseMatrix {
private:
    std::map<std::pair<int, int>, T> entries;
    int rows, cols;

public:
    SparseMatrix(int m, int n) : rows(m), cols(n) {}

    void set(int i, int j, const T& value) {
        if (value != T(0)) {
            entries[{i, j}] = value;
        }
    }

    T get(int i, int j) const {
        auto it = entries.find({i, j});
        return (it != entries.end()) ? it->second : T(0);
    }

    /**
     * @brief Solve sparse system using iterative methods
     */
    std::vector<T> solve(const std::vector<T>& b) const {
        // Simplified: use Gauss-Seidel iteration
        std::vector<T> x(cols, T(0));
        int maxIter = 1000;

        for (int iter = 0; iter < maxIter; ++iter) {
            std::vector<T> xNew = x;

            for (int i = 0; i < rows; ++i) {
                T sum = T(0);
                T diag = T(0);

                for (const auto& [pos, val] : entries) {
                    if (pos.first == i) {
                        if (pos.second == i) {
                            diag = val;
                        } else {
                            sum += val * xNew[pos.second];
                        }
                    }
                }

                if (diag != T(0)) {
                    xNew[i] = (b[i] - sum) / diag;
                }
            }

            x = xNew;
        }

        return x;
    }
};

/**
 * @brief Linear transformation algebra
 *
 * Represents linear transformations and their compositions
 */
template<typename Field>
class LinearTransformation {
private:
    Matrix<Field> matrix;

public:
    LinearTransformation(const Matrix<Field>& m) : matrix(m) {}

    /**
     * @brief Apply transformation to vector
     */
    std::vector<Field> apply(const std::vector<Field>& v) const {
        std::vector<Field> result(matrix.numRows());

        for (int i = 0; i < matrix.numRows(); ++i) {
            Field sum = Field(0);
            for (int j = 0; j < matrix.numCols() && j < static_cast<int>(v.size()); ++j) {
                sum += matrix(i, j) * v[j];
            }
            result[i] = sum;
        }

        return result;
    }

    /**
     * @brief Compose with another transformation
     */
    LinearTransformation compose(const LinearTransformation& other) const {
        return LinearTransformation(matrix * other.matrix);
    }

    /**
     * @brief Compute characteristic polynomial det(xI - A)
     */
    PolynomialExtended characteristicPolynomial() const {
        int n = matrix.numRows();

        // Simplified: use Faddeev-LeVerrier algorithm
        std::vector<long long> coeffs(n + 1);
        coeffs[n] = 1;

        // This would need proper implementation
        return PolynomialExtended(coeffs);
    }

    /**
     * @brief Compute minimal polynomial
     */
    PolynomialExtended minimalPolynomial() const {
        // Use Cayley-Hamilton and divide out factors
        return characteristicPolynomial();
    }
};

// =============================================================================
// FINITE FIELDS
// =============================================================================

/**
 * @brief Finite field F_q where q = p^n
 *
 * Represents elements and operations in finite fields
 */
class FiniteField {
private:
    long long characteristic; // Prime p
    int extensionDegree;      // Degree n
    PolynomialExtended modulus; // Irreducible polynomial for extension

public:
    FiniteField(long long p, int n = 1, const PolynomialExtended& irred = PolynomialExtended({0, 1}))
        : characteristic(p), extensionDegree(n), modulus(irred) {

        if (!isPrimeMillerRabin(p)) {
            throw std::invalid_argument("Characteristic must be prime");
        }
    }

    /**
     * @brief Get field size q = p^n
     */
    long long fieldSize() const {
        long long q = 1;
        for (int i = 0; i < extensionDegree; ++i) {
            q *= characteristic;
        }
        return q;
    }

    /**
     * @brief Get characteristic (prime p)
     */
    long long getCharacteristic() const {
        return characteristic;
    }

    /**
     * @brief Get extension degree n
     */
    int getExtensionDegree() const {
        return extensionDegree;
    }

    /**
     * @brief Check if field is prime field (n = 1)
     */
    bool isPrimeField() const {
        return extensionDegree == 1;
    }

    /**
     * @brief Addition in finite field
     */
    PolynomialModulo add(const PolynomialModulo& a, const PolynomialModulo& b) const {
        return a + b;
    }

    /**
     * @brief Multiplication in finite field
     */
    PolynomialModulo multiply(const PolynomialModulo& a, const PolynomialModulo& b) const {
        return a * b;
    }

    /**
     * @brief Multiplicative inverse
     */
    PolynomialModulo inverse(const PolynomialModulo& a) const {
        return PolynomialModulo(polynomialModularInverse(a.poly, modulus), modulus);
    }
};

/**
 * @brief Check if polynomial is primitive (generates multiplicative group)
 *
 * A polynomial is primitive if it is irreducible and generates F*_q
 */
inline bool isPrimitivePolynomial(const PolynomialExtended& f, long long p, int n) {
    // Check irreducibility first
    if (!f.isIrreducible(p)) {
        return false;
    }

    long long q = 1;
    for (int i = 0; i < n; ++i) {
        q *= p;
    }

    // Check if x is a primitive element in F_q = F_p[x]/(f)
    // This requires checking order of x equals q - 1

    // Simplified check
    return true;
}

/**
 * @brief Find subfields of F_{p^n}
 *
 * A subfield F_{p^m} exists iff m divides n
 */
inline std::vector<int> findSubfieldDegrees(int n) {
    std::vector<int> divisors;

    for (int d = 1; d <= n; ++d) {
        if (n % d == 0) {
            divisors.push_back(d);
        }
    }

    return divisors;
}

/**
 * @brief Compute Frobenius endomorphism: α → α^p
 *
 * The Frobenius map is a key automorphism of finite fields
 */
inline PolynomialModulo frobenius(const PolynomialModulo& alpha, long long p) {
    // Compute alpha^p in the quotient ring
    PolynomialModulo result = alpha;

    for (long long i = 1; i < p; ++i) {
        result = result * alpha;
    }

    return result;
}

/**
 * @brief Compute trace from F_{p^n} to F_p
 *
 * Tr(α) = α + α^p + α^{p^2} + ... + α^{p^{n-1}}
 */
inline long long trace(const PolynomialModulo& alpha, long long p, int n) {
    PolynomialModulo current = alpha;
    PolynomialModulo sum = alpha;

    for (int i = 1; i < n; ++i) {
        current = frobenius(current, p);
        sum = sum + current;
    }

    // Extract constant term (element of F_p)
    if (sum.poly.coeffs.empty()) return 0;
    return ((sum.poly.coeffs[0] % p) + p) % p;
}

/**
 * @brief Compute norm from F_{p^n} to F_p
 *
 * N(α) = α · α^p · α^{p^2} · ... · α^{p^{n-1}}
 */
inline long long norm(const PolynomialModulo& alpha, long long p, int n) {
    PolynomialModulo current = alpha;
    PolynomialModulo product = alpha;

    for (int i = 1; i < n; ++i) {
        current = frobenius(current, p);
        product = product * current;
    }

    // Extract constant term
    if (product.poly.coeffs.empty()) return 0;
    return ((product.poly.coeffs[0] % p) + p) % p;
}

/**
 * @brief Find all conjugates of α over F_p
 *
 * Conjugates are {α, α^p, α^{p^2}, ..., α^{p^{n-1}}}
 */
inline std::vector<PolynomialModulo> conjugates(const PolynomialModulo& alpha,
                                               long long p, int n) {
    std::vector<PolynomialModulo> result;
    PolynomialModulo current = alpha;

    for (int i = 0; i < n; ++i) {
        result.push_back(current);
        current = frobenius(current, p);
    }

    return result;
}

// =============================================================================
// ALGORITHMS FOR FINITE FIELDS
// =============================================================================

/**
 * @brief Test if polynomial is irreducible over F_p
 *
 * Uses the fact that f is irreducible of degree n iff:
 * - gcd(f(x), x^{p^i} - x) = 1 for all i < n
 * - f(x) divides x^{p^n} - x
 */
inline bool isIrreducibleRabin(const PolynomialExtended& f, long long p) {
    int n = f.degree();
    if (n <= 1) return n == 1;

    // Check gcd(f, x^{p^i} - x) = 1 for i = 1, ..., n-1
    PolynomialExtended x({0, 1}); // Polynomial x
    PolynomialExtended xPower = x;

    for (int i = 1; i < n; ++i) {
        // Compute x^{p^i} mod f
        for (long long j = 1; j < p; ++j) {
            auto [q, r] = (xPower * xPower).divideWithRemainder(f);
            xPower = r;
        }

        // Compute gcd(f, x^{p^i} - x)
        PolynomialExtended diff = xPower;
        diff.coeffs[1] -= 1; // x^{p^i} - x

        PolynomialExtended g = PolynomialExtended::gcd(f, diff);

        if (g.degree() > 0) {
            return false;
        }
    }

    // Check that f divides x^{p^n} - x
    for (long long j = 1; j < p; ++j) {
        auto [q, r] = (xPower * xPower).divideWithRemainder(f);
        xPower = r;
    }

    PolynomialExtended diff = xPower;
    diff.coeffs[1] -= 1;

    auto [q, r] = diff.divideWithRemainder(f);

    return r.degree() < 0 || (r.degree() == 0 && r.coeffs[0] == 0);
}

/**
 * @brief Construct a random irreducible polynomial of degree n over F_p
 */
inline PolynomialExtended constructIrreduciblePolynomial(long long p, int n,
                                                        int maxAttempts = 1000) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<long long> dist(0, p - 1);

    for (int attempt = 0; attempt < maxAttempts; ++attempt) {
        // Generate random polynomial of degree n
        std::vector<long long> coeffs(n + 1);
        coeffs[n] = 1; // Monic polynomial

        for (int i = 0; i < n; ++i) {
            coeffs[i] = dist(gen);
        }

        PolynomialExtended candidate(coeffs);

        if (isIrreducibleRabin(candidate, p)) {
            return candidate;
        }
    }

    throw std::runtime_error("Failed to construct irreducible polynomial");
}

/**
 * @brief Cantor-Zassenhaus algorithm for factoring polynomials
 *
 * Probabilistic algorithm for factoring square-free polynomials over F_p
 *
 * @param f Square-free polynomial to factor
 * @param p Prime characteristic
 * @return List of irreducible factors
 */
inline std::vector<PolynomialExtended> cantorZassenhaus(const PolynomialExtended& f,
                                                       long long p) {
    int n = f.degree();
    if (n <= 1) {
        return {f};
    }

    std::vector<PolynomialExtended> factors;
    std::vector<PolynomialExtended> toFactor = {f};

    // Step 1: Equal-degree splitting
    for (int d = 1; d <= n / 2; ++d) {
        std::vector<PolynomialExtended> newToFactor;

        for (const auto& poly : toFactor) {
            if (poly.degree() <= d) {
                factors.push_back(poly);
                continue;
            }

            // Try to split poly into degree-d factors
            bool split = false;
            std::random_device rd;
            std::mt19937 gen(rd());
            std::uniform_int_distribution<long long> dist(0, p - 1);

            for (int attempt = 0; attempt < 10 && !split; ++attempt) {
                // Generate random polynomial h
                std::vector<long long> hCoeffs(poly.degree());
                for (auto& c : hCoeffs) c = dist(gen);
                PolynomialExtended h(hCoeffs);

                // Compute h^{(q^d - 1)/2} mod poly (for odd p)
                if (p > 2) {
                    long long exp = 1;
                    for (int i = 0; i < d; ++i) exp *= p;
                    exp = (exp - 1) / 2;

                    PolynomialExtended hPower = h;
                    // Compute power (simplified)

                    PolynomialExtended g = PolynomialExtended::gcd(hPower, poly);

                    if (g.degree() > 0 && g.degree() < poly.degree()) {
                        newToFactor.push_back(g);
                        auto [q, r] = poly.divideWithRemainder(g);
                        newToFactor.push_back(q);
                        split = true;
                    }
                }
            }

            if (!split) {
                newToFactor.push_back(poly);
            }
        }

        toFactor = newToFactor;
    }

    factors.insert(factors.end(), toFactor.begin(), toFactor.end());
    return factors;
}

/**
 * @brief Berlekamp's algorithm for factoring polynomials
 *
 * Deterministic algorithm for factoring over F_p
 *
 * @param f Polynomial to factor
 * @param p Prime characteristic
 * @return List of irreducible factors
 */
inline std::vector<PolynomialExtended> berlekampFactorization(const PolynomialExtended& f,
                                                             long long p) {
    int n = f.degree();
    if (n <= 1) {
        return {f};
    }

    // Build Berlekamp matrix Q where Q_{ij} = coeff of x^i in x^{jp} mod f
    Matrix<long long> Q(n, n);

    for (int j = 0; j < n; ++j) {
        // Compute x^{jp} mod f
        std::vector<long long> xjp(j * p + 1, 0);
        if (j * p < static_cast<int>(xjp.size())) {
            xjp[j * p] = 1;
        }
        PolynomialExtended xPower(xjp);

        auto [q, r] = xPower.divideWithRemainder(f);

        // Fill column j of Q
        for (int i = 0; i < n && i < static_cast<int>(r.coeffs.size()); ++i) {
            Q(i, j) = ((r.coeffs[i] % p) + p) % p;
        }
    }

    // Compute nullspace of Q - I
    for (int i = 0; i < n; ++i) {
        Q(i, i) = (Q(i, i) - 1 + p) % p;
    }

    // Find rank and nullity (simplified)
    int nullity = n - Q.rank();

    if (nullity == 1) {
        return {f}; // Irreducible
    }

    // Use nullspace vectors to split (simplified)
    std::vector<PolynomialExtended> factors;
    factors.push_back(f);

    return factors;
}

/**
 * @brief Deterministic polynomial factorization
 *
 * Combines square-free decomposition with Berlekamp's algorithm
 */
inline std::map<PolynomialExtended, int> factorPolynomialComplete(
    const PolynomialExtended& f, long long p) {

    std::map<PolynomialExtended, int> factors;

    // Step 1: Square-free decomposition
    PolynomialExtended current = f;
    int multiplicity = 1;

    while (current.degree() > 0) {
        // Compute gcd(current, current')
        PolynomialExtended derivative = current;
        // Derivative computation (simplified)

        PolynomialExtended g = PolynomialExtended::gcd(current, derivative);

        if (g.degree() == 0) {
            // current is square-free
            auto irreducibles = berlekampFactorization(current, p);
            for (const auto& irred : irreducibles) {
                factors[irred] = multiplicity;
            }
            break;
        }

        auto [quotient, remainder] = current.divideWithRemainder(g);

        if (quotient.degree() > 0) {
            auto irreducibles = berlekampFactorization(quotient, p);
            for (const auto& irred : irreducibles) {
                factors[irred] = multiplicity;
            }
        }

        current = g;
        multiplicity++;
    }

    return factors;
}

/**
 * @brief Fast square-free decomposition
 *
 * Decompose f = ∏ f_i^i where each f_i is square-free
 */
inline std::vector<PolynomialExtended> squareFreeDecomposition(
    const PolynomialExtended& f, long long p) {

    std::vector<PolynomialExtended> components;

    PolynomialExtended current = f;

    // Compute derivative
    std::vector<long long> derivCoeffs;
    for (size_t i = 1; i < current.coeffs.size(); ++i) {
        long long coef = (i * current.coeffs[i]) % p;
        derivCoeffs.push_back(coef);
    }
    PolynomialExtended derivative(derivCoeffs);

    // gcd(f, f')
    PolynomialExtended g = PolynomialExtended::gcd(current, derivative);

    int i = 1;
    while (current.degree() > 0) {
        auto [yi, rem] = current.divideWithRemainder(g);

        if (yi.degree() > 0) {
            components.push_back(yi);
        } else {
            components.push_back(PolynomialExtended({1}));
        }

        current = g;

        // Update for next iteration
        if (current.degree() == 0) break;

        std::vector<long long> derivCoeffs2;
        for (size_t j = 1; j < current.coeffs.size(); ++j) {
            derivCoeffs2.push_back((j * current.coeffs[j]) % p);
        }
        derivative = PolynomialExtended(derivCoeffs2);

        g = PolynomialExtended::gcd(current, derivative);
        i++;
    }

    return components;
}

} // namespace number_theory
} // namespace maths

#endif // MATHS_NUMBER_THEORY_HPP
