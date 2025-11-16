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

} // namespace number_theory
} // namespace maths

#endif // MATHS_NUMBER_THEORY_HPP
