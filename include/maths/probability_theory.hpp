#ifndef MATHS_PROBABILITY_THEORY_HPP
#define MATHS_PROBABILITY_THEORY_HPP

#include <vector>
#include <set>
#include <map>
#include <functional>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <random>
#include <limits>

namespace maths {
namespace probability_theory {

using Outcome = int;
using Event = std::set<Outcome>;
using SigmaAlgebra = std::set<Event>;

// ============================================================================
// PROBABILITY SPACES
// ============================================================================

class ProbabilitySpace {
protected:
    Event sample_space_;
    SigmaAlgebra sigma_algebra_;
    std::map<Event, double> probability_measure_;

public:
    ProbabilitySpace(const Event& Omega) : sample_space_(Omega) {
        sigma_algebra_.insert(Event());
        sigma_algebra_.insert(Omega);
        probability_measure_[Event()] = 0.0;
        probability_measure_[Omega] = 1.0;
    }

    void addEvent(const Event& E, double prob) {
        if (prob < 0.0 || prob > 1.0) {
            throw std::invalid_argument("Probability must be in [0,1]");
        }
        sigma_algebra_.insert(E);
        probability_measure_[E] = prob;

        Event complement;
        std::set_difference(sample_space_.begin(), sample_space_.end(),
                          E.begin(), E.end(),
                          std::inserter(complement, complement.begin()));
        sigma_algebra_.insert(complement);
        probability_measure_[complement] = 1.0 - prob;
    }

    double P(const Event& E) const {
        auto it = probability_measure_.find(E);
        if (it != probability_measure_.end()) {
            return it->second;
        }
        return 0.0;
    }

    Event complement(const Event& E) const {
        Event comp;
        std::set_difference(sample_space_.begin(), sample_space_.end(),
                          E.begin(), E.end(),
                          std::inserter(comp, comp.begin()));
        return comp;
    }

    Event unionEvent(const Event& A, const Event& B) const {
        Event result;
        std::set_union(A.begin(), A.end(), B.begin(), B.end(),
                      std::inserter(result, result.begin()));
        return result;
    }

    Event intersection(const Event& A, const Event& B) const {
        Event result;
        std::set_intersection(A.begin(), A.end(), B.begin(), B.end(),
                            std::inserter(result, result.begin()));
        return result;
    }

    double conditionalProbability(const Event& A, const Event& B) const {
        double P_B = P(B);
        if (P_B < 1e-15) return 0.0;
        return P(intersection(A, B)) / P_B;
    }

    bool areIndependent(const Event& A, const Event& B, double tol = 1e-10) const {
        double P_A = P(A);
        double P_B = P(B);
        double P_AB = P(intersection(A, B));
        return std::abs(P_AB - P_A * P_B) < tol;
    }
};

class DiscreteProbabilitySpace : public ProbabilitySpace {
public:
    DiscreteProbabilitySpace(const std::vector<double>& probabilities)
        : ProbabilitySpace(Event()) {

        double sum = std::accumulate(probabilities.begin(), probabilities.end(), 0.0);
        if (std::abs(sum - 1.0) > 1e-10) {
            throw std::invalid_argument("Probabilities must sum to 1");
        }

        for (size_t i = 0; i < probabilities.size(); ++i) {
            sample_space_.insert(i);
            Event singleton = {static_cast<Outcome>(i)};
            sigma_algebra_.insert(singleton);
            probability_measure_[singleton] = probabilities[i];
        }
        probability_measure_[sample_space_] = 1.0;
    }
};

// ============================================================================
// RANDOM VARIABLES
// ============================================================================

class RandomVariable {
protected:
    ProbabilitySpace* space_;
    std::function<double(Outcome)> X_;

public:
    RandomVariable(ProbabilitySpace* space, std::function<double(Outcome)> X)
        : space_(space), X_(X) {}

    double operator()(Outcome omega) const {
        return X_(omega);
    }

    double expectation(const std::vector<Outcome>& outcomes,
                      const std::vector<double>& probabilities) const {
        double E = 0.0;
        for (size_t i = 0; i < outcomes.size(); ++i) {
            E += X_(outcomes[i]) * probabilities[i];
        }
        return E;
    }

    double variance(const std::vector<Outcome>& outcomes,
                   const std::vector<double>& probabilities) const {
        double E_X = expectation(outcomes, probabilities);
        double E_X2 = 0.0;
        for (size_t i = 0; i < outcomes.size(); ++i) {
            double val = X_(outcomes[i]);
            E_X2 += val * val * probabilities[i];
        }
        return E_X2 - E_X * E_X;
    }

    double standardDeviation(const std::vector<Outcome>& outcomes,
                            const std::vector<double>& probabilities) const {
        return std::sqrt(variance(outcomes, probabilities));
    }

    ProbabilitySpace* space() const {
        return space_;
    }

    std::vector<double> cumulativeDistribution(const std::vector<Outcome>& outcomes,
                                              const std::vector<double>& probabilities,
                                              const std::vector<double>& x_values) const {
        std::vector<double> F_X;
        for (double x : x_values) {
            double prob = 0.0;
            for (size_t i = 0; i < outcomes.size(); ++i) {
                if (X_(outcomes[i]) <= x) {
                    prob += probabilities[i];
                }
            }
            F_X.push_back(prob);
        }
        return F_X;
    }

    double moment(int k, const std::vector<Outcome>& outcomes,
                 const std::vector<double>& probabilities) const {
        double M_k = 0.0;
        for (size_t i = 0; i < outcomes.size(); ++i) {
            M_k += std::pow(X_(outcomes[i]), k) * probabilities[i];
        }
        return M_k;
    }

    double centralMoment(int k, const std::vector<Outcome>& outcomes,
                        const std::vector<double>& probabilities) const {
        double mu = expectation(outcomes, probabilities);
        double M_k = 0.0;
        for (size_t i = 0; i < outcomes.size(); ++i) {
            M_k += std::pow(X_(outcomes[i]) - mu, k) * probabilities[i];
        }
        return M_k;
    }
};

// ============================================================================
// CONVERGENCE OF RANDOM VARIABLES
// ============================================================================

class Convergence {
public:
    static bool convergesInProbability(
        const std::vector<RandomVariable>& sequence,
        const RandomVariable& limit,
        const std::vector<Outcome>& outcomes,
        const std::vector<double>& probabilities,
        double epsilon = 1e-6) {

        int n = sequence.size();
        if (n == 0) return false;

        const auto& X_n = sequence.back();
        double prob_diff = 0.0;
        for (size_t i = 0; i < outcomes.size(); ++i) {
            if (std::abs(X_n(outcomes[i]) - limit(outcomes[i])) > epsilon) {
                prob_diff += probabilities[i];
            }
        }
        return prob_diff < epsilon;
    }

    static bool convergesAlmostSurely(
        const std::vector<RandomVariable>& sequence,
        const RandomVariable& limit,
        const std::vector<Outcome>& outcomes,
        double epsilon = 1e-6) {

        for (const auto& omega : outcomes) {
            bool converges_at_omega = true;
            for (size_t n = sequence.size() - 10; n < sequence.size(); ++n) {
                if (std::abs(sequence[n](omega) - limit(omega)) > epsilon) {
                    converges_at_omega = false;
                    break;
                }
            }
            if (!converges_at_omega) return false;
        }
        return true;
    }

    static bool convergesInDistribution(
        const std::vector<RandomVariable>& sequence,
        const RandomVariable& limit,
        const std::vector<Outcome>& outcomes,
        const std::vector<double>& probabilities,
        const std::vector<double>& test_points,
        double tol = 1e-6) {

        auto F_X = limit.cumulativeDistribution(outcomes, probabilities, test_points);
        auto F_Xn = sequence.back().cumulativeDistribution(outcomes, probabilities, test_points);

        for (size_t i = 0; i < test_points.size(); ++i) {
            if (std::abs(F_Xn[i] - F_X[i]) > tol) {
                return false;
            }
        }
        return true;
    }
};

// ============================================================================
// LAW OF LARGE NUMBERS
// ============================================================================

class LawOfLargeNumbers {
public:
    static bool verifyWeakLLN(
        const std::vector<RandomVariable>& iid_sequence,
        const std::vector<Outcome>& outcomes,
        const std::vector<double>& probabilities,
        double epsilon = 1e-6) {

        double mu = iid_sequence[0].expectation(outcomes, probabilities);

        auto sample_mean = [&](int n) {
            return [&, n](Outcome omega) {
                double sum = 0.0;
                for (int i = 0; i < n && i < static_cast<int>(iid_sequence.size()); ++i) {
                    sum += iid_sequence[i](omega);
                }
                return sum / n;
            };
        };

        int n = iid_sequence.size();
        RandomVariable X_bar(iid_sequence[0].space(), sample_mean(n));

        double prob_deviation = 0.0;
        for (size_t i = 0; i < outcomes.size(); ++i) {
            if (std::abs(X_bar(outcomes[i]) - mu) > epsilon) {
                prob_deviation += probabilities[i];
            }
        }

        double variance = iid_sequence[0].variance(outcomes, probabilities);
        double chebyshev_bound = variance / (n * epsilon * epsilon);

        return prob_deviation <= chebyshev_bound + 1e-10;
    }

    static std::vector<double> sampleMeans(
        const std::vector<std::vector<double>>& samples) {

        std::vector<double> means;
        int cumulative_n = 0;
        double cumulative_sum = 0.0;

        for (const auto& sample : samples) {
            for (double x : sample) {
                cumulative_sum += x;
                ++cumulative_n;
                means.push_back(cumulative_sum / cumulative_n);
            }
        }
        return means;
    }
};

// ============================================================================
// CENTRAL LIMIT THEOREM
// ============================================================================

class CentralLimitTheorem {
public:
    static std::vector<double> standardizedSums(
        const std::vector<std::vector<double>>& iid_samples,
        double mu, double sigma) {

        std::vector<double> Z_n;
        for (size_t n = 1; n < iid_samples.size(); ++n) {
            double sum = 0.0;
            for (size_t i = 0; i <= n; ++i) {
                sum += iid_samples[i][0];
            }
            double Z = (sum - (n + 1) * mu) / (sigma * std::sqrt(n + 1));
            Z_n.push_back(Z);
        }
        return Z_n;
    }

    static bool verifyCLT(const std::vector<double>& standardized_sums,
                         int n_bins = 20, double tol = 0.1) {

        if (standardized_sums.empty()) return false;

        double z_min = *std::min_element(standardized_sums.begin(), standardized_sums.end());
        double z_max = *std::max_element(standardized_sums.begin(), standardized_sums.end());
        double bin_width = (z_max - z_min) / n_bins;

        std::vector<int> histogram(n_bins, 0);
        for (double z : standardized_sums) {
            int bin = std::min(static_cast<int>((z - z_min) / bin_width), n_bins - 1);
            histogram[bin]++;
        }

        auto normal_pdf = [](double z) {
            return std::exp(-0.5 * z * z) / std::sqrt(2.0 * M_PI);
        };

        double chi_squared = 0.0;
        for (int i = 0; i < n_bins; ++i) {
            double z_center = z_min + (i + 0.5) * bin_width;
            double expected = normal_pdf(z_center) * bin_width * standardized_sums.size();
            if (expected > 0) {
                chi_squared += std::pow(histogram[i] - expected, 2) / expected;
            }
        }

        double critical_value = 30.0;
        return chi_squared < critical_value;
    }
};

// ============================================================================
// MARTINGALES
// ============================================================================

class Martingale {
private:
    std::vector<RandomVariable> sequence_;
    std::vector<std::vector<Event>> filtration_;

public:
    Martingale(const std::vector<RandomVariable>& seq,
              const std::vector<std::vector<Event>>& filt)
        : sequence_(seq), filtration_(filt) {}

    bool isMartingale(const std::vector<Outcome>& outcomes,
                     const std::vector<double>& probabilities,
                     double tol = 1e-10) const {

        for (size_t n = 0; n < sequence_.size() - 1; ++n) {
            for (const auto& F_n_event : filtration_[n]) {
                double E_Xn = 0.0;
                double E_Xn1_given_Fn = 0.0;
                double P_Fn = 0.0;

                for (size_t i = 0; i < outcomes.size(); ++i) {
                    if (F_n_event.count(outcomes[i]) > 0) {
                        E_Xn += sequence_[n](outcomes[i]) * probabilities[i];
                        E_Xn1_given_Fn += sequence_[n+1](outcomes[i]) * probabilities[i];
                        P_Fn += probabilities[i];
                    }
                }

                if (P_Fn > 1e-15) {
                    E_Xn /= P_Fn;
                    E_Xn1_given_Fn /= P_Fn;

                    if (std::abs(E_Xn1_given_Fn - E_Xn) > tol) {
                        return false;
                    }
                }
            }
        }
        return true;
    }

    bool isSubmartingale(const std::vector<Outcome>& outcomes,
                        const std::vector<double>& probabilities,
                        double tol = 1e-10) const {

        for (size_t n = 0; n < sequence_.size() - 1; ++n) {
            double E_Xn = sequence_[n].expectation(outcomes, probabilities);
            double E_Xn1 = sequence_[n+1].expectation(outcomes, probabilities);

            if (E_Xn1 < E_Xn - tol) {
                return false;
            }
        }
        return true;
    }

    bool isSupermartingale(const std::vector<Outcome>& outcomes,
                          const std::vector<double>& probabilities,
                          double tol = 1e-10) const {

        for (size_t n = 0; n < sequence_.size() - 1; ++n) {
            double E_Xn = sequence_[n].expectation(outcomes, probabilities);
            double E_Xn1 = sequence_[n+1].expectation(outcomes, probabilities);

            if (E_Xn1 > E_Xn + tol) {
                return false;
            }
        }
        return true;
    }

    std::vector<double> doobDecomposition(
        const std::vector<Outcome>& outcomes,
        const std::vector<double>& probabilities) const {

        std::vector<double> martingale_part;
        std::vector<double> predictable_part;

        martingale_part.push_back(sequence_[0].expectation(outcomes, probabilities));
        predictable_part.push_back(0.0);

        for (size_t n = 1; n < sequence_.size(); ++n) {
            double E_Xn = sequence_[n].expectation(outcomes, probabilities);
            double E_Xn_minus_1 = sequence_[n-1].expectation(outcomes, probabilities);

            predictable_part.push_back(E_Xn - E_Xn_minus_1 + predictable_part.back());
            martingale_part.push_back(E_Xn - predictable_part.back());
        }

        return martingale_part;
    }
};

// ============================================================================
// STOPPING TIMES
// ============================================================================

class StoppingTime {
private:
    std::function<int(const std::vector<Outcome>&)> tau_;

public:
    StoppingTime(std::function<int(const std::vector<Outcome>&)> tau) : tau_(tau) {}

    int operator()(const std::vector<Outcome>& path) const {
        return tau_(path);
    }

    bool isStoppingTime(const std::vector<std::vector<Event>>& filtration,
                       const std::vector<Outcome>& sample_outcomes) const {

        std::vector<Outcome> path;
        for (size_t n = 0; n < filtration.size(); ++n) {
            path.push_back(sample_outcomes[n % sample_outcomes.size()]);
            int tau_value = tau_(path);

            if (tau_value == static_cast<int>(n)) {
                bool in_filtration = false;
                for (const auto& F_n_event : filtration[n]) {
                    if (F_n_event.count(sample_outcomes[n % sample_outcomes.size()]) > 0) {
                        in_filtration = true;
                        break;
                    }
                }
                if (!in_filtration) return false;
            }
        }
        return true;
    }

    double optionalStoppingExpectation(
        const Martingale& M,
        const std::vector<Outcome>& outcomes,
        const std::vector<double>& probabilities) const {

        double E_X_tau = 0.0;
        for (size_t i = 0; i < outcomes.size(); ++i) {
            std::vector<Outcome> path = {outcomes[i]};
            int tau_i = tau_(path);

            if (tau_i >= 0 && tau_i < 1) {
                E_X_tau += outcomes[i] * probabilities[i];
            }
        }
        return E_X_tau;
    }
};

// ============================================================================
// CHARACTERISTIC FUNCTIONS
// ============================================================================

class CharacteristicFunction {
public:
    static std::complex<double> compute(
        const RandomVariable& X,
        double t,
        const std::vector<Outcome>& outcomes,
        const std::vector<double>& probabilities) {

        std::complex<double> phi_t(0, 0);
        for (size_t i = 0; i < outcomes.size(); ++i) {
            double x = X(outcomes[i]);
            std::complex<double> exp_itx(std::cos(t * x), std::sin(t * x));
            phi_t += exp_itx * probabilities[i];
        }
        return phi_t;
    }

    static std::vector<std::complex<double>> characteristic(
        const RandomVariable& X,
        const std::vector<double>& t_values,
        const std::vector<Outcome>& outcomes,
        const std::vector<double>& probabilities) {

        std::vector<std::complex<double>> phi;
        for (double t : t_values) {
            phi.push_back(compute(X, t, outcomes, probabilities));
        }
        return phi;
    }

    static bool levyContinuityTheorem(
        const std::vector<std::complex<double>>& phi_n,
        const std::vector<std::complex<double>>& phi,
        double tol = 1e-6) {

        if (phi_n.size() != phi.size()) return false;

        for (size_t i = 0; i < phi.size(); ++i) {
            if (std::abs(phi_n[i] - phi[i]) > tol) {
                return false;
            }
        }
        return true;
    }
};

// ============================================================================
// BROWNIAN MOTION (Discrete Approximation)
// ============================================================================

class BrownianMotion {
private:
    std::mt19937 rng_;
    std::normal_distribution<double> normal_;

public:
    BrownianMotion(unsigned seed = std::random_device{}())
        : rng_(seed), normal_(0.0, 1.0) {}

    std::vector<double> path(double T, int n_steps) {
        std::vector<double> B(n_steps + 1, 0.0);
        double dt = T / n_steps;
        double sqrt_dt = std::sqrt(dt);

        for (int i = 1; i <= n_steps; ++i) {
            B[i] = B[i-1] + sqrt_dt * normal_(rng_);
        }
        return B;
    }

    std::vector<std::vector<double>> multiplePaths(double T, int n_steps, int n_paths) {
        std::vector<std::vector<double>> paths;
        for (int i = 0; i < n_paths; ++i) {
            paths.push_back(path(T, n_steps));
        }
        return paths;
    }

    bool verifyMartingaleProperty(const std::vector<double>& path, double tol = 0.5) {
        double expected = 0.0;
        for (size_t i = 0; i < path.size(); ++i) {
            expected += path[i];
        }
        expected /= path.size();
        return std::abs(expected) < tol;
    }

    bool verifyQuadraticVariation(double T, int n_steps, int n_realizations,
                                  double tol = 0.2) {
        double total_QV = 0.0;

        for (int r = 0; r < n_realizations; ++r) {
            auto B = path(T, n_steps);
            double QV = 0.0;
            for (size_t i = 1; i < B.size(); ++i) {
                QV += std::pow(B[i] - B[i-1], 2);
            }
            total_QV += QV;
        }

        double avg_QV = total_QV / n_realizations;
        return std::abs(avg_QV - T) < tol;
    }
};

} // namespace probability_theory
} // namespace maths

#endif // MATHS_PROBABILITY_THEORY_HPP
