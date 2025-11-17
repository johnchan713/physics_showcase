#ifndef MATHS_ACTUARIAL_LIFE_TABLES_HPP
#define MATHS_ACTUARIAL_LIFE_TABLES_HPP

#include <cmath>
#include <vector>
#include <string>
#include <stdexcept>
#include <algorithm>
#include <functional>

/**
 * @file life_tables.hpp
 * @brief Actuarial life tables and survival analysis
 *
 * Implements:
 * - Mortality and survival functions
 * - Life expectancy calculations
 * - Actuarial present values
 * - Life annuities and life insurance
 * - Premium calculation principles
 * - Loss distributions
 */

namespace maths::actuarial {

/**
 * @class MortalityFunctions
 * @brief Fundamental mortality and survival functions
 */
class MortalityFunctions {
public:
    /**
     * @brief Survival function: S(x) = P(T > x)
     *
     * Probability of surviving to age x
     */
    static double survivalFunction(double x, std::function<double(double)> hazard) {
        // S(x) = exp(-∫₀ˣ μ(t)dt) where μ is hazard/force of mortality
        // For simple implementation, approximate with discrete sum
        double integral = 0.0;
        double dt = 0.01;
        int steps = static_cast<int>(x / dt);

        for (int i = 0; i < steps; ++i) {
            integral += hazard(i * dt) * dt;
        }

        return std::exp(-integral);
    }

    /**
     * @brief Probability density function of death
     */
    static double deathPDF(double x, std::function<double(double)> hazard) {
        // f(x) = S(x) · μ(x)
        return survivalFunction(x, hazard) * hazard(x);
    }

    /**
     * @brief n-year survival probability: ₙp_x = P(T_x > n)
     *
     * Probability that (x) survives n more years
     */
    static double nYearSurvival(double x, double n, std::function<double(double)> hazard) {
        // ₙp_x = S(x+n) / S(x)
        if (n < 0.0) {
            throw std::invalid_argument("n must be non-negative");
        }

        double S_x = survivalFunction(x, hazard);
        double S_xn = survivalFunction(x + n, hazard);

        if (S_x < 1e-10) return 0.0;
        return S_xn / S_x;
    }

    /**
     * @brief n-year death probability: ₙq_x = P(T_x ≤ n)
     */
    static double nYearDeath(double x, double n, std::function<double(double)> hazard) {
        return 1.0 - nYearSurvival(x, n, hazard);
    }

    /**
     * @brief Deferred death probability: ₘ|ₙq_x
     *
     * Probability that (x) survives m years but dies within next n years
     */
    static double deferredDeath(double x, double m, double n,
                                std::function<double(double)> hazard) {
        // ₘ|ₙq_x = ₘp_x · ₙq_{x+m}
        double prob_survive_m = nYearSurvival(x, m, hazard);
        double prob_die_next_n = nYearDeath(x + m, n, hazard);
        return prob_survive_m * prob_die_next_n;
    }

    /**
     * @brief Fundamental relationships
     */
    static std::string relationships() {
        return "Fundamental Mortality Relationships:\n"
               "\n"
               "Survival function: S(x) = P(T > x)\n"
               "PDF of death: f(x) = -S'(x) = S(x)·μ(x)\n"
               "Hazard/force of mortality: μ(x) = f(x)/S(x) = -S'(x)/S(x)\n"
               "\n"
               "From life table at age x:\n"
               "- lₓ = number alive at age x\n"
               "- dₓ = lₓ - lₓ₊₁ = deaths between x and x+1\n"
               "- qₓ = dₓ/lₓ = probability of death within 1 year\n"
               "- pₓ = 1 - qₓ = probability of survival 1 year\n"
               "\n"
               "General relationships:\n"
               "- ₙp_x = lₓ₊ₙ/lₓ\n"
               "- ₙq_x = 1 - ₙp_x\n"
               "- ₘ₊ₙp_x = ₘp_x · ₙp_{x+m}\n"
               "- ₘ|ₙq_x = ₘp_x · ₙq_{x+m}\n"
               "\n"
               "Standard mortality laws:\n"
               "- Gompertz: μ(x) = Be^(cx) (exponential increase)\n"
               "- Makeham: μ(x) = A + Be^(cx) (adds constant)\n"
               "- Weibull: μ(x) = kx^(n-1) (power law)";
    }
};

/**
 * @class LifeExpectancy
 * @brief Complete and curtate life expectancy
 */
class LifeExpectancy {
public:
    /**
     * @brief Complete expectation of life: e̊ₓ = E[T_x]
     *
     * Expected future lifetime at age x
     */
    static double completeExpectation(double x, std::function<double(double)> hazard,
                                      double max_age = 120.0) {
        // e̊ₓ = ∫₀^∞ ₜp_x dt
        double expectation = 0.0;
        double dt = 0.1;
        int steps = static_cast<int>((max_age - x) / dt);

        for (int i = 0; i < steps; ++i) {
            double t = i * dt;
            expectation += MortalityFunctions::nYearSurvival(x, t, hazard) * dt;
        }

        return expectation;
    }

    /**
     * @brief Curtate expectation of life: eₓ = E[K_x]
     *
     * Expected number of complete future years
     * K_x = floor(T_x) = curtate future lifetime
     */
    static double curtateExpectation(double x, std::function<double(double)> hazard,
                                     double max_age = 120.0) {
        // eₓ = Σ(k=1 to ∞) ₖp_x
        double expectation = 0.0;
        int max_years = static_cast<int>(max_age - x);

        for (int k = 1; k <= max_years; ++k) {
            expectation += MortalityFunctions::nYearSurvival(x, k, hazard);
        }

        return expectation;
    }

    /**
     * @brief Temporary life expectancy: e̊ₓ:n̅|
     *
     * Expected lifetime in next n years (capped at n)
     */
    static double temporaryExpectation(double x, double n,
                                       std::function<double(double)> hazard) {
        double expectation = 0.0;
        double dt = 0.1;
        int steps = static_cast<int>(n / dt);

        for (int i = 0; i < steps; ++i) {
            double t = i * dt;
            expectation += MortalityFunctions::nYearSurvival(x, t, hazard) * dt;
        }

        return expectation;
    }
};

/**
 * @class ActuarialPresentValue
 * @brief Present values of life-contingent cashflows
 */
class ActuarialPresentValue {
public:
    /**
     * @brief Whole life insurance: Āₓ
     *
     * Present value of $1 paid at moment of death
     * Āₓ = E[v^T_x] where v = e^(-δ)
     */
    static double wholeLifeInsurance(double x, double delta,
                                     std::function<double(double)> hazard,
                                     double max_age = 120.0) {
        // Āₓ = ∫₀^∞ e^(-δt) · ₜp_x · μ(x+t) dt
        double apv = 0.0;
        double dt = 0.1;
        int steps = static_cast<int>((max_age - x) / dt);

        for (int i = 0; i < steps; ++i) {
            double t = i * dt;
            double discount = std::exp(-delta * t);
            double survival = MortalityFunctions::nYearSurvival(x, t, hazard);
            double force = hazard(x + t);

            apv += discount * survival * force * dt;
        }

        return apv;
    }

    /**
     * @brief Term life insurance: Āₓ:n̅|
     *
     * $1 paid if death occurs within n years
     */
    static double termLifeInsurance(double x, double n, double delta,
                                    std::function<double(double)> hazard) {
        double apv = 0.0;
        double dt = 0.1;
        int steps = static_cast<int>(n / dt);

        for (int i = 0; i < steps; ++i) {
            double t = i * dt;
            double discount = std::exp(-delta * t);
            double survival = MortalityFunctions::nYearSurvival(x, t, hazard);
            double force = hazard(x + t);

            apv += discount * survival * force * dt;
        }

        return apv;
    }

    /**
     * @brief Endowment insurance: Āₓ:n̅|
     *
     * $1 paid at death if within n years, or at time n if alive
     */
    static double endowmentInsurance(double x, double n, double delta,
                                     std::function<double(double)> hazard) {
        // Term insurance + pure endowment
        double term = termLifeInsurance(x, n, delta, hazard);
        double pure_endow = std::exp(-delta * n) *
                           MortalityFunctions::nYearSurvival(x, n, hazard);
        return term + pure_endow;
    }

    /**
     * @brief Whole life annuity: äₓ
     *
     * Present value of $1/year paid while alive (continuous)
     */
    static double wholeLifeAnnuity(double x, double delta,
                                   std::function<double(double)> hazard,
                                   double max_age = 120.0) {
        // äₓ = ∫₀^∞ e^(-δt) · ₜp_x dt
        double apv = 0.0;
        double dt = 0.1;
        int steps = static_cast<int>((max_age - x) / dt);

        for (int i = 0; i < steps; ++i) {
            double t = i * dt;
            apv += std::exp(-delta * t) *
                   MortalityFunctions::nYearSurvival(x, t, hazard) * dt;
        }

        return apv;
    }

    /**
     * @brief Temporary annuity: äₓ:n̅|
     *
     * $1/year for n years or until death, whichever is first
     */
    static double temporaryAnnuity(double x, double n, double delta,
                                   std::function<double(double)> hazard) {
        double apv = 0.0;
        double dt = 0.1;
        int steps = static_cast<int>(n / dt);

        for (int i = 0; i < steps; ++i) {
            double t = i * dt;
            apv += std::exp(-delta * t) *
                   MortalityFunctions::nYearSurvival(x, t, hazard) * dt;
        }

        return apv;
    }

    /**
     * @brief Deferred annuity: ₘ|äₓ
     *
     * Annuity starting m years from now
     */
    static double deferredAnnuity(double x, double m, double delta,
                                  std::function<double(double)> hazard,
                                  double max_age = 120.0) {
        // ₘ|äₓ = ₘp_x · v^m · äₓ₊ₘ
        double prob_survive = MortalityFunctions::nYearSurvival(x, m, hazard);
        double discount = std::exp(-delta * m);
        double annuity_value = wholeLifeAnnuity(x + m, delta, hazard, max_age);

        return prob_survive * discount * annuity_value;
    }

    /**
     * @brief Actuarial notation and relationships
     */
    static std::string notation() {
        return "Actuarial Present Value Notation:\n"
               "\n"
               "Insurance (lump sum at death):\n"
               "- Āₓ = whole life insurance\n"
               "- Āₓ:n̅| = n-year term insurance\n"
               "- Āₓ:n̅| = n-year endowment insurance\n"
               "- Aₓ = discrete version (paid end of year of death)\n"
               "\n"
               "Annuities (regular payments while alive):\n"
               "- äₓ = continuous whole life annuity-due\n"
               "- äₓ:n̅| = n-year temporary annuity\n"
               "- aₓ = discrete annuity (payments at year-end)\n"
               "\n"
               "Key relationships:\n"
               "- äₓ = (1 - Āₓ)/δ (insurance-annuity identity)\n"
               "- Āₓ:n̅| = Āₓ:n̅| + v^n · ₙp_x (endowment = term + pure)\n"
               "- Variance: ²Āₓ - (Āₓ)² where ²Āₓ uses v^(2T_x)\n"
               "\n"
               "Premium principles:\n"
               "- Equivalence: Premium PV = Benefit PV\n"
               "- Net premium: P = Āₓ / äₓ\n"
               "- Gross premium: includes expenses, profit loading";
    }
};

/**
 * @class PremiumCalculation
 * @brief Insurance premium calculation principles
 */
class PremiumCalculation {
public:
    /**
     * @brief Net single premium
     *
     * Lump sum paid at inception
     */
    static double netSinglePremium(double benefit,
                                   double insurance_apv) {
        return benefit * insurance_apv;
    }

    /**
     * @brief Net level premium
     *
     * Equal payments while alive (using equivalence principle)
     */
    static double netLevelPremium(double benefit,
                                  double insurance_apv,
                                  double annuity_apv) {
        if (annuity_apv < 1e-10) {
            throw std::runtime_error("Annuity APV too small");
        }
        return benefit * insurance_apv / annuity_apv;
    }

    /**
     * @brief Gross premium with loadings
     */
    static double grossPremium(double net_premium,
                              double expense_loading = 0.0,
                              double profit_loading = 0.0) {
        // Gross = Net + Expenses + Profit
        return net_premium * (1.0 + expense_loading + profit_loading);
    }

    /**
     * @brief Premium principles
     */
    static std::string principles() {
        return "Premium Calculation Principles:\n"
               "\n"
               "1. Equivalence Principle:\n"
               "   E[PV(Premiums)] = E[PV(Benefits)]\n"
               "   Most common, results in fair premium\n"
               "\n"
               "2. Expected Value Principle:\n"
               "   Premium = (1 + θ)·E[Loss]\n"
               "   where θ is safety loading\n"
               "\n"
               "3. Variance Principle:\n"
               "   Premium = E[Loss] + α·Var[Loss]\n"
               "   Accounts for risk via variance\n"
               "\n"
               "4. Standard Deviation Principle:\n"
               "   Premium = E[Loss] + β·SD[Loss]\n"
               "\n"
               "5. Percentile Principle:\n"
               "   Premium = VaR_p[Loss] (pth percentile)\n"
               "   E.g., p = 0.95 for 95% confidence\n"
               "\n"
               "Components of gross premium:\n"
               "- Net premium (equivalence)\n"
               "- Expense loading (acquisition, admin)\n"
               "- Profit loading (company profit margin)\n"
               "- Safety margin (adverse deviation)";
    }
};

/**
 * @class Reserves
 * @brief Policy reserves and accounting
 */
class Reserves {
public:
    /**
     * @brief Prospective reserve at time t
     *
     * Reserve = E[PV(Future Benefits) - PV(Future Premiums)]
     */
    static std::string prospectiveMethod() {
        return "Prospective Reserve Method:\n"
               "\n"
               "Reserve at time t:\n"
               "  ₜV = E[PV(Future Benefits) - PV(Future Premiums) | Aₜ]\n"
               "\n"
               "For whole life policy issued at age x:\n"
               "  ₜV = Āₓ₊ₜ - P·äₓ₊ₜ\n"
               "\n"
               "where P = Āₓ/äₓ is net level premium.\n"
               "\n"
               "Properties:\n"
               "- At issue: ₀V = 0 (equivalence principle)\n"
               "- Increases over time as mortality risk rises\n"
               "- At maturity: reserve = face amount\n"
               "\n"
               "Reserves are:\n"
               "- Liability on insurer's balance sheet\n"
               "- Assets backing future obligations\n"
               "- Required by regulation";
    }

    /**
     * @brief Retrospective reserve
     */
    static std::string retrospectiveMethod() {
        return "Retrospective Reserve Method:\n"
               "\n"
               "Reserve = Accumulated Value of:\n"
               "  (Past Premiums) - (Past Claims)\n"
               "\n"
               "For survivor:\n"
               "  ₜV = [P·s̈ₓ:t̅| - (1 - ₜpₓ)] / ₜpₓ\n"
               "\n"
               "where s̈ₓ:t̅| = accumulated annuity value.\n"
               "\n"
               "Equivalence: Prospective = Retrospective\n"
               "(fundamental theorem of reserves)\n"
               "\n"
               "Thiele's differential equation:\n"
               "  d/dt[ₜV] = δ·ₜV + P - μₓ₊ₜ·(1 - ₜV)\n"
               "\n"
               "Interpretation:\n"
               "- δ·ₜV: interest earned\n"
               "- P: premium income\n"
               "- μₓ₊ₜ·(1 - ₜV): expected death claim";
    }
};

} // namespace maths::actuarial

#endif // MATHS_ACTUARIAL_LIFE_TABLES_HPP
