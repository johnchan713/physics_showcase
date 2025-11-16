#ifndef PHYSICS_CONDENSED_MATTER_HPP
#define PHYSICS_CONDENSED_MATTER_HPP

#include <vector>
#include <complex>
#include <cmath>
#include <functional>
#include <algorithm>
#include <numeric>

namespace physics {
namespace condensed_matter {

using Complex = std::complex<double>;

constexpr double k_B = 1.380649e-23;
constexpr double hbar = 1.054571817e-34;
constexpr double e = 1.602176634e-19;
constexpr double m_e = 9.1093837015e-31;

// ============================================================================
// BAND STRUCTURE
// ============================================================================

class BandStructure {
public:
    virtual double energy(const std::vector<double>& k) const = 0;
    virtual ~BandStructure() = default;

    double groupVelocity(const std::vector<double>& k, int direction,
                        double dk = 1e-6) const {
        std::vector<double> k_plus = k, k_minus = k;
        k_plus[direction] += dk;
        k_minus[direction] -= dk;
        return (energy(k_plus) - energy(k_minus)) / (2.0 * hbar * dk);
    }

    double effectiveMass(const std::vector<double>& k, int direction,
                        double dk = 1e-6) const {
        std::vector<double> k_plus = k, k_minus = k;
        k_plus[direction] += dk;
        k_minus[direction] -= dk;
        double d2E = (energy(k_plus) - 2.0 * energy(k) + energy(k_minus)) / (dk * dk);
        return hbar * hbar / d2E;
    }

    std::vector<double> fermiSurface(double E_F, double k_max, int n_points = 100) const {
        std::vector<double> fermi_surface;
        for (int i = 0; i < n_points; ++i) {
            double k = -k_max + 2.0 * k_max * i / n_points;
            std::vector<double> k_vec = {k, 0, 0};
            if (std::abs(energy(k_vec) - E_F) < 0.1 * E_F) {
                fermi_surface.push_back(k);
            }
        }
        return fermi_surface;
    }
};

class FreElectronBand : public BandStructure {
public:
    double energy(const std::vector<double>& k) const override {
        double k2 = 0.0;
        for (double ki : k) {
            k2 += ki * ki;
        }
        return hbar * hbar * k2 / (2.0 * m_e);
    }
};

class TightBindingBand : public BandStructure {
private:
    double t_;
    double a_;

public:
    TightBindingBand(double hopping, double lattice_const)
        : t_(hopping), a_(lattice_const) {}

    double energy(const std::vector<double>& k) const override {
        double E = 0.0;
        for (size_t i = 0; i < k.size() && i < 3; ++i) {
            E += -2.0 * t_ * std::cos(k[i] * a_);
        }
        return E;
    }
};

class KronigPenneyModel : public BandStructure {
private:
    double V0_;
    double a_;

public:
    KronigPenneyModel(double barrier_height, double lattice_const)
        : V0_(barrier_height), a_(lattice_const) {}

    double energy(const std::vector<double>& k) const override {
        double ka = k[0] * a_;
        double alpha = std::sqrt(2.0 * m_e * V0_) / hbar;
        double cos_ka = std::cos(ka);

        return hbar * hbar * k[0] * k[0] / (2.0 * m_e)
             + V0_ * (1.0 - cos_ka) / (2.0 + alpha * alpha * a_ * a_);
    }
};

// ============================================================================
// DENSITY OF STATES
// ============================================================================

class DensityOfStates {
public:
    static double freElectron3D(double E) {
        if (E < 0) return 0.0;
        return (1.0 / (2.0 * M_PI * M_PI)) * std::pow(2.0 * m_e / (hbar * hbar), 1.5)
             * std::sqrt(E);
    }

    static double freElectron2D(double E) {
        if (E < 0) return 0.0;
        return m_e / (M_PI * hbar * hbar);
    }

    static double freElectron1D(double E) {
        if (E <= 0) return 0.0;
        return std::sqrt(m_e / (2.0 * E)) / (M_PI * hbar);
    }

    static double integrate(std::function<double(double)> g, double E_min,
                          double E_max, int n_points = 1000) {
        double sum = 0.0;
        double dE = (E_max - E_min) / n_points;
        for (int i = 0; i < n_points; ++i) {
            double E = E_min + (i + 0.5) * dE;
            sum += g(E) * dE;
        }
        return sum;
    }

    static double fermiLevel(std::function<double(double)> g, double n_electrons,
                            double E_max, int n_points = 1000) {
        double E_F = 0.0;
        double dE = E_max / n_points;

        for (int i = 0; i < n_points; ++i) {
            double E = i * dE;
            double n = integrate(g, 0.0, E, n_points);
            if (n >= n_electrons) {
                E_F = E;
                break;
            }
        }
        return E_F;
    }
};

// ============================================================================
// BCS SUPERCONDUCTIVITY
// ============================================================================

class BCSTheory {
private:
    double T_c_;
    double Delta_0_;
    double E_F_;

public:
    BCSTheory(double critical_temp, double fermi_energy)
        : T_c_(critical_temp), E_F_(fermi_energy) {
        Delta_0_ = 1.764 * k_B * T_c_;
    }

    double energyGap(double T) const {
        if (T >= T_c_) return 0.0;

        double t = T / T_c_;
        return Delta_0_ * std::sqrt(1.0 - std::pow(t, 4));
    }

    double quasiparticleEnergy(double epsilon, double T) const {
        double Delta = energyGap(T);
        return std::sqrt(epsilon * epsilon + Delta * Delta);
    }

    double coherenceLength(double T) const {
        if (T >= T_c_) return 0.0;

        double Delta = energyGap(T);
        double v_F = std::sqrt(2.0 * E_F_ / m_e);
        return hbar * v_F / (M_PI * Delta);
    }

    double penetrationDepth(double n_s, double T) const {
        if (T >= T_c_) return std::numeric_limits<double>::infinity();

        return std::sqrt(m_e / (mu_0 * n_s * e * e));
    }

    double cooperPairBindingEnergy() const {
        return 2.0 * Delta_0_;
    }

    double criticalField(double T) const {
        if (T >= T_c_) return 0.0;

        double Delta = energyGap(T);
        return Delta * std::sqrt(2.0) / (1.602e-19);
    }

    double specificHeat(double T) const {
        if (T >= T_c_) {
            return normalSpecificHeat(T);
        }

        double Delta = energyGap(T);
        return 1.43 * normalSpecificHeat(T_c_) * std::exp(-Delta / (k_B * T));
    }

    double normalSpecificHeat(double T) const {
        return (M_PI * M_PI / 3.0) * k_B * k_B * T;
    }

    double condensationEnergy() const {
        return 0.5 * Delta_0_ * Delta_0_ / k_B;
    }

private:
    static constexpr double mu_0 = 1.25663706212e-6;
};

// ============================================================================
// LONDON EQUATIONS
// ============================================================================

class LondonTheory {
private:
    double lambda_L_;
    double n_s_;

public:
    LondonTheory(double penetration_depth, double superfluid_density)
        : lambda_L_(penetration_depth), n_s_(superfluid_density) {}

    std::vector<double> superfluidVelocity(const std::vector<double>& A) const {
        std::vector<double> v_s(A.size());
        for (size_t i = 0; i < A.size(); ++i) {
            v_s[i] = -(e / (m_e)) * A[i];
        }
        return v_s;
    }

    std::vector<double> supercurrent(const std::vector<double>& A) const {
        auto v_s = superfluidVelocity(A);
        std::vector<double> j_s(v_s.size());
        for (size_t i = 0; i < v_s.size(); ++i) {
            j_s[i] = -n_s_ * e * v_s[i];
        }
        return j_s;
    }

    double magneticFieldPenetration(double x) const {
        return std::exp(-x / lambda_L_);
    }

    double meissnerEffect(double B_external, double x) const {
        return B_external * std::exp(-x / lambda_L_);
    }

    double fluxQuantum() const {
        return M_PI * hbar / e;
    }

    int fluxoidQuantization(double total_flux) const {
        double Phi_0 = fluxQuantum();
        return static_cast<int>(std::round(total_flux / Phi_0));
    }
};

// ============================================================================
// JOSEPHSON EFFECT
// ============================================================================

class JosephsonEffect {
private:
    double I_c_;
    double V_c_;

public:
    JosephsonEffect(double critical_current, double gap_voltage)
        : I_c_(critical_current), V_c_(gap_voltage) {}

    double dcCurrent(double phase_diff) const {
        return I_c_ * std::sin(phase_diff);
    }

    double acJosephson(double V, double t) const {
        double omega_J = 2.0 * e * V / hbar;
        return I_c_ * std::sin(omega_J * t);
    }

    double josephsonFrequency(double V) const {
        return 2.0 * e * V / hbar;
    }

    double inverseFrohlichEffect(double omega) const {
        return hbar * omega / (2.0 * e);
    }

    double shapiroSteps(int n, double omega_rf) const {
        return n * hbar * omega_rf / (2.0 * e);
    }

    double criticalCurrent() const {
        return I_c_;
    }

    double junctionEnergy(double phase) const {
        double Phi_0 = M_PI * hbar / e;
        return -(Phi_0 * I_c_ / (2.0 * M_PI)) * std::cos(phase);
    }
};

// ============================================================================
// FERMI LIQUID THEORY
// ============================================================================

class FermiLiquidTheory {
private:
    double m_star_;
    double m_;

public:
    FermiLiquidTheory(double effective_mass, double bare_mass)
        : m_star_(effective_mass), m_(bare_mass) {}

    double massEnhancement() const {
        return m_star_ / m_;
    }

    double quasiparticleEnergy(double k, double k_F) const {
        double v_F = hbar * k_F / m_star_;
        return v_F * (k - k_F);
    }

    double quasiparticleLifetime(double E, double T) const {
        double E_F = hbar * hbar * 1e10 / (2.0 * m_star_);
        if (std::abs(E) < k_B * T) {
            return hbar / (k_B * T);
        }
        return hbar * E_F / (E * E);
    }

    double specificHeat(double T, double g_0) const {
        return (m_star_ / m_) * g_0 * M_PI * M_PI * k_B * k_B * T / 3.0;
    }

    double compressibility(double n, double k_F) const {
        double E_F = hbar * hbar * k_F * k_F / (2.0 * m_star_);
        return 3.0 * n / (2.0 * E_F);
    }

    double spinSusceptibility(double chi_0) const {
        return (m_star_ / m_) * chi_0;
    }
};

// ============================================================================
// PHONONS
// ============================================================================

class PhononDispersion {
private:
    double omega_D_;
    double v_s_;

public:
    PhononDispersion(double debye_freq, double sound_velocity)
        : omega_D_(debye_freq), v_s_(v_s_) {}

    double acousticBranch(double k) const {
        return v_s_ * k;
    }

    double opticalBranch(double k, double omega_0) const {
        return omega_0;
    }

    double debyeFrequency() const {
        return omega_D_;
    }

    double debyeTemperature() const {
        return hbar * omega_D_ / k_B;
    }

    double boseEinsteinDistribution(double omega, double T) const {
        double x = hbar * omega / (k_B * T);
        if (x > 100) return 0.0;
        return 1.0 / (std::exp(x) - 1.0);
    }

    double phononHeatCapacity(double T, double T_D) const {
        if (T < T_D / 50.0) {
            return 12.0 * M_PI * M_PI * M_PI * M_PI * k_B * std::pow(T / T_D, 3) / 5.0;
        }
        return 3.0 * k_B;
    }

    double thermalConductivity(double T, double mean_free_path) const {
        double c_v = phononHeatCapacity(T, debyeTemperature());
        return (1.0 / 3.0) * c_v * v_s_ * mean_free_path;
    }
};

// ============================================================================
// QUANTUM HALL EFFECT
// ============================================================================

class QuantumHallEffect {
private:
    double B_;
    int filling_factor_;

public:
    QuantumHallEffect(double magnetic_field, int nu)
        : B_(magnetic_field), filling_factor_(nu) {}

    double hallConductance() const {
        return filling_factor_ * e * e / (2.0 * M_PI * hbar);
    }

    double hallResistance() const {
        double sigma_H = hallConductance();
        return 1.0 / sigma_H;
    }

    double vonKlitzingConstant() const {
        return 2.0 * M_PI * hbar / (e * e);
    }

    double landauLevelEnergy(int n) const {
        double omega_c = e * B_ / m_e;
        return hbar * omega_c * (n + 0.5);
    }

    double magneticLength() const {
        return std::sqrt(hbar / (e * B_));
    }

    double cyclotronFrequency() const {
        return e * B_ / m_e;
    }

    int fillingFactor() const {
        return filling_factor_;
    }

    double laughlinWavefunction(const std::vector<Complex>& z, int m) const {
        Complex psi = 1.0;
        for (size_t i = 0; i < z.size(); ++i) {
            for (size_t j = i + 1; j < z.size(); ++j) {
                psi *= std::pow(z[i] - z[j], m);
            }
        }
        for (const auto& zi : z) {
            psi *= std::exp(-0.25 * std::norm(zi));
        }
        return std::abs(psi);
    }
};

// ============================================================================
// LANDAU THEORY OF PHASE TRANSITIONS
// ============================================================================

class LandauTheory {
private:
    double T_c_;

public:
    LandauTheory(double critical_temp) : T_c_(critical_temp) {}

    double freeEnergy(double m, double T, double h, double a, double b, double c) const {
        double t = (T - T_c_) / T_c_;
        return a * t * m * m + b * m * m * m * m + c * m * m * m * m * m * m - h * m;
    }

    double orderParameter(double T, double b) const {
        if (T >= T_c_) return 0.0;

        double t = (T - T_c_) / T_c_;
        return std::sqrt(-t / (2.0 * b));
    }

    double susceptibility(double T, double a) const {
        double t = (T - T_c_) / T_c_;
        if (T > T_c_) {
            return 1.0 / (2.0 * a * t);
        }
        return 1.0 / (4.0 * a * std::abs(t));
    }

    double correlationLength(double T, double xi_0) const {
        double t = std::abs((T - T_c_) / T_c_);
        return xi_0 / std::sqrt(t);
    }

    double criticalExponentBeta() const {
        return 0.5;
    }

    double criticalExponentGamma() const {
        return 1.0;
    }

    double criticalExponentNu() const {
        return 0.5;
    }
};

// ============================================================================
// HUBBARD MODEL
// ============================================================================

class HubbardModel {
private:
    double t_;
    double U_;
    int L_;

public:
    HubbardModel(double hopping, double interaction, int lattice_size)
        : t_(hopping), U_(interaction), L_(lattice_size) {}

    double kineticEnergy(const std::vector<int>& occupation) const {
        double E_kin = 0.0;
        for (size_t i = 0; i < occupation.size() - 1; ++i) {
            E_kin += -t_ * occupation[i] * occupation[i + 1];
        }
        return E_kin;
    }

    double interactionEnergy(const std::vector<int>& n_up,
                            const std::vector<int>& n_down) const {
        double E_int = 0.0;
        for (size_t i = 0; i < n_up.size(); ++i) {
            E_int += U_ * n_up[i] * n_down[i];
        }
        return E_int;
    }

    double totalEnergy(const std::vector<int>& n_up,
                      const std::vector<int>& n_down) const {
        std::vector<int> n_total(n_up.size());
        for (size_t i = 0; i < n_up.size(); ++i) {
            n_total[i] = n_up[i] + n_down[i];
        }
        return kineticEnergy(n_total) + interactionEnergy(n_up, n_down);
    }

    bool isMottInsulator(double filling) const {
        return (std::abs(filling - 1.0) < 0.1) && (U_ / t_ > 4.0);
    }

    double bandGap() const {
        if (U_ / t_ > 4.0) {
            return U_ - 8.0 * t_;
        }
        return 0.0;
    }
};

} // namespace condensed_matter
} // namespace physics

#endif // PHYSICS_CONDENSED_MATTER_HPP
