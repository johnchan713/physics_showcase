#ifndef PHYSICS_STATISTICAL_MODELS_HPP
#define PHYSICS_STATISTICAL_MODELS_HPP

#include <cmath>
#include <vector>
#include <stdexcept>
#include <map>
#include <set>
#include <algorithm>
#include <functional>
#include <queue>
#include <numeric>

/**
 * @file statistical_models.hpp
 * @brief Comprehensive implementation of statistical and graphical models
 *
 * This module implements:
 * - Ising model (1D, 2D, 3D)
 * - Markov Random Fields (MRF) on graphs
 * - Finite models (finite state machines, automata)
 * - Tree models (graphical models, belief propagation)
 * - Potts model
 * - Spin glass models
 * - Graph partitioning and clustering
 *
 * All calculations use standard statistical physics conventions.
 */

namespace physics {
namespace statistical_models {

// ============================================================================
// CONSTANTS
// ============================================================================

namespace constants {
    constexpr double BOLTZMANN_CONSTANT = 1.380649e-23;  // J/K
    constexpr double PI = 3.141592653589793;
}

// ============================================================================
// ISING MODEL - 1D
// ============================================================================

/**
 * @class IsingModel1D
 * @brief One-dimensional Ising model implementation
 *
 * Hamiltonian: H = -J·Σ s_i·s_{i+1} - h·Σ s_i
 * where s_i ∈ {-1, +1}, J is coupling constant, h is external field
 */
class IsingModel1D {
private:
    std::vector<int> spins;  // Spin configuration: +1 or -1
    double coupling;          // J: coupling constant
    double externalField;     // h: external magnetic field
    double temperature;       // T: temperature in Kelvin
    bool periodicBoundary;    // Periodic boundary conditions

public:
    /**
     * @brief Constructor for 1D Ising model
     *
     * @param size Number of spins
     * @param J Coupling constant (positive = ferromagnetic)
     * @param h External field
     * @param T Temperature (K)
     * @param periodic Use periodic boundary conditions
     */
    IsingModel1D(int size, double J, double h, double T, bool periodic = true)
        : spins(size, 1), coupling(J), externalField(h),
          temperature(T), periodicBoundary(periodic) {
        if (size <= 0) {
            throw std::invalid_argument("Size must be positive");
        }
        if (temperature <= 0.0) {
            throw std::invalid_argument("Temperature must be positive");
        }
    }

    /**
     * @brief Set spin at position i
     *
     * @param i Position
     * @param spin Spin value (+1 or -1)
     */
    void setSpin(int i, int spin) {
        if (i < 0 || i >= static_cast<int>(spins.size())) {
            throw std::out_of_range("Index out of range");
        }
        if (spin != 1 && spin != -1) {
            throw std::invalid_argument("Spin must be +1 or -1");
        }
        spins[i] = spin;
    }

    /**
     * @brief Get spin at position i
     */
    int getSpin(int i) const {
        if (i < 0 || i >= static_cast<int>(spins.size())) {
            throw std::out_of_range("Index out of range");
        }
        return spins[i];
    }

    /**
     * @brief Initialize spins randomly
     *
     * @param upProbability Probability of spin up (default 0.5)
     */
    void randomize(double upProbability = 0.5) {
        for (auto& spin : spins) {
            spin = (static_cast<double>(std::rand()) / RAND_MAX < upProbability) ? 1 : -1;
        }
    }

    /**
     * @brief Calculate total energy of configuration
     *
     * E = -J·Σ s_i·s_{i+1} - h·Σ s_i
     *
     * @return Total energy
     */
    double energy() const {
        double E = 0.0;
        int N = spins.size();

        // Nearest-neighbor interaction
        for (int i = 0; i < N - 1; ++i) {
            E -= coupling * spins[i] * spins[i + 1];
        }

        // Periodic boundary condition
        if (periodicBoundary && N > 1) {
            E -= coupling * spins[N - 1] * spins[0];
        }

        // External field contribution
        for (int spin : spins) {
            E -= externalField * spin;
        }

        return E;
    }

    /**
     * @brief Calculate local energy at site i
     *
     * E_i = -s_i·(J·(s_{i-1} + s_{i+1}) + h)
     *
     * @param i Site index
     * @return Local energy
     */
    double localEnergy(int i) const {
        if (i < 0 || i >= static_cast<int>(spins.size())) {
            throw std::out_of_range("Index out of range");
        }

        int N = spins.size();
        double neighborSum = 0.0;

        // Left neighbor
        if (i > 0) {
            neighborSum += spins[i - 1];
        } else if (periodicBoundary) {
            neighborSum += spins[N - 1];
        }

        // Right neighbor
        if (i < N - 1) {
            neighborSum += spins[i + 1];
        } else if (periodicBoundary) {
            neighborSum += spins[0];
        }

        return -spins[i] * (coupling * neighborSum + externalField);
    }

    /**
     * @brief Calculate energy change if spin i is flipped
     *
     * ΔE = 2·s_i·(J·(s_{i-1} + s_{i+1}) + h)
     *
     * @param i Site index
     * @return Energy change
     */
    double energyChangeFlip(int i) const {
        return -2.0 * localEnergy(i);
    }

    /**
     * @brief Calculate total magnetization
     *
     * M = Σ s_i
     *
     * @return Total magnetization
     */
    double magnetization() const {
        return std::accumulate(spins.begin(), spins.end(), 0.0);
    }

    /**
     * @brief Calculate magnetization per spin
     *
     * m = M/N = (1/N)·Σ s_i
     *
     * @return Magnetization per spin
     */
    double magnetizationPerSpin() const {
        return magnetization() / spins.size();
    }

    /**
     * @brief Perform single Monte Carlo step (Metropolis algorithm)
     *
     * @param i Site to attempt flip
     * @return true if flip was accepted
     */
    bool metropolisStep(int i) {
        double deltaE = energyChangeFlip(i);

        // Accept if energy decreases
        if (deltaE <= 0.0) {
            spins[i] = -spins[i];
            return true;
        }

        // Accept with Boltzmann probability
        double probability = std::exp(-deltaE / (constants::BOLTZMANN_CONSTANT * temperature));
        if (static_cast<double>(std::rand()) / RAND_MAX < probability) {
            spins[i] = -spins[i];
            return true;
        }

        return false;
    }

    /**
     * @brief Perform one Monte Carlo sweep (attempt to flip all spins)
     *
     * @return Number of accepted flips
     */
    int monteCarlSweep() {
        int accepted = 0;
        for (size_t i = 0; i < spins.size(); ++i) {
            if (metropolisStep(i)) {
                accepted++;
            }
        }
        return accepted;
    }

    /**
     * @brief Calculate partition function (exact for small systems)
     *
     * Z = Σ exp(-βE)
     *
     * Warning: Exponential complexity O(2^N)
     *
     * @return Partition function
     */
    double partitionFunction() const {
        int N = spins.size();
        if (N > 20) {
            throw std::runtime_error("System too large for exact partition function");
        }

        double beta = 1.0 / (constants::BOLTZMANN_CONSTANT * temperature);
        double Z = 0.0;

        // Iterate over all 2^N configurations
        for (long long config = 0; config < (1LL << N); ++config) {
            // Create temporary configuration
            std::vector<int> tempSpins(N);
            for (int i = 0; i < N; ++i) {
                tempSpins[i] = (config & (1LL << i)) ? 1 : -1;
            }

            // Calculate energy for this configuration
            double E = 0.0;
            for (int i = 0; i < N - 1; ++i) {
                E -= coupling * tempSpins[i] * tempSpins[i + 1];
            }
            if (periodicBoundary && N > 1) {
                E -= coupling * tempSpins[N - 1] * tempSpins[0];
            }
            for (int spin : tempSpins) {
                E -= externalField * spin;
            }

            Z += std::exp(-beta * E);
        }

        return Z;
    }

    /**
     * @brief Calculate free energy
     *
     * F = -kT·ln(Z)
     *
     * @return Free energy
     */
    double freeEnergy() const {
        double Z = partitionFunction();
        return -constants::BOLTZMANN_CONSTANT * temperature * std::log(Z);
    }

    /**
     * @brief Get system size
     */
    int size() const {
        return spins.size();
    }

    /**
     * @brief Get temperature
     */
    double getTemperature() const {
        return temperature;
    }

    /**
     * @brief Set temperature
     */
    void setTemperature(double T) {
        if (T <= 0.0) {
            throw std::invalid_argument("Temperature must be positive");
        }
        temperature = T;
    }
};

// ============================================================================
// ISING MODEL - 2D
// ============================================================================

/**
 * @class IsingModel2D
 * @brief Two-dimensional Ising model on square lattice
 *
 * Hamiltonian: H = -J·Σ_{<i,j>} s_i·s_j - h·Σ_i s_i
 * where <i,j> denotes nearest neighbors on square lattice
 */
class IsingModel2D {
private:
    std::vector<std::vector<int>> spins;  // 2D spin configuration
    int rows, cols;
    double coupling;
    double externalField;
    double temperature;
    bool periodicBoundary;

public:
    /**
     * @brief Constructor for 2D Ising model
     *
     * @param nRows Number of rows
     * @param nCols Number of columns
     * @param J Coupling constant
     * @param h External field
     * @param T Temperature (K)
     * @param periodic Use periodic boundary conditions
     */
    IsingModel2D(int nRows, int nCols, double J, double h, double T, bool periodic = true)
        : rows(nRows), cols(nCols), coupling(J), externalField(h),
          temperature(T), periodicBoundary(periodic) {
        if (nRows <= 0 || nCols <= 0) {
            throw std::invalid_argument("Dimensions must be positive");
        }
        if (temperature <= 0.0) {
            throw std::invalid_argument("Temperature must be positive");
        }

        spins.resize(rows, std::vector<int>(cols, 1));
    }

    /**
     * @brief Set spin at position (i, j)
     */
    void setSpin(int i, int j, int spin) {
        if (i < 0 || i >= rows || j < 0 || j >= cols) {
            throw std::out_of_range("Index out of range");
        }
        if (spin != 1 && spin != -1) {
            throw std::invalid_argument("Spin must be +1 or -1");
        }
        spins[i][j] = spin;
    }

    /**
     * @brief Get spin at position (i, j)
     */
    int getSpin(int i, int j) const {
        if (i < 0 || i >= rows || j < 0 || j >= cols) {
            throw std::out_of_range("Index out of range");
        }
        return spins[i][j];
    }

    /**
     * @brief Initialize spins randomly
     */
    void randomize(double upProbability = 0.5) {
        for (auto& row : spins) {
            for (auto& spin : row) {
                spin = (static_cast<double>(std::rand()) / RAND_MAX < upProbability) ? 1 : -1;
            }
        }
    }

    /**
     * @brief Calculate total energy
     *
     * E = -J·Σ_{<i,j>} s_i·s_j - h·Σ_i s_i
     *
     * @return Total energy
     */
    double energy() const {
        double E = 0.0;

        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < cols; ++j) {
                // Right neighbor
                int jRight = (j + 1) % cols;
                if (periodicBoundary || j + 1 < cols) {
                    E -= coupling * spins[i][j] * spins[i][jRight];
                }

                // Bottom neighbor
                int iBottom = (i + 1) % rows;
                if (periodicBoundary || i + 1 < rows) {
                    E -= coupling * spins[i][j] * spins[iBottom][j];
                }

                // External field
                E -= externalField * spins[i][j];
            }
        }

        return E;
    }

    /**
     * @brief Calculate energy change if spin at (i,j) is flipped
     *
     * @param i Row index
     * @param j Column index
     * @return Energy change
     */
    double energyChangeFlip(int i, int j) const {
        if (i < 0 || i >= rows || j < 0 || j >= cols) {
            throw std::out_of_range("Index out of range");
        }

        int neighborSum = 0;

        // Up neighbor
        int iUp = (i - 1 + rows) % rows;
        if (periodicBoundary || i > 0) {
            neighborSum += spins[iUp][j];
        }

        // Down neighbor
        int iDown = (i + 1) % rows;
        if (periodicBoundary || i < rows - 1) {
            neighborSum += spins[iDown][j];
        }

        // Left neighbor
        int jLeft = (j - 1 + cols) % cols;
        if (periodicBoundary || j > 0) {
            neighborSum += spins[i][jLeft];
        }

        // Right neighbor
        int jRight = (j + 1) % cols;
        if (periodicBoundary || j < cols - 1) {
            neighborSum += spins[i][jRight];
        }

        return 2.0 * spins[i][j] * (coupling * neighborSum + externalField);
    }

    /**
     * @brief Calculate total magnetization
     */
    double magnetization() const {
        double M = 0.0;
        for (const auto& row : spins) {
            M += std::accumulate(row.begin(), row.end(), 0.0);
        }
        return M;
    }

    /**
     * @brief Calculate magnetization per spin
     */
    double magnetizationPerSpin() const {
        return magnetization() / (rows * cols);
    }

    /**
     * @brief Perform single Metropolis step at site (i,j)
     */
    bool metropolisStep(int i, int j) {
        double deltaE = energyChangeFlip(i, j);

        if (deltaE <= 0.0) {
            spins[i][j] = -spins[i][j];
            return true;
        }

        double probability = std::exp(-deltaE / (constants::BOLTZMANN_CONSTANT * temperature));
        if (static_cast<double>(std::rand()) / RAND_MAX < probability) {
            spins[i][j] = -spins[i][j];
            return true;
        }

        return false;
    }

    /**
     * @brief Perform one Monte Carlo sweep
     */
    int monteCarloSweep() {
        int accepted = 0;
        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < cols; ++j) {
                if (metropolisStep(i, j)) {
                    accepted++;
                }
            }
        }
        return accepted;
    }

    /**
     * @brief Calculate critical temperature (Onsager solution)
     *
     * T_c = 2J / (k_B·ln(1 + √2)) ≈ 2.269·J/k_B
     *
     * Valid for 2D square lattice with zero external field
     *
     * @return Critical temperature (K)
     */
    static double criticalTemperature(double J) {
        return 2.0 * J / (constants::BOLTZMANN_CONSTANT * std::log(1.0 + std::sqrt(2.0)));
    }

    /**
     * @brief Get dimensions
     */
    std::pair<int, int> dimensions() const {
        return {rows, cols};
    }

    /**
     * @brief Get temperature
     */
    double getTemperature() const {
        return temperature;
    }

    /**
     * @brief Set temperature
     */
    void setTemperature(double T) {
        if (T <= 0.0) {
            throw std::invalid_argument("Temperature must be positive");
        }
        temperature = T;
    }
};

// ============================================================================
// MARKOV RANDOM FIELD ON GRAPHS
// ============================================================================

/**
 * @class MarkovRandomField
 * @brief Markov Random Field on arbitrary graph structure
 *
 * Energy: E(x) = Σ_i φ_i(x_i) + Σ_{i,j} ψ_{ij}(x_i, x_j)
 * Probability: P(x) = (1/Z)·exp(-E(x)/T)
 */
class MarkovRandomField {
private:
    int numNodes;
    int numStates;  // Number of states per node
    std::vector<int> nodeStates;  // Current state of each node
    std::vector<std::vector<int>> adjList;  // Adjacency list

    // Node potentials: φ_i(x_i)
    std::vector<std::vector<double>> nodePotentials;

    // Edge potentials: ψ_{ij}(x_i, x_j)
    std::map<std::pair<int, int>, std::vector<std::vector<double>>> edgePotentials;

    double temperature;

public:
    /**
     * @brief Constructor for MRF
     *
     * @param nodes Number of nodes
     * @param states Number of states per node
     * @param T Temperature parameter
     */
    MarkovRandomField(int nodes, int states, double T = 1.0)
        : numNodes(nodes), numStates(states), temperature(T) {
        if (nodes <= 0 || states <= 0) {
            throw std::invalid_argument("Nodes and states must be positive");
        }
        if (temperature <= 0.0) {
            throw std::invalid_argument("Temperature must be positive");
        }

        nodeStates.resize(nodes, 0);
        adjList.resize(nodes);
        nodePotentials.resize(nodes, std::vector<double>(states, 0.0));
    }

    /**
     * @brief Add edge to graph
     *
     * @param i First node
     * @param j Second node
     */
    void addEdge(int i, int j) {
        if (i < 0 || i >= numNodes || j < 0 || j >= numNodes) {
            throw std::out_of_range("Node index out of range");
        }
        if (i == j) {
            throw std::invalid_argument("Self-loops not allowed");
        }

        adjList[i].push_back(j);
        adjList[j].push_back(i);

        // Initialize edge potential to zero
        auto key = std::make_pair(std::min(i, j), std::max(i, j));
        if (edgePotentials.find(key) == edgePotentials.end()) {
            edgePotentials[key] = std::vector<std::vector<double>>(
                numStates, std::vector<double>(numStates, 0.0)
            );
        }
    }

    /**
     * @brief Set node potential
     *
     * @param node Node index
     * @param state State index
     * @param potential Potential value
     */
    void setNodePotential(int node, int state, double potential) {
        if (node < 0 || node >= numNodes) {
            throw std::out_of_range("Node index out of range");
        }
        if (state < 0 || state >= numStates) {
            throw std::out_of_range("State index out of range");
        }
        nodePotentials[node][state] = potential;
    }

    /**
     * @brief Set edge potential
     *
     * @param i First node
     * @param j Second node
     * @param stateI State of node i
     * @param stateJ State of node j
     * @param potential Potential value
     */
    void setEdgePotential(int i, int j, int stateI, int stateJ, double potential) {
        if (i < 0 || i >= numNodes || j < 0 || j >= numNodes) {
            throw std::out_of_range("Node index out of range");
        }
        if (stateI < 0 || stateI >= numStates || stateJ < 0 || stateJ >= numStates) {
            throw std::out_of_range("State index out of range");
        }

        auto key = std::make_pair(std::min(i, j), std::max(i, j));
        if (edgePotentials.find(key) == edgePotentials.end()) {
            throw std::invalid_argument("Edge does not exist");
        }

        if (i < j) {
            edgePotentials[key][stateI][stateJ] = potential;
        } else {
            edgePotentials[key][stateJ][stateI] = potential;
        }
    }

    /**
     * @brief Set state of node
     */
    void setNodeState(int node, int state) {
        if (node < 0 || node >= numNodes) {
            throw std::out_of_range("Node index out of range");
        }
        if (state < 0 || state >= numStates) {
            throw std::out_of_range("State index out of range");
        }
        nodeStates[node] = state;
    }

    /**
     * @brief Get state of node
     */
    int getNodeState(int node) const {
        if (node < 0 || node >= numNodes) {
            throw std::out_of_range("Node index out of range");
        }
        return nodeStates[node];
    }

    /**
     * @brief Calculate energy of current configuration
     *
     * E = Σ_i φ_i(x_i) + Σ_{i,j} ψ_{ij}(x_i, x_j)
     *
     * @return Energy
     */
    double energy() const {
        double E = 0.0;

        // Node potentials
        for (int i = 0; i < numNodes; ++i) {
            E += nodePotentials[i][nodeStates[i]];
        }

        // Edge potentials
        for (const auto& [edge, potential] : edgePotentials) {
            int i = edge.first;
            int j = edge.second;
            E += potential[nodeStates[i]][nodeStates[j]];
        }

        return E;
    }

    /**
     * @brief Calculate energy change if node changes state
     *
     * @param node Node to change
     * @param newState New state
     * @return Energy change
     */
    double energyChangeState(int node, int newState) const {
        if (node < 0 || node >= numNodes) {
            throw std::out_of_range("Node index out of range");
        }
        if (newState < 0 || newState >= numStates) {
            throw std::out_of_range("State index out of range");
        }

        int oldState = nodeStates[node];
        if (oldState == newState) {
            return 0.0;
        }

        double deltaE = 0.0;

        // Change in node potential
        deltaE += nodePotentials[node][newState] - nodePotentials[node][oldState];

        // Change in edge potentials
        for (int neighbor : adjList[node]) {
            auto key = std::make_pair(std::min(node, neighbor), std::max(node, neighbor));
            if (edgePotentials.find(key) != edgePotentials.end()) {
                int neighborState = nodeStates[neighbor];

                if (node < neighbor) {
                    deltaE += edgePotentials.at(key)[newState][neighborState] -
                              edgePotentials.at(key)[oldState][neighborState];
                } else {
                    deltaE += edgePotentials.at(key)[neighborState][newState] -
                              edgePotentials.at(key)[neighborState][oldState];
                }
            }
        }

        return deltaE;
    }

    /**
     * @brief Perform Gibbs sampling step at node
     *
     * Sample new state from conditional distribution:
     * P(x_i | x_{-i}) ∝ exp(-E_i(x_i)/T)
     *
     * @param node Node to sample
     */
    void gibbsSamplingStep(int node) {
        std::vector<double> probabilities(numStates);
        double maxLogProb = -std::numeric_limits<double>::infinity();

        // Calculate unnormalized probabilities for each state
        for (int state = 0; state < numStates; ++state) {
            double E = nodePotentials[node][state];

            // Add edge contributions
            for (int neighbor : adjList[node]) {
                auto key = std::make_pair(std::min(node, neighbor), std::max(node, neighbor));
                if (edgePotentials.find(key) != edgePotentials.end()) {
                    int neighborState = nodeStates[neighbor];

                    if (node < neighbor) {
                        E += edgePotentials.at(key)[state][neighborState];
                    } else {
                        E += edgePotentials.at(key)[neighborState][state];
                    }
                }
            }

            double logProb = -E / temperature;
            if (logProb > maxLogProb) {
                maxLogProb = logProb;
            }
            probabilities[state] = logProb;
        }

        // Normalize probabilities using log-sum-exp trick
        double sumProb = 0.0;
        for (int state = 0; state < numStates; ++state) {
            probabilities[state] = std::exp(probabilities[state] - maxLogProb);
            sumProb += probabilities[state];
        }

        for (double& prob : probabilities) {
            prob /= sumProb;
        }

        // Sample from distribution
        double r = static_cast<double>(std::rand()) / RAND_MAX;
        double cumulative = 0.0;
        for (int state = 0; state < numStates; ++state) {
            cumulative += probabilities[state];
            if (r <= cumulative) {
                nodeStates[node] = state;
                break;
            }
        }
    }

    /**
     * @brief Perform one Gibbs sampling sweep (sample all nodes)
     */
    void gibbsSamplingSweep() {
        for (int i = 0; i < numNodes; ++i) {
            gibbsSamplingStep(i);
        }
    }

    /**
     * @brief Get number of nodes
     */
    int getNumNodes() const {
        return numNodes;
    }

    /**
     * @brief Get number of states
     */
    int getNumStates() const {
        return numStates;
    }

    /**
     * @brief Get adjacency list for node
     */
    const std::vector<int>& getNeighbors(int node) const {
        if (node < 0 || node >= numNodes) {
            throw std::out_of_range("Node index out of range");
        }
        return adjList[node];
    }
};

// ============================================================================
// FINITE STATE MODELS
// ============================================================================

/**
 * @class FiniteStateModel
 * @brief Finite state machine with probabilistic transitions
 */
class FiniteStateModel {
private:
    int numStates;
    int currentState;
    std::vector<std::vector<double>> transitionMatrix;  // P(j|i)
    std::vector<double> emissionProbabilities;  // Current emission probabilities

public:
    /**
     * @brief Constructor
     *
     * @param states Number of states
     * @param initialState Initial state
     */
    FiniteStateModel(int states, int initialState = 0)
        : numStates(states), currentState(initialState) {
        if (states <= 0) {
            throw std::invalid_argument("Number of states must be positive");
        }
        if (initialState < 0 || initialState >= states) {
            throw std::out_of_range("Initial state out of range");
        }

        transitionMatrix.resize(states, std::vector<double>(states, 0.0));
        emissionProbabilities.resize(states, 0.0);
    }

    /**
     * @brief Set transition probability from state i to state j
     *
     * @param i From state
     * @param j To state
     * @param probability Transition probability
     */
    void setTransition(int i, int j, double probability) {
        if (i < 0 || i >= numStates || j < 0 || j >= numStates) {
            throw std::out_of_range("State index out of range");
        }
        if (probability < 0.0 || probability > 1.0) {
            throw std::invalid_argument("Probability must be in [0, 1]");
        }
        transitionMatrix[i][j] = probability;
    }

    /**
     * @brief Get transition probability
     */
    double getTransition(int i, int j) const {
        if (i < 0 || i >= numStates || j < 0 || j >= numStates) {
            throw std::out_of_range("State index out of range");
        }
        return transitionMatrix[i][j];
    }

    /**
     * @brief Normalize transition probabilities from state i
     */
    void normalizeTransitions(int i) {
        if (i < 0 || i >= numStates) {
            throw std::out_of_range("State index out of range");
        }

        double sum = 0.0;
        for (double prob : transitionMatrix[i]) {
            sum += prob;
        }

        if (sum > 1e-10) {
            for (double& prob : transitionMatrix[i]) {
                prob /= sum;
            }
        }
    }

    /**
     * @brief Sample next state from transition probabilities
     *
     * @return Next state
     */
    int sampleNextState() {
        double r = static_cast<double>(std::rand()) / RAND_MAX;
        double cumulative = 0.0;

        for (int j = 0; j < numStates; ++j) {
            cumulative += transitionMatrix[currentState][j];
            if (r <= cumulative) {
                currentState = j;
                return j;
            }
        }

        // Fallback (shouldn't happen if properly normalized)
        return currentState;
    }

    /**
     * @brief Calculate stationary distribution (eigenvector with eigenvalue 1)
     *
     * Simplified power iteration method
     *
     * @param iterations Number of iterations
     * @return Stationary distribution
     */
    std::vector<double> stationaryDistribution(int iterations = 1000) const {
        std::vector<double> distribution(numStates, 1.0 / numStates);

        for (int iter = 0; iter < iterations; ++iter) {
            std::vector<double> newDist(numStates, 0.0);

            for (int j = 0; j < numStates; ++j) {
                for (int i = 0; i < numStates; ++i) {
                    newDist[j] += distribution[i] * transitionMatrix[i][j];
                }
            }

            distribution = newDist;
        }

        return distribution;
    }

    /**
     * @brief Calculate n-step transition probability
     *
     * P^(n)(i, j) = (T^n)_{ij}
     *
     * @param i From state
     * @param j To state
     * @param n Number of steps
     * @return n-step transition probability
     */
    double nStepTransition(int i, int j, int n) const {
        if (n < 0) {
            throw std::invalid_argument("Number of steps must be non-negative");
        }
        if (n == 0) {
            return (i == j) ? 1.0 : 0.0;
        }

        // Matrix power using repeated multiplication
        std::vector<std::vector<double>> result = transitionMatrix;

        for (int step = 1; step < n; ++step) {
            std::vector<std::vector<double>> temp(numStates, std::vector<double>(numStates, 0.0));

            for (int row = 0; row < numStates; ++row) {
                for (int col = 0; col < numStates; ++col) {
                    for (int k = 0; k < numStates; ++k) {
                        temp[row][col] += result[row][k] * transitionMatrix[k][col];
                    }
                }
            }

            result = temp;
        }

        return result[i][j];
    }

    /**
     * @brief Get current state
     */
    int getCurrentState() const {
        return currentState;
    }

    /**
     * @brief Set current state
     */
    void setCurrentState(int state) {
        if (state < 0 || state >= numStates) {
            throw std::out_of_range("State out of range");
        }
        currentState = state;
    }

    /**
     * @brief Get number of states
     */
    int getNumStates() const {
        return numStates;
    }
};

// ============================================================================
// TREE MODEL (GRAPHICAL MODEL)
// ============================================================================

/**
 * @class TreeModel
 * @brief Tree-structured graphical model for belief propagation
 */
class TreeModel {
private:
    int numNodes;
    int numStates;
    std::vector<int> parent;  // Parent of each node (-1 for root)
    std::vector<std::vector<int>> children;  // Children of each node

    // Node potentials
    std::vector<std::vector<double>> nodePotentials;

    // Edge potentials: parent -> child
    std::vector<std::vector<std::vector<double>>> edgePotentials;

    // Messages: upward and downward
    std::vector<std::vector<double>> upwardMessages;
    std::vector<std::vector<double>> downwardMessages;

    // Marginals
    std::vector<std::vector<double>> marginals;

public:
    /**
     * @brief Constructor
     *
     * @param nodes Number of nodes
     * @param states Number of states per node
     */
    TreeModel(int nodes, int states)
        : numNodes(nodes), numStates(states) {
        if (nodes <= 0 || states <= 0) {
            throw std::invalid_argument("Nodes and states must be positive");
        }

        parent.resize(nodes, -1);
        children.resize(nodes);
        nodePotentials.resize(nodes, std::vector<double>(states, 1.0));
        edgePotentials.resize(nodes, std::vector<std::vector<double>>(
            states, std::vector<double>(states, 1.0)
        ));
        upwardMessages.resize(nodes, std::vector<double>(states, 1.0));
        downwardMessages.resize(nodes, std::vector<double>(states, 1.0));
        marginals.resize(nodes, std::vector<double>(states, 0.0));
    }

    /**
     * @brief Set parent-child relationship
     *
     * @param child Child node
     * @param par Parent node
     */
    void setParent(int child, int par) {
        if (child < 0 || child >= numNodes || par < -1 || par >= numNodes) {
            throw std::out_of_range("Node index out of range");
        }
        if (child == par) {
            throw std::invalid_argument("Node cannot be its own parent");
        }

        parent[child] = par;
        if (par >= 0) {
            children[par].push_back(child);
        }
    }

    /**
     * @brief Set node potential
     */
    void setNodePotential(int node, int state, double potential) {
        if (node < 0 || node >= numNodes) {
            throw std::out_of_range("Node index out of range");
        }
        if (state < 0 || state >= numStates) {
            throw std::out_of_range("State index out of range");
        }
        nodePotentials[node][state] = potential;
    }

    /**
     * @brief Set edge potential (parent state -> child state)
     */
    void setEdgePotential(int child, int parentState, int childState, double potential) {
        if (child < 0 || child >= numNodes) {
            throw std::out_of_range("Node index out of range");
        }
        if (parentState < 0 || parentState >= numStates ||
            childState < 0 || childState >= numStates) {
            throw std::out_of_range("State index out of range");
        }
        edgePotentials[child][parentState][childState] = potential;
    }

    /**
     * @brief Perform upward message passing (leaf to root)
     *
     * m_{i→par(i)}(x_par) = Σ_{x_i} ψ(x_par, x_i)·φ_i(x_i)·∏_{j∈ch(i)} m_{j→i}(x_i)
     */
    void upwardPass(int node) {
        // Recursively compute messages from children first
        for (int child : children[node]) {
            upwardPass(child);
        }

        // Compute message to parent (if not root)
        if (parent[node] >= 0) {
            std::fill(upwardMessages[node].begin(), upwardMessages[node].end(), 0.0);

            for (int parentState = 0; parentState < numStates; ++parentState) {
                for (int nodeState = 0; nodeState < numStates; ++nodeState) {
                    double product = edgePotentials[node][parentState][nodeState] *
                                   nodePotentials[node][nodeState];

                    // Multiply messages from children
                    for (int child : children[node]) {
                        product *= upwardMessages[child][nodeState];
                    }

                    upwardMessages[node][parentState] += product;
                }
            }
        }
    }

    /**
     * @brief Perform downward message passing (root to leaf)
     */
    void downwardPass(int node) {
        // Compute message from parent (if not root)
        if (parent[node] >= 0) {
            int par = parent[node];
            std::fill(downwardMessages[node].begin(), downwardMessages[node].end(), 0.0);

            for (int nodeState = 0; nodeState < numStates; ++nodeState) {
                for (int parentState = 0; parentState < numStates; ++parentState) {
                    double product = edgePotentials[node][parentState][nodeState] *
                                   nodePotentials[par][parentState] *
                                   downwardMessages[par][parentState];

                    // Multiply messages from siblings
                    for (int sibling : children[par]) {
                        if (sibling != node) {
                            product *= upwardMessages[sibling][parentState];
                        }
                    }

                    downwardMessages[node][nodeState] += product;
                }
            }
        } else {
            // Root node
            std::fill(downwardMessages[node].begin(), downwardMessages[node].end(), 1.0);
        }

        // Recursively process children
        for (int child : children[node]) {
            downwardPass(child);
        }
    }

    /**
     * @brief Compute marginal probabilities for all nodes
     *
     * @param root Root node
     */
    void computeMarginals(int root = 0) {
        // Upward pass
        upwardPass(root);

        // Downward pass
        downwardPass(root);

        // Compute marginals
        for (int node = 0; node < numNodes; ++node) {
            double sum = 0.0;

            for (int state = 0; state < numStates; ++state) {
                marginals[node][state] = nodePotentials[node][state] *
                                        downwardMessages[node][state];

                // Multiply messages from children
                for (int child : children[node]) {
                    marginals[node][state] *= upwardMessages[child][state];
                }

                sum += marginals[node][state];
            }

            // Normalize
            if (sum > 1e-10) {
                for (int state = 0; state < numStates; ++state) {
                    marginals[node][state] /= sum;
                }
            }
        }
    }

    /**
     * @brief Get marginal probability of state at node
     */
    double getMarginal(int node, int state) const {
        if (node < 0 || node >= numNodes) {
            throw std::out_of_range("Node index out of range");
        }
        if (state < 0 || state >= numStates) {
            throw std::out_of_range("State index out of range");
        }
        return marginals[node][state];
    }

    /**
     * @brief Get most likely state for node (MAP estimate)
     */
    int getMAPState(int node) const {
        if (node < 0 || node >= numNodes) {
            throw std::out_of_range("Node index out of range");
        }

        int maxState = 0;
        double maxProb = marginals[node][0];

        for (int state = 1; state < numStates; ++state) {
            if (marginals[node][state] > maxProb) {
                maxProb = marginals[node][state];
                maxState = state;
            }
        }

        return maxState;
    }

    /**
     * @brief Get number of nodes
     */
    int getNumNodes() const {
        return numNodes;
    }

    /**
     * @brief Get number of states
     */
    int getNumStates() const {
        return numStates;
    }
};

} // namespace statistical_models
} // namespace physics

#endif // PHYSICS_STATISTICAL_MODELS_HPP
