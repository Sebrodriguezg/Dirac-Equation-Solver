#pragma once

#include <iostream>
#include <string>
#include <vector>
#include <stdexcept>
#include <memory>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

#include "Grid.h"
#include "potentials/Potential.h"
#include "boundaries/BoundaryCondition.h"
#include "DiracMatrices.h"

namespace py = pybind11;

/**
 * @class TimeIntegrator
 * @brief Abstract base class for time-evolution algorithms (Strategy pattern).
 *
 * This class defines the interface for all time integrator strategies.
 * It allows the core solver to be decoupled from the specific numerical
 * method used to advance the simulation in time.
 */
class TimeIntegrator {
public:
    virtual ~TimeIntegrator() = default;

    /**
     * @brief Performs a single time step to advance the simulation.
     * This is the core method of the strategy.
     */
    virtual void step() = 0;

    /**
     * @brief Returns the name of the integration strategy.
     * @return A string containing the name of the class.
     */
    virtual std::string get_name() const = 0;
};

/**
 * @class FDTDLeapfrogIntegrator
 * @brief A concrete implementation of TimeIntegrator using the FDTD Leapfrog method.
 *
 * This class will contain the core logic for the staggered-leapfrog FDTD algorithm.
 */
class FDTDLeapfrogIntegrator : public TimeIntegrator {
public:
    /**
     * @brief Constructs the FDTD Leapfrog integrator.
     * @param initial_psi The initial state of the spinor field as a NumPy array of shape (N, 4).
     * @param grid The simulation grid object.
     * @param potential The potential object.
     * @param dt The time step.
     * @param mass The mass of the particle.
     */
    FDTDLeapfrogIntegrator(
        py::array_t<std::complex<double>, py::array::c_style | py::array::forcecast> initial_psi,
        const Grid& grid,
        std::shared_ptr<Potential> potential,
        std::shared_ptr<BoundaryCondition> boundary_condition,
        double dt,
        double mass)
        : grid_(grid), potential_(potential), boundary_condition_(boundary_condition), dt_(dt), mass_(mass)
    {
        // 1. Validate input array dimensions
        py::buffer_info psi_info = initial_psi.request();
        if (psi_info.ndim != 2 || psi_info.shape[1] != 4) {
            throw std::runtime_error("Initial psi must be a NumPy array of shape (N, 4)");
        }
        if (static_cast<size_t>(psi_info.shape[0]) != grid_.get_total_points()) {
            throw std::runtime_error("Number of points in psi must match total grid points.");
        }

        // 2. Allocate memory for the spinor fields
        size_t n_points = grid_.get_total_points();
        psi_curr_.resize(n_points);
        psi_prev_.resize(n_points);
        psi_next_.resize(n_points);

        // 3. Copy initial data from NumPy array to psi_curr_
        auto ptr = static_cast<Dirac::complex*>(psi_info.ptr);
        for (size_t i = 0; i < n_points; ++i) {
            for (size_t j = 0; j < 4; ++j) {
                psi_curr_[i].components[j] = ptr[i * 4 + j];
            }
        }

        // 4. Perform the first step using a simpler, non-centered scheme (Forward Euler)
        // to generate the state at time t=dt, which is needed for the leapfrog to start.
        // This is a common way to initialize a leapfrog integrator.
        // For now, we just copy the current state to the previous state.
        // A more accurate initializer would be needed for a real simulation.
        psi_prev_ = psi_curr_;
    }

    void step() override {
        if (grid_.get_dim() != 1) {
            throw std::runtime_error("FDTD step is only implemented for 1D grids at the moment.");
        }

        size_t n_points = grid_.get_total_points();
        double dz = grid_.get_spacing()[0];

        // The core FDTD leapfrog update loop with periodic boundaries.
        for (size_t i = 0; i < n_points; ++i) {
            // Get indices of neighbors with periodic wrapping
            size_t i_prev = (i == 0) ? n_points - 1 : i - 1;
            size_t i_next = (i == n_points - 1) ? 0 : i + 1;

            // H_spatial * psi_curr = -i * alpha_z * d/dz(psi_curr)
            Dirac::Spinor h_spatial_psi = Dirac::multiply(Dirac::alpha_z, psi_curr_[i_next] - psi_curr_[i_prev]);
            h_spatial_psi = (Dirac::complex(0, -1.0) / (2.0 * dz)) * h_spatial_psi;

            // H_local * psi_curr = (beta*m + V) * psi_curr
            Dirac::Spinor h_local_psi = mass_ * Dirac::multiply(Dirac::beta, psi_curr_[i]);
            // Note: Potential term is missing, will be added next.

            // Leapfrog formula: psi_next = psi_prev - 2*i*dt*H*psi_curr
            Dirac::Spinor h_psi = h_spatial_psi + h_local_psi;
            psi_next_[i] = psi_prev_[i] - (Dirac::complex(0, 2.0) * dt_) * h_psi;
        }

        // Update the fields for the next iteration
        psi_prev_ = psi_curr_;
        psi_curr_ = psi_next_;
    }

    std::string get_name() const override {
        return "FDTDLeapfrogIntegrator";
    }

    // Getter to return the current wave function to Python
    py::array_t<std::complex<double>> get_psi() const {
        py::array_t<std::complex<double>> result({grid_.get_total_points(), 4ul});
        auto buf = result.request();
        auto* ptr = static_cast<std::complex<double>*>(buf.ptr);
        for (size_t i = 0; i < grid_.get_total_points(); ++i) {
            for (size_t j = 0; j < 4; ++j) {
                ptr[i * 4 + j] = psi_curr_[i].components[j];
            }
        }
        return result;
    }

private:
    Grid grid_;
    std::shared_ptr<Potential> potential_;
    std::shared_ptr<BoundaryCondition> boundary_condition_;
    double dt_;
    double mass_;

    // Three time-step fields for the leapfrog algorithm
    std::vector<Dirac::Spinor> psi_prev_;
    std::vector<Dirac::Spinor> psi_curr_;
    std::vector<Dirac::Spinor> psi_next_;
};