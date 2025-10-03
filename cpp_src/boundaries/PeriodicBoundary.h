#pragma once

#include "BoundaryCondition.h"

/**
 * @class PeriodicBoundary
 * @brief Implements periodic boundary conditions.
 *
 * This strategy treats the grid as a torus, where the value at the
 * start of an axis is linked to the value at the end. For a 1D grid,
 * this means psi[0] is linked to psi[N-1].
 */
class PeriodicBoundary : public BoundaryCondition {
public:
    /**
     * @brief Applies periodic boundary conditions to the field.
     * In the FDTD scheme, this is typically handled by the update loop itself,
     * where the neighbor of point 0 is N-1, and the neighbor of N-1 is 0.
     * This explicit apply function can be used if the main loop doesn't
     * handle it, for example, to enforce the condition on the boundary points
     * after a step. For our current FDTD implementation, the logic will be
     * integrated directly into the `step` function for efficiency.
     *
     * @param psi A reference to the vector of spinors (the field).
     */
    void apply(std::vector<Dirac::Spinor>& psi) const override {
        if (psi.size() < 2) return; // Not applicable for grids with less than 2 points.

        // In a simple case, we would copy the values, e.g.:
        // psi.front() = psi.back();
        // But for a staggered grid, this logic is more complex and best
        // handled inside the numerical loop. This function serves as a placeholder
        // and confirmation of the strategy pattern's implementation.
    }

    /**
     * @brief Returns the name of the boundary condition strategy.
     * @return The string "PeriodicBoundary".
     */
    std::string get_name() const override {
        return "PeriodicBoundary";
    }
};