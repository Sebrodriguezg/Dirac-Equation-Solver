#pragma once

#include <vector>
#include <string>
#include <memory>
#include "../DiracMatrices.h" // For Dirac::Spinor

/**
 * @class BoundaryCondition
 * @brief Abstract base class for boundary conditions (Strategy pattern).
 *
 * This class defines the interface for all boundary condition strategies.
 * It allows the core solver to be decoupled from the specific method
 * used to handle the edges of the simulation grid.
 */
class BoundaryCondition {
public:
    virtual ~BoundaryCondition() = default;

    /**
     * @brief Applies the boundary condition to a given spinor field.
     * @param psi A reference to the vector of spinors representing the field.
     */
    virtual void apply(std::vector<Dirac::Spinor>& psi) const = 0;

    /**
     * @brief Returns the name of the boundary condition strategy.
     * @return A string containing the name of the class.
     */
    virtual std::string get_name() const = 0;
};