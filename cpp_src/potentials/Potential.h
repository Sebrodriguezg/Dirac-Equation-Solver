#pragma once

#include <vector>
#include <string>
#include <memory>

/**
 * @class Potential
 * @brief Abstract base class for all scalar potentials (Template pattern).
 *
 * This class defines the interface for any potential V(r) to be used in
 * the Dirac equation. It allows different potential functions to be
 * implemented and used interchangeably by the solver.
 */
class Potential {
public:
    virtual ~Potential() = default;

    /**
     * @brief Evaluates the potential at a given spatial coordinate.
     * @param r A vector representing the spatial coordinate (e.g., [x, y, z]).
     * @return The scalar value of the potential at that point.
     */
    virtual double evaluate(const std::vector<double>& r) const = 0;

    /**
     * @brief Returns the name of the potential.
     * @return A string containing the name of the class.
     */
    virtual std::string get_name() const = 0;
};