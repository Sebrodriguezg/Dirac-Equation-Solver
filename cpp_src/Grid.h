#pragma once

#include <vector>
#include <numeric>
#include <stdexcept>
#include <cmath>

/**
 * @class Grid
 * @brief Represents a spatial grid (1D, 2D, or 3D) for the simulation.
 *
 * This class stores the parameters of the simulation grid, such as its
 * dimensions, number of points, and cell spacing. It is a C++ counterpart
 * to the Python Grid class.
 */
class Grid {
public:
    /**
     * @brief Constructs a Grid object.
     * @param shape A vector containing the number of points in each dimension (e.g., {nx, ny}).
     * @param spacing A vector containing the cell spacing in each dimension (e.g., {dx, dy}).
     */
    Grid(std::vector<size_t> shape, std::vector<double> spacing)
        : shape_(shape), spacing_(spacing), dim_(shape.size()) {
        if (shape.size() != spacing.size()) {
            throw std::invalid_argument("Shape and spacing must have the same number of dimensions.");
        }
        if (dim_ < 1 || dim_ > 3) {
            throw std::invalid_argument("Grid must be 1D, 2D, or 3D.");
        }
    }

    // Getters for grid properties
    size_t get_dim() const { return dim_; }
    const std::vector<size_t>& get_shape() const { return shape_; }
    const std::vector<double>& get_spacing() const { return spacing_; }

    /**
     * @brief Calculates the total number of points in the grid.
     * @return The product of the number of points in each dimension.
     */
    size_t get_total_points() const {
        size_t total = 1;
        for (size_t n : shape_) {
            total *= n;
        }
        return total;
    }

private:
    std::vector<size_t> shape_;  // Number of points in each dimension
    std::vector<double> spacing_;// Cell spacing in each dimension
    size_t dim_;                 // Number of dimensions
};