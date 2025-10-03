#pragma once

#include "Potential.h"

/**
 * @class FreeParticle
 * @brief A concrete implementation of Potential for a free particle (V=0).
 *
 * This class represents the case where there is no external potential.
 * It fulfills the requirement FR-3 for a free particle.
 */
class FreeParticle : public Potential {
public:
    /**
     * @brief Evaluates the potential, which is always zero for a free particle.
     * @param r The spatial coordinate (ignored).
     * @return Always returns 0.0.
     */
    double evaluate(const std::vector<double>& r) const override {
        // For a free particle, the potential is always zero.
        (void)r; // Mark 'r' as unused to prevent compiler warnings.
        return 0.0;
    }

    /**
     * @brief Returns the name of the potential.
     * @return The string "FreeParticle".
     */
    std::string get_name() const override {
        return "FreeParticle";
    }
};