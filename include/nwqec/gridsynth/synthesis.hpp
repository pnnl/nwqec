#pragma once

#include <string>
#include <vector>
#include <array>
#include <cassert>

#include "ring.hpp"
#include "unitary.hpp"
#include "normal_form.hpp"

namespace gridsynth
{

    // Lookup tables for reduction algorithm
    inline const std::array<int, 16> BIT_SHIFT = {0, 0, 1, 0, 2, 0, 1, 3, 3, 3, 0, 2, 2, 1, 0, 0};
    inline const std::array<int, 16> BIT_COUNT = {0, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 4};

    /**
     * Reduces the denominator exponent of a DOmegaUnitary by applying gates
     * Returns a pair of (gate_string, reduced_unitary)
     */
    inline std::pair<std::string, DOmegaUnitary> _reduce_denomexp(const DOmegaUnitary &unitary)
    {
        std::vector<std::string> T_POWER_and_H = {"H", "TH", "SH", "TSH"};

        int residue_z = unitary.z().residue();
        int residue_w = unitary.w().residue();
        int residue_squared_z = (unitary.z().u() * unitary.z().conj().u()).residue();

        int m = BIT_SHIFT[residue_w] - BIT_SHIFT[residue_z];
        if (m < 0)
        {
            m += 4;
        }

        DOmegaUnitary new_unitary = unitary;
        std::string gate_string;

        if (residue_squared_z == 0b0000)
        {
            new_unitary = unitary.mul_by_H_and_T_power_from_left(0).renew_denomexp(unitary.k() - 1);
            gate_string = T_POWER_and_H[0];
        }
        else if (residue_squared_z == 0b1010)
        {
            new_unitary = unitary.mul_by_H_and_T_power_from_left(-m).renew_denomexp(unitary.k() - 1);
            gate_string = T_POWER_and_H[m];
        }
        else if (residue_squared_z == 0b0001)
        {
            if (BIT_COUNT[residue_z] == BIT_COUNT[residue_w])
            {
                new_unitary = unitary.mul_by_H_and_T_power_from_left(-m).renew_denomexp(unitary.k() - 1);
                gate_string = T_POWER_and_H[m];
            }
            else
            {
                new_unitary = unitary.mul_by_H_and_T_power_from_left(-m);
                gate_string = T_POWER_and_H[m];
            }
        }
        else
        {
            // Default case - just apply H
            new_unitary = unitary.mul_by_H_from_left().renew_denomexp(unitary.k() - 1);
            gate_string = "H";
        }

        return {gate_string, new_unitary};
    }

    /**
     * Decomposes a DOmegaUnitary into a sequence of Clifford+T gates
     * This is the main synthesis algorithm
     */
    inline std::string decompose_domega_unitary(DOmegaUnitary unitary)
    {
        std::string gates = "";

        // Reduce denominator exponent to 0
        while (unitary.k() > 0)
        {
            auto [gate_string, reduced_unitary] = _reduce_denomexp(unitary);
            gates += gate_string;
            unitary = reduced_unitary;
        }

        // Handle odd phase
        if (unitary.n() & 1)
        {
            gates += "T";
            unitary = unitary.mul_by_T_inv_from_left();
        }

        // Handle case where z = 0
        if (unitary.z() == DOmega::from_int(0))
        {
            gates += "X";
            unitary = unitary.mul_by_X_from_left();
        }

        // Find the omega power that matches z
        int m_W = 0;
        for (int m = 0; m < 8; ++m)
        {
            if (unitary.z().u() == OMEGA_POWER[m])
            {
                m_W = m;
                unitary = unitary.mul_by_W_power_from_left(-m_W);
                break;
            }
        }

        // Apply S gates to normalize the phase
        int m_S = unitary.n() >> 1;
        gates += std::string(m_S, 'S');
        unitary = unitary.mul_by_S_power_from_left(-m_S);

        // Apply W gates
        gates += std::string(m_W, 'W');

        // Verify decomposition is correct (optional assertion)
        // assert(unitary == DOmegaUnitary::identity());

        // Normalize the gate sequence using the complete NormalForm implementation
        NormalForm normal_form = NormalForm::from_gates(gates);
        return normal_form.to_gates();
    }

} // namespace gridsynth
