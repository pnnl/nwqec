#pragma once

#include <complex>
#include <vector>
#include <string>
#include <array>
#include <cmath>
#include <chrono>
#include <iostream>
#include <algorithm>
#include <sstream>
#include <memory>
#include <optional>
#include "nwqec/core/constants.hpp"

#include "nwqec/gridsynth/ring.hpp"
#include "nwqec/gridsynth/grid_op.hpp"
#include "nwqec/gridsynth/region.hpp"
#include "nwqec/gridsynth/epsilon_region.hpp"
#include "nwqec/gridsynth/unitary.hpp"
#include "nwqec/gridsynth/to_upright.hpp"
#include "nwqec/gridsynth/tdgp.hpp"
#include "nwqec/gridsynth/diophantine.hpp"
#include "nwqec/gridsynth/synthesis.hpp"
#include "nwqec/gridsynth/mymath.hpp"
#include "nwqec/gridsynth/types.hpp"

namespace gridsynth
{

    /**
     * Calculate error between target rotation and synthesized unitary
     */
    inline std::string error(std::string theta, const std::string &gates)
    {
        // Compute Rz(theta) and U coordinates directly to avoid building std::complex matrix
        Float theta_f = Float(theta) / Float(2.0);
        Float cos_theta = cos(theta_f);
        Float sin_theta = sin(theta_f);

        DOmegaUnitary Uu = DOmegaUnitary::from_gates(gates);
        auto M = Uu.to_matrix(); // 2x2 DOmega
        const Float inv_scale = Float(1.0) / pow_sqrt2(Uu.k());
        const Float sqrt2_over_2 = SQRT2 / Float(2.0);

        // Extract U entries
        Float u00r, u00i, u01r, u01i, u10r, u10i, u11r, u11i;
        M[0][0].coords_into_with(inv_scale, sqrt2_over_2, u00r, u00i);
        M[0][1].coords_into_with(inv_scale, sqrt2_over_2, u01r, u01i);
        M[1][0].coords_into_with(inv_scale, sqrt2_over_2, u10r, u10i);
        M[1][1].coords_into_with(inv_scale, sqrt2_over_2, u11r, u11i);

        // E = U - Rz
        Float A = u00r - cos_theta;
        Float B = u00i + sin_theta;
        Float C = u11r - cos_theta;
        Float D = u11i - sin_theta;
        Float E = u01r;
        Float F = u01i;
        Float G = u10r;
        Float H = u10i;

        // det(E) = e00*e11 - e01*e10
        Float det_re = (A * C - B * D) - (E * G - F * H);
        Float det_im = (A * D + B * C) - (E * H + F * G);
        Float abs_det = sqrt(det_re * det_re + det_im * det_im);
        return sqrt(abs_det).to_string();
    }

    /**
     * Main gridsynth algorithm - finds a DOmegaUnitary approximation
     *
     * @param theta Target rotation angle
     * @param epsilon Error tolerance
     * @param diophantine_timeout_ms Timeout for diophantine solving in milliseconds
     * @param factoring_timeout_ms Timeout for factoring in milliseconds
     * @param verbose Enable verbose output
     * @param measure_time Enable timing measurements
     * @return DOmegaUnitary approximation
     */
    inline DOmegaUnitary gridsynth(
        Float theta,
        Float epsilon,
        int diophantine_timeout_ms = NWQEC::DEFAULT_DIOPHANTINE_TIMEOUT_MS,
        int factoring_timeout_ms = NWQEC::DEFAULT_FACTORING_TIMEOUT_MS,
        bool verbose = false,
        bool measure_time = false)
    {

        double time_of_solve_TDGP = 0.0;
        double time_of_diophantine_dyadic = 0.0;
        double time_of_to_upright = 0.0;

        // Create proper EpsilonRegion and UnitDisk objects
        EpsilonRegion epsilon_region(theta, epsilon);
        UnitDisk unit_disk;

        // Transform to upright position
        auto start = std::chrono::high_resolution_clock::now();

        ToUpright transformer(epsilon_region.ellipse(), unit_disk.ellipse());
        transformer.run(verbose);
        auto transformed = transformer.get_result();

        // exit(0); // (was present previously; keeping commented out for safety)

        if (measure_time)
        {
            auto end = std::chrono::high_resolution_clock::now();
            auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
            time_of_to_upright += duration.count() / 1000.0;
        }

        if (verbose)
        {
            std::cout << "------------------" << std::endl;
        }

        // fatten y intervals by width * 1e-4
        Float epsilon_factor = Float(1e-4);
        Interval bboxA_y_fattened = transformed.bboxA.I_y().fatten(transformed.bboxA.I_y().width() * epsilon_factor);
        Interval bboxB_y_fattened = transformed.bboxB.I_y().fatten(transformed.bboxB.I_y().width() * epsilon_factor);
        GridOp opG_inv = transformed.opG.inv();

        int num_diophantine_calls = 0;

        TdgpSolver gp_solver(epsilon_region,
                             unit_disk,
                             opG_inv,
                             transformed.bboxA, transformed.bboxB,
                             bboxA_y_fattened, bboxB_y_fattened);
        Integer k = 0;

        while (true) // Use infinite loop like Python version
        {
            // Solve TDGP
            if (measure_time)
            {
                start = std::chrono::high_resolution_clock::now();
            }

            auto sol = gp_solver.solve(k, verbose);

            if (measure_time)
            {
                auto end = std::chrono::high_resolution_clock::now();
                auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
                time_of_solve_TDGP += duration.count() / 1000.0;
                start = std::chrono::high_resolution_clock::now();
            }

            // Try each solution from TDGP
            for (const DOmega &z : sol)
            {
                if ((z * z.conj()).residue() == 0)
                {
                    continue;
                }

                DRootTwo xi = DRootTwo(1) - DRootTwo::fromDOmega(z.conj() * z);
                std::optional<DOmega> w_opt = diophantine_dyadic(xi, diophantine_timeout_ms, factoring_timeout_ms);
                num_diophantine_calls++;
                if (w_opt.has_value())
                {

                    DOmega z_reduced = z.reduce_denomexp();
                    DOmega w_reduced = w_opt->reduce_denomexp();

                    // Align denominator exponents
                    if (z_reduced.k() > w_reduced.k())
                    {
                        w_reduced = w_reduced.renew_denomexp(z_reduced.k());
                    }
                    else if (z_reduced.k() < w_reduced.k())
                    {
                        z_reduced = z_reduced.renew_denomexp(w_reduced.k());
                    }

                    DOmegaUnitary u_approx(DOmega::from_int(0), DOmega::from_int(0), 0); // Initialize with dummy values
                    if ((z_reduced + w_reduced).reduce_denomexp().k() < z_reduced.k())
                    {
                        u_approx = DOmegaUnitary(z_reduced, w_reduced, 0);
                    }
                    else
                    {
                        u_approx = DOmegaUnitary(z_reduced, w_reduced.mul_by_omega(), 0);
                    }

                    if (measure_time)
                    {
                        auto end = std::chrono::high_resolution_clock::now();
                        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
                        time_of_diophantine_dyadic += duration.count() / 1000.0;

                        std::cout << "time of to_upright: " << time_of_to_upright << " ms" << std::endl;
                        std::cout << "time of solve_TDGP: " << time_of_solve_TDGP << " ms" << std::endl;
                        std::cout << "time of diophantine(" << num_diophantine_calls << "): " << time_of_diophantine_dyadic << " ms" << std::endl;
                    }

                    if (verbose)
                    {
                        std::cout << "z=" << z_reduced.to_string()
                                  << ", w=" << w_reduced.to_string() << std::endl;
                        std::cout << "------------------" << std::endl;
                    }

                    return u_approx;
                }
            }

            if (measure_time)
            {
                auto end = std::chrono::high_resolution_clock::now();
                auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
                time_of_diophantine_dyadic += duration.count() / 1000.0;
            }

            k++;
        }
    }

    /**
     * Main gridsynth algorithm - returns gate sequence
     *
     * @param theta Target rotation angle
     * @param epsilon Error tolerance
     * @param diophantine_timeout_ms Timeout for diophantine solving in milliseconds
     * @param factoring_timeout_ms Timeout for factoring in milliseconds
     * @param verbose Enable verbose output
     * @param measure_time Enable timing measurements
     * @return Vector of gate strings
     */
    inline std::string gridsynth_gates(
        std::string theta,
        std::string epsilon,
        int diophantine_timeout_ms = NWQEC::DEFAULT_DIOPHANTINE_TIMEOUT_MS,
        int factoring_timeout_ms = NWQEC::DEFAULT_FACTORING_TIMEOUT_MS,
        bool verbose = false,
        bool measure_time = false)
    {
        auto start_total = std::chrono::high_resolution_clock::now();

        DOmegaUnitary u_approx = gridsynth(
            Float(theta), Float(epsilon),
            diophantine_timeout_ms, factoring_timeout_ms,
            verbose, measure_time);

        std::string gates_str = decompose_domega_unitary(u_approx);

        if (measure_time)
        {
            auto end = std::chrono::high_resolution_clock::now();
            auto total_duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start_total);

            std::cout << "Gridsynth time: " << total_duration.count() / 1000.0 << " ms" << std::endl;
        }

        return gates_str;
    }

} // namespace gridsynth
