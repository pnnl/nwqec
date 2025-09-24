#pragma once

#include <cmath>
#include <algorithm>
#include <iostream>

#include "grid_op.hpp"
#include "region.hpp"
#include "ring.hpp"
#include "mymath.hpp"

/**
 * Object-oriented version of to_upright.py - class for transforming ellipse pairs to "upright" position
 * for optimal quantum circuit synthesis
 */
namespace gridsynth
{

    /**
     * Apply a reduction step: new_opG * ellipse_pair, opG_l, new_opG * opG_r, False
     */
    inline void _reduction(EllipsePair &ellipse_pair,
                           GridOp &opG_r, const GridOp &new_opG)
    {
        ellipse_pair.apply_grid_op(new_opG);
        opG_r = new_opG * opG_r;
    }

    /**
     * Shift ellipse pair by lambda^n scaling
     */
    void _shift_ellipse_pair(EllipsePair &ellipse_pair, Integer n)
    {
        // Get lambda^n and lambda^(-n)
        ZRootTwo lambda_n = LAMBDA.pow(n);
        ZRootTwo lambda_inv_n = lambda_n.inv();

        Float lambda_n_real = lambda_n.to_real();
        Float lambda_inv_n_real = lambda_inv_n.to_real();

        // Get mutable copies of the ellipses
        Ellipse A = ellipse_pair.e1();
        Ellipse B = ellipse_pair.e2();

        // Apply scaling: A.a *= lambda^(-n), A.d *= lambda^n
        A.scale_a(lambda_inv_n_real);
        A.scale_d(lambda_n_real);

        // Apply scaling: B.a *= lambda^n, B.d *= lambda^(-n)
        B.scale_a(lambda_n_real);
        B.scale_d(lambda_inv_n_real);

        // If n is odd, flip B.b sign
        if (n.is_odd())
        {
            B.flip_b();
        }

        ellipse_pair.set_e1(A);
        ellipse_pair.set_e2(B);
    }

    /**
     * Main step lemma - decides which transformation to apply
     */
    void _step_lemma(EllipsePair &ellipse_pair, GridOp &opG_l,
                     GridOp &opG_r, bool &end, bool verbose = false)
    {
        const Ellipse &A = ellipse_pair.e1();
        const Ellipse &B = ellipse_pair.e2();

        // Precompute values to avoid repeated calls
        Float ellipse_pair_bias = ellipse_pair.bias();
        Float ellipse_pair_skew = ellipse_pair.skew();
        Float bias_A = A.bias();
        Float bias_B = B.bias();

        if (verbose)
        {
            std::cout << "-----\n";
            std::cout << "skew: " << ellipse_pair_skew
                      << ", bias: " << ellipse_pair_bias << "\n";
            std::cout << "bias(A): " << bias_A << ", bias(B): " << bias_B
                      << ", sign(A.b):" << (A.b() >= 0 ? "+" : "-")
                      << ", sign(B.b):" << (B.b() >= 0 ? "+" : "-") << "\n";
            std::cout << "-----\n";
        }

        // Z operation: if B.b < 0
        if (B.b() < 0)
        {
            if (verbose)
                std::cout << "Z\n";
            GridOp OP_Z(ZOmega(0, 0, 0, 1), ZOmega(0, -1, 0, 0));
            _reduction(ellipse_pair, opG_r, OP_Z);
            end = false;
            return;
        }

        // X operation: if A.bias * B.bias < 1
        else if (bias_A * bias_B < 1)
        {
            if (verbose)
                std::cout << "X\n";
            GridOp OP_X(ZOmega(0, 1, 0, 0), ZOmega(0, 0, 0, 1));
            _reduction(ellipse_pair, opG_r, OP_X);
            end = false;
            return;
        }

        // S operation: extreme bias values
        if (ellipse_pair_bias > 33.971 || ellipse_pair_bias < 0.029437)
        {
            Float lambda_real = LAMBDA.to_real();
            Integer n = round_to_integer(log(ellipse_pair_bias) / log(lambda_real) / 8);
            if (verbose)
                std::cout << "S (n=" << n << ")\n";
            GridOp OP_S(ZOmega(-1, 0, 1, 1), ZOmega(1, -1, 1, 0));
            _reduction(ellipse_pair, opG_r, OP_S.pow(n));
            end = false;
            return;
        }

        // Check if we're done: skew <= 15
        if (ellipse_pair_skew <= 15)
        {
            end = true;
            return;
        }

        // Sigma operation: moderate bias values
        if (ellipse_pair_bias > 5.8285 || ellipse_pair_bias < 0.17157)
        {
            Float lambda_real = LAMBDA.to_real();
            Integer n = round_to_integer(log(ellipse_pair_bias) / log(lambda_real) / 4);
            if (verbose)
                std::cout << "sigma (n=" << n << ")\n";

            _shift_ellipse_pair(ellipse_pair, n);

            if (n >= 0)
            {
                GridOp OP_SIGMA_L = GridOp(ZOmega(-1, 0, 1, 1), ZOmega(0, 1, 0, 0)).pow(n);
                GridOp OP_SIGMA_R = GridOp(ZOmega(0, 0, 0, 1), ZOmega(1, -1, 1, 0)).pow(n);
                opG_l = opG_l * OP_SIGMA_L;
                opG_r = OP_SIGMA_R * opG_r;
            }
            else
            {
                GridOp OP_SIGMA_L = GridOp(ZOmega(-1, 0, 1, -1), ZOmega(0, 1, 0, 0)).pow(-n);
                GridOp OP_SIGMA_R = GridOp(ZOmega(0, 0, 0, 1), ZOmega(1, 1, 1, 0)).pow(-n);
                opG_l = opG_l * OP_SIGMA_L;
                opG_r = OP_SIGMA_R * opG_r;
            }
            end = false;
            return;
        }

        // R operation: both biases in moderate range
        if (0.24410 <= bias_A && bias_A <= 4.0968 &&
            0.24410 <= bias_B && bias_B <= 4.0968)
        {
            if (verbose)
                std::cout << "R\n";
            GridOp OP_R(ZOmega(0, 0, 1, 0), ZOmega(1, 0, 0, 0));
            _reduction(ellipse_pair, opG_r, OP_R);
            end = false;
            return;
        }

        // K operation: A.b >= 0 and A.bias <= 1.6969
        if (A.b() >= 0 && bias_A <= 1.6969)
        {
            if (verbose)
                std::cout << "K\n";
            GridOp OP_K(ZOmega(-1, -1, 0, 0), ZOmega(0, -1, 1, 0));
            _reduction(ellipse_pair, opG_r, OP_K);
            end = false;
            return;
        }

        // K_conj_sq2 operation: A.b >= 0 and B.bias <= 1.6969
        if (A.b() >= 0 && bias_B <= 1.6969)
        {
            if (verbose)
                std::cout << "K_conj_sq2\n";
            GridOp OP_K_conj_sq2(ZOmega(1, -1, 0, 0), ZOmega(0, -1, -1, 0));
            _reduction(ellipse_pair, opG_r, OP_K_conj_sq2);
            end = false;
            return;
        }

        // A operation: A.b >= 0
        if (A.b() >= 0)
        {
            Integer n = max(Integer(1), floorsqrt(min(bias_A, bias_B) / 4));
            if (verbose)
                std::cout << "A (n=" << n << ")\n";
            GridOp OP_A_n(ZOmega(0, 0, 0, 1), ZOmega(0, 1, 0, 2 * n));
            _reduction(ellipse_pair, opG_r, OP_A_n);
            end = false;
            return;
        }

        // B operation: fallback case
        {
            Integer n = max(Integer(1), floorsqrt(min(bias_A, bias_B) / 2));
            if (verbose)
                std::cout << "B (n=" << n << ")\n";
            GridOp OP_B_n(ZOmega(0, 0, 0, 1), ZOmega(n, 1, -n, 0));
            _reduction(ellipse_pair, opG_r, OP_B_n);
            end = false;
            return;
        }
    }

    // Result structure for to_upright_set_pair
    struct UprightResult
    {
        GridOp opG;
        Rectangle bboxA;
        Rectangle bboxB;

        UprightResult(const GridOp &op, const Rectangle &bA, const Rectangle &bB)
            : opG(op), bboxA(bA), bboxB(bB) {}
    };

    /**
     * Object-oriented class for transforming ellipse pairs to upright position
     */
    class ToUpright
    {
    private:
        Ellipse original_setA_;
        Ellipse original_setB_;
        Ellipse setA_;
        Ellipse setB_;
        EllipsePair original_pair_;
        EllipsePair current_pair_;
        GridOp opG_l_;
        GridOp opG_r_;
        bool done_;

    public:
        /**
         * Constructor: Initialize with copies of setA and setB
         */
        ToUpright(const Ellipse &setA, const Ellipse &setB)
            : original_setA_(setA), original_setB_(setB),
              setA_(setA.normalize()), setB_(setB.normalize()),
              original_pair_(original_setA_, original_setB_), current_pair_(setA_, setB_),
              opG_l_(GridOp(ZOmega(0, 0, 0, 1), ZOmega(0, 1, 0, 0))), // Identity
              opG_r_(GridOp(ZOmega(0, 0, 0, 1), ZOmega(0, 1, 0, 0))), // Identity
              done_(false)
        {
        }

        /**
         * Run the full transformation until done
         */
        void run(bool verbose = false)
        {
            while (!done_)
            {
                if (done_)
                    return;

                _step_lemma(current_pair_, opG_l_, opG_r_, done_, verbose);
            }
        }

        /**
         * Get the final result
         */
        UprightResult get_result()
        {
            GridOp opG = opG_l_ * opG_r_;
            original_pair_.apply_grid_op(opG);

            Ellipse ellipseA_upright = original_pair_.e1();
            Ellipse ellipseB_upright = original_pair_.e2();

            Rectangle bboxA = ellipseA_upright.bbox();
            Rectangle bboxB = ellipseB_upright.bbox();

            return UprightResult(opG, bboxA, bboxB);
        }
    };

} // namespace gridsynth
