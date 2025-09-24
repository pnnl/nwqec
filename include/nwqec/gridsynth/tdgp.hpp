#pragma once

#include <vector>
#include <optional>
#include <iostream>
#include <algorithm>
#include <memory>
// profiler includes removed

#include "odgp.hpp"
#include "ring.hpp"
#include "region.hpp"
#include "grid_op.hpp"
#include "types.hpp"

namespace gridsynth
{

    // Object-oriented TDGP solver that reuses internal buffers and avoids
    // excessive temporary object creation.
    class TdgpSolver
    {
    public:
        TdgpSolver(const ConvexSet &setA,
                   const ConvexSet &setB,
                   const GridOp &opG_inv,
                   const Rectangle &bboxA,
                   const Rectangle &bboxB,
                   const Interval &bboxA_y_fattened,
                   const Interval &bboxB_y_fattened)
            : setA_(setA.clone()),
              setB_(setB.clone()),
              opG_inv_(opG_inv),
              bboxA_(bboxA),
              bboxB_(bboxB),
              bboxA_y_fattened_(bboxA_y_fattened),
              bboxB_y_fattened_(bboxB_y_fattened)
        {
        }

        // Solve for a given scale k and return solutions.
        std::vector<DOmega> solve(Integer k, bool verbose = false)
        {
            // Precompute invariants for this k
            const Float inv_scale_k = Float(1.0) / pow_sqrt2(k);
            const Float inv_scale_kp1 = Float(1.0) / pow_sqrt2(k + 1);
            const Float sqrt2_over_2 = SQRT2 / Float(2.0);

            // x- and y-direction ODGP solves
            dr2_x_ = odgp_.solve_scaled(bboxA_.I_x(), bboxB_.I_x(), k + 1);
            dr2_y_ = odgp_.solve_scaled(bboxA_y_fattened_, bboxB_y_fattened_, k + 1);

            suf_.clear();
            if (!dr2_x_.empty() && !dr2_y_.empty())
            {
                const DRootTwo &alpha0 = dr2_x_[0];
                const DRootTwo dx = DRootTwo::power_of_inv_sqrt2(k);

                // Precompute v and its conjugate coordinates once per k
                const DOmega v_common = opG_inv_ * DOmega::from_droottwo_vector(dx, DRootTwo::from_int(0), k);
                // Reuse local arrays to reduce temporaries
                Float v_re, v_im;
                v_common.coords_into_with(inv_scale_k, sqrt2_over_2, v_re, v_im);
                const std::array<Float, 2> v_coords = {v_re, v_im};
                Float v_conj_re, v_conj_im;
                DOmega v_conj = v_common.conj_sq2();
                v_conj.coords_into_with(inv_scale_k, sqrt2_over_2, v_conj_re, v_conj_im);
                const std::array<Float, 2> v_conj_sq2_coords = {v_conj_re, v_conj_im};

                for (const auto &beta : dr2_y_)
                {
                    DOmega z0 = opG_inv_ * DOmega::from_droottwo_vector(alpha0, beta, k + 1);

                    Float z0_re, z0_im;
                    z0.coords_into_with(inv_scale_kp1, sqrt2_over_2, z0_re, z0_im);
                    const std::array<Float, 2> z0_coords = {z0_re, z0_im};
                    auto t_A_opt = setA_->intersect(z0_coords, v_coords);

                    Float z0_conj_re, z0_conj_im;
                    DOmega z0_conj = z0.conj_sq2();
                    z0_conj.coords_into_with(inv_scale_kp1, sqrt2_over_2, z0_conj_re, z0_conj_im);
                    const std::array<Float, 2> z0_conj_sq2_coords = {z0_conj_re, z0_conj_im};
                    auto t_B_opt = setB_->intersect(z0_conj_sq2_coords, v_conj_sq2_coords);

                    if (!t_A_opt.has_value() || !t_B_opt.has_value())
                        continue;

                    auto [tA_l, tA_r] = *t_A_opt;
                    auto [tB_l, tB_r] = *t_B_opt;

                    Interval intA(tA_l, tA_r);
                    Interval intB(tB_l, tB_r);

                    DRootTwo parity = (beta - alpha0).mul_by_sqrt2_power_renewing_denomexp(k);

                    const Float intA_width = intA.width();
                    const Float intB_width = intB.width();
                    static const Float TEN("10.0");
                    const Float two_pow_k = Float((Integer(1) << k));
                    const Float dtA = TEN / max(TEN, two_pow_k * intB_width);
                    const Float dtB = TEN / max(TEN, two_pow_k * intA_width);

                    intA = intA.fatten(dtA);
                    intB = intB.fatten(dtB);

                    dr2_t_ = odgp_.solve_scaled_with_parity(intA, intB, 1, parity);
                    for (const auto &alpha : dr2_t_)
                    {
                        DRootTwo new_alpha = alpha * dx + alpha0;
                        suf_.push_back(DOmega::from_droottwo_vector(new_alpha, beta, k));
                    }
                }
            }

            // Transform and filter
            std::vector<DOmega> out;
            if (!suf_.empty())
            {
                out.reserve(suf_.size());
                for (const auto &z : suf_)
                {
                    DOmega z_tr = opG_inv_ * z;
                    Float zr, zi;
                    z_tr.coords_into_with(inv_scale_k, sqrt2_over_2, zr, zi);
                    const std::array<Float, 2> z_coords = {zr, zi};
                    Float zcr, zci;
                    DOmega z_tr_conj = z_tr.conj_sq2();
                    z_tr_conj.coords_into_with(inv_scale_k, sqrt2_over_2, zcr, zci);
                    const std::array<Float, 2> z_conj_sq2_coords = {zcr, zci};
                    if (setA_->inside(z_coords) && setB_->inside(z_conj_sq2_coords))
                        out.push_back(std::move(z_tr));
                }
            }

            if (verbose)
            {
                std::cout << "k=" << k
                          << " size of sol: " << out.size()
                          << ", size of candidates: " << suf_.size() << "\n";
            }

            return out;
        }

        // Monotone-increasing stepping from k=0 upwards.
        std::vector<DOmega> solve_next(bool verbose = false)
        {
            ++step_k_;
            return solve(step_k_, verbose);
        }

        // Quick check utility
        bool verify_solution(const DOmega &candidate) const
        {
            const Float inv_scale = Float(1.0) / pow_sqrt2(candidate.k());
            const Float sqrt2_over_2 = SQRT2 / Float(2.0);
            Float zr, zi;
            candidate.coords_into_with(inv_scale, sqrt2_over_2, zr, zi);
            const std::array<Float, 2> z_coords = {zr, zi};
            Float zcr, zci;
            DOmega cand_conj = candidate.conj_sq2();
            cand_conj.coords_into_with(inv_scale, sqrt2_over_2, zcr, zci);
            const std::array<Float, 2> z_conj_sq2_coords = {zcr, zci};
            return setA_->inside(z_coords) && setB_->inside(z_conj_sq2_coords);
        }

        // profiler removed

    private:
        std::unique_ptr<ConvexSet> setA_;
        std::unique_ptr<ConvexSet> setB_;
        GridOp opG_inv_;
        Rectangle bboxA_;
        Rectangle bboxB_;
        Interval bboxA_y_fattened_;
        Interval bboxB_y_fattened_;

        // Single ODGP solver instance
        OdgpSolver odgp_;

        // step counter for solve_next
        Integer step_k_ = -1;

        // Reused buffers
        std::vector<DRootTwo> dr2_x_;
        std::vector<DRootTwo> dr2_y_;
        std::vector<DRootTwo> dr2_t_;
        std::vector<DOmega> suf_;
        std::vector<DOmega> tmp_;

        // profiler removed
    };

} // namespace gridsynth
