#pragma once

#include <vector>
#include <cmath>
#include <mpfr.h>

#include "mymath.hpp"
#include "ring.hpp"
#include "region.hpp"
#include "types.hpp"

namespace gridsynth
{

    // Object-oriented ODGP solver that reuses internal buffers to minimize
    // allocations and temporary objects across calls.
    class OdgpSolver
    {
    public:
        OdgpSolver() = default;
        ~OdgpSolver()
        {
            if (scratch_init_)
            {
                mpfr_clear(mp_a_);
                mpfr_clear(mp_b_);
                mpfr_clear(mp_tmp_);
                mpfr_clear(mp_real_);
                mpfr_clear(mp_conj_);
                mpfr_clear(mp_tmp2_);
                mpfr_clear(mp_q_);
                mpfr_clear(mp_R0_);
                mpfr_clear(mp_Rslope_);
                mpfr_clear(mp_C0_);
                mpfr_clear(mp_Cslope_);
                mpfr_clear(mp_rs_step_);
                mpfr_clear(mp_cs_step_);
                mpfr_clear(mp_I_l_);
                mpfr_clear(mp_I_r_);
                mpfr_clear(mp_J_l_);
                mpfr_clear(mp_J_r_);
            }
        }

        // profiler removed

        // Solve the Orthogonal Diophantine Grid Problem: find all beta in Z[√2]
        // with beta ∈ I and beta* ∈ J.
        std::vector<ZRootTwo> solve(const Interval &I, const Interval &J)
        {
            std::vector<ZRootTwo> out;
            if (I.width() < 0 || J.width() < 0)
                return out;

            // Initial shift (alpha) to reduce search region
            static const Float TWO(2.0);
            static const Float FOUR(4.0);
            Integer a = floor_to_integer((I.l() + J.l()) / TWO);
            Integer b = floor_to_integer(SQRT2 * (I.l() - J.l()) / FOUR);
            ZRootTwo alpha(a, b);

            Interval shifted_I = I - alpha.to_real();
            Interval shifted_J = J - alpha.conj_sq2().to_real();

            // Hoist original interval endpoints into local mpfr scratch for locality
            ensure_scratch();
            {
                GMPFloat I_l_f = I.l();
                GMPFloat I_r_f = I.r();
                GMPFloat J_l_f = J.l();
                GMPFloat J_r_f = J.r();
                mpfr_set(mp_I_l_, I_l_f.get_mpfr(), MPFR_RNDN);
                mpfr_set(mp_I_r_, I_r_f.get_mpfr(), MPFR_RNDN);
                mpfr_set(mp_J_l_, J_l_f.get_mpfr(), MPFR_RNDN);
                mpfr_set(mp_J_r_, J_r_f.get_mpfr(), MPFR_RNDN);
            }

            ZRootTwo M = ZRootTwo::from_int(1);
            bool conj_flag = false;
            solve_internal_emit(shifted_I, shifted_J, M, conj_flag, alpha, I, J, out);
            return out;
        }

        // Solve with parity constraint
        std::vector<ZRootTwo> solve_with_parity(const Interval &I, const Interval &J, const ZRootTwo &beta)
        {
            Integer p = beta.parity();
            static const Float TWO(2.0);
            Interval scaled_I = (I + (-static_cast<Float>(p))) * (SQRT2 / TWO);
            Interval scaled_J = (J + (-static_cast<Float>(p))) * (-SQRT2 / TWO);

            tmp_z_.clear();
            auto base = solve(scaled_I, scaled_J);

            std::vector<ZRootTwo> out;
            out.reserve(base.size());
            for (const auto &alpha : base)
                out.push_back(alpha * ZRootTwo(0, 1) + ZRootTwo::from_int(p));
            return out;
        }

        // Solve scaled ODGP: returns DRootTwo with denominator exponent k
        std::vector<DRootTwo> solve_scaled(const Interval &I, const Interval &J, Integer k)
        {
            Float scale = pow_sqrt2(k);
            Interval scaled_I = I * scale;
            Interval scaled_J(Float(0.0), Float(0.0));
            if (k & 1)
                scaled_J = J * (-scale);
            else
                scaled_J = J * scale;

            auto sol = solve(scaled_I, scaled_J);
            std::vector<DRootTwo> out;
            out.reserve(sol.size());
            for (const auto &alpha : sol)
                out.emplace_back(alpha, k);
            return out;
        }

        // Scaled ODGP with parity
        std::vector<DRootTwo> solve_scaled_with_parity(const Interval &I, const Interval &J, Integer k, const DRootTwo &beta)
        {
            if (k == 0)
            {
                ZRootTwo beta_z = beta.renew_denomexp(0).alpha();
                auto sol = solve_with_parity(I, J, beta_z);
                std::vector<DRootTwo> out;
                out.reserve(sol.size());
                for (const auto &a : sol)
                    out.push_back(DRootTwo::from_zroottwo(a));
                return out;
            }

            Integer p = beta.renew_denomexp(k).parity();
            DRootTwo offset = (p == 0) ? DRootTwo::from_int(0) : DRootTwo::power_of_inv_sqrt2(k);
            Interval shifted_I = I + (-offset.to_real());
            Interval shifted_J = J + (-offset.conj_sq2().to_real());

            auto sol = solve_scaled(shifted_I, shifted_J, k - 1);
            std::vector<DRootTwo> out;
            out.reserve(sol.size());
            for (const auto &a : sol)
                out.push_back(a + offset);
            return out;
        }

    private:
        void ensure_scratch()
        {
            if (!scratch_init_)
            {
                mpfr_init2(mp_a_, GMPFloat::get_default_precision());
                mpfr_init2(mp_b_, GMPFloat::get_default_precision());
                mpfr_init2(mp_tmp_, GMPFloat::get_default_precision());
                mpfr_init2(mp_real_, GMPFloat::get_default_precision());
                mpfr_init2(mp_conj_, GMPFloat::get_default_precision());
                mpfr_init2(mp_tmp2_, GMPFloat::get_default_precision());
                mpfr_init2(mp_q_, GMPFloat::get_default_precision());
                mpfr_init2(mp_R0_, GMPFloat::get_default_precision());
                mpfr_init2(mp_Rslope_, GMPFloat::get_default_precision());
                mpfr_init2(mp_C0_, GMPFloat::get_default_precision());
                mpfr_init2(mp_Cslope_, GMPFloat::get_default_precision());
                mpfr_init2(mp_rs_step_, GMPFloat::get_default_precision());
                mpfr_init2(mp_cs_step_, GMPFloat::get_default_precision());
                mpfr_init2(mp_I_l_, GMPFloat::get_default_precision());
                mpfr_init2(mp_I_r_, GMPFloat::get_default_precision());
                mpfr_init2(mp_J_l_, GMPFloat::get_default_precision());
                mpfr_init2(mp_J_r_, GMPFloat::get_default_precision());
                scratch_init_ = true;
            }
        }
        // Recursive internal solver. Appends solutions to tmp_z_ after applying
        // accumulated transform (conjugation flag and right-multiplication by M).
        void solve_internal(const Interval &I, const Interval &J, const ZRootTwo &M, bool conj_flag)
        {
                        if (I.width() < 0 || J.width() < 0)
                return;

            if (I.width() > 0 && J.width() <= 0)
            {
                // Swap intervals and toggle conjugation of the base solutions
                solve_internal(J, I, M, !conj_flag);
                return;
            }

            Integer n = 0;
            if (J.width() > 0)
            {
                static const Float LAMBDA_REAL = LAMBDA.to_real();
                auto [n_val, _] = floorlog(J.width(), LAMBDA_REAL);
                n = n_val;
            }

            if (n == 0)
            {
                static const Float TWO(2.0);
                Integer a_min = ceil_to_integer((I.l() + J.l()) / TWO);
                Integer a_max = floor_to_integer((I.r() + J.r()) / TWO);

                for (Integer a = a_min; a <= a_max; ++a)
                {
                    Integer b_min = ceil_to_integer(SQRT2 * (a - J.r()) / TWO);
                    Integer b_max = floor_to_integer(SQRT2 * (a - J.l()) / TWO);
                    for (Integer b = b_min; b <= b_max; ++b)
                    {
                        ZRootTwo beta(a, b);
                        if (conj_flag)
                            beta = beta.conj_sq2();
                        tmp_z_.push_back(beta * M);
                    }
                }
            }
            else
            {
                const auto &lp = get_lambda_powers(n);
                const ZRootTwo &lambda_inv_n = lp.lambda_inv_n;

                Interval scaled_I = I * lp.lambda_n_real;
                Interval scaled_J = J * lp.lambda_conj_n_real;

                ZRootTwo newM = M * lambda_inv_n;
                solve_internal(scaled_I, scaled_J, newM, conj_flag);
            }
        }

        // Recursive internal solver that emits fully checked candidates directly into 'out'.
        void solve_internal_emit(const Interval &I, const Interval &J,
                                 const ZRootTwo &M, bool conj_flag,
                                 const ZRootTwo &alpha,
                                 const Interval &orig_I, const Interval &orig_J,
                                 std::vector<ZRootTwo> &out)
        {
                        if (I.width() < 0 || J.width() < 0)
                return;

            if (I.width() > 0 && J.width() <= 0)
            {
                solve_internal_emit(J, I, M, !conj_flag, alpha, orig_I, orig_J, out);
                return;
            }

            Integer n = 0;
            if (J.width() > 0)
            {
                static const Float LAMBDA_REAL = LAMBDA.to_real();
                auto [n_val, _] = floorlog(J.width(), LAMBDA_REAL);
                n = n_val;
            }

            if (n == 0)
            {
                static const Float TWO("2.0");
                Integer a_min = ceil_to_integer((I.l() + J.l()) / TWO);
                Integer a_max = floor_to_integer((I.r() + J.r()) / TWO);

                // Use hoisted mpfr endpoints cached in solve()
                const mpfr_t &I_l = mp_I_l_;
                const mpfr_t &I_r = mp_I_r_;
                const mpfr_t &J_l = mp_J_l_;
                const mpfr_t &J_r = mp_J_r_;

                // Ensure scratch mpfr variables are initialized once
                ensure_scratch();

                // Cache integer coefficients of M and alpha
                const Integer Ma = M.a();
                const Integer Mb = M.b();
                const Integer Alpa = alpha.a();
                const Integer Alpb = alpha.b();

                const Integer two_Mb = Mb * Integer(2);
                const Integer step = conj_flag ? Integer(-1) : Integer(1);
                // Precompute slopes independent of 'a'
                mpfr_set_z(mp_Rslope_, two_Mb.get_mpz_t(), MPFR_RNDN);
                mpfr_mul_z(mp_tmp_, SQRT2.get_mpfr(), Ma.get_mpz_t(), MPFR_RNDN);
                mpfr_add(mp_Rslope_, mp_Rslope_, mp_tmp_, MPFR_RNDN);
                mpfr_set_z(mp_Cslope_, two_Mb.get_mpz_t(), MPFR_RNDN);
                mpfr_mul_z(mp_tmp_, SQRT2.get_mpfr(), Ma.get_mpz_t(), MPFR_RNDN);
                mpfr_sub(mp_Cslope_, mp_Cslope_, mp_tmp_, MPFR_RNDN);
                if (step > 0)
                {
                    mpfr_set(mp_rs_step_, mp_Rslope_, MPFR_RNDN);
                    mpfr_set(mp_cs_step_, mp_Cslope_, MPFR_RNDN);
                }
                else
                {
                    mpfr_neg(mp_rs_step_, mp_Rslope_, MPFR_RNDN);
                    mpfr_neg(mp_cs_step_, mp_Cslope_, MPFR_RNDN);
                }

                for (Integer a = a_min; a <= a_max; ++a)
                {
                    Integer b_min = ceil_to_integer(SQRT2 * (a - J.r()) / TWO);
                    Integer b_max = floor_to_integer(SQRT2 * (a - J.l()) / TWO);
                    if (b_max < b_min) continue;

                    const Integer base_a = a * Ma;
                    const Integer base_b = a * Mb;

                    // Refine b-range using linear constraints with bb
                    // R(bb) = (base_a+Alpa) + (base_b+Alpb)*sqrt2 + bb*(2*Mb + Ma*sqrt2)
                    // C(bb) = (base_a+Alpa) - (base_b+Alpb)*sqrt2 + bb*(2*Mb - Ma*sqrt2)
                    ensure_scratch();
                    GMPInteger BaA = base_a + Alpa;
                    GMPInteger BbB = base_b + Alpb;

                    // Build offsets
                    mpfr_set_z(mp_R0_, BaA.get_mpz_t(), MPFR_RNDN);
                    mpfr_mul_z(mp_tmp_, SQRT2.get_mpfr(), BbB.get_mpz_t(), MPFR_RNDN);
                    mpfr_add(mp_R0_, mp_R0_, mp_tmp_, MPFR_RNDN);

                    mpfr_set_z(mp_C0_, BaA.get_mpz_t(), MPFR_RNDN);
                    mpfr_mul_z(mp_tmp_, SQRT2.get_mpfr(), BbB.get_mpz_t(), MPFR_RNDN);
                    mpfr_sub(mp_C0_, mp_C0_, mp_tmp_, MPFR_RNDN);

                    // slopes already computed outside loop

                    auto intersect_bounds = [&](const mpfr_t &L, const mpfr_t &R, const mpfr_t &slope, const mpfr_t &offset, GMPInteger &cur_min, GMPInteger &cur_max) {
                        // |slope| <= tol => offset must be inside [L,R]
                        mpfr_abs(mp_tmp2_, slope, MPFR_RNDN);
                        if (mpfr_cmp_d(mp_tmp2_, 1e-40) <= 0) {
                            if (!(mpfr_cmp(offset, L) >= 0 && mpfr_cmp(offset, R) <= 0)) {
                                cur_max = cur_min - Integer(1); // empty
                            }
                            return;
                        }
                        // Compute q_low and q_high
                        // ql = ceil((L-offset)/slope), qh = floor((R-offset)/slope)
                        mpfr_sub(mp_q_, L, offset, MPFR_RNDN);
                        mpfr_div(mp_q_, mp_q_, slope, MPFR_RNDN);
                        GMPInteger ql;
                        mpfr_get_z(ql.get_mpz_t(), mp_q_, MPFR_RNDU);

                        mpfr_sub(mp_q_, R, offset, MPFR_RNDN);
                        mpfr_div(mp_q_, mp_q_, slope, MPFR_RNDN);
                        GMPInteger qh;
                        mpfr_get_z(qh.get_mpz_t(), mp_q_, MPFR_RNDD);

                        if (mpfr_sgn(slope) < 0) {
                            GMPInteger t = ql; ql = qh; qh = t;
                        }
                        if (ql > cur_min) cur_min = ql;
                        if (qh < cur_max) cur_max = qh;
                    };

                    // Work in bb space initially
                    GMPInteger bb_min = conj_flag ? -b_max : b_min;
                    GMPInteger bb_max = conj_flag ? -b_min : b_max;
                    intersect_bounds(I_l, I_r, mp_Rslope_, mp_R0_, bb_min, bb_max);
                    if (bb_max < bb_min) continue;
                    intersect_bounds(J_l, J_r, mp_Cslope_, mp_C0_, bb_min, bb_max);
                    if (bb_max < bb_min) continue;

                    // Map back to b range and clamp to original
                    if (conj_flag) {
                        // b = -bb
                        GMPInteger nb_min = -bb_max;
                        GMPInteger nb_max = -bb_min;
                        bb_min = nb_min; bb_max = nb_max;
                    }
                    if (bb_min > b_min) b_min = bb_min;
                    if (bb_max < b_max) b_max = bb_max;
                    if (b_max < b_min) continue;
                    Integer b = b_min;
                    Integer bb = conj_flag ? -b : b;
                    Integer prod_a = base_a + two_Mb * bb;
                    Integer prod_b = base_b + Ma * bb;

                    // Initialize exact MPFR linear forms at bb0
                    // real = R0 + bb*Rslope; conj = C0 + bb*Cslope
                    mpfr_mul_z(mp_tmp_, mp_Rslope_, bb.get_mpz_t(), MPFR_RNDN);
                    mpfr_add(mp_real_, mp_R0_, mp_tmp_, MPFR_RNDN);
                    mpfr_mul_z(mp_tmp_, mp_Cslope_, bb.get_mpz_t(), MPFR_RNDN);
                    mpfr_add(mp_conj_, mp_C0_, mp_tmp_, MPFR_RNDN);

                    // Prepare per-step slope (handles conj via step sign)
                    mpfr_t rs_step, cs_step;
                    mpfr_init2(rs_step, GMPFloat::get_default_precision());
                    mpfr_init2(cs_step, GMPFloat::get_default_precision());
                    if (step > 0) {
                        mpfr_set(rs_step, mp_Rslope_, MPFR_RNDN);
                        mpfr_set(cs_step, mp_Cslope_, MPFR_RNDN);
                    } else {
                        mpfr_neg(rs_step, mp_Rslope_, MPFR_RNDN);
                        mpfr_neg(cs_step, mp_Cslope_, MPFR_RNDN);
                    }

                    for (; b <= b_max; ++b)
                    {
                        // candidate coefficients with alpha
                        Integer cand_a = prod_a + Alpa;
                        Integer cand_b = prod_b + Alpb;

                        // mp_real_/mp_conj_ are currently at bb; just test
                        if (mpfr_cmp(mp_real_, I_l) >= 0 && mpfr_cmp(mp_real_, I_r) <= 0 &&
                            mpfr_cmp(mp_conj_, J_l) >= 0 && mpfr_cmp(mp_conj_, J_r) <= 0)
                        {
                            out.emplace_back(cand_a, cand_b);
                        }

                        // advance to next b
                        bb += step;
                        prod_a += two_Mb * step;
                        prod_b += Ma * step;
                        mpfr_add(mp_real_, mp_real_, rs_step, MPFR_RNDN);
                        mpfr_add(mp_conj_, mp_conj_, cs_step, MPFR_RNDN);
                    }

                    mpfr_clear(rs_step);
                    mpfr_clear(cs_step);
                }
            }
            else
            {
                const auto &lp = get_lambda_powers(n);
                Interval scaled_I = I * lp.lambda_n_real;
                Interval scaled_J = J * lp.lambda_conj_n_real;
                ZRootTwo newM = M * lp.lambda_inv_n;
                solve_internal_emit(scaled_I, scaled_J, newM, conj_flag, alpha, orig_I, orig_J, out);
            }
        }

        // Cache for powers of LAMBDA and its conjugates (and their real values)
        struct LambdaPowTriplet
        {
            ZRootTwo lambda_n;
            ZRootTwo lambda_conj_n;
            ZRootTwo lambda_inv_n;
            Float lambda_n_real;
            Float lambda_conj_n_real;
        };

        const LambdaPowTriplet &get_lambda_powers(Integer n)
        {
            size_t idx = static_cast<size_t>(n);
            // static thread-local cache shared across instances
            static thread_local std::vector<LambdaPowTriplet> cache;
            if (cache.empty())
            {
                LambdaPowTriplet base{ZRootTwo::from_int(1), ZRootTwo::from_int(1), ZRootTwo::from_int(1), Float(1.0), Float(1.0)};
                cache.push_back(base);
                // Pre-seed a small window to improve first-use locality
                static const ZRootTwo LAMBDA_CONJ = LAMBDA.conj_sq2();
                static const ZRootTwo LAMBDA_INV = LAMBDA.inv();
                for (int i = 0; i < 8; ++i)
                {
                    const auto &prev = cache.back();
                    LambdaPowTriplet next{prev.lambda_n * LAMBDA,
                                          prev.lambda_conj_n * LAMBDA_CONJ,
                                          prev.lambda_inv_n * LAMBDA_INV,
                                          Float(0.0), Float(0.0)};
                    next.lambda_n_real = next.lambda_n.to_real();
                    next.lambda_conj_n_real = next.lambda_conj_n.to_real();
                    cache.push_back(std::move(next));
                }
            }
            if (idx >= cache.size())
            {
                static const ZRootTwo LAMBDA_CONJ = LAMBDA.conj_sq2();
                static const ZRootTwo LAMBDA_INV = LAMBDA.inv();
                while (cache.size() <= idx)
                {
                    const auto &prev = cache.back();
                    LambdaPowTriplet next{prev.lambda_n * LAMBDA,
                                          prev.lambda_conj_n * LAMBDA_CONJ,
                                          prev.lambda_inv_n * LAMBDA_INV,
                                          Float(0.0), Float(0.0)};
                    next.lambda_n_real = next.lambda_n.to_real();
                    next.lambda_conj_n_real = next.lambda_conj_n.to_real();
                    cache.push_back(std::move(next));
                }
            }
            return cache[idx];
        }

        std::vector<ZRootTwo> tmp_z_;
        std::vector<DRootTwo> tmp_d_;

        // mpfr scratch for tight loops
        bool scratch_init_ = false;
        mpfr_t mp_a_, mp_b_, mp_tmp_, mp_real_, mp_conj_;
        mpfr_t mp_tmp2_, mp_q_, mp_R0_, mp_Rslope_, mp_C0_, mp_Cslope_;
        mpfr_t mp_rs_step_, mp_cs_step_;
        mpfr_t mp_I_l_, mp_I_r_, mp_J_l_, mp_J_r_;
    };

} // namespace gridsynth
