/*
  Copyright (C) 2024-26 Jeffrey Pullin
  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 or 3 of the License
  (at your option).

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   A copy of the GNU General Public License is available at
   http://www.r-project.org/Licenses/
*/

#ifndef LMM_SC_INT_H
#define LMM_SC_INT_H

#include <Eigen/Dense>
#include <algorithm>
#include <vector>

class LMM_SC_INT {

    private:
        const Eigen::Ref<Eigen::MatrixXd> X;
        const Eigen::Ref<Eigen::VectorXd> y;
        const Eigen::Ref<Eigen::VectorXd> x;
        const Eigen::Ref<Eigen::VectorXd> ns;
        const Eigen::Ref<Eigen::MatrixXd> XTX;
        const std::vector<Eigen::MatrixXd>& ZiTZis;
        const std::vector<Eigen::MatrixXd>& ZiTXis;
        size_t N;
        size_t c;
        size_t n;
        double tol = 1e-5;
        int max_iter = 50;

	public:

        Eigen::VectorXd cum_ns;
        Eigen::MatrixXd XTX_inv;
        std::vector<Eigen::MatrixXd> Mi_invs;
        Eigen::MatrixXd XtSigma_invX;
        Eigen::MatrixXd XtSigma_invX_inv;
        Eigen::VectorXd XtSigma_invy;

        Eigen::VectorXd XTy;
        double yTy;
        std::vector<Eigen::VectorXd> ZiTyis;
        std::vector<Eigen::VectorXd> ZiTeis;
        std::vector<Eigen::MatrixXd> ZiTSigma_invZis;
        std::vector<Eigen::MatrixXd> ZiTSigma_invXis;
        std::vector<Eigen::VectorXd> ZiTSigma_inveis;
        std::vector<Eigen::VectorXd> jis;

        Eigen::VectorXd y_out;
        Eigen::VectorXd mu_out;
        Eigen::MatrixXd XtWX_inv;
        Eigen::VectorXd Xty_res;
        Eigen::MatrixXd XtWZ;
        Eigen::VectorXd ZtDy_res;
        Eigen::MatrixXd XtWDZ;
        Eigen::VectorXd Zty_res;
        Eigen::VectorXd ZtSigma_invZ_diag;
        Eigen::VectorXd ZtDSigma_invDZ_diag;
        Eigen::VectorXd ZtDSigma_invZ_diag;
        Eigen::MatrixXd ZtSigma_invX;
        Eigen::MatrixXd ZtDSigma_invX;
        Eigen::VectorXd d_out;
        Eigen::VectorXd dw_out;
        Eigen::VectorXd dwd_out;

        double sigma2;
        double tau;
        double tau_sq;
        Eigen::VectorXd theta;
        Eigen::VectorXd theta_prev;
        Eigen::VectorXd beta;
        Eigen::VectorXd beta_prev;
        Eigen::MatrixXd G;

        bool lmm_converged;
        int iter;
        double step_size;

        void compute_cum_ns() {
            cum_ns = Eigen::VectorXd::Zero(n);
            cum_ns(0) = 0;
            for (size_t i = 1; i < n; ++i) {
                cum_ns(i) = cum_ns(i - 1) + ns(i - 1);
            }
        }

        void update_G() {
            G(0, 0) = theta(1);
            G(0, 1) = theta(2);
            G(1, 0) = theta(2);
            G(1, 1) = theta(3);
        }

        void compute_sufficient_stats() {
            XTy = X.transpose() * y;
            yTy = y.squaredNorm();

            ZiTyis.resize(n);
            for (size_t i = 0; i < n; ++i) {
                const auto x_i = x.segment(cum_ns(i), ns(i));
                const auto y_i = y.segment(cum_ns(i), ns(i));

                Eigen::VectorXd ZiTyi(2);
                ZiTyi(0) = y_i.sum();
                ZiTyi(1) = x_i.dot(y_i);
                ZiTyis[i] = ZiTyi;
            }
        }

        void update_ZiTeis() {
            ZiTeis.resize(n);
            for (size_t i = 0; i < n; ++i) {
                ZiTeis[i] = ZiTyis[i] - ZiTXis[i] * beta;
            }
        }

        void update_ZiTSigma_inveis() {
            ZiTSigma_inveis.resize(n);
            for (size_t i = 0; i < n; ++i) {
                ZiTSigma_inveis[i] = tau * ZiTeis[i] - tau_sq * (ZiTZis[i] * Mi_invs[i] * ZiTeis[i]);
            }
        }

        void update_jis() {
            jis.resize(n);
            for (size_t i = 0; i < n; ++i) {
                jis[i] = tau * ZiTSigma_inveis[i] - tau_sq * (ZiTZis[i] * Mi_invs[i] * ZiTSigma_inveis[i]);
            }
        }

        void update_Mi_invs() {
            Mi_invs.clear();
            Eigen::MatrixXd G_inv = G.inverse();
            for (size_t i = 0; i < n; ++i) {
                Mi_invs.push_back((G_inv + tau * ZiTZis[i]).inverse());
            }
        }

        void update_ZiTSigma_invZis() {
            ZiTSigma_invZis.resize(n);
            for (size_t i = 0; i < n; ++i) {
                ZiTSigma_invZis[i] = tau * ZiTZis[i]
                    - tau_sq * (ZiTZis[i] * Mi_invs[i] * ZiTZis[i]);
            }
        }

        void update_ZiTSigma_invXis() {
            ZiTSigma_invXis.resize(n);
            for (size_t i = 0; i < n; ++i) {
                ZiTSigma_invXis[i] = tau * ZiTXis[i]
                    - tau_sq * (ZiTZis[i] * Mi_invs[i] * ZiTXis[i]);
            }
        }

        void update_XtSigma_invX_inv() {
            XtSigma_invX = tau * XTX;
            for (size_t i = 0; i < n; ++i) {
                XtSigma_invX.noalias() -=
                    tau_sq * (ZiTXis[i].transpose() * Mi_invs[i] * ZiTXis[i]);
            }
            XtSigma_invX_inv = XtSigma_invX.inverse();
        }

        void update_XtSigma_invy() {
            XtSigma_invy = tau * XTy;
            for (size_t i = 0; i < n; ++i) {
                XtSigma_invy.noalias() -=
                    tau_sq * (ZiTXis[i].transpose() * Mi_invs[i] * ZiTyis[i]);
            }
        }

        double compute_tr_Sigma_inv() {
            double tr_MZiTZi = 0;
            for (size_t i = 0; i < n; ++i) {
                tr_MZiTZi += (Mi_invs[i] * ZiTZis[i]).trace();
            }
            return tau * N - tau_sq * tr_MZiTZi;
        }

        void update_beta() {
            beta_prev = beta;
            beta = XtSigma_invX_inv * XtSigma_invy;
        }

        void check_psd_boundary() {
            const double eps = 1e-8;
            bool at_boundary = theta.hasNaN() || theta.size() < 4;
            if (!at_boundary) {
                const double t0 = theta(1);
                const double t1 = theta(2);
                const double t2 = theta(3);
                if (t0 <= eps || t2 <= eps) {
                    at_boundary = true;
                } else {
                    const double det = t0 * t2 - t1 * t1;
                    at_boundary = (det <= eps * t0 * t2);
                }
            }
            if (at_boundary) {
                lmm_converged = false;
            }
        }

        void update_theta() {
            theta_prev = theta;

            Eigen::VectorXd U = Eigen::VectorXd::Zero(4);
            Eigen::MatrixXd AI = Eigen::MatrixXd::Zero(4, 4);

            double eTe = yTy - 2.0 * beta.dot(XTy) + beta.dot(XTX * beta);
            const double tau3 = tau_sq * tau;
            const double tau4 = tau_sq * tau_sq;
            double rTr = tau_sq * eTe;
            for (size_t i = 0; i < n; ++i) {
                const Eigen::VectorXd M_Zite = Mi_invs[i] * ZiTeis[i];
                rTr -= 2.0 * tau3 * ZiTeis[i].dot(M_Zite);
                rTr += tau4 * M_Zite.dot(ZiTZis[i] * M_Zite);
            }
            double tr_Sigma_inv = compute_tr_Sigma_inv();

            Eigen::MatrixXd XtSigma_inv2X = tau_sq * XTX;
            for (size_t i = 0; i < n; ++i) {
                const Eigen::MatrixXd M_ZiTX = Mi_invs[i] * ZiTXis[i];
                XtSigma_inv2X.noalias() -= 2.0 * tau3 * (ZiTXis[i].transpose() * M_ZiTX);
                XtSigma_inv2X.noalias() += tau4 * (M_ZiTX.transpose() * ZiTZis[i] * M_ZiTX);
            }

            Eigen::VectorXd u_e = Eigen::VectorXd::Zero(c);
            double sum_hMh = 0;
            for (size_t i = 0; i < n; ++i) {
                const Eigen::VectorXd Mh = Mi_invs[i] * ZiTSigma_inveis[i];
                u_e.noalias() -= tau_sq * (ZiTXis[i].transpose() * Mh);
                sum_hMh += ZiTSigma_inveis[i].dot(Mh);
            }

            double tr_P = tr_Sigma_inv - (XtSigma_inv2X * XtSigma_invX_inv).trace();

            U(0) = 0.5 * (rTr - tr_P);
            AI(0, 0) = 0.5 * (
                tau * rTr - tau_sq * sum_hMh - u_e.dot(XtSigma_invX_inv * u_e)
            );

            Eigen::Matrix2d re_H[3];
            re_H[0] << 1.0, 0.0, 0.0, 0.0;
            re_H[1] << 0.0, 1.0, 1.0, 0.0;
            re_H[2] << 0.0, 0.0, 0.0, 1.0;

            Eigen::VectorXd re_quad_terms = Eigen::VectorXd::Zero(3);
            Eigen::VectorXd re_traces = Eigen::VectorXd::Zero(3);
            Eigen::VectorXd uks[3];
            for (int k = 0; k < 3; ++k) {
                Eigen::MatrixXd sum_H = Eigen::MatrixXd::Zero(c, c);
                uks[k] = Eigen::VectorXd::Zero(c);
                double ai_term = 0;
                double ai_e = 0;
                double ai_kl[3] = {0, 0, 0};
                for (size_t i = 0; i < n; ++i) {
                    const Eigen::VectorXd Hh = re_H[k] * ZiTSigma_inveis[i];
                    re_quad_terms(k) += ZiTSigma_inveis[i].dot(Hh);
                    re_traces(k) += (re_H[k] * ZiTSigma_invZis[i]).trace();
                    sum_H.noalias() +=
                        ZiTSigma_invXis[i].transpose() * re_H[k] * ZiTSigma_invXis[i];
                    uks[k].noalias() += ZiTSigma_invXis[i].transpose() * Hh;
                    ai_term += Hh.dot(ZiTSigma_invZis[i] * Hh);
                    ai_e += Hh.dot(jis[i]);
                    for (int l = 0; l < k; ++l) {
                        ai_kl[l] += Hh.dot(ZiTSigma_invZis[i] * (re_H[l] * ZiTSigma_inveis[i]));
                    }
                }
                re_traces(k) -= (XtSigma_invX_inv * sum_H).trace();
                U(k + 1) = 0.5 * (re_quad_terms(k) - re_traces(k));
                AI(k + 1, k + 1) = 0.5 * (
                    ai_term - uks[k].dot(XtSigma_invX_inv * uks[k])
                );
                AI(0, k + 1) = 0.5 * (
                    ai_e - u_e.dot(XtSigma_invX_inv * uks[k])
                );
                AI(k + 1, 0) = AI(0, k + 1);
                for (int l = 0; l < k; ++l) {
                    AI(k + 1, l + 1) = 0.5 * (
                        ai_kl[l] - uks[k].dot(XtSigma_invX_inv * uks[l])
                    );
                    AI(l + 1, k + 1) = AI(k + 1, l + 1);
                }
            }

            Eigen::VectorXd step = AI.inverse() * U;

            double ss = step_size;
            theta = theta_prev + ss * step;
            while (
                (theta(0) < 0.0) ||
                (theta(1) < 0.0) ||
                (theta(3) < 0.0) ||
                (theta(1) * theta(3) - theta(2) * theta(2) < 0.0)
            ) {
                ss *= 0.5;
                theta = theta_prev + ss * step;
                if (ss < 1e-10) {
                    theta = theta_prev;
                    break;
                }
            }

            sigma2 = theta(0);
            tau = 1.0 / sigma2;
            tau_sq = tau * tau;
        }

        void update_step_size() {
            if ((iter + 1) % 10 == 0) {
                step_size = 0.9 * step_size;
            }
        }

        void check_converge() {
            Eigen::VectorXd tol_vec_beta = Eigen::VectorXd::Constant(beta.size(), tol);
            double diff1 = ((beta - beta_prev).cwiseAbs().cwiseQuotient(
                (beta).cwiseAbs() + (beta_prev).cwiseAbs() + tol_vec_beta)).maxCoeff();
            Eigen::VectorXd tol_vec_theta = Eigen::VectorXd::Constant(theta.size(), tol);
            double diff2 = ((theta - theta_prev).cwiseAbs().cwiseQuotient(
                (theta).cwiseAbs() + (theta_prev).cwiseAbs() + tol_vec_theta)).maxCoeff();
            lmm_converged = (2 * std::max(diff1, diff2)) < tol;
        }

        void compute_output() {
            y_out = Eigen::VectorXd::Zero(n);
            ZtDy_res = Eigen::VectorXd::Zero(n);
            ZtSigma_invZ_diag = Eigen::VectorXd::Zero(n);
            ZtDSigma_invDZ_diag = Eigen::VectorXd::Zero(n);
            ZtDSigma_invZ_diag = Eigen::VectorXd::Zero(n);
            d_out = Eigen::VectorXd::Zero(n);
            dwd_out = Eigen::VectorXd::Zero(n);
            ZtSigma_invX = Eigen::MatrixXd::Zero(n, c);
            ZtDSigma_invX = Eigen::MatrixXd::Zero(n, c);
            XtWZ = Eigen::MatrixXd::Zero(c, n);
            XtWDZ = Eigen::MatrixXd::Zero(c, n);

            for (size_t i = 0; i < n; ++i) {
                y_out(i) = ZiTSigma_inveis[i](0);
                ZtDy_res(i) = ZiTSigma_inveis[i](1);
                ZtSigma_invZ_diag(i) = ZiTSigma_invZis[i](0, 0);
                ZtDSigma_invZ_diag(i) = ZiTSigma_invZis[i](0, 1);
                ZtDSigma_invDZ_diag(i) = ZiTSigma_invZis[i](1, 1);
                d_out(i) = ZiTZis[i](0, 1);
                dwd_out(i) = ZiTZis[i](1, 1);
                ZtSigma_invX.row(i) = ZiTSigma_invXis[i].row(0);
                ZtDSigma_invX.row(i) = ZiTSigma_invXis[i].row(1);
            }

            mu_out = ns;
            Zty_res = y_out;
            dw_out = d_out;
            // P-score correction: t = X^T Σ^{-1} Z g, using (X^T Σ^{-1} X)^{-1}.
            XtWZ = ZtSigma_invX.transpose();
            XtWDZ = ZtDSigma_invX.transpose();
            XtWX_inv = XtSigma_invX_inv;
            Xty_res = XtSigma_invy - XtSigma_invX * beta;
        }

        LMM_SC_INT(
            const Eigen::Ref<Eigen::MatrixXd> X_,
            const Eigen::Ref<Eigen::VectorXd> y_,
            const Eigen::Ref<Eigen::VectorXd> x_,
            const Eigen::Ref<Eigen::VectorXd> ns_,
            const Eigen::Ref<Eigen::MatrixXd> XTX_,
            const std::vector<Eigen::MatrixXd>& ZiTZis_,
            const std::vector<Eigen::MatrixXd>& ZiTXis_
        ) :
            X(X_),
            y(y_),
            x(x_),
            ns(ns_),
            XTX(XTX_),
            ZiTZis(ZiTZis_),
            ZiTXis(ZiTXis_)
        {
            N = X.rows();
            c = X.cols();
            n = ns.size();
            init_params();
        };

        void init_params() {
            compute_cum_ns();
            compute_sufficient_stats();
            XTX_inv = XTX.inverse();
            beta = XTX_inv * XTy;
            update_ZiTeis();
            sigma2 = 1;
            tau = 1.0 / sigma2;
            tau_sq = tau * tau;
            G = Eigen::MatrixXd::Zero(2, 2);
            theta = Eigen::VectorXd::Zero(4);
            theta(0) = sigma2;
            theta(1) = 1.0;
            theta(2) = 0.5;
            theta(3) = 1.0;
            update_G();
            update_Mi_invs();
            update_ZiTSigma_invZis();
            update_ZiTSigma_invXis();
            update_ZiTSigma_inveis();
            update_jis();
            iter = 0;
            lmm_converged = false;
            step_size = 1;
        }

        void fit() {

            while (iter < max_iter) {

                update_Mi_invs();
                update_ZiTSigma_invZis();
                update_ZiTSigma_invXis();
                update_XtSigma_invX_inv();
                update_XtSigma_invy();
                update_beta();
                update_ZiTeis();
                update_ZiTSigma_inveis();
                update_jis();
                update_theta();
                update_G();
                check_converge();
                if (lmm_converged) {
                    break;
                }
                update_step_size();
                iter += 1;
            }

            if (theta.hasNaN() || (beta.hasNaN())) {
                lmm_converged = false;
            }
            check_psd_boundary();

            update_Mi_invs();
            update_ZiTSigma_invZis();
            update_ZiTSigma_invXis();
            update_XtSigma_invX_inv();
            update_XtSigma_invy();
            update_beta();
            update_ZiTeis();
            update_ZiTSigma_inveis();
            compute_output();
        }
};

#endif
