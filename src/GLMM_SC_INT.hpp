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

#ifndef GLMM_SC_INT_H
#define GLMM_SC_INT_H

#include <Eigen/Dense>
#include <limits>
#include <vector>
#include "GLM.hpp"
#include "Family.hpp"

class GLMM_SC_INT {

    private:
        const Eigen::Ref<Eigen::MatrixXd> X;
        const Eigen::Ref<Eigen::VectorXd> y;
        const Eigen::Ref<Eigen::MatrixXd> X_int;
        const Eigen::Ref<Eigen::VectorXd> offset;
        std::unique_ptr<Family> family;
        const Eigen::Ref<Eigen::VectorXd> ns;
        size_t N;
        size_t c;
        size_t n;
        size_t K;
        size_t m;
        double tol = 1e-5;
        int max_iter = 50;

	public:

        Eigen::VectorXd cum_ns;
        Eigen::VectorXd eta;
        Eigen::VectorXd mu;
        Eigen::VectorXd u;
        Eigen::VectorXd w;
        Eigen::MatrixXd XtSigma_invX_inv;
        Eigen::MatrixXd XSigma_inv;
        Eigen::VectorXd y_tilde;
        Eigen::VectorXd Py_tilde;
        Eigen::MatrixXd Sigma_invX;
        std::vector<Eigen::MatrixXd> Xis;
        std::vector<Eigen::MatrixXd> Zis;
        std::vector<Eigen::VectorXd> wis;
        std::vector<Eigen::MatrixXd> xis;
        std::vector<Eigen::MatrixXd> M_invs;
        Eigen::MatrixXd ZtSigma_invX;
        std::vector<Eigen::MatrixXd> ZtAkSigma_invX;

        // Output.
        Eigen::VectorXd y_out;
        Eigen::VectorXd mu_out;
        Eigen::MatrixXd XtWX_inv;
        Eigen::VectorXd Xty_res;
        Eigen::MatrixXd XtWZ;
        Eigen::VectorXd ZtSigma_invZ_diag;
        Eigen::MatrixXd ZtASigma_invAZ;
        Eigen::MatrixXd ZtAy_res;
        std::vector<Eigen::MatrixXd> XtWAkZ;

        // Parameters
        Eigen::VectorXd beta;
        Eigen::VectorXd beta_prev;
        Eigen::MatrixXd G;
        Eigen::VectorXd tau;
        Eigen::VectorXd tau_prev;

        // Control parameters.
        bool glmm_converged;
        int iter;
        double step_size;
        double rel_diff_beta;
        double rel_diff_tau;

        void update_y_tilde() {
            Eigen::VectorXd mu_eta_vec = family->mu_eta(mu);
            y_tilde = (eta - offset) + (y - mu).cwiseProduct(mu_eta_vec);
        }

        void update_eta() {
            eta = X * beta + offset + u;
        }

        void update_mu() {
            mu = family->invlink(eta);
        }

        void update_w() {
            Eigen::VectorXd v_vec = family->var(mu);
            Eigen::VectorXd mu_eta_vec = family->mu_eta(mu);
            w = (v_vec.array() * mu_eta_vec.array().square()).inverse().matrix();
        }

        void compute_cum_ns() {
            cum_ns = Eigen::VectorXd::Zero(n);
            cum_ns(0) = 0;
            for (size_t i = 1; i < n; ++i) {
                cum_ns(i) = cum_ns(i - 1) + ns(i - 1);
            }
        }

        void update_G() {
            G.setZero();
            G(0, 0) = tau(0);
            for (size_t k = 1; k < m; ++k) {
                G(k, k) = tau(1);
                G(0, k) = tau(2);
                G(k, 0) = tau(2);
            }
        }

        bool tau_feasible() const {
            if (tau(0) < 0.0 || tau(1) < 0.0) {
                return false;
            }
            const double k = static_cast<double>(K);
            return (tau(0) * tau(1) - k * tau(2) * tau(2)) >= 0.0;
        }

        void update_M_invs() {
            Eigen::MatrixXd G_inv = G.inverse();
            for (size_t i = 0; i < n; ++i) {
                Eigen::MatrixXd WZ = Zis[i];
                WZ.array().colwise() *= wis[i].array();
                Eigen::MatrixXd ZiTWiZi = Zis[i].transpose() * WZ;
                M_invs[i] = (G_inv + ZiTWiZi).inverse();
            }
        }

        Eigen::VectorXd Sigma_inv_x(Eigen::VectorXd y) {
            Eigen::VectorXd out = Eigen::VectorXd::Zero(N);
            for (size_t i = 0; i < n; ++i) {
                const auto y_i = y.segment(cum_ns(i), ns(i));
                const Eigen::VectorXd& w_i = wis[i];

                Eigen::VectorXd ZtWy = Zis[i].transpose() * w_i.cwiseProduct(y_i);
                Eigen::VectorXd v = M_invs[i] * ZtWy;

                out.segment(cum_ns(i), ns(i)).array() =
                    w_i.array() * (y_i.array() - (Zis[i] * v).array());
            }
            return out;
        }

        void compute_Zis() {
            Zis.clear();
            Zis.reserve(n);
            for (size_t i = 0; i < n; ++i) {
                Eigen::MatrixXd Zi(static_cast<Eigen::Index>(ns(i)), static_cast<Eigen::Index>(m));
                Zi.col(0) = Eigen::VectorXd::Ones(ns(i));
                if (K > 0) {
                    Zi.rightCols(K) = xis[i];
                }
                Zis.push_back(Zi);
            }
        }

        void compute_Xis() {
            Xis.clear();
            Xis.reserve(n);
            for (size_t i = 0; i < n; ++i) {
                Xis.push_back(X.block(cum_ns(i), 0, ns(i), X.cols()));
            }
        }
        
        void compute_wis() {
            wis.clear();
            wis.reserve(n);
            for (size_t i = 0; i < n; ++i) {
                wis.push_back(w.segment(cum_ns(i), ns(i)));
            }
        }

        void compute_xis() {
            xis.clear();
            xis.reserve(n);
            for (size_t i = 0; i < n; ++i) {
                xis.push_back(X_int.block(cum_ns(i), 0, ns(i), K));
            }
        }

        void update_Sigma_invX() {
            for (size_t i = 0; i < n; ++i) {
                auto block = Sigma_invX.block(cum_ns(i), 0, ns(i), c);
                const Eigen::VectorXd& w_i = wis[i];

                block = Xis[i];
                block.array().colwise() *= w_i.array();

                Eigen::MatrixXd WZ = Zis[i];
                WZ.array().colwise() *= w_i.array();
                Eigen::MatrixXd A = M_invs[i] * (Zis[i].transpose() * block);
                block.noalias() -= WZ * A;
            }
        }

        void update_XtSigma_invX_inv() {
            XtSigma_invX_inv = (X.transpose() * Sigma_invX).inverse();
        }
        
        Eigen::VectorXd P_x(Eigen::VectorXd x) {
            Eigen::VectorXd a = Sigma_inv_x(x);
            Eigen::VectorXd tmp = XtSigma_invX_inv * (Sigma_invX.transpose() * x);
            return a - Sigma_invX * tmp;
        }

        Eigen::VectorXd Omega_x(Eigen::VectorXd x) {
            Eigen::VectorXd out = Eigen::VectorXd::Zero(N);
            for (size_t i = 0; i < n; ++i) {
                double sum_i = x.segment(cum_ns(i), ns(i)).sum();
                out.segment(cum_ns(i), ns(i)).setConstant(sum_i);
            }
            return out;
        }

        Eigen::VectorXd Delta_x(Eigen::VectorXd x) {
            Eigen::VectorXd out = Eigen::VectorXd::Zero(N);
            for (size_t i = 0; i < n; ++i) {
                Eigen::VectorXd x_seg = x.segment(cum_ns(i), ns(i));
                out.segment(cum_ns(i), ns(i)) = xis[i] * (xis[i].transpose() * x_seg);
            }
            return out;
        }

        Eigen::VectorXd Pi_x(Eigen::VectorXd x) {
            Eigen::VectorXd out = Eigen::VectorXd::Zero(N);
            for (size_t i = 0; i < n; ++i) {
                Eigen::VectorXd x_seg = x.segment(cum_ns(i), ns(i));
                Eigen::VectorXd xsum = xis[i].rowwise().sum();
                out.segment(cum_ns(i), ns(i)) =
                    xsum * x_seg.sum() + Eigen::VectorXd::Constant(ns(i), xsum.dot(x_seg));
            }
            return out;
        }

        Eigen::VectorXd ZtGZ_x(Eigen::VectorXd x) {
            Eigen::VectorXd out = Eigen::VectorXd::Zero(N);
            for (size_t i = 0; i < n; ++i) {
                Eigen::VectorXd Ztv = Zis[i].transpose() * x.segment(cum_ns(i), ns(i));
                out.segment(cum_ns(i), ns(i)) = Zis[i] * (G * Ztv);
            }
            return out;
        }

        Eigen::VectorXd collapse_vec(Eigen::VectorXd x) {
            Eigen::VectorXd out = Eigen::VectorXd::Zero(n);
            for (size_t i = 0; i < n; ++i) {
                out(i) = x.segment(cum_ns(i), ns(i)).sum();
            }
            return out;
        }

        void update_ZtAkSigma_invX() {
            ZtAkSigma_invX.assign(m, Eigen::MatrixXd::Zero(n, c));
            for (size_t i = 0; i < n; ++i) {
                auto block = Sigma_invX.block(cum_ns(i), 0, ns(i), c);
                ZtAkSigma_invX[0].row(i) = block.colwise().sum();
                for (size_t k = 0; k < K; ++k) {
                    ZtAkSigma_invX[k + 1].row(i) = xis[i].col(k).transpose() * block;
                }
            }
            ZtSigma_invX = ZtAkSigma_invX[0];
        }

        double compute_trPOmega() {

            double tr_Sigma_invOmega = 0;
            double tr_Sigma_invRest = 0;

            for (size_t i = 0; i < n; ++i) {
                double wi_sum = wis[i].sum();
                Eigen::VectorXd ZiTWi1 = Zis[i].transpose() * wis[i];
                tr_Sigma_invOmega += wi_sum - ZiTWi1.dot(M_invs[i] * ZiTWi1);
            }

            Eigen::MatrixXd tmp = XtSigma_invX_inv * ZtSigma_invX.transpose();
            tr_Sigma_invRest = (ZtSigma_invX * tmp).trace();

            return tr_Sigma_invOmega - tr_Sigma_invRest; 
        }

        double compute_trPDelta() {
            double tr_Sigma_invDelta = 0;
            double tr_Sigma_invRest = 0;

            for (size_t i = 0; i < n; ++i) {
                tr_Sigma_invDelta += wis[i].dot(xis[i].rowwise().squaredNorm());

                Eigen::MatrixXd WXint = xis[i];
                WXint.array().colwise() *= wis[i].array();
                Eigen::MatrixXd ZtWXint = Zis[i].transpose() * WXint;
                tr_Sigma_invDelta -= (ZtWXint.transpose() * M_invs[i] * ZtWXint).trace();
            }

            for (size_t k = 0; k < K; ++k) {
                const Eigen::MatrixXd& Dtk = ZtAkSigma_invX[k + 1];
                Eigen::MatrixXd tmp = XtSigma_invX_inv * Dtk.transpose();
                tr_Sigma_invRest += (Dtk * tmp).trace();
            }

            return tr_Sigma_invDelta - tr_Sigma_invRest; 
        }

        double compute_trPPi() {
            double tr_Sigma_invPi = 0;
            double tr_Sigma_invRest = 0;

            for (size_t i = 0; i < n; ++i) {
                Eigen::VectorXd xsum = xis[i].rowwise().sum();
                tr_Sigma_invPi += 2.0 * wis[i].dot(xsum);
                Eigen::VectorXd ZiTWixsum = Zis[i].transpose() * wis[i].cwiseProduct(xsum);
                Eigen::VectorXd ZiTWi1 = Zis[i].transpose() * wis[i];
                tr_Sigma_invPi -= 2.0 * ZiTWixsum.dot(M_invs[i] * ZiTWi1);
            }

            Eigen::MatrixXd DtSigma_invX = Eigen::MatrixXd::Zero(static_cast<Eigen::Index>(n), static_cast<Eigen::Index>(c));
            for (size_t k = 0; k < K; ++k) {
                DtSigma_invX += ZtAkSigma_invX[k + 1];
            }
            Eigen::MatrixXd tmp = XtSigma_invX_inv * DtSigma_invX.transpose();
            tr_Sigma_invRest = 2.0 * (ZtSigma_invX * tmp).trace();

            return tr_Sigma_invPi - tr_Sigma_invRest;
        }

        void update_beta() {
            beta_prev = beta;
            beta = XtSigma_invX_inv * (Sigma_invX.transpose() * y_tilde);
        }

        void update_Py_tilde() {
            Py_tilde = Sigma_inv_x(y_tilde);
            Py_tilde.noalias() -= Sigma_invX * beta;
        }

        void update_u() {
            u = ZtGZ_x(Py_tilde);
        }

        void update_tau() {
            tau_prev = tau;

            Eigen::VectorXd U = Eigen::VectorXd::Zero(3);
            Eigen::MatrixXd AI = Eigen::MatrixXd::Zero(3, 3);

            Eigen::VectorXd OmegaPy_tilde = Omega_x(Py_tilde);
            Eigen::VectorXd DeltaPy_tilde = Delta_x(Py_tilde);
            Eigen::VectorXd PiPy_tilde = Pi_x(Py_tilde);

            U(0) = Py_tilde.dot(OmegaPy_tilde) - compute_trPOmega();
            U(1) = Py_tilde.dot(DeltaPy_tilde) - compute_trPDelta();
            U(2) = Py_tilde.dot(PiPy_tilde) - compute_trPPi();

            Eigen::VectorXd POmegaPy_tilde = P_x(OmegaPy_tilde);
            Eigen::VectorXd PDeltaPy_tilde = P_x(DeltaPy_tilde);
            Eigen::VectorXd PPiPy_tilde = P_x(PiPy_tilde);

            AI(0, 0) = OmegaPy_tilde.dot(POmegaPy_tilde);
            AI(1, 1) = DeltaPy_tilde.dot(PDeltaPy_tilde);
            AI(2, 2) = PiPy_tilde.dot(PPiPy_tilde);
            AI(0, 1) = OmegaPy_tilde.dot(PDeltaPy_tilde);
            AI(0, 2) = OmegaPy_tilde.dot(PPiPy_tilde);
            AI(1, 2) = DeltaPy_tilde.dot(PPiPy_tilde);
            AI(1, 0) = AI(0, 1);
            AI(2, 0) = AI(0, 2);
            AI(2, 1) = AI(1, 2);

            Eigen::VectorXd step = AI.inverse() * U;
            
            double ss = step_size;
            tau = tau_prev + ss * step;
            while (!tau_feasible()) {
                ss *= 0.5;
                tau = tau_prev + ss * step;
                if (ss < 1e-10) {
                    tau = tau_prev;
                    break;
                }
            }
    
        }

        void update_step_size() {
            if ((iter + 1) % 10 == 0) {
                step_size = 0.9 * step_size;
            }
        }

        void check_converge() {
            Eigen::VectorXd tol_vec_beta = Eigen::VectorXd::Constant(beta.size(), tol);
            rel_diff_beta = ((beta - beta_prev).cwiseAbs().cwiseQuotient(
                (beta).cwiseAbs() + (beta_prev).cwiseAbs() + tol_vec_beta)).maxCoeff();
            Eigen::VectorXd tol_vec_tau = Eigen::VectorXd::Constant(tau.size(), tol);
            rel_diff_tau = ((tau - tau_prev).cwiseAbs().cwiseQuotient(
                (tau).cwiseAbs() + (tau_prev).cwiseAbs() + tol_vec_tau)).maxCoeff();
            glmm_converged = (2 * std::max(rel_diff_beta, rel_diff_tau)) < tol;
        }

        void check_psd_boundary() {
            const double eps = 1e-8;
            bool at_boundary = tau.hasNaN() || tau.size() < 3;
            if (!at_boundary) {
                at_boundary = (tau(0) <= eps) || (tau(1) <= eps) || !tau_feasible();
            }
            if (at_boundary) {
                glmm_converged = false;
            }
        }

        void compute_output() {

            mu_out = collapse_vec(mu);
            Eigen::VectorXd y_res = (y.array() - mu.array());
            y_out = collapse_vec(y_res);
    
            XtWX_inv = (X.transpose() * w.asDiagonal() * X).inverse();
            Xty_res = X.transpose() * y_res;

            XtWAkZ.assign(m, Eigen::MatrixXd::Zero(c, n));
            ZtAy_res = Eigen::MatrixXd::Zero(n, m);
            ZtASigma_invAZ = Eigen::MatrixXd::Zero(n, m * m);

            for (size_t i = 0; i < n; ++i) {
                const Eigen::Index start = static_cast<Eigen::Index>(cum_ns(i));
                const Eigen::Index len = static_cast<Eigen::Index>(ns(i));
                const Eigen::VectorXd y_res_i = y_res.segment(start, len);
                const Eigen::VectorXd& w_i = wis[i];

                XtWAkZ[0].col(i) = Xis[i].transpose() * w_i;
                ZtAy_res(i, 0) = y_res_i.sum();
                for (size_t k = 0; k < K; ++k) {
                    Eigen::VectorXd wxk = xis[i].col(k).cwiseProduct(w_i);
                    XtWAkZ[k + 1].col(i) = Xis[i].transpose() * wxk;
                    ZtAy_res(i, k + 1) = xis[i].col(k).dot(y_res_i);
                }

                Eigen::MatrixXd WZ = Zis[i];
                WZ.array().colwise() *= w_i.array();
                Eigen::MatrixXd A = Zis[i].transpose() * WZ;
                Eigen::MatrixXd ZiTSigmaInvZi = A - A * M_invs[i] * A;
                for (size_t r = 0; r < m; ++r) {
                    for (size_t col = 0; col < m; ++col) {
                        ZtASigma_invAZ(i, r * m + col) = ZiTSigmaInvZi(r, col);
                    }
                }
            }

            XtWZ = XtWAkZ[0];
            ZtSigma_invZ_diag = ZtASigma_invAZ.col(0);
        }

        GLMM_SC_INT(
            const Eigen::Ref<Eigen::MatrixXd> X_, 
            const Eigen::Ref<Eigen::VectorXd> y_, 
            const Eigen::Ref<Eigen::MatrixXd> X_int_,
            const Eigen::Ref<Eigen::VectorXd> offset_,
            std::unique_ptr<Family> family_, 
            const Eigen::Ref<Eigen::VectorXd> ns_
        ) : 
            X(X_),
            y(y_),
            X_int(X_int_),
            offset(offset_),
            family(std::move(family_)),
            ns(ns_)
        {   
            N = X.rows();
            c = X.cols();
            n = ns.size();
            K = static_cast<size_t>(X_int.cols());
            m = K + 1;
            init_params();
        };

        void init_params() {
            auto poisson = std::unique_ptr<Family>(new Poisson());
            GLM glm(X, y, offset, std::move(poisson));
            glm.fit();
            beta = glm.beta;
            mu = glm.mu;
            compute_cum_ns();
            compute_Xis();
            compute_xis();
            compute_Zis();
            u = Eigen::VectorXd::Zero(N);
            update_eta();
            update_w();
            compute_wis();
            update_y_tilde();
            G = Eigen::MatrixXd::Zero(m, m);
            tau = Eigen::VectorXd::Zero(3);
            tau(0) = 1.0;
            tau(1) = (K == 1) ? 1.0 : 1e-4;
            tau(2) = (K == 1) ? 0.5 : 0.0;
            update_G();
            M_invs.resize(n);
            update_M_invs();
            Sigma_invX.resize(N, c);
            Py_tilde.resize(N);
            iter = 0;
            glmm_converged = false;
            step_size = 1;
            rel_diff_beta = 0.0;
            rel_diff_tau = 0.0;
        }

        void fit() {

            while (iter < max_iter) {
           
                update_M_invs();
                update_Sigma_invX();
                update_XtSigma_invX_inv();
                update_ZtAkSigma_invX();
                update_beta();
                update_Py_tilde();
                update_u();
                update_eta();
                update_mu();
                update_tau();
                update_G();
                update_w();
                compute_wis();
                update_y_tilde();
                check_converge();
                if (glmm_converged) {
                    break;
                }
                update_step_size();
                iter += 1;
            }

            if (tau.hasNaN() || (beta.hasNaN())) {
                glmm_converged = false;
            }
            check_psd_boundary();

            update_M_invs();
            update_Sigma_invX();
            update_XtSigma_invX_inv();
            update_ZtAkSigma_invX();
            compute_output();
        }
};

#endif
