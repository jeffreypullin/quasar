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
#include <iostream>
#include "GLM.hpp"
#include "Family.hpp"

class GLMM_SC_INT {

    private:
        const Eigen::Ref<Eigen::MatrixXd> X;
        const Eigen::Ref<Eigen::VectorXd> y;
        const Eigen::Ref<Eigen::VectorXd> x;
        const Eigen::Ref<Eigen::VectorXd> offset;
        std::unique_ptr<Family> family;
        const Eigen::Ref<Eigen::VectorXd> ns;
        size_t N;
        size_t c;
        size_t n;
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
        std::vector<Eigen::VectorXd> xis;
        std::vector<Eigen::Matrix2d> M_invs;
        Eigen::MatrixXd ZtSigma_invX;
        Eigen::MatrixXd ZtDSigma_invX;

        // Output.
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
        Eigen::VectorXd d_out;
        Eigen::VectorXd dw_out;
        Eigen::VectorXd dwd_out;

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

        void update_y_tilde() {
            Eigen::VectorXd mu_eta_vec = family->mu_eta(mu);
            y_tilde = (eta - offset) + (y - mu).cwiseProduct(mu_eta_vec);
        }

        void update_eta() {
            eta = X * beta + offset + u;
            for (size_t i = 0; i < N; i++) {
                if (eta(i) < -30) eta(i) = -30.0;
                if (eta(i) >  30) eta(i) =  30.0;
            }
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
            G(0, 0) = tau(0);
            G(0, 1) = tau(1);
            G(1, 0) = tau(1);
            G(1, 1) = tau(2);
        }

        void update_M_invs() {
            Eigen::Matrix2d G_inv = G.inverse();
            for (size_t i = 0; i < n; ++i) {
                Eigen::Matrix2d ZiTWiZi = Zis[i].transpose() * wis[i].asDiagonal() * Zis[i];
                M_invs[i] = (G_inv + ZiTWiZi).inverse();
            }
        }

        Eigen::VectorXd Sigma_inv_x(Eigen::VectorXd y) {
            Eigen::VectorXd out = Eigen::VectorXd::Zero(N);
            for (size_t i = 0; i < n; ++i) {
                const auto y_i = y.segment(cum_ns(i), ns(i));
                const Eigen::VectorXd& w_i = wis[i];
                const Eigen::VectorXd& x_i = xis[i];

                Eigen::Vector2d ZtWy;
                ZtWy(0) = w_i.dot(y_i);
                ZtWy(1) = (w_i.array() * x_i.array() * y_i.array()).sum();
                Eigen::Vector2d v = M_invs[i] * ZtWy;

                out.segment(cum_ns(i), ns(i)).array() =
                    w_i.array() * (y_i.array() - v(0) - v(1) * x_i.array());
            }
            return out;
        }

        void compute_Zis() {
            Zis.clear();
            for (size_t i = 0; i < n; ++i) {
                Eigen::MatrixXd Zi(static_cast<Eigen::Index>(ns(i)), 2);
                Zi.col(0) = Eigen::VectorXd::Ones(ns(i));
                Zi.col(1) = xis[i];
                Zis.push_back(Zi);
            }
        }

        void compute_Xis() {
            Xis.clear();
            for (size_t i = 0; i < n; ++i) {
                Xis.push_back(X.block(cum_ns(i), 0, ns(i), X.cols()));
            }
        }
        
        void compute_wis() {
            wis.clear();
            for (size_t i = 0; i < n; ++i) {
                wis.push_back(w.segment(cum_ns(i), ns(i)));
            }
        }

        void compute_xis() {
            xis.clear();
            for (size_t i = 0; i < n; ++i) {
                xis.push_back(x.segment(cum_ns(i), ns(i)));
            }
        }

        void update_Sigma_invX() {
            for (size_t i = 0; i < n; ++i) {
                auto block = Sigma_invX.block(cum_ns(i), 0, ns(i), c);
                const Eigen::VectorXd& w_i = wis[i];

                block = Xis[i];
                block.array().colwise() *= w_i.array();

                Eigen::MatrixXd WZ(static_cast<Eigen::Index>(ns(i)), 2);
                WZ.col(0) = w_i;
                WZ.col(1) = w_i.cwiseProduct(xis[i]);
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
                Eigen::VectorXd tmp = Eigen::VectorXd::Constant(ns(i), sum_i);
                out.segment(cum_ns(i), ns(i)) = tmp;
            }
            return out;
        }

        Eigen::VectorXd Delta_x(Eigen::VectorXd x) {
            Eigen::VectorXd out = Eigen::VectorXd::Zero(N);
            for (size_t i = 0; i < n; ++i) {
                Eigen::VectorXd x_seg = x.segment(cum_ns(i), ns(i));
                out.segment(cum_ns(i), ns(i)) = xis[i] * xis[i].dot(x_seg);
            }
            return out;
        }

        Eigen::VectorXd Pi_x(Eigen::VectorXd x) {
            Eigen::VectorXd out = Eigen::VectorXd::Zero(N);
            for (size_t i = 0; i < n; ++i) {
                Eigen::VectorXd x_seg = x.segment(cum_ns(i), ns(i));
                double x_seg_sum = x_seg.sum();
                Eigen::VectorXd tmp = Eigen::VectorXd::Constant(ns(i), xis[i].dot(x_seg));
                out.segment(cum_ns(i), ns(i)) = xis[i] * x_seg_sum + tmp;
            }
            return out;
        }

        Eigen::VectorXd ZtGZ_x(Eigen::VectorXd x) {
            Eigen::VectorXd out = Eigen::VectorXd::Zero(N);
            for (size_t i = 0; i < n; ++i) {
                const auto v_i = x.segment(cum_ns(i), ns(i));
                const double s = v_i.sum();
                const double d = xis[i].dot(v_i);
                out.segment(cum_ns(i), ns(i)).array() =
                    (G(0, 0) * s + G(0, 1) * d) +
                    (G(0, 1) * s + G(1, 1) * d) * xis[i].array();
            }
            return out;
        }

        Eigen::MatrixXd Omega_X(Eigen::MatrixXd X) {
            Eigen::MatrixXd out = Eigen::MatrixXd::Zero(X.rows(), X.cols());
            for (int j = 0; j < X.cols(); ++j) {
                out.col(j) = Omega_x(X.col(j));
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

        void update_ZtSigma_invX() {
            ZtSigma_invX = Eigen::MatrixXd::Zero(n, c);
            for (size_t i = 0; i < n; ++i) {
                ZtSigma_invX.row(i) = Sigma_invX.block(cum_ns(i), 0, ns(i), c).colwise().sum();
            }
        }

        void update_ZtDSigma_invX() {
            ZtDSigma_invX = Eigen::MatrixXd::Zero(n, c);
            for (size_t i = 0; i < n; ++i) {
                ZtDSigma_invX.row(i) =
                    xis[i].transpose() * Sigma_invX.block(cum_ns(i), 0, ns(i), c);
            }
        }

        Eigen::MatrixXd compute_DtSigma_invX() {
            Eigen::MatrixXd DtSigma_invX = Eigen::MatrixXd::Zero(n, c);
            for (size_t i = 0; i < n; ++i) {
                DtSigma_invX.row(i) = xis[i].transpose() * Sigma_invX.block(cum_ns(i), 0, ns(i), c);
            }
            return DtSigma_invX;
        }

        double compute_trPOmega() {

            double tr_Sigma_invOmega = 0;
            double tr_Sigma_invRest = 0;

            for (size_t i = 0; i < n; ++i) {
                double wi_sum = wis[i].sum();
                Eigen::VectorXd ZiTWi1 = Zis[i].transpose() * wis[i];
                tr_Sigma_invOmega += wi_sum - ZiTWi1.dot(M_invs[i] * ZiTWi1);
            }

            Eigen::MatrixXd tmp;
            tmp = (XtSigma_invX_inv) * ZtSigma_invX.transpose();
            tr_Sigma_invRest = (ZtSigma_invX * tmp).trace();

            return tr_Sigma_invOmega - tr_Sigma_invRest; 
        }

        double compute_trPDelta() {
            double tr_Sigma_invDelta = 0;
            double tr_Sigma_invRest = 0;

            for (size_t i = 0; i < n; ++i) {
                tr_Sigma_invDelta += (wis[i].array() * xis[i].array().square()).sum();

                Eigen::VectorXd ZiTWixi = Zis[i].transpose() * wis[i].cwiseProduct(xis[i]);
                tr_Sigma_invDelta -= ZiTWixi.dot(M_invs[i] * ZiTWixi);
            }

            Eigen::MatrixXd DtSigma_invX = compute_DtSigma_invX();
            Eigen::MatrixXd tmp = (XtSigma_invX_inv) * DtSigma_invX.transpose();
            tr_Sigma_invRest = (DtSigma_invX * tmp).trace();

            return tr_Sigma_invDelta - tr_Sigma_invRest; 
        }

        double compute_trPPi() {
            double tr_Sigma_invPi = 0;
            double tr_Sigma_invRest = 0;

            for (size_t i = 0; i < n; ++i) {
                double wixi_sum = wis[i].dot(xis[i]);
                Eigen::VectorXd ZiTWixi = Zis[i].transpose() * wis[i].cwiseProduct(xis[i]);
                Eigen::VectorXd ones = Eigen::VectorXd::Ones(wis[i].size());
                Eigen::VectorXd ZiTWi1 = Zis[i].transpose() * (wis[i].cwiseProduct(ones));
           
                tr_Sigma_invPi += 2 * wixi_sum - 2 * ZiTWixi.dot(M_invs[i] * ZiTWi1);
            }

            Eigen::MatrixXd DtSigma_invX = compute_DtSigma_invX();
            Eigen::MatrixXd tmp = XtSigma_invX_inv * DtSigma_invX.transpose();
            tr_Sigma_invRest = 2 * (ZtSigma_invX * tmp).trace();

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
            U(1) = Py_tilde.dot(PiPy_tilde) - compute_trPPi();
            U(2) = Py_tilde.dot(DeltaPy_tilde) - compute_trPDelta();

            Eigen::VectorXd POmegaPy_tilde = P_x(OmegaPy_tilde);
            Eigen::VectorXd PPiPy_tilde = P_x(PiPy_tilde);
            Eigen::VectorXd PDeltaPy_tilde = P_x(DeltaPy_tilde);

            AI(0, 0) = (OmegaPy_tilde.transpose() * POmegaPy_tilde);
            AI(1, 1) = (PiPy_tilde.transpose() * PPiPy_tilde);
            AI(2, 2) = (DeltaPy_tilde.transpose() * PDeltaPy_tilde);
            AI(0, 1) = (OmegaPy_tilde.transpose() * PPiPy_tilde);
            AI(0, 2) = (OmegaPy_tilde.transpose() * PDeltaPy_tilde);
            AI(1, 2) = (PiPy_tilde.transpose() * PDeltaPy_tilde);
            AI(1, 0) = AI(0, 1);
            AI(2, 0) = AI(0, 2);
            AI(2, 1) = AI(1, 2);

            Eigen::VectorXd step = AI.inverse() * U;
            
            double ss = step_size;
            tau = tau_prev + ss * step;
            while (
                (tau(0) < 0.0) ||
                (tau(2) < 0.0) ||
                (tau(0) * tau(2) - tau(1) * tau(1) < 0.0)
            ) {
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
            double diff1 = ((beta - beta_prev).cwiseAbs().cwiseQuotient(
                (beta).cwiseAbs() + (beta_prev).cwiseAbs() + tol_vec_beta)).maxCoeff();
            Eigen::VectorXd tol_vec_tau = Eigen::VectorXd::Constant(tau.size(), tol);
            double diff2 = ((tau - tau_prev).cwiseAbs().cwiseQuotient(
                (tau).cwiseAbs() + (tau_prev).cwiseAbs() + tol_vec_tau)).maxCoeff();
            glmm_converged = (2 * std::max(diff1, diff2)) < tol;
        }

        void check_psd_boundary() {
            const double eps = 1e-8;
            bool at_boundary = tau.hasNaN() || tau.size() < 3;
            if (!at_boundary) {
                const double t0 = tau(0);
                const double t1 = tau(1);
                const double t2 = tau(2);
                if (t0 <= eps || t2 <= eps) {
                    at_boundary = true;
                } else {
                    const double det = t0 * t2 - t1 * t1;
                    at_boundary = (det <= eps * t0 * t2);
                }
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
            XtWZ = Eigen::MatrixXd::Zero(c, n);
            for (size_t i = 0; i < n; ++i) {
                XtWZ.col(i) = Xis[i].transpose() * wis[i];
            }
            ZtDy_res = collapse_vec(x.cwiseProduct(y_res));
            XtWDZ = Eigen::MatrixXd::Zero(c, n);
            for (size_t i = 0; i < n; ++i) {
                XtWDZ.col(i) = Xis[i].transpose() * xis[i].cwiseProduct(wis[i]);
            }
            Zty_res = collapse_vec(y_res);
            d_out = collapse_vec(x);
            dw_out = collapse_vec(x.cwiseProduct(w));
            dwd_out = collapse_vec(x.cwiseProduct(w).cwiseProduct(x));

            ZtSigma_invZ_diag = Eigen::VectorXd::Zero(n);
            ZtDSigma_invDZ_diag = Eigen::VectorXd::Zero(n);
            ZtDSigma_invZ_diag = Eigen::VectorXd::Zero(n);
            for (size_t i = 0; i < n; ++i) {
                Eigen::VectorXd ZiTWi1 = Zis[i].transpose() * wis[i];
                ZtSigma_invZ_diag(i) = wis[i].sum() - ZiTWi1.dot(M_invs[i] * ZiTWi1);

                Eigen::VectorXd ZiTWixi = Zis[i].transpose() * wis[i].cwiseProduct(xis[i]);
                double tmp1 = (wis[i].array() * xis[i].array().square()).sum();
                double tmp2 = ZiTWixi.dot(M_invs[i] * ZiTWixi);
                ZtDSigma_invDZ_diag(i) = tmp1 - tmp2;
                ZtDSigma_invZ_diag(i) = wis[i].dot(xis[i]) - ZiTWixi.dot(M_invs[i] * ZiTWi1);
            }
        }

        GLMM_SC_INT(
            const Eigen::Ref<Eigen::MatrixXd> X_, 
            const Eigen::Ref<Eigen::VectorXd> y_, 
            const Eigen::Ref<Eigen::VectorXd> x_,
            const Eigen::Ref<Eigen::VectorXd> offset_,
            std::unique_ptr<Family> family_, 
            const Eigen::Ref<Eigen::VectorXd> ns_
        ) : 
            X(X_),
            y(y_),
            x(x_),
            offset(offset_),
            family(std::move(family_)),
            ns(ns_)
        {   
            N = X.rows();
            c = X.cols();
            n = ns.size();
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
            G = Eigen::MatrixXd::Zero(2, 2);
            tau = Eigen::VectorXd::Zero(3);
            tau(0) = 1.0;
            tau(1) = 0.5;
            tau(2) = 1.0;
            update_G();
            M_invs.resize(n);
            update_M_invs();
            Sigma_invX.resize(N, c);
            Py_tilde.resize(N);
            iter = 0;
            glmm_converged = false;
            step_size = 1;
        }

        void fit() {

            while (iter < max_iter) {
           
                update_M_invs();
                update_Sigma_invX();
                update_XtSigma_invX_inv();
                update_ZtSigma_invX();
                update_ZtDSigma_invX();
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
            update_ZtSigma_invX();
            update_ZtDSigma_invX();
            compute_output();
        }
};

#endif