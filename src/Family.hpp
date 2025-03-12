/*
 * MIT License
 *
 * Copyright (c) 2024 Jeffrey Pullin
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

#ifndef FAMILY_H
#define FAMILY_H

#include <Eigen/Dense>
#include <iostream>

class Family {

    protected:
        std::string name;

    public: 
        virtual ~Family() = default;
        virtual Eigen::VectorXd link(Eigen::VectorXd mu) = 0;
        virtual Eigen::VectorXd invlink(Eigen::VectorXd eta) = 0;
        virtual Eigen::VectorXd init(Eigen::VectorXd y) = 0;
        virtual Eigen::VectorXd var(Eigen::VectorXd mu) = 0;
        virtual Eigen::VectorXd mu_eta(Eigen::VectorXd mu) = 0;
        virtual double dev(Eigen::VectorXd y, Eigen::VectorXd mu) = 0;
};

class Poisson : public Family {

    public: 
        Poisson() { name = "poisson"; }

        Eigen::VectorXd link(Eigen::VectorXd mu) {
            return mu.array().log();
        };
        
        Eigen::VectorXd invlink(Eigen::VectorXd eta) {
            return eta.array().exp();
        };

        Eigen::VectorXd init(Eigen::VectorXd y) {
             if ((y.array() < 0).any()) {
                throw std::invalid_argument("Poisson family requires non-negative values");
            }
            return y.array() + 0.1;
        };

        Eigen::VectorXd mu_eta(Eigen::VectorXd mu) {
            return mu.array().inverse();
        };

        Eigen::VectorXd var(Eigen::VectorXd mu) {
            return mu;
        };

        double dev(Eigen::VectorXd y, Eigen::VectorXd mu) {
            Eigen::VectorXd r(y.size());
            for (int i = 0; i < y.size(); ++i) {
                if (y(i) == 0)  {
                    r(i) = mu(i);
                } else {
                    r(i) = y(i) * log(y(i) / mu(i)) - (y(i) - mu(i));
                }
           }
           return 2 * r.sum();
        };
};

class NegativeBinomial : public Family {

    private:
        double phi;

    public:
        NegativeBinomial(double phi_) : phi(phi_) { 
            name = "negative_binomial"; 
            if (phi_ <= 0) {
                throw std::invalid_argument("NegativeBinomial family requires positive theta parameter");
            }
        }

        Eigen::VectorXd link(Eigen::VectorXd mu) {
            return mu.array().log();
        };
        
        Eigen::VectorXd invlink(Eigen::VectorXd eta) {
            return eta.array().exp();
        };

        Eigen::VectorXd init(Eigen::VectorXd y) {
            if ((y.array() < 0).any()) {
                throw std::invalid_argument("NegativeBinomial family requires non-negative values");
            }
            return y.array() + 0.1;
        };

        Eigen::VectorXd mu_eta(Eigen::VectorXd mu) {
            return mu.array().inverse();
        };

        Eigen::VectorXd var(Eigen::VectorXd mu) {
            return mu.array() + phi * mu.array().square();
        };

        double dev(Eigen::VectorXd y, Eigen::VectorXd mu) {
            Eigen::VectorXd r(y.size());
            double theta = 1 / phi;
            for (int i = 0; i < y.size(); ++i) {
                if (y(i) == 0)  {
                    r(i) = -(y(i) + theta) * std::log((y(i) + theta)/(mu(i) + theta));
                } else {
                    r(i) = y(i) * std::log(y(i) / mu(i)) - (y(i) + theta) * std::log((y(i) + theta)/(mu(i) + theta));
                }
            }
            return 2 * r.sum();
        };

};


#endif
