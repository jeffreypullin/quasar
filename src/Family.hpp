/*
  Copyright (C) 2024-25 Jeffrey Pullin

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
