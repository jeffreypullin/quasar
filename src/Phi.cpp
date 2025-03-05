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

#include <Eigen/Dense>
#include <boost/math/special_functions/digamma.hpp>
#include <boost/math/special_functions/trigamma.hpp>

double theta_score(Eigen::VectorXd y, Eigen::VectorXd mu, double theta) {
    
    double score = 0.0;
    for (int i = 0; i < y.size(); ++i) {
        score += boost::math::digamma(theta + y(i)) 
                 - boost::math::digamma(theta) 
                 + std::log(theta) 
                 + 1.0 
                 - std::log(theta + mu(i)) 
                 - (y(i) + theta) / (mu(i) + theta);
    }
    return score;
}

double theta_info(Eigen::VectorXd y, Eigen::VectorXd mu, double theta) {
    
    double info = 0.0;
    for (int i = 0; i < y.size(); ++i) {
        info -= boost::math::trigamma(theta + y(i)) 
                + boost::math::trigamma(theta) 
                - 1.0 / theta 
                + 2.0 / (mu(i) + theta) 
                - (y(i) + theta) / std::pow(mu(i) + theta, 2);
    }
    return info;
}

double estimate_phi_ml(Eigen::VectorXd y, Eigen::VectorXd mu) {

    // TODO: Make these parameters.
    int max_iter = 10;
    double tol = 1e-5;
    
    double phi;
    double theta;

    int iter = 0;
    double delta = 1;

    while (std::abs(delta) > tol) {

        theta = std::abs(theta);
        delta = theta_score(y, mu, theta) / theta_info(y, mu, theta);
        theta+= delta;

        // Too many iterations. 
        if (iter > max_iter) {
            break;
        }
        iter++;
    }

    phi = 1 / theta;
    return phi;
}

