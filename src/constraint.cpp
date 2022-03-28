// [[Rcpp::depends(BH)]]
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppParallel)]]
// [[Rcpp::depends(StanHeaders)]]

// [[Rcpp::plugins(cpp14)]]

#include <boost/math/tools/promotion.hpp>
#include <stan/math.hpp>
#include <RcppEigen.h>
#include <limits>
#include <cmath>
#include <string>
#include <unordered_map>

using namespace Rcpp;
using stan::math::abs;
using stan::math::pow;
using stan::math::sqrt;
using stan::math::log;
using boost::math::tools::promote_args;

template <typename T>
Eigen::Matrix<typename promote_args<T, double>::type, Eigen::Dynamic, 1>
constraint(const Eigen::VectorXd& y, const Eigen::Matrix<
  T, Eigen::Dynamic, 1> theta) {

  // Extract the parameters.
  const T w = theta(0);
  const T a = theta(1);
  const T b = theta(2);
  const T g = theta(3);
  const T n = theta(4);
  const T wbar = w / (1.0 - b);

  const double delta = pow(10, std::numeric_limits<double>::min_exponent10);
  typedef typename promote_args<T, double>::type TYPE;
  TYPE cond_1 = -1.0 * (a + g);
  TYPE cond_2 = 0.0;
  TYPE dphi = 0.0;

  double d = 0.0;

  for (auto t = 1; t < y.size(); t++) {
    d = 0;
    if (y(t-1) < 0) d = 1;
    dphi = (a + g*d) * (n + 1) * pow(y(t), 4);
    dphi /= pow((n - 2) * wbar + y(t)*y(t), 2);
    dphi += b;
    cond_2 += log(abs(dphi));
  }

  // Main condition
  cond_2 /= y.size();
  cond_2 += delta;

  // Result
  Eigen::Matrix<TYPE, Eigen::Dynamic, 1> res(2);
  res(0) = cond_1; res(1) = cond_2;

  return res;
}

//' Non-Linear Inequality Constraints
//'
//' @details This function evaluates the left-hand side of the constraints
//' \enumerate{
//'   \item \eqn{-(\alpha + \gamma) \le 0} and
//'   \item \eqn{\frac{1}{n} \sum_{t=1}^{n} \log \Lambda_{t} \left( \theta \right) + \delta.}
//'  }
//'
//' @param y numeric vector of observations \eqn{y_0, y_1, ...}
//' @param theta numeric vector of parameter values
//'
//' @returns Returns a list that contains the value of the left-hand side
//' (`"objective"`) of each constraint and the Jacobian ("jacobian") of the
//' constraint.
//'
// [[Rcpp::export]]
List cnstr(const Eigen::VectorXd y, const Eigen::VectorXd theta) {
  // Declare variables that store function value and gradient.
  Eigen::VectorXd fx;
  Eigen::MatrixXd jac_fx;

  // Gradient calculation.
  stan::math::jacobian(
    [&y](auto par) {
      return constraint(y, par);
    },
    theta, fx, jac_fx);

  return List::create(Named("objective") = fx, Named("jacobian") = jac_fx);
}



