// [[Rcpp::depends(BH)]]
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppParallel)]]
// [[Rcpp::depends(StanHeaders)]]

// [[Rcpp::plugins(cpp14)]]

#include <boost/math/tools/promotion.hpp>
#include <stan/math.hpp>
#include <RcppEigen.h>
#include <cmath>
#include <string>
#include <unordered_map>

using namespace Rcpp;
using stan::math::multiply;
using stan::math::divide;
using stan::math::pow;
using stan::math::sqrt;
using stan::math::log;
using stan::math::student_t_lpdf;
using boost::math::tools::promote_args;

// Student's t log-density i.t.o. variance.
template <typename T1, typename T2>
typename promote_args<T1, T2>::type log_pdf(const T1 x, const T2 f,
                                            const T2 nu) {
  // Define the promoted type
  typedef typename promote_args<T1, T2>::type TYPE;
  // Calculate the scale parameter such that V[Y] = f
  const TYPE phi = sqrt(divide(nu - 2.0, nu));
  const TYPE sigma = multiply(sqrt(f), phi);
  // Calculate the log-density of Y~T with df = nu and V[Y] = f.
  return student_t_lpdf(x, nu, 0.0, sigma);
}

// score, i.t.o. epsilon
template <typename T1, typename T2>
typename promote_args<T1, T2>::type qscore(const T1 e_t, const T2 n) {
  const auto sqrd_dev = e_t * e_t;
  const auto numer = (n + 1.0) * sqrd_dev;
  const auto denom = (n - 2.0) + sqrd_dev;
  return numer / denom;
}

// recursion
template <typename T>
std::unordered_map<std::string, Eigen::Matrix<T, Eigen::Dynamic, 1>>
 recursion_map(const Eigen::VectorXd& y, const double f_0,
               const Eigen::Matrix<T, Eigen::Dynamic, 1> theta) {
  // Extract the parameters.
  const T w = theta(0);
  const T a = theta(1);
  const T b = theta(2);
  const T g = theta(3);
  const T n = theta(4);

  // Prepare the vectors to store results.
  Eigen::Matrix<T, Eigen::Dynamic, 1> ll(y.size());
  ll.fill(0.0);
  Eigen::Matrix<T, Eigen::Dynamic, 1> f(y.size());
  f.fill(w);
  f(0) = f_0;
  Eigen::Matrix<T, Eigen::Dynamic, 1> e(y.size());
  e.fill(y(0) / sqrt(f_0));
  Eigen::Matrix<T, Eigen::Dynamic, 1> u(y.size());
  u.fill(qscore(e(0), n));
  double d = 0.0;

  for (auto t = 1; t < y.size(); t++) {
    d = y(t-1) < 0 ? 1 : 0;
    f(t) += ((a + g*d) * u(t-1) + b) * f(t-1);
    e(t) = y(t) / sqrt(f(t));
    u(t) = qscore(e(t), n);
    ll(t) = log_pdf(y(t), f(t), n);
  }

  std::unordered_map<std::string, Eigen::Matrix<T, Eigen::Dynamic, 1>> res = {
    {"f", f}, {"llik_t", ll}, {"s_t", u}, {"e_t", e}};

  return res;
}

//' Beta-t-GARCH Recursion/Filter
//'
//' @details
//'   Calculates the values of
//'   \itemize{
//'     \item each observation's contribution to the log-likelihood,
//'     \item the conditional variance,
//'     \item score innovation/residual and
//'     \item model residuals.
//'   }
//'   End-users should rather use the \code{BetaTGARCHfit} class's \code{filter}
//'   method.
//'
//' @param y numeric vector of observations \eqn{y_0, y_1, ...}
//' @param f_0 initial value of the conditional variance
//' @param theta numeric vector of parameter values
//'
//' @returns
//'   Returns a list of vectors. The elements of the list are (in order):
//'   \enumerate{
//'     \item each observation's contribution to the log-likelihood (`"llik_t"`),
//'     \item the conditional variances (`"f_t"`),
//'     \item score residuals (`"s_t"`) and
//'     \item model residuals (`"residuals"`).
//'   }
//'
// [[Rcpp::export]]
List recursion_lst(const Eigen::VectorXd& y, const double f_0,
                   const Eigen::VectorXd theta) {
  auto res = recursion_map(y, f_0, theta);
  return List::create(
    Named("llik_t") = res["llik_t"],
    Named("f_t") = res["f"],
    Named("s_t") = res["s_t"],
    Named("residuals") = res["e_t"]);
}

//' Negative log-likelihood and its gradient.
//'
//' @details
//'   Calculates the negative log-likelihood of the Beta-GARCH(1, 1) model
//'   with leverage and its gradient w.r.t. to vector of parameters.
//'
//' @param y numeric vector of observations \eqn{y_0, y_1, ...}
//' @param f_0 initial value of the conditional variance
//' @param theta numeric vector of parameter values
//'
//' @returns
//'   Returns a list that contains the negative log-likelihood (`"objective"`)
//'   and the gradient (`"gradient"`) of the negative log-likelihood w.r.t. the
//'   model parameters.
//'
// [[Rcpp::export]]
List nll(const Eigen::VectorXd y, const double f_0,
         const Eigen::VectorXd theta) {
  // Declare variables that store function value and gradient.
  double fx;
  Eigen::VectorXd grad_fx;

  // Gradient calculation.
  stan::math::gradient(
    [&y, &f_0](auto par) {
      auto res = recursion_map(y, f_0, par);
      auto llik_t = res["llik_t"];
      auto nll = -1.0 * stan::math::sum(llik_t);
      return nll;
    },
    theta, fx, grad_fx);

  return List::create(Named("objective") = fx, Named("gradient") = grad_fx);
}

//' Simulate Values From Beta-t-GARCH Model Given a Vector of Innovations
//'
//' @details Simulates observations from a Beta-t-GARCH model.
//'
//' @param e numeric vector of zero-mean and unit-variance innovations
//' @param f_0 initial value of the conditional variance
//' @param theta numeric vector of parameter values
//'
//' @returns
//'   Returns a list of vectors. The elements of the list are (in order):
//'   \itemize{
//'    \item the simulated observations (`"y_t"`),
//'    \item the conditional variance (`"f_t"`),
//'    \item score innovations (`"s_t"`),
//'    \item the contributions to the log-likelihood (`"llik_t"`) and
//'    \item the original series of innovations (`"innovations"`).
//'   }
//'
// [[Rcpp::export]]
List simulate_lst(const Eigen::VectorXd& e, const double f_0,
                  const Eigen::VectorXd theta) {
  // Extract the parameters.
  const double w = theta(0);
  const double a = theta(1);
  const double b = theta(2);
  const double g = theta(3);
  const double n = theta(4);

  // Prepare the vectors to store results.
  Eigen::VectorXd y(e.size());
  Eigen::VectorXd f(e.size());
  Eigen::VectorXd u(e.size());
  Eigen::VectorXd ll(e.size());

  // Initial values
  f.fill(w);
  f(0) = f_0;
  u.fill(qscore(e(0), n));

  double d = 0.0;

  for (auto t = 1; t < y.size(); t++) {
    d = y(t-1) < 0 ? 1 : 0;
    f(t) += ((a + g*d) * u(t-1) + b) * f(t-1);
    ll(t) = log_pdf(y(t), f(t), n);
    y(t) = e(t) * sqrt(f(t));
    u(t) = qscore(e(t), n);
  }

  auto res = List::create(Named("y_t") = y, Named("f_t") = f, Named("s_t") = u,
                          Named("llik_t") = ll, Named("innovations") = e);
  return res;
}
