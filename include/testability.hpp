#ifndef TESTABILITY_HPP
#define TESTABILITY_HPP

#include "config.hpp"

#include <vector>
#include <string>
#include <functional>
#include <cmath>
#include <utility>
#include <type_traits>

#include <boost/math/tools/roots.hpp>
#include <boost/math/special_functions/binomial.hpp>
#include <boost/assert.hpp>

#include <range/v3/view/transform.hpp>
#include <range/v3/view/filter.hpp>
#include <range/v3/view/remove_if.hpp>
#include <range/v3/view/iota.hpp>
#include <range/v3/view/reverse.hpp>
#include <range/v3/view/map.hpp>
#include <range/v3/view/adjacent_remove_if.hpp>
#include <range/v3/view/zip.hpp>
#include <range/v3/numeric/inner_product.hpp>
#include <range/v3/algorithm/for_each.hpp>
#include <range/v3/algorithm/max.hpp>
#include <range/v3/algorithm/sort.hpp>
#include <range/v3/algorithm/max_element.hpp>
#include <range/v3/algorithm/find_if.hpp>
#include <range/v3/algorithm/find_if_not.hpp>
#include <range/v3/begin_end.hpp>
#include <range/v3/size.hpp>

using Frequency = unsigned;

namespace thesis
{
  // In: n1 = D1.size, n2 = D2.size, f = frequency
  double min_p_(const unsigned n1, const unsigned n2, const Frequency f);

  template<typename Algorithm>
  unsigned
  one_pass_(Algorithm alg, const unsigned n1, const unsigned n2, const double alpha);

  template<typename Algorithm>
  unsigned
  lamp_dec_(Algorithm m, const unsigned n1, const unsigned n2, const double alpha);

  template<typename Algorithm>
  unsigned
  early_term_(Algorithm m_et, const unsigned n1, unsigned n2, const double alpha);

  template<typename Algorithm>
  unsigned
  bis_leap_(Algorithm m_et, const unsigned n1, const unsigned n2, const double alpha);


} // namespace thesis

/**************************************************************************************/

namespace th = thesis;
// In: n1 = D1.size, n2 = D2.size, f = frequency. Requires[f <= n1 + n2]
// Out: (double) phi(f) 
double th::min_p_(const unsigned n1, const unsigned n2, const Frequency f)
{

  using namespace ranges;

  auto N = n1 + n2;
  BOOST_ASSERT(f < N); // precondition
  BOOST_ASSERT(n1 > 0);
  auto op2 = [](const auto a, const auto b) { // correctly compute n/m
    return static_cast<double>(a) / static_cast<double>(b);};

  auto init_v = static_cast<double>(n1) / static_cast<double>(N); // init value

  auto base1 = f > n1 ? 1u : n1 - f + 1;
  auto base2 = f > n1 ? N - n1 + 1 : N - f + 1;

  auto r1 = view::iota(base1,n1); // {base1, base1+1, ..., n1}
  auto r2 = view::iota(base2,N); // {base2, base2+1, ..., N}
  return inner_product(r1, r2, init_v, std::multiplies<>{}, op2); // calculate base1/base2 x (base1 + 1)/(base2 + 1) x ... x n1/N


}

template<typename Algorithm>
unsigned
th::one_pass_(Algorithm alg, const unsigned n1, const unsigned n2, const double alpha) 
{
  using namespace ranges;

  const auto phi = [=](const auto f){return th::min_p_(n1, n2, f);}; // minimum p-value function
  auto freq_rng = view::iota(1u,n1) | view::remove_if([=](const auto f){
      return phi(f) > alpha;}); 
  const auto min_freq = *begin(freq_rng); // min. admissible frequency for mining at significance level: alpha


  auto freqs = alg(min_freq); BOOST_ASSERT(size(freqs) > 0); //run the mining algorithm
  const auto cmp = [](const auto a, const auto b){return a > b;};
  sort(freqs,cmp); // then sort (increasing order)

  auto map_ = view::zip(view::iota(1ul, size(freqs)), freqs) |
    view::adjacent_remove_if([](const auto& a, const auto& b){
	return std::get<1>(a) == std::get<1>(b);}); // enumerate the mapping (freq, m(freq))

  const auto pred = [=](const auto& p){
    return phi(std::get<1>(p)) * std::get<0>(p) > alpha;}; // search for the SIGMA root frequency

  return std::get<1>(*find_if(map_, pred)) + 1; // and return it
}

template<typename Algorithm>
unsigned
th::lamp_dec_(Algorithm m, const unsigned n1, const unsigned n2, const double alpha) 
{
  using namespace ranges;
  const auto freq_rng = view::iota(1u,n1) | view::reverse; // {n1, n1 - 1, ..., 1}
  const auto phi = [=](const auto f){return th::min_p_(n1, n2, f);}; // minimum p-value function
  const auto pred = [=](const auto f){
    const auto m_ = m(f); BOOST_ASSERT(m_ > 0); // <-- very important!
    return m_ * phi(f) > alpha;}; 

  // find the first value f for which m(f) * phi(f) <= alpha is satisfied and return it  
  return *find_if(freq_rng, pred) + 1; 
}

template<typename Algorithm>
unsigned
th::early_term_(Algorithm m_et, const unsigned n1, const unsigned n2, const double alpha){
  using namespace ranges;
  const auto phi = [=](const auto f){return th::min_p_(n1, n2, f);};
  const auto pred = [=](const auto f){
    const auto pv = phi(f);
    auto m =  m_et(f, pv); BOOST_ASSERT(m > 0); // <-- very important!
    return m * pv <= alpha;
  }; // continue the meaning process until m * pv <= alpha is satisfied

  // min. admissible range of frequency  
  auto freq_rng = view::iota(1,n1)
    | view::remove_if([=](const auto f){return phi(f) > alpha;});
  
  BOOST_ASSERT(!pred(*begin(freq_rng))); // <-- starting frequency too large

  // return the first freq that satisfies the predicate pred  
  return *find_if(freq_rng, pred);
}

template<typename Algorithm>
unsigned
th::bis_leap_(Algorithm m_et, const unsigned n1, const unsigned n2, const double alpha) 
{
  using namespace ranges;
  using namespace boost::math;
    
  const auto phi = [=](const auto f){return th::min_p_(n1, n2, f);};
  const auto g = [=](const auto f){ // function whose root is to be found
    const auto pv = phi(f);
    const auto m = m_et(f, pv); BOOST_ASSERT(m > 0); // <-- very important!
    return m * pv - alpha;};
  auto freq_rng = view::iota(1u,n1)| view::remove_if([=](const auto f){
      return phi(f) > alpha;}); // min. admissible frequency range

  // setup [freq_min, freq_max] interval
  const auto min = static_cast<double>(*begin(freq_rng));
  const auto max = static_cast<double>(n1);

  // termination condition freq_max - freq_min == 1
  const auto tol = [](const auto a, const auto b){
    return b - a <= 1;};

  // apply the bisection method and return freq_min 
  const auto res = tools::bisect(g, min, max, tol);
  return std::floor(res.first);
}
// }   

#endif
