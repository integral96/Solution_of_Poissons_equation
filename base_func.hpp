#pragma once

#include <boost/mpl/bool.hpp>
#include <boost/mpl/at.hpp>
#include <boost/mpl/or.hpp>
#include <boost/mpl/if.hpp>
#include <boost/mpl/vector.hpp>
#include <boost/mpl/vector_c.hpp>
#include <type_traits>
#include <boost/type_traits/enable_if.hpp>
#include <boost/any.hpp>
#include <boost/hana.hpp>
#include <boost/proto/proto.hpp>

#include <array>
#include <random>
#include <ctime>
#include <complex>
#include <chrono>
#include <thread>

#define _SILENCE_ALL_CXX17_DEPRECATION_WARNINGS

namespace mpl = boost::mpl;
namespace hana = boost::hana;
namespace proto = boost::proto;
using boost::hana::literals::operator ""_c;

#define BASE_FUNCTION_COUNT 10

///Factorial

template<size_t N>
struct factorial {
    static constexpr size_t value = N*factorial<N - 1>::value;
};
template<>
struct factorial<0> {
    static constexpr size_t value = 1;
};

///Met func
///
// is vector for multiply vfnrbx on vector
template <typename F, size_t N, typename Vector>
struct is_vector : boost::mpl::false_ {};
template <typename F, size_t N>
struct is_vector<F, N, boost::array<F, N>> : boost::mpl::true_ {
    static constexpr size_t size = N;
};
template <typename F, size_t N>
struct is_vector<F, N, std::array<F, N>> : boost::mpl::true_ {
    static constexpr size_t size = N;
};

template<typename T> struct is_int : boost::mpl::false_ {};
template<> struct is_int<int> : boost::mpl::true_ {};
template<> struct is_int<unsigned> : boost::mpl::true_ {};
template<> struct is_int<size_t> : boost::mpl::true_ {};

template<typename T> struct is_complex_d : boost::mpl::false_ {};
template<> struct is_complex_d<std::complex<double>> : boost::mpl::true_ {};
template<typename T> struct is_complex_i : boost::mpl::false_ {};
template<> struct is_complex_i<std::complex<int>> : boost::mpl::true_ {};

///
///Вычисление степени
///

template<int N, typename T>
typename std::enable_if_t<(N < 0), T> pow_(const T& x) {
    return T(1)/pow_<-N>(x);
}
template<int N, typename T>
typename std::enable_if_t<(N == 0), T> pow_(const T& x) {
    return T(1);
}
template<int N, typename T>
typename std::enable_if_t<(N > 0) && (N%2 == 0), T> pow_(const T& x) {
    T p = pow_<N / 2>(x);
    return p*p;
}
template<int N, typename T>
typename std::enable_if_t<(N > 0) && (N%2 == 1), T> pow_(const T& x) {
    return pow_<N - 1>(x)*x;
}



/*!
 * struct meta_loop
 */
template <size_t N, size_t I, class Closure>
typename std::enable_if_t<(I == N)> is_meta_loop(Closure& closure) {}

template <size_t N, size_t I, class Closure>
typename std::enable_if_t<(I < N)> is_meta_loop(Closure& closure) {
    closure.template apply<I>();
    is_meta_loop<N, I + 1>(closure);
}
template <size_t N, class Closure>
void meta_loop(Closure& closure) {
    is_meta_loop<N, 0>(closure);
}
template <size_t N, class Closure>
void meta_loopUV(Closure& closure) {
    is_meta_loop<N, 1>(closure);
}
template <size_t N, size_t K, class Closure>
void meta_loop_KN(Closure& closure) {
    is_meta_loop<N, K>(closure);
}
///++
///
/*!
 * struct meta_loop_inv
 */
template <int N, int I, class Closure>
typename std::enable_if_t<(I == -1)> is_meta_loop_inv(Closure& closure) {}

template <int N, int I, class Closure>
typename std::enable_if_t<(I >= 0)> is_meta_loop_inv(Closure& closure) {
    closure.template apply<I>();
    is_meta_loop_inv<0, I - 1>(closure);
}
template <int N, class Closure>
void meta_loop_inv(Closure& closure) {
    is_meta_loop_inv<0, N>(closure);
}

/////Calculate Binom

template<size_t N, size_t K>
struct BC {
    static constexpr size_t value = factorial<N>::value / factorial<K>::value / factorial<N - K>::value;
};
/*!
 * struct abstract_sum
 */
template<class Closure>
struct abstract_sum_closures {
    typedef typename Closure::value_type value_type;
    abstract_sum_closures(Closure &closure) :  closure(closure), result(value_type()){}

    template<unsigned I>
    void apply(){
        result += closure.template value<I>();
    }
    Closure &closure;
    value_type result;
};

template<unsigned N, class Closure>
typename Closure::value_type abstract_sums(Closure &closure) {
    abstract_sum_closures<Closure> my_closure(closure);
    meta_loop<N>(my_closure);
    return my_closure.result;
}

/*!
 * struct abstract_subtract
 */
template<class Closure>
struct abstract_subtract_closures {
    typedef typename Closure::value_type value_type;
    abstract_subtract_closures(Closure &closure) :  closure(closure), result(value_type()){}

    template<unsigned I>
    void apply(){
        result -= closure.template value<I>();
    }
    Closure &closure;
    value_type result;
};

template<unsigned N, class Closure>
typename Closure::value_type abstract_subtract(Closure &closure) {
    abstract_subtract_closures<Closure> my_closure(closure);
    meta_loop<N>(my_closure);
    return my_closure.result;
}
/*!
 * struct abstract_mult
 */
template<class Closure>
struct abstract_multiple_closures {
    using value_type = typename Closure::value_type;
    abstract_multiple_closures(Closure &closure) : closure(closure), result(value_type(1)){}
    template<size_t I>
    void apply(){
        result *= closure.template value<I>();
    }
    Closure &closure;
    value_type result;
};
template<size_t K, class Closure>
typename Closure::value_type abstract_multiple(Closure &closure) {
    abstract_multiple_closures<Closure> my_closure(closure);
    meta_loop<K>(my_closure);
    return my_closure.result;
}

/*!
 * struct abstract_divide
 */
template<class Closure>
struct abstract_divide_closures {
    typedef typename Closure::value_type value_type;
    abstract_divide_closures(Closure &closure) :  closure(closure), result(value_type(1)){}

    template<unsigned I>
    void apply(){
        result /= closure.template value<I>();
    }
    Closure &closure;
    value_type result;
};

template<unsigned N, class Closure>
typename Closure::value_type abstract_divide(Closure &closure) {
    abstract_subtract_closures<Closure> my_closure(closure);
    meta_loop<N>(my_closure);
    return my_closure.result;
}

////Math Base && Expression && variable
///
///Random init
template<typename T, typename Matrix, typename = boost::enable_if_t<(Matrix::dimension > 0)>>
inline void gen_rand_matrix(Matrix& A, T min, T max) {
    std::time_t now = std::time(0);
    std::mt19937 gen{static_cast<std::uint32_t>(now)};
        if constexpr(Matrix::dimension == 4 && is_int<T>::value) {
            std::uniform_int_distribution<> dist{min, max};
            for(size_t i = 0; i < Matrix::row; ++i)
                for(size_t j = 0; j < Matrix::col; ++j)
                    for(size_t k = 0; k < Matrix::lin; ++k)
                        for(size_t l = 0; l < Matrix::vol; ++l)
                            A(i, j, k, l) = dist(gen);
        }
        if constexpr(Matrix::dimension == 4 && !is_int<T>::value) {
            std::uniform_real_distribution<> dist{min, max};
            for(size_t i = 0; i < Matrix::row; ++i)
                for(size_t j = 0; j < Matrix::col; ++j)
                    for(size_t k = 0; k < Matrix::lin; ++k)
                        for(size_t l = 0; l < Matrix::vol; ++l)
                            A(i, j, k, l) = dist(gen);
        }
        if constexpr(Matrix::dimension == 3 && is_int<T>::value) {
            std::uniform_int_distribution<> dist{min, max};
            for(size_t i = 0; i < Matrix::row; ++i)
                for(size_t j = 0; j < Matrix::col; ++j)
                    for(size_t k = 0; k < Matrix::lin; ++k)
                        A(i, j, k) = dist(gen);
        }
        if constexpr(Matrix::dimension == 3 && !is_int<T>::value) {
            std::uniform_real_distribution<> dist{min, max};
            for(size_t i = 0; i < Matrix::row; ++i)
                for(size_t j = 0; j < Matrix::col; ++j)
                    for(size_t k = 0; k < Matrix::lin; ++k)
                        A(i, j, k) = dist(gen);
        }
        if constexpr(Matrix::dimension == 2 && is_int<T>::value) {
            std::uniform_int_distribution<> dist{min, max};
            for(size_t i = 0; i < Matrix::row; ++i)
                for(size_t j = 0; j < Matrix::col; ++j)
                        A(i, j) = dist(gen);
        }
        if constexpr(Matrix::dimension == 2 && !is_int<T>::value) {
            std::uniform_real_distribution<> dist{min, max};
            for(size_t i = 0; i < Matrix::row; ++i)
                for(size_t j = 0; j < Matrix::col; ++j)
                        A(i, j) = dist(gen);
        }
}
template<size_t N, typename T, typename Array, typename = boost::enable_if_t<std::is_same_v<Array, boost::array<T, N>>>>
inline void gen_rand_array(Array& A, T min, T max, int n = 0) {
    auto start = std::chrono::system_clock::now();
    std::this_thread::sleep_for(std::chrono::milliseconds(n));
    auto end = std::chrono::system_clock::now();
    std::time_t end_time = std::chrono::system_clock::to_time_t(end);
    std::mt19937 gen{static_cast<std::uint32_t>(end_time)};
    std::uniform_real_distribution<> dist{min, max};
    for(size_t i = 0; i < N; ++i)
            A[i] = dist(gen);
}

