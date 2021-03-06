#pragma once

#include <omp.h>

#include <iostream>
#include <tuple>
#include <algorithm>
#include <cmath>
#include <vector>
#include <future>

#include <boost/any.hpp>
#include <boost/unordered_map.hpp>
#include <boost/functional/hash.hpp>
#include <boost/bind/bind.hpp>

#include "matrix.hpp"

using namespace boost::placeholders;

template<typename T, size_t N>
std::array<T, N> operator + (const std::array<T, N>& A, const std::array<T, N>& B) {
    std::array<T, N> tmp;
    for(size_t i = 0; i < N; ++i ) tmp[i] = A[i] + B[i];
    return tmp;
}
template<typename T, size_t N>
boost::array<T, N> operator + (const boost::array<T, N>& A, const boost::array<T, N>& B) {
    boost::array<T, N> tmp;
    for(size_t i = 0; i < N; ++i ) tmp[i] = A[i] + B[i];
    return tmp;
}

template<size_t M, size_t N, typename T = double>
struct OMEGA {
protected:
    typedef T value_type;
    value_type H1, H2;
public:
    constexpr OMEGA(value_type l1, value_type l2) : H1(l1/M), H2(l2/N) { }
};

template<size_t M, size_t N, typename Func_G, typename Func_Phi>
class PUASSON_SOLVER : OMEGA<M, N> {
private:
    typedef OMEGA<M, N> super;
    typedef typename OMEGA<M, N>::value_type value_type;
    typedef std::array<value_type, M> array_M;
    typedef typename boost::unordered_map<size_t, std::array<value_type, M>> vector_map;
    typedef typename boost::unordered_map<size_t, matrix_2<M, N, value_type>> matrix_map;
    typedef matrix_2<M, N, value_type> Matrix;
public:
    vector_map YY_1;
    vector_map FF_1;

    Matrix CC;
    Matrix CC_invert;

    matrix_map ALPHA;
    vector_map BETA;
    Matrix EDIN;

    array_M vec1;
    array_M vec2;
    array_M vec3;

    value_type H1, H2;
    value_type stp_H1, stp_H2;

private:
    Func_G   func_g;
    Func_Phi func_phi;

private:
    struct SOLV_ALPHA  {
        SOLV_ALPHA(PUASSON_SOLVER<M, N, Func_G, Func_Phi>& solv) : solv(solv) {}
        template<size_t K>
        void apply(){
            if constexpr(K == 0) solv.ALPHA[K] = solv.CC_invert;
            solv.ALPHA[K + 1] = ((solv.CC - solv.EDIN * solv.ALPHA[K])^(solv.CC - solv.EDIN * solv.ALPHA[K]))*solv.EDIN;
        }
    private:
        PUASSON_SOLVER<M, N, Func_G, Func_Phi>& solv;
    };
    struct SOLV_BETA {
        SOLV_BETA(PUASSON_SOLVER<M, N, Func_G, Func_Phi>& solv) : solv(solv) {}
        template<size_t K>
        void apply(){
            if constexpr(K == 0) solv.BETA[K] = solv.CC_invert*solv.FF_1[0];
            solv.BETA[K + 1] = ((solv.CC - solv.EDIN * solv.ALPHA[K])^(solv.CC - solv.EDIN * solv.ALPHA[K])) *
                    (solv.FF_1[K] + solv.BETA[K]);
        }
    private:
        PUASSON_SOLVER<M, N, Func_G, Func_Phi>& solv;
    };
    struct SOLV_YY {
        SOLV_YY(PUASSON_SOLVER<M, N, Func_G, Func_Phi>& solv) : solv(solv) {}
        template<int K>
        void apply(){
            if constexpr(K == M) solv.YY_1[K] = solv.BETA[K + 1];
            solv.YY_1[K] = solv.ALPHA[K + 1]*solv.YY_1[K + 1] + solv.BETA[K + 1];
        }
    private:
        PUASSON_SOLVER<M, N, Func_G, Func_Phi>& solv;
    };
    /////////////////////////////////////
    template<int J>
    struct Border_conditions_1 {
        Border_conditions_1(PUASSON_SOLVER<M, N, Func_G, Func_Phi>& solv) : solv(solv) {}
    template<int I>
    void apply(){
        if constexpr(I == 0) {
            solv.vec1[I] = (solv.H2*(solv.func_phi(solv.stp_H1*I, solv.stp_H2*J) +
                                          1/solv.H1*solv.func_g(solv.stp_H1*(I - 1), solv.stp_H2*J)));
        } else if constexpr(I == M) {
            solv.vec1[I] = (solv.H2*(solv.func_phi(solv.stp_H1*I, solv.stp_H2*J) +
                                                   1/solv.H1*solv.func_g(solv.stp_H1*(I + 1), solv.stp_H2*J)));
        } else {
            solv.vec1[I] = (solv.H2*solv.func_phi(solv.stp_H1*I, solv.stp_H2*J));
        }
    }
    private:
        PUASSON_SOLVER<M, N, Func_G, Func_Phi>& solv;
    };
    ///////////////////////////////////////
    template<int J>
    struct Border_conditions_2 {
        Border_conditions_2(PUASSON_SOLVER<M, N, Func_G, Func_Phi>& solv) : solv(solv) {}
    template<int I>
    void apply(){
        solv.vec2[I] = solv.func_g(solv.stp_H1*I, solv.stp_H2*0.);
    }
    private:
        PUASSON_SOLVER<M, N, Func_G, Func_Phi>& solv;
    };
    ///////////////////////////////////////
    template<int J>
    struct Border_conditions_3 {
        Border_conditions_3(PUASSON_SOLVER<M, N, Func_G, Func_Phi>& solv) : solv(solv) {}
    template<int I>
    void apply(){
        solv.vec3[I] = solv.func_g(solv.stp_H1*I, solv.stp_H2*N);
    }
    private:
        PUASSON_SOLVER<M, N, Func_G, Func_Phi>& solv;
    };
    ///////////////////////////////////////
    struct Border_conditions {
        Border_conditions(PUASSON_SOLVER<M, N, Func_G, Func_Phi>& solv) : solv(solv) {
        }
    template<int J>
    void apply(){
        if constexpr(J == 0) {
            Border_conditions_2<J> closure(solv);
            meta_loopUV<M + 1>(closure);
            solv.FF_1.insert(std::make_pair(J, solv.vec2));
        } if constexpr(J == N) {
            Border_conditions_3<J> closure(solv);
            meta_loopUV<M + 1>(closure);
            solv.FF_1.insert(std::make_pair(J, solv.vec3));
        } else {
            Border_conditions_1<J> closure(solv);
            meta_loop<M + 1>(closure);
            solv.FF_1.insert(std::make_pair(J, solv.vec1));
        }
    }
    private:
        PUASSON_SOLVER<M, N, Func_G, Func_Phi>& solv;
    };
public:
    constexpr PUASSON_SOLVER(const Func_G& func_g_, const Func_Phi& func_phi_, value_type l1, value_type l2) :
        func_g(func_g_), func_phi(func_phi_), super(l1, l2)
    {
        H1 = super::H1*super::H1;
        H2 = super::H2*super::H2;
        stp_H1 = super::H1;
        stp_H2 = super::H2;
        value_type alpha = H2/H1;
        EDIN.identity_matrix();
        CC.three_diag_matrix(alpha);
        CC_invert = CC^CC;
        //
//            Border_conditions

            struct Border_conditions closure_cond(*this);
            meta_loop<M + 1>(closure_cond);
            //
            // Solver
                struct SOLV_ALPHA closure_alpha(*this);
                meta_loop<M + 1>(closure_alpha);
                struct SOLV_BETA closure_beta(*this);
                meta_loop<M + 1>(closure_beta);
                struct SOLV_YY closure_YY(*this);
                meta_loop_inv<M + 1>(closure_YY);
    }

    const std::array<value_type, M>& YY_(size_t i) {
        return YY_1[i];
    }

};

