#pragma once

#include "base_func.hpp"

#include <iostream>
#include <iomanip>
#include <algorithm>
#include <variant>
#include <cmath>



template<typename Iter>
constexpr bool is_input_iterator =
        hana::is_convertible<typename std::iterator_traits<Iter>::iterator_category, std::input_iterator_tag>{};
template<typename Iter>
constexpr bool is_random_iterator =
        hana::is_convertible<typename std::iterator_traits<Iter>::iterator_category, std::random_access_iterator_tag>{};

//constexpr auto is_3d_hana = hana::metafunction_class<matrix_3<N1, N2, N3, int>>;

struct MatrixGrammar : proto::or_<
        proto::terminal<proto::_>,
        proto::plus<MatrixGrammar, MatrixGrammar>,
        proto::minus<MatrixGrammar, MatrixGrammar>,
        proto::multiplies<MatrixGrammar, MatrixGrammar>

        >
{};
template<typename Expr>
struct MatrixExpr;
struct MatrixDomain : proto::domain<proto::generator<MatrixExpr>, MatrixGrammar> {};
struct SubscriptCntxt : proto::callable_context<const SubscriptCntxt> {
    typedef double result_type;
    int i, j;
    SubscriptCntxt(int i, int j) : i(i), j(j) {}
    template<class Matrix>
    double operator()(proto::tag::terminal, const Matrix& mtrx) const {
        return mtrx(i, j);
    }
    template<typename E1, typename E2>
    double operator()(proto::tag::plus, const E1& e1, const E2& e2) const {
        return proto::eval(e1, *this) + proto::eval(e2, *this);
    }
    template<typename E1, typename E2>
    double operator()(proto::tag::minus, const E1& e1, const E2& e2) const {
        return proto::eval(e1, *this) - proto::eval(e2, *this);
    }
    template<typename E1, typename E2>
    double operator()(proto::tag::multiplies, const E1& e1, const E2& e2) const {
        return proto::eval(e1, *this) * proto::eval(e2, *this);
    }
};
template<typename Expr>
struct MatrixExpr : proto::extends<Expr, MatrixExpr<Expr>, MatrixDomain> {
    explicit MatrixExpr(const Expr& e)
        : proto::extends<Expr, MatrixExpr<Expr>, MatrixDomain>(e) {
    }

    typename proto::result_of::eval< Expr, SubscriptCntxt>::type
    operator ()(int i, int j) const {
        const SubscriptCntxt ctx(i, j);
        return proto::eval(*this, ctx);
    }
};

///2-х мерная матрица.
///
///

template<size_t N, size_t M, typename T>
using Matrix2D_Type = std::array<std::array<T, M>, N>;

template<size_t N, size_t M, typename T>
struct matrix_2 : public Matrix2D_Type<N, M, T> {
    static constexpr size_t dimension = 2;
    static constexpr size_t row = N;
    static constexpr size_t col = M;
    typedef Matrix2D_Type<N, M, T> type;
    typedef typename type::iterator iterator;
    typedef typename type::const_iterator  constant_iterator;
    typedef T value_type;
    type MTRX;
public:
    matrix_2() : MTRX() {}

    matrix_2(std::initializer_list<std::initializer_list<value_type>> listlist) : MTRX() {
        auto rows = listlist.size();
        auto cols = listlist.begin()->size();
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                MTRX[i][j] = ((listlist.begin()+i)->begin())[j];
            }
        }
    }
    void fill(T x) {
        for (size_t i = 0; i < N; ++i)
            for (size_t j = 0; j < M; ++j) MTRX[i][j] = x;
    }
    matrix_2<N, M, T>& identity_matrix() {
        calc_identity_matrix(*this);
        return *this;
    }
    void three_diag_matrix(value_type alpha) {
        calc_three_diag_matrix(*this, alpha);
    }
    T& operator () (size_t i, size_t j) {
        return MTRX[i][j];
    }
    T const& operator () (size_t i, size_t j) const {
        return MTRX[i][j];
    }

    T& at (size_t i, size_t j) {
        return MTRX[i][j];
    }
    matrix_2<N, M, T>& operator = (matrix_2<N, M, T> const& other) {
        MTRX = other.MTRX;
        return *this;
    }
    template<typename Expr>
    matrix_2<N, M, T>& operator = ( const Expr& expr ) {
        for(int i = 0; i < N; ++i) {
            for(int j = 0; j < M; ++j) {
                const SubscriptCntxt ctx(i, j);
                MTRX[i][j] = proto::eval(proto::as_expr<MatrixDomain>(expr), ctx);
            }
        }
        return *this;
    }
    friend matrix_2<N, M, T> operator - (const matrix_2<N, M, T>& other1, const matrix_2<N, M, T>& other2) {
        BOOST_STATIC_ASSERT(N == M);
        matrix_2<N, M, T> tmp;
        calc_minus_matrix(tmp, other1, other2);
        return tmp;
    }
    friend matrix_2<N, M, T> operator + (const matrix_2<N, M, T>& other1, const matrix_2<N, M, T>& other2) {
        BOOST_STATIC_ASSERT(N == M);
        matrix_2<N, M, T> tmp;
        calc_plus_matrix(tmp, other1, other2);
        return tmp;
    }
    friend matrix_2<N, N, T> operator * (const matrix_2<N, M, T>& other1, const matrix_2<M, N, T>& other2) {
        matrix_2<N, N, T> tmp;
        calc_multy_matrix(tmp, other1, other2);
        return tmp;
    }
    template<class Vector>
    friend std::enable_if_t<(is_vector<T, N, Vector>::value), Vector> operator * (const matrix_2<N, M, T>& other1, const Vector& vector) {
        Vector tmp;
        calc_multy_vector(tmp, other1, vector);
        return tmp;
    }
    template<typename F>
    friend std::enable_if_t<std::is_convertible<T, F>::value, matrix_2<N, M, T>> operator / (matrix_2<N, M, T>& other1, const F scalar) {
        calc_divide_on_scalar(other1, scalar);
        return other1;
    }
    template<typename F>
    friend std::enable_if_t<std::is_convertible<T, F>::value, matrix_2<N, M, T>> operator * (matrix_2<N, M, T>& other1, const F scalar) {
        calc_multiply_on_scalar(other1, scalar);
        return other1;
    }
    matrix_2<N, M, T> operator ^(const matrix_2<N, M, T>& ORIG) { ///Invers matrix
        BOOST_STATIC_ASSERT(N == M);
        matrix_2<N, M, T> RESULT;
        matrix_2<N, 2*M, T> TMP;
        calc_multy_2N(TMP, ORIG);
        calc_invers_tmp_matrix(TMP);
        calc_invers_matrix(TMP);
        calc_multy_toN(RESULT, TMP);
        return RESULT;
    }
    friend std::ostream& operator << (std::ostream& os, const matrix_2<N, M, T>& A){
        for (size_t i = 0; i < N; ++i) {
            for (size_t j = 0; j < M; ++j) {
//                std::visit([&](auto&& arg) {
//                    os << arg << "\t";
//                },
                os << std::round(A.MTRX[i][j]*10000)/10000<< "\t";
            }
            os << std::endl;
        }
        return os;
    }
};


///З-х мерная матрица.
///
template<size_t N, size_t M, size_t K, typename T>
using Matrix3D_Type = std::array<std::array<std::array<T, K>, M>, N>;

template<size_t N, size_t M, size_t K, typename T>
struct matrix_3 : public Matrix3D_Type<N, M, K, T> {
    static constexpr size_t dimension = 3;
    static constexpr size_t row = N;
    static constexpr size_t col = M;
    static constexpr size_t lin = K;
    typedef Matrix3D_Type<N, M, K, T> super;
    typedef typename super::iterator iterator;
    typedef typename super::const_iterator  constant_iterator;
    typedef T value_type;
private:
    super MTRX;
public:
    matrix_3() : MTRX() {}

    iterator begin() {
        return MTRX.begin();
    }
    iterator end() {
        return MTRX.end();
    }
    constant_iterator cbegin() {
        return MTRX.cbegin();
    }
    constant_iterator cend() {
        return MTRX.cend();
    }
    T& operator () (size_t i, size_t j, size_t k) {
        return MTRX[i][j][k];
    }
    T const& operator () (size_t i, size_t j, size_t k) const {
        return MTRX[i][j][k];
    }
    T& at (size_t i, size_t j, size_t k) {
        return MTRX[i][j][k];
    }
    matrix_3<N, M, K, T>& operator = (matrix_3<N, M, K, T> const& other) {
        MTRX = other.MTRX;
        return *this;
    }
    constexpr bool operator == (matrix_3<N, M, K, T> const& other) const {
        return MTRX == other.MTRX;
    }
    constexpr bool operator != (matrix_3<N, M, K, T> const& other) const {
        return MTRX != other.MTRX;
    }
    friend std::ostream& operator << (std::ostream& os, const matrix_3<N, M, K, T>& A){
        for(const auto& x : A.MTRX) {
            for(const auto& y : x) {
                for(const auto& z : y) os << z << "\t";
                os << std::endl;
            } os << std::endl;
        } os << std::endl;
        return os;
    }
};

///4-х мерная матрица.
///

template<size_t N, size_t M, size_t K, size_t L, typename T>
using Matrix4D_Type = std::array<std::array<std::array<std::array<T, L>, K>, M>, N>;

template<size_t N, size_t M, size_t K, size_t L, typename T>
struct matrix_4 : public Matrix4D_Type<N, M, K, L, T> {
    static constexpr size_t dimension = 4;
    static constexpr size_t row = N;
    static constexpr size_t col = M;
    static constexpr size_t lin = K;
    static constexpr size_t vol = L;
    typedef Matrix4D_Type<N, M, K, L, T> super;
    typedef typename super::iterator iterator;
    typedef typename super::const_iterator  constant_iterator;
    typedef T value_type;
private:
    super MTRX;
public:
    matrix_4() : MTRX() {}
    iterator begin() {
        return MTRX.begin();
    }
    iterator end() {
        return MTRX.end();
    }
    constant_iterator cbegin() {
        return MTRX.cbegin();
    }
    constant_iterator cend() {
        return MTRX.cend();
    }
    T& operator () (size_t i, size_t j, size_t k, size_t l) {
        return MTRX[i][j][k][l];
    }
    T const& operator () (size_t i, size_t j, size_t k, size_t l) const {
        return MTRX[i][j][k][l];
    }
    T& at (size_t i, size_t j, size_t k, size_t l) {
        return MTRX[i][j][k][l];
    }
    matrix_4<N, M, K, L, T>& operator = (matrix_4<N, M, K, L, T> const& other) {
        MTRX = other.MTRX;
        return *this;
    }
    friend std::ostream& operator << (std::ostream& os, const matrix_4<N, M, K, L, T>& A){
        for(const auto& x : A.MTRX) {
            for(const auto& y : x) {
                for(const auto& z : y) {
                    for(const auto& f : z) os << f << "\t";
                    os << std::endl;
                } os << std::endl;
            } os << std::endl;
        } os << std::endl;
        return os;
    }
};

///Minus matrix

template <size_t i, typename Matrix>
struct matrix_minus_1
{
    typedef typename Matrix::value_type value_type;
private:
    const Matrix& A;
    const Matrix& B;
    Matrix& RESULT;
public:
    matrix_minus_1(Matrix& RESULT_, const Matrix& A_, const Matrix& B_) : RESULT(RESULT_), A(A_), B(B_)  {}

    template <size_t j>
    void apply() {
        RESULT(i, j) = A(i, j) - B(i, j);
    }
};
template <typename Matrix>
struct matrix_minus_2
{
    typedef typename Matrix::value_type value_type;
private:
    const Matrix& A;
    const Matrix& B;
    Matrix& RESULT;
public:
    matrix_minus_2(Matrix& RESULT_, const Matrix& A_, const Matrix& B_) : RESULT(RESULT_), A(A_), B(B_)  {}

    template <size_t i>
    void apply() {
        matrix_minus_1<i, Matrix> closure(RESULT, A, B);
        meta_loop<Matrix::row>(closure);
    }
};
template <typename Matrix>
inline void calc_minus_matrix(Matrix& RESULT, const Matrix& A, const Matrix& B) {
    matrix_minus_2<Matrix> closure(RESULT, A, B);
    meta_loop<Matrix::col>(closure);
}
////close operator minus
///
///Plus matrix

template <size_t i, typename Matrix>
struct matrix_plus_1
{
    typedef typename Matrix::value_type value_type;
private:
    const Matrix& A;
    const Matrix& B;
    Matrix& RESULT;
public:
    matrix_plus_1(Matrix& RESULT_, const Matrix& A_, const Matrix& B_) : RESULT(RESULT_), A(A_), B(B_)  {}

    template <size_t j>
    void apply() {
        RESULT(i, j) = A(i, j) + B(i, j);
    }
};
template <typename Matrix>
struct matrix_plus_2
{
    typedef typename Matrix::value_type value_type;
private:
    const Matrix& A;
    const Matrix& B;
    Matrix& RESULT;
public:
    matrix_plus_2(Matrix& RESULT_, const Matrix& A_, const Matrix& B_) : RESULT(RESULT_), A(A_), B(B_)  {}

    template <size_t i>
    void apply() {
        matrix_plus_1<i, Matrix> closure(RESULT, A, B);
        meta_loop<Matrix::row>(closure);
    }
};
template <typename Matrix>
void calc_plus_matrix(Matrix& RESULT, const Matrix& A, const Matrix& B) {
    matrix_plus_2<Matrix> closure(RESULT, A, B);
    meta_loop<Matrix::col>(closure);
}
////close operator plus
///
///DIVIDE matrix on scalar

template <typename T, size_t i, typename Matrix>
struct matrix_divide_1
{
private:
    Matrix& RESULT;
    T scalar;
public:
    matrix_divide_1(Matrix& RESULT_, T scalar_) : RESULT(RESULT_), scalar(scalar_)  {}

    template <size_t j>
    void apply() {
        RESULT(i, j) /= scalar;
    }
};
template <typename T, typename Matrix>
struct matrix_divide_2
{
private:
    Matrix& RESULT;
    T scalar;
public:
    matrix_divide_2(Matrix& RESULT_, T scalar_) : RESULT(RESULT_), scalar(scalar_)   {}

    template <size_t i>
    void apply() {
        matrix_divide_1<T, i, Matrix> closure(RESULT, scalar);
        meta_loop<Matrix::row>(closure);
    }
};
template <typename T, typename Matrix>
void calc_divide_on_scalar(Matrix& RESULT, T scalar) {
    matrix_divide_2<T, Matrix> closure(RESULT, scalar);
    meta_loop<Matrix::col>(closure);
}
////close operator devide on scalar
///
///multiply matrix on scalar

template <typename T, size_t i, typename Matrix>
struct matrix_multiply_1
{
private:
    Matrix& RESULT;
    T scalar;
public:
    matrix_multiply_1(Matrix& RESULT_, T scalar_) : RESULT(RESULT_), scalar(scalar_)  {}

    template <size_t j>
    void apply() {
        RESULT(i, j) *= scalar;
    }
};
template <typename T, typename Matrix>
struct matrix_multiply_2
{
private:
    Matrix& RESULT;
    T scalar;
public:
    matrix_multiply_2(Matrix& RESULT_, T scalar_) : RESULT(RESULT_), scalar(scalar_)   {}

    template <size_t i>
    void apply() {
        matrix_multiply_1<T, i, Matrix> closure(RESULT, scalar);
        meta_loop<Matrix::row>(closure);
    }
};
template <typename T, typename Matrix>
void calc_multiply_on_scalar(Matrix& RESULT, T scalar) {
    matrix_divide_2<T, Matrix> closure(RESULT, scalar);
    meta_loop<Matrix::col>(closure);
}
////close operator multiply on scalar
///
///
///
///identity matrix

template <size_t i, typename Matrix>
struct identity_matrix_1
{
    typedef typename Matrix::value_type value_type;
private:
    Matrix& RESULT;
public:
    identity_matrix_1(Matrix& RESULT_) : RESULT(RESULT_)  {}

    template <size_t j>
    void apply() {
        if constexpr (i == j) RESULT(i, j) = value_type(1);
        else RESULT(i, j) = value_type(0);
    }
};
template <typename Matrix>
struct identity_matrix_2
{
    typedef typename Matrix::value_type value_type;
private:
    Matrix& RESULT;
public:
    identity_matrix_2(Matrix& RESULT_) : RESULT(RESULT_)  {}

    template <size_t i>
    void apply() {
        identity_matrix_1<i, Matrix> closure(RESULT);
        meta_loop<Matrix::row>(closure);
    }
};
template <typename Matrix>
inline void calc_identity_matrix(Matrix& RESULT) {
    identity_matrix_2<Matrix> closure(RESULT);
    meta_loop<Matrix::col>(closure);
}
////close identity_matrix
/*!
 *  Умножение матриц
 */



template <size_t I, size_t J, typename MatrixA, typename MatrixB>
struct matrix_prod_closure
{
    typedef typename MatrixA::value_type value_type;
private:
    const MatrixA& A;
    const MatrixB& B;
public:
    matrix_prod_closure(const MatrixA& A, const MatrixB& B) : A(A), B(B) {}
    template<size_t K>
    value_type value() const {
        return A(I, K) * B(K, J);
    }
};

template <size_t i, typename MatrixC, typename MatrixA, typename MatrixB>
struct A_m_1
{
private:
    const MatrixA& A;
    const MatrixB& B;
    MatrixC& C;

public:
    A_m_1(MatrixC& C, const MatrixA& A, const MatrixB& B): C(C), A(A), B(B){}

    template <size_t j>
    void apply() {
        matrix_prod_closure<i, j, MatrixA, MatrixB> closure(A, B);
        C(i, j) = abstract_sums<MatrixA::col>(closure);
    }
};

template <typename MatrixC, typename MatrixA, typename MatrixB>
struct CALC_A
{
private:
    const MatrixA& A;
    const MatrixB& B;
    MatrixC& C;
public:
    CALC_A(MatrixC& C, const MatrixA& A, const MatrixB& B): C(C), A(A), B(B){}

    template <size_t i>
    void apply() {
        A_m_1<i, MatrixC, MatrixA, MatrixB> closure(C, A, B);
        meta_loop<MatrixA::row>(closure);
    }
};
template <typename MatrixC, typename MatrixA, typename MatrixB>
inline void calc_multy_matrix(MatrixC& C, const MatrixA& A, const MatrixB& B) {
    CALC_A<MatrixC, MatrixA, MatrixB> closure(C, A, B);
    meta_loop<MatrixC::row>(closure);
}

/*!
 *  Умножение матриц на вектор
 */

template <size_t I, typename Matrix, typename Vector>
struct vector_prod_closure
{
    typedef typename Matrix::value_type value_type;
private:
    const Matrix& A;
    const Vector& B;
public:
    vector_prod_closure(const Matrix& A, const Vector& B) : A(A), B(B) {}
    template<size_t K>
    value_type value() const {
        return A(I, K) * B[K];
    }
};

template <typename Matrix, typename Vector>
struct A_Vector_1
{
private:
    const Matrix& A;
    const Vector& B;
    Vector& C;

public:
    A_Vector_1(Vector& C, const Matrix& A, const Vector& B): C(C), A(A), B(B){}

    template <size_t I>
    void apply() {
        vector_prod_closure<I, Matrix, Vector> closure(A, B);
        C[I] = abstract_sums<Matrix::row>(closure);
    }
};

template <typename Matrix, typename Vector>
void calc_multy_vector(Vector& C, const Matrix& A, const Vector& B) {
    A_Vector_1<Matrix, Vector> closure(C, A, B);
    meta_loop<Matrix::row>(closure);
}

/*!
 *  Норма матриц
 */

template <size_t I, typename Matrix>
struct norm_closure
{
    typedef typename Matrix::value_type value_type;
private:
    const Matrix& A;
public:
    norm_closure(const Matrix& A) : A(A) {}
    template<size_t J>
    value_type value() const {
        return std::abs(A(I, J));
    }
};

template <typename Matrix>
struct norm_closure_1
{
private:
    const Matrix& A;
    typedef typename Matrix::value_type value_type;
public:
    value_type C;
public:
    norm_closure_1(const Matrix& A): A(A) {}

    template <size_t I>
    void apply() {
        norm_closure<I, Matrix> closure(A);
        value_type tmp;
        tmp = abstract_sums<Matrix::row>(closure);
        if (C < tmp) C = tmp;
    }
};

template <typename Matrix>
typename Matrix::value_type calc_norma(const Matrix& A) {
    norm_closure_1<Matrix> closure(A);
    meta_loop<Matrix::row>(closure);
    return closure.C;
}


/*!
 *  Обратная матрица
 */

template <size_t i, typename Matrix, typename MatrixA>
struct matrix_2N
{
    typedef typename Matrix::value_type value_type;
private:
    const Matrix& A;
    MatrixA& C;

public:
    matrix_2N(MatrixA& C, const Matrix& A): C(C), A(A) {}

    template <size_t j>
    void apply() {
        if constexpr (j < Matrix::row) {
            C(i, j) = A(i, j);
        }
        else if constexpr (i == (j - Matrix::row))
            C(i, j) = value_type(1);
        else
            C(i, j) = value_type(0);
    }
};
template <typename Matrix, typename MatrixA>
struct matrix_2N_i
{
private:
    const Matrix& A;
    MatrixA& C;
public:
    matrix_2N_i(MatrixA& C, const Matrix& A): C(C), A(A) {}

    template <size_t i>
    void apply() {
        matrix_2N<i, Matrix, MatrixA> closure(C, A);
        meta_loop<MatrixA::col>(closure);
    }
};
template <typename Matrix, typename MatrixA>
void calc_multy_2N(MatrixA& C, const Matrix& A) {
    matrix_2N_i<Matrix, MatrixA> closure(C, A);
    meta_loop<MatrixA::row>(closure);
}

template <size_t i, size_t j, typename Matrix>
struct matrix_invers_1
{
    typedef typename Matrix::value_type value_type;
private:
    Matrix& C;
    value_type ratio;
public:
    matrix_invers_1(Matrix& C, value_type ratio_) : C(C), ratio(ratio_)  {}

    template <size_t k>
    void apply() {
        C(j, k) -= ratio * C(i, k);
    }
};

template <size_t i, typename Matrix>
struct matrix_invers_2
{
    typedef typename Matrix::value_type value_type;
private:
    Matrix& C;

public:
    matrix_invers_2(Matrix& C): C(C) {}

    template <size_t j>
    void apply() {
        if constexpr (i != j) {
            value_type ratio = C(j, i)/C(i, i);
            matrix_invers_1<i, j, Matrix> closure(C, ratio);
            meta_loop<Matrix::col>(closure);
        }
    }
};

template <typename Matrix>
struct matrix_invers
{
    typedef typename Matrix::value_type value_type;
private:

    Matrix& C;
public:
    matrix_invers(Matrix& C): C(C) {}

    template <size_t i>
    void apply() {
        matrix_invers_2<i, Matrix> closure(C);
        meta_loop<Matrix::row>(closure);
    }
};
template <typename Matrix>
void calc_invers_tmp_matrix(Matrix& C) {
    matrix_invers<Matrix> closure(C);
    meta_loop<Matrix::row>(closure);
}
////////////////==============

template <size_t i, typename Matrix>
struct matrix_invers_b
{
    typedef typename Matrix::value_type value_type;
private:
    Matrix& C;
    value_type tmp;
public:
    matrix_invers_b(Matrix& C, value_type tmp_) : C(C), tmp(tmp_) {}

    template <size_t j>
    void apply() {
        C(i, j) /= tmp;
    }
};
template <typename Matrix>
struct matrix_invers_c
{
private:
    Matrix& C;
public:
    matrix_invers_c(Matrix& C) : C(C) {}

    template <size_t i>
    void apply() {
        matrix_invers_b<i, Matrix> closure(C, C(i, i));
        meta_loop<Matrix::col>(closure);
    }
};
template <typename Matrix>
void calc_invers_matrix(Matrix& C) {
    matrix_invers_c<Matrix> closure(C);
    meta_loop<Matrix::row>(closure);
}
//////////Result invers Matrix
///

template <size_t i, typename Matrix, typename MatrixA>
struct matrix_toN
{
    typedef typename Matrix::value_type value_type;
private:
    const MatrixA& TMP;
    Matrix& RES;

public:
    matrix_toN(Matrix& RES_, const MatrixA& TMP_): RES(RES_), TMP(TMP_) {}

    template <size_t j>
    void apply() {
        if constexpr (j < Matrix::row) {
            RES(i, j) = TMP(i, j + Matrix::row);
        }
    }
};

template <typename Matrix, typename MatrixA>
struct matrix_toN_i
{
private:
    const MatrixA& TMP;
    Matrix& RES;
public:
    matrix_toN_i(Matrix& RES_, const MatrixA& TMP_): RES(RES_), TMP(TMP_) {}

    template <size_t i>
    void apply() {
        matrix_toN<i, Matrix, MatrixA> closure(RES, TMP);
        meta_loop<MatrixA::row>(closure);
    }
};
template <typename Matrix, typename MatrixA>
void calc_multy_toN(Matrix& RES, const MatrixA& TMP) {
    matrix_toN_i<Matrix, MatrixA> closure(RES, TMP);
    meta_loop<MatrixA::row>(closure);
}

template <typename Matrix, typename MatrixA>
void Result_invers_Matrix(Matrix& TMP, MatrixA& RESULT, const MatrixA& ORIG) {
    calc_multy_2N(TMP, ORIG);
    calc_invers_tmp_matrix(TMP);
    calc_invers_matrix(TMP);
    calc_multy_toN(RESULT, TMP);
}
//// close invers matrix
///
///
///three diag matrix

template <size_t i, typename Matrix>
struct three_diag_matrix_1
{
    typedef typename Matrix::value_type value_type;
private:
    Matrix& RESULT;
    value_type alpha;
public:
    three_diag_matrix_1(Matrix& RESULT_, value_type alpha_) : RESULT(RESULT_), alpha(alpha_)  {}

    template <size_t j>
    void apply() {
        if constexpr(i == j) RESULT(i, j) = value_type(2)*(value_type(1) + alpha);
        else if constexpr(i == j + 1) {
            RESULT(i, j) = - alpha;
            RESULT(j, i) = RESULT(i, j);
        }
        else RESULT(i, j) = value_type(0);
    }
};
template <typename Matrix>
struct three_diag_matrix_2
{
    typedef typename Matrix::value_type value_type;
private:
    Matrix& RESULT;
    value_type alpha;
public:
    three_diag_matrix_2(Matrix& RESULT_, value_type alpha_) : RESULT(RESULT_), alpha(alpha_)    {}

    template <size_t i>
    void apply() {
        three_diag_matrix_1<i, Matrix> closure(RESULT, alpha);
        meta_loop<Matrix::row>(closure);
    }
};
template <typename Matrix>
void calc_three_diag_matrix(Matrix& RESULT, typename Matrix::value_type alpha) {
    three_diag_matrix_2<Matrix> closure(RESULT, alpha);
    meta_loop<Matrix::col>(closure);
}

//#endif // MATRIX_HPP
