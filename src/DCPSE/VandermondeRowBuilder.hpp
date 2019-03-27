//
// Created by tommaso on 22/03/19.
//

#ifndef OPENFPM_PDATA_VANDERMONDEROW_HPP
#define OPENFPM_PDATA_VANDERMONDEROW_HPP

#include "../../openfpm_numerics/src/DMatrix/EMatrix.hpp"
#include "MonomialBasis.hpp"

template <unsigned int dim, typename T>
class VandermondeRowBuilder
{
private:
    MonomialBasis<dim> monomialBasis;

public:
    VandermondeRowBuilder(const MonomialBasis<dim> &monomialBasis) : monomialBasis(monomialBasis) {}

    void buildRow(EMatrix<T, Eigen::Dynamic, Eigen::Dynamic> &M, unsigned int row, Point<dim, T> x, T eps);

private:
    T computePower(Monomial<dim> mbe, Point<dim, T> x);
};

template<unsigned int dim, typename T>
void VandermondeRowBuilder<dim, T>::buildRow(EMatrix<T, Eigen::Dynamic, Eigen::Dynamic> &M, unsigned int row, Point<dim, T> x, T eps)
{
    unsigned int col = 0;
    for (auto& basisElement : monomialBasis.getElements())
    {
        auto & mbe = monomialBasis.getElement(col);
        M(row, col) = computePower(mbe, x);
        M(row, col) /= pow(eps, mbe.order());
        ++col;
    }
}

template<unsigned int dim, typename T>
T VandermondeRowBuilder<dim, T>::computePower(Monomial<dim> mbe, Point<dim, T> x)
{
    T res = 1;
    for (int i = 0; i < dim; ++i)
    {
        res *= pow(x.value(i), mbe.getExponent(i));
    }
    return res;
}

#endif //OPENFPM_PDATA_VANDERMONDEROW_HPP
