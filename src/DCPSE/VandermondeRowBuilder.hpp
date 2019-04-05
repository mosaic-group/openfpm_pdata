//
// Created by tommaso on 22/03/19.
//

#ifndef OPENFPM_PDATA_VANDERMONDEROW_HPP
#define OPENFPM_PDATA_VANDERMONDEROW_HPP

#include "MonomialBasis.hpp"

template <unsigned int dim, typename T>
class VandermondeRowBuilder
{
private:
    const MonomialBasis<dim> monomialBasis;

public:
    VandermondeRowBuilder(const MonomialBasis<dim> &monomialBasis) : monomialBasis(monomialBasis) {}

    template <typename MatrixType>
    void buildRow(MatrixType &M, unsigned int row, Point<dim, T> x, T eps);
};

template<unsigned int dim, typename T>
template <typename MatrixType>
void VandermondeRowBuilder<dim, T>::buildRow(MatrixType &M, unsigned int row, Point<dim, T> x, T eps)
{
    unsigned int col = 0;
    for (auto& basisElement : monomialBasis.getElements())
    {
        Monomial<dim> m = monomialBasis.getElement(col);
        M(row, col) = m.evaluate(x);
        M(row, col) /= pow(eps, m.order());
        ++col;
    }
}

#endif //OPENFPM_PDATA_VANDERMONDEROW_HPP
