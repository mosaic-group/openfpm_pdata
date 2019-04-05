//
// Created by tommaso on 28/03/19.
//

#ifndef OPENFPM_PDATA_DCPSERHS_HPP
#define OPENFPM_PDATA_DCPSERHS_HPP

#include "MonomialBasis.hpp"

template<unsigned int dim>
class DcpseRhs
{
private:
    const Point<dim, unsigned int> differentialSignature;
    const MonomialBasis<dim> derivatives;
    int sign;
public:
    DcpseRhs(const MonomialBasis<dim> &monomialBasis, const Point<dim, unsigned int> &differentialSignature);

    template<typename T, typename MatrixType>
    MatrixType &getVector(MatrixType &b);
};

// Definitions below

template<unsigned int dim>
DcpseRhs<dim>::DcpseRhs(const MonomialBasis<dim> &monomialBasis,
                        const Point<dim, unsigned int> &differentialSignature)
        : differentialSignature(differentialSignature), derivatives(monomialBasis.getDerivative(differentialSignature))
{
    unsigned int order = (Monomial<dim>(differentialSignature)).order();
    if (order % 2 == 0)
    {
        sign = 1;
    } else
    {
        sign = -1;
    }
}

template<unsigned int dim>
template<typename T, typename MatrixType>
MatrixType &DcpseRhs<dim>::getVector(MatrixType &b)
{
    // The given vector V should have the right dimensions
    assert(b.cols() == 1);
    assert(b.rows() == derivatives.size());
    for (unsigned int i = 0; i < derivatives.size(); ++i)
    {
        const Monomial<dim> dm = derivatives.getElement(i);
        b(i, 0) = sign * dm.evaluate(Point<dim, T>(0));
    }
    return b;
}

#endif //OPENFPM_PDATA_DCPSERHS_HPP
