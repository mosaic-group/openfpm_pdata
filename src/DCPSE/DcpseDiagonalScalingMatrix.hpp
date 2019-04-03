//
// Created by tommaso on 29/03/19.
//

#ifndef OPENFPM_PDATA_DCPSEDIAGONALSCALINGMATRIX_HPP
#define OPENFPM_PDATA_DCPSEDIAGONALSCALINGMATRIX_HPP

#include "MonomialBasis.hpp"
#include "Support.hpp"

template <unsigned int dim>
class DcpseDiagonalScalingMatrix
{
private:
    const MonomialBasis<dim> monomialBasis;

public:
    DcpseDiagonalScalingMatrix(const MonomialBasis<dim> &monomialBasis) : monomialBasis(monomialBasis) {}

    template <typename T, typename MatrixType, typename Prop>
    void buildMatrix(MatrixType &M, Support<dim, T, Prop> support, T eps);
};

template<unsigned int dim>
template<typename T, typename MatrixType, typename Prop>
void DcpseDiagonalScalingMatrix<dim>::buildMatrix(MatrixType &M, Support<dim, T, Prop> support, T eps)
{
    // Check that all the dimension constraints are met
    assert(support.size() > monomialBasis.size());
    assert(M.rows() == support.size());
    assert(M.cols() == support.size());

    // Fill the diagonal matrix
    M.setZero(); // Make sure the rest of the matrix is zero!
    int i = 0;
    for (const auto& pt : support.getOffsets())
    {
        M(i,i) = exp(- norm2(pt) / (2.0 * eps * eps));
        ++i;
    }
}

#endif //OPENFPM_PDATA_DCPSEDIAGONALSCALINGMATRIX_HPP
