//
// Created by tommaso on 21/03/19.
//

#ifndef OPENFPM_PDATA_VANDERMONDE_HPP
#define OPENFPM_PDATA_VANDERMONDE_HPP

#include "../../openfpm_numerics/src/DMatrix/EMatrix.hpp"
#include "MonomialBasis.hpp"
#include "VandermondeRowBuilder.hpp"

template<unsigned int dim, typename T, typename MatrixType>
class Vandermonde
{
private:
    const Point<dim, T> point;
    std::vector<Point<dim, T>> offsets;
    const MonomialBasis<dim> monomialBasis;
    T eps;

public:
    Vandermonde(const Point<dim, T> &point, const std::vector<Point<dim, T>> &neighbours,
                const MonomialBasis<dim> &monomialBasis);

    MatrixType & getMatrix(MatrixType &M);

private:
    void computeOffsets(const std::vector<Point<dim, T>> &neighbours);

    void computeEps(T factor);

    static T computeAbsSum(Point<dim, T> &x);
};

template<unsigned int dim, typename T, typename MatrixType>
Vandermonde<dim, T, MatrixType>::Vandermonde(const Point<dim, T> &point, const std::vector<Point<dim, T>> &neighbours,
                                 const MonomialBasis<dim> &monomialBasis) : point(point),
                                                                            monomialBasis(monomialBasis)
{
    // Compute the offsets from point to all neighbours (and store them)
    computeOffsets(neighbours);
    // Compute eps for this point
    computeEps(2);
}

template<unsigned int dim, typename T, typename MatrixType>
void Vandermonde<dim, T, MatrixType>::computeOffsets(const std::vector<Point<dim, T>> &neighbours)
{
    for (auto &other : neighbours)
    {
        Point<dim, T> curOffset(point);
        curOffset -= other;
        offsets.push_back(curOffset);
    }
}

template<unsigned int dim, typename T, typename MatrixType>
MatrixType & Vandermonde<dim, T, MatrixType>::getMatrix(MatrixType &M)
{
    // Build the Vandermonde matrix, row-by-row
    VandermondeRowBuilder<dim, T> vrb(monomialBasis);
    unsigned int row = 0;
    for (auto &offset : offsets)
    {
        vrb.buildRow(M, row, offset, eps);
        ++row;
    }
    return M;
}

template<unsigned int dim, typename T, typename MatrixType>
void Vandermonde<dim, T, MatrixType>::computeEps(T factor)
{
    T avgNeighbourSpacing = 0;
    for (auto &offset : offsets)
    {
        avgNeighbourSpacing += computeAbsSum(offset);
    }
    avgNeighbourSpacing /= offsets.size();
    eps = factor * avgNeighbourSpacing;
    std::cout << "eps=" << eps << std::endl;
}

template<unsigned int dim, typename T, typename MatrixType>
T Vandermonde<dim, T, MatrixType>::computeAbsSum(Point<dim, T> &x)
{
    T absSum = 0;
    for (const auto elem : x.asArray())
    {
        absSum += abs(elem);
    }
    return absSum;
}

#endif //OPENFPM_PDATA_VANDERMONDE_HPP
