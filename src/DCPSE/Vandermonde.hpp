//
// Created by tommaso on 21/03/19.
//

#ifndef OPENFPM_PDATA_VANDERMONDE_HPP
#define OPENFPM_PDATA_VANDERMONDE_HPP

#include "MonomialBasis.hpp"
#include "VandermondeRowBuilder.hpp"
#include "Support.hpp"

template<unsigned int dim, typename T, typename MatrixType>
class Vandermonde
{
private:
    const Point<dim, T> point;
    const std::vector<Point<dim, T>> offsets;
    const MonomialBasis<dim> monomialBasis;
    T eps;

public:
    Vandermonde(const Point<dim, T> &point, const std::vector<Point<dim, T>> &neighbours,
                const MonomialBasis<dim> &monomialBasis);

    template<typename Prop>
    Vandermonde(const Support<dim, T, Prop> &support,
                const MonomialBasis<dim> &monomialBasis);

    MatrixType &getMatrix(MatrixType &M);

    T getEps();

private:
    static std::vector<Point<dim, T>>
    computeOffsets(const Point<dim, T> &point, const std::vector<Point<dim, T>> &neighbours);

    void computeEps(T factor);

    static T computeAbsSum(const Point<dim, T> &x);

    void initialize();
};

template<unsigned int dim, typename T, typename MatrixType>
void Vandermonde<dim, T, MatrixType>::initialize()
{
    // First check that the number of points given is enough for building the Vandermonde matrix
    if (offsets.size() < monomialBasis.size())
    {
        ACTION_ON_ERROR(std::length_error("Not enough neighbour points passed for Vandermonde matrix construction!"));
    }
    // Compute eps for this point
    computeEps(2);
}

template<unsigned int dim, typename T, typename MatrixType>
Vandermonde<dim, T, MatrixType>::Vandermonde(const Point<dim, T> &point, const std::vector<Point<dim, T>> &neighbours,
                                             const MonomialBasis<dim> &monomialBasis)
        : point(point),
          offsets(computeOffsets(point, neighbours)),
          monomialBasis(monomialBasis)
{
    initialize();
}

template<unsigned int dim, typename T, typename MatrixType>
template<typename Prop>
Vandermonde<dim, T, MatrixType>::Vandermonde(const Support<dim, T, Prop> &support,
                                             const MonomialBasis<dim> &monomialBasis)
        : point(support.getReferencePoint()),
          offsets(support.getOffsets()),
          monomialBasis(monomialBasis)
{
    initialize();
}


template<unsigned int dim, typename T, typename MatrixType>
std::vector<Point<dim, T>>
Vandermonde<dim, T, MatrixType>::computeOffsets(const Point<dim, T> &point,
                                                const std::vector<Point<dim, T>> &neighbours)
{
    std::vector<Point<dim, T>> offsets;
    for (auto &other : neighbours)
    {
        Point<dim, T> curOffset(point);
        curOffset -= other;
        offsets.push_back(curOffset);
    }
    return offsets;
}

template<unsigned int dim, typename T, typename MatrixType>
MatrixType &Vandermonde<dim, T, MatrixType>::getMatrix(MatrixType &M)
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
    assert(eps != 0);
}

template<unsigned int dim, typename T, typename MatrixType>
T Vandermonde<dim, T, MatrixType>::computeAbsSum(const Point<dim, T> &x)
{
    T absSum = 0;
    for (unsigned int i = 0; i < dim; ++i)
    {
        absSum += fabs(x.value(i));
    }
    return absSum;
}

template<unsigned int dim, typename T, typename MatrixType>
T Vandermonde<dim, T, MatrixType>::getEps()
{
    return eps;
}

#endif //OPENFPM_PDATA_VANDERMONDE_HPP
