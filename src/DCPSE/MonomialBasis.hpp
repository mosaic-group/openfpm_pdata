//
// Created by tommaso on 20/03/19.
//

#ifndef OPENFPM_PDATA_MONOMIALBASIS_H
#define OPENFPM_PDATA_MONOMIALBASIS_H

#include <vector>
#include <Grid/grid_sm.hpp>
#include <Grid/iterators/grid_key_dx_iterator_sub_bc.hpp>
#include "Monomial.hpp"

template<unsigned int dim>
class MonomialBasis
{
private:
    std::vector<Monomial<dim>> basis;

public:
    MonomialBasis(const std::vector<unsigned int> &degrees, unsigned int convergenceOrder);

    MonomialBasis(unsigned int degrees[dim], unsigned int convergenceOrder);

//    explicit MonomialBasis(Point<dim, unsigned int> degrees, unsigned int convergenceOrder);

    explicit MonomialBasis(const std::vector<Monomial<dim>> &basis) : basis(basis) {}

    MonomialBasis(const MonomialBasis &other);

    MonomialBasis &operator=(const MonomialBasis &other);

    unsigned int size() const;

    const Monomial<dim> &getElement(unsigned int i) const;

    Monomial<dim> &getElement(unsigned int i);

    const std::vector<Monomial<dim>> &getElements() const;

    MonomialBasis<dim> getDerivative(Point<dim, unsigned int> differentialOrder) const;

    bool operator==(const MonomialBasis &other) const;

    template<typename charT, typename traits>
    friend std::basic_ostream<charT, traits> &
    operator<<(std::basic_ostream<charT, traits> &lhs, MonomialBasis<dim> const &rhs)
    {
        lhs << "MonomialBasis: size=" << rhs.size() << ", elements={ ";
        for (const auto &el : rhs.getElements())
        {
            lhs << "(" << el << ") ";
        }
        lhs << "}" << std::endl;
        return lhs;
    }

private:
    void generateBasis(std::vector<unsigned int> m, unsigned int r);
};

//// Definitions below

template<unsigned int dim>
MonomialBasis<dim>::MonomialBasis(const std::vector<unsigned int> &degrees, unsigned int convergenceOrder)
{
    generateBasis(degrees, convergenceOrder);
}

template<unsigned int dim>
MonomialBasis<dim>::MonomialBasis(unsigned int *degrees, unsigned int convergenceOrder)
        : MonomialBasis(std::vector<unsigned int>(degrees, degrees + dim), convergenceOrder) {}

template<unsigned int dim>
MonomialBasis<dim>::MonomialBasis(const MonomialBasis &other)
{
    basis = other.basis; // Here it works because both std::vector and Monomial perform a deep copy.
}

template<unsigned int dim>
MonomialBasis<dim> &MonomialBasis<dim>::operator=(const MonomialBasis &other)
{
    basis = other.basis; // Here it works because both std::vector and Monomial perform a deep copy.
    return *this;
}

template<unsigned int dim>
unsigned int MonomialBasis<dim>::size() const
{
    return basis.size();
}

template<unsigned int dim>
const Monomial<dim> &MonomialBasis<dim>::getElement(unsigned int i) const
{
    return basis[i];
}

template<unsigned int dim>
Monomial<dim> &MonomialBasis<dim>::getElement(unsigned int i)
{
    return basis[i];
}

template<unsigned int dim>
void MonomialBasis<dim>::generateBasis(std::vector<unsigned int> m, unsigned int r)
{
    // Compute the vector of actual dimensions to iterate over
    // NOTE: each index can go up to sum(m)+r
    unsigned int mSum = std::accumulate(m.begin(), m.end(), 0U);
    unsigned int orderLimit = mSum + r;
    size_t dimensions[dim];
    std::fill(dimensions, dimensions + dim, orderLimit);

    // Now initialize grid with appropriate size, then start-stop points and boundary conditions for the iterator
    grid_sm<dim, void> grid(dimensions);

    long int startV[dim] = {}; // 0-initialized
    grid_key_dx<dim, long int> start(startV);
    grid_key_dx<dim, long int> stop(dimensions);

    size_t bc[dim];
    std::fill(bc, bc + dim, NON_PERIODIC);

    grid_key_dx_iterator_sub_bc<dim> it(grid, start, stop, bc);

    // Finally compute alpha_min
//    unsigned char alphaMin = static_cast<unsigned char>(!(mSum % 2)); // if mSum is even, alpha_min must be 1
    unsigned char alphaMin = 0; // we want to always have 1 in the basis

    while (it.isNext())
    {
        Point<dim, long int> p = it.get().get_k();
        Monomial<dim> candidateBasisElement(p);
        // Filter out the elements which don't fullfil the theoretical condition for being in the vandermonde matrix
        if (candidateBasisElement.order() < orderLimit && candidateBasisElement.order() >= alphaMin)
        {
            basis.push_back(candidateBasisElement);
        }
        ++it;
    }
}

template<unsigned int dim>
const std::vector<Monomial<dim>> &MonomialBasis<dim>::getElements() const
{
    return basis;
}

template<unsigned int dim>
MonomialBasis<dim> MonomialBasis<dim>::getDerivative(const Point<dim, unsigned int> differentialOrder) const
{
    std::vector<Monomial<dim>> derivatives;
    for (const auto &monomial : getElements())
    {
        derivatives.push_back(monomial.getDerivative(differentialOrder));
    }
    return MonomialBasis<dim>(derivatives);
}

template<unsigned int dim>
bool MonomialBasis<dim>::operator==(const MonomialBasis &other) const
{
    return basis == other.basis;
}

//template<unsigned int dim>
//MonomialBasis<dim>::MonomialBasis(Point<dim, unsigned int> degrees, unsigned int convergenceOrder)
//        : MonomialBasis(degrees.asArray(), convergenceOrder) {}

#endif //OPENFPM_PDATA_MONOMIALBASIS_H
