//
// Created by tommaso on 21/03/19.
//

#ifndef OPENFPM_PDATA_MONOMIALBASISELEMENT_H
#define OPENFPM_PDATA_MONOMIALBASISELEMENT_H


#include <vector>
#include <ostream>
#include "Space/Shape/Point.hpp"

template<unsigned int dim>
class MonomialBasisElement
{
//    The base idea is that this should behave like a pair <unsigned int, std::vector<unsigned int>>
private:
    unsigned int sum = 0;
    Point<dim, unsigned int> exponents;

public:
    MonomialBasisElement();

    explicit MonomialBasisElement(const Point<dim, unsigned int> &other);

    explicit MonomialBasisElement(const Point<dim, long int> &other);

    explicit MonomialBasisElement(const unsigned int other[dim]);

    MonomialBasisElement(const MonomialBasisElement<dim> &other) : MonomialBasisElement(other.exponents) {}

    MonomialBasisElement<dim> &operator=(const MonomialBasisElement<dim> &other);

    bool operator==(const MonomialBasisElement<dim> &other) const;

    unsigned int order();

    unsigned int getExponent(unsigned int i);

    void setExponent(unsigned int i, unsigned int value);

    template<typename charT, typename traits>
    friend std::basic_ostream<charT, traits> &
    operator<<(std::basic_ostream<charT, traits> &lhs, MonomialBasisElement<dim> const &rhs)
    {
        return lhs << rhs.exponents.toString();
    }

private:
    void updateSum();
};

////// Definitions below

template<unsigned int dim>
MonomialBasisElement<dim>::MonomialBasisElement()
{
    exponents.zero();
    sum = 0;
}

template<unsigned int dim>
MonomialBasisElement<dim>::MonomialBasisElement(const Point<dim, unsigned int> &other)
        : exponents(other)
{
    updateSum();
}

template<unsigned int dim>
MonomialBasisElement<dim>::MonomialBasisElement(const Point<dim, long int> &other)
{
    for (size_t i = 0; i < other.nvals; ++i)
    {
        exponents.get(i) = other.value(i);
    }
    updateSum();
}

template<unsigned int dim>
MonomialBasisElement<dim>::MonomialBasisElement(const unsigned int other[dim])
        : MonomialBasisElement(Point<3, unsigned int>(other)) {}

template<unsigned int dim>
MonomialBasisElement<dim> &MonomialBasisElement<dim>::operator=(const MonomialBasisElement<dim> &other)
{
    exponents = other.exponents;
    sum = other.sum;
    return *this;
}

template<unsigned int dim>
void MonomialBasisElement<dim>::updateSum()
{
    unsigned int partialSum = 0;
    for (unsigned int i = 0; i < dim; ++i)
    {
        partialSum += exponents.value(i);
    }
    sum = partialSum;
}

template<unsigned int dim>
unsigned int MonomialBasisElement<dim>::order()
{
    return sum;
}

template<unsigned int dim>
unsigned int MonomialBasisElement<dim>::getExponent(unsigned int i)
{
    return exponents.value(i);
}

template<unsigned int dim>
void MonomialBasisElement<dim>::setExponent(unsigned int i, unsigned int value)
{
    exponents.get(i) = value;
    updateSum();
}

template<unsigned int dim>
bool MonomialBasisElement<dim>::operator==
        (const MonomialBasisElement<dim> &other) const
{
    return exponents == other.exponents;
}


#endif //OPENFPM_PDATA_MONOMIALBASISELEMENT_H
