//
// Created by tommaso on 21/03/19.
//

#ifndef OPENFPM_PDATA_MONOMIALBASISELEMENT_H
#define OPENFPM_PDATA_MONOMIALBASISELEMENT_H


#include <vector>
#include <ostream>
#include "Space/Shape/Point.hpp"

template<unsigned int dim>
class Monomial
{
//    The base idea is that this should behave like a pair <unsigned int, std::vector<unsigned int>>
private:
    unsigned int sum = 0;
    Point<dim, unsigned int> exponents;

public:
    Monomial();

    explicit Monomial(const Point<dim, unsigned int> &other);

    explicit Monomial(const Point<dim, long int> &other);

    explicit Monomial(const unsigned int other[dim]);

    Monomial(const Monomial<dim> &other) : Monomial(other.exponents) {}

    Monomial<dim> &operator=(const Monomial<dim> &other);

    bool operator==(const Monomial<dim> &other) const;

    unsigned int order() const;

    unsigned int getExponent(unsigned int i) const;

    void setExponent(unsigned int i, unsigned int value);

    template<typename T>
    T evaluate(const Point<dim, T> x) const;

    template<typename T>
    T evaluate(const T x[dim]) const;

    Monomial<dim> getDerivative(const Point<dim, unsigned int> differentialOrder) const;

    template<typename charT, typename traits>
    friend std::basic_ostream<charT, traits> &
    operator<<(std::basic_ostream<charT, traits> &lhs, Monomial<dim> const &rhs)
    {
        return lhs << rhs.exponents.toString();
    }

private:
    void updateSum();
};

////// Definitions below

template<unsigned int dim>
Monomial<dim>::Monomial()
{
    exponents.zero();
    sum = 0;
}

template<unsigned int dim>
Monomial<dim>::Monomial(const Point<dim, unsigned int> &other)
        : exponents(other)
{
    updateSum();
}

template<unsigned int dim>
Monomial<dim>::Monomial(const Point<dim, long int> &other)
{
    for (size_t i = 0; i < other.nvals; ++i)
    {
        exponents.get(i) = other.value(i);
    }
    updateSum();
}

template<unsigned int dim>
Monomial<dim>::Monomial(const unsigned int other[dim])
        : Monomial(Point<3, unsigned int>(other)) {}

template<unsigned int dim>
Monomial<dim> &Monomial<dim>::operator=(const Monomial<dim> &other)
{
    exponents = other.exponents;
    sum = other.sum;
    return *this;
}

template<unsigned int dim>
void Monomial<dim>::updateSum()
{
    unsigned int partialSum = 0;
    for (unsigned int i = 0; i < dim; ++i)
    {
        partialSum += exponents.value(i);
    }
    sum = partialSum;
}

template<unsigned int dim>
unsigned int Monomial<dim>::order() const
{
    return sum;
}

template<unsigned int dim>
unsigned int Monomial<dim>::getExponent(unsigned int i) const
{
    return exponents.value(i);
}

template<unsigned int dim>
void Monomial<dim>::setExponent(unsigned int i, unsigned int value)
{
    exponents.get(i) = value;
    updateSum();
}

template<unsigned int dim>
bool Monomial<dim>::operator==
        (const Monomial<dim> &other) const
{
    return exponents == other.exponents;
}

template<unsigned int dim>
template<typename T>
T Monomial<dim>::evaluate(const Point<dim, T> x) const
{
    T res = 1;
    for (unsigned int i = 0; i < dim; ++i)
    {
        res *= pow(x.value(i), getExponent(i));
    }
    return res;
}

template<unsigned int dim>
Monomial<dim> Monomial<dim>::getDerivative(const Point<dim, unsigned int> differentialOrder) const
{
    Point<dim, unsigned int> res(exponents);
    for (unsigned int i = 0; i < dim; ++i)
    {
        res.get(i) = static_cast<unsigned int>(
                std::max(static_cast<int>(res.value(i)) - static_cast<int>(differentialOrder.value(i)), 0)
        );
    }
    return Monomial(res);
}

template<unsigned int dim>
template<typename T>
T Monomial<dim>::evaluate(const T x[dim]) const
{
    return evaluate(Point<dim, T>(x));
}


#endif //OPENFPM_PDATA_MONOMIALBASISELEMENT_H
