//
// Created by tommaso on 25/03/19.
//

#ifndef OPENFPM_PDATA_SUPPORTBUILDER_HPP
#define OPENFPM_PDATA_SUPPORTBUILDER_HPP

// This is to automatically get the right support (set of particles) for doing DCPSE on a given particle.
// todo: This could be enhanced to support an algorithm for choosing the support in a smart way (to lower condition num)

#include <Space/Shape/Point.hpp>
#include <Vector/vector_dist.hpp>
#include "Support.hpp"

template<unsigned int dim, typename T, typename Prop>
class SupportBuilder
{
private:
    vector_dist<dim, T, Prop> &domain;
    CellList<dim, T, Mem_fast<HeapMemory, local_index_>> cellList;
    const Point<dim, unsigned int> differentialSignature;

public:
    SupportBuilder(vector_dist<dim, T, Prop> &domain, Point<dim, unsigned int> differentialSignature, T rCut);

    SupportBuilder(vector_dist<dim, T, Prop> &domain, unsigned int differentialSignature[dim], T rCut);

    Support<dim, T, Prop> getSupport(vector_dist_iterator itPoint, unsigned int requiredSize);

private:
    size_t getCellLinId(const grid_key_dx<dim> &cellKey);

    size_t getNumElementsInCell(const grid_key_dx<dim> &cellKey);

    size_t getNumElementsInSetOfCells(const std::set<grid_key_dx<dim>> &set);

    void enlargeSetOfCellsUntilSize(std::set<grid_key_dx<dim>> &set, unsigned int requiredSize);

    std::vector<size_t> getPointsInSetOfCells(std::set<grid_key_dx<dim>> set);

    bool isCellKeyInBounds(grid_key_dx<dim> key);
};

// Method definitions below

template<unsigned int dim, typename T, typename Prop>
SupportBuilder<dim, T, Prop>::SupportBuilder(vector_dist<dim, T, Prop> &domain, const Point<dim, unsigned int> differentialSignature,
                               T rCut) : domain(domain), differentialSignature(differentialSignature)
{
    cellList = domain.template getCellList<CellList<dim, T, Mem_fast<HeapMemory, local_index_>>>(rCut);
}

template<unsigned int dim, typename T, typename Prop>
Support<dim, T, Prop> SupportBuilder<dim, T, Prop>::getSupport(vector_dist_iterator itPoint, unsigned int requiredSize)
{
    // Get spatial position from point iterator
    vect_dist_key_dx p = itPoint.get();
    Point<dim, T> pos = domain.getPos(p.getKey());

    // Get cell containing current point and add it to the set of cell keys
    grid_key_dx<dim> curCellKey = cellList.getCellGrid(pos); // Here get the key of the cell where the current point is
    std::set<grid_key_dx<dim>> supportCells;
    supportCells.insert(curCellKey);

    // Make sure to consider a set of cells providing enough points for the support
    enlargeSetOfCellsUntilSize(supportCells, requiredSize + 1); // NOTE: this +1 is because we then remove the point itself

    // Now return all the points from the support into a vector
    std::vector<size_t> supportKeys = getPointsInSetOfCells(supportCells);
    std::remove(supportKeys.begin(), supportKeys.end(), p.getKey());
    return Support<dim, T, Prop>(domain, p.getKey(), supportKeys);
}

template<unsigned int dim, typename T, typename Prop>
size_t SupportBuilder<dim, T, Prop>::getNumElementsInCell(const grid_key_dx<dim> &cellKey)
{
    const size_t curCellId = getCellLinId(cellKey);
    size_t numElements = cellList.getNelements(curCellId);
    return numElements;
}

template<unsigned int dim, typename T, typename Prop>
size_t SupportBuilder<dim, T, Prop>::getNumElementsInSetOfCells(const std::set<grid_key_dx<dim>> &set)
{
    size_t tot = 0;
    for (const auto cell : set)
    {
        tot += getNumElementsInCell(cell);
    }
    return tot;
}

template<unsigned int dim, typename T, typename Prop>
void SupportBuilder<dim, T, Prop>::enlargeSetOfCellsUntilSize(std::set<grid_key_dx<dim>> &set, unsigned int requiredSize)
{
    while (getNumElementsInSetOfCells(set) < requiredSize)
    {
        auto tmpSet = set;
        for (const auto el : tmpSet)
        {
            for (unsigned int i = 0; i < dim; ++i)
            {
                const auto pOneEl = el.move(i, +1);
                const auto mOneEl = el.move(i, -1);
                if (isCellKeyInBounds(pOneEl))
                {
                    set.insert(pOneEl);
                }
                if (isCellKeyInBounds(mOneEl))
                {
                    set.insert(mOneEl);
                }
            }
        }
    }
}

template<unsigned int dim, typename T, typename Prop>
size_t SupportBuilder<dim, T, Prop>::getCellLinId(const grid_key_dx<dim> &cellKey)
{
    mem_id id = cellList.getGrid().LinId(cellKey);
    return static_cast<size_t>(id);
}

template<unsigned int dim, typename T, typename Prop>
std::vector<size_t> SupportBuilder<dim, T, Prop>::getPointsInSetOfCells(std::set<grid_key_dx<dim>> set)
{
    std::vector<size_t> points;
    for (const auto cellKey : set)
    {
        const size_t cellLinId = getCellLinId(cellKey);
        const size_t elemsInCell = getNumElementsInCell(cellKey);
        for (size_t k = 0; k < elemsInCell; ++k)
        {
            size_t el = cellList.get(cellLinId, k);
//            Point<dim, T> pos = domain.getPos(el);
            points.push_back(el);
        }
    }
    return points;
}

template<unsigned int dim, typename T, typename Prop>
SupportBuilder<dim, T, Prop>::SupportBuilder(vector_dist<dim, T, Prop> &domain, unsigned int *differentialSignature, T rCut)
        : SupportBuilder(domain, Point<dim, unsigned int>(differentialSignature), rCut) {}

template<unsigned int dim, typename T, typename Prop>
bool SupportBuilder<dim, T, Prop>::isCellKeyInBounds(grid_key_dx<dim> key)
{
    const size_t *cellGridSize = cellList.getGrid().getSize();
    for (size_t i = 0; i < dim; ++i)
    {
        if (key.value(i) < 0 || key.value(i) >= cellGridSize[i])
        {
            return false;
        }
    }
    return true;
}

#endif //OPENFPM_PDATA_SUPPORTBUILDER_HPP
