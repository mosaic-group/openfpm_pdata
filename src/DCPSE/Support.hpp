//
// Created by tommaso on 29/03/19.
//

#ifndef OPENFPM_PDATA_SUPPORT_HPP
#define OPENFPM_PDATA_SUPPORT_HPP

#include <Space/Shape/Point.hpp>
#include <Vector/vector_dist.hpp>

template<unsigned int dim, typename T, typename Prop>
class Support
{
    // This class is basically a wrapper around a point and a set of offsets.
    // Offsets are stored as they are the data that is mostly required in DCPSE, so we
    // pre-compute and store them, while the positions can be re-computed everytime is
    // necessary (it should almost never be the case) (todo: check if this is required)

private:
    const vector_dist<dim, T, Prop> &domain;
    const size_t referencePointKey;
    const std::vector<size_t> keys;
    const std::vector<Point<dim, T>> offsets;

public:
    Support() {};

    Support(const vector_dist<dim, T, Prop> &domain, const size_t &referencePoint, const std::vector<size_t> &keys)
            : domain(domain),
              referencePointKey(referencePoint),
              keys(keys),
              offsets(computeOffsets(referencePoint, keys)) {}

    Support(const Support<dim, T, Prop> &other);

    size_t size();

    const Point<dim, T> getReferencePoint() const;

    const size_t getReferencePointKey() const;

    const std::vector<size_t> &getKeys() const;

    const std::vector<Point<dim, T>> &getOffsets() const;

private:
    std::vector<Point<dim, T>>
    computeOffsets(const size_t referencePoint, const std::vector<size_t> &keys);
};

template<unsigned int dim, typename T, typename Prop>
std::vector<Point<dim, T>>
Support<dim, T, Prop>::computeOffsets(const size_t referencePoint, const std::vector<size_t> &keys)
{
    std::vector<Point<dim, T>> offsets;
    for (auto &otherK : keys)
    {
        Point<dim, T> curOffset(domain.getPos(referencePoint));
        curOffset -= domain.getPos(otherK);
        offsets.push_back(curOffset);
    }
    return offsets;
}

template<unsigned int dim, typename T, typename Prop>
const Point<dim, T> Support<dim, T, Prop>::getReferencePoint() const
{
    return Point<dim, T>(domain.getPos(referencePointKey));
}

template<unsigned int dim, typename T, typename Prop>
const std::vector<Point<dim, T>> &Support<dim, T, Prop>::getOffsets() const
{
    return offsets;
}

template<unsigned int dim, typename T, typename Prop>
size_t Support<dim, T, Prop>::size()
{
    return offsets.size();
}

template<unsigned int dim, typename T, typename Prop>
Support<dim, T, Prop>::Support(const Support<dim, T, Prop> &other)
        : domain(other.domain),
          referencePointKey(other.referencePointKey),
          keys(other.keys),
          offsets(other.offsets) {}

template<unsigned int dim, typename T, typename Prop>
const size_t Support<dim, T, Prop>::getReferencePointKey() const
{
    return referencePointKey;
}

template<unsigned int dim, typename T, typename Prop>
const std::vector<size_t> &Support<dim, T, Prop>::getKeys() const
{
    return keys;
}

#endif //OPENFPM_PDATA_SUPPORT_HPP
