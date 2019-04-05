//
// Created by tommaso on 29/03/19.
//

#ifndef OPENFPM_PDATA_DCPSE_HPP
#define OPENFPM_PDATA_DCPSE_HPP

#ifdef HAVE_EIGEN

#include "Vector/vector_dist.hpp"
#include "MonomialBasis.hpp"
#include "../../openfpm_numerics/src/DMatrix/EMatrix.hpp"
#include "SupportBuilder.hpp"
#include "Support.hpp"
#include "Vandermonde.hpp"
#include "DcpseDiagonalScalingMatrix.hpp"
#include "DcpseRhs.hpp"

template<unsigned int dim, typename T, typename... list>
class Dcpse
{
    // This works in this way:
    // 1) User constructs this by giving a domain of points (where one of the properties is the value of our f),
    //    the signature of the differential operator and the error order bound.
    // 2) The machinery for assembling and solving the linear system for coefficients starts...
    // 3) The user can then call an evaluate(point) method to get the evaluation of the differential operator
    //    on the given point.
private:
    const Point<dim, unsigned int> differentialSignature;
    const unsigned int differentialOrder;
    const MonomialBasis<dim> monomialBasis;
    std::vector<EMatrix<T, Eigen::Dynamic, 1>> localCoefficients; // Each MPI rank has just access to the local ones
    std::vector<Support<dim, T, aggregate<list...>>> localSupports; // Each MPI rank has just access to the local ones
    std::vector<T> localEps; // Each MPI rank has just access to the local ones

public:
    // Here we require the first element of the aggregate to be:
    // 1) the value of the function f on the point
    Dcpse(vector_dist<dim, T, aggregate<list...>> &particles,
          Point<dim, unsigned int> differentialSignature,
          unsigned int convergenceOrder,
          T rCut,
          T supportSizeFactor = 1);

    /**
     * Computes the value of the differential operator on all the particles,
     * using the f values stored at the fValuePos position in the aggregate
     * and storing the resulting Df values at the DfValuePos position in the aggregate.
     * @tparam fValuePos Position in the aggregate of the f values to use.
     * @tparam DfValuePos Position in the aggregate of the Df values to store.
     * @param particles The set of particles to iterate over.
     */
    template<unsigned int fValuePos, unsigned int DfValuePos>
    void computeDifferentialOperator(vector_dist<dim, T, aggregate<list...>> &particles);

private:
    void initializeAdaptive(vector_dist<dim, T, aggregate<list...>> &particles,
                            unsigned int convergenceOrder,
                            T rCut);

    void initializeStaticSize(vector_dist<dim, T, aggregate<list...>> &particles,
                              unsigned int convergenceOrder,
                              T rCut,
                              T supportSizeFactor);

    T computeKernel(Point<dim, T> x, EMatrix<T, Eigen::Dynamic, 1> a) const;

    T conditionNumber(const EMatrix<T, -1, -1> &V, T condTOL) const;
};

template<unsigned int dim, typename T, typename... list>
Dcpse<dim, T, list...>::Dcpse(vector_dist<dim, T, aggregate<list...>> &particles,
                              Point<dim, unsigned int> differentialSignature,
                              unsigned int convergenceOrder, T rCut, T supportSizeFactor) :
        differentialSignature(differentialSignature),
        differentialOrder(Monomial<dim>(differentialSignature).order()),
        monomialBasis(differentialSignature.asArray(), convergenceOrder)
{
    if (supportSizeFactor < 1)
    {
        initializeAdaptive(particles, convergenceOrder, rCut);
    } else
    {
        initializeStaticSize(particles, convergenceOrder, rCut, supportSizeFactor);
    }
}

template<unsigned int dim, typename T, typename ... list>
T Dcpse<dim, T, list...>::conditionNumber(const EMatrix<T, -1, -1> &V, T condTOL) const
{
    Eigen::JacobiSVD<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>> svd(V);
    T cond = svd.singularValues()(0)
             / svd.singularValues()(svd.singularValues().size() - 1);
    if (cond > condTOL)
    {
        std::cout
                << "WARNING: cond(V) = " << cond
                << " is greater than TOL = " << condTOL
                << ",  numPoints(V) = " << V.rows()
                << std::endl; // debug
    }
    return cond;
}

template<unsigned int dim, typename T, typename ... list>
template<unsigned int fValuePos, unsigned int DfValuePos>
void Dcpse<dim, T, list...>::computeDifferentialOperator(vector_dist<dim, T, aggregate<list...>> &particles)
{
    char sign = 1;
    if (differentialOrder % 2 == 0)
    {
        sign = -1;
    }

    auto it = particles.getDomainIterator();
    auto coefficientsIt = localCoefficients.begin();
    auto supportsIt = localSupports.begin();
    auto epsIt = localEps.begin();
    while (it.isNext())
    {
        double eps = *epsIt;

        T Dfxp = 0;
        Support<dim, T, aggregate<list...>> support = *supportsIt;
        size_t xpK = support.getReferencePointKey();
        Point<dim, T> xp = support.getReferencePoint();
        T fxp = sign * particles.template getProp<fValuePos>(xpK);
        for (auto &xqK : support.getKeys())
        {
            Point<dim, T> xq = particles.getPos(xqK);
            T fxq = particles.template getProp<fValuePos>(xqK);
            Point<dim, T> normalizedArg = (xp - xq) / eps;
            EMatrix<T, Eigen::Dynamic, 1> &a = *coefficientsIt;
            Dfxp += (fxq + fxp) * computeKernel(normalizedArg, a);
        }
        Dfxp /= pow(eps, differentialOrder);
        //
        T trueDfxp = particles.template getProp<2>(xpK);
        // Store Dfxp in the right position
        particles.template getProp<DfValuePos>(xpK) = Dfxp;
        //
        ++it;
        ++coefficientsIt;
        ++supportsIt;
        ++epsIt;
    }
}

template<unsigned int dim, typename T, typename ... list>
T Dcpse<dim, T, list...>::computeKernel(Point<dim, T> x, EMatrix<T, Eigen::Dynamic, 1> a) const
{
    T res = 0;
    unsigned int counter = 0;
    for (const Monomial<dim> &m : monomialBasis.getElements())
    {
        T coeff = a(counter);
        T mbValue = m.evaluate(x);
        T expFactor = exp(-norm2(x));
        res += coeff * mbValue * expFactor;
        ++counter;
    }
    return res;
}

template<unsigned int dim, typename T, typename... list>
void Dcpse<dim, T, list...>::initializeAdaptive(vector_dist<dim, T, aggregate<list...>> &particles,
                                                unsigned int convergenceOrder, T rCut)
{
    SupportBuilder<dim, T, aggregate<list...>>
            supportBuilder(particles, differentialSignature, rCut);
    unsigned int requiredSupportSize = monomialBasis.size();

    auto it = particles.getDomainIterator();
    while (it.isNext())
    {
        const T condVTOL = 1e3;

        // Get the points in the support of the DCPSE kernel and store the support for reuse
        Support<dim, T, aggregate<list...>> support = supportBuilder.getSupport(it, requiredSupportSize);
        EMatrix<T, Eigen::Dynamic, Eigen::Dynamic> V(support.size(), monomialBasis.size());

        // Vandermonde matrix computation
        Vandermonde<dim, T, EMatrix<T, Eigen::Dynamic, Eigen::Dynamic>>
                vandermonde(support, monomialBasis);
        vandermonde.getMatrix(V);

        T condV = conditionNumber(V, condVTOL);
        T eps = vandermonde.getEps();

        if (condV > condVTOL)
        {
            requiredSupportSize *= 2;
            std::cout
                    << "INFO: Increasing, requiredSupportSize = " << requiredSupportSize
                    << std::endl; // debug
            continue;
        } else
        {
            requiredSupportSize = monomialBasis.size();
        }

        localSupports.push_back(support);
        localEps.push_back(eps);
        // Compute the diagonal matrix E
        DcpseDiagonalScalingMatrix<dim> diagonalScalingMatrix(monomialBasis);
        EMatrix<T, Eigen::Dynamic, Eigen::Dynamic> E(support.size(), support.size());
        diagonalScalingMatrix.buildMatrix(E, support, eps);
        // Compute intermediate matrix B
        EMatrix<T, Eigen::Dynamic, Eigen::Dynamic> B = E * V;
        // Compute matrix A
        EMatrix<T, Eigen::Dynamic, Eigen::Dynamic> A = B.transpose() * B;
        // Compute RHS vector b
        DcpseRhs<dim> rhs(monomialBasis, differentialSignature);
        EMatrix<T, Eigen::Dynamic, 1> b(monomialBasis.size(), 1);
        rhs.template getVector<T>(b);
        // Get the vector where to store the coefficients...
        EMatrix<T, Eigen::Dynamic, 1> a(monomialBasis.size(), 1);
        // ...solve the linear system...
        a = A.colPivHouseholderQr().solve(b);
        // ...and store the solution for later reuse
        localCoefficients.push_back(a);
        //
        ++it;
    }
}

template<unsigned int dim, typename T, typename... list>
void Dcpse<dim, T, list...>::initializeStaticSize(vector_dist<dim, T, aggregate<list...>> &particles,
                                                  unsigned int convergenceOrder, T rCut, T supportSizeFactor)
{
    SupportBuilder<dim, T, aggregate<list...>>
            supportBuilder(particles, differentialSignature, rCut);
    unsigned int requiredSupportSize = monomialBasis.size() * supportSizeFactor;

    auto it = particles.getDomainIterator();
    while (it.isNext())
    {
        // Get the points in the support of the DCPSE kernel and store the support for reuse
        Support<dim, T, aggregate<list...>> support = supportBuilder.getSupport(it, requiredSupportSize);
        EMatrix<T, Eigen::Dynamic, Eigen::Dynamic> V(support.size(), monomialBasis.size());

        // Vandermonde matrix computation
        Vandermonde<dim, T, EMatrix<T, Eigen::Dynamic, Eigen::Dynamic>>
                vandermonde(support, monomialBasis);
        vandermonde.getMatrix(V);

        T eps = vandermonde.getEps();

        localSupports.push_back(support);
        localEps.push_back(eps);
        // Compute the diagonal matrix E
        DcpseDiagonalScalingMatrix<dim> diagonalScalingMatrix(monomialBasis);
        EMatrix<T, Eigen::Dynamic, Eigen::Dynamic> E(support.size(), support.size());
        diagonalScalingMatrix.buildMatrix(E, support, eps);
        // Compute intermediate matrix B
        EMatrix<T, Eigen::Dynamic, Eigen::Dynamic> B = E * V;
        // Compute matrix A
        EMatrix<T, Eigen::Dynamic, Eigen::Dynamic> A = B.transpose() * B;
        // Compute RHS vector b
        DcpseRhs<dim> rhs(monomialBasis, differentialSignature);
        EMatrix<T, Eigen::Dynamic, 1> b(monomialBasis.size(), 1);
        rhs.template getVector<T>(b);
        // Get the vector where to store the coefficients...
        EMatrix<T, Eigen::Dynamic, 1> a(monomialBasis.size(), 1);
        // ...solve the linear system...
        a = A.colPivHouseholderQr().solve(b);
        // ...and store the solution for later reuse
        localCoefficients.push_back(a);
        //
        ++it;
    }
}

#endif // HAVE_EIGEN

#endif //OPENFPM_PDATA_DCPSE_HPP
