/*
 * LB_Model.hpp
 *
 *  Created on: Jan 18, 2017
 *      Author: i-bird
 */

#ifndef SRC_DLB_LB_MODEL_HPP_
#define SRC_DLB_LB_MODEL_HPP_

/*! \brief Linear model
 *
 * The linear model count each particle as weight one
 *
 */
struct ModelLin
{
	size_t factor = 1;

	ModelLin(size_t factor)
	:factor(factor)
	{}

	ModelLin()	{}

	template<typename Decomposition, typename vector> inline void addComputation(Decomposition & dec, const vector & vd, size_t v, size_t p)
	{
		dec.addComputationCost(v, 1);
	}

	template<typename Decomposition> inline void applyModel(Decomposition & dec, size_t v)
	{
		dec.setSubSubDomainComputationCost(v, dec.getSubSubDomainComputationCost(v));
	}

	double distributionTol()
	{
		return 1.01;
	}
};

/*! \brief Linear model
 *
 * The linear model count each particle as weight one
 *
 */
struct ModelSquare
{
	size_t factor = 1;

	template<typename Decomposition, typename vector> inline void addComputation(Decomposition & dec, const vector & vd, size_t v, size_t p)
	{
		dec.addComputationCost(v, 1);
	}

	template<typename Decomposition> inline void applyModel(Decomposition & dec, size_t v)
	{
		dec.setSubSubDomainComputationCost(v, dec.getSubSubDomainComputationCost(v) * dec.getSubSubDomainComputationCost(v));
	}

	double distributionTol()
	{
		return 1.01;
	}
};


#endif /* SRC_DLB_LB_MODEL_HPP_ */
