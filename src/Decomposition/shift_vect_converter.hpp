/*
 * shift_vect_converter.hpp
 *
 *  Created on: Feb 8, 2018
 *      Author: i-bird
 */

#ifndef SRC_DECOMPOSITION_SHIFT_VECT_CONVERTER_HPP_
#define SRC_DECOMPOSITION_SHIFT_VECT_CONVERTER_HPP_

#include "Space/Shape/HyperCube.hpp"

/*! \brief in case of high dimensions shift vector converter
 *
 * In case of high-dimensions the number of shift vectors explode, this class
 * handle such case
 *
 */
template<unsigned int dim, typename T, typename Memory, template<typename> class layout_base>
class shift_vect_converter
{
	//! Indicate which indexes are non_periodic
	size_t red_shift_v[dim];

	// indexes
	size_t tmp[dim];

	// Dimension
	int dim_r = 0;

	/*! \brief Here we generare the shift vectors for the low dimension case
	 *
	 * \param domain box that describe the domain
	 *
	 */
	void generateShiftVectors_ld(const Box<dim,T> & domain, size_t (& bc)[dim],
			                     openfpm::vector<Point<dim,T>,Memory,layout_base> & shifts)
	{
		shifts.resize(openfpm::math::pow(3,dim));

		HyperCube<dim> hyp;

		for (long int i = dim ; i >= 0 ; i--)
		{
			std::vector<comb<dim>> cmbs = hyp.getCombinations_R(i);

			for (size_t j = 0 ; j < cmbs.size() ; j++)
			{
				for (size_t k = 0 ; k < dim ; k++)
				{
					switch (cmbs[j][k])
					{
					case 1:
						shifts.get(cmbs[j].lin()).template get<0>()[k] = -(domain.getHigh(k) - domain.getLow(k));
						break;
					case 0:
						shifts.get(cmbs[j].lin()).template get<0>()[k] = 0;
						break;
					case -1:
						shifts.get(cmbs[j].lin()).template get<0>()[k] = (domain.getHigh(k) - domain.getLow(k));
						break;
					}
				}
			}
		}
	}

	/*! \brief Here we generare the shift vectors for the high dimension case
	 *
	 * \param domain box that describe the domain
	 *
	 */
	void generateShiftVectors_hd(const Box<dim,T> & domain, size_t (& bc)[dim],
			                     openfpm::vector<Point<dim,T>,Memory,layout_base> & shifts)
	{
		// get the indexes of the free degree of freedom
		for (size_t i = 0 ; i < dim ; i++)
		{
			if (bc[i] == PERIODIC)
			{
				red_shift_v[dim_r] = i;
				dim_r++;
			}
		}

		HyperCube<dim> hyp;

		// precalculate the nuber of shift vectors
		size_t nsv = 0;
		for (long int i = dim-1 ; i >= 0 ; i--)
		{nsv += hyp.getCombinations_R_bc(i,bc).size();}
		shifts.resize(nsv+1);

		for (long int i = dim-1 ; i >= 0 ; i--)
		{
			std::vector<comb<dim>> cmbs = hyp.getCombinations_R_bc(i,bc);

			for (size_t j = 0 ; j < cmbs.size() ; j++)
			{
				size_t lin_cmb = linId_hd(cmbs[j]);

				for (size_t k = 0 ; k < dim ; k++)
				{
					switch (cmbs[j][k])
					{
					case 1:
						shifts.get(lin_cmb).template get<0>()[k] = -(domain.getHigh(k) - domain.getLow(k));
						break;
					case 0:
						shifts.get(lin_cmb).template get<0>()[k] = 0;
						break;
					case -1:
						shifts.get(lin_cmb).template get<0>()[k] = (domain.getHigh(k) - domain.getLow(k));
						break;
					}
				}
			}
		}
	}

public:

	/*! \brief Here we generare the shift vectors for the low dimension case
	 *
	 * \param domain box that describe the domain
	 *
	 */
	void generateShiftVectors(const Box<dim,T> & domain, size_t (& bc)[dim],
			                  openfpm::vector<Point<dim,T>,Memory,layout_base> & shifts)
	{
		if (dim < 10)
		{generateShiftVectors_ld(domain,bc,shifts);}
		else
		{generateShiftVectors_hd(domain,bc,shifts);}
	}

	/*! \brief Initialize
	 *
	 * \param bc boundary conditions
	 *
	 */
	void Initialize(size_t (& bc)[dim])
	{
		// get the indexes of the free degree of freedom
		for (size_t i = 0 ; i < dim ; i++)
		{
			if (bc[i] == PERIODIC)
			{
				red_shift_v[dim] = i;
				dim_r++;
			}
		}
	}

	/*! \brief linearize the combination in case of high dimension
	 *
	 * \param cmb combination
	 *
	 */
	size_t linId_hd(const comb<dim> & cmb)
	{
		size_t cul = 1;
		size_t lin = 0;
		for (long int i = 0 ; i < dim_r ; i++)
		{
			lin += cul*(cmb.c[red_shift_v[i]] + 1);
			cul *= 3;
		}

		return lin;
	}

	/*! \brief linearize the combination in case of low dimensions
	 *
	 * \param cmb combination
	 *
	 */
	inline size_t linId_ld(const comb<dim> & cmb)
	{
		return cmb.lin();
	}

	/*! \brief linearize the combination in case of high dimensions
	 *
	 * \param cmb combination
	 *
	 */
	inline size_t linId(const comb<dim> & cmb)
	{
		if (dim < 10)
		{return linId_ld(cmb);}

		return linId_hd(cmb);
	}

};


#endif /* SRC_DECOMPOSITION_SHIFT_VECT_CONVERTER_HPP_ */
