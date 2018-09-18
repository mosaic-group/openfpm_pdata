/*
 * vector_dist_util_unit_tests.hpp
 *
 *  Created on: Feb 14, 2018
 *      Author: i-bird
 */

#ifndef SRC_VECTOR_TESTS_VECTOR_DIST_UTIL_UNIT_TESTS_HPP_
#define SRC_VECTOR_TESTS_VECTOR_DIST_UTIL_UNIT_TESTS_HPP_


/*! \brief Count local and non local
 *
 * \param vd distributed vector
 * \param it iterator
 * \param bc boundary conditions
 * \param box domain box
 * \param dom_ext domain + ghost box
 * \param l_cnt local particles counter
 * \param nl_cnt non local particles counter
 * \param n_out out of domain + ghost particles counter
 *
 */
template<unsigned int dim,typename vector_dist>
inline void count_local_n_local(vector_dist & vd,
								vector_dist_iterator & it,
								size_t (& bc)[dim] ,
								Box<dim,typename vector_dist::stype> & box,
								Box<dim,typename vector_dist::stype> & dom_ext,
								size_t & l_cnt,
								size_t & nl_cnt,
								size_t & n_out)
{
	auto & ct = vd.getDecomposition();

	while (it.isNext())
	{
		auto key = it.get();
		// Check if it is in the domain
		if (box.isInsideNP(vd.getPos(key)) == true)
		{
			Point<dim,typename vector_dist::stype> xp = vd.getPos(key);

			// Check if local
			if (ct.isLocalBC(xp,bc) == true)
				l_cnt++;
			else
				nl_cnt++;
		}
		else
		{
			nl_cnt++;
		}

		Point<dim,typename vector_dist::stype> xp = vd.getPos(key);

		// Check that all particles are inside the Domain + Ghost part
		if (dom_ext.isInside(xp) == false)
				n_out++;

		++it;
	}
}



#endif /* SRC_VECTOR_TESTS_VECTOR_DIST_UTIL_UNIT_TESTS_HPP_ */
