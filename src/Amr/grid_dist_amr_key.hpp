/*
 * grid_dist_amr_key.hpp
 *
 *  Created on: Sep 23, 2017
 *      Author: i-bird
 */

#ifndef SRC_AMR_GRID_DIST_AMR_KEY_HPP_
#define SRC_AMR_GRID_DIST_AMR_KEY_HPP_

/*! \brief Amr grid distributed key
 *
 * \tparam dim dimensionality
 *
 */
template<unsigned int dim>
class grid_dist_amr_key
{
	//! actual level
	size_t lvl;

	//! actual position in the distributed grid
	grid_dist_key_dx<dim> key;


public:

	/*! \constructor
	 *
	 * \param lvl level
	 * \param key distributed grid key
	 * \param offsets to move between levels
	 *
	 */
	inline grid_dist_amr_key(size_t lvl,
			          grid_dist_key_dx<dim> key)
	:lvl(lvl),key(key)
	{}

	/*! \brief Return the grid key
	 *
	 * \return the distributed key
	 *
	 */
	inline const grid_dist_key_dx<dim> & getKey() const
	{
		return key;
	}

	/*! \brief Return the grid key (as reference)
	 *
	 * \return the distributed key
	 *
	 */
	inline grid_dist_key_dx<dim> & getKeyRef()
	{
		return key;
	}


	/*! \brief Return the level
	 *
	 * \return the level
	 *
	 */
	inline size_t getLvl() const
	{
		return lvl;
	}

	/*! \brief Return the level
	 *
	 * \param lvl level to set
	 *
	 */
	inline void setLvl(size_t lvl)
	{
		this->lvl = lvl;
	}

	/*! \brief Create a new key moving the old one
	 *
	 * \param s dimension id
	 * \param s number of steps
	 *
	 * \return new key
	 *
	 */
	inline grid_dist_amr_key<dim> moveSpace(size_t d,size_t s)
	{
		return grid_dist_amr_key<dim>(lvl,key.move(d,s));
	}
};




#endif /* SRC_AMR_GRID_DIST_AMR_KEY_HPP_ */
