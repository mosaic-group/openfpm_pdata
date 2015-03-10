/*
 * Vector.hpp
 *
 *  Created on: Mar 5, 2015
 *      Author: Pietro Incardona
 */

#ifndef VECTOR_HPP_
#define VECTOR_HPP_

#include "Space/space.hpp"

#define NO_ID false
#define ID true

/*! \brief Distributed vector
 *
 *
 *
 */

template<typename space, typename prop, typename Box, typename Decomposition , typename Memory=HeapMemory, bool with_id=false>
class vector_dist
{
private:

	//! Space Decomposition
	Decomposition dec;

	// Space for space position
	grid_dist_id<1,space,Decomposition,Memory> pos;

	// Space for properties
	grid_dist_id<1,prop,Decomposition,Memory> prp;

	// Virtual cluster
	Vcluster & v_cl;

public:

	/*! \brief Constructor
	 *
	 * \param number of elements
	 *
	 */
	vector_dist(size_t np)
	:dec(Decomposition(*global_v_cluster)),v_cl(*global_v_cluster)
	{
		// resize the position vector
		pos.resize(np);

		// resize the properties vector
		prp.resize(np);

		// Create a valid decomposition of the space
		// Get the number of processor and calculate the number of sub-domain
		// for decomposition
		size_t n_proc = v_cl.getProcessingUnits();
		size_t n_sub = n_proc * SUB_UNIT_FACTOR;

		// Calculate the maximum number (before merging) of sub-domain on
		// each dimension
		size_t div[space::size];
		for (int i = 0 ; i < space::size ; i++)
		{div[i] = round_big_2(pow(n_sub,1.0/space::size));}

		// Create the sub-domains
		dec.setParameters(div);
	}

	/*! \brief Get position of an object
	 *
	 * \param vec_key vector element
	 *
	 */
	template<unsigned int id> auto getPos(size_t vec_key) -> decltype(pos.template get<id>(vec_key))
	{
		return pos.template get<id>(vec_key);
	}

	/*! \brief Get the property of the object
	 *
	 * \param vec_key vector element
	 *
	 */
	template<unsigned int id> auto getProp(size_t vec_key) -> decltype(prp.template get<id>(vec_key))
	{
		return prp.template get<id>(vec_key);
	}

	/*! \brief It communicate the particle to the respective processor
	 *
	 */
	void map()
	{
		// allocate n vector with n = number of processors
//		boost::shared_ptr<openfpm::vector<space>> (new openfpm::vector<space>[v_cl.getProcessingUnits()]);

		// allocate n vector with n = number of processors
//		boost::shared_ptr<openfpm::vector<prop>> (new openfpm::vector<space>[v_cl.getProcessingUnits()]);

		// Contain the map of the processor should communicate
		openfpm::vector<unsigned char> p_map;

		// Contain the processor id of each particle (basically where they have to go)
		openfpm::vector<size_t> lbl_p(pos.size());

		// It contain the list of the processors it should to communicate
		openfpm::vector<size_t> p_list;

		auto it = pos.getIterator();

		// Label all the particles it the processor id where they should go
		while (it.isNext())
		{
			auto key = it.get();

			size_t p_id = dec.processorID(pos.get_o(key));

			lbl_p.get(key) = p_id;

			// It has to communicate
			p_map.get(p_id) = 1;

			++it;
		}
	}


};


#endif /* VECTOR_HPP_ */
