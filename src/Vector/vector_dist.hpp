/*
 * Vector.hpp
 *
 *  Created on: Mar 5, 2015
 *      Author: Pietro Incardona
 */

#ifndef VECTOR_HPP_
#define VECTOR_HPP_

#include "VCluster.hpp"
#include "Space/Shape/Point.hpp"
#include "Vector/vector_dist_iterator.hpp"
#include "Space/Shape/Box.hpp"
#include "Vector/vector_dist_key.hpp"

#define NO_ID false
#define ID true

/*! \brief Distributed vector
 *
 */

template<typename point, typename prop, typename Box, typename Decomposition , typename Memory=HeapMemory, bool with_id=false>
class vector_dist
{
private:

	//! Space Decomposition
	Decomposition dec;

	// Particle position vector for each sub-domain the last one is the unassigned particles vector
	Vcluster_object_array<openfpm::vector<point>> v_pos;

	// Particle properties vector for each sub-domain the last one is the unassigned particles vector
	Vcluster_object_array<openfpm::vector<prop>> v_prp;

	// Virtual cluster
	Vcluster & v_cl;

public:

	/*! \brief Constructor
	 *
	 * \param number of elements
	 *
	 */
	vector_dist(size_t np, Box box)
	:dec(Decomposition(*global_v_cluster)),v_cl(*global_v_cluster)
	{
		// Allocate unassigned particles vectors
		v_pos = v_cl.template allocate<openfpm::vector<point>>(1);
		v_prp = v_cl.template allocate<openfpm::vector<prop>>(1);

		// resize the position vector
		v_pos.get(0).resize(np);

		// resize the properties vector
		v_prp.get(0).resize(np);

		// Create a valid decomposition of the space
		// Get the number of processor and calculate the number of sub-domain
		// for decomposition
		size_t n_proc = v_cl.getProcessingUnits();
		size_t n_sub = n_proc * SUB_UNIT_FACTOR;

		// Calculate the maximum number (before merging) of sub-domain on
		// each dimension
		size_t div[point::dims];
		for (int i = 0 ; i < point::dims ; i++)
		{div[i] = openfpm::math::round_big_2(pow(n_sub,1.0/point::dims));}

		// Create the sub-domains
		dec.setParameters(div,box);
	}

	/*! \brief Get position of an object
	 *
	 * \param vec_key vector element
	 *
	 */
	template<unsigned int id> inline auto getPos(vect_dist_key_dx vec_key) -> decltype(v_pos.get(vec_key.getSub()).template get<id>(vec_key.getKey()))
	{
		return v_pos.get(vec_key.getSub()).template get<id>(vec_key.getKey());
	}

	/*! \brief Get the property of the object
	 *
	 * \param vec_key vector element
	 *
	 */
	template<unsigned int id> inline auto getProp(vect_dist_key_dx vec_key) -> decltype(v_prp.get(vec_key.v_c).template get<id>(vec_key.key))
	{
		return v_prp.get(vec_key.v_c).template get<id>(vec_key.key);
	}

	/*! \brief It communicate the particle to the respective processor
	 *
	 */
	void map()
	{
		// Unassigned particle vector, is always the last one
		size_t up_v = v_pos.size()-1;

		// Contain the map of the processor should communicate
		openfpm::vector<unsigned char> p_map;

		// Contain the processor id of each particle (basically where they have to go)
		openfpm::vector<size_t> lbl_p(v_pos.size());

		// It contain the list of the processors it should to communicate
		openfpm::vector<size_t> p_list;

		auto it = v_pos.get(up_v).getIterator();

		// Label all the particles with the processor id where they should go
		while (it.isNext())
		{
			auto key = it.get();

			size_t p_id = dec.processorID(v_pos.get(up_v).get(key));

			lbl_p.get(key) = p_id;

			// It has to communicate
			p_map.get(p_id) = 1;

			++it;
		}
	}

	/*! \brief Get the iterator across the position of the particles
	 *
	 * \return an iterator
	 *
	 */
	vector_dist_iterator<openfpm::vector<point>> getIterator()
	{
		return vector_dist_iterator<openfpm::vector<point>>(v_pos);
	}

	/*! \brief Get the iterator across the properties of the particles
	 *
	 * \return an iterator
	 *
	 */
	vector_dist_iterator<openfpm::vector<point>> getPropIterator()
	{
		return vector_dist_iterator<openfpm::vector<prop>>(v_prp);
	}

	/*! \brief Get the decomposition
	 *
	 * \return
	 *
	 */
	Decomposition & getDecomposition()
	{
		return dec;
	}

};


#endif /* VECTOR_HPP_ */
