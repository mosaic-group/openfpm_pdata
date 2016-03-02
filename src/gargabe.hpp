/*
 * gargabe.hpp
 *
 *  Created on: Jan 13, 2015
 *      Author: i-bird
 */

#ifndef GARGABE_HPP_
#define GARGABE_HPP_



	template <unsigned int j, unsigned int i, typename Graph> void optimize(size_t start_p, Graph & graph)
	{
		// We assume that Graph is the rapresentation of a cartesian graph
		// this mean that the direction d is at the child d

		// Create an Hyper-cube

		HyperCube<dim> hyp;

		// Get the number of wavefronts

		size_t n_wf = hyp.getNumberOfElements_R(0);

		// Get the number of intersecting wavefront



		// Get the number of sub-dimensional common wavefront
		// basically are a list of all the subdomain common to two or more

		// Create n_wf wavefront queue

		openfpm::vector<wavefront> v_w;
		v.reserve(n_wf);

		// direction of expansion

		size_t domain_id = 0;
		int exp_dir = 0;
		bool can_expand = true;

		// while is possible to expand

		while (can_expand)
		{
			// for each direction of expansion expand the wavefront

			for (int d = 0 ; d < n_wf ; d++)
			{
				// get the wavefront at direction d

				openfpm::vector<size_t> & wf_d = v_w.get<wavefront::domains>(d);

				// flag to indicate if the wavefront can expand

				bool w_can_expand = true;

				// for each subdomain

				for (size_t sub = 0 ; sub < wf_d.size() ; sub++)
				{
					// check if the adjacent domain in direction d exist
					// and is of the same id

					// get the starting subdomain
					size_t sub_w = wf_d.get<0>(sub);

					// we get the processor id of the neighborhood sub-domain on direction d
					size_t exp_p = graph.getChild(sub_w,d).get<j>();

					// we check if it is the same processor id
					if (exp_p != domain_id)
					{
						w_can_expand = false;
					}
				}

				// if we can expand the wavefront expand it
				if (w_can_expand == true)
				{
					// for each subdomain
					for (size_t sub = 0 ; sub < wf_d.size() ; sub++)
					{
						// update the position of the wavefront
						wf_d.get<0>(sub) = wf_d.get<0>(sub) + gh.stride(d);
					}

					// here we add sub-domains to all the other queues
					// get the face of the hyper-cube

					SubHyperCube<dim,dim-1> sub_hyp = hyp.getSubHyperCube(d);

					std::vector<comb<dim>> q_comb = sub_hyp.getCombinations_R(dim-2);
				}
			}
		}

		// For each point in the Hyper-cube check if we can move the wave front


	}

#ifndef PARALLEL_DECOMPOSITION
//		CreateSubspaces();
#endif

#ifndef USE_METIS_GP

		// Here we do not use METIS
		// Distribute the divided domains

		// Get the number of processing units
		size_t Np = v_cl.getProcessingUnits();

		// Get the ID of this processing unit
		// and push the subspace is taking this
		// processing unit

		for (size_t p_id = v_cl.getProcessUnitID(); p_id < Np ; p_id += Np)
			id_sub.push_back(p_id);
#else


#endif



<<<<<<< HEAD
		/////////////// DEBUG /////////////////////

		// get the decomposition
		auto & dec = g_dist.getDecomposition();

		Vcluster & v_cl = *global_v_cluster;

		// check the consistency of the decomposition
		val = dec.check_consistency();
		BOOST_REQUIRE_EQUAL(val,true);

		// for each local volume
		// Get the number of local grid needed
		size_t n_grid = dec.getNLocalHyperCube();

		size_t vol = 0;

		openfpm::vector<Box<2,size_t>> v_b;

		// Allocate the grids
		for (size_t i = 0 ; i < n_grid ; i++)
		{
			// Get the local hyper-cube
			SpaceBox<2,float> sub = dec.getLocalHyperCube(i);

			Box<2,size_t> g_box = g_dist.getCellDecomposer().convertDomainSpaceIntoGridUnits(sub);
			v_b.add(g_box);

			vol += g_box.getVolumeKey();
		}

		v_cl.reduce(vol);
		v_cl.execute();

		BOOST_REQUIRE_EQUAL(vol,k*k);

		/////////////////////////////////////


		// 3D test

	//	g_dist.write("");

	/*	auto g_it = g_dist.getIteratorBulk();

		auto g_it_halo = g_dist.getHalo();

		// Let try to solve the poisson equation d2(u) = f with f = 1 and computation
		// comunication overlap (100 Jacobi iteration)

		for (int i = 0 ; i < 100 ; i++)
		{
			g_dist.ghost_get();

			// Compute the bulk

			jacobi_iteration(g_it);

			g_dist.ghost_sync();

			// Compute the halo

			jacobi_iteration(g_it_halo);
		}*/


		BOOST_AUTO_TEST_CASE( grid_dist_id_poisson_test_use)
		{
			// grid size
		/*	size_t sz[2] = {1024,1024};

			// Distributed grid with id decomposition

			grid_dist_id<2, scalar<float>, CartDecomposition<2,size_t>> g_dist(sz);

			// Create the grid on memory

			g_dist.Create();*/

		/*	auto g_it = g_dist.getIteratorBulk();

			auto g_it_halo = g_dist.getHalo();

			// Let try to solve the poisson equation d2(u) = f with f = 1 and computation
			// comunication overlap (100 Jacobi iteration)

			for (int i = 0 ; i < 100 ; i++)
			{
				g_dist.ghost_get();

				// Compute the bulk

				jacobi_iteration(g_it);

				g_dist.ghost_sync();

				// Compute the halo

				jacobi_iteration(g_it_halo);
			}*/
		}

		template<typename iterator> void jacobi_iteration(iterator g_it, grid_dist_id<2, float, scalar<float>, CartDecomposition<2,float>> & g_dist)
		{
			// scalar
			typedef scalar<float> S;

			// iterator

			while(g_it.isNext())
			{
				// Jacobi update

				auto pos = g_it.get();

				g_dist.template get<S::ele>(pos) = (g_dist.template get<S::ele>(pos.move(0,1)) +
			                             g_dist.template get<S::ele>(pos.move(0,-1)) +
			                             g_dist.template get<S::ele>(pos.move(1,1)) +
			                             g_dist.template get<S::ele>(pos.move(1,-1)) / 4.0);

				++g_it;
			}
		}

=======

		/*
		 * CartDecomposition.cpp
		 *
		 *  Created on: Aug 15, 2014
		 *      Author: Pietro Incardona
		 */

		#include "CartDecomposition.hpp"



		/*! \brief The the bulk part of the data set, or the data that does not depend
		 *  from the ghosts layers
		 *
		 * The the bulk part of the data set, or the data that does not depend from the
		 *  ghosts layers
		 *
		 */

		/*template<typename T> T CartDecomposition<T>::getBulk(T data)
		{
			// for each element in data

			for (size_t i = 0; i < data.size() ; i++)
			{
				if (localSpace.isInside())
			}

		}

		template<typename T> T CartDecomposition<T>::getInternal()
		{

		}*/

		/*! \brief Check if is border or bulk
		 *
		 * \param neighboorhood define the neighboorhood of all the points
		 * \return true if border, false if bulk
		 *
		 */

		bool borderOrBulk(neighborhood & nb)
		{
			device::grid<1,size_t> nbr = nb.next();

			// check the neighborhood

			// get neighborhood iterator

			grid_key_dx_iterator<dim> iterator_nbr = nbr.getIterator();

			while (iterator_nbr.hasNext())
			{
				grid_key_dx key_nbr = iterator_nbr.next();

				// check if the neighboorhood is internal

				if(subspace.isBound(data.template get<Point::x>(key_nbr)) == false)
				{
					// it is border

					return true;

					ret.bord.push_back(key);
					break;
				}
			}

			return false;
		}

		/*! \brief This function divide the data set into bulk, border, external and internal part
		 *
		 * \tparam dim dimensionality of the structure storing your data
		 *         (example if they are in 3D grid, has to be 3)
		 * \tparam T type of object we are dividing
		 * \tparam device type of layout selected
		 * \param data 1-dimensional grid of point
		 * \param nb define the neighborhood of all the points
		 * \return a structure with the set of objects divided
		 *
		 */

		template<unsigned int dim, typename T, template<typename> class layout, typename Memory, template<unsigned int, typename> class Domain, template<typename, typename, typename> class data_s>
		dataDiv<T> CartDecomposition<dim,T,layout>::divide(device::grid<1,Point<dim,T>> & data, neighborhood & nb)
		{
			//! allocate the 3 subset

			dataDiv<T> ret;

			ret.bord = new boost::shared_ptr<T>(new T());
			ret.inte = new boost::shared_ptr<T>(new T());
			ret.ext = new boost::shared_ptr<T>(new T());

			//! get grid iterator

			grid_key_dx_iterator<dim> iterator = data.getIterator();

			//! we iterate trough all the set of objects

			while (iterator.hasNext())
			{
				grid_key_dx<dim> key = iterator.next();

				//! Check if the object is inside the subspace

				if (subspace.isBound(data.template get<Point<3,T>::x>(key)))
				{
					//! Check if the neighborhood is inside the subspace

					if (borderOrBulk(nb) == true)
					{
						// It is border

						ret.bord.push_back(key);
					}
					else
					{
						// It is bulk

						ret.bulk.push_back(key);
					}
				}
				else
				{
					//! it is external

					ret.ext.push_back(key);
				}
			}
		}


>>>>>>> Jenkin script for taurus


/*! \brief Allocate a set of objects
 *
 * \tparam obj
 * \param n number of object
 *
 * \return an object representing an array of objects
 *
 */
/*	template <typename obj> Vcluster_object_array<obj> allocate(size_t n)
{
	// Vcluster object array
	Vcluster_object_array<obj> vo;

	// resize the array
	vo.resize(n);

	// Create the object on memory and return a Vcluster_object_array
	return vo;
}*/


/*template<typename T>
class Vcluster_object_array : public VObject
{
	std::vector<T> objects;

public:*/

	/*! \brief Constructor of object array
	 *
	 */
/*	Vcluster_object_array()
	{

	}*/

	/*! \brief Return the size of the objects array
	 *
	 * \return the size of the array
	 *
	 */
/*	size_t size() const
	{
		return objects.size();
	}*/

	/*! \brief Return the element i
	 *
	 * \return a reference to the object i
	 *
	 */

/*	T & get(unsigned int i)
	{
		return objects[i];
	}*/

	/*! \brief Return the element i
	 *
	 * \return a reference to the object i
	 *
	 */
/*	const T & get(unsigned int i) const
	{
		return objects[i];
	}*/

	/*! \brief Check if this Object is an array
	 *
	 * \return true, it is an array
	 *
	 */
/*	bool isArray()
	{
		return true;
	}*/

	/*! \brief Destroy the object
	 *
	 */
/*	virtual void destroy()
	{
		// Destroy the objects
		objects.clear();
	}*/

	/*! \brief Get the size of the memory needed to pack the object
	 *
	 * \return the size of the message to pack the object
	 *
	 */
/*	size_t packObjectSize()
	{
		size_t message = 0;

		// Destroy each objects
		for (size_t i = 0 ; i < objects.size() ; i++)
		{
			message += objects[i].packObjectSize();
		}

		return message;
	}*/


	/*! \brief Get the size of the memory needed to pack the object
	 *
	 * \param Memory where to write the packed object
	 *
	 * \return the size of the message to pack the object
	 *
	 */
/*	size_t packObject(void * mem)
	{
		// Pointer is zero
		size_t ptr = 0;
		unsigned char * m = (unsigned char *)mem;

		// pack each object
		for (size_t i = 0 ; i < objects.size() ; i++)
		{
			ptr += objects[i].packObject(&m[ptr]);
		}

#ifdef DEBUG
		if (ptr != packObjectSize())
		{
			std::cerr << "Error " << __FILE__ << " " << __LINE__ << " the pack object size does not match the message" << "\n";
		}
#endif

		return ptr;
	}*/

	/*! \brief Calculate the size to pack an object in the array
	 *
	 * \param array object index
	 *
	 */
/*	size_t packObjectInArraySize(size_t i)
	{
		return objects[i].packObjectSize();
	}*/

	/*! \brief pack the object in the array (the message produced can be used to move one)
	 * object from one processor to another
	 *
	 * \param i index of the object to pack
	 * \param p Memory of the packed object message
	 *
	 */
/*	size_t packObjectInArray(size_t i, void * p)
	{
		return objects[i].packObject(p);
	}*/

	/*! \brief Destroy an object from the array
	 *
	 * \param i object to destroy
	 *
	 */
/*	void destroy(size_t i)
	{
		objects.erase(objects.begin() + i);
	}*/

	/*! \brief Return the object j in the array
	 *
	 * \param j element j
	 *
	 */
/*	T & operator[](size_t j)
	{
		return objects[j];
	}*/

	/*! \brief Return the object j in the array
	 *
	 * \param j element j
	 *
	 */
/*	const T & operator[](size_t j) const
	{
		return objects[j];
	}*/

	/*! \brief Resize the array
	 *
	 * \param size
	 *
	 */
/*	void resize(size_t n)
	{
		objects.resize(n);
	}
};*/

/*! \brief VObject
 *
 * Any object produced by the Virtual cluster (MUST) inherit this class
 *
 */

/*class VObject
{
public:

	// Check if this Object is an array
	virtual bool isArray() = 0;

	// destroy the object
	virtual void destroy() = 0;

	// get the size of the memory needed to pack the object
	virtual size_t packObjectSize() = 0;

	// pack the object
	virtual size_t packObject(void *) = 0;

	// get the size of the memory needed to pack the object in the array
	virtual size_t packObjectInArraySize(size_t i) = 0;

	// pack the object in the array (the message produced can be used to move one)
	// object from one processor to another
	virtual size_t packObjectInArray(size_t i, void * p) = 0;

	// destroy an element from the array
	virtual void destroy(size_t n) = 0;
};*/

#endif /* GARGABE_HPP_ */
