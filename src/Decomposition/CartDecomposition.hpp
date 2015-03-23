/*
 * CartDecomposition.hpp
 *
 *  Created on: Aug 15, 2014
 *      Author: Pietro Incardona
 */

#ifndef CARTDECOMPOSITION_HPP
#define CARTDECOMPOSITION_HPP

#include "config.h"
#include "Decomposition.hpp"
#include "Vector/map_vector.hpp"
#include <vector>
#include "global_const.hpp"
#include <initializer_list>
#include "SubdomainGraphNodes.hpp"
#include "metis_util.hpp"
#include "dec_optimizer.hpp"
#include "Space/Shape/Box.hpp"
#include "Space/Shape/Point.hpp"

/**
 * \brief This class decompose a space into subspaces
 *
 * This class decompose a space into regular hyper-cube subspaces, and give the possibilities to
 * select one subspace
 *
 * \tparam dim is the dimensionality of the physical domain we are going to decompose.
 * \tparam T type of the space we decompose, Real, Integer, Complex ...
 * \tparam layout to use
 * \tparam Memory Memory factory used to allocate memory
 * \tparam Domain Structure that contain the information of your physical domain
 * \tparam data type of structure that store the sub-domain decomposition can be an openfpm structure like
 *        vector, ...
 *
 * \note if PARALLEL_DECOMPOSITION macro is defined a parallel decomposition algorithm is used, basically
 *       each processor does not recompute the same decomposition
 *
 */

template<unsigned int dim, typename T, template<typename> class device_l=openfpm::device_cpu, typename Memory=HeapMemory, template<unsigned int, typename> class Domain=Box, template<typename, typename, typename, typename, unsigned int> class data_s = openfpm::vector>
class CartDecomposition
{
public:
	//! Type of the domain we are going to decompose
	typedef T domain_type;

	//! It simplify to access the SpaceBox element
	typedef SpaceBox<dim,T> Box;

private:

	//! This is the key type toaccess  data_s, for example in the case of vector
	//! acc_key is size_t
	typedef typename data_s<SpaceBox<dim,T>,device_l<SpaceBox<dim,T>>,Memory,openfpm::vector_grow_policy_default,openfpm::vect_isel<SpaceBox<dim,T>>::value >::access_key acc_key;

	//! Subspace selected
	//! access_key in case of grid is just the set of the index to access the grid
	std::vector<acc_key> id_sub;

	//! the margin of the sub-domain selected
	SpaceBox<dim,T> sub_domain;

	//! the set of all local sub-domain as vector
	data_s<SpaceBox<dim,T>,device_l<SpaceBox<dim,T>>,Memory,openfpm::vector_grow_policy_default, openfpm::vect_isel<SpaceBox<dim,T>>::value > sub_domains;

	//! Structure that contain for each sub-domain box the processor id
	//! exist for efficient global communication
	openfpm::vector<size_t> fine_s;

	//! number of total sub-domain
	size_t N_tot;

	//! number of sub-domain on each dimension
	size_t div[dim];

	//! rectangular domain to decompose
	Domain<dim,T> domain;

	//! Box Spacing
	T spacing[dim];

	//! Runtime virtual cluster machine
	Vcluster & v_cl;

	/*! \brief Create internally the decomposition
	 *
     * \param v_cl Virtual cluster, used internally to handle or pipeline communication
	 *
	 */
	void CreateDecomposition(Vcluster & v_cl)
	{
		// Calculate the total number of box and and the spacing
		// on each direction

		N_tot = 1;

		// Get the box containing the domain
		SpaceBox<dim,T> bs = domain.getBox();

		for (unsigned int i = 0; i < dim ; i++)
		{
			// Calculate the spacing

			spacing[i] = (bs.getHigh(i) - bs.getLow(i)) / div[i];
			N_tot *= div[i];
		}

		// Here we use METIS

		// Create a cartesian grid graph
		CartesianGraphFactory<dim,Graph_CSR<nm_part_v,nm_part_e>> g_factory_part;

		// Processor graph

		Graph_CSR<nm_part_v,nm_part_e> gp = g_factory_part.template construct<NO_EDGE,T,dim-1>(div,domain);

		// Get the number of processing units
		size_t Np = v_cl.getProcessingUnits();

		// Get the processor id
		long int p_id = v_cl.getProcessUnitID();

		// Convert the graph to metis
		Metis<Graph_CSR<nm_part_v,nm_part_e>> met(gp,Np);

		// decompose

		met.decompose<nm_part_v::id>();

		// fill the structure that store the processor id for each sub-domain

		fine_s.resize(N_tot);

		// Optimize the decomposition creating bigger spaces
		// And reducing Ghost over-stress

		dec_optimizer<dim,Graph_CSR<nm_part_v,nm_part_e>> d_o(gp,div);

		// set of Boxes produced by the decomposition optimizer
		openfpm::vector<::Box<dim,size_t>> loc_box;

		// optimize the decomposition
		d_o.template optimize<nm_part_v::sub_id,nm_part_v::id>(gp,p_id,loc_box);

		// convert into sub-domain
		for (size_t s = 0 ; s < loc_box.size() ; s++)
		{
			SpaceBox<dim,T> sub_d(loc_box.get(s));

			// re-scale with spacing
			sub_d.spacing(spacing);

			// add the sub-domain
			sub_domains.add(sub_d);
		}
	}

	/*! \brief Create the subspaces that decompose your domain
	 *
	 * Create the subspaces that decompose your domain
	 *
	 */

	void CreateSubspaces()
	{
		// Create a grid where each point is a space

		grid<3,void> g(div);

		// create a grid_key_dx iterator

		grid_key_dx_iterator<dim> gk_it(g);

		// Divide the space into subspaces

		while (gk_it.isNext())
		{
			//! iterate through all subspaces
			grid_key_dx<dim> key = gk_it.get();

			//! Create a new subspace

			SpaceBox<dim,T> tmp;

			//! fill with the Margin of the box

			for (int i = 0 ; i < dim ; i++)
			{
				tmp.setHigh(i,(key.get(i)+1)*spacing[i]);
				tmp.setLow(i,key.get(i)*spacing[i]);
			}

			//! add the space box

			sub_domains.add(tmp);

			// add the iterator

			++gk_it;
		}
	}

public:

	/*! \brief Cartesian decomposition copy constructor
	 *
     * \param v_cl Virtual cluster, used internally to handle or pipeline communication
	 *
	 */
	CartDecomposition(CartDecomposition<dim,T,device_l,Memory,Domain,data_s> && cd)
	:sub_domain(cd.sub_domain),N_tot(cd.N_tot),domain(cd.domain),v_cl(cd.v_cl)
	{
		//! Subspace selected
		//! access_key in case of grid is just the set of the index to access the grid
		id_sub.swap(cd.id_sub);

		//! the set of all local sub-domain as vector
		sub_domains.swap(cd.sub_domains);

		for (int i = 0 ; i < dim ; i++)
		{
			this->div[i] = div[dim];

			//! Box Spacing
			this->spacing[i] = spacing[i];
		}
	}

	/*! \brief Cartesian decomposition constructor
	 *
     * \param v_cl Virtual cluster, used internally to handle or pipeline communication
	 *
	 */
	CartDecomposition(Vcluster & v_cl)
	:id_sub(0),N_tot(0),v_cl(v_cl)
	{}

	/*! \brief Cartesian decomposition constructor, it divide the space in boxes
	 *
	 * \param dec is a vector that store how to divide on each dimension
	 * \param domain is the domain to divide
	 * \param v_cl are information of the cluster runtime machine
	 *
	 */
	CartDecomposition(std::vector<size_t> dec, Domain<dim,T> domain, Vcluster & v_cl)
	:id_sub(0),div(dec),domain(domain),v_cl(v_cl)
	{
		// Create the decomposition

		CreateDecomposition(v_cl);
	}

	//! Cartesian decomposition destructor
	~CartDecomposition()
	{}

	/*! \brief processorID return in which processor the particle should go
	 *
	 * \return processorID
	 *
	 */

	template<typename Mem> size_t inline processorID(encapc<1, Point<dim,T>, Mem> p)
	{
		size_t pid = 0;

		for (size_t i = 0 ; i < dim ; i++)
		{
			pid += p.template get<Point<dim,T>::x>()[i];
		}

		return pid;
	}

	/*! \brief processorID return in which processor the particle should go
	 *
	 * \return processorID
	 *
	 */

	size_t inline processorID(T (&p)[dim])
	{
		size_t pid = 0;

		for (size_t i = 0 ; i < dim ; i++)
		{
			pid += p[i];
		}

		return pid;
	}

	/*! \brief Set the parameter of the decomposition
	 *
     * \param div_ std::vector storing into how many domain to decompose on each dimension
     * \param domain_ domain to decompose
	 *
	 */
	void setParameters(std::vector<size_t> div_, Domain<dim,T> domain_)
	{
		// Set the decomposition parameters

		div = div_;
		domain = domain_;

		//! Create the decomposition

		CreateDecomposition(v_cl);
	}

	/*! \brief Set the parameter of the decomposition
	 *
     * \param div_ std::vector storing into how many domain to decompose on each dimension
     * \param domain_ domain to decompose
	 *
	 */
	void setParameters(size_t div_[dim], Domain<dim,T> domain_)
	{
		// Set the decomposition parameters

		for (int i = 0 ; i < dim ; i++)
			div[i] = div_[i];

		domain = domain_;

		//! Create the decomposition

		CreateDecomposition(v_cl);
	}

	/*! \brief Get the number of local local hyper-cubes or sub-domains
	 *
	 * \return the number of sub-domains
	 *
	 */
	size_t getNLocalHyperCube()
	{
		return sub_domains.size();
	}

	/*! The the bulk part of the data set, or the data that
	 * does not depend from the ghosts layers
	 *
	 * \return the bulk of your data
	 *
	 */
	T getBulk()
	{

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

//	dataDiv<T> CartDecomposition<dim,T,layout>::divide(layout::grid<1,Point<dim,T>> & data, neighborhood & nb);

	/*! The the internal part of the data set, or the data that
	 * are inside the local space
	 *
	 * \return the internal part of your data
	 *
	 */
	T getInternal()
	{

	}

	/*! Get the internal part of the dataset, or the data that
	 * depend from the ghost layers
	 *
	 * \return the ghost part of your data
	 *
	 */

	T getBorder()
	{

	}

	/*! Get the external part of the dataset, or the data that
	 * are outside localSpace including ghost
	 *
	 * \return the external part of your data
	 *
	 */
	T getExternal()
	{

	}

	/*! \brief Get the number of one set of hyper-cube enclosing one particular
	 *         subspace, the hyper-cube enclose your space, even if one box is enough
	 *         can be more that one to increase occupancy
	 *
     * In case of Cartesian decomposition it just return 1, each subspace
	 * has one hyper-cube, and occupancy 1
	 *
	 * \param id of the subspace
	 * \return the number of hyper-cube enclosing your space
	 *
	 */
	size_t getNHyperCube(size_t id)
	{
		return 1;
	}

	/*! \brief Get the hyper-cube margins id_c has to be 0
	 *
	 * Get the hyper-cube margins id_c has to be 0, each subspace
	 * has one hyper-cube
	 *
	 * \param id of the subspace
	 * \param id_c
	 * \return The specified hyper-cube space
	 *
	 */
	SpaceBox<dim,T> & getHyperCubeMargins(size_t id, size_t id_c)
	{
#ifdef DEBUG
		// Check if this subspace exist
		if (id >= N_tot)
		{
			std::cerr << "Error CartDecomposition: id > N_tot";
		}
		else if (id_c > 0)
		{
			// Each subspace is an hyper-cube so return error if id_c > 0
			std::cerr << "Error CartDecomposition: id_c > 0";
		}
#endif

		return sub_domains.get<Object>(id);
	}

	/*! \brief Get the total number of Hyper-cube
	 *
	 * Get the total number of Hyper-cube
	 *
	 * \return The total number of hyper-cube
	 *
	 */

	size_t getNHyperCube()
	{
		return N_tot;
	}

	/*! \brief produce an hyper-cube approximation of the space decomposition
	 *
	 */

	void hyperCube()
	{
	}

	/*! \brief Select the local space
	 *
	 * Select the local space
	 *
	 * \param sub select the sub-space
	 *
	 */
	void setSpace(size_t sub)
	{
		id_sub.push_back(sub);
	}


	/*! \brief Get the local grids
	 *
	 * Get the local grids
	 *
	 * \return the local grids
	 *
	 */

	auto getLocalHyperCubes() -> decltype(sub_domains) &
	{
		return sub_domains;
	}

	/*! \brief Get the local hyper-cubes
	 *
	 * Get the local hyper-cubes
	 *
	 * \param lc is the id of the space
	 * \return the local hyper-cube
	 *
	 */

	SpaceBox<dim,T> getLocalHyperCube(size_t lc)
	{
		// Create a space box
		SpaceBox<dim,T> sp;

		// fill the space box

		for (size_t k = 0 ; k < dim ; k++)
		{
			// create the SpaceBox Low and High
			sp.setLow(k,sub_domains.template get<Box::p1>(lc)[k]);
			sp.setHigh(k,sub_domains.template get<Box::p2>(lc)[k]);
		}

		return sp;
	}

	/*! \brief Return the structure that store the physical domain
	 *
	 * Return the structure that store the physical domain
	 *
	 * \return The physical domain
	 *
	 */

	Domain<dim,T> & getDomain()
	{
		return domain;
	}

	/*! \brief Check if the particle is local
	 *
	 * \param p object position
	 *
	 * \return true if it is local
	 *
	 */
	template<typename Mem> bool isLocal(encapc<1, Point<dim,T>, Mem> p)
	{
		return processorID<Mem>() == v_cl.getProcessUnitID();
	}

	/*! \brief Check if the particle is local
	 *
	 * \param p object position
	 *
	 * \return true if it is local
	 *
	 */
	bool isLocal(T (&pos)[dim])
	{
		return processorID(pos) == v_cl.getProcessUnitID();
	}
};


#endif
