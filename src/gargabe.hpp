/*
 * gargabe.hpp
 *
 *  Created on: Jan 13, 2015
 *      Author: i-bird
 */

#ifdef GARGABE_HPP_
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









/*! \brief Impose an operator
 *
 * This function impose an operator on a particular grid region to produce the system
 *
 * Ax = b
 *
 * ## Stokes equation, lid driven cavity with one splipping wall
 *
 * \param op Operator to impose (A term)
 * \param num right hand side of the term (b term)
 * \param id Equation id in the system that we are imposing
 * \param it_d iterator that define where you want to impose
 *
 */
template<typename T> void impose(const T & op , typename Sys_eqs::stype num ,long int id ,grid_dist_iterator_sub<Sys_eqs::dims,typename g_map_type::d_grid> it_d, bool skip_first = false)
{
	//////////////////////// DEBUG /////////////////

	SparseMatrix<double,int> Al;
	Al.load("debug_matrix_single_processor");

	// Construct the map 3 processors 1 processors

	std::unordered_map<size_t,size_t> map_row;

	auto it2 = g_map.getDomainGhostIterator();
	auto ginfo = g_map.getGridInfoVoid();

	while (it2.isNext())
	{
		auto key = it2.get();
		auto key_g = g_map.getGKey(key);
		key_g += pd.getKP1();

		// To linearize must be positive
		bool is_negative = false;
		for (size_t i = 0 ; i < Sys_eqs::dims ; i++)
		{
			if (key_g.get(i) < 0)
				is_negative = true;
		}

		if (is_negative == true)
		{
			++it2;
			continue;
		}

		// Carefull g map is extended, so the original (0,0) is shifted in g_map by

		if (g_map.template get<0>(key) == 7)
		{
			int debug = 0;
			debug++;
		}

		map_row[g_map.template get<0>(key)] = ginfo.LinId(key_g);

		++it2;
	}

	////////////////////////////////////////////////

	Vcluster & v_cl = *global_v_cluster;

	openfpm::vector<triplet> & trpl = A.getMatrixTriplets();

	auto it = it_d;
	grid_sm<Sys_eqs::dims,void> gs = g_map.getGridInfoVoid();

	std::unordered_map<long int,float> cols;

	// resize b if needed
	b.resize(Sys_eqs::nvar * g_map.size());

	bool is_first = skip_first;

	// iterate all the grid points
	while (it.isNext())
	{
		if (is_first == true && v_cl.getProcessUnitID() == 0)
		{
			++it;
			is_first = false;
			continue;
		}
		else
			is_first = false;

		// get the position
		auto key = it.get();

		// Calculate the non-zero colums
		T::value(g_map,key,gs,spacing,cols,1.0);

		//////////// DEBUG //////////////////

		auto g_calc_pos = g_map.getGKey(key);
		g_calc_pos += pd.getKP1();

		/////////////////////////////////////

		// create the triplet

		for ( auto it = cols.begin(); it != cols.end(); ++it )
		{
			trpl.add();
			trpl.last().row() = g_map.template get<0>(key)*Sys_eqs::nvar + id;
			trpl.last().col() = it->first;
			trpl.last().value() = it->second;

			///////////// DEBUG ///////////////////////

			auto ginfo = g_map.getGridInfoVoid();

			size_t r = (trpl.last().row() / Sys_eqs::nvar);
			size_t r_rest = (trpl.last().row() % Sys_eqs::nvar);
			size_t c = (trpl.last().col() / Sys_eqs::nvar);
			size_t c_rest = (trpl.last().col() % Sys_eqs::nvar);
			double val = trpl.last().value();

			// Transform

			size_t rf = map_row[r] * 3 + r_rest;
			size_t cf = map_row[c] * 3 + c_rest;

			auto position_row = ginfo.InvLinId(rf / 3);
			auto position_col = ginfo.InvLinId(cf / 3);

			double valf = Al.getValue(rf,cf);

			if (val != valf)
			{
				int debug = 0;
				debug++;
			}

			///////////////////////////////////////////

//				std::cout << "(" << trpl.last().row() << "," << trpl.last().col() << "," << trpl.last().value() << ")" << "\n";
		}

		b(g_map.template get<0>(key)*Sys_eqs::nvar + id) = num;

		cols.clear();
//			std::cout << "\n";

		// if SE_CLASS1 is defined check the position
#ifdef SE_CLASS1
//			T::position(key,gs,s_pos);
#endif

		++row;
		++row_b;
		++it;
	}
}

typename Sys_eqs::SparseMatrix_type A;

/*! \brief produce the Matrix
 *
 *  \return the Sparse matrix produced
 *
 */
typename Sys_eqs::SparseMatrix_type & getA()
{
#ifdef SE_CLASS1
	consistency();
#endif
	A.resize(g_map.size()*Sys_eqs::nvar,g_map.size()*Sys_eqs::nvar);

	///////////////// DEBUG SAVE //////////////////

//		A.save("debug_matrix_single_processor");

	////////////////////////////////////////////////

	return A;

}


typename Sys_eqs::SparseMatrix_type A;

/*! \brief produce the Matrix
 *
 *  \return the Sparse matrix produced
 *
 */
typename Sys_eqs::SparseMatrix_type & getA()
{
#ifdef SE_CLASS1
	consistency();
#endif
	A.resize(g_map.size()*Sys_eqs::nvar,g_map.size()*Sys_eqs::nvar);

	///////////////// DEBUG SAVE //////////////////

//		A.save("debug_matrix_single_processor");

	////////////////////////////////////////////////

	return A;

}


/*! \brief produce the B vector
 *
 *  \return the vector produced
 *
 */
typename Sys_eqs::Vector_type & getB()
{
#ifdef SE_CLASS1
	consistency();
#endif

	// size of the matrix
//		B.resize(g_map.size()*Sys_eqs::nvar);

	// copy the vector
//		for (size_t i = 0; i < row_b; i++)
//			B.insert(i,b.get(i));

	return b;
}
};




/*! \brief Given an external ghost box, I have an internal ghost box with the same id this function link them
 *
 *
 */
void link_ebox_with_ibox()
{
/*
#ifdef SE_CLASS1

	// No box must be unlinked
	for (size_t i = 0 ; i < proc_int_box.size() ; i++)
	{
		for (size_t j = 0 ; j < proc_int_box.get(i).ibx.size() ; j++)
			proc_int_box.get(i).ibx.get(j).link = -1;

		for (size_t j = 0 ; j < proc_int_box.get(i).ebx.size() ; j++)
			proc_int_box.get(i).ebx.get(j).link= -1;
	}
#endif

	// Get the number of near processors
	for (size_t i = 0 ; i < proc_int_box.size() ; i++)
	{
		std::unordered_map<size_t,std::pair<size_t,size_t>> from_id_to_ibox;
		std::unordered_map<size_t,std::pair<size_t,size_t>> from_id_to_ebox;

		for (size_t j = 0 ; j < getProcessorNIGhost(i) ; j++)
		{
			std::pair<size_t,size_t> & ele = from_id_to_ibox[getProcessorIGhostId(i,j)];
			ele.first = i;
			ele.second = j;
		}

		for (size_t j = 0 ; j < getProcessorNEGhost(i) ; j++)
		{
			std::pair<size_t,size_t> & ele = from_id_to_ebox[getProcessorEGhostId(i,j)];

			ele.first = i;
			ele.second = j;
		}

		// iterate across all the ibox

		for ( auto it = from_id_to_ibox.begin(); it != from_id_to_ibox.end(); ++it )
		{
			auto ite = from_id_to_ebox.find(it->first);

			if(ite == from_id_to_ebox.end())
				std::cerr << __FILE__ << ":" << __LINE__ << " error exist an internal ghost box that does not have an external ghost box" << std::endl;

			if (ite->first != it->first)
				std::cerr << __FILE__ << ":" << __LINE__ << " error exist an internal ghost box with inconsistent information about its origin" << std::endl;

			proc_int_box.get(i).ibx.get(it->second.second).link = ite->second.second;
			proc_int_box.get(i).ebx.get(ite->second.second).link = it->second.second;
		}
	}

#ifdef SE_CLASS1

	// No box must be unlinked
	for (size_t i = 0 ; i < proc_int_box.size() ; i++)
	{
		for (size_t j = 0 ; j < proc_int_box.get(i).ibx.size() ; j++)
		{
			if (proc_int_box.get(i).ibx.get(j).link == -1)
				std::cerr << __FILE__ << ":" << __LINE__ << " error detected unlinked internal ghost box" << std::endl;
		}

		for (size_t j = 0 ; j < proc_int_box.get(i).ebx.size() ; j++)
		{
			if (proc_int_box.get(i).ibx.get(j).link == -1)
				std::cerr << __FILE__ << ":" << __LINE__ << " error detected unlinked external ghost box" << std::endl;
		}
	}
#endif*/


	/*		for (size_t i = 0 ; i < this->getNNProcessors() ; i++)
			{
				for (size_t j = 0 ; j < this->getProcessorNIGhost(i) ; j++)
				{
					size_t id_i = this->getProcessorIGhostId(i,j);
					long int link = this->getProcessorIGhostLink(i,j);

					if (link == -1)
						return false;

					size_t id_e = this->getProcessorEGhostId(i,link);

					if (id_i != id_e)
						return false;
				}
			}*/
}





/////////////////////////////// Fixing  IG BOX not clear if it is really needed /////////////////

/*! \brief Fix the destination box based on the source box
 *
 * in case of periodic grids external ghost box and internal ghost box can miss-match
 * in size if the external ghost box is outside the domain, or more practically
 * if internal and external ghost boxes are linked by periodicity.
 * The two boxes has been calculated in two different way and round-off problem can happen
 * In this call we fix such problem maching the received ghost box to the external ghost box
 *
 * \param bs source box
 * \param dom_i domain from where the source box has been created
 * \param bd destination box
 * \param cmb sector of the destination box
 *
 */
inline bool fix_box_ig(Box<dim,size_t> & bs, Box<dim,long int> & dom_i, const Box<dim,size_t> & bd, comb<dim> & cmb)
{
	// Each dimension must match
	for (size_t k = 0 ; k < dim ; k++)
	{
		size_t iw = bs.getHigh(k) - bs.getLow(k);
		size_t ew = bd.getHigh(k) - bd.getLow(k);

		if (iw != ew)
		{
			std::cout << "Fixing internal external" << std::endl;

			Box<dim,size_t> & bst = bs;

			if (cmb.c[k] == -1)
				bst.setHigh(k,bd.getHigh(k) - (iw - ew));
			else if (cmb.c[k] == 1)
				bst.setLow(k,bs.getLow(k) + (iw - ew));
			else
				return false;

			// points in direction k of the domain
			long int dom_ext = dom_i.getHigh(k) - dom_i.getLow(k);
			// points in direction k of the internal ghost box
			long int ext_ibox = bst.getHigh(k) - bst.getLow(k);

			// internal ghost box cannot be bigger than the domain
			// notify the failure in fixing
			if (dom_ext < ext_ibox)
				return false;

			bs = bst;
		}
	}

	return true;
}

/////////////// GHOST LOCAL FIX


bool ret = fix_box_ig(bx_src,gdb_ext.get(i).Dbox,bx_dst,loc_eg_box.get(sub_id_dst).bid.get(k).cmb);

if (ret == false)
	std::cerr << "ERROR FAIL TO FIX " << std::endl;


/////////////////////


/*! \brief Fix the internal and external ghost box to be consistent
 *
 * in case of periodic grids external ghost box and internal ghost box can miss-match
 * in size if the external ghost box is outside the domain, or more practically
 * if internal and external ghost boxes are linked by periodicity.
 * The two boxes has been calculated in two different way and round-off problem can happen
 * In this call we fix such problem maching each processor communicate its calculate external
 * ghost boxes out of the boundary of the domain the receiving processor fix the size of the
 * connected internal ghost box
 *
 */
inline void fix_ie_g_box()
{
	if (init_fix_ie_g_box == true)	return;

	comb<dim> zero;
	zero.zero();

	// Here we collect all the external ghost box in the sector different from 0 that this processor has

	openfpm::vector<size_t> prc;
	openfpm::vector<size_t> prc_recv;
	openfpm::vector<size_t> sz_recv;
	openfpm::vector<openfpm::vector<Box_fix<dim>>> box_ext_send(dec.getNNProcessors());
	openfpm::vector<openfpm::vector<Box_fix<dim>>> box_ext_recv;

	// It contain the map g_id as key, and the pair, processor id, box-id
	std::unordered_map<long int,std::pair<long int,long int>> iglist;

	// Here we create list of all the internal ghost box linked with an external ghost box
	// by periodicity
	for(size_t i = 0 ; i < dec.getNNProcessors() ; i++)
	{
		for (size_t j = 0 ; j < ig_box.get(i).bid.size() ; j++)
		{
			if (ig_box.get(i).bid.get(j).cmb != zero)
			{
				auto & ele = iglist[ig_box.get(i).bid.get(j).g_id];
				ele.first = i;
				ele.second = j;
			}
		}
	}

	for(size_t i = 0 ; i < dec.getNNProcessors() ; i++)
	{
		for (size_t j = 0 ; j < eg_box.get(i).bid.size() ; j++)
		{
			if (eg_box.get(i).bid.get(j).cmb != zero)
			{
				box_ext_send.get(i).add();
				box_ext_send.get(i).last().bx = eg_box.get(i).bid.get(j).l_e_box;
				box_ext_send.get(i).last().g_id = eg_box.get(i).bid.get(j).g_id;
			}
		}
		prc.add(dec.IDtoProc(i));
	}

	v_cl.SSendRecv(box_ext_send,box_ext_recv,prc,prc_recv,sz_recv);

	// Received the external boxes we do fixation for each processor
	for (size_t i = 0 ; i < box_ext_recv.size() ; i++)
	{
		// For each received external ghost box
		for (size_t j = 0 ; j < box_ext_recv.get(i).size() ; j++)
		{
			// ig box linked
			size_t proc_id = dec.ProctoID(prc_recv.get(i));

			auto it = g_id_to_internal_ghost_box.get(proc_id).find(box_ext_recv.get(i).get(j).g_id);

#ifdef SE_CLASS1

			if (it == g_id_to_internal_ghost_box.get(proc_id).end())
			{
				std::cerr << __FILE__ << ":" << __LINE__ << " warning unlinked external ghost box" << std::endl;
				continue;
			}

#endif

			size_t link = it->second;

			Box<dim,size_t> & box_i = ig_box.get(proc_id).bid.get(link).box;

			// local Sub-domain from where this internal ghost box is calculated
			Box<dim,long int> & box_sub_i = gdb_ext.get(ig_box.get(proc_id).bid.get(link).sub).Dbox;

			comb<dim> cmb = ig_box.get(proc_id).bid.get(link).cmb;

			// the fixing can fail
			// if it fail put the ig_box into a list
			// The fix can fail (for example) if the external ghost box require 7 point on x
			// but the domain has 6 point, in this case we cannot correct the internal ghost box
			bool ret = fix_box_ig(box_i,box_sub_i,box_ext_recv.get(i).get(j).bx,cmb);

			if (ret == false)
				std::cerr << __FILE__ << ":" << __LINE__ << " and inconsistency between internal and external ghost boxes has been detected. The fix is not possible please change your ghost size (by a small amount) on the order of 10^-5 if you use float 10^-14 if you use double"  << std::endl;

			// Invalidate the ig_box in the list
			auto & ele = iglist[box_ext_recv.get(i).get(j).g_id];
			ele.first = -1;
			ele.second = -1;
		}
	}

	// Here we check if all the internal ghost box has been explored
	// if one internal ghost box has not been explored, it been that, there is not
	// corresponding external ghost box on the other side. so we invalidate

	for ( auto it = iglist.begin(); it != iglist.end(); ++it )
	{
		// If has not been explored invalidate, there is not external ghost
		if (it->second.first != -1)
		{
			size_t a = it->second.first;
			size_t b = it->second.second;
			ig_box.get(a).bid.get(b).box.invalidate();
		}
	}
}


//////////////////////////////////////////////////////////////

// Fix the exteenal and internal ghost boxes in ghost get
fix_ie_g_box();

//////////////////////


#endif /* GARGABE_HPP_ */
