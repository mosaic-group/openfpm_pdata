/*
 * Vector.hpp
 *
 *  Created on: Mar 5, 2015
 *      Author: Pietro Incardona
 */

#ifndef VECTOR_HPP_
#define VECTOR_HPP_

#include "config.h"
#include "HDF5_wr/HDF5_wr.hpp"
#include "VCluster/VCluster.hpp"
#include "Space/Shape/Point.hpp"
#include "Vector/Iterators/vector_dist_iterator.hpp"
#include "Space/Shape/Box.hpp"
#include "Vector/vector_dist_key.hpp"
#include "memory/PtrMemory.hpp"
#include "NN/CellList/CellList.hpp"
#include "NN/CellList/CellListFast_gen.hpp"
#include "util/common.hpp"
#include "util/object_util.hpp"
#include "memory/ExtPreAlloc.hpp"
#include "CSVWriter/CSVWriter.hpp"
#include "VTKWriter/VTKWriter.hpp"
#include "Decomposition/common.hpp"
#include "Grid/Iterators/grid_dist_id_iterator_dec.hpp"
#include "Grid/grid_key_dx_iterator_hilbert.hpp"
#include "Vector/vector_dist_ofb.hpp"
#include "Decomposition/CartDecomposition.hpp"
#include "data_type/aggregate.hpp"
#include "NN/VerletList/VerletList.hpp"
#include "vector_dist_comm.hpp"
#include "DLB/LB_Model.hpp"
#include "Vector/vector_map_iterator.hpp"
#include "NN/CellList/ParticleIt_Cells.hpp"
#include "NN/CellList/ProcKeys.hpp"
#include "Vector/vector_dist_kernel.hpp"
#include "NN/CellList/cuda/CellList_gpu.hpp"

#define DEC_GRAN(gr) ((size_t)gr << 32)

#define VECTOR_DIST_ERROR_OBJECT std::runtime_error("Runtime vector distributed error");

#ifdef SE_CLASS3
#include "se_class3_vector.hpp"
#endif

#ifdef SE_CLASS3
	#define SE_CLASS3_VDIST_CONSTRUCTOR ,se3(getDecomposition(),*this)
#else
	#define SE_CLASS3_VDIST_CONSTRUCTOR
#endif


#define NO_ID false
#define ID true

// Perform a ghost get or a ghost put
#define GET	1
#define PUT 2

// Write the particles with ghost
#define NO_GHOST 0
#define WITH_GHOST 2

#define GCL_NON_SYMMETRIC 0
#define GCL_SYMMETRIC 1
#define GCL_HILBERT 2


//! General function t get a cell-list
template<unsigned int dim, typename St, typename CellL, typename Vector, unsigned int impl>
struct gcl
{
	/*! \brief Get the Cell list based on the type
	 *
	 * \param vd Distributed vector
	 * \param r_cut Cut-off radius
	 * \param g Ghost
	 *
	 * \return the constructed cell-list
	 *
	 */
	static inline CellL get(Vector & vd, const St & r_cut, const Ghost<dim,St> & g)
	{
		return vd.template getCellList<CellL>(r_cut);
	}
};

//! General function t get a cell-list
template<unsigned int dim, typename St, typename CellL, typename Vector>
struct gcl<dim,St,CellL,Vector,GCL_HILBERT>
{
	/*! \brief Get the Cell list based on the type
	 *
	 * \param vd Distributed vector
	 * \param r_cut Cut-off radius
	 * \param g Ghost
	 *
	 * \return the constructed cell-list
	 *
	 */
	static inline CellL get(Vector & vd, const St & r_cut, const Ghost<dim,St> & g)
	{
		return vd.getCellList_hilb(r_cut,g);
	}
};

//! General function t get a cell-list
template<unsigned int dim, typename St, typename CellL, typename Vector>
struct gcl<dim,St,CellL,Vector,GCL_SYMMETRIC>
{
	/*! \brief Get the Cell list based on the type
	 *
	 * \param vd Distributed vector
	 * \param r_cut Cut-off radius
	 * \param g Ghost
	 *
	 * \return the constructed cell-list
	 *
	 */
	static inline CellL get(Vector & vd, const St & r_cut, const Ghost<dim,St> & g)
	{
		return vd.getCellListSym(r_cut);
	}
};

/////////////////// GCL anisotropic ///////////////////////////////////////////

//! General function t get a cell-list
template<unsigned int dim, typename St, typename CellL, typename Vector, unsigned int impl>
struct gcl_An
{
	/*! \brief Get the Cell list based on the type
	 *
	 * \param vd Distributed vector
	 * \param r_cut Cut-off radius
	 * \param g Ghost
	 *
	 * \return the constructed cell-list
	 *
	 */
	static inline CellL get(Vector & vd, const size_t (& div)[dim], const size_t (& pad)[dim], const Ghost<dim,St> & g)
	{
		return vd.template getCellListSym<CellL>(div,pad);
	}
};

#define CELL_MEMFAST(dim,St) CellList_gen<dim, St, Process_keys_lin, Mem_fast<>, shift<dim, St> >
#define CELL_MEMBAL(dim,St) CellList_gen<dim, St, Process_keys_lin, Mem_bal<>, shift<dim, St> >
#define CELL_MEMMW(dim,St) CellList_gen<dim, St, Process_keys_lin, Mem_mw<>, shift<dim, St> >

#define CELL_MEMFAST_HILB(dim,St) CellList_gen<dim, St, Process_keys_hilb, Mem_fast<>, shift<dim, St> >
#define CELL_MEMBAL_HILB(dim,St) CellList_gen<dim, St, Process_keys_hilb, Mem_bal<>, shift<dim, St> >
#define CELL_MEMMW_HILB(dim,St) CellList_gen<dim, St, Process_keys_hilb, Mem_mw<>, shift<dim, St> >

#define VERLET_MEMFAST(dim,St) VerletList<dim,St,Mem_fast<>,shift<dim,St> >
#define VERLET_MEMBAL(dim,St)  VerletList<dim,St,Mem_bal<>,shift<dim,St> >
#define VERLET_MEMMW(dim,St)   VerletList<dim,St,Mem_mw<>,shift<dim,St> >

#define VERLET_MEMFAST_INT(dim,St) VerletList<dim,St,Mem_fast<HeapMemory,unsigned int>,shift<dim,St> >
#define VERLET_MEMBAL_INT(dim,St)  VerletList<dim,St,Mem_bal<unsigned int>,shift<dim,St> >
#define VERLET_MEMMW_INT(dim,St)   VerletList<dim,St,Mem_mw<unsigned int>,shift<dim,St> >

enum reorder_opt
{
	NO_REORDER = 0,
	HILBERT = 1,
	LINEAR = 2
};

/*! \brief Distributed vector
 *
 * This class reppresent a distributed vector, the distribution of the structure
 * is based on the positional information of the elements the vector store
 *
 * ## Create a vector of random elements on each processor 2D
 * \snippet vector_dist_unit_test.hpp Create a vector of random elements on each processor 2D
 *
 * ## Create a vector of random elements on each processor 3D
 * \snippet vector_dist_unit_test.hpp Create a vector of random elements on each processor 3D
 *
 * ## Create a vector of elements distributed on a grid like way
 * \snippet vector_dist_unit_test.hpp Create a vector of elements distributed on a grid like way
 *
 * ## Redistribute the particles and sync the ghost properties
 * \snippet vector_dist_unit_test.hpp Redistribute the particles and sync the ghost properties
 *
 * \tparam dim Dimensionality of the space where the elements lives
 * \tparam St type of space float, double ...
 * \tparam prop properties the vector element store in OpenFPM data structure format
 * \tparam Decomposition Decomposition strategy to use CartDecomposition ...
 * \tparam Memory Memory pool where store the information HeapMemory ...
 * \tparam Memory layout
 *
 */

template<unsigned int dim,
         typename St,
         typename prop,
         typename Decomposition = CartDecomposition<dim,St>,
         typename Memory = HeapMemory,
         template<typename> class layout_base = memory_traits_lin>
class vector_dist : public vector_dist_comm<dim,St,prop,Decomposition,Memory,layout_base>
{
public:

	//! Self type
	typedef vector_dist<dim,St,prop,Decomposition,Memory,layout_base> self;

	//! property object
	typedef prop value_type;

private:

	//! Ghost marker, all the particle with id > g_m are ghost all with g_m < are real particle
	size_t g_m = 0;

	//! Particle position vector, (It has 2 elements) the first has real particles assigned to a processor
	//! the second element contain unassigned particles
	openfpm::vector<Point<dim, St>,Memory,typename layout_base<Point<dim,St>>::type,layout_base> v_pos;

	//! Particle properties vector, (It has 2 elements) the first has real particles assigned to a processor
	//! the second element contain unassigned particles
	openfpm::vector<prop,Memory,typename layout_base<prop>::type,layout_base> v_prp;

#ifdef CUDA_GPU

	openfpm::vector<prop,Memory,typename layout_base<prop>::type,layout_base> v_prp_out;

	openfpm::vector<Point<dim, St>,Memory,typename layout_base<Point<dim,St>>::type,layout_base> v_pos_out;

#endif

	//! Virtual cluster
	Vcluster<Memory> & v_cl;

	//! option used to create this vector
	size_t opt = 0;

	//! Name of the properties
	openfpm::vector<std::string> prp_names;

#ifdef SE_CLASS3

	se_class3_vector<prop::max_prop,dim,St,Decomposition,self> se3;

#endif

	/*! \brief Initialize the structures
	 *
	 * \param np number of particles
	 *
	 */
	void init_structures(size_t np)
	{
		// convert to a local number of elements
		size_t p_np = np / v_cl.getProcessingUnits();

		// Get non divisible part
		size_t r = np % v_cl.getProcessingUnits();

		// Distribute the remain particles
		if (v_cl.getProcessUnitID() < r)
			p_np++;

		// resize the position vector
		v_pos.resize(p_np);

		// resize the properties vector
		v_prp.resize(p_np);

		g_m = p_np;
	}

	/*! \brief Check if the parameters describe a valid vector. In case it does not report an error
	 *
	 * \param box Box to check
	 *
	 */
	void check_parameters(Box<dim,St> & box)
	{
		// if the box is not valid return an error
		if (box.isValid() == false)
		{
			std::cerr << __FILE__ << ":" << __LINE__ << "  Error the domain is not valid " << box.toString() << std::endl;
			ACTION_ON_ERROR(VECTOR_DIST_ERROR_OBJECT);
		}

	}

	/*! \brief It check that the r_cut is not bugger than the ghost
	 *
	 * \param r_cut cut-off radius
	 *
	 */
	void check_ghost_compatible_rcut(St r_cut)
	{
		for (size_t i = 0 ; i < dim ; i++)
		{
			if (fabs(getDecomposition().getGhost().getLow(i)) < r_cut)
			{
				std::cerr << __FILE__ << ":" << __LINE__ << " Error the cut off radius " << r_cut << " is bigger that the ghost layer on the dimension " << i << " lower=" << getDecomposition().getGhost().getLow(i) << std::endl;
				ACTION_ON_ERROR(VECTOR_DIST_ERROR_OBJECT);
			}
		}
	}

	/*! \brief Reorder based on hilbert space filling curve
	 *
	 * \param v_pos_dest reordered vector of position
	 * \param v_prp_dest reordered vector of properties
	 * \param m order of the space filling curve
	 * \param cell_list cell-list
	 *
	 */
	template<typename CellL, typename sfc_it>
	void reorder_sfc(openfpm::vector<Point<dim,St>> & v_pos_dest,
						 openfpm::vector<prop> & v_prp_dest,
						 sfc_it & h_it,
						 CellL & cell_list)
	{
		v_pos_dest.resize(v_pos.size());
		v_prp_dest.resize(v_prp.size());

		//Index for v_pos_dest
		size_t count = 0;

		grid_key_dx<dim> ksum;

		for (size_t i = 0; i < dim ; i++)
		{ksum.set_d(i,cell_list.getPadding(i));}

		while (h_it.isNext())
		{
		  auto key = h_it.get();
		  key += ksum;

		  size_t lin = cell_list.getGrid().LinId(key);

		  // for each particle in the Cell "lin"
		  for (size_t i = 0; i < cell_list.getNelements(lin); i++)
		  {
			  //reorder
			  auto v = cell_list.get(lin,i);
			  v_pos_dest.get(count) = v_pos.get(v);
			  v_prp_dest.get(count) = v_prp.get(v);

			  count++;
		  }
		  ++h_it;
		}
	}

public:

	//! space type
	typedef St stype;

	//! dimensions of space
	static const unsigned int dims = dim;

	/*! \brief Operator= for distributed vector
	 *
	 * \param v vector to copy
	 *
	 * \return itself
	 *
	 */
	vector_dist<dim,St,prop,Decomposition,Memory,layout_base> &
	operator=(const vector_dist<dim,St,prop,Decomposition,Memory,layout_base> & v)
	{
		static_cast<vector_dist_comm<dim,St,prop,Decomposition,Memory,layout_base> *>(this)->operator=(static_cast<vector_dist_comm<dim,St,prop,Decomposition,Memory,layout_base>>(v));

		g_m = v.g_m;
		v_pos = v.v_pos;
		v_prp = v.v_prp;

#ifdef SE_CLASS3
		se3 = v.se3;
#endif

		opt = v.opt;

		return *this;
	}

	/*! \brief Operator= for distributed vector
	 *
	 * \param v vector to copy
	 *
	 * \return itself
	 *
	 */
	vector_dist<dim,St,prop,Decomposition,Memory,layout_base> &
	operator=(vector_dist<dim,St,prop,Decomposition,Memory,layout_base> && v)
	{
		static_cast<vector_dist_comm<dim,St,prop,Decomposition,Memory,layout_base> *>(this)->operator=(static_cast<vector_dist_comm<dim,St,prop,Decomposition,Memory,layout_base> >(v));

		g_m = v.g_m;
		v_pos.swap(v.v_pos);
		v_prp.swap(v.v_prp);

#ifdef SE_CLASS3
		se3 = v.se3;
#endif

		opt = v.opt;

		return *this;
	}


	/*! \brief Copy Constructor
	 *
	 * \param v vector to copy
	 *
	 */
	vector_dist(const vector_dist<dim,St,prop,Decomposition,Memory,layout_base> & v)
	:vector_dist_comm<dim,St,prop,Decomposition,Memory,layout_base>(v.getDecomposition()),v_cl(v.v_cl) SE_CLASS3_VDIST_CONSTRUCTOR
	{
#ifdef SE_CLASS2
		check_new(this,8,VECTOR_DIST_EVENT,4);
#endif

		this->operator=(v);
	}

	/*! \brief Copy constructor
	 *
	 * \param v vector to copy
	 *
	 */
	vector_dist(vector_dist<dim,St,prop,Decomposition,Memory,layout_base> && v) noexcept
	:v_cl(v.v_cl) SE_CLASS3_VDIST_CONSTRUCTOR
	{
#ifdef SE_CLASS2
		check_new(this,8,VECTOR_DIST_EVENT,4);
#endif

		this->operator=(v);

#ifdef SE_CLASS3
		se3.Initialize();
#endif
	}

	/*! \brief Constructor with predefined decomposition
	 *
	 * \param dec is the decomposition
	 * \param np number of particles
	 *
	 */
	vector_dist(const Decomposition & dec, size_t np) :
	vector_dist_comm<dim,St,prop,Decomposition,Memory,layout_base>(dec), v_cl(create_vcluster<Memory>()) SE_CLASS3_VDIST_CONSTRUCTOR
	{
#ifdef SE_CLASS2
		check_new(this,8,VECTOR_DIST_EVENT,4);
#endif

		init_structures(np);

#ifdef SE_CLASS3
		se3.Initialize();
#endif
	}


	/*! \brief Constructor of a distributed vector
	 *
	 * \param np number of elements
	 * \param box domain where the vector of elements live
	 * \param bc boundary conditions
	 * \param g Ghost margins
	 * \param opt [Optional] additional options. BIND_DEC_TO_GHOST Bind the decomposition to be multiple of the
	 *          ghost size. This is required if we want to use symmetric to eliminate
	 *          ghost communications.
	 * \param gdist [Optional] override the default distribution grid
	 *
	 */
	vector_dist(size_t np, Box<dim, St> box, const size_t (&bc)[dim], const Ghost<dim, St> & g, size_t opt = 0, const grid_sm<dim,void> & gdist = grid_sm<dim,void>())
	:v_cl(create_vcluster<Memory>()),opt(opt) SE_CLASS3_VDIST_CONSTRUCTOR
	{
#ifdef SE_CLASS2
		check_new(this,8,VECTOR_DIST_EVENT,4);
#endif

		if (opt >> 32 != 0)
			this->setDecompositionGranularity(opt >> 32);

		check_parameters(box);

		init_structures(np);
		this->init_decomposition(box,bc,g,opt,gdist);

#ifdef SE_CLASS3
		se3.Initialize();
#endif
	}

	~vector_dist()
	{
#ifdef SE_CLASS2
		check_delete(this);
#endif
	}

	/*! \brief remove all the elements
	 *
	 *
	 */
	void clear()
	{
		resize(0);
	}

	/*! \brief return the local size of the vector
	 *
	 * \return local size
	 *
	 */
	size_t size_local() const
	{
		return g_m;
	}

	/*! \brief return the local size of the vector
	 *
	 * \return local size
	 *
	 */
	size_t size_local_with_ghost() const
	{
		return v_pos.size();
	}

#ifndef ONLY_READWRITE_GETTER

	/*! \brief Get the position of an element
	 *
	 * see the vector_dist iterator usage to get an element key
	 *
	 * \param vec_key element
	 *
	 * \return the position of the element in space
	 *
	 */
	inline auto getPos(vect_dist_key_dx vec_key) -> decltype(v_pos.template get<0>(vec_key.getKey()))
	{
#ifdef SE_CLASS3
		check_for_pos_nan_inf<prop::max_prop_real,prop::max_prop>(*this,vec_key.getKey());
#endif

		return v_pos.template get<0>(vec_key.getKey());
	}

	/*! \brief Get the position of an element
	 *
	 * see the vector_dist iterator usage to get an element key
	 *
	 * \param vec_key element
	 *
	 * \return the position of the element in space
	 *
	 */
	inline auto getPos(vect_dist_key_dx vec_key) const -> decltype(v_pos.template get<0>(vec_key.getKey()))
	{
#ifdef SE_CLASS3
		check_for_pos_nan_inf<prop::max_prop_real,prop::max_prop>(*this,vec_key.getKey());
#endif
		return v_pos.template get<0>(vec_key.getKey());
	}

	/*! \brief Get the position of an element
	 *
	 * see the vector_dist iterator usage to get an element key
	 *
	 * \param vec_key element
	 *
	 * \return the position of the element in space
	 *
	 */
	inline auto getPos(size_t vec_key) -> decltype(v_pos.template get<0>(vec_key))
	{
#ifdef SE_CLASS3
		check_for_pos_nan_inf<prop::max_prop_real,prop::max_prop>(*this,vec_key);
#endif
		return v_pos.template get<0>(vec_key);
	}

	/*! \brief Get the position of an element
	 *
	 * see the vector_dist iterator usage to get an element key
	 *
	 * \param vec_key element
	 *
	 * \return the position of the element in space
	 *
	 */
	inline auto getPos(size_t vec_key) const -> decltype(v_pos.template get<0>(vec_key))
	{
#ifdef SE_CLASS3
		check_for_pos_nan_inf<prop::max_prop_real,prop::max_prop>(*this,vec_key);
#endif
		return v_pos.template get<0>(vec_key);
	}

	/*! \brief Get the property of an element
	 *
	 * see the vector_dist iterator usage to get an element key
	 *
	 * \tparam id property id
	 * \param vec_key vector element
	 *
	 * \return return the selected property of the vector element
	 *
	 */
	template<unsigned int id> inline auto getProp(vect_dist_key_dx vec_key) -> decltype(v_prp.template get<id>(vec_key.getKey()))
	{
#ifdef SE_CLASS3
		check_for_prop_nan_inf<id,prop::max_prop+SE3_STATUS>(*this,vec_key.getKey());
#endif
		return v_prp.template get<id>(vec_key.getKey());
	}

	/*! \brief Get the property of an element
	 *
	 * see the vector_dist iterator usage to get an element key
	 *
	 * \tparam id property id
	 * \param vec_key vector element
	 *
	 * \return return the selected property of the vector element
	 *
	 */
	template<unsigned int id> inline auto getProp(vect_dist_key_dx vec_key) const -> decltype(v_prp.template get<id>(vec_key.getKey()))
	{
#ifdef SE_CLASS3
		check_for_prop_nan_inf<id,prop::max_prop+SE3_STATUS>(*this,vec_key.getKey());
#endif
		return v_prp.template get<id>(vec_key.getKey());
	}

	/*! \brief Get the property of an element
	 *
	 * see the vector_dist iterator usage to get an element key
	 *
	 * \tparam id property id
	 * \param vec_key vector element
	 *
	 * \return return the selected property of the vector element
	 *
	 */
	template<unsigned int id> inline auto getProp(size_t vec_key) -> decltype(v_prp.template get<id>(vec_key))
	{
#ifdef SE_CLASS3
		check_for_prop_nan_inf<id,prop::max_prop+SE3_STATUS>(*this,vec_key);
#endif
		return v_prp.template get<id>(vec_key);
	}

	/*! \brief Get the property of an element
	 *
	 * see the vector_dist iterator usage to get an element key
	 *
	 * \tparam id property id
	 * \param vec_key vector element
	 *
	 * \return return the selected property of the vector element
	 *
	 */
	template<unsigned int id> inline auto getProp(size_t vec_key) const -> decltype(v_prp.template get<id>(vec_key))
	{
#ifdef SE_CLASS3
		check_for_prop_nan_inf<id,prop::max_prop+SE3_STATUS>(*this,vec_key);
#endif
		return v_prp.template get<id>(vec_key);
	}

#endif

///////////////////// Read and write with no check

	/*! \brief Get the position of an element
	 *
	 * see the vector_dist iterator usage to get an element key
	 *
	 * \param vec_key element
	 *
	 * \return the position of the element in space
	 *
	 */
	inline auto getPosNC(vect_dist_key_dx vec_key) -> decltype(v_pos.template get<0>(vec_key.getKey()))
	{
		return v_pos.template get<0>(vec_key.getKey());
	}

	/*! \brief Get the position of an element
	 *
	 * see the vector_dist iterator usage to get an element key
	 *
	 * \param vec_key element
	 *
	 * \return the position of the element in space
	 *
	 */
	inline auto getPosNC(vect_dist_key_dx vec_key) const -> decltype(v_pos.template get<0>(vec_key.getKey()))
	{
		return v_pos.template get<0>(vec_key.getKey());
	}

	/*! \brief Get the position of an element
	 *
	 * see the vector_dist iterator usage to get an element key
	 *
	 * \param vec_key element
	 *
	 * \return the position of the element in space
	 *
	 */
	inline auto getPosNC(size_t vec_key) -> decltype(v_pos.template get<0>(vec_key))
	{
		return v_pos.template get<0>(vec_key);
	}

	/*! \brief Get the position of an element
	 *
	 * see the vector_dist iterator usage to get an element key
	 *
	 * \param vec_key element
	 *
	 * \return the position of the element in space
	 *
	 */
	inline auto getPosNC(size_t vec_key) const -> decltype(v_pos.template get<0>(vec_key))
	{
		return v_pos.template get<0>(vec_key);
	}

	/*! \brief Get the property of an element
	 *
	 * see the vector_dist iterator usage to get an element key
	 *
	 * \tparam id property id
	 * \param vec_key vector element
	 *
	 * \return return the selected property of the vector element
	 *
	 */
	template<unsigned int id> inline auto getPropNC(vect_dist_key_dx vec_key) -> decltype(v_prp.template get<id>(vec_key.getKey()))
	{
		return v_prp.template get<id>(vec_key.getKey());
	}

	/*! \brief Get the property of an element
	 *
	 * see the vector_dist iterator usage to get an element key
	 *
	 * \tparam id property id
	 * \param vec_key vector element
	 *
	 * \return return the selected property of the vector element
	 *
	 */
	template<unsigned int id> inline auto getPropNC(vect_dist_key_dx vec_key) const -> decltype(v_prp.template get<id>(vec_key.getKey()))
	{
		return v_prp.template get<id>(vec_key.getKey());
	}

	/*! \brief Get the property of an element
	 *
	 * see the vector_dist iterator usage to get an element key
	 *
	 * \tparam id property id
	 * \param vec_key vector element
	 *
	 * \return return the selected property of the vector element
	 *
	 */
	template<unsigned int id> inline auto getPropNC(size_t vec_key) -> decltype(v_prp.template get<id>(vec_key))
	{
		return v_prp.template get<id>(vec_key);
	}

	/*! \brief Get the property of an element
	 *
	 * see the vector_dist iterator usage to get an element key
	 *
	 * \tparam id property id
	 * \param vec_key vector element
	 *
	 * \return return the selected property of the vector element
	 *
	 */
	template<unsigned int id> inline auto getPropNC(size_t vec_key) const -> decltype(v_prp.template get<id>(vec_key))
	{
		return v_prp.template get<id>(vec_key);
	}

///////////////////// Read and Write function

	/*! \brief Get the position of an element
	 *
	 * see the vector_dist iterator usage to get an element key
	 *
	 * \param vec_key element
	 *
	 * \return the position of the element in space
	 *
	 */
	inline auto getPosWrite(vect_dist_key_dx vec_key) -> decltype(v_pos.template get<0>(vec_key.getKey()))
	{
#ifdef SE_CLASS3
		se3.template write<prop::max_prop_real>(*this,vec_key.getKey());
#endif

		return v_pos.template get<0>(vec_key.getKey());
	}

	/*! \brief Get the position of an element
	 *
	 * see the vector_dist iterator usage to get an element key
	 *
	 * \param vec_key element
	 *
	 * \return the position of the element in space
	 *
	 */
	inline auto getPosRead(vect_dist_key_dx vec_key) const -> decltype(v_pos.template get<0>(vec_key.getKey()))
	{
#ifdef SE_CLASS3
		se3.template read<prop::max_prop_real>(*this,vec_key.getKey());
#endif

		return v_pos.template get<0>(vec_key.getKey());
	}

	/*! \brief Get the property of an element
	 *
	 * see the vector_dist iterator usage to get an element key
	 *
	 * \tparam id property id
	 * \param vec_key vector element
	 *
	 * \return return the selected property of the vector element
	 *
	 */
	template<unsigned int id> inline auto getPropWrite(vect_dist_key_dx vec_key) -> decltype(v_prp.template get<id>(vec_key.getKey()))
	{
#ifdef SE_CLASS3
		se3.template write<id>(*this,vec_key.getKey());
#endif

		return v_prp.template get<id>(vec_key.getKey());
	}

	/*! \brief Get the property of an element
	 *
	 * see the vector_dist iterator usage to get an element key
	 *
	 * \tparam id property id
	 * \param vec_key vector element
	 *
	 * \return return the selected property of the vector element
	 *
	 */
	template<unsigned int id> inline auto getPropRead(vect_dist_key_dx vec_key) const -> decltype(v_prp.template get<id>(vec_key.getKey()))
	{
#ifdef SE_CLASS3
		se3.template read<id>(*this,vec_key.getKey());
#endif

		return v_prp.template get<id>(vec_key.getKey());
	}

//////////////////////////////////////////////

	/*! \brief Add local particle
	 *
	 * It add a local particle, with "local" we mean in this processor
	 * the particle can be also created out of the processor domain, in this
	 * case a call to map is required. Added particles are always created at the
	 * end and can be accessed with getLastPos and getLastProp
	 *
	 */
	void add()
	{
		v_prp.insert(g_m);
		v_pos.insert(g_m);

		g_m++;

#ifdef SE_CLASS3
		for (size_t i = 0 ; i < prop::max_prop_real+1 ; i++)
			v_prp.template get<prop::max_prop_real>(g_m-1)[i] = UNINITIALIZED;
#endif
	}

#ifndef ONLY_READWRITE_GETTER

	/*! \brief Get the position of the last element
	 *
	 * \return the position of the element in space
	 *
	 */
	inline auto getLastPos() -> decltype(v_pos.template get<0>(0))
	{
		return v_pos.template get<0>(g_m - 1);
	}

	/*! \brief Get the property of the last element
	 *
	 * \tparam id property id
	 *
	 * \return return the selected property of the vector element
	 *
	 */
	template<unsigned int id> inline auto getLastProp() -> decltype(v_prp.template get<id>(0))
	{
		return v_prp.template get<id>(g_m - 1);
	}

#endif

////////////////////////////// READ AND WRITE VERSION //////////

	/*! \brief Get the position of the last element
	 *
	 * \return the position of the element in space
	 *
	 */
	inline auto getLastPosRead() -> decltype(v_pos.template get<0>(0))
	{
#ifdef SE_CLASS3
		se3.template read<prop::max_prop_real>(*this,g_m-1);
#endif

		return v_pos.template get<0>(g_m - 1);
	}

	/*! \brief Get the property of the last element
	 *
	 * \tparam id property id
	 *
	 * \return return the selected property of the vector element
	 *
	 */
	template<unsigned int id> inline auto getLastPropRead() -> decltype(v_prp.template get<id>(0))
	{
#ifdef SE_CLASS3
		se3.read<id>(*this,g_m-1);
#endif

		return v_prp.template get<id>(g_m - 1);
	}


	/*! \brief Get the position of the last element
	 *
	 * \return the position of the element in space
	 *
	 */
	inline auto getLastPosWrite() -> decltype(v_pos.template get<0>(0))
	{
#ifdef SE_CLASS3
		se3.template write<prop::max_prop_real>(*this,g_m-1);
#endif

		return v_pos.template get<0>(g_m - 1);
	}

	/*! \brief Get the property of the last element
	 *
	 * \tparam id property id
	 *
	 * \return return the selected property of the vector element
	 *
	 */
	template<unsigned int id> inline auto getLastPropWrite() -> decltype(v_prp.template get<id>(0))
	{
#ifdef SE_CLASS3
		se3.template write<id>(*this,g_m-1);
#endif

		return v_prp.template get<id>(g_m - 1);
	}

////////////////////////////////////////////////////////////////

	/*! \brief Construct a cell list symmetric based on a cut of radius
	 *
	 * \tparam CellL CellList type to construct
	 *
	 * \param r_cut interation radius, or size of each cell
	 *
	 * \return the Cell list
	 *
	 */
	template<typename CellL = CellList<dim, St, Mem_fast<>, shift<dim, St> > > CellL getCellListSym(St r_cut)
	{
#ifdef SE_CLASS1
		if (!(opt & BIND_DEC_TO_GHOST))
		{
			if (getDecomposition().getGhost().getLow(dim-1) == 0.0)
			{
				std::cerr << __FILE__ << ":" << __LINE__ << " Error the vector has been constructed without BIND_DEC_TO_GHOST, If you construct a vector without BIND_DEC_TO_GHOST the ghost must be full without reductions " << std::endl;
				ACTION_ON_ERROR(VECTOR_DIST_ERROR_OBJECT);
			}
		}
#endif

		// Cell list
		CellL cell_list;

		size_t pad = 0;
		CellDecomposer_sm<dim,St,shift<dim,St>> cd_sm;
		cl_param_calculateSym(getDecomposition().getDomain(),cd_sm,getDecomposition().getGhost(),r_cut,pad);

		// Processor bounding box
		Box<dim, St> pbox = getDecomposition().getProcessorBounds();

		// Ghost padding extension
		Ghost<dim,size_t> g_ext(0);
		cell_list.Initialize(cd_sm,pbox,pad);
		cell_list.set_ndec(getDecomposition().get_ndec());

		updateCellListSym(cell_list);

		return cell_list;
	}

	/*! \brief Construct a cell list symmetric based on a cut of radius
	 *
	 * \tparam CellL CellList type to construct
	 *
	 * \param r_cut interation radius, or size of each cell
	 *
	 * \return the Cell list
	 *
	 */
	template<typename CellL = CellList<dim, St, Mem_fast<>, shift<dim, St> > >
	CellL getCellListSym(const size_t (& div)[dim],
						 const size_t (& pad)[dim])
	{
#ifdef SE_CLASS1
		if (!(opt & BIND_DEC_TO_GHOST))
		{
			if (getDecomposition().getGhost().getLow(dim-1) == 0.0)
			{
				std::cerr << __FILE__ << ":" << __LINE__ << " Error the vector has been constructed without BIND_DEC_TO_GHOST, If you construct a vector without BIND_DEC_TO_GHOST the ghost must be full without reductions " << std::endl;
				ACTION_ON_ERROR(VECTOR_DIST_ERROR_OBJECT);
			}
		}
#endif

		size_t pad_max = pad[0];
		for (size_t i = 1 ; i < dim ; i++)
		{if (pad[i] > pad_max)	{pad_max = pad[i];}}

		// Cell list
		CellL cell_list;

		CellDecomposer_sm<dim,St,shift<dim,St>> cd_sm;
		cd_sm.setDimensions(getDecomposition().getDomain(),div,pad_max);

		// Processor bounding box
		Box<dim, St> pbox = getDecomposition().getProcessorBounds();

		// Ghost padding extension
		Ghost<dim,size_t> g_ext(0);
		cell_list.Initialize(cd_sm,pbox,pad_max);
		cell_list.set_ndec(getDecomposition().get_ndec());

		updateCellListSym(cell_list);

		return cell_list;
	}

	/*! \brief Construct a cell list starting from the stored particles
	 *
	 * \tparam CellL CellList type to construct
	 *
	 * \param r_cut interation radius, or size of each cell
	 * \param no_se3 avoid SE_CLASS3 checking
	 *
	 * \return the Cell list
	 *
	 */
	template<typename CellL = CellList_gen<dim, St, Process_keys_lin, Mem_fast<>, shift<dim, St> > >
	CellL getCellList(St r_cut, bool no_se3 = false)
	{
#ifdef SE_CLASS3
		if (no_se3 == false)
		{se3.getNN();}
#endif
#ifdef SE_CLASS1
		check_ghost_compatible_rcut(r_cut);
#endif

		// Get ghost and anlarge by 1%
		Ghost<dim,St> g = getDecomposition().getGhost();
		g.magnify(1.013);

		return getCellList<CellL>(r_cut, g,no_se3);
	}

#ifdef CUDA_GPU

	/*! \brief Construct a cell list starting from the stored particles, this version of cell-list can be offloaded with
	 *         the function toGPU
	 *
	 * \param r_cut interation radius, or size of each cell
	 *
	 * \return the Cell list
	 *
	 */
	CellList_gpu<dim,St,CudaMemory,shift_only<dim, St>> getCellListGPU(St r_cut, bool no_se3 = false)
	{
#ifdef SE_CLASS3
		if (no_se3 == false)
		{se3.getNN();}
#endif
#ifdef SE_CLASS1
		check_ghost_compatible_rcut(r_cut);
#endif

		// Get ghost and anlarge by 1%
		Ghost<dim,St> g = getDecomposition().getGhost();
		g.magnify(1.013);

		return getCellListGPU(r_cut, g,no_se3);
	}

	/*! \brief Construct a cell list starting from the stored particles
	 *
	 * It differ from the get getCellList for an additional parameter, in case the
	 * domain + ghost is not big enough to contain additional padding particles, a Cell list
	 * with bigger space can be created
	 * (padding particles in general are particles added by the user out of the domains)
	 *
	 * \tparam CellL CellList type to construct
	 *
	 * \param r_cut interation radius, or size of each cell
	 * \param enlarge In case of padding particles the cell list must be enlarged, like a ghost this parameter say how much must be enlarged
	 *
	 * \return the CellList
	 *
	 */
	CellList_gpu<dim,St,CudaMemory,shift_only<dim, St>> getCellListGPU(St r_cut, const Ghost<dim, St> & enlarge, bool no_se3 = false)
	{
#ifdef SE_CLASS3
		if (no_se3 == false)
		{se3.getNN();}
#endif

		// Division array
		size_t div[dim];

		// get the processor bounding box
		Box<dim, St> pbox = getDecomposition().getProcessorBounds();

		// Processor bounding box
		cl_param_calculate(pbox, div, r_cut, enlarge);

		CellList_gpu<dim,St,CudaMemory,shift_only<dim, St>> cell_list(pbox,div);

		v_prp_out.resize(v_pos.size());
		v_pos_out.resize(v_pos.size());

		cell_list.template construct<decltype(v_pos),decltype(v_prp)>(v_pos,v_pos_out,v_prp,v_prp_out);

		return cell_list;
	}

#endif

	/*! \brief Construct an hilbert cell list starting from the stored particles
	 *
	 * \tparam CellL CellList type to construct
	 *
	 * \param r_cut interation radius, or size of each cell
	 *
	 * \return the Cell list
	 *
	 */
	template<typename CellL = CellList_gen<dim, St, Process_keys_hilb, Mem_fast<>, shift<dim, St> > >
	CellL getCellList_hilb(St r_cut)
	{
#ifdef SE_CLASS3
		se3.getNN();
#endif
#ifdef SE_CLASS1
		check_ghost_compatible_rcut(r_cut);
#endif

		// Get ghost and anlarge by 1%
		Ghost<dim,St> g = getDecomposition().getGhost();
		g.magnify(1.013);

		return getCellList_hilb(r_cut, g);
	}

	/*! \brief Update a cell list using the stored particles
	 *
	 * \tparam CellL CellList type to construct
	 *
	 * \param cell_list Cell list to update
	 * \param no_se3 avoid se class 3 checking
	 *
	 */
	template<typename CellL> void updateCellList(CellL & cell_list, bool no_se3 = false)
	{
#ifdef SE_CLASS3
		if (no_se3 == false)
		{se3.getNN();}
#endif

		// This function assume equal spacing in all directions
		// but in the worst case we take the maximum
		St r_cut = 0;
		for (size_t i = 0 ; i < dim ; i++)
		{r_cut = std::max(r_cut,cell_list.getCellBox().getHigh(i));}

		// Here we have to check that the Cell-list has been constructed
		// from the same decomposition
		bool to_reconstruct = cell_list.get_ndec() != getDecomposition().get_ndec();

		if (to_reconstruct == false)
		{
			populate_cell_list(v_pos,cell_list,g_m,CL_NON_SYMMETRIC);

			cell_list.set_gm(g_m);
		}
		else
		{
			CellL cli_tmp = gcl<dim,St,CellL,self,GCL_NON_SYMMETRIC>::get(*this,r_cut,getDecomposition().getGhost());

			cell_list.swap(cli_tmp);
		}
	}

	/*! \brief Update a cell list using the stored particles
	 *
	 * \tparam CellL CellList type to construct
	 *
	 * \param cell_list Cell list to update
	 *
	 */
	template<typename CellL = CellList<dim, St, Mem_fast<>, shift<dim, St> > >
	void updateCellListSym(CellL & cell_list)
	{
#ifdef SE_CLASS3
		se3.getNN();
#endif

		// Here we have to check that the Cell-list has been constructed
		// from the same decomposition
		bool to_reconstruct = cell_list.get_ndec() != getDecomposition().get_ndec();

		if (to_reconstruct == false)
		{
			populate_cell_list(v_pos,cell_list,g_m,CL_SYMMETRIC);

			cell_list.set_gm(g_m);
		}
		else
		{
			CellL cli_tmp = gcl_An<dim,St,CellL,self,GCL_SYMMETRIC>::get(*this,
																		 cell_list.getDivWP(),
																		 cell_list.getPadding(),
																		 getDecomposition().getGhost());

			cell_list.swap(cli_tmp);
		}
	}

	/*! \brief Construct a cell list starting from the stored particles
	 *
	 * It differ from the get getCellList for an additional parameter, in case the
	 * domain + ghost is not big enough to contain additional padding particles, a Cell list
	 * with bigger space can be created
	 * (padding particles in general are particles added by the user out of the domains)
	 *
	 * \tparam CellL CellList type to construct
	 *
	 * \param r_cut interation radius, or size of each cell
	 * \param enlarge In case of padding particles the cell list must be enlarged, like a ghost this parameter say how much must be enlarged
	 * \param no_se3 avoid se_class3 cheking default false
	 *
	 * \return the CellList
	 *
	 */
	template<typename CellL = CellList_gen<dim, St, Process_keys_lin, Mem_fast<>, shift<dim, St> > >
	CellL getCellList(St r_cut, const Ghost<dim, St> & enlarge, bool no_se3 = false)
	{
#ifdef SE_CLASS3
		if (no_se3 == false)
		{se3.getNN();}
#endif

		CellL cell_list;

		// Division array
		size_t div[dim];

		// get the processor bounding box
		Box<dim, St> pbox = getDecomposition().getProcessorBounds();

		// Processor bounding box
		cl_param_calculate(pbox, div, r_cut, enlarge);

		cell_list.Initialize(pbox, div);
		cell_list.set_gm(g_m);
		cell_list.set_ndec(getDecomposition().get_ndec());

		updateCellList(cell_list,no_se3);

		return cell_list;
	}

	/*! \brief Construct an hilbert cell list starting from the stored particles
	 *
	 * It differ from the get getCellList for an additional parameter, in case the
	 * domain + ghost is not big enough to contain additional padding particles, a Cell list
	 * with bigger space can be created
	 * (padding particles in general are particles added by the user out of the domains)
	 *
	 * \tparam CellL CellList type to construct
	 *
	 * \param r_cut interation radius, or size of each cell
	 * \param enlarge In case of padding particles the cell list must be enlarged, like a ghost this parameter say how much must be enlarged
	 *
	 * \return The Cell-list
	 *
	 */
	template<typename CellL = CellList_gen<dim, St, Process_keys_hilb, Mem_fast<>, shift<dim, St> > > CellL getCellList_hilb(St r_cut, const Ghost<dim, St> & enlarge)
	{
#ifdef SE_CLASS3
		se3.getNN();
#endif

		CellL cell_list;

		// Division array
		size_t div[dim];

		// get the processor bounding box
		Box<dim, St> pbox = getDecomposition().getProcessorBounds();

		// Processor bounding box
		cl_param_calculate(pbox,div, r_cut, enlarge);

		cell_list.Initialize(pbox, div);
		cell_list.set_gm(g_m);
		cell_list.set_ndec(getDecomposition().get_ndec());

		updateCellList(cell_list);

		return cell_list;
	}

	/*! \brief for each particle get the symmetric verlet list
	 *
	 * \param r_cut cut-off radius
	 *
	 * \return the verlet list
	 *
	 */
	template <typename VerletL = VerletList<dim,St,Mem_fast<>,shift<dim,St> >>
	VerletL getVerletSym(St r_cut)
	{
#ifdef SE_CLASS3
		se3.getNN();
#endif

		VerletL ver;

		// Processor bounding box
		Box<dim, St> pbox = getDecomposition().getProcessorBounds();

		ver.InitializeSym(getDecomposition().getDomain(),pbox,getDecomposition().getGhost(),r_cut,v_pos,g_m);

		ver.set_ndec(getDecomposition().get_ndec());

		return ver;
	}

	/*! \brief for each particle get the symmetric verlet list
	 *
	 * \param r_cut cut-off radius
	 *
	 * \return the verlet list
	 *
	 */
	template <typename VerletL = VerletList<dim,St,Mem_fast<>,shift<dim,St> >>
	VerletL getVerletCrs(St r_cut)
	{
#ifdef SE_CLASS1
		if (!(opt & BIND_DEC_TO_GHOST))
		{
			std::cerr << __FILE__ << ":" << __LINE__ << " Error the vector has been constructed without BIND_DEC_TO_GHOST, getVerletCrs require the vector to be constructed with BIND_DEC_TO_GHOST option " << std::endl;
			ACTION_ON_ERROR(VECTOR_DIST_ERROR_OBJECT);
		}
#endif

#ifdef SE_CLASS3
		se3.getNN();
#endif

		VerletL ver;

		// Processor bounding box
		Box<dim, St> pbox = getDecomposition().getProcessorBounds();

		// Initialize the verlet list
		ver.InitializeCrs(getDecomposition().getDomain(),pbox,getDecomposition().getGhost(),r_cut,v_pos,g_m);

		// Get the internal cell list
		auto & NN = ver.getInternalCellList();

		// Shift
		grid_key_dx<dim> shift;

		// Add padding
		for (size_t i = 0 ; i < dim ; i++)
			shift.set_d(i,NN.getPadding(i));

		grid_sm<dim,void> gs = NN.getInternalGrid();

		getDecomposition().setNNParameters(shift,gs);

		ver.createVerletCrs(r_cut,g_m,v_pos,
				            getDecomposition().getCRSDomainCells(),
							getDecomposition().getCRSAnomDomainCells());

		ver.set_ndec(getDecomposition().get_ndec());

		return ver;
	}

	/*! \brief for each particle get the verlet list
	 *
	 * \param r_cut cut-off radius
	 *
	 * \return a VerletList object
	 *
	 */
	template <typename VerletL = VerletList<dim,St,Mem_fast<>,shift<dim,St> >>
	VerletL getVerlet(St r_cut)
	{
#ifdef SE_CLASS3
		se3.getNN();
#endif

		VerletL ver;

		// get the processor bounding box
		Box<dim, St> bt = getDecomposition().getProcessorBounds();

		// Get the ghost
		Ghost<dim,St> g = getDecomposition().getGhost();
		g.magnify(1.013);

		// enlarge the box where the Verlet is defined
		bt.enlarge(g);

		ver.Initialize(bt,getDecomposition().getProcessorBounds(),r_cut,v_pos,g_m,VL_NON_SYMMETRIC);

		ver.set_ndec(getDecomposition().get_ndec());

		return ver;
	}

	/*! \brief for each particle get the verlet list
	 *
	 * \param r_cut cut-off radius
	 * \param ver Verlet to update
	 * \param r_cut cutoff radius
	 * \param opt option like VL_SYMMETRIC and VL_NON_SYMMETRIC or VL_CRS_SYMMETRIC
	 *
	 */
	template<typename Mem_type> void updateVerlet(VerletList<dim,St,Mem_type,shift<dim,St> > & ver, St r_cut, size_t opt = VL_NON_SYMMETRIC)
	{
#ifdef SE_CLASS3
		se3.getNN();
#endif

		if (opt == VL_SYMMETRIC)
		{
			auto & NN = ver.getInternalCellList();

			// Here we have to check that the Box defined by the Cell-list is the same as the domain box of this
			// processor. if it is not like that we have to completely reconstruct from stratch
			bool to_reconstruct = NN.get_ndec() != getDecomposition().get_ndec();

			if (to_reconstruct == false)
				ver.update(getDecomposition().getDomain(),r_cut,v_pos,g_m, opt);
			else
			{
				VerletList<dim,St,Mem_type,shift<dim,St> > ver_tmp;

				ver_tmp = getVerlet<VerletList<dim,St,Mem_type,shift<dim,St> >>(r_cut);
				ver.swap(ver);
			}
		}
		else if (opt == VL_CRS_SYMMETRIC)
		{
#ifdef SE_CLASS1
			if ((opt & BIND_DEC_TO_GHOST))
			{
				std::cerr << __FILE__ << ":" << __LINE__ << " Error the vector has been constructed without BIND_DEC_TO_GHOST, updateVerlet with the option VL_CRS_SYMMETRIC require the vector to be constructed with BIND_DEC_TO_GHOST option " << std::endl;
				ACTION_ON_ERROR(VECTOR_DIST_ERROR_OBJECT);
			}
#endif

			auto & NN = ver.getInternalCellList();

			// Here we have to check that the Box defined by the Cell-list is the same as the domain box of this
			// processor. if it is not like that we have to completely reconstruct from stratch
			bool to_reconstruct = NN.get_ndec() != getDecomposition().get_ndec();

			if (to_reconstruct == false)
			{
				// Shift
				grid_key_dx<dim> shift;

				// Add padding
				for (size_t i = 0 ; i < dim ; i++)
					shift.set_d(i,NN.getPadding(i));

				grid_sm<dim,void> gs = NN.getInternalGrid();

				getDecomposition().setNNParameters(shift,gs);

				ver.updateCrs(getDecomposition().getDomain(),r_cut,v_pos,g_m,
						      getDecomposition().getCRSDomainCells(),
							  getDecomposition().getCRSAnomDomainCells());
			}
			else
			{
				VerletList<dim,St,Mem_type,shift<dim,St> > ver_tmp;

				ver_tmp = getVerletCrs<VerletList<dim,St,Mem_type,shift<dim,St> >>(r_cut);
				ver.swap(ver_tmp);
			}
		}
		else
		{
			auto & NN = ver.getInternalCellList();

			// Here we have to check that the Box defined by the Cell-list is the same as the domain box of this
			// processor. if it is not like that we have to completely reconstruct from stratch
			bool to_reconstruct = NN.get_ndec() != getDecomposition().get_ndec();

			if (to_reconstruct == false)
				ver.update(getDecomposition().getDomain(),r_cut,v_pos,g_m, opt);
			else
			{
				VerletList<dim,St,Mem_type,shift<dim,St> > ver_tmp;

				ver_tmp = getVerlet<VerletList<dim,St,Mem_type,shift<dim,St> >>(r_cut);
				ver.swap(ver_tmp);
			}
		}
	}


	/*! \brief Construct a cell list starting from the stored particles and reorder a vector according to the Hilberts curve
	 *
	 * \tparam CellL CellList type to construct
	 *
	 * \param m an order of a hilbert curve
	 *
	 */
	template<typename CellL=CellList_gen<dim,St,Process_keys_lin,Mem_bal<>,shift<dim,St> > >
	void reorder (int32_t m, reorder_opt opt = reorder_opt::HILBERT)
	{
		reorder<CellL>(m,getDecomposition().getGhost(),opt);
	}


	/*! \brief Construct a cell list starting from the stored particles and reorder a vector according to the Hilberts curve
	 *
	 *
	 *It differs from the reorder(m) for an additional parameter, in case the
	 * domain + ghost is not big enough to contain additional padding particles, a Cell list
	 * with bigger space can be created
	 * (padding particles in general are particles added by the user out of the domains)
	 *
	 * \param m order of a curve
	 * \param enlarge In case of padding particles the cell list must be enlarged, like a ghost this parameter say how much must be enlarged
	 *
	 */
	template<typename CellL=CellList_gen<dim,St,Process_keys_lin,Mem_bal<>,shift<dim,St> > >
	void reorder(int32_t m, const Ghost<dim,St> & enlarge, reorder_opt opt = reorder_opt::HILBERT)
	{
		// reset the ghost part
		v_pos.resize(g_m);
		v_prp.resize(g_m);


		CellL cell_list;

		// calculate the parameters of the cell list

		// get the processor bounding box
		Box<dim,St> pbox = getDecomposition().getProcessorBounds();
		// extend by the ghost
		pbox.enlarge(enlarge);

		size_t div[dim];

		// Calculate the division array and the cell box
		for (size_t i = 0 ; i < dim ; i++)
		{
			div[i] = 1 << m;
		}

		cell_list.Initialize(pbox,div);
		cell_list.set_gm(g_m);

		// for each particle add the particle to the cell list

		auto it = getIterator();

		while (it.isNext())
		{
			auto key = it.get();

			Point<dim,St> xp = this->getPos(key);

			cell_list.add(xp,key.getKey());

			++it;
		}

		// Use cell_list to reorder v_pos

		//destination vector
		openfpm::vector<Point<dim,St>> v_pos_dest;
		openfpm::vector<prop> v_prp_dest;

		if (opt == reorder_opt::HILBERT)
		{
			grid_key_dx_iterator_hilbert<dim> h_it(m);

			reorder_sfc<CellL,grid_key_dx_iterator_hilbert<dim>>(v_pos_dest,v_prp_dest,h_it,cell_list);
		}
		else if (opt == reorder_opt::LINEAR)
		{
			grid_sm<dim,void> gs(div);
			grid_key_dx_iterator<dim> h_it(gs);

			reorder_sfc<CellL,grid_key_dx_iterator<dim>>(v_pos_dest,v_prp_dest,h_it,cell_list);
		}
		else
		{
			// We do nothing, we second swap nullify the first
			v_pos.swap(v_pos_dest);
			v_prp.swap(v_prp_dest);
		}

		v_pos.swap(v_pos_dest);
		v_prp.swap(v_prp_dest);
	}

	/*! \brief It return the number of particles contained by the previous processors
	 *
	 * \warning It only work with the initial decomposition
	 *
	 * Given 1000 particles and 3 processors, you will get
	 *
	 * * Processor 0: 0
	 * * Processor 1: 334
	 * * Processor 2: 667
	 *
	 * \param np initial number of particles
	 *
	 * \return number of particles contained by the previous processors
	 *
	 */
	size_t init_size_accum(size_t np)
	{
		size_t accum = 0;

		// convert to a local number of elements
		size_t p_np = np / v_cl.getProcessingUnits();

		// Get non divisible part
		size_t r = np % v_cl.getProcessingUnits();

		accum = p_np * v_cl.getProcessUnitID();

		// Distribute the remain particles
		if (v_cl.getProcessUnitID() <= r)
			accum += v_cl.getProcessUnitID();
		else
			accum += r;

		return accum;
	}

	/*! \brief Get an iterator that traverse domain and ghost particles
	 *
	 * \return an iterator
	 *
	 */
	vector_dist_iterator getIterator()
	{
#ifdef SE_CLASS3
		se3.getIterator();
#endif
		return vector_dist_iterator(0, v_pos.size());
	}

	/*! /brief Get a grid Iterator
	 *
	 * Usefull function to place particles on a grid or grid-like (grid + noise)
	 *
	 * \param sz size of the grid
	 *
	 * \return a Grid iterator
	 *
	 */
	inline grid_dist_id_iterator_dec<Decomposition> getGridIterator(const size_t (&sz)[dim])
	{
		grid_key_dx<dim> start;
		grid_key_dx<dim> stop;
		for (size_t i = 0; i < dim; i++)
		{
			start.set_d(i, 0);
			stop.set_d(i, sz[i] - 1);
		}

		grid_dist_id_iterator_dec<Decomposition> it_dec(getDecomposition(), sz, start, stop);
		return it_dec;
	}

	/*! \brief Get the iterator across the position of the ghost particles
	 *
	 * \return an iterator
	 *
	 */
	vector_dist_iterator getGhostIterator() const
	{
#ifdef SE_CLASS3
		se3.getIterator();
#endif

		return vector_dist_iterator(g_m, v_pos.size());
	}

	/*! \brief Get the iterator across the position of the ghost particles
	 *
	 * \return an iterator
	 *
	 */
	vector_dist_iterator getGhostIterator_no_se3() const
	{
		return vector_dist_iterator(g_m, v_pos.size());
	}

	/*! \brief Get an iterator that traverse the particles in the domain
	 *         using a cell list
	 *
	 * \param NN Cell-list
	 *
	 * \return an iterator over the particles
	 *
	 */
	template<typename CellList> ParticleIt_Cells<dim,CellList>
	getDomainIteratorCells(CellList & NN)
	{
#ifdef SE_CLASS3
		se3.getIterator();
#endif

		// Shift
		grid_key_dx<dim> shift;

		// Add padding
		for (size_t i = 0 ; i < dim ; i++)
			shift.set_d(i,NN.getPadding(i));

		grid_sm<dim,void> gs = NN.getInternalGrid();

		getDecomposition().setNNParameters(shift,gs);

		return ParticleIt_Cells<dim,CellList>(NN,getDecomposition().getDomainCells(),g_m);
	}

	/*! \brief Get an iterator that traverse the particles in the domain
	 *
	 * \return an iterator
	 *
	 */
	vector_dist_iterator getDomainIterator() const
	{
#ifdef SE_CLASS3
		se3.getIterator();
#endif

		return vector_dist_iterator(0, g_m);
	}

#ifdef CUDA_GPU

	/*! \brief Get an iterator that traverse the particles in the domain
	 *
	 * \return an iterator
	 *
	 */
	ite_gpu<1> getDomainIteratorGPU(size_t n_thr = 1024) const
	{
#ifdef SE_CLASS3
		se3.getIterator();
#endif

		return v_pos.getGPUIteratorTo(n_thr);
	}

#endif

	/*! \brief Get an iterator that traverse the particles in the domain
	 *
	 * \return an iterator
	 *
	 */
	vector_dist_iterator getDomainIterator_no_se3() const
	{
		return vector_dist_iterator(0, g_m);
	}

	/*! \brief Get an iterator that traverse the particles in the domain
	 *
	 * \return an iterator
	 *
	 */
	vector_dist_iterator getDomainAndGhostIterator() const
	{
#ifdef SE_CLASS3
		se3.getIterator();
#endif

		return vector_dist_iterator(0, v_pos.size());
	}

	/*! \brief Get an iterator that traverse the particles in the domain
	 *
	 * \return an iterator
	 *
	 */
	vector_dist_iterator getDomainAndGhostIterator_no_se3() const
	{
		return vector_dist_iterator(0, v_pos.size());
	}

	/*! \brief Get the decomposition
	 *
	 * \return
	 *
	 */
	inline Decomposition & getDecomposition()
	{
		return vector_dist_comm<dim,St,prop,Decomposition,Memory,layout_base>::getDecomposition();
	}

	/*! \brief Get the decomposition
	 *
	 * \return
	 *
	 */
	inline const Decomposition & getDecomposition() const
	{
		return vector_dist_comm<dim,St,prop,Decomposition,Memory,layout_base>::getDecomposition();
	}

	/*! \brief It move all the particles that does not belong to the local processor to the respective processor
	 *
	 * \tparam out of bound policy it specify what to do when the particles are detected out of bound
	 *
	 * In general this function is called after moving the particles to move the
	 * elements out the local processor. Or just after initialization if each processor
	 * contain non local particles
	 *
	 * \tparam prp properties to communicate
	 *
	 * \param opt options
	 *
	 */
	template<unsigned int ... prp> void map_list(size_t opt = NONE)
	{
#ifdef SE_CLASS3
		se3.map_pre();
#endif

		this->template map_list_<prp...>(v_pos,v_prp,g_m,opt);

#ifdef SE_CLASS3
		se3.map_post();
#endif
	}


	/*! \brief It move all the particles that does not belong to the local processor to the respective processor
	 *
	 * \tparam out of bound policy it specify what to do when the particles are detected out of bound
	 *
	 * In general this function is called after moving the particles to move the
	 * elements out the local processor. Or just after initialization if each processor
	 * contain non local particles
	 *
	 * \param opt options
	 *
	 */
	template<typename obp = KillParticle> void map(size_t opt = NONE)
	{
#ifdef SE_CLASS3
		se3.map_pre();
#endif

		this->template map_<obp>(v_pos,v_prp,g_m,opt);

#ifdef SE_CLASS3
		se3.map_post();
#endif
	}

	/*! \brief It synchronize the properties and position of the ghost particles
	 *
	 * \tparam prp list of properties to get synchronize
	 *
	 * \param opt options WITH_POSITION, it send also the positional information of the particles
	 *
	 */
	template<int ... prp> inline void ghost_get(size_t opt = WITH_POSITION)
	{
#ifdef SE_CLASS1
		if (getDecomposition().getProcessorBounds().isValid() == false && size_local() != 0)
		{
			std::cerr << __FILE__ << ":" << __LINE__ << " Error the processor " << v_cl.getProcessUnitID() << " has particles, but is supposed to be unloaded" << std::endl;
			ACTION_ON_ERROR(VECTOR_DIST_ERROR_OBJECT);
		}
#endif

#ifdef SE_CLASS3
		se3.template ghost_get_pre<prp...>(opt);
#endif

		this->template ghost_get_<prp...>(v_pos,v_prp,g_m,opt);

#ifdef SE_CLASS3

		this->template ghost_get_<prop::max_prop_real>(v_pos,v_prp,g_m,opt | KEEP_PROPERTIES);

		se3.template ghost_get_post<prp...>(opt);
#endif
	}

	/*! \brief It synchronize the properties and position of the ghost particles
	 *
	 * \tparam op which kind of operation to apply
	 * \tparam prp list of properties to get synchronize
	 *
	 * \param opt_ options. It is an optional parameter.
	 *             If set to NO_CHANGE_ELEMENTS the communication has lower latencies. This option has some usage limitations, so please refere to the examples
	 *             for further explanations
	 *
	 *
	 */
	template<template<typename,typename> class op, int ... prp> inline void ghost_put(size_t opt_ = NONE)
	{
#ifdef SE_CLASS3
		se3.template ghost_put<prp...>();
#endif
		this->template ghost_put_<op,prp...>(v_pos,v_prp,g_m,opt_);
	}

	/*! \brief Remove a set of elements from the distributed vector
	 *
	 * \warning keys must be sorted
	 *
	 * \param keys vector of elements to eliminate
	 * \param start from where to eliminate
	 *
	 */
	void remove(openfpm::vector<size_t> & keys, size_t start = 0)
	{
		v_pos.remove(keys, start);
		v_prp.remove(keys, start);

		g_m -= keys.size();
	}

	/*! \brief Remove one element from the distributed vector
	 *
	 * \param key remove one element from the vector
	 *
	 */
	void remove(size_t key)
	{
		v_pos.remove(key);
		v_prp.remove(key);

		g_m--;
	}

	/*! \brief Add the computation cost on the decomposition coming
	 * from the particles
	 *
	 * \param md Model to use
	 * \param vd external vector to add for the computational cost
	 *
	 */
	template <typename Model=ModelLin>inline void addComputationCosts(const self & vd, Model md=Model())
	{
		CellDecomposer_sm<dim, St, shift<dim,St>> cdsm;

		Decomposition & dec = getDecomposition();

		cdsm.setDimensions(dec.getDomain(), dec.getDistGrid().getSize(), 0);

		auto it = vd.getDomainIterator();

		while (it.isNext())
		{
			size_t v = cdsm.getCell(vd.getPos(it.get()));

			md.addComputation(dec,vd,v,it.get().getKey());

			++it;
		}
	}

	/*! \brief Add the computation cost on the decomposition coming
	 * from the particles
	 *
	 * \param md Model to use
	 * \param ts It is an optional parameter approximately should be the number of ghost get between two
	 *           rebalancing at first decomposition this number can be ignored (default = 1) because not used
	 *
	 */
	template <typename Model=ModelLin> void finalizeComputationCosts(Model md=Model(), size_t ts = 1)
	{
		Decomposition & dec = getDecomposition();
		auto & dist = getDecomposition().getDistribution();

		dec.computeCommunicationAndMigrationCosts(ts);

		// Go throught all the sub-sub-domains and apply the model

		for (size_t i = 0 ; i < dist.getNOwnerSubSubDomains(); i++)
		{md.applyModel(dec,dist.getOwnerSubSubDomain(i));}

		dist.setDistTol(md.distributionTol());
	}

	/*! \brief Initialize the computational cost
	 *
	 */
	void initializeComputationCosts()
	{
		Decomposition & dec = getDecomposition();
		auto & dist = getDecomposition().getDistribution();

		for (size_t i = 0; i < dist.getNOwnerSubSubDomains() ; i++)
		{dec.setSubSubDomainComputationCost(dist.getOwnerSubSubDomain(i) , 1);}
	}

	/*! \brief Add the computation cost on the decomposition coming
	 * from the particles
	 *
	 * \param md Model to use
	 * \param ts It is an optional parameter approximately should be the number of ghost get between two
	 *           rebalancing at first decomposition this number can be ignored (default = 1) because not used
	 *
	 */
	template <typename Model=ModelLin>inline void addComputationCosts(Model md=Model(), size_t ts = 1)
	{
		initializeComputationCosts();

		addComputationCosts(*this,md);

		finalizeComputationCosts(md,ts);
	}

	/*! \brief Save the distributed vector on HDF5 file
	 *
	 * \param filename file where to save
	 *
	 */
	inline void save(const std::string & filename) const
	{
		HDF5_writer<VECTOR_DIST> h5s;

		h5s.save(filename,v_pos,v_prp);
	}

	/*! \brief Load the distributed vector from an HDF5 file
	 *
	 * \param filename file from where to load
	 *
	 */
	inline void load(const std::string & filename)
	{
		HDF5_reader<VECTOR_DIST> h5l;

		h5l.load(filename,v_pos,v_prp,g_m);
	}

	/*! \brief Output particle position and properties
	 *
	 * \param out output filename
	 * \param opt VTK_WRITER, CSV_WRITER, it is also possible to choose the format for  VTK
	 *            FORMAT_BINARY. (the default is ASCII format)
	 *
	 * \return true if the file has been written without error
	 *
	 */
	inline bool write(std::string out, int opt = VTK_WRITER)
	{

		if ((opt & 0x0FFF0000) == CSV_WRITER)
		{
			// CSVWriter test
			CSVWriter<openfpm::vector<Point<dim, St>,Memory,typename layout_base<Point<dim,St>>::type,layout_base>,
			          openfpm::vector<prop,Memory,typename layout_base<prop>::type,layout_base> > csv_writer;

			std::string output = std::to_string(out + "_" + std::to_string(v_cl.getProcessUnitID()) + std::to_string(".csv"));

			// Write the CSV
			return csv_writer.write(output,v_pos,v_prp);
		}
		else
		{
			file_type ft = file_type::ASCII;

			if (opt & FORMAT_BINARY)
				ft = file_type::BINARY;

			// VTKWriter for a set of points
			VTKWriter<boost::mpl::pair<openfpm::vector<Point<dim, St>,Memory,typename layout_base<Point<dim,St>>::type,layout_base>,
									   openfpm::vector<prop,Memory,typename layout_base<prop>::type,layout_base>>,
			                           VECTOR_POINTS> vtk_writer;
			vtk_writer.add(v_pos,v_prp,g_m);

			std::string output = std::to_string(out + "_" + std::to_string(v_cl.getProcessUnitID()) + std::to_string(".vtk"));

			// Write the VTK file
			return vtk_writer.write(output,prp_names,"particles",ft);
		}
	}

	/*! \brief Delete the particles on the ghost
	 *
	 *
	 */
	void deleteGhost()
	{
		v_pos.resize(g_m);
		v_prp.resize(g_m);
	}

	/*! \brief Resize the vector (locally)
	 *
	 * \warning It automatically delete the ghosts
	 *
	 * \param rs
	 *
	 */
	void resize(size_t rs)
	{
		deleteGhost();

		v_pos.resize(rs);
		v_prp.resize(rs);

		g_m = rs;
	}

	/*! \brief Output particle position and properties
	 *
	 * \param out output
	 * \param iteration (we can append the number at the end of the file_name)
	 * \param opt VTK_WRITER, CSV_WRITER, it is also possible to choose the format for  VTK
	 *            FORMAT_BINARY. (the default is ASCII format)
	 *
	 * \return if the file has been written correctly
	 *
	 */
	inline bool write_frame(std::string out, size_t iteration, int opt = VTK_WRITER)
	{
		if ((opt & 0x0FFF0000) == CSV_WRITER)
		{
			// CSVWriter test
			CSVWriter<openfpm::vector<Point<dim, St>>, openfpm::vector<prop> > csv_writer;

			std::string output = std::to_string(out + "_" + std::to_string(v_cl.getProcessUnitID()) + "_" + std::to_string(iteration) + std::to_string(".csv"));

			// Write the CSV
			return csv_writer.write(output, v_pos, v_prp);
		}
		else
		{
			file_type ft = file_type::ASCII;

			if (opt & FORMAT_BINARY)
				ft = file_type::BINARY;

			// VTKWriter for a set of points
			VTKWriter<boost::mpl::pair<openfpm::vector<Point<dim,St>>, openfpm::vector<prop>>, VECTOR_POINTS> vtk_writer;
			vtk_writer.add(v_pos,v_prp,g_m);

			std::string output = std::to_string(out + "_" + std::to_string(v_cl.getProcessUnitID()) + "_" + std::to_string(iteration) + std::to_string(".vtk"));

			// Write the VTK file
			return vtk_writer.write(output,prp_names,"particles",ft);
		}
	}

	/*! \brief Get the Celllist parameters
	 *
	 * \param r_cut spacing of the cell-list
	 * \param div division required for the cell-list
	 * \param box where the Cell list must be defined (In general Processor domain + Ghost)
	 * \param enlarge Optionally a request to make the space a littler bit larger than Processor domain + Ghost
	 *        keeping the cell list consistent with the requests
	 *
	 */
	void getCellListParams(St r_cut, size_t (&div)[dim],Box<dim, St> & box, Ghost<dim,St> enlarge = Ghost<dim,St>(0.0))
	{
		// get the processor bounding box
		Box<dim, St> pbox = getDecomposition().getProcessorBounds();

		// enlarge the processor bounding box by the ghost
		Ghost<dim,St> g = getDecomposition().getGhost();
		pbox.enlarge(g);

		cl_param_calculate(pbox, div,r_cut,enlarge);

		// output the fixed domain
		box = pbox;
	}

	/*! \brief It return the id of structure in the allocation list
	 *
	 * \see print_alloc and SE_CLASS2
	 *
	 * \return the id
	 *
	 */
	long int who()
	{
#ifdef SE_CLASS2
		return check_whoami(this,8);
#else
		return -1;
#endif
	}

	/*! \brief Get the Virtual Cluster machine
	 *
	 * \return the Virtual cluster machine
	 *
	 */

	Vcluster<> & getVC()
	{
#ifdef SE_CLASS2
		check_valid(this,8);
#endif
		return v_cl;
	}

	/*! \brief return the position vector of all the particles
	 *
	 * \return the particle position vector
	 *
	 */
	const openfpm::vector<Point<dim,St>> & getPosVector() const
	{
		return v_pos;
	}

	/*! \brief return the position vector of all the particles
	 *
	 * \return the particle position vector
	 *
	 */
	openfpm::vector<Point<dim,St>> & getPosVector()
	{
		return v_pos;
	}

	/*! \brief return the property vector of all the particles
	 *
	 * \return the particle property vector
	 *
	 */
	const openfpm::vector<prop> & getPropVector() const
	{
		return v_prp;
	}

	/*! \brief return the property vector of all the particles
	 *
	 * \return the particle property vector
	 *
	 */
	openfpm::vector<prop> & getPropVector()
	{
		return v_prp;
	}

	/*! \brief It return the sum of the particles in the previous processors
	 *
	 * \return the particles number
	 *
	 */
	size_t accum()
	{
		openfpm::vector<size_t> accu;

		size_t sz = size_local();

		v_cl.allGather(sz,accu);
		v_cl.execute();

		sz = 0;

		for (size_t i = 0 ; i < v_cl.getProcessUnitID() ; i++)
		{sz += accu.get(i);}

		return sz;
	}

	/*! \brief Get a special particle iterator able to iterate across particles using
	 *         symmetric crossing scheme
	 *
	 * \param NN Cell-List neighborhood
	 *
	 * \return Particle iterator
	 *
	 */
	template<typename cli> ParticleItCRS_Cells<dim,cli> getParticleIteratorCRS_Cell(cli & NN)
	{
#ifdef SE_CLASS3
		se3.getIterator();
#endif

#ifdef SE_CLASS1
		if (!(opt & BIND_DEC_TO_GHOST))
		{
			std::cerr << __FILE__ << ":" << __LINE__ << " Error the vector has been constructed without BIND_DEC_TO_GHOST, getParticleIteratorCRS_Cell require the vector to be constructed with BIND_DEC_TO_GHOST option " << std::endl;
			ACTION_ON_ERROR(VECTOR_DIST_ERROR_OBJECT);
		}
#endif

		// Shift
		grid_key_dx<dim> shift;

		// Add padding
		for (size_t i = 0 ; i < dim ; i++)
			shift.set_d(i,NN.getPadding(i));

		grid_sm<dim,void> gs = NN.getInternalGrid();

		getDecomposition().setNNParameters(shift,gs);

		// First we check that
		return ParticleItCRS_Cells<dim,cli>(NN,getDecomposition().getCRSDomainCells(),
				                               getDecomposition().getCRSAnomDomainCells(),
											   NN.getNNc_sym());
	}

	/*! \brief Set the properties names
	 *
	 * It is useful to specify name for the properties in vtk writers
	 *
	 * \param names set of properties names
	 *
	 */
	void setPropNames(const openfpm::vector<std::string> & names)
	{
		prp_names = names;
	}

	/*! \brief Get a special particle iterator able to iterate across particles using
	 *         symmetric crossing scheme
	 *
	 * \param NN Verlet list neighborhood
	 *
	 * \return Particle iterator
	 *
	 */
	template<typename vrl> openfpm::vector_key_iterator_seq<typename vrl::Mem_type_type::local_index_type> getParticleIteratorCRS(vrl & NN)
	{
#ifdef SE_CLASS1
		if (!(opt & BIND_DEC_TO_GHOST))
		{
			std::cerr << __FILE__ << ":" << __LINE__ << " Error the vector has been constructed without BIND_DEC_TO_GHOST, getParticleIteratorCRS_Cell require the vector to be constructed with BIND_DEC_TO_GHOST option " << std::endl;
			ACTION_ON_ERROR(VECTOR_DIST_ERROR_OBJECT);
		}
#endif

		// First we check that
		return openfpm::vector_key_iterator_seq<typename vrl::Mem_type_type::local_index_type>(NN.getParticleSeq());
	}

	/*! \brief Return from which cell we have to start in case of CRS interation
	 *         scheme
	 *
	 * \param NN cell-list
	 *
	 * \return The starting cell point
	 *
	 */
	template<typename Celllist> grid_key_dx<dim> getCRSStart(Celllist & NN)
	{
		return NN.getStartDomainCell();
	}

	/*! \brief Return from which cell we have to stop in case of CRS interation
	 *         scheme
	 *
	 * \param NN cell-list
	 *
	 * \return The stop cell point
	 *
	 */
	template<typename Celllist> grid_key_dx<dim> getCRSStop(Celllist & NN)
	{
		grid_key_dx<dim> key = NN.getStopDomainCell();

		for (size_t i = 0 ; i < dim ; i++)
			key.set_d(i,key.get(i) + 1);
		return key;
	}


#ifdef CUDA_GPU

		/*! \brief Convert the grid into a data-structure compatible for computing into GPU
		 *
		 *  The object created can be considered like a reference of the original
		 *
		 * \return an usable vector in the kernel
		 *
		 */
		template<unsigned int ... prp> vector_dist_ker<dim,St,prop> toKernel()
		{
			vector_dist_ker<dim,St,prop> v(v_pos.toKernel(), v_prp.toKernel());

			return v;
		}

		/*! \brief Convert the grid into a data-structure compatible for computing into GPU
		 *
		 *  In comparison with toGPU return a version sorted better for coalesced memory
		 *
		 * \return an usable vector in the kernel
		 *
		 */
		template<unsigned int ... prp> vector_dist_ker<dim,St,prop> toKernel_sorted()
		{
			vector_dist_ker<dim,St,prop> v(v_pos_out.toKernel(), v_prp_out.toKernel());

			return v;
		}

		/*! \brief Move the memory from the device to host memory
		 *
		 * \tparam property to move use POS_PROP for position property
		 *
		 */
		template<unsigned int ... prp> void deviceToHostProp()
		{
			v_prp.template deviceToHost<prp ...>();
		}

		/*! \brief Move the memory from the device to host memory
		 *
		 * \tparam property to move use POS_PROP for position property
		 *
		 */
		void deviceToHostPos()
		{
			v_pos.template deviceToHost<0>();
		}

		/*! \brief Move the memory from the device to host memory
		 *
		 * \tparam property to move use POS_PROP for position property
		 *
		 */
		template<unsigned int ... prp> void hostToDeviceProp()
		{
			v_prp.template hostToDevice<prp ...>();
		}

		/*! \brief Move the memory from the device to host memory
		 *
		 * \tparam property to move use POS_PROP for position property
		 *
		 */
		void hostToDevicePos()
		{
			v_pos.template hostToDevice<0>();
		}

#endif


#ifdef SE_CLASS3

	se_class3_vector<prop::max_prop,dim,St,Decomposition,self> & get_se_class3()
	{
		return se3;
	}

#endif
};


template<unsigned int dim, typename St, typename prop, typename Decomposition = CartDecomposition<dim,St,CudaMemory,memory_traits_inte>> using vector_dist_gpu = vector_dist<dim,St,prop,Decomposition,CudaMemory,memory_traits_inte>;
template<unsigned int dim, typename St, typename prop, typename Decomposition = CartDecomposition<dim,St,HeapMemory,memory_traits_inte>> using vector_dist_soa = vector_dist<dim,St,prop,Decomposition,HeapMemory,memory_traits_inte>;

#endif /* VECTOR_HPP_ */
