/*
 * Vector.hpp
 *
 *  Created on: Mar 5, 2015
 *      Author: Pietro Incardona
 */

#ifndef VECTOR_HPP_
#define VECTOR_HPP_

#include "config.h"
#include "util/cuda_util.hpp"
#include "HDF5_wr/HDF5_wr.hpp"
#include "VCluster/VCluster.hpp"
#include "Space/Shape/Point.hpp"
#include "Vector/Iterators/vector_dist_iterator.hpp"
#include "Space/Shape/Box.hpp"
#include "Vector/vector_dist_key.hpp"
#include "memory/PtrMemory.hpp"
#include "NN/CellList/CellList.hpp"
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
#include "NN/CellList/SFCKeys.hpp"
#include "Vector/vector_dist_kernel.hpp"
#include "NN/CellList/cuda/CellList_gpu.hpp"
#include "lib/pdata.hpp"
#include "cuda/vector_dist_operators_list_ker.hpp"
#include "util/PathsAndFiles.hpp"

#include <type_traits>

#define DEC_GRAN(gr) ((size_t)gr << 32)

#ifdef CUDA_GPU
template<unsigned int dim,typename St> using CELLLIST_GPU_SPARSE = CellList_gpu<dim,St,CudaMemory,shift_only<dim, St>,true>;
#endif

#define VECTOR_DIST_ERROR_OBJECT std::runtime_error("Runtime vector distributed error");

#ifdef SE_CLASS3
#include "se_class3_vector.hpp"
#endif

#ifdef SE_CLASS3
	#define SE_CLASS3_VDIST_CONSTRUCTOR ,se3(getDecomposition(),*this)
#else
	#define SE_CLASS3_VDIST_CONSTRUCTOR
#endif

// Perform a ghost get or a ghost put
constexpr int GET = 1;
constexpr int PUT = 2;

// Write the particles with ghost
constexpr int NO_GHOST = 0;
constexpr int WITH_GHOST = 2;

enum reorder_opt
{
	NO_REORDER = 0,
	HILBERT = 1,
	LINEAR = 2
};


template<typename vector_type>
struct decrement_memory
{
	//! vector
	vector_type & v;


	/*! \brief constructor
	 *
	 *
	 * \param src encapsulated object1
	 * \param dst encapsulated object2
	 *
	 */
	inline decrement_memory(vector_type & v)
	:v(v)
	{};


	/*!  \brief It call the copy function for each property
	 *
	 * \param t each member
	 *
	 */
	template<typename T>
	inline void operator()(T& t) const
	{
		v.template getMemory<T::value>().decRef();
	}
};

/*! \brief Distributed vector
 *
 * This class represent a distributed vector, the distribution of the structure
 * is based on the positional information of the elements the vector store
 *
 * ## Create a vector of random elements on each processor 2D
 * \snippet Vector/tests/vector_dist_unit_test.cpp Create a vector of random elements on each processor 2D
 *
 * ## Create a vector of random elements on each processor 3D
 * \snippet Vector/tests/vector_dist_unit_test.cpp Create a vector of random elements on each processor 3D
 *
 * ## Create a vector of elements distributed on a grid like way
 * \snippet Vector/tests/vector_dist_unit_test.cpp Create a vector of elements distributed on a grid like way
 *
 * ## Redistribute the particles and sync the ghost properties
 * \snippet Vector/tests/vector_dist_unit_test.cpp Redistribute the particles and sync the ghost properties
 *
 * ## Create a gpu distributed vector [St = float or double]
 * \snippet Vector/cuda/vector_dist_gpu_unit_tests.cu Create a gpu vector
 *
 * ## Fill a GPU vector_dist on CPU and move the information to GPU and redistribute [St = float or double]
 * \snippet Vector/cuda/vector_dist_gpu_unit_tests.cu Fill gpu vector and move to GPU
 *
 * ## Fill the ghost on GPU
 * \snippet Vector/cuda/vector_dist_gpu_unit_tests.cu Fill the ghost on GPU
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
         template<typename> class layout_base = memory_traits_lin,
		 typename vector_dist_pos = openfpm::vector<Point<dim, St>,Memory,layout_base>,
		 typename vector_dist_prop = openfpm::vector<prop,Memory,layout_base> >
class vector_dist : public vector_dist_comm<dim,St,prop,Decomposition,Memory,layout_base>,
#ifdef CUDA_GPU
					private vector_dist_ker_list<vector_dist_ker<dim,St,prop,layout_base>>
#else
					private vector_dist_ker_list<int>
#endif
{

public:

	//! Self type
	typedef vector_dist<dim,St,prop,Decomposition,Memory,layout_base> self;

	//! template parameters typedefs
	static const unsigned int dims = dim;
	typedef St stype;
	typedef prop value_type;
	typedef Decomposition Decomposition_type;
	typedef Memory Memory_type;

private:

	//! Ghost marker, all the particle with id > ghostMarker are ghost all with ghostMarker < are real particle
	size_t ghostMarker = 0;

	//! Particle position vector, (It has 2 elements) the first has real particles assigned to a processor
	//! the second element contain unassigned particles
	vector_dist_pos vPos;

	//! Used when cell list is constructed with CL_GPU_REORDER_POSITION. Stores a copy of vPos
	//! Stores a sorted copy of vPos in between updateCellListGPU and restoreOrder
	vector_dist_pos vPosReordered;

	//! Particle properties vector, (It has 2 elements) the first has real particles assigned to a processor
	//! the second element contain unassigned particles
	vector_dist_prop vPrp;

	//! Used when cell list is constructed with CL_GPU_REORDER_PROPERTY. Stores a copy of vPrp
	//! Stores a sorted copy of vPos in between updateCellListGPU and restoreOrder
	vector_dist_prop vPrpReordered;

	//! option used to create this vector
	size_t opt = 0;

	//! Name of the properties
	openfpm::vector<std::string> prp_names;

#ifdef SE_CLASS3

	se_class3_vector<prop::max_prop,dim,St,Decomposition,self> se3;

#endif

#ifdef SE_CLASS1
	 int map_ctr=0;
	 //ghost check to be done.
     int ghostget_ctr=0;
#endif

	/*! \brief Initialize the structures
	 *
	 * \param np number of particles
	 *
	 */
	void init_structures(size_t np)
	{
		Vcluster<Memory> & v_cl = create_vcluster<Memory>();

		// convert to a local number of elements
		size_t p_np = np / v_cl.getProcessingUnits();

		// Get non divisible part
		size_t r = np % v_cl.getProcessingUnits();

		// Distribute the remain particles
		if (v_cl.getProcessUnitID() < r)
			p_np++;

		// resize the position vector
		vPos.resize(p_np);

		// resize the properties vector
		vPrp.resize(p_np);

		ghostMarker = p_np;
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
	 * \param cellList cell-list
	 *
	 */
	template<typename CellL, typename sfc_it>
	void reorder_sfc(openfpm::vector<Point<dim,St>> & v_pos_dest,
						 openfpm::vector<prop> & v_prp_dest,
						 sfc_it & h_it,
						 CellL & cellList)
	{
		v_pos_dest.resize(vPos.size());
		v_prp_dest.resize(vPrp.size());

		//Index for v_pos_dest
		size_t count = 0;

		grid_key_dx<dim> ksum;

		for (size_t i = 0; i < dim ; i++)
		{ksum.set_d(i,cellList.getPadding(i));}

		while (h_it.isNext())
		{
		  auto key = h_it.get();
		  key += ksum;

		  size_t lin = cellList.getGrid().LinId(key);

		  // for each particle in the Cell "lin"
		  for (size_t i = 0; i < cellList.getNelements(lin); i++)
		  {
			  //reorder
			  auto v = cellList.get(lin,i);
			  v_pos_dest.get(count) = vPos.get(v);
			  v_prp_dest.get(count) = vPrp.get(v);

			  count++;
		  }
		  ++h_it;
		}
	}

public:
	typedef decltype(vPos) internal_position_vector_type;

	typedef CellList<dim, St, Mem_fast<>, shift<dim, St>, internal_position_vector_type > CellList_type;

	//! yes I am vector dist
	typedef int yes_i_am_vector_dist;

	//! yes I am vector subset dist
	typedef std::integral_constant<bool,false> is_it_a_subset;


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

		ghostMarker = v.ghostMarker;
		vPos = v.vPos;
		vPrp = v.vPrp;

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

		ghostMarker = v.ghostMarker;
		vPos.swap(v.vPos);
		vPrp.swap(v.vPrp);

#ifdef SE_CLASS3
		se3 = v.se3;
#endif

		opt = v.opt;

		return *this;
	}

	// default constructor (structure contain garbage)
	vector_dist()
	{}


	/*! \brief Copy Constructor
	 *
	 * \param v vector to copy
	 *
	 */
	vector_dist(const vector_dist<dim,St,prop,Decomposition,Memory,layout_base> & v)
	:vector_dist_comm<dim,St,prop,Decomposition,Memory,layout_base>(v.getDecomposition()) SE_CLASS3_VDIST_CONSTRUCTOR
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
	vector_dist_comm<dim,St,prop,Decomposition,Memory,layout_base>(dec) SE_CLASS3_VDIST_CONSTRUCTOR
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
	vector_dist(size_t np, Box<dim, St> box, const size_t (&bc)[dim], const Ghost<dim, St> & g, const grid_sm<dim,void> & gdist)
	:opt(0) SE_CLASS3_VDIST_CONSTRUCTOR
	{
		if (opt >> 32 != 0)
		{this->setDecompositionGranularity(opt >> 32);}

		check_parameters(box);

		init_structures(np);

		this->init_decomposition_gr_cell(box,bc,g,opt,gdist);


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
	:opt(opt) SE_CLASS3_VDIST_CONSTRUCTOR
	{
#ifdef SE_CLASS2
		check_new(this,8,VECTOR_DIST_EVENT,4);
#endif

		if (opt >> 32 != 0)
		{this->setDecompositionGranularity(opt >> 32);}

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

	/*! Set reference counter to one
	 *
	 */
	void setReferenceCounterToOne()
	{
		for (int i = 0 ; i < vPos.template getMemory<0>().ref() - 1; i++)
		{
			vPos.template getMemory<0>().decRef();
		}

		for (int i = 0 ; i < vPrp.template getMemory<0>().ref() - 1; i++)
		{
			decrement_memory<decltype(vPrp)> m(vPrp);

			boost::mpl::for_each_ref<boost::mpl::range_c<int,0,decltype(vPrp)::value_type::max_prop>>(m);
		}
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
		return ghostMarker;
	}

	/*! \brief return the local size of the vector
	 *
	 * \return local size
	 *
	 */
	size_t size_local_with_ghost() const
	{
		return vPos.size();
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
	inline auto getPos(vect_dist_key_dx vec_key) -> decltype(vPos.template get<0>(vec_key.getKey()))
	{
#ifdef SE_CLASS3
		check_for_pos_nan_inf<prop::max_prop_real,prop::max_prop>(*this,vec_key.getKey());
#endif

		return vPos.template get<0>(vec_key.getKey());
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
	inline auto getPos(vect_dist_key_dx vec_key) const -> decltype(vPos.template get<0>(vec_key.getKey()))
	{
#ifdef SE_CLASS3
		check_for_pos_nan_inf<prop::max_prop_real,prop::max_prop>(*this,vec_key.getKey());
#endif
		return vPos.template get<0>(vec_key.getKey());
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
	inline auto getPos(size_t vec_key) -> decltype(vPos.template get<0>(vec_key))
	{
#ifdef SE_CLASS3
		check_for_pos_nan_inf<prop::max_prop_real,prop::max_prop>(*this,vec_key);
#endif
		return vPos.template get<0>(vec_key);
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
	inline auto getPos(size_t vec_key) const -> decltype(vPos.template get<0>(vec_key))
	{
#ifdef SE_CLASS3
		check_for_pos_nan_inf<prop::max_prop_real,prop::max_prop>(*this,vec_key);
#endif
		return vPos.template get<0>(vec_key);
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
	template<unsigned int id> inline auto getProp(vect_dist_key_dx vec_key) -> decltype(vPrp.template get<id>(vec_key.getKey()))
	{
#ifdef SE_CLASS3
		check_for_prop_nan_inf<id,prop::max_prop+SE3_STATUS>(*this,vec_key.getKey());
#endif
		return vPrp.template get<id>(vec_key.getKey());
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
	template<unsigned int id> inline auto getProp(vect_dist_key_dx vec_key) const -> decltype(vPrp.template get<id>(vec_key.getKey()))
	{
#ifdef SE_CLASS3
		check_for_prop_nan_inf<id,prop::max_prop+SE3_STATUS>(*this,vec_key.getKey());
#endif
		return vPrp.template get<id>(vec_key.getKey());
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
	template<unsigned int id> inline auto getProp(size_t vec_key) -> decltype(vPrp.template get<id>(vec_key))
	{
#ifdef SE_CLASS3
		check_for_prop_nan_inf<id,prop::max_prop+SE3_STATUS>(*this,vec_key);
#endif
		return vPrp.template get<id>(vec_key);
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
	template<unsigned int id> inline auto getProp(size_t vec_key) const -> decltype(vPrp.template get<id>(vec_key))
	{
#ifdef SE_CLASS3
		check_for_prop_nan_inf<id,prop::max_prop+SE3_STATUS>(*this,vec_key);
#endif
		return vPrp.template get<id>(vec_key);
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
	inline auto getPosNC(vect_dist_key_dx vec_key) -> decltype(vPos.template get<0>(vec_key.getKey()))
	{
		return vPos.template get<0>(vec_key.getKey());
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
	inline auto getPosNC(vect_dist_key_dx vec_key) const -> decltype(vPos.template get<0>(vec_key.getKey()))
	{
		return vPos.template get<0>(vec_key.getKey());
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
	inline auto getPosNC(size_t vec_key) -> decltype(vPos.template get<0>(vec_key))
	{
		return vPos.template get<0>(vec_key);
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
	inline auto getPosNC(size_t vec_key) const -> decltype(vPos.template get<0>(vec_key))
	{
		return vPos.template get<0>(vec_key);
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
	template<unsigned int id> inline auto getPropNC(vect_dist_key_dx vec_key) -> decltype(vPrp.template get<id>(vec_key.getKey()))
	{
		return vPrp.template get<id>(vec_key.getKey());
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
	template<unsigned int id> inline auto getPropNC(vect_dist_key_dx vec_key) const -> decltype(vPrp.template get<id>(vec_key.getKey()))
	{
		return vPrp.template get<id>(vec_key.getKey());
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
	template<unsigned int id> inline auto getPropNC(size_t vec_key) -> decltype(vPrp.template get<id>(vec_key))
	{
		return vPrp.template get<id>(vec_key);
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
	template<unsigned int id> inline auto getPropNC(size_t vec_key) const -> decltype(vPrp.template get<id>(vec_key))
	{
		return vPrp.template get<id>(vec_key);
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
	inline auto getPosWrite(vect_dist_key_dx vec_key) -> decltype(vPos.template get<0>(vec_key.getKey()))
	{
#ifdef SE_CLASS3
		se3.template write<prop::max_prop_real>(*this,vec_key.getKey());
#endif

		return vPos.template get<0>(vec_key.getKey());
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
	inline auto getPosRead(vect_dist_key_dx vec_key) const -> decltype(vPos.template get<0>(vec_key.getKey()))
	{
#ifdef SE_CLASS3
		se3.template read<prop::max_prop_real>(*this,vec_key.getKey());
#endif

		return vPos.template get<0>(vec_key.getKey());
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
	template<unsigned int id> inline auto getPropWrite(vect_dist_key_dx vec_key) -> decltype(vPrp.template get<id>(vec_key.getKey()))
	{
#ifdef SE_CLASS3
		se3.template write<id>(*this,vec_key.getKey());
#endif

		return vPrp.template get<id>(vec_key.getKey());
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
	template<unsigned int id> inline auto getPropRead(vect_dist_key_dx vec_key) const -> decltype(vPrp.template get<id>(vec_key.getKey()))
	{
#ifdef SE_CLASS3
		se3.template read<id>(*this,vec_key.getKey());
#endif

		return vPrp.template get<id>(vec_key.getKey());
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
		vPrp.insert(ghostMarker);
		vPos.insert(ghostMarker);

		ghostMarker++;

#ifdef SE_CLASS3
		for (size_t i = 0 ; i < prop::max_prop_real+1 ; i++)
			vPrp.template get<prop::max_prop_real>(ghostMarker-1)[i] = UNINITIALIZED;
#endif
	}

    //////////////////////////////////////////////

	/*! \brief Add at the END of local and ghost particle
	 *
	 * It add a local particle at the end of local and ghost, with "local" we mean in this processor
	 * the particle can be also created out of the processor domain, in this
	 * case a call to map is required. Added particles are always created at the
	 * end and can be accessed with getLastPos and getLastProp
	 *
	 */
	void appendLocal()
	{
		vPos.add();
	}

#ifndef ONLY_READWRITE_GETTER

	/*! \brief Get the position of the last element
	 *
	 * \return the position of the element in space
	 *
	 */
	inline auto getLastPos() -> decltype(vPos.template get<0>(0))
	{
		return vPos.template get<0>(ghostMarker - 1);
	}


    /*! \brief Get the position of the last element after ghost
	 *
	 * \return the position of the element in space
	 *
	 */
	inline auto getLastPosEnd() -> decltype(vPos.template get<0>(0))
	{
		return vPos.template get<0>(vPos.size() - 1);
	}

	/*! \brief Get the property of the last element
	 *
	 * \tparam id property id
	 *
	 * \return return the selected property of the vector element
	 *
	 */
	template<unsigned int id> inline auto getLastProp() -> decltype(vPrp.template get<id>(0))
	{
		return vPrp.template get<id>(ghostMarker - 1);
	}

#endif

////////////////////////////// READ AND WRITE VERSION //////////

	/*! \brief Get the position of the last element
	 *
	 * \return the position of the element in space
	 *
	 */
	inline auto getLastPosRead() -> decltype(vPos.template get<0>(0))
	{
#ifdef SE_CLASS3
		se3.template read<prop::max_prop_real>(*this,ghostMarker-1);
#endif

		return vPos.template get<0>(ghostMarker - 1);
	}

	/*! \brief Get the property of the last element
	 *
	 * \tparam id property id
	 *
	 * \return return the selected property of the vector element
	 *
	 */
	template<unsigned int id> inline auto getLastPropRead() -> decltype(vPrp.template get<id>(0))
	{
#ifdef SE_CLASS3
		se3.read<id>(*this,ghostMarker-1);
#endif

		return vPrp.template get<id>(ghostMarker - 1);
	}


	/*! \brief Get the position of the last element
	 *
	 * \return the position of the element in space
	 *
	 */
	inline auto getLastPosWrite() -> decltype(vPos.template get<0>(0))
	{
#ifdef SE_CLASS3
		se3.template write<prop::max_prop_real>(*this,ghostMarker-1);
#endif

		return vPos.template get<0>(ghostMarker - 1);
	}

	/*! \brief Get the property of the last element
	 *
	 * \tparam id property id
	 *
	 * \return return the selected property of the vector element
	 *
	 */
	template<unsigned int id> inline auto getLastPropWrite() -> decltype(vPrp.template get<id>(0))
	{
#ifdef SE_CLASS3
		se3.template write<id>(*this,ghostMarker-1);
#endif

		return vPrp.template get<id>(ghostMarker - 1);
	}

////////////////////////////////////////////////////////////////

	/*! \brief Construct a cell list symmetric based on a cut-off radius
	 *
	 * \tparam CellList_type CellList type to construct
	 *
	 * \param r_cut interation radius, or size of each cell
	 *
	 * \return the Cell list
	 *
	 */
	template<typename CellList_type = CellList<dim, St, Mem_fast<>, shift<dim, St>,internal_position_vector_type > >
	CellList_type getCellListSym(St r_cut, size_t opt = CL_LINEAR_CELL_KEYS)
	{
		return getCellList(r_cut, opt | CL_SYMMETRIC);
	}

	/*! \brief Construct a local symmetric cell list based on a cut-off radius
	 *
	 * \tparam CellList_type CellList type to construct
	 *
	 * \param r_cut interation radius, or size of each cell
	 *
	 * \return the Cell list
	 *
	 */
	template<typename CellList_type = CellList<dim, St, Mem_fast<>, shift<dim, St>,internal_position_vector_type > >
	CellList_type getCellListSymLocal(St r_cut, size_t opt = CL_LINEAR_CELL_KEYS)
	{
		return getCellList(r_cut, opt | CL_LOCAL_SYMMETRIC);
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
	template<typename CellList_type = CellList<dim, St, Mem_fast<>, shift<dim, St>, decltype(vPos) > >
	CellList_type getCellList(St r_cut, size_t opt = CL_NON_SYMMETRIC | CL_LINEAR_CELL_KEYS, bool no_se3 = false, float ghostEnlargeFactor = 1.013)
	{
		if (opt & CL_SYMMETRIC || opt & CL_LOCAL_SYMMETRIC) {
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
			size_t pad = 0;
			CellDecomposer_sm<dim,St,shift<dim,St>> cd_sm;
			cl_param_calculateSym(getDecomposition().getDomain(),cd_sm,getDecomposition().getGhost(),r_cut,pad);

			// Processor bounding box
			Box<dim, St> pbox = getDecomposition().getProcessorBounds();

			CellList_type cellList;

			cellList.setOpt(opt);
			cellList.Initialize(cd_sm,pbox,pad);
			cellList.set_ndec(getDecomposition().get_ndec());
			cellList.setGhostMarker(ghostMarker);

			cellList.fill(vPos, vPrp, ghostMarker);

			return cellList;
		}

		// CL_NON_SYMMETRIC
		else {
#ifdef SE_CLASS3
			if (no_se3 == false)
			{se3.getNN();}
#endif
#ifdef SE_CLASS1
			check_ghost_compatible_rcut(r_cut);
#endif
			// Get ghost and anlarge by 1%
			Ghost<dim,St> ghostEnlarge = getDecomposition().getGhost();
			ghostEnlarge.magnify(ghostEnlargeFactor);

			// Division array
			size_t div[dim];

			// get the processor bounding box
			Box<dim, St> pbox = getDecomposition().getProcessorBounds();

			// Processor bounding box
			cl_param_calculate(pbox, div, r_cut, ghostEnlarge);

			CellList_type cellList;

			cellList.setOpt(opt);
			cellList.Initialize(pbox, div);
			cellList.setGhostMarker(ghostMarker);
			cellList.set_ndec(getDecomposition().get_ndec());

			cellList.fill(vPos, vPrp, ghostMarker);

			return cellList;
		}
	}


	/*! \brief Construct a symmetric cell list based on a number of divisions and padding
	 *
	 * \tparam CellL CellList type to construct
	 *
	 * \param div number of cells per dimension
	 * \param number of padding layers per dimension
	 *
	 * \return the Cell list
	 *
	 */
	template<typename CellL = CellList<dim, St, Mem_fast<>, shift<dim, St> > >
	CellL getCellListSym(
		const size_t (& div)[dim],
		const size_t (& pad)[dim],
		size_t opt = CL_LINEAR_CELL_KEYS)
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
			if (pad[i] > pad_max) pad_max = pad[i];


		CellDecomposer_sm<dim,St,shift<dim,St>> cd_sm;
		cd_sm.setDimensions(getDecomposition().getDomain(),div,pad_max);

		// Processor bounding box
		Box<dim, St> pbox = getDecomposition().getProcessorBounds();

		CellL cellList;

		cellList.setOpt(opt | CL_SYMMETRIC);
		cellList.Initialize(cd_sm,pbox,pad_max);
		cellList.set_ndec(getDecomposition().get_ndec());

		cellList.fill(vPos, vPrp, ghostMarker);

		return cellList;
	}


#ifdef CUDA_GPU

	/*! \brief Construct a cell list starting from the stored particles
	 *
	 * \param r_cut interation radius, or size of each cell
	 * \param NNIteratorBox sets the number of neighborhood cell layers used to iterate through by getNNIteratorBox()
	 *
	 * \return the Cell list
	 *
	 */
	template<typename CellType = CellList_gpu<dim,St,CudaMemory,shift_only<dim, St>>>
	CellType getCellListGPU(St r_cut, size_t opt = CL_NON_SYMMETRIC, size_t NNIteratorBox = 1, bool no_se3 = false, float ghostEnlargeFactor = 1.013)
	{
#ifdef SE_CLASS3
		if (no_se3 == false)
		{se3.getNN();}
#endif
#ifdef SE_CLASS1
		check_ghost_compatible_rcut(r_cut);
#endif

		// Get ghost and anlarge by 1%
		Ghost<dim,St> ghostEnlarge = getDecomposition().getGhost();
		ghostEnlarge.magnify(ghostEnlargeFactor);

		Vcluster<Memory> & v_cl = create_vcluster<Memory>();

		// Division array
		size_t div[dim];

		// get the processor bounding box
		Box<dim, St> pbox = getDecomposition().getProcessorBounds();

		// Processor bounding box
		cl_param_calculate(pbox, div, r_cut, ghostEnlarge);

		CellType cellList(pbox,div);

		// getNNIteratorBox has to be set here compared to getNNIteratorRadius (could be set later)
		// As sparse cell list on gpu uses it in construct()
		// If not set, the dafault value of 1 would have been used by construct()
		cellList.setOpt(opt);
		cellList.setBoxNN(NNIteratorBox);
		cellList.set_ndec(getDecomposition().get_ndec());
		cellList.setGhostMarker(ghostMarker);

		return cellList;

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
	template<typename CellL = CellList<dim, St, Mem_fast<>, shift<dim, St> > >
	CellL getCellList_hilb(St r_cut, size_t opt = CL_NON_SYMMETRIC)
	{
		return getCellList(r_cut, opt | CL_HILBERT_CELL_KEYS);
	}

	/*! \brief Update a cell list using the stored particles
	 *
	 * \tparam CellL CellList type to construct
	 *
	 * \param cellList Cell list to update
	 * \param no_se3 avoid se class 3 checking
	 *
	 */
	template<unsigned int ... prp,typename CellL>
	void updateCellList(CellL & cellList, bool no_se3 = false)
	{
#ifdef SE_CLASS3
		if (no_se3 == false)
		{se3.getNN();}
#endif

		// Here we have to check that the Cell-list has been constructed
		// from the same decomposition
		bool to_reconstruct = cellList.get_ndec() != getDecomposition().get_ndec();

		if (to_reconstruct == false)
		{
			cellList.fill(vPos, vPrp, ghostMarker);
			cellList.setGhostMarker(ghostMarker);
		}
		else
		{
			// This function assume equal spacing in all directions
			// but in the worst case we take the maximum
			St r_cut = cellList.getCellBox().getRcut();

			if (cellList.getOpt() & CL_SYMMETRIC) {
				CellL cellListTmp = getCellListSym<CellL>(cellList.getDivWP(), cellList.getPadding());
				cellList.swap(cellListTmp);
			}

			else {
				CellL cellListTmp = getCellList<CellL>(r_cut, cellList.getOpt());
				cellList.swap(cellListTmp);
			}

		}
	}


#ifdef CUDA_GPU

	/*! \brief Update a cell list on GPU using the stored particles
	 *
	 * \tparam CellL CellList type to construct
	 *
	 * \param cellList Cell list to update
	 * \param no_se3 avoid se class 3 checking
	 *
	 */
	template<unsigned int ... prp, typename CellList_type>
	void updateCellListGPU(CellList_type & cellList, bool no_se3 = false)
	{
#ifdef SE_CLASS3
		if (no_se3 == false)
		{se3.getNN();}
#endif
		Vcluster<Memory> & v_cl = create_vcluster<Memory>();

		if (cellList.getOpt() & CL_GPU_REORDER_POSITION)
			if (vPosReordered.size() != vPos.size()) vPosReordered.resize(vPos.size());

		if (cellList.getOpt() & CL_GPU_REORDER_PROPERTY)
			if (vPrpReordered.size() != vPrp.size()) vPrpReordered.resize(vPrp.size());

		// Here we have to check that the Cell-list has been constructed
		// from the same decomposition
		bool to_reconstruct = cellList.get_ndec() != getDecomposition().get_ndec();

		if (to_reconstruct == false)
		{
			if (cellList.getOpt() & CL_GPU_REORDER_POSITION || cellList.getOpt() & CL_GPU_REORDER_PROPERTY) {

				cellList.template construct<decltype(vPos),decltype(vPrp),prp ...>(
					vPos,
					vPrp,
					vPosReordered,
					vPrpReordered,
					v_cl.getGpuContext(),
					ghostMarker,
					0,
					vPos.size());

				if (cellList.getOpt() & CL_GPU_REORDER_POSITION) vPos.swap(vPosReordered);
				if (cellList.getOpt() & CL_GPU_REORDER_PROPERTY) vPrp.swap(vPrpReordered);
			}
			else
				cellList.construct(vPos,vPrp,v_cl.getGpuContext(),ghostMarker,0,vPos.size());
		}

		else
		{
			// This function assume equal spacing in all directions
			// but in the worst case we take the maximum
			St r_cut = cellList.getCellBox().getRcut();

			CellList_type cellListTmp = getCellListGPU<CellList_type>(r_cut, cellList.getOpt(), cellList.getBoxNN());

			if (cellList.getOpt() & CL_GPU_REORDER_POSITION || cellList.getOpt() & CL_GPU_REORDER_PROPERTY) {
				cellListTmp.template construct<decltype(vPos),decltype(vPrp),prp ...>(
					vPos,
					vPrp,
					vPosReordered,
					vPrpReordered,
					v_cl.getGpuContext(),
					ghostMarker,
					0,
					vPos.size());

				if (cellList.getOpt() & CL_GPU_REORDER_POSITION) vPos.swap(vPosReordered);
				if (cellList.getOpt() & CL_GPU_REORDER_PROPERTY) vPrp.swap(vPrpReordered);
			}
			else
				cellListTmp.construct(vPos,vPrp,v_cl.getGpuContext(),ghostMarker,0,vPos.size());

			cellList.swap(cellListTmp);
		}
	}

#endif

	/*! \brief for each particle get the symmetric verlet list
	 *
	 * \param r_cut cut-off radius
	 *
	 * \return the verlet list
	 *
	 */
	template <typename VerletList_type = VerletList<dim,St,VL_SYMMETRIC,Mem_fast<>,shift<dim,St>,decltype(vPos)>>
	VerletList_type getVerletSym(St r_cut)
	{
		return getVerlet<VL_SYMMETRIC, VerletList_type>(r_cut);
	}

	/*! \brief for each particle get the symmetric verlet list
	 *
	 * \param r_cut cut-off radius
	 *
	 * \return the verlet list
	 *
	 */
	template <typename VerletList_type = VerletList<dim,St,VL_CRS_SYMMETRIC,Mem_fast<>,shift<dim,St>,decltype(vPos)>>
	VerletList_type getVerletCrs(St r_cut)
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

		VerletList_type verletList;

		// Processor bounding box
		Box<dim, St> pbox = getDecomposition().getProcessorBounds();

		// Initialize the verlet list
		verletList.initializeCrs(getDecomposition().getDomain(),pbox,getDecomposition().getGhost(),r_cut,vPos,ghostMarker);

		// Get the internal cell list
		auto & cellList = verletList.getInternalCellList();

		// Shift
		grid_key_dx<dim> shift;

		// Add padding
		for (size_t i = 0 ; i < dim ; i++)
			shift.set_d(i,cellList.getPadding(i));

		grid_sm<dim,void> gs = cellList.getInternalGrid();

		getDecomposition().setNNParameters(shift,gs);

		verletList.fillCRSSymmetric(vPos,r_cut,ghostMarker,
			getDecomposition().getCRSDomainCells(),
			getDecomposition().getCRSAnomDomainCells()
		);

		verletList.set_ndec(getDecomposition().get_ndec());

		return verletList;
	}

	/*! \brief Get Verlet list with unique cut-off radius for every particle
	 *
	 * \param rCuts list of cut-off radii for every particle in vPos
	 *
	 * \return the verlet list
	 *
	 */
	template <typename VerletList_type = VerletList<dim,St,VL_NON_SYMMETRIC|VL_ADAPTIVE_RCUT,Mem_fast<>,shift<dim,St>,decltype(vPos)>>
	VerletList_type getVerletAdaptRCut(openfpm::vector<St> &rCuts)
	{
#ifdef SE_CLASS3
		se3.getNN();
#endif

		VerletList_type verletList;

		// get the processor bounding box
		Ghost<dim,St> g = getDecomposition().getGhost();
		g.magnify(1.013);

		Box<dim, St> pbox = getDecomposition().getProcessorBounds();
		// enlarge the box where the Verlet is defined
		pbox.enlarge(g);

		verletList.InitializeNonSymmAdaptive(pbox,rCuts,vPos,ghostMarker);

		return verletList;
	}

	/*! \brief for each particle get the verlet list
	 *
	 * \param r_cut cut-off radius
	 * \param nMax Max number of neighboring particles added per p (sorted by distance). By default all are included
	 *	Should be supplied additionally with option flag VL_NMAX_NEIGHBOR to VerletList_type template type
	 *
	 * \return a VerletList object
	 *
	 */
	template <unsigned int optVerlet=VL_NON_SYMMETRIC, typename VerletList_type = VerletList<dim,St,optVerlet,Mem_fast<>,shift<dim,St>,decltype(vPos)>>
	VerletList_type getVerlet(St r_cut, size_t neighborMaxNum = 0)
	{
#ifdef SE_CLASS3
		se3.getNN();
#endif

		VerletList_type verletList;

		if (neighborMaxNum)
			verletList.setNeighborMaxNum(neighborMaxNum);

		// get the processor bounding box
		Box<dim, St> pbox = getDecomposition().getProcessorBounds();

		if (verletList.getOpt() & VL_SYMMETRIC)
		{
			verletList.InitializeSym(getDecomposition().getDomain(),pbox,getDecomposition().getGhost(),r_cut,vPos,ghostMarker);
		}

		else
		{
			Ghost<dim,St> g = getDecomposition().getGhost();
			g.magnify(1.013);
			// enlarge the box where the Verlet is defined
			pbox.enlarge(g);

			verletList.Initialize(pbox,r_cut,vPos,ghostMarker);
		}

		verletList.set_ndec(getDecomposition().get_ndec());

		return verletList;
	}

	/*! \brief for each particle get the verlet list
	 *
	 * \param r_cut cut-off radius
	 * \param verletList Verlet to update
	 *
	 */
	template<unsigned int opt, typename Mem_type>
	void updateVerlet(VerletList<dim,St,opt,Mem_type,shift<dim,St> > & verletList, St r_cut)
	{
#ifdef SE_CLASS3
		se3.getNN();
#endif
		if ((opt & VL_SYMMETRIC) || (opt & VL_NON_SYMMETRIC))
		{
			auto & cellList = verletList.getInternalCellList();

			// Here we have to check that the Box defined by the Cell-list is the same as the domain box of this
			// processor. if it is not like that we have to completely reconstruct from stratch
			bool to_reconstruct = cellList.get_ndec() != getDecomposition().get_ndec();

			if (to_reconstruct == false)
			{
				verletList.update(getDecomposition().getDomain(),r_cut,vPos,ghostMarker);
			}

			else
			{
				VerletList<dim,St,opt,Mem_type,shift<dim,St>> ver_tmp = getVerlet<opt, VerletList<dim,St,opt,Mem_type,shift<dim,St>>>(r_cut);
				verletList.swap(ver_tmp);
			}
		}
		else if (opt & VL_CRS_SYMMETRIC)
		{
#ifdef SE_CLASS1
			if ((this->opt & BIND_DEC_TO_GHOST))
			{
				std::cerr << __FILE__ << ":" << __LINE__ << " Error the vector has been constructed without BIND_DEC_TO_GHOST, updateVerlet with the option VL_CRS_SYMMETRIC require the vector to be constructed with BIND_DEC_TO_GHOST option " << std::endl;
				ACTION_ON_ERROR(VECTOR_DIST_ERROR_OBJECT);
			}
#endif

			auto & cellList = verletList.getInternalCellList();

			// Here we have to check that the Box defined by the Cell-list is the same as the domain box of this
			// processor. if it is not like that we have to completely reconstruct from stratch
			bool to_reconstruct = cellList.get_ndec() != getDecomposition().get_ndec();

			if (to_reconstruct == false)
			{
				// Shift
				grid_key_dx<dim> shift;

				// Add padding
				for (size_t i = 0 ; i < dim ; i++)
					shift.set_d(i,cellList.getPadding(i));

				grid_sm<dim,void> gs = cellList.getInternalGrid();

				getDecomposition().setNNParameters(shift,gs);

				verletList.updateCrs(vPos,r_cut,ghostMarker,
					getDecomposition().getCRSDomainCells(),
					getDecomposition().getCRSAnomDomainCells()
				);
			}
			else
			{
				VerletList<dim,St,opt,Mem_type,shift<dim,St>> ver_tmp = getVerletCrs<VerletList<dim,St,opt,Mem_type,shift<dim,St> >>(r_cut);
				verletList.swap(ver_tmp);
			}
		}
	}

	/*! \brief Update non-symmetric adaptive r-cut Verlet list
	 *
	 * \param verletList Verlet to update
	 * \param rCuts list of cut-off radii for every particle in pos
	 *
	 */
	template<unsigned int opt, typename Mem_type>
	void updateVerletAdaptRCut(VerletList<dim,St,opt,Mem_type,shift<dim,St> > & verletList, openfpm::vector<St> &rCuts)
	{
#ifdef SE_CLASS3
		se3.getNN();
#endif
		// in this mode the Verlet list doesn't depend on the decomposition counter
		// has to be fully reconstructed on update
		Box<dim, St> pbox = getDecomposition().getProcessorBounds();
		verletList.InitializeNonSymmAdaptive(pbox,rCuts,vPos,ghostMarker);
	}


	/*! \brief Construct a cell list starting from the stored particles and reorder a vector according to the Hilberts curve
	 *
	 * \tparam CellL CellList type to construct
	 *
	 * \param m an order of a hilbert curve
	 *
	 */
	template<typename CellL=CellList<dim,St,Mem_bal<>,shift<dim,St> > >
	void reorder (int32_t m, reorder_opt opt = reorder_opt::HILBERT)
	{
		reorder<CellL>(m,getDecomposition().getGhost(),opt);
	}


	/*! \brief Construct a cell list starting from the stored particles and reorder a vector according to the Hilberts curve
	 *
	 * \warning it kill the ghost and invalidate cell-lists
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
	template<typename CellL=CellList<dim,St,Mem_bal<>,shift<dim,St> > >
	void reorder(int32_t m, const Ghost<dim,St> & enlarge, reorder_opt opt = reorder_opt::HILBERT)
	{
		// reset the ghost part
		vPos.resize(ghostMarker);
		vPrp.resize(ghostMarker);


		CellL cellList;

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

		cellList.Initialize(pbox,div);
		cellList.setGhostMarker(ghostMarker);

		// for each particle add the particle to the cell list

		auto it = getIterator();

		while (it.isNext())
		{
			auto key = it.get();

			Point<dim,St> xp = this->getPos(key);

			cellList.add(xp,key.getKey());

			++it;
		}

		// Use cellList to reorder vPos

		//destination vector
		openfpm::vector<Point<dim,St>> v_pos_dest;
		openfpm::vector<prop> v_prp_dest;

		if (opt == reorder_opt::HILBERT)
		{
			grid_key_dx_iterator_hilbert<dim> h_it(m);

			reorder_sfc<CellL,grid_key_dx_iterator_hilbert<dim>>(v_pos_dest,v_prp_dest,h_it,cellList);
		}
		else if (opt == reorder_opt::LINEAR)
		{
			grid_sm<dim,void> gs(div);
			grid_key_dx_iterator<dim> h_it(gs);

			reorder_sfc<CellL,grid_key_dx_iterator<dim>>(v_pos_dest,v_prp_dest,h_it,cellList);
		}
		else
		{
			// We do nothing, we second swap nullify the first
			vPos.swap(v_pos_dest);
			vPrp.swap(v_prp_dest);
		}

		vPos.swap(v_pos_dest);
		vPrp.swap(v_prp_dest);
	}

	/*! \brief Construct a cell list starting from the stored particles and reorder a vector according to the Hilberts curve
	 *
	 * \warning it kill the ghost and invalidate cell-lists
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
	template<typename CellL=CellList<dim,St,Mem_bal<>,shift<dim,St> > >
	void reorder_rcut(St r_cut)
	{
		// reset the ghost part
		vPos.resize(ghostMarker);
		vPrp.resize(ghostMarker);

		auto cellList = getCellList<CellL>(r_cut);

		// Use cellList to reorder vPos

		//destination vector
		openfpm::vector<Point<dim,St>> vPosReorder;
		openfpm::vector<prop> vPrpReorder;

		size_t div[dim];
		for (size_t i = 0 ; i < dim ; i++)
		{div[i] = cellList.getGrid().size(i) - 2*cellList.getPadding()[i];}

		grid_sm<dim,void> gs(div);
		grid_key_dx_iterator<dim> h_it(gs);

		reorder_sfc<CellL,grid_key_dx_iterator<dim>>(vPosReorder,vPrpReorder,h_it,cellList);

		vPos.swap(vPosReorder);
		vPrp.swap(vPrpReorder);
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
		Vcluster<Memory> & v_cl = create_vcluster<Memory>();

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
		return vector_dist_iterator(0, vPos.size());
	}

	/*! \brief Get an iterator that traverse domain and ghost particles
	 *
	 * \param start particle
	 * \param stop particle
	 *
	 * \return an iterator
	 *
	 */
	vector_dist_iterator getIterator(size_t start, size_t stop)
	{
#ifdef SE_CLASS3
		se3.getIterator();
#endif
		return vector_dist_iterator(start, stop);
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

		return vector_dist_iterator(ghostMarker, vPos.size());
	}

	/*! \brief Get the iterator across the position of the ghost particles
	 *
	 * \return an iterator
	 *
	 */
	vector_dist_iterator getGhostIterator_no_se3() const
	{
		return vector_dist_iterator(ghostMarker, vPos.size());
	}

	/*! \brief Get an iterator that traverse the particles in the domain
	 *         using a cell list
	 *
	 * \param cellList Cell-list
	 *
	 * \return an iterator over the particles
	 *
	 */
	template<typename CellList> ParticleIt_Cells<dim,CellList>
	getDomainIteratorCells(CellList & cellList)
	{
#ifdef SE_CLASS3
		se3.getIterator();
#endif

		// Shift
		grid_key_dx<dim> shift;

		// Add padding
		for (size_t i = 0 ; i < dim ; i++)
			shift.set_d(i,cellList.getPadding(i));

		grid_sm<dim,void> gs = cellList.getInternalGrid();

		getDecomposition().setNNParameters(shift,gs);

		return ParticleIt_Cells<dim,CellList>(cellList,getDecomposition().getDomainCells(),ghostMarker);
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

		return vector_dist_iterator(0, ghostMarker);
	}

#ifdef CUDA_GPU

	/*! \brief Get an iterator that traverse the particles in the domain
	 *
	 * \return an iterator
	 *
	 */
	ite_gpu<1> getDomainIteratorGPU(size_t n_thr = default_kernel_wg_threads_) const
	{
#ifdef SE_CLASS3
		se3.getIterator();
#endif

		return vPos.getGPUIteratorTo(ghostMarker-1,n_thr);
	}

	/*! \brief Get an iterator that traverse the particles in the domain
	 *
	 * \return an iterator
	 *
	 */
	ite_gpu<1> getDomainAndGhostIteratorGPU(size_t n_thr = default_kernel_wg_threads_) const
	{
#ifdef SE_CLASS3
		se3.getIterator();
#endif

		return vPos.getGPUIteratorTo(vPos.size()-1,n_thr);
	}

	/*! \brief print a vector type property
	 *
	 *
	 * \tparam property
	 *
	 */
	template<unsigned int prp>
	void debugPrintVector()
	{
		this->vPrp.template deviceToHost<prp>();

		auto it = this->getDomainIterator();

		while(it.isNext())
		{
			auto p = it.get();

			for (size_t i = 0 ; i < std::extent<typename boost::mpl::at<typename prop::type,boost::mpl::int_<prp>>::type>::value ; i++)
			{
				std::cout << vPrp.template get<prp>(p.getKey())[i] << "   ";
			}

			std::cout << std::endl;

			++it;
		}
	}

	/*! \brief print a scalar type property
	 *
	 *
	 * \tparam property
	 *
	 */
	template<unsigned int prp>
	void debugPrintScalar()
	{
		this->vPrp.template deviceToHost<prp>();

		auto it = this->getDomainIterator();

		while(it.isNext())
		{
			auto p = it.get();

			std::cout << vPrp.template get<prp>(p.getKey()) << "   " << std::endl;

			++it;
		}
	}

#endif

#ifdef CUDA_GPU

        /*! \brief Get an iterator that traverse the particles in the domain
         *
         * \return an iterator
         *
         */
        auto getDomainIteratorDevice(size_t n_thr = default_kernel_wg_threads_) const -> decltype(this->getDomainIteratorGPU(n_thr))
        {
                return this->getDomainIteratorGPU(n_thr);
        }


#else

        /*! \brief Get an iterator that traverse the particles in the domain
         *
         * \return an iterator
         *
         */
        auto getDomainIteratorDevice(size_t n_thr = default_kernel_wg_threads_) const -> decltype(this->getDomainIterator())
        {
                return this->getDomainIterator();
        }


#endif

	/*! \brief Get an iterator that traverse the particles in the domain
	 *
	 * \return an iterator
	 *
	 */
	vector_dist_iterator getDomainIterator_no_se3() const
	{
		return vector_dist_iterator(0, ghostMarker);
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

		return vector_dist_iterator(0, vPos.size());
	}

	/*! \brief Get an iterator that traverse the particles in the domain
	 *
	 * \return an iterator
	 *
	 */
	vector_dist_iterator getDomainAndGhostIterator_no_se3() const
	{
		return vector_dist_iterator(0, vPos.size());
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

		this->template map_list_<prp...>(vPos,vPrp,ghostMarker,opt);

#ifdef CUDA_GPU
		this->update(this->toKernel());
#endif

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
#ifdef SE_CLASS1
	    map_ctr++;
#endif

		this->template map_<obp>(vPos,vPrp,ghostMarker,opt);

#ifdef CUDA_GPU
		this->update(this->toKernel());
#endif

#ifdef SE_CLASS3
		se3.map_post();
#endif
	}
#ifdef SE_CLASS1
    int getMapCtr() const
    {
	    return map_ctr;
    }
#endif


	/*! \brief Stub does not do anything
	*
	*/
	void ghost_get_subset()
	{
    /*  #ifdef SE_CLASS1
        This is not a ghost get on subset.
       std::cerr<<__FILE__<<":"<<__LINE__<<":You Used a ghost_get on a subset. This does not do anything. Please use ghostget on the entire set.";
    #endif */
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
		Vcluster<Memory> & v_cl = create_vcluster<Memory>();

		if (getDecomposition().getProcessorBounds().isValid() == false && size_local() != 0)
		{
			std::cerr << __FILE__ << ":" << __LINE__ << " Error the processor " << v_cl.getProcessUnitID() << " has particles, but is supposed to be unloaded" << std::endl;
			ACTION_ON_ERROR(VECTOR_DIST_ERROR_OBJECT);
		}
#endif

#ifdef SE_CLASS3
		se3.template ghost_get_pre<prp...>(opt);
#endif

		this->template ghost_get_<GHOST_SYNC,prp...>(vPos,vPrp,ghostMarker,opt);

#ifdef CUDA_GPU
		this->update(this->toKernel());
#endif

#ifdef SE_CLASS3

		this->template ghost_get_<prop::max_prop_real>(vPos,vPrp,ghostMarker,opt | KEEP_PROPERTIES);

		se3.template ghost_get_post<prp...>(opt);
#endif
	}


	/*! \brief It synchronize the properties and position of the ghost particles
	 *
	 * \tparam prp list of properties to get synchronize
	 *
	 * \param opt options WITH_POSITION, it send also the positional information of the particles
	 *
	 */
	template<int ... prp> inline void Ighost_get(size_t opt = WITH_POSITION)
	{
#ifdef SE_CLASS1
		Vcluster<Memory> & v_cl = create_vcluster<Memory>();

		if (getDecomposition().getProcessorBounds().isValid() == false && size_local() != 0)
		{
			std::cerr << __FILE__ << ":" << __LINE__ << " Error the processor " << v_cl.getProcessUnitID() << " has particles, but is supposed to be unloaded" << std::endl;
			ACTION_ON_ERROR(VECTOR_DIST_ERROR_OBJECT);
		}
#endif

#ifdef SE_CLASS3
		se3.template ghost_get_pre<prp...>(opt);
#endif

		this->template ghost_get_<GHOST_ASYNC,prp...>(vPos,vPrp,ghostMarker,opt);
	}

	/*! \brief It synchronize the properties and position of the ghost particles
	 *
	 * \tparam prp list of properties to get synchronize
	 *
	 * \param opt options WITH_POSITION, it send also the positional information of the particles
	 *
	 */
	template<int ... prp> inline void ghost_wait(size_t opt = WITH_POSITION)
	{
#ifdef SE_CLASS1
		Vcluster<Memory> & v_cl = create_vcluster<Memory>();

		if (getDecomposition().getProcessorBounds().isValid() == false && size_local() != 0)
		{
			std::cerr << __FILE__ << ":" << __LINE__ << " Error the processor " << v_cl.getProcessUnitID() << " has particles, but is supposed to be unloaded" << std::endl;
			ACTION_ON_ERROR(VECTOR_DIST_ERROR_OBJECT);
		}
#endif

		this->template ghost_wait_<prp...>(vPos,vPrp,ghostMarker,opt);

#ifdef CUDA_GPU
		this->update(this->toKernel());
#endif

#ifdef SE_CLASS3

		this->template ghost_get_<prop::max_prop_real>(vPos,vPrp,ghostMarker,opt | KEEP_PROPERTIES);

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
		this->template ghost_put_<op,prp...>(vPos,vPrp,ghostMarker,opt_);
	}

	/*! \brief Remove a set of elements from the distributed vector
	 *
	 * \param keys std::set of elements to eliminate.
	 * 	Values inside std::set are sorted by definition of the container
	 *
	 */
	void remove(std::set<size_t> & keys)
	{
		openfpm::vector<size_t> v_keys; v_keys.reserve(keys.size());

		for (auto it = keys.begin(); it != keys.end(); ++it)
			v_keys.add(*it);

		// keys are sorted and unique
		this->remove(v_keys, 0);
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
		vPos.remove(keys, start);
		vPrp.remove(keys, start);

		ghostMarker -= keys.size();
	}

	/*! \brief Remove a set of elements from the distributed vector
	 *
	 * \warning keys must be sorted
	 *
	 * \param keys vector of elements to eliminate
	 * \param start from where to eliminate
	 *
	 */
	void remove(openfpm::vector<aggregate<int>> & keys, size_t start = 0)
	{
		vPos.remove(keys, start);
		vPrp.remove(keys, start);

		ghostMarker -= keys.size();
	}

	/*! \brief Remove one element from the distributed vector
	 *
	 * \param key remove one element from the vector
	 *
	 */
	void remove(size_t key)
	{
		vPos.remove(key);
		vPrp.remove(key);

		ghostMarker--;
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
			Point<dim,St> p = vd.getPos(it.get());
			size_t v = cdsm.getCell(p);

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

		h5s.save(filename,vPos,vPrp);
	}

	/*! \brief Load the distributed vector from an HDF5 file
	 *
	 * \param filename file from where to load
	 *
	 */
	inline void load(const std::string & filename)
	{
		HDF5_reader<VECTOR_DIST> h5l;

		h5l.load(filename,vPos,vPrp,ghostMarker);
	}

	/*! \brief Reserve space for the internal vectors
	 *
	 * \param
	 *
	 */
	void setCapacity(unsigned int ns)
	{
		vPos.reserve(ns);
		vPrp.reserve(ns);
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
	inline bool write(std::string out ,int opt = VTK_WRITER)
	{
		return write(out,"",opt);
	}

	/*! \brief Output particle position and properties
	 *
	 * \param out output filename
	 * \param meta_info meta information example ("time = 1.234" add the information time to the VTK file)
	 * \param opt VTK_WRITER, CSV_WRITER, it is also possible to choose the format for  VTK
	 *            FORMAT_BINARY. (the default is ASCII format)
	 *
	 * \return true if the file has been written without error
	 *
	 */
	inline bool write(std::string out, std::string meta_info ,int opt = VTK_WRITER)
	{
		Vcluster<Memory> & v_cl = create_vcluster<Memory>();

		if ((opt & 0x0FFF0000) == CSV_WRITER)
		{
			// CSVWriter test
			CSVWriter<vector_dist_pos,
			          vector_dist_prop > csv_writer;

			std::string output = std::to_string(out + "_" + std::to_string(v_cl.getProcessUnitID()) + std::to_string(".csv"));

			// Write the CSV
			return csv_writer.write(output,vPos,vPrp);
		}
		else
		{
			file_type ft = file_type::ASCII;

			if (opt & FORMAT_BINARY)
				ft = file_type::BINARY;

			// VTKWriter for a set of points
			VTKWriter<boost::mpl::pair<vector_dist_pos,
									   vector_dist_prop>,
			                           VECTOR_POINTS> vtk_writer;
			vtk_writer.add(vPos,vPrp,ghostMarker);

			std::string output = std::to_string(out + "_" + std::to_string(v_cl.getProcessUnitID()) + std::to_string(".vtp"));
            //Create Directory for VTP files and write the PVTP metadata
            if(v_cl.rank()==0)
            {
                create_directory_if_not_exist("VTPDATA",1);
                vtk_writer.write_pvtp(out,prp_names,v_cl.size());
            }
            v_cl.barrier();
			// Write the VTK file
			bool ret=vtk_writer.write(output,prp_names,"particles",meta_info,ft);
			return ret;
		}
	}

	/*! \brief Delete the particles on the ghost
	 *
	 *
	 */
	void deleteGhost()
	{
		vPos.resize(ghostMarker);
		vPrp.resize(ghostMarker);
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

		vPos.resize(rs);
		vPrp.resize(rs);

		ghostMarker = rs;

#ifdef CUDA_GPU
		this->update(this->toKernel());
#endif
	}

    /*! \brief Resize the vector at the end of the ghost (locally)
	 *
	 * \warning It doesn't delete the ghosts
	 *
	 * \param rs
	 *
	 */
	void discardLocalAppend(size_t rs)
	{
		vPos.resize(rs);
#ifdef CUDA_GPU
		this->update(this->toKernel());
#endif
	}

	/*! \brief Output particle position and properties
	 *
	 * \param out output
	 * \param iteration (we can append the number at the end of the file_name)
	 * \param meta_info meta information example ("time = 1.234" add the information time to the VTK file)
	 * \param opt VTK_WRITER, CSV_WRITER, it is also possible to choose the format for  VTK
	 *            FORMAT_BINARY. (the default is ASCII format)
	 *
	 * \return if the file has been written correctly
	 *
	 */
	inline bool write_frame(std::string out, size_t iteration, int opt = VTK_WRITER)
	{
		return write_frame(out,iteration,"",opt);
	}

	/*! \brief Output particle position and properties
	 *
	 * \param out output
	 * \param iteration (we can append the number at the end of the file_name)
	 * \param meta_info meta information example ("time = 1.234" add the information time to the VTK file)
	 * \param opt VTK_WRITER, CSV_WRITER, it is also possible to choose the format for  VTK
	 *            FORMAT_BINARY. (the default is ASCII format)
	 *
	 * \return if the file has been written correctly
	 *
	 */
	inline bool write_frame(std::string out, size_t iteration, std::string meta_info, int opt = VTK_WRITER)
	{
		Vcluster<Memory> & v_cl = create_vcluster<Memory>();

		if ((opt & 0x0FFF0000) == CSV_WRITER)
		{
			// CSVWriter test
			CSVWriter<vector_dist_pos,
					  vector_dist_prop > csv_writer;

			std::string output = std::to_string(out + "_" + std::to_string(v_cl.getProcessUnitID()) + "_" + std::to_string(iteration) + std::to_string(".csv"));

			// Write the CSV
			return csv_writer.write(output, vPos, vPrp);
		}
		else
		{
			file_type ft = file_type::ASCII;

			if (opt & FORMAT_BINARY)
				ft = file_type::BINARY;

			// VTKWriter for a set of points
			VTKWriter<boost::mpl::pair<vector_dist_pos,
									   vector_dist_prop>, VECTOR_POINTS> vtk_writer;
			vtk_writer.add(vPos,vPrp,ghostMarker);

			std::string output = std::to_string(out + "_" + std::to_string(v_cl.getProcessUnitID()) + "_" + std::to_string(iteration) + std::to_string(".vtp"));

            //Create Directory for VTP files and write the PVTP metadata
            if(v_cl.rank()==0)
            {
                create_directory_if_not_exist("VTPDATA",1);
                vtk_writer.write_pvtp(out,prp_names,v_cl.size(),iteration);
            }
            v_cl.barrier();

			// Write the VTK file
			bool ret=vtk_writer.write(output,prp_names,"particles",meta_info,ft);

            return ret;
		}
	}

    /*! \brief Output particle position and properties and add a time stamp to pvtp
	 *
	 * \param out output
	 * \param iteration (we can append the number at the end of the file_name)
	 * \param time = 1.234 to add to the information time to the PVTP file)
	 * \param opt VTK_WRITER, CSV_WRITER, it is also possible to choose the format for  VTK
	 *            FORMAT_BINARY. (the default is ASCII format)
	 *
	 * \return if the file has been written correctly
	 *
	 */
	inline bool write_frame(std::string out, size_t iteration, double time, int opt = VTK_WRITER)
	{
		Vcluster<Memory> & v_cl = create_vcluster<Memory>();

		if ((opt & 0x0FFF0000) == CSV_WRITER)
		{
			// CSVWriter test
			CSVWriter<vector_dist_pos, vector_dist_prop> csv_writer;

			std::string output = std::to_string(out + "_" + std::to_string(v_cl.getProcessUnitID()) + "_" + std::to_string(iteration) + std::to_string(".csv"));

			// Write the CSV
			return csv_writer.write(output, vPos, vPrp);
		}
		else
		{
			file_type ft = file_type::ASCII;

			if (opt & FORMAT_BINARY)
				ft = file_type::BINARY;

			// VTKWriter for a set of points
			VTKWriter<boost::mpl::pair<vector_dist_pos, vector_dist_prop>, VECTOR_POINTS> vtk_writer;
			vtk_writer.add(vPos,vPrp,ghostMarker);

			std::string output = std::to_string(out + "_" + std::to_string(v_cl.getProcessUnitID()) + "_" + std::to_string(iteration) + std::to_string(".vtp"));

            //Create Directory for VTP files and write the PVTP metadata
            if(v_cl.rank()==0)
            {
                create_directory_if_not_exist("VTPDATA",1);
                vtk_writer.write_pvtp(out,prp_names,v_cl.size(),iteration,time);
            }
            v_cl.barrier();
			// Write the VTK file
			bool ret=vtk_writer.write(output,prp_names,"particles","",ft);
            return ret;
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

	/*! \brief Get the Virtual Cluster machine
	 *
	 * \return the Virtual cluster machine
	 *
	 */

	Vcluster<Memory> & getVC()
	{
#ifdef SE_CLASS2
		check_valid(this,8);
#endif
		return create_vcluster<Memory>();;
	}

	/*! \brief return the position vector of all the particles
	 *
	 * \return the particle position vector
	 *
	 */
	const vector_dist_pos & getPosVector() const
	{
		return vPos;
	}

	/*! \brief return the position vector of all the particles
	 *
	 * \return the particle position vector
	 *
	 */
	vector_dist_pos & getPosVector()
	{
		return vPos;
	}

	/*! \brief return the property vector of all the particles
	 *
	 * \return the particle property vector
	 *
	 */
	const vector_dist_prop & getPropVector() const
	{
		return vPrp;
	}

	/*! \brief return the property vector of all the particles
	 *
	 * \return the particle property vector
	 *
	 */
	vector_dist_prop & getPropVector()
	{
		return vPrp;
	}

	/*! \brief It return the sum of the particles in the previous processors
	 *
	 * \return the particles number
	 *
	 */
	size_t accum()
	{
		Vcluster<Memory> & v_cl = create_vcluster<Memory>();

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
	 * \param cellList Cell-List neighborhood
	 *
	 * \return Particle iterator
	 *
	 */
	template<typename cli> ParticleItCRS_Cells<dim,cli,decltype(vPos)> getParticleIteratorCRS_Cell(cli & cellList)
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
			shift.set_d(i,cellList.getPadding(i));

		grid_sm<dim,void> gs = cellList.getInternalGrid();

		getDecomposition().setNNParameters(shift,gs);

		// First we check that
		return ParticleItCRS_Cells<dim,cli,decltype(vPos)>(
			cellList,getDecomposition().getCRSDomainCells(),
			getDecomposition().getCRSAnomDomainCells(),
			cellList.getNNc_sym()
		);
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

    /*! \brief Get the properties names
     *
     * It is useful to get name for the properties in vtk writers
     *
     */
    openfpm::vector<std::string> &  getPropNames()
    {
        return prp_names;
    }


    /*! \brief Get a special particle iterator able to iterate across particles using
	 *         symmetric crossing scheme
	 *
	 * \param cellList Verlet list neighborhood
	 *
	 * \return Particle iterator
	 *
	 */
	template<typename VerletList_type>
	openfpm::vector_key_iterator_seq<typename VerletList_type::Mem_type_type::local_index_type> getParticleIteratorCRS(VerletList_type & cellList)
	{
#ifdef SE_CLASS1
		if (!(opt & BIND_DEC_TO_GHOST))
		{
			std::cerr << __FILE__ << ":" << __LINE__ << " Error the vector has been constructed without BIND_DEC_TO_GHOST, getParticleIteratorCRS_Cell require the vector to be constructed with BIND_DEC_TO_GHOST option " << std::endl;
			ACTION_ON_ERROR(VECTOR_DIST_ERROR_OBJECT);
		}
#endif

		// First we check that
		return openfpm::vector_key_iterator_seq<typename VerletList_type::Mem_type_type::local_index_type>(cellList.getParticleSeq());
	}

	/*! \brief Return from which cell we have to start in case of CRS interation
	 *         scheme
	 *
	 * \param cellList cell-list
	 *
	 * \return The starting cell point
	 *
	 */
	template<typename Celllist> grid_key_dx<dim> getCRSStart(Celllist & cellList)
	{
		return cellList.getStartDomainCell();
	}

	/*! \brief Return from which cell we have to stop in case of CRS interation
	 *         scheme
	 *
	 * \param cellList cell-list
	 *
	 * \return The stop cell point
	 *
	 */
	template<typename Celllist> grid_key_dx<dim> getCRSStop(Celllist & cellList)
	{
		grid_key_dx<dim> key = cellList.getStopDomainCell();

		for (size_t i = 0 ; i < dim ; i++)
			key.set_d(i,key.get(i) + 1);
		return key;
	}

	/*! \brief Indicate that this class is not a subset
	 *
	 * \return false
	 *
	 */
	bool isSubset() const
	{
		return false;
	}

#ifdef CUDA_GPU

		/*! \brief Convert the grid into a data-structure compatible for computing into GPU
		 *
		 *  The object created can be considered like a reference of the original
		 *
		 * \return an usable vector in the kernel
		 *
		 */
		template<unsigned int ... prp> vector_dist_ker<dim,St,prop,layout_base> toKernel()
		{
			vector_dist_ker<dim,St,prop,layout_base> v(ghostMarker,vPos.toKernel(), vPrp.toKernel());

			return v;
		}

		/*! \brief Return the internal vector_dist_ker_list structure
		 *
		 *
		 *
		 * \return
		 */
		vector_dist_ker_list<vector_dist_ker<dim,St,prop,layout_base>> & private_get_vector_dist_ker_list()
		{
			return *this;
		}

		/*! \brief Move the memory from the device to host memory
		 *
		 * \tparam property to move use POS_PROP for position property
		 *
		 */
		template<unsigned int ... prp> void deviceToHostProp()
		{
			vPrp.template deviceToHost<prp ...>();
		}

		/*! \brief Move the memory from the device to host memory
		 *
		 * \tparam property to move use POS_PROP for position property
		 *
		 * \param start point
		 * \param stop point (included)
		 *
		 */
		template<unsigned int ... prp> void deviceToHostProp(size_t start, size_t stop)
		{
			vPrp.template deviceToHost<prp ...>(start,stop);
		}

		/*! \brief Move the memory from the device to host memory
		 *
		 * \tparam property to move use POS_PROP for position property
		 *
		 */
		void deviceToHostPos()
		{
			vPos.template deviceToHost<0>();
		}

		/*! \brief Move the memory from the device to host memory
		 *
		 * \tparam property to move use POS_PROP for position property
		 *
		 */
		template<unsigned int ... prp> void hostToDeviceProp()
		{
			vPrp.template hostToDevice<prp ...>();
		}

		/*! \brief Move the memory from the device to host memory
		 *
		 * \tparam property to move use POS_PROP for position property
		 *
		 */
		void hostToDevicePos()
		{
			vPos.template hostToDevice<0>();
		}

		void setGhostMarker(size_t ghostMarker)
		{
			this->ghostMarker = ghostMarker;
		}

		/*! \brief this function restores order after sorting done by cell list (GPU-only)
		 *
		 * \param cellList Cell-list to use to restore the order
		 *
		 */
		template<unsigned int ...prp, typename CellList_type>
		void restoreOrder(CellList_type & cellList)
		{
			if (cellList.getOpt() & CL_GPU_REORDER_POSITION) vPos.swap(vPosReordered);
			if (cellList.getOpt() & CL_GPU_REORDER_PROPERTY) vPrp.swap(vPrpReordered);

			cellList.template restoreOrder<vector_dist_pos, vector_dist_prop, prp...>(vPosReordered, vPrpReordered, vPos, vPrp);
		}

        /*! \brief This function compare if the host and device buffer position match up to some tolerance
         *
         * \tparam prp property to check
         *
         * \param tol tollerance absolute
         *
         */
        bool compareHostAndDevicePos(St tol, St near  = -1.0, bool silent = false)
        {
			return compare_host_device<Point<dim,St>,0>::compare(vPos,tol,near,silent);
        }


        /*! \brief This function compare if the host and device buffer position match up to some tolerance
         *
         * \tparam prp property to check
         *
         * \param tol tollerance absolute
         *
         */
        template<unsigned int prp>
        bool compareHostAndDeviceProp(St tol, St near  = -1.0, bool silent = false)
        {
			return compare_host_device<typename boost::mpl::at<typename prop::type,
												boost::mpl::int_<prp> >::type,prp>::compare(vPrp,tol,near,silent);
        }

#else

		/*! \brief Move the memory from the device to host memory
		 *
		 * \tparam property to move use POS_PROP for position property
		 *
		 */
		template<unsigned int ... prp> void deviceToHostProp()
		{}

		/*! \brief Move the memory from the device to host memory
		 *
		 * \tparam property to move use POS_PROP for position property
		 *
		 */
		void deviceToHostPos()
		{}

		/*! \brief Move the memory from the device to host memory
		 *
		 * \tparam property to move use POS_PROP for position property
		 *
		 */
		template<unsigned int ... prp> void hostToDeviceProp()
		{}

		/*! \brief Move the memory from the device to host memory
		 *
		 * \tparam property to move use POS_PROP for position property
		 *
		 */
		void hostToDevicePos()
		{}

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
template<unsigned int dim, typename St, typename prop, typename Decomposition = CartDecomposition<dim,St,CudaMemory,memory_traits_inte>> using vector_dist_dev = vector_dist<dim,St,prop,Decomposition,CudaMemory,memory_traits_inte>;

#endif /* VECTOR_HPP_ */
