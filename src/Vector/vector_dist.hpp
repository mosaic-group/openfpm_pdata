/*
 * Vector.hpp
 *
 *  Created on: Mar 5, 2015
 *      Author: Pietro Incardona
 */

#ifndef VECTOR_HPP_
#define VECTOR_HPP_

#include "../../openfpm_data/src/NN/CellList/CellListFast_gen.hpp"
#include "../../openfpm_data/src/NN/CellList/CellListFast_gen.hpp"
#include "HDF5_XdmfWriter/HDF5_XdmfWriter.hpp"
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
#include "NN/CellList/ProcKeys.hpp"

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

#define DEC_GRAN(gr) ((size_t)gr << 32)

// Perform a ghost get or a ghost put
#define GET	1
#define PUT 2

// Write the particles with ghost
#define NO_GHOST 0
#define WITH_GHOST 2

//! General function t get a cell-list
template<unsigned int dim, typename St, typename CellL, typename Vector>
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
		return vd.getCellList(r_cut);
	}
};

//! General function t get a cell-list
template<unsigned int dim, typename St, typename Vector, typename Mem_type>
struct gcl<dim,St,CellList_gen<dim, St, Process_keys_hilb,Mem_type, shift<dim, St> >,Vector>
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
	static inline CellList_gen<dim, St, Process_keys_hilb, Mem_type, shift<dim, St> > get(Vector & vd, const St & r_cut, const Ghost<dim,St> & g)
	{
		return vd.getCellList_hilb(r_cut,g);
	}
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
 *
 */

template<unsigned int dim, typename St, typename prop, typename Decomposition = CartDecomposition<dim,St>, typename Memory = HeapMemory>
class vector_dist : public vector_dist_comm<dim,St,prop,Decomposition,Memory>
{
public:

	//! Self type
	typedef vector_dist<dim,St,prop,Decomposition,Memory> self;

	//! property object
	typedef prop value_type;

private:

	//! Ghost marker, all the particle with id > g_m are ghost all with g_m < are real particle
	size_t g_m = 0;

	//! Particle position vector, (It has 2 elements) the first has real particles assigned to a processor
	//! the second element contain unassigned particles
	openfpm::vector<Point<dim, St>> v_pos;

	//! Particle properties vector, (It has 2 elements) the first has real particles assigned to a processor
	//! the second element contain unassigned particles
	openfpm::vector<prop> v_prp;

	//! Virtual cluster
	Vcluster & v_cl;

	//! option used to create this vector
	size_t opt = 0;

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
	vector_dist<dim,St,prop,Decomposition,Memory> & operator=(const vector_dist<dim,St,prop,Decomposition,Memory> & v)
	{
		static_cast<vector_dist_comm<dim,St,prop,Decomposition,Memory> *>(this)->operator=(static_cast<vector_dist_comm<dim,St,prop,Decomposition,Memory>>(v));

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
	vector_dist<dim,St,prop,Decomposition,Memory> & operator=(vector_dist<dim,St,prop,Decomposition,Memory> && v)
	{
		static_cast<vector_dist_comm<dim,St,prop,Decomposition,Memory> *>(this)->operator=(static_cast<vector_dist_comm<dim,St,prop,Decomposition,Memory> >(v));

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
	vector_dist(const vector_dist<dim,St,prop,Decomposition,Memory> & v)
	:vector_dist_comm<dim,St,prop,Decomposition,Memory>(v.getDecomposition()),v_cl(v.v_cl) SE_CLASS3_VDIST_CONSTRUCTOR
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
	vector_dist(vector_dist<dim,St,prop,Decomposition,Memory> && v) noexcept
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
	vector_dist_comm<dim,St,prop,Decomposition,Memory>(dec), v_cl(create_vcluster()) SE_CLASS3_VDIST_CONSTRUCTOR
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
	 * \param opt additional options. BIND_DEC_TO_GHOST Bind the decomposition to be multiple of the
	 *          ghost size. This is required if we want to use symmetric to eliminate
	 *          ghost communications.
	 *
	 */
	vector_dist(size_t np, Box<dim, St> box, const size_t (&bc)[dim], const Ghost<dim, St> & g, size_t opt = 0, const grid_sm<dim,void> & gdist = grid_sm<dim,void>())
	:v_cl(create_vcluster()),opt(opt) SE_CLASS3_VDIST_CONSTRUCTOR
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
	template<unsigned int id> inline auto getProp(vect_dist_key_dx vec_key) const -> const decltype(v_prp.template get<id>(vec_key.getKey()))
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
	template<unsigned int id> inline auto getProp(size_t vec_key) -> decltype(v_prp.template get<id>(vec_key))
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
	template<unsigned int id> inline auto getProp(size_t vec_key) const -> const decltype(v_prp.template get<id>(vec_key))
	{
		return v_prp.template get<id>(vec_key);
	}

#endif

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
	template<unsigned int id> inline auto getPropRead(vect_dist_key_dx vec_key) const -> const decltype(v_prp.template get<id>(vec_key.getKey()))
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

////////////////////////////// READ AND WRITE VERSION //////////

	/*! \brief Get the position of the last element
	 *
	 * \return the position of the element in space
	 *
	 */
	inline auto getLastPosRead() -> decltype(v_pos.template get<0>(0))
	{
#ifdef SE_CLASS3
		se3.read<prop::max_prop_real>(*this,g_m-1);
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
		se3.write<prop::max_prop_real>(*this,g_m-1);
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
		se3.write<id>(*this,g_m-1);
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
	template<typename CellL = CellList<dim, St, Mem_fast<dim,St>, shift<dim, St> > > CellL getCellListSym(St r_cut)
	{
#ifdef SE_CLASS1
		if ((opt & BIND_DEC_TO_GHOST))
		{
			std::cerr << __FILE__ << ":" << __LINE__ << " error to get symmetric cell-list you must construct the vector with the option BIND_DEC_TO_GHOST " << std::endl;
			ACTION_ON_ERROR(VECTOR_DIST_ERROR_OBJECT);
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
	template<typename CellL = CellList<dim, St, FAST, shift<dim, St> > > CellL getCellListSymNoBind(St r_cut)
	{
		return getCellList(r_cut);
	}


	/*! \brief Construct a cell list starting from the stored particles
	 *
	 * \tparam CellL CellList type to construct
	 *
	 * \param r_cut interation radius, or size of each cell
	 *
	 * \return the Cell list
	 *
	 */
	template<typename CellL = CellList_gen<dim, St, Process_keys_lin, Mem_fast<dim,St>, shift<dim, St> > > CellL getCellList(St r_cut)
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

		return getCellList(r_cut, g);
	}

	/*! \brief Construct an hilbert cell list starting from the stored particles
	 *
	 * \tparam CellL CellList type to construct
	 *
	 * \param r_cut interation radius, or size of each cell
	 *
	 * \return the Cell list
	 *
	 */
	template<typename CellL = CellList_gen<dim, St, Process_keys_hilb, Mem_fast<dim,St>, shift<dim, St> > > CellL getCellList_hilb(St r_cut)
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
	 *
	 */
	template<typename CellL> void updateCellList(CellL & cell_list)
	{
#ifdef SE_CLASS3
		se3.getNN();
#endif

		// This function assume equal spacing in all directions
		// but in the worst case we take the maximum
		St r_cut = 0;
		for (size_t i = 0 ; i < dim ; i++)
			r_cut = std::max(r_cut,cell_list.getCellBox().getHigh(i));

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
			CellL cli_tmp = gcl<dim,St,CellL,self>::get(*this,r_cut,getDecomposition().getGhost());

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
	template<typename CellL = CellList<dim, St, Mem_fast<dim,St>, shift<dim, St> > > void updateCellListSym(CellL & cell_list)
	{
#ifdef SE_CLASS3
		se3.getNN();
#endif

		// This function assume equal spacing in all directions
		// but in the worst case we take the maximum
		St r_cut = 0;
		for (size_t i = 0 ; i < dim ; i++)
			r_cut = std::max(r_cut,cell_list.getCellBox().getHigh(i));

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
			CellL cli_tmp = gcl<dim,St,CellL,self>::get(*this,r_cut,getDecomposition().getGhost());

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
	 *
	 * \return the CellList
	 *
	 */
	template<typename CellL = CellList_gen<dim, St, Process_keys_lin, Mem_fast<dim,St>, shift<dim, St> > > CellL getCellList(St r_cut, const Ghost<dim, St> & enlarge)
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
		cl_param_calculate(pbox, div, r_cut, enlarge);

		cell_list.Initialize(pbox, div, g_m);
		cell_list.set_ndec(getDecomposition().get_ndec());

		updateCellList(cell_list);

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
	template<typename CellL = CellList_gen<dim, St, Process_keys_hilb, Mem_fast<dim,St>, shift<dim, St> > > CellL getCellList_hilb(St r_cut, const Ghost<dim, St> & enlarge)
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

		cell_list.Initialize(pbox, div, g_m);
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
	VerletList<dim,St,FAST,shift<dim,St> > getVerletSym(St r_cut)
	{
#ifdef SE_CLASS3
		se3.getNN();
#endif

		VerletList<dim,St,FAST,shift<dim,St>> ver;

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
	VerletList<dim,St,FAST,shift<dim,St> > getVerletCrs(St r_cut)
	{
#ifdef SE_CLASS3
		se3.getNN();
#endif

		VerletList<dim,St,FAST,shift<dim,St>> ver;

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
	VerletList<dim,St,FAST,shift<dim,St> > getVerlet(St r_cut)
	{
#ifdef SE_CLASS3
		se3.getNN();
#endif

		VerletList<dim,St,FAST,shift<dim,St>> ver;

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
	void updateVerlet(VerletList<dim,St,FAST,shift<dim,St> > & ver, St r_cut, size_t opt = VL_NON_SYMMETRIC)
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
				VerletList<dim,St,FAST,shift<dim,St> > ver_tmp;

				ver_tmp = getVerlet(r_cut);
				ver.swap(ver);
			}
		}
		else if (opt == VL_CRS_SYMMETRIC)
		{
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
				VerletList<dim,St,FAST,shift<dim,St> > ver_tmp;

				ver_tmp = getVerletCrs(r_cut);
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
				VerletList<dim,St,FAST,shift<dim,St> > ver_tmp;

				ver_tmp = getVerlet(r_cut);
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
	template<typename CellL=CellList_gen<dim,St,Process_keys_lin,Mem_fast<dim,St>,shift<dim,St> > > void reorder (int32_t m)
	{
		reorder(m,getDecomposition().getGhost());
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
	template<typename CellL=CellList_gen<dim,St,Process_keys_lin,Mem_fast<dim,St>,shift<dim,St> > > void reorder(int32_t m, const Ghost<dim,St> & enlarge)
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

		cell_list.Initialize(pbox,div,g_m);

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

		v_pos_dest.resize(v_pos.size());
		v_prp_dest.resize(v_prp.size());

		//hilberts curve iterator
		grid_key_dx_iterator_hilbert<dim> h_it(m);

		//Index for v_pos_dest
		size_t count = 0;

		grid_key_dx<dim> ksum;

		for (size_t i = 0; i < dim ; i++)
			ksum.set_d(i,cell_list.getPadding(i));

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
	 *
	 * \return an iterator
	 *
	 */
	template<typename CellList> ParticleIt_Cells<dim,CellList> getDomainIteratorCells(CellList & NN)
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
		return vector_dist_comm<dim,St,prop,Decomposition,Memory>::getDecomposition();
	}

	/*! \brief Get the decomposition
	 *
	 * \return
	 *
	 */
	inline const Decomposition & getDecomposition() const
	{
		return vector_dist_comm<dim,St,prop,Decomposition,Memory>::getDecomposition();
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
	 *
	 */
	template<unsigned int ... prp> void map_list()
	{
		this->template map_list_<prp...>(v_pos,v_prp,g_m);
	}


	/*! \brief It move all the particles that does not belong to the local processor to the respective processor
	 *
	 * \tparam out of bound policy it specify what to do when the particles are detected out of bound
	 *
	 * In general this function is called after moving the particles to move the
	 * elements out the local processor. Or just after initialization if each processor
	 * contain non local particles
	 *
	 *
	 */
	template<typename obp = KillParticle> void map()
	{
#ifdef SE_CLASS3
		se3.map_pre();
#endif

		this->template map_<obp>(v_pos,v_prp,g_m);

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
	 * \param ts It is an optional parameter approximately should be the number of ghost get between two
	 *           rebalancing at first decomposition this number can be ignored (default = 1) because not used
	 *
	 */
	template <typename Model=ModelLin>inline void addComputationCosts(Model md=Model(), size_t ts = 1)
	{
		CellDecomposer_sm<dim, St, shift<dim,St>> cdsm;

		Decomposition & dec = getDecomposition();
		auto & dist = getDecomposition().getDistribution();

		cdsm.setDimensions(dec.getDomain(), dec.getDistGrid().getSize(), 0);

		for (size_t i = 0; i < dist.getNOwnerSubSubDomains() ; i++)
			dec.setSubSubDomainComputationCost(dist.getOwnerSubSubDomain(i) , 1);

		auto it = getDomainIterator();

		while (it.isNext())
		{
			size_t v = cdsm.getCell(this->getPos(it.get()));

			md.addComputation(dec,*this,v,it.get().getKey());

			++it;
		}

		dec.computeCommunicationAndMigrationCosts(ts);

		// Go throught all the sub-sub-domains and apply the model

		for (size_t i = 0 ; i < dist.getNOwnerSubSubDomains(); i++)
			md.applyModel(dec,dist.getOwnerSubSubDomain(i));

		dist.setDistTol(md.distributionTol());
	}

	inline void save(const std::string & filename) const
	{
		//Pack_request vector
		size_t req = 0;

		//std::cout << "V_pos.size() before save: " << v_pos.size() << std::endl;

		//Pack request
		Packer<decltype(v_pos),HeapMemory>::packRequest(v_pos,req);
		Packer<decltype(v_prp),HeapMemory>::packRequest(v_prp,req);

		std::cout << "Req: " << req << std::endl;

		// allocate the memory
		HeapMemory pmem;
		//pmem.allocate(req);
		ExtPreAlloc<HeapMemory> & mem = *(new ExtPreAlloc<HeapMemory>(req,pmem));
		mem.incRef();

		//Packing

		Pack_stat sts;

		Packer<decltype(v_pos),HeapMemory>::pack(mem,v_pos,sts);
		Packer<decltype(v_prp),HeapMemory>::pack(mem,v_prp,sts);

	    /*****************************************************************
	     * Create a new file with default creation and access properties.*
	     * Then create a dataset and write data to it and close the file *
	     * and dataset.                                                  *
	     *****************************************************************/

		int mpi_rank = v_cl.getProcessUnitID();
		int mpi_size = v_cl.getProcessingUnits();

		MPI_Comm comm = v_cl.getMPIComm();
		MPI_Info info  = MPI_INFO_NULL;
/*
	    //Initialize MPI
	    MPI_Comm_size(comm, &mpi_size);
	    MPI_Comm_rank(comm, &mpi_rank);
*/
	    // Set up file access property list with parallel I/O access

		hid_t plist_id = H5Pcreate(H5P_FILE_ACCESS);
		H5Pset_fapl_mpio(plist_id, comm, info);

		// Create a new file collectively and release property list identifier.
		hid_t file = H5Fcreate (filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
		H5Pclose(plist_id);

		size_t sz = pmem.size();
		//std::cout << "Pmem.size: " << pmem.size() << std::endl;
		openfpm::vector<size_t> sz_others;
		v_cl.allGather(sz,sz_others);
		v_cl.execute();

		size_t sum = 0;

		for (size_t i = 0; i < sz_others.size(); i++)
			sum += sz_others.get(i);

		//Size for data space in file
		hsize_t fdim[1] = {sum};

		//Size for data space in file
		hsize_t fdim2[1] = {(size_t)mpi_size};

		//Create data space in file
		hid_t file_dataspace_id = H5Screate_simple(1, fdim, NULL);

		//Create data space in file
		hid_t file_dataspace_id_2 = H5Screate_simple(1, fdim2, NULL);

		//Size for data space in memory
		hsize_t mdim[1] = {pmem.size()};

		//Create data space in memory
		hid_t mem_dataspace_id = H5Screate_simple(1, mdim, NULL);

		if (mpi_rank == 0)
			std::cout << "Total object size: " << sum << std::endl;

		//Create data set in file
		hid_t file_dataset = H5Dcreate (file, "vector_dist", H5T_NATIVE_CHAR, file_dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

		//Create data set 2 in file
		hid_t file_dataset_2 = H5Dcreate (file, "metadata", H5T_NATIVE_INT, file_dataspace_id_2, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

	    //H5Pclose(plist_id);
	    H5Sclose(file_dataspace_id);
	    H5Sclose(file_dataspace_id_2);

	    hsize_t block[1] = {pmem.size()};

	    //hsize_t stride[1] = {1};

		hsize_t count[1] = {1};

	    hsize_t offset[1] = {0};

	    for (int i = 0; i < mpi_rank; i++)
	    {
	    	if (mpi_rank == 0)
				offset[0] = 0;
	    	else
	    		offset[0] += sz_others.get(i);
	    }

	    std::cout << "MPI rank: " << mpi_rank << ", MPI size: " << mpi_size << ", Offset: " << offset[0] << ", Block: " << block[0] << std::endl;

	    int metadata[mpi_size];

	    for (int i = 0; i < mpi_size; i++)
	    	metadata[i] = sz_others.get(i);

	    //Select hyperslab in the file.
	    file_dataspace_id = H5Dget_space(file_dataset);
	    H5Sselect_hyperslab(file_dataspace_id, H5S_SELECT_SET, offset, NULL, count, block);

	    file_dataspace_id_2 = H5Dget_space(file_dataset_2);


	    //Create property list for collective dataset write.
	    plist_id = H5Pcreate(H5P_DATASET_XFER);
	    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);

		//Write a data set to a file
		H5Dwrite(file_dataset, H5T_NATIVE_CHAR, mem_dataspace_id, file_dataspace_id, plist_id, (const char *)pmem.getPointer());

		//Write a data set 2 to a file
		H5Dwrite(file_dataset_2, H5T_NATIVE_INT, H5S_ALL, file_dataspace_id_2, plist_id, metadata);


	    //Close/release resources.
	    H5Dclose(file_dataset);
	    H5Sclose(file_dataspace_id);
	    H5Dclose(file_dataset_2);
	    H5Sclose(file_dataspace_id_2);
	    H5Sclose(mem_dataspace_id);
	    H5Pclose(plist_id);
	    H5Fclose(file);
	}

	void load_block(long int bid, hssize_t mpi_size_old, int * metadata_out, openfpm::vector<size_t> metadata_accum, hid_t plist_id, hid_t dataset_2)
	{
/*	  	if (mpi_size >= mpi_size_old)
	  	{
			if (mpi_rank >= mpi_size_old)
				block[0] = 0;
			else
				block[0] = {(size_t)metadata_out[mpi_rank]};
	  	}
	  	else
	  	{
	  		int x = mpi_size_old/mpi_size;
	  		int shift = mpi_rank*x;
	  		for (int i = 0; i < x; i++)
	  		{
	  			//block0.get(mpi_rank).add(metadata_out[shift]);
	  			block[0] += metadata_out[shift];
	  			shift++;
	  		}
	  		int y = mpi_size_old%mpi_size;
	  		if (mpi_rank < y)
	  		{
				block_add[0] += metadata_out[mpi_size*x+mpi_rank];
				//block_add0.get(mpi_rank).add(metadata_out[mpi_size*x+mpi_rank]);
	  		}
	  	}*/

//		std::cout << "BID: " << bid << std::endl;

		hsize_t offset[1];
		hsize_t block[1];

		if (bid < mpi_size_old && bid != -1)
		{
			offset[0] = metadata_accum.get(bid);
			block[0] = metadata_out[bid];
		}
		else
		{
			offset[0] = 0;
			block[0] = 0;
		}

//		std::cout << "Offset: " << offset[0] << "; Block: " << block[0]<<  std::endl;
//	    hsize_t offset_add[1] = {0};

/*	    if (mpi_size >= mpi_size_old)
		{
			if (mpi_rank >= mpi_size_old)
				offset[0] = 0;
			else
			{
				for (int i = 0; i < mpi_rank; i++)
				offset[0] += metadata_out[i];
			}
		}
	    else
	    {
	  		int x = mpi_size_old/mpi_size;
	  		int shift = mpi_rank*x;

	  		for (int i = 0; i < shift; i++)
	  		{
	  			offset[0] += metadata_out[i];
	  			//offset0.get(mpi_rank).add(metadata_out[i]);
	  		}

	  		int y = mpi_size_old%mpi_size;
	  		if (mpi_rank < y)
	  		{
	  			for (int i = 0; i < mpi_size*x + mpi_rank; i++)
	  			{
	  				offset_add[0] += metadata_out[i];
	  				//offset_add0.get(mpi_rank).add(metadata_out[i]);
	  			}
	  		}
	    }*/

	    //hsize_t stride[1] = {1};
	    hsize_t count[1] = {1};

	    //std::cout << "LOAD: MPI rank: " << mpi_rank << ", MPI size: " << mpi_size << ", Offset: " << offset[0] << ", Offset_add: " << offset_add[0] << ", Block: " << block[0] << ", Block_add: " << block_add[0] << std::endl;
/*
	    std::cout << "LOAD: MPI rank: " << mpi_rank << ", MPI size: " << mpi_size << std::endl;
	    for (size_t i = 0; i < offset0.get(mpi_rank).size(); i++)
	    	std::cout << ", Offset: " << offset0.get(mpi_rank).get(i) << std::endl;
		for (size_t i = 0; i < offset_add0.get(mpi_rank).size(); i++)
			std::cout << ", Offset_add: " << offset_add0.get(mpi_rank).get(i) << std::endl;
		for (size_t i = 0; i < block0.get(mpi_rank).size(); i++)
			std::cout << ", Block: " << block0.get(mpi_rank).get(i) << std::endl;
		for (size_t i = 0; i < block_add0.get(mpi_rank).size(); i++)
			std::cout << ", Block_add: " << block_add0.get(mpi_rank).get(i) << std::endl;
*/

		//Select file dataspace
		hid_t file_dataspace_id_2 = H5Dget_space(dataset_2);

        H5Sselect_hyperslab(file_dataspace_id_2, H5S_SELECT_SET, offset, NULL, count, block);

		//Select file dataspace
/*		hid_t file_dataspace_id_3 = H5Dget_space(dataset_2);

        H5Sselect_hyperslab(file_dataspace_id_3, H5S_SELECT_SET, offset_add, NULL, count, block_add);*/

        hsize_t mdim_2[1] = {block[0]};
//        hsize_t mdim_3[1] = {block_add[0]};


		//Size for data space in memory

		/*if (mpi_rank >= mpi_size_old)
			mdim_2[0] = 0;
		else
			mdim_2[0] = metadata_out[mpi_rank];*/

		//Create data space in memory
		hid_t mem_dataspace_id_2 = H5Screate_simple(1, mdim_2, NULL);
//		hid_t mem_dataspace_id_3 = H5Screate_simple(1, mdim_3, NULL);

		//if (mpi_rank == 0)
/*		{
			hssize_t size2;

			size2 = H5Sget_select_npoints (mem_dataspace_id_2);
			printf ("\nLOAD: memspace_id_2 size: %llu\n", size2);
			size2 = H5Sget_select_npoints (file_dataspace_id_2);
			printf ("LOAD: dataspace_id_2 size: %llu\n", size2);
		}*/
/*
		if (mpi_rank == 0)
		{
			hssize_t size2;

			size2 = H5Sget_select_npoints (mem_dataspace_id_3);
			printf ("\nLOAD: memspace_id_3 size: %llu\n", size2);
			size2 = H5Sget_select_npoints (file_dataspace_id_3);
			printf ("LOAD: dataspace_id_3 size: %llu\n", size2);
		}
*/
	size_t sum = 0;

		for (int i = 0; i < mpi_size_old; i++)
		{
			sum += metadata_out[i];
		}

//		std::cout << "LOAD: sum: " << sum << std::endl;

		// allocate the memory
		HeapMemory pmem;
//		HeapMemory pmem2;
		//pmem.allocate(req);
		ExtPreAlloc<HeapMemory> & mem = *(new ExtPreAlloc<HeapMemory>(block[0],pmem));
		mem.incRef();
//		ExtPreAlloc<HeapMemory> & mem2 = *(new ExtPreAlloc<HeapMemory>(block_add[0],pmem2));
//		mem2.incRef();

	  	// Read the dataset.
	    H5Dread(dataset_2, H5T_NATIVE_CHAR, mem_dataspace_id_2, file_dataspace_id_2, plist_id, (char *)mem.getPointer());

	    // Read the dataset.
//		H5Dread(dataset_2, H5T_NATIVE_CHAR, mem_dataspace_id_3, file_dataspace_id_3, plist_id, (char *)mem2.getPointer());

		mem.allocate(pmem.size());
//		mem2.allocate(pmem2.size());
//		std::cout << "Mem.size(): " << mem.size() << " = " << block[0] << std::endl;

		Unpack_stat ps;

		openfpm::vector<Point<dim, St>> v_pos_unp;

		openfpm::vector<prop> v_prp_unp;

		Unpacker<decltype(v_pos_unp),HeapMemory>::unpack(mem,v_pos_unp,ps,1);
		Unpacker<decltype(v_prp_unp),HeapMemory>::unpack(mem,v_prp_unp,ps,1);

//		Unpack_stat ps2;


//		Unpacker<decltype(v_pos),HeapMemory>::unpack(mem2,v_pos_unp,ps2,1);
//		Unpacker<decltype(v_prp),HeapMemory>::unpack(mem2,v_prp_unp,ps2,1);

//		std::cout << "V_pos.size(): " << v_pos.size() << std::endl;
//		std::cout << "V_pos_unp.size(): " << v_pos_unp.size() << std::endl;

		mem.decRef();
		delete &mem;

		for (size_t i = 0; i < v_pos_unp.size(); i++)
			v_pos.add(v_pos_unp.get(i));

		for (size_t i = 0; i < v_prp_unp.size(); i++)
			v_prp.add(v_prp_unp.get(i));

		g_m = v_pos.size();
	}

	inline void load(const std::string & filename)
	{
		v_pos.clear();
		v_prp.clear();

		g_m = 0;

		MPI_Comm comm = v_cl.getMPIComm();
		MPI_Info info  = MPI_INFO_NULL;

		int mpi_rank = v_cl.getProcessUnitID();
		//int mpi_size = v_cl.getProcessingUnits();

		// Set up file access property list with parallel I/O access
		hid_t plist_id = H5Pcreate(H5P_FILE_ACCESS);
		H5Pset_fapl_mpio(plist_id, comm, info);

		//Open a file
	    hid_t file = H5Fopen (filename.c_str(), H5F_ACC_RDONLY, plist_id);
	    H5Pclose(plist_id);

	    //Open dataset
	    hid_t dataset = H5Dopen (file, "metadata", H5P_DEFAULT);

	    //Create property list for collective dataset read
	  	plist_id = H5Pcreate(H5P_DATASET_XFER);
	  	H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);

		//Select file dataspace
		hid_t file_dataspace_id = H5Dget_space(dataset);

		hssize_t mpi_size_old = H5Sget_select_npoints (file_dataspace_id);
/*
		if (mpi_rank == 0)
			printf ("\nOld MPI size: %llu\n", mpi_size_old);
*/
	  	//Where to read metadata
	  	int metadata_out[mpi_size_old];

	  	for (int i = 0; i < mpi_size_old; i++)
	  	{
	  		metadata_out[i] = 0;
	  	}

		//Size for data space in memory
		hsize_t mdim[1] = {(size_t)mpi_size_old};

		//Create data space in memory
		hid_t mem_dataspace_id = H5Screate_simple(1, mdim, NULL);
/*
		if (mpi_rank == 0)
		{
			hssize_t size;

			size = H5Sget_select_npoints (mem_dataspace_id);
			printf ("LOAD: memspace_id size: %llu\n", size);
			size = H5Sget_select_npoints (file_dataspace_id);
			printf ("LOAD: dataspace_id size: %llu\n", size);
		}
*/
	  	// Read the dataset.
	    H5Dread(dataset, H5T_NATIVE_INT, mem_dataspace_id, file_dataspace_id, plist_id, metadata_out);
/*
		if (mpi_rank == 0)
		{
			std::cout << "Metadata_out[]: ";
			for (int i = 0; i < mpi_size_old; i++)
			{
				std::cout << metadata_out[i] << " ";
			}
			std::cout << " " << std::endl;
		}
*/
	    openfpm::vector<size_t> metadata_accum;
	    metadata_accum.resize(mpi_size_old);

	    metadata_accum.get(0) = 0;
	    for (int i = 1 ; i < mpi_size_old ; i++)
	    	metadata_accum.get(i) = metadata_accum.get(i-1) + metadata_out[i-1];

	    //Open dataset
	    hid_t dataset_2 = H5Dopen (file, "vector_dist", H5P_DEFAULT);

	    //Create property list for collective dataset read
	  	plist_id = H5Pcreate(H5P_DATASET_XFER);
	  	H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);

/*
	  	openfpm::vector<openfpm::vector<hsize_t>> block0;
	  	openfpm::vector<openfpm::vector<hsize_t>> block_add0;
	  	openfpm::vector<openfpm::vector<hsize_t>> offset0;
	  	openfpm::vector<openfpm::vector<hsize_t>> offset_add0;

	  	block0.resize(mpi_size);
	  	offset0.resize(mpi_size);
	  	block_add0.resize(mpi_size);
		offset_add0.resize(mpi_size);
*/
	  	openfpm::vector<size_t> n_block;
	  	n_block.resize(v_cl.getProcessingUnits());


	  	for(size_t i = 0 ; i < n_block.size() ; i++)
	  		n_block.get(i) = mpi_size_old / v_cl.getProcessingUnits();

	  	size_t rest_block = mpi_size_old % v_cl.getProcessingUnits();

	  	//std::cout << "MPI size old: " << mpi_size_old << std::endl;
	  	//std::cout << "MPI size: " << v_cl.getProcessingUnits() << std::endl;


	  //	std::cout << "Rest block: " << rest_block << std::endl;

	  	size_t max_block;

	  	if (rest_block != 0)
	  		max_block = n_block.get(0) + 1;
	  	else
	  		max_block = n_block.get(0);

	  	//for(size_t i = 0 ; i < n_block.size() ; i++)
	  	for(size_t i = 0 ; i < rest_block ; i++)
	  		n_block.get(i) += 1;


	  	//for(size_t i = 0 ; i < n_block.size() ; i++)
	  		//std::cout << "n_block.get(i): " << n_block.get(i) << std::endl;

	  	size_t start_block = 0;
	  	size_t stop_block = 0;


	  	if (v_cl.getProcessUnitID() != 0)
	  	{
			for(size_t i = 0 ; i < v_cl.getProcessUnitID() ; i++)
				start_block += n_block.get(i);
	  	}

	  	stop_block = start_block + n_block.get(v_cl.getProcessUnitID());

	  	//std::cout << "ID: " << v_cl.getProcessUnitID() << "; Start block: " << start_block << "; " << "Stop block: " << stop_block << std::endl;

	  	if (mpi_rank >= mpi_size_old)
	  		load_block(start_block,mpi_size_old,metadata_out,metadata_accum,plist_id,dataset_2);
	  	else
	  	{
	  		size_t n_bl = 0;
	  		size_t lb = start_block;
			for ( ; lb < stop_block ; lb++, n_bl++)
				load_block(lb,mpi_size_old,metadata_out,metadata_accum,plist_id,dataset_2);

			if (n_bl < max_block)
				load_block(-1,mpi_size_old,metadata_out,metadata_accum,plist_id,dataset_2);
	  	}

	    // Close the dataset.
	    H5Dclose(dataset);
	    H5Dclose(dataset_2);
	    // Close the file.
	    H5Fclose(file);
	    H5Pclose(plist_id);

		//std::cout << "V_pos.size() after merge: " << v_pos.size() << std::endl;

		// Map particles
		//map();

		//std::cout << "V_pos.size() after merge and map: " << v_pos.size() << std::endl;
	}

	/*! \brief Output particle position and properties
	 *
	 * \param out output
	 * \param opt VTK_WRITER or CSV_WRITER
	 *
	 * \return true if the file has been written without error
	 *
	 */
	inline bool write(std::string out, int opt = VTK_WRITER)
	{

		if ((opt & 0x0FFF0000) == CSV_WRITER)
		{
			// CSVWriter test
			CSVWriter<openfpm::vector<Point<dim,St>>, openfpm::vector<prop> > csv_writer;

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
			VTKWriter<boost::mpl::pair<openfpm::vector<Point<dim,St>>, openfpm::vector<prop>>, VECTOR_POINTS> vtk_writer;
			vtk_writer.add(v_pos,v_prp,g_m);

			std::string output = std::to_string(out + "_" + std::to_string(v_cl.getProcessUnitID()) + std::to_string(".vtk"));

			// Write the VTK file
			return vtk_writer.write(output,"particles",ft);
		}

		return false;
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
	 * \param opt NO_GHOST or WITH_GHOST
	 *
	 * \return if the file has been written correctly
	 *
	 */
	inline bool write(std::string out, size_t iteration, int opt = NO_GHOST)
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
			return vtk_writer.write(output,"particles",ft);
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

	Vcluster & getVC()
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
			sz += accu.get(i);

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

	/*! \brief Get a special particle iterator able to iterate across particles using
	 *         symmetric crossing scheme
	 *
	 * \param NN Verlet list neighborhood
	 *
	 * \return Particle iterator
	 *
	 */
	template<typename vrl> openfpm::vector_key_iterator_seq<typename vrl::local_index_t> getParticleIteratorCRS(vrl & NN)
	{
		// First we check that
		return openfpm::vector_key_iterator_seq<typename vrl::local_index_t>(NN.getParticleSeq());
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

#ifdef SE_CLASS3

	se_class3_vector<prop::max_prop,dim,St,Decomposition,self> & get_se_class3()
	{
		return se3;
	}

#endif
};


#endif /* VECTOR_HPP_ */
