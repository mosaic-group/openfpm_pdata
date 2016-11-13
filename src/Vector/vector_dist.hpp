/*
 * Vector.hpp
 *
 *  Created on: Mar 5, 2015
 *      Author: Pietro Incardona
 */

#ifndef VECTOR_HPP_
#define VECTOR_HPP_

#include "HDF5_XdmfWriter/HDF5_XdmfWriter.hpp"
#include "VCluster.hpp"
#include "Space/Shape/Point.hpp"
#include "Vector/vector_dist_iterator.hpp"
#include "Space/Shape/Box.hpp"
#include "Vector/vector_dist_key.hpp"
#include "memory/PtrMemory.hpp"
#include "NN/CellList/CellList.hpp"
#include "NN/CellList/CellListFast_hilb.hpp"
#include "util/common.hpp"
#include "util/object_util.hpp"
#include "memory/ExtPreAlloc.hpp"
#include "CSVWriter/CSVWriter.hpp"
#include "VTKWriter/VTKWriter.hpp"
#include "Decomposition/common.hpp"
#include "Grid/grid_dist_id_iterator_dec.hpp"
#include "Grid/grid_key_dx_iterator_hilbert.hpp"
#include "Vector/vector_dist_ofb.hpp"
#include "Decomposition/CartDecomposition.hpp"
#include "data_type/aggregate.hpp"
#include "NN/VerletList/VerletList.hpp"
#include "vector_dist_comm.hpp"

#define NO_ID false
#define ID true

// Perform a ghost get or a ghost put
#define GET	1
#define PUT 2

// Write the particles with ghost
#define NO_GHOST 0
#define WITH_GHOST 2

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

		return *this;
	}

	/*! \brief Copy Constructor
	 *
	 * \param v vector to copy
	 *
	 */
	vector_dist(const vector_dist<dim,St,prop,Decomposition,Memory> & v)
	:vector_dist_comm<dim,St,prop,Decomposition,Memory>(v.getDecomposition()),v_cl(v.v_cl)
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
	:v_cl(v.v_cl)
	{
#ifdef SE_CLASS2
		check_new(this,8,VECTOR_DIST_EVENT,4);
#endif

		this->operator=(v);
	}

	/*! \brief Constructor with predefined decomposition
	 *
	 * \param dec is the decomposition
	 * \param np number of particles
	 *
	 */
	vector_dist(const Decomposition & dec, size_t np) :
	vector_dist_comm<dim,St,prop,Decomposition,Memory>(dec), v_cl(create_vcluster())
	{
#ifdef SE_CLASS2
		check_new(this,8,VECTOR_DIST_EVENT,4);
#endif

		init_structures(np);
	}


	/*! \brief Constructor
	 *
	 * \param np number of elements
	 * \param box domain where the vector of elements live
	 * \param bc boundary conditions
	 * \param g Ghost margins
	 * \param opt additional options.
	 *        * BIND_DEC_TO_GHOST Bind the decomposition to be multiple of the
	 *          ghost size. This is required if we want to use symmetric to eliminate
	 *          ghost communications.
	 *
	 */
	vector_dist(size_t np, Box<dim, St> box, const size_t (&bc)[dim], const Ghost<dim, St> & g, size_t opt = 0)
	:v_cl(create_vcluster())
	{
#ifdef SE_CLASS2
		check_new(this,8,VECTOR_DIST_EVENT,4);
#endif

		init_structures(np);
		this->init_decomposition(box,bc,g,opt);
	}

	~vector_dist()
	{
#ifdef SE_CLASS2
		check_delete(this);
#endif
	}

	/*! \brief return the local size of the vector
	 *
	 * \return local size
	 *
	 */
	size_t size_local()
	{
		return g_m;
	}

	/*! \brief return the local size of the vector
	 *
	 * \return local size
	 *
	 */
	size_t size_local_with_ghost()
	{
		return v_pos.size();
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

	/*! \brief Construct a cell list symmetric based on a cut of radius
	 *
	 * \tparam CellL CellList type to construct
	 *
	 * \param r_cut interation radius, or size of each cell
	 *
	 * \return the Cell list
	 *
	 */
	template<typename CellL = CellList<dim, St, FAST, shift<dim, St> > > CellL getCellListSym(St r_cut)
	{
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

		updateCellList(cell_list);

		return cell_list;
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
	template<typename CellL = CellList<dim, St, FAST, shift<dim, St> > > CellL getCellList(St r_cut)
	{
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
	template<typename CellL = CellList_hilb<dim, St, FAST, shift<dim, St> > > CellL getCellList_hilb(St r_cut)
	{
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
	template<typename CellL = CellList<dim, St, FAST, shift<dim, St> > > void updateCellList(CellL & cell_list)
	{
		populate_cell_list(v_pos,cell_list,g_m,CL_NON_SYMMETRIC);

		cell_list.set_gm(g_m);
	}

	/*! \brief Update a cell list using the stored particles
	 *
	 * \tparam CellL CellList type to construct
	 *
	 * \param cell_list Cell list to update
	 *
	 */
	template<typename CellL = CellList<dim, St, FAST, shift<dim, St> > > void updateCellListSym(CellL & cell_list)
	{
		populate_cell_list(v_pos,cell_list,g_m,CL_SYMMETRIC);

		cell_list.set_gm(g_m);
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
	template<typename CellL = CellList<dim, St, FAST, shift<dim, St> > > CellL getCellList(St r_cut, const Ghost<dim, St> & enlarge)
	{
		CellL cell_list;

		// Division array
		size_t div[dim];

		// get the processor bounding box
		Box<dim, St> pbox = getDecomposition().getProcessorBounds();

		// Processor bounding box
		cl_param_calculate(pbox, div, r_cut, enlarge);

		cell_list.Initialize(pbox, div);

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
	template<typename CellL = CellList_hilb<dim, St, FAST, shift<dim, St> > > CellL getCellList_hilb(St r_cut, const Ghost<dim, St> & enlarge)
	{
		CellL cell_list;

		// Division array
		size_t div[dim];

		// get the processor bounding box
		Box<dim, St> pbox = getDecomposition().getProcessorBounds();

		// Processor bounding box
		cl_param_calculate(pbox,div, r_cut, enlarge);

		cell_list.Initialize(pbox, div, g_m);

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
		VerletList<dim,St,FAST,shift<dim,St>> ver;

		// Processor bounding box
		Box<dim, St> pbox = getDecomposition().getProcessorBounds();

		ver.InitializeSym(getDecomposition().getDomain(),pbox,getDecomposition().getGhost(),r_cut,v_pos,g_m);

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
		VerletList<dim,St,FAST,shift<dim,St>> ver;

		// get the processor bounding box
		Box<dim, St> bt = getDecomposition().getProcessorBounds();

		// Get the ghost
		Ghost<dim,St> g = getDecomposition().getGhost();
		g.magnify(1.013);

		// enlarge the box where the Verlet is defined
		bt.enlarge(g);

		ver.Initialize(bt,getDecomposition().getProcessorBounds(),r_cut,v_pos,g_m,VL_NON_SYMMETRIC);

		return ver;
	}

	/*! \brief for each particle get the verlet list
	 *
	 * \param r_cut cut-off radius
	 * \param ver Verlet to update
	 * \param r_cut cutoff radius
	 * \param opt option like VL_SYMMETRIC and VL_NON_SYMMETRIC
	 *
	 */
	void updateVerlet(VerletList<dim,St,FAST,shift<dim,St> > & ver, St r_cut, size_t opt = VL_NON_SYMMETRIC)
	{
		ver.update(getDecomposition().getDomain(),r_cut,v_pos,g_m, opt);
	}

	/*! \brief for each particle get the verlet list
	 *
	 * \param verlet output verlet list for each particle
	 * \param r_cut cut-off radius
	 *
	 * \deprecated
	 *
	 */
	void getVerletDeprecated(openfpm::vector<openfpm::vector<size_t>> & verlet, St r_cut)
	{
		// resize verlet to store the number of particles
		verlet.resize(size_local());

		// get the cell-list
		auto cl = getCellList(r_cut);

		// square of the cutting radius
		St r_cut2 = r_cut * r_cut;

		// iterate the particles
		auto it_p = this->getDomainIterator();
		while (it_p.isNext())
		{
			// key
			vect_dist_key_dx key = it_p.get();

			// Get the position of the particles
			Point<dim, St> p = this->getPos(key);

			// Clear the neighborhood of the particle
			verlet.get(key.getKey()).clear();

			// Get the neighborhood of the particle
			auto NN = cl.template getNNIterator<NO_CHECK>(cl.getCell(p));
			while (NN.isNext())
			{
				auto nnp = NN.get();

				// p != q
				if (nnp == key.getKey())
				{
					++NN;
					continue;
				}

				Point<dim, St> q = this->getPos(nnp);

				if (p.distance2(q) < r_cut2)
					verlet.get(key.getKey()).add(nnp);

				// Next particle
				++NN;
			}

			// next particle
			++it_p;
		}
	}

	/*! \brief Construct a cell list starting from the stored particles and reorder a vector according to the Hilberts curve
	 *
	 * \tparam CellL CellList type to construct
	 *
	 * \param m an order of a hilbert curve
	 *
	 *
	 *
	 */
	template<typename CellL=CellList<dim,St,FAST,shift<dim,St> > > void reorder (int32_t m)
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
	template<typename CellL=CellList<dim,St,FAST,shift<dim,St> > > void reorder(int32_t m, const Ghost<dim,St> & enlarge)
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

		// for each particle add the particle to the cell list

		auto it = getIterator();

		while (it.isNext())
		{
			auto key = it.get();

			cell_list.add(this->getPos(key),key.getKey());

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
		return vector_dist_iterator(g_m, v_pos.size());
	}

	/*! \brief Get an iterator that traverse the particles in the domain
	 *
	 * \return an iterator
	 *
	 */
	vector_dist_iterator getDomainIterator() const
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
		this->template map_<obp>(v_pos,v_prp,g_m);
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
		this->template ghost_get_<prp...>(v_pos,v_prp,g_m,opt);
	}

	/*! \brief It synchronize the properties and position of the ghost particles
	 *
	 * \tparam op which kind of operation to apply
	 * \tparam prp list of properties to get synchronize
	 *
	 *
	 */
	template<template<typename,typename> class op, int ... prp> inline void ghost_put()
	{
		this->template ghost_put_<op,prp...>(v_pos,v_prp,g_m);
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

	/*! \brief Add the computation cost on the decomposition comming from the particles
	 *
	 */
	inline void addComputationCosts()
	{
		CellDecomposer_sm<dim, St> cdsm;

		Decomposition & dec = getDecomposition();

		cdsm.setDimensions(dec.getDomain(), dec.getGrid().getSize(), 0);

		for (size_t i = 0; i < getDecomposition().getNSubSubDomains(); i++)
		{
			dec.setSubSubDomainComputationCost(i, 1);
		}

		auto it = getDomainIterator();

		while (it.isNext())
		{
			size_t v = cdsm.getCell(this->getPos(it.get()));

			dec.addComputationCost(v, 1);

			++it;
		}

	}

	/*! \brief Output particle position and properties
	 *
	 * \param out output
	 * \param opt VTK_WRITER or CSV_WRITER
	 *
	 * \return true if the file has been written without error
	 *
	 */
	inline bool write(std::string out, int opt = NO_GHOST | VTK_WRITER )
	{

		if ((opt & 0xFFFF0000) == CSV_WRITER)
		{
			// CSVWriter test
			CSVWriter<openfpm::vector<Point<dim,St>>, openfpm::vector<prop> > csv_writer;

			std::string output = std::to_string(out + "_" + std::to_string(v_cl.getProcessUnitID()) + std::to_string(".csv"));

			// Write the CSV
			return csv_writer.write(output,v_pos,v_prp);
		}
		else if ((opt & 0xFFFF0000) == VTK_WRITER)
		{
			// VTKWriter for a set of points
			VTKWriter<boost::mpl::pair<openfpm::vector<Point<dim,St>>, openfpm::vector<prop>>, VECTOR_POINTS> vtk_writer;
			vtk_writer.add(v_pos,v_prp,g_m);

			std::string output = std::to_string(out + "_" + std::to_string(v_cl.getProcessUnitID()) + std::to_string(".vtk"));

			// Write the VTK file
			return vtk_writer.write(output);
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
		if ((opt & 0xFFFF0000) == CSV_WRITER)
		{
			// CSVWriter test
			CSVWriter<openfpm::vector<Point<dim, St>>, openfpm::vector<prop> > csv_writer;

			std::string output = std::to_string(out + "_" + std::to_string(v_cl.getProcessUnitID()) + "_" + std::to_string(iteration) + std::to_string(".csv"));

			// Write the CSV
			return csv_writer.write(output, v_pos, v_prp);
		}
		else
		{
			// VTKWriter for a set of points
			VTKWriter<boost::mpl::pair<openfpm::vector<Point<dim,St>>, openfpm::vector<prop>>, VECTOR_POINTS> vtk_writer;
			vtk_writer.add(v_pos,v_prp,g_m);

			std::string output = std::to_string(out + "_" + std::to_string(v_cl.getProcessUnitID()) + "_" + std::to_string(iteration) + std::to_string(".vtk"));

			// Write the VTK file
			return vtk_writer.write(output);
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
};


#endif /* VECTOR_HPP_ */
