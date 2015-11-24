/*
 * CartDecomposition.hpp
 *
 *  Created on: Oct 07, 2015
 *      Author: Antonio Leo
 */

#ifndef CARTDECOMPOSITION_HPP
#define CARTDECOMPOSITION_HPP

#include "config.h"
#include <cmath>
#include "VCluster.hpp"
#include "Graph/CartesianGraphFactory.hpp"
#include "Decomposition.hpp"
#include "Vector/map_vector.hpp"
#include <vector>
#include <initializer_list>
#include "SubdomainGraphNodes.hpp"
#include "dec_optimizer.hpp"
#include "Space/Shape/Box.hpp"
#include "Space/Shape/Point.hpp"
#include "NN/CellList/CellDecomposer.hpp"
#include <unordered_map>
#include "NN/CellList/CellList.hpp"
#include "Space/Ghost.hpp"
#include "common.hpp"
#include "ie_loc_ghost.hpp"
#include "ie_ghost.hpp"
#include "nn_processor.hpp"
#include "GraphMLWriter.hpp"
#include "ParMetisDistribution.hpp"
#include "MetisDistribution.hpp"
#include "DLB.hpp"

#define CARTDEC_ERROR 2000lu

// Macro that decide what to do in case of error
#ifdef STOP_ON_ERROR
#define ACTION_ON_ERROR() exit(1);
#elif defined(THROW_ON_ERROR)
#define ACTION_ON_ERROR() throw CARTDEC_ERROR;
#else
#define ACTION_ON_ERROR()
#endif

/**
 * \brief This class decompose a space into subspaces
 *
 * \tparam dim is the dimensionality of the physical domain we are going to decompose.
 * \tparam T type of the space we decompose, Real, Integer, Complex ...
 * \tparam Memory Memory factory used to allocate memory
 * \tparam Domain Structure that contain the information of your physical domain
 *
 * Given an N-dimensional space, this class decompose the space into a Cartesian grid of small
 * sub-sub-domain. At each sub-sub-domain is assigned  an id that identify which processor is
 * going to take care of that part of space (in general the space assigned to a processor is
 * simply connected), a second step merge several sub-sub-domain with same id into bigger region
 *  sub-domain with the id. Each sub-domain has an extended space called ghost part
 *
 * Assuming that VCluster.getProcessUnitID(), equivalent to the MPI processor rank, return the processor local
 * processor id, we define
 *
 * * local processor: processor rank
 * * local sub-domain: sub-domain given to the local processor
 * * external ghost box: (or ghost box) are the boxes that compose the ghost space of the processor, or the
 *   boxes produced expanding every local sub-domain by the ghost extension and intersecting with the sub-domain
 *   of the other processors
 * * Near processors are the processors adjacent to the local processor, where with adjacent we mean all the processor
 *   that has a non-zero intersection with the ghost part of the local processor, or all the processors that
 *   produce non-zero external boxes with the local processor, or all the processor that should communicate
 *   in case of ghost data synchronization
 * * internal ghost box: is the part of ghost of the near processor that intersect the space of the
 *       processor, or the boxes produced expanding the sub-domain of the near processors with the local sub-domain
 * * Near processor sub-domain: is a sub-domain that live in the a near (or contiguous) processor
 * * Near processor list: the list of all the near processor of the local processor (each processor has a list
 *                        of the near processor)
 * * Local ghosts interal or external are all the ghosts that does not involve inter-processor communications
 *
 * \see calculateGhostBoxes() for a visualization of internal and external ghost boxes
 *
 * ### Create a Cartesian decomposition object on a Box space, distribute, calculate internal and external ghost boxes
 * \snippet CartDecomposition_unit_test.hpp Create CartDecomposition
 *
 */

template<unsigned int dim, typename T, typename Memory = HeapMemory,
		template<unsigned int, typename > class Domain = Box, typename Distribution = ParMetisDistribution<dim, T>>
class CartDecomposition: public ie_loc_ghost<dim, T>, public nn_prcs<dim, T>, public ie_ghost<dim, T>
{

public:

	//! Type of the domain we are going to decompose
	typedef T domain_type;

	//! It simplify to access the SpaceBox element
	typedef SpaceBox<dim, T> Box;

private:

	//! This is the key type to access  data_s, for example in the case of vector
	//! acc_key is size_t
	typedef typename openfpm::vector<SpaceBox<dim, T>, Memory, openfpm::vector_grow_policy_default,
			openfpm::vect_isel<SpaceBox<dim, T>>::value>::access_key acc_key;

	//! the set of all local sub-domain as vector
	openfpm::vector<SpaceBox<dim, T>> sub_domains;

	//! for each sub-domain, contain the list of the neighborhood processors
	openfpm::vector<openfpm::vector<long unsigned int> > box_nn_processor;

	//! Structure that contain for each sub-sub-domain box the processor id
	//! exist for efficient global communication
	openfpm::vector<size_t> fine_s;

	//! Structure that store the cartesian grid information
	grid_sm<dim, void> gr;

	//! Structure that decompose your structure into cell without creating them
	//! useful to convert positions to CellId or sub-domain id in this case
	CellDecomposer_sm<dim, T> cd;

	//! rectangular domain to decompose
	Domain<dim, T> domain;

	//! Box Spacing
	T spacing[dim];

	//! Runtime virtual cluster machine
	Vcluster & v_cl;

	//! Create ditribution
	Distribution dist;

	//! Cell-list that store the geometrical information of the local internal ghost boxes
	CellList<dim, T, FAST> lgeo_cell;

	/*! \brief Constructor, it decompose and distribute the sub-domains across the processors
	 *
	 * \param v_cl Virtual cluster, used internally for communications
	 *
	 */
	void CreateDecomposition(Vcluster & v_cl)
	{
#ifdef SE_CLASS1
		if (&v_cl == NULL)
		{
			std::cerr << __FILE__ << ":" << __LINE__ << " error VCluster instance is null, check that you ever initialized it \n";
			ACTION_ON_ERROR()
		}
#endif

		int p_id = v_cl.getProcessUnitID();

		// Calculate the total number of box and and the spacing
		// on each direction
		// Get the box containing the domain
		SpaceBox<dim, T> bs = domain.getBox();

		for (unsigned int i = 0; i < dim; i++)
		{
			// Calculate the spacing
			spacing[i] = (bs.getHigh(i) - bs.getLow(i)) / gr.size(i);
		}

		// fill the structure that store the processor id for each sub-domain
		fine_s.resize(gr.size());

		// Optimize the decomposition creating bigger spaces
		// And reducing Ghost over-stress
		dec_optimizer<dim, Graph_CSR<nm_v, nm_e>> d_o(dist.getGraph(), gr.getSize());

		// set of Boxes produced by the decomposition optimizer
		openfpm::vector<::Box<dim, size_t>> loc_box;

		// optimize the decomposition
		d_o.template optimize<nm_v::sub_id, nm_v::proc_id>(dist.getGraph(), p_id, loc_box, box_nn_processor);

		// Initialize ss_box and bbox
		if (loc_box.size() >= 0)
		{
			SpaceBox<dim, size_t> sub_dc = loc_box.get(0);
			SpaceBox<dim, T> sub_d(sub_dc);
			sub_d.mul(spacing);
			sub_d.expand(spacing);

			// Fixing sub-domains to cover all the domain

			// Fixing sub_d
			// if (loc_box) is a the boundary we have to ensure that the box span the full
			// domain (avoiding rounding off error)
			for (size_t i = 0; i < dim; i++)
			{
				if (sub_dc.getHigh(i) == cd.getGrid().size(i) - 1)
				{
					sub_d.setHigh(i, domain.getHigh(i));
				}
			}

			// add the sub-domain
			sub_domains.add(sub_d);

			ss_box = sub_d;
			ss_box -= ss_box.getP1();
			bbox = sub_d;
		}

		// convert into sub-domain
		for (size_t s = 1; s < loc_box.size(); s++)
		{
			SpaceBox<dim, size_t> sub_dc = loc_box.get(s);
			SpaceBox<dim, T> sub_d(sub_dc);

			// re-scale and add spacing (the end is the starting point of the next domain + spacing)
			sub_d.mul(spacing);
			sub_d.expand(spacing);

			// Fixing sub-domains to cover all the domain

			// Fixing sub_d
			// if (loc_box) is a the boundary we have to ensure that the box span the full
			// domain (avoiding rounding off error)
			for (size_t i = 0; i < dim; i++)
			{
				if (sub_dc.getHigh(i) == cd.getGrid().size(i) - 1)
				{
					sub_d.setHigh(i, domain.getHigh(i));
				}
			}

			// add the sub-domain
			sub_domains.add(sub_d);

			// Calculate the bound box
			bbox.enclose(sub_d);

			// Create the smallest box contained in all sub-domain
			ss_box.contained(sub_d);
		}

		nn_prcs<dim, T>::create(box_nn_processor, sub_domains);

		// fill fine_s structure
		// fine_s structure contain the processor id for each sub-sub-domain
		// with sub-sub-domain we mean the sub-domain decomposition before
		// running dec_optimizer (before merging sub-domains)
		auto it = dist.getGraph().getVertexIterator();

		while (it.isNext())
		{
			size_t key = it.get();

			// fill with the fine decomposition
			fine_s.get(key) = dist.getGraph().template vertex_p<nm_v::proc_id>(key);

			++it;
		}

		// Get the smallest sub-division on each direction
		::Box<dim, T> unit = getSmallestSubdivision();
		// Get the processor bounding Box
		::Box<dim, T> bound = getProcessorBounds();

		// calculate the sub-divisions
		size_t div[dim];
		for (size_t i = 0; i < dim; i++)
			div[i] = (size_t) ((bound.getHigh(i) - bound.getLow(i)) / unit.getHigh(i));

		// Create shift
		Point<dim, T> orig;

		// p1 point of the Processor bound box is the shift
		for (size_t i = 0; i < dim; i++)
			orig.get(i) = bound.getLow(i);

		// Initialize the geo_cell structure
		ie_ghost<dim, T>::Initialize_geo_cell(domain, div, orig);
		lgeo_cell.Initialize(domain, div, orig);

	}

	/*! \brief Calculate communication and migration costs
	 *
	 * \param gh_s ghost thickness
	 * \param ts how many timesteps have passed since last calculation, used to approximate the cost
	 */
	void computeCommunicationAndMigrationCosts(size_t ts)
	{

		size_t p_id = v_cl.getProcessUnitID();
		float migration;

		SpaceBox<dim, T> cellBox = cd.getCellBox();
		float b_s = cellBox.getHigh(0);
		float gh_s = ghost.getHigh(0);

		// compute the gh_area for 2 dim case
		float gh_v = (gh_s * b_s);

		// multiply for sub-sub-domain side for each domain
		for (int i = 2; i < dim; i++)
			gh_v *= b_s;

		size_t norm = (size_t) (1.0 / gh_v);

		migration = pow(b_s, dim);

		size_t prev = 0;

		for (size_t i = 0; i < dist.getNSubSubDomains(); i++)
		{
			dist.setMigrationCost(i, norm * migration * dist.getVertexWeight(i));

			for (size_t s = 0; s < dist.getNSubSubDomainNeighbors(i); s++)
			{
				dist.setCommunicationCost(prev + s, 1 * dist.getVertexWeight(i) * ts);
			}
			prev += dist.getNSubSubDomainNeighbors(i);
		}
	}

	// Save the ghost boundaries
	Ghost<dim, T> ghost;

	/*! \brief Create the subspaces that decompose your domain
	 *
	 */
	void CreateSubspaces()
	{
		// Create a grid where each point is a space
		grid_sm<dim, void> g(div);

		// create a grid_key_dx iterator
		grid_key_dx_iterator<dim> gk_it(g);

		// Divide the space into subspaces
		while (gk_it.isNext())
		{
			//! iterate through all subspaces
			grid_key_dx<dim> key = gk_it.get();

			//! Create a new subspace
			SpaceBox<dim, T> tmp;

			//! fill with the Margin of the box
			for (int i = 0; i < dim; i++)
			{
				tmp.setHigh(i, (key.get(i) + 1) * spacing[i]);
				tmp.setLow(i, key.get(i) * spacing[i]);
			}

			//! add the space box
			sub_domains.add(tmp);

			// add the iterator
			++gk_it;
		}
	}

	// Heap memory receiver
	HeapMemory hp_recv;

	// vector v_proc
	openfpm::vector<size_t> v_proc;

	// Receive counter
	size_t recv_cnt;

	// reference counter of the object in case is shared between object
	long int ref_cnt;

public:

	/*! \brief Cart decomposition constructor
	 *
	 * \param v_cl Virtual cluster, used internally to handle or pipeline communication
	 *
	 */
	CartDecomposition(Vcluster & v_cl) :
			nn_prcs<dim, T>(v_cl), v_cl(v_cl), dist(v_cl)
	{
		// Reset the box to zero
		bbox.zero();
	}

	//! Cart decomposition destructor
	~CartDecomposition()
	{
	}

	/*! \brief class to select the returned id by ghost_processorID
	 *
	 */
	class box_id
	{
	public:
		/*! \brief Return the box id
		 *
		 * \param p structure containing the id informations
		 * \param b_id box_id
		 *
		 * \return box id
		 *
		 */
		inline static size_t id(p_box<dim, T> & p, size_t b_id)
		{
			return b_id;
		}
	};

	/*! \brief class to select the returned id by ghost_processorID
	 *
	 */
	class processor_id
	{
	public:
		/*! \brief Return the processor id
		 *
		 * \param p structure containing the id informations
		 * \param b_id box_id
		 *
		 * \return processor id
		 *
		 */
		inline static size_t id(p_box<dim, T> & p, size_t b_id)
		{
			return p.proc;
		}
	};

	/*! \brief class to select the returned id by ghost_processorID
	 *
	 */
	class lc_processor_id
	{
	public:
		/*! \brief Return the near processor id
		 *
		 * \param p structure containing the id informations
		 * \param b_id box_id
		 *
		 * \return local processor id
		 *
		 */
		inline static size_t id(p_box<dim, T> & p, size_t b_id)
		{
			return p.lc_proc;
		}
	};

	/*! It calculate the internal ghost boxes
	 *
	 * Example: Processor 10 calculate
	 * B8_0 B9_0 B9_1 and B5_0
	 *
	 *
	 *
	 \verbatim

	 +----------------------------------------------------+
	 |                                                    |
	 |                 Processor 8                        |
	 |                 Sub-domain 0                       +-----------------------------------+
	 |                                                    |                                   |
	 |                                                    |                                   |
	 ++--------------+---+---------------------------+----+        Processor 9                |
	 |              |   |     B8_0                  |    |        Subdomain 0                |
	 |              +------------------------------------+                                   |
	 |              |   |                           |    |                                   |
	 |              |   |  XXXXXXXXXXXXX XX         |B9_0|                                   |
	 |              | B |  X Processor 10 X         |    |                                   |
	 | Processor 5  | 5 |  X Sub-domain 0 X         |    |                                   |
	 | Subdomain 0  | _ |  X              X         +----------------------------------------+
	 |              | 0 |  XXXXXXXXXXXXXXXX         |    |                                   |
	 |              |   |                           |    |                                   |
	 |              |   |                           |    |        Processor 9                |
	 |              |   |                           |B9_1|        Subdomain 1                |
	 |              |   |                           |    |                                   |
	 |              |   |                           |    |                                   |
	 |              |   |                           |    |                                   |
	 +--------------+---+---------------------------+----+                                   |
	 |                                   |
	 +-----------------------------------+

	 \endverbatim

	 and also
	 G8_0 G9_0 G9_1 G5_0 (External ghost boxes)

	 \verbatim

	 +----------------------------------------------------+
	 |                                                    |
	 |                 Processor 8                        |
	 |                 Sub-domain 0                       +-----------------------------------+
	 |           +---------------------------------------------+                              |
	 |           |         G8_0                           |    |                              |
	 ++--------------+------------------------------------+    |   Processor 9                |
	 |          |   |                                    |    |   Subdomain 0                |
	 |          |   |                                    |G9_0|                              |
	 |          |   |                                    |    |                              |
	 |          |   |      XXXXXXXXXXXXX XX              |    |                              |
	 |          |   |      X Processor 10 X              |    |                              |
	 | Processor|5  |      X Sub-domain 0 X              |    |                              |
	 | Subdomain|0  |      X              X              +-----------------------------------+
	 |          |   |      XXXXXXXXXXXXXXXX              |    |                              |
	 |          | G |                                    |    |                              |
	 |          | 5 |                                    |    |   Processor 9                |
	 |          | | |                                    |    |   Subdomain 1                |
	 |          | 0 |                                    |G9_1|                              |
	 |          |   |                                    |    |                              |
	 |          |   |                                    |    |                              |
	 +--------------+------------------------------------+    |                              |
	 |                                        |    |                              |
	 +----------------------------------------+----+------------------------------+


	 \endverbatim

	 *
	 *
	 *
	 * \param ghost margins for each dimensions (p1 negative part) (p2 positive part)
	 *
	 *
	 \verbatim
	 ^ p2[1]
	 |
	 |
	 +----+----+
	 |         |
	 |         |
	 p1[0]<-----+         +----> p2[0]
	 |         |
	 |         |
	 +----+----+
	 |
	 v  p1[1]

	 \endverbatim

	 *
	 *
	 */
	void calculateGhostBoxes()
	{
#ifdef DEBUG
		// the ghost margins are assumed to be smaller
		// than one sub-domain

		for (size_t i = 0; i < dim; i++)
		{
			if (ghost.template getLow(i) >= domain.template getHigh(i) / gr.size(i)
					|| ghost.template getHigh(i) >= domain.template getHigh(i) / gr.size(i))
			{
				std::cerr << "Error " << __FILE__ << ":" << __LINE__ << " : Ghost are bigger than one domain" << "\n";
			}
		}
#endif

		// Intersect all the local sub-domains with the sub-domains of the contiguous processors

		// create the internal structures that store ghost information
		ie_ghost<dim, T>::create_box_nn_processor_ext(v_cl, ghost, sub_domains, box_nn_processor, *this);
		ie_ghost<dim, T>::create_box_nn_processor_int(v_cl, ghost, sub_domains, box_nn_processor, *this);

		// ebox must come after ibox (in this case)

		ie_loc_ghost<dim, T>::create_loc_ghost_ibox(ghost, sub_domains);
		ie_loc_ghost<dim, T>::create_loc_ghost_ebox(ghost, sub_domains);

		// get the smallest sub-domain dimension on each direction
		for (size_t i = 0; i < dim; i++)
		{
			if (ghost.template getLow(i) >= ss_box.getHigh(i)
					|| ghost.template getHigh(i) >= domain.template getHigh(i) / gr.size(i))
			{
				std::cerr << "Error " << __FILE__ << ":" << __LINE__ << " : Ghost are bigger than one domain" << "\n";
			}
		}
	}

	/*! \brief The default grid size
	 *
	 *  The default grid is always an isotropic grid that adapt with the number of processors,
	 *  it define in how many cell it will be divided the space for a particular required minimum
	 *  number of sub-domain
	 *
	 */
	static size_t getDefaultGrid(size_t n_sub)
	{
		// Calculate the number of sub-sub-domain on
		// each dimension
		return openfpm::math::round_big_2(pow(n_sub, 1.0 / dim));
	}

	/*! \brief Given a point return in which processor the particle should go
	 *
	 * \return processorID
	 *
	 */
	template<typename Mem> size_t inline processorID(encapc<1, Point<dim, T>, Mem> p)
	{
		return fine_s.get(cd.getCell(p));
	}

	// Smallest subdivision on each direction
	::Box<dim, T> ss_box;

	/*! \brief Get the smallest subdivision of the domain on each direction
	 *
	 * \return a box p1 is set to zero
	 *
	 */
	const ::Box<dim, T> & getSmallestSubdivision()
	{
		return ss_box;
	}

	/*! \brief Given a point return in which processor the particle should go
	 *
	 * \return processorID
	 *
	 */

	size_t inline processorID(const T (&p)[dim]) const
	{
		return fine_s.get(cd.getCell(p));
	}

	/*! \brief Set the parameter of the decomposition
	 *
	 * \param div_ storing into how many domain to decompose on each dimension
	 * \param domain_ domain to decompose
	 *
	 */
	void setParameters(const size_t (&div_)[dim], Domain<dim, T> domain_, Ghost<dim, T> ghost = Ghost<dim, T>())
	{
		// set the ghost
		this->ghost = ghost;

		// Set the decomposition parameters
		gr.setDimensions(div_);
		domain = domain_;
		cd.setDimensions(domain, div_, 0);

		// init distribution
		dist.init(gr, domain);

	}

	/*! \brief Start decomposition
	 *
	 */
	void decompose()
	{
		computeCommunicationAndMigrationCosts(1);

		dist.decompose();

		CreateDecomposition(v_cl);
	}

	/*! \brief Refine the decomposition, available only for ParMetis distribution, for Metis it is a null call
	 *
	 */
	void rebalance()
	{
		computeCommunicationAndMigrationCosts(1);

		dist.refine();
	}

	/*! \brief Refine the decomposition, available only for ParMetis distribution, for Metis it is a null call
	 *
	 */
	void rebalance(DLB & dlb)
	{
		if (dlb.rebalanceNeeded())
		{
			computeCommunicationAndMigrationCosts(dlb.getNTimeStepSinceDLB());
			dist.refine();
		}
	}

	/*! \brief function that return the position of the cell in the space
	 *
	 * \param id vertex id
	 * \param pos vector that will contain x, y, z
	 *
	 */
	inline void getSubSubDomainPosition(size_t id, openfpm::vector<real_t> &pos)
	{
		dist.getVertexPosition(id, pos);
	}

	/*! \brief function that set the weight of the vertex
	 *
	 * \param id vertex id
	 *
	 */
	inline void setSubSubDomainComputationCost(size_t id, size_t weight)
	{
		dist.setVertexWeight(id, weight);
	}

	/*! \brief Get the number of local sub-domains
	 *
	 * \return the number of sub-domains
	 *
	 */
	size_t getNLocalHyperCube()
	{
		return sub_domains.size();
	}

	/*! \brief Get the local sub-domain
	 *
	 * \param i (each local processor can have more than one sub-domain)
	 * \return the sub-domain
	 *
	 */
	SpaceBox<dim, T> getLocalHyperCube(size_t lc)
	{
		// Create a space box
		SpaceBox<dim, T> sp;

		// fill the space box

		for (size_t k = 0; k < dim; k++)
		{
			// create the SpaceBox Low and High
			sp.setLow(k, sub_domains.template get<Box::p1>(lc)[k]);
			sp.setHigh(k, sub_domains.template get<Box::p2>(lc)[k]);
		}

		return sp;
	}

	/*! \brief Get the local sub-domain with ghost extension
	 *
	 * \param i (each local processor can have more than one sub-domain)
	 * \return the sub-domain
	 *
	 */
	SpaceBox<dim, T> getSubDomainWithGhost(size_t lc)
	{
		// Create a space box
		SpaceBox<dim, T> sp = sub_domains.get(lc);

		// enlarge with ghost
		sp.enlarge(ghost);

		return sp;
	}

	/*! \brief Return the structure that store the physical domain
	 *
	 * \return The physical domain
	 *
	 */
	Domain<dim, T> & getDomain()
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
	template<typename Mem> bool isLocal(const encapc<1, Point<dim, T>, Mem> p) const
	{
		return processorID<Mem>(p) == v_cl.getProcessUnitID();
	}

	/*! \brief Check if the particle is local
	 *
	 * \param p object position
	 *
	 * \return true if it is local
	 *
	 */
	bool isLocal(const T (&pos)[dim]) const
	{
		return processorID(pos) == v_cl.getProcessUnitID();
	}

	::Box<dim, T> bbox;

	/*! \brief Return the bounding box containing union of all the sub-domains for the local processor
	 *
	 * \return The bounding box
	 *
	 */
	::Box<dim, T> & getProcessorBounds()
	{
		return bbox;
	}

	////////////// Functions to get decomposition information ///////////////

	/*! \brief Write the decomposition as VTK file
	 *
	 * The function generate several files
	 *
	 * * subdomains_X.vtk domain for the local processor (X) as union of sub-domain
	 * * subdomains_adjacent_X.vtk sub-domains adjacent to the local processor (X)
	 * * internal_ghost_X.vtk Internal ghost boxes for the local processor (X)
	 * * external_ghost_X.vtk External ghost boxes for the local processor (X)
	 * * local_internal_ghost_X.vtk internal local ghost boxes for the local processor (X)
	 * * local_external_ghost_X.vtk external local ghost boxes for the local processor (X)
	 *
	 * where X is the local processor rank
	 *
	 * \param output directory where to write the files
	 *
	 */
	bool write(std::string output) const
	{
		//! subdomains_X.vtk domain for the local processor (X) as union of sub-domain
		VTKWriter<openfpm::vector<::SpaceBox<dim, T>>, VECTOR_BOX> vtk_box1;
		vtk_box1.add(sub_domains);
		vtk_box1.write(
				output + std::string("subdomains_") + std::to_string(v_cl.getProcessUnitID()) + std::string(".vtk"));

		nn_prcs<dim, T>::write(output);
		ie_ghost<dim, T>::write(output, v_cl.getProcessUnitID());
		ie_loc_ghost<dim, T>::write(output, v_cl.getProcessUnitID());

		return true;
	}

	/*! \brief function to check the consistency of the information of the decomposition
	 *
	 * \return false if is inconsistent
	 *
	 */
	bool check_consistency()
	{
		if (ie_loc_ghost<dim, T>::check_consistency(getNLocalHyperCube()) == false)
			return false;

		return true;
	}

	/*! \brief Print subdomains, external and internal ghost boxes
	 *
	 */
	void debugPrint()
	{
		std::cout << "Subdomains\n";
		for (size_t p = 0; p < sub_domains.size(); p++)
		{
			std::cout << ::SpaceBox<dim, T>(sub_domains.get(p)).toString() << "\n";
		}

		std::cout << "External ghost box\n";

		for (size_t p = 0; p<nn_prcs < dim, T>::getNNProcessors(); p++)
		{
			for (size_t i = 0; i<ie_ghost < dim, T>::getProcessorNEGhost(p); i++)
			{
				std::cout << ie_ghost<dim, T>::getProcessorEGhostBox(p, i).toString() << "   prc="
						<< nn_prcs<dim, T>::IDtoProc(p) << "   id=" << ie_ghost<dim, T>::getProcessorEGhostId(p, i)
						<< "\n";
			}
		}

		std::cout << "Internal ghost box\n";

		for (size_t p = 0; p<nn_prcs < dim, T>::getNNProcessors(); p++)
		{
			for (size_t i = 0; i<ie_ghost < dim, T>::getProcessorNIGhost(p); i++)
			{
				std::cout << ie_ghost<dim, T>::getProcessorIGhostBox(p, i).toString() << "   prc="
						<< nn_prcs<dim, T>::IDtoProc(p) << "   id=" << ie_ghost<dim, T>::getProcessorIGhostId(p, i)
						<< "\n";
			}
		}
	}

	/*! \brief Print current graph and save it to file with name test_graph_[id]
	 *
	 * \param id to attach to the filename
	 *
	 */
	void printCurrentDecomposition(int id)
	{
		if (v_cl.getProcessUnitID() == 0)
		{
			dist.printCurrentDecomposition(id);
		}
	}

	//! Increment the reference counter
	void incRef()
	{
		ref_cnt++;
	}

	//! Decrement the reference counter
	void decRef()
	{
		ref_cnt--;
	}

	//! Return the reference counter
	long int ref()
	{
		return ref_cnt;
	}
};

#endif
