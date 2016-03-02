/*
 * CartDecomposition.hpp
 *
 *  Created on: Oct 07, 2015
 *      Author: Pietro Incardona, Antonio Leo
 */

#ifndef CARTDECOMPOSITION_HPP
#define CARTDECOMPOSITION_HPP

#include "config.h"
#include <cmath>
#include "VCluster.hpp"
#include "Graph/CartesianGraphFactory.hpp"
#include "Graph/DistCartesianGraphFactory.hpp"
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
#include "GraphMLWriter/GraphMLWriter.hpp"
#include "Distribution/ParMetisDistribution.hpp"
#include "Distribution/DistParMetisDistribution.hpp"
#include "Distribution/MetisDistribution.hpp"
#include "DLB/DLB.hpp"
#include "util/se_util.hpp"
#include "util/mathutil.hpp"

#define CARTDEC_ERROR 2000lu

/**
 * \brief This class decompose a space into subspaces
 *
 * \tparam dim is the dimensionality of the physical domain we are going to decompose.
 * \tparam T type of the space we decompose, Real, Integer, Complex ...
 * \tparam Memory Memory factory used to allocate memory
 * \tparam Distribution type of distribution, can be ParMetisDistribution or MetisDistribution
 *
 * Given an N-dimensional space, this class decompose the space into a Cartesian grid of small
 * sub-sub-domain. To each sub-sub-domain is assigned an id that identify at which processor is
 * assigned (in general the union of all the sub-sub-domain assigned to a processor is
 * simply connected space), a second step merge several sub-sub-domain with same id into bigger region
 *  sub-domain. Each sub-domain has an extended space called ghost part
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

template<unsigned int dim, typename T, typename Memory = HeapMemory, typename Distribution = ParMetisDistribution<dim, T>>
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
	typedef typename openfpm::vector<SpaceBox<dim, T>, Memory, openfpm::vector_grow_policy_default, openfpm::vect_isel<SpaceBox<dim, T>>::value>::access_key acc_key;

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
	::Box<dim,T> domain;

	//! Box Spacing
	T spacing[dim];

	//! Runtime virtual cluster machine
	Vcluster & v_cl;

	//! Create distribution
	Distribution dist;

	// Smallest subdivision on each direction
	::Box<dim,T> ss_box;

	::Box<dim,T> bbox;

	// Heap memory receiver
	HeapMemory hp_recv;

	// vector v_proc
	openfpm::vector<size_t> v_proc;

	// reference counter of the object in case is shared between object
	long int ref_cnt;

	// ghost info
	Ghost<dim,T> ghost;

	// Boundary condition info
	size_t bc[dim];

	/*! \brief Constructor, it decompose and distribute the sub-domains across the processors
	 *
	 * \param v_cl Virtual cluster, used internally for communications
	 *
	 */
	void createSubdomains(Vcluster & v_cl, const size_t (& bc)[dim])
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
		d_o.template optimize<nm_v::sub_id, nm_v::proc_id>(dist.getGraph(), p_id, loc_box, box_nn_processor,bc);

		// Initialize ss_box and bbox
		if (loc_box.size() >= 0)
		{
			SpaceBox<dim, size_t> sub_dc = loc_box.get(0);
			SpaceBox<dim, T> sub_d(sub_dc);
			sub_d.mul(spacing);
			sub_d.expand(spacing);
			sub_d += domain.getP1();

			// we add the

			// Fixing sub-domains to cover all the domain

			// Fixing sub_d
			// if (loc_box) is at the boundary we have to ensure that the box span the full
			// domain (avoiding rounding off error)
			for (size_t i = 0; i < dim; i++)
			{
				if (sub_dc.getHigh(i) == cd.getGrid().size(i) - 1)
					sub_d.setHigh(i, domain.getHigh(i));

				if (sub_dc.getLow(i) == 0)
					sub_d.setLow(i,domain.getLow(i));
			}

			// add the sub-domain
			sub_domains.add(sub_d);

			ss_box = sub_d;
			ss_box -= ss_box.getP1();
			bbox = sub_d;
		}

/*		if (loc_box.size())
		bbox.zero();
		ss_box = domain;*/

		// convert into sub-domain
		for (size_t s = 1; s < loc_box.size(); s++)
		{
			SpaceBox<dim, size_t> sub_dc = loc_box.get(s);
			SpaceBox<dim, T> sub_d(sub_dc);

			// re-scale and add spacing (the end is the starting point of the next domain + spacing)
			sub_d.mul(spacing);
			sub_d.expand(spacing);
			sub_d += domain.getP1();

			// Fixing sub-domains to cover all the domain

			// Fixing sub_d
			// if (loc_box) is a the boundary we have to ensure that the box span the full
			// domain (avoiding rounding off error)
			for (size_t i = 0; i < dim; i++)
			{
				if (sub_dc.getHigh(i) == cd.getGrid().size(i) - 1)
					sub_d.setHigh(i, domain.getHigh(i));

				if (sub_dc.getLow(i) == 0)
					sub_d.setLow(i,domain.getLow(i));
			}

			// add the sub-domain
			sub_domains.add(sub_d);

			// Calculate the bound box
			bbox.enclose(sub_d);

			// Create the smallest box contained in all sub-domain
			ss_box.contained(sub_d);
		}

		nn_prcs<dim,T>::create(box_nn_processor, sub_domains);
		nn_prcs<dim,T>::refine_ss_box(ss_box);
		nn_prcs<dim,T>::applyBC(domain,ghost,bc);

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

		Initialize_geo_cell_lists();
	}

	/*! \brief Initialize geo_cell lists
	 *
	 *
	 *
	 */
	void Initialize_geo_cell_lists()
	{
		// Get the smallest sub-division on each direction
		::Box<dim, T> unit = getSmallestSubdivision();
		// Get the processor bounding Box
		::Box<dim,T> bound = getProcessorBounds();
		// Not necessary, but I prefer
		bound.enlarge(ghost);

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
		ie_ghost<dim,T>::Initialize_geo_cell(bound,div,orig);

		// Initialize shift vectors
		ie_ghost<dim,T>::generateShiftVectors(domain);
	}

	/*! \brief Calculate communication and migration costs
	 *
	 * \param ts how many timesteps have passed since last calculation, used to approximate the cost
	 */
	void computeCommunicationAndMigrationCosts(size_t ts)
	{
		float migration = 0;

		SpaceBox<dim, T> cellBox = cd.getCellBox();
		float b_s = cellBox.getHigh(0);
		float gh_s = ghost.getHigh(0);

		// compute the gh_area for 2 dim case
		float gh_v = (gh_s * b_s);

		// multiply for sub-sub-domain side for each domain
		for (size_t i = 2; i < dim; i++)
			gh_v *= b_s;

		size_t norm = (size_t) (1.0 / gh_v);

		migration = pow(b_s, dim);

		size_t prev = 0;

		for (size_t i = 0; i < dist.getNSubSubDomains(); i++)
		{
			dist.setMigrationCost(i, norm * migration * dist.getVertexWeight(i));

			for (size_t s = 0; s < dist.getNSubSubDomainNeighbors(i); s++)
			{
				dist.setCommunicationCost(i, s, 1 * dist.getVertexWeight(i) * ts);
			}
			prev += dist.getNSubSubDomainNeighbors(i);
		}
	}

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


	/*! \brief It copy the sub-domains into another CartesianDecomposition object extending them
	 *
	 * \see duplicate (in case of extended domain)
	 *
	 * \param cart Cartesian decomposition object
	 * \param box Extended domain
	 *
	 */
	void extend_subdomains(CartDecomposition<dim,T> & cart, const ::Box<dim,T> & ext_dom) const
	{
		// Box
		typedef ::Box<dim,T> b;

		cart.bbox = ext_dom;
		cart.ss_box = ext_dom;

		for (size_t i = 0 ; i < sub_domains.size() ; i++)
		{
			::Box<dim,T> box;

			// Calculate the extended box
			for (size_t j = 0 ; j < dim ; j++)
			{
				if (sub_domains.template get<b::p1>(i)[j] == domain.getLow(j))
					box.setLow(j,ext_dom.getLow(j));
				else
					box.setLow(j,sub_domains.template get<b::p1>(i)[j]);

				if (sub_domains.template get<b::p2>(i)[j] == domain.getHigh(j))
					box.setHigh(j,ext_dom.getHigh(j));
				else
					box.setHigh(j,sub_domains.template get<b::p2>(i)[j]);
			}

			// add the subdomain
			cart.sub_domains.add(box);

			// Calculate the bound box
			cart.bbox.enclose(box);

			// Create the smallest box contained in all sub-domain
			cart.ss_box.contained(box);
		}
	}

	/*! \brief Extend the fines for the new Cartesian decomposition
	 *
	 * \param new_fines extended fine_s
	 * \param old_fines old fine_s
	 *
	 */
	void extend_fines(CartDecomposition<dim,T> & cart) const
	{
		// Extension, first we calculate the extensions of the new domain compared
		// to the old one in cell units (each cell unit is a sub-sub-domain)
		::Box<dim,size_t> ext;
		// Extension of the new fines structure
		::Box<dim,size_t> n_fines_ext;
		// Extension of the old fines structure
		::Box<dim,size_t> o_fines_ext;

		size_t sz_new[dim];
		size_t sz_old[dim];

		for (size_t i = 0; i < dim ; i++)
		{
			size_t p1 = (domain.getLow(i) - this->domain.getLow(i)) / cd.getCellBox().getHigh(i) + 1;
			size_t p2 = (domain.getLow(i) - this->domain.getLow(i)) / cd.getCellBox().getHigh(i) + 1;

			ext.setLow(i,p1);
			ext.setHigh(i,p2);
			sz_new[i] = p1+p2+cd.getGrid().size(i);
			sz_old[i] = cd.getGrid().size(i);
		}

		grid_sm<dim,void> info_new(sz_new);
		grid_sm<dim,void> info_old(sz_old);

		// resize the new fines
		cart.fine_s.resize(info_new.size());

		// we create an iterator that iterate across the full new fines
		grid_key_dx_iterator<dim> fines_t(info_new);

		while (fines_t.isNext())
		{
			auto key = fines_t.get();

			// new_fines is bigger than old_fines structure
			// out of bound key must be adjusted
			// The adjustment produce a natural extension
			// a representation can be seen in the figure of
			// CartDecomposition duplicate function with extended domains

			grid_key_dx<dim> key_old;
			for (size_t i = 0 ; i < dim ; i++)
			{
				key_old.set_d(i,(long int)key.get(i) - ext.getLow(i));
				if (key_old.get(i) < 0)
					key_old.set_d(i,0);
				else if(key_old.get(i) >= (long int)info_old.size(i) )
					key_old.set_d(i,info_old.size(i)-1);
			}

			cart.fine_s.get(info_new.LinId(key)) = fine_s.get(info_old.LinId(key_old));

			++fines_t;
		}

		cart.gr.setDimensions(sz_new);

		// the new extended CellDecomposer must be consistent with the old cellDecomposer.
		cart.cd.setDimensions(cd,ext);
	}

public:

	static constexpr int dims = dim;

	typedef T stype;

	//! Increment the reference counter
	void incRef()
	{ref_cnt++;}

	//! Decrement the reference counter
	void decRef()
	{ref_cnt--;}

	//! Return the reference counter
	long int ref()
	{
		return ref_cnt;
	}

	/*! \brief Cartesian decomposition constructor
	 *
	 * \param v_cl Virtual cluster, used internally to handle or pipeline communication
	 *
	 */
	CartDecomposition(Vcluster & v_cl) :
			nn_prcs<dim, T>(v_cl), v_cl(v_cl), dist(v_cl),ref_cnt(0)
	{
		// Reset the box to zero
		bbox.zero();
	}

	/*! \brief Cartesian decomposition copy constructor
	 *
     * \param cart object to copy
	 *
	 */
	CartDecomposition(const CartDecomposition<dim,T,Memory> & cart)
	:nn_prcs<dim,T>(cart.v_cl),v_cl(cart.v_cl),dist(v_cl),ref_cnt(0)
	{
		this->operator=(cart);
	}

	/*! \brief Cartesian decomposition copy constructor
	 *
     * \param cart object to copy
	 *
	 */
	CartDecomposition(CartDecomposition<dim,T,Memory> && cart)
	:nn_prcs<dim,T>(cart.v_cl),v_cl(cart.v_cl),dist(v_cl),ref_cnt(0)
	{
		this->operator=(cart);
	}

	//! Cartesian decomposition destructor
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

	/*! \brief class to select the returned id by ghost_processorID
	 *
	 */
	class shift_id
	{
	public:
		/*! \brief Return the shift id
		 *
		 * \param p structure containing the id informations
		 * \param b_id box_id
		 *
		 * \return shift_id id
		 *
		 */
		inline static size_t id(p_box<dim,T> & p, size_t b_id)
		{
			return p.shift_id;
		}
	};

	/*! \brief Apply boundary condition to the point
	 *
	 * \param p Point to apply the boundary condition
	 *
	 */
	void applyPointBC(float (& pt)[dim]) const
	{
		for (size_t i = 0 ; i < dim ; i++)
		{
			if (bc[i] == PERIODIC)
				pt[i] = openfpm::math::periodic_l(pt[i],domain.getHigh(i),domain.getLow(i));
		}
	}

	/*! \brief Apply boundary condition to the point
	 *
	 * \param p Point to apply the boundary condition
	 *
	 */
	void applyPointBC(Point<dim,T> & pt) const
	{
		for (size_t i = 0 ; i < dim ; i++)
		{
			if (bc[i] == PERIODIC)
				pt.get(i) = openfpm::math::periodic_l(pt.get(i),domain.getHigh(i),domain.getLow(i));
		}
	}

	/*! \brief Apply boundary condition to the point
	 *
	 * \param encapsulated object
	 *
	 */
	template<typename Mem> void applyPointBC(encapc<1,Point<dim,T>,Mem> && pt) const
	{
		for (size_t i = 0 ; i < dim ; i++)
		{
			if (bc[i] == PERIODIC)
				pt.template get<0>()[i] = openfpm::math::periodic_l(pt.template get<0>()[i],domain.getHigh(i),domain.getLow(i));
		}
	}

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
|                 Sub+domain 0                       +-----------------------------------+
|                                                    |                                   |
|                                                    |                                   |
++--------------+---+---------------------------+----+        Processor 9                |
 |              |   |     B8_0                  |    |        Subdomain 0                |
 |              +------------------------------------+                                   |
 |              |   |                           |    |                                   |
 |              |   |                           |B9_0|                                   |
 |              | B |    Local processor        |    |                                   |
 | Processor 5  | 5 |    Subdomain 0            |    |                                   |
 | Subdomain 0  | _ |                           +----------------------------------------+
 |              | 0 |                           |    |                                   |
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

      +----------------------------------------------------+
      |                 Processor 8                        |
      |                 Subdomain 0                        +-----------------------------------+
      |                                                    |                                   |
      |           +---------------------------------------------+                              |
      |           |         G8_0                           |    |                              |
+-----+---------------+------------------------------------+    |   Processor 9                |
|                 |   |                                    |    |   Subdomain 0                |
|                 |   |                                    |G9_0|                              |
|                 |   |                                    |    |                              |
|                 |   |                                    |    |                              |
|                 |   |        Local processor             |    |                              |
|  Processor 5    |   |        Sub+domain 0                |    |                              |
|  Subdomain 0    |   |                                    +-----------------------------------+
|                 |   |                                    |    |                              |
|                 | G |                                    |    |                              |
|                 | 5 |                                    |    |   Processor 9                |
|                 | | |                                    |    |   Subdomain 1                |
|                 | 0 |                                    |G9_1|                              |
|                 |   |                                    |    |                              |
|                 |   |                                    |    |                              |
+---------------------+------------------------------------+    |                              |
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
			if (fabs(ghost.template getLow(i)) >= ss_box.getHigh(i) || ghost.template getHigh(i) >= ss_box.getHigh(i))
			{
				std::cerr << "Error " << __FILE__ << ":" << __LINE__  << " : Ghost are bigger than one sub-domain" << "\n";
			}
		}
#endif

		// Intersect all the local sub-domains with the sub-domains of the contiguous processors

		// create the internal structures that store ghost information
		ie_ghost<dim, T>::create_box_nn_processor_ext(v_cl, ghost, sub_domains, box_nn_processor, *this);
		ie_ghost<dim, T>::create_box_nn_processor_int(v_cl, ghost, sub_domains, box_nn_processor, *this);

		ie_loc_ghost<dim,T>::create(sub_domains,domain,ghost,bc);

		// get the smallest sub-domain dimension on each direction
		for (size_t i = 0; i < dim; i++)
		{
			if (fabs(ghost.template getLow(i)) >= ss_box.getHigh(i) || ghost.template getHigh(i) >= ss_box.getHigh(i))
			{
				std::cerr << "Error " << __FILE__ << ":" << __LINE__  << " : Ghost are bigger than one sub-domain" << "\n";
			}
		}
	}

	/*! \brief It create another object that contain the same decomposition information but with different ghost boxes
	 *
	 * \param g ghost
	 *
	 * \return a duplicated decomposition with different ghost boxes
	 *
	 */
	CartDecomposition<dim,T,Memory> duplicate(const Ghost<dim,T> & g) const
	{
		CartDecomposition<dim,T,Memory> cart(v_cl);

		cart.box_nn_processor = box_nn_processor;
		cart.sub_domains = sub_domains;
		cart.fine_s = fine_s;

		cart.gr = gr;
		cart.cd = cd;
		cart.domain = domain;
		std::copy(spacing,spacing+3,cart.spacing);

		//! Runtime virtual cluster
		cart.v_cl = v_cl;

		cart.bbox = bbox;
		cart.ss_box = ss_box;
		cart.ghost = g;

		cart.dist = dist;

		for (size_t i = 0 ; i < dim ; i++)
			cart.bc[i] = bc[i];

		(static_cast<nn_prcs<dim,T> &>(cart)).create(box_nn_processor, sub_domains);
		(static_cast<nn_prcs<dim,T> &>(cart)).applyBC(domain,ghost,bc);

		cart.Initialize_geo_cell_lists();
		cart.calculateGhostBoxes();

		return cart;
	}

	/*! \brief It create another object that contain the same decomposition information but with different ghost boxes and an extended domain
	 *
	 * The domain extension is produced extending the boxes at the border like in figure
	 *
	 * \verbatim
	 *
+--------------^--------^----------^----------+
|              |        |          |          |
|        A     |    E   |     F    |    N     |
|    +-----------------------------------+---->
|    |         |        |          |     |    |
|  A |   A     |        |     F    |     |    |
|    |         |        |          |     |    |
|    |         |    E   +----------+  N  |  N |
<--------------+        |          |     |    |
|    |         |        |          |     |    |
|    |         |        |     G    |     |    |
|    |         |        |          +---------->
|  B |   B     |        +----------+     |    |
|    |         +--------+          |  M  |  M |
|    |         |        |     H    |     |    |
|    |         |        +-----+----+---------->
<--------------+    D   |     |          |    |
|    |         |        |  I  |     L    |  L |
|  C |   C     |        |     |          |    |
|    |         |        |     |          |    |
|    +-----------------------------------+    |
|              |        |     |               |
|        C     |    D   |  I  |     L         |
+--------------v--------v-----v---------------+

	 *
	 * \endverbatim
	 *
	 * \param g ghost
	 * \param domain extended domain (MUST be extended)
	 *
	 * \return a duplicated decomposition with different ghost boxes and an extended domain
	 *
	 */
	CartDecomposition<dim,T,Memory> duplicate(const Ghost<dim,T> & g, const ::Box<dim,T> & ext_domain) const
	{
		CartDecomposition<dim,T,Memory> cart(v_cl);

		cart.box_nn_processor = box_nn_processor;

		// Calculate new sub-domains for extended domain
		extend_subdomains(cart,ext_domain);

		// Calculate fine_s structure for the extended domain
		// update the cell decomposer and gr
		extend_fines(cart);

		// Get the old sub-sub-domain grid extension

		cart.domain = ext_domain;

		// spacing does not change
		std::copy(spacing,spacing+3,cart.spacing);

		//! Runtime virtual cluster
		cart.v_cl = v_cl;

		cart.ghost = g;
		cart.dist = dist;

		for (size_t i = 0 ; i < dim ; i++)
			cart.bc[i] = bc[i];

		(static_cast<nn_prcs<dim,T> &>(cart)).create(cart.box_nn_processor, cart.sub_domains);
		(static_cast<nn_prcs<dim,T> &>(cart)).applyBC(ext_domain,ghost,bc);

		cart.Initialize_geo_cell_lists();
		cart.calculateGhostBoxes();

		return cart;
	}

	/*! \brief It create another object that contain the same information and act in the same way
	 *
	 * \return a duplicated decomposition
	 *
	 */
	CartDecomposition<dim,T,Memory> duplicate() const
	{
		CartDecomposition<dim,T,Memory> cart(v_cl);

		(static_cast<ie_loc_ghost<dim,T>*>(&cart))->operator=(static_cast<ie_loc_ghost<dim,T>>(*this));
		(static_cast<nn_prcs<dim,T>*>(&cart))->operator=(static_cast<nn_prcs<dim,T>>(*this));
		(static_cast<ie_ghost<dim,T>*>(&cart))->operator=(static_cast<ie_ghost<dim,T>>(*this));

		cart.sub_domains = sub_domains;
		cart.box_nn_processor = box_nn_processor;
		cart.fine_s = fine_s;
		cart.gr = gr;
		cart.cd = cd;
		cart.domain = domain;
		std::copy(spacing,spacing+3,cart.spacing);

		//! Runtime virtual cluster
		cart.v_cl = v_cl;

		cart.ghost = ghost;

		cart.bbox = bbox;
		cart.ss_box = ss_box;

		for (size_t i = 0 ; i < dim ; i++)
			cart.bc[i] = this->bc[i];

		return cart;
	}

	/*! \brief Copy the element
	 *
	 * \param cart element to copy
	 *
	 */
	CartDecomposition<dim,T,Memory> & operator=(const CartDecomposition & cart)
	{
		static_cast<ie_loc_ghost<dim,T>*>(this)->operator=(static_cast<ie_loc_ghost<dim,T>>(cart));
		static_cast<nn_prcs<dim,T>*>(this)->operator=(static_cast<nn_prcs<dim,T>>(cart));
		static_cast<ie_ghost<dim,T>*>(this)->operator=(static_cast<ie_ghost<dim,T>>(cart));

		sub_domains = cart.sub_domains;
		box_nn_processor = cart.box_nn_processor;
		fine_s = cart.fine_s;
		gr = cart.gr;
		cd = cart.cd;
		domain = cart.domain;
		std::copy(cart.spacing,cart.spacing+3,spacing);

		//! Runtime virtual cluster
		v_cl = cart.v_cl;

		ghost = cart.ghost;

		bbox = cart.bbox;
		ss_box = cart.ss_box;

		for (size_t i = 0 ; i < dim ; i++)
			bc[i] = cart.bc[i];

		return *this;
	}

	/*! \brief Copy the element, move semantic
	 *
	 * \param cart element to copy
	 *
	 */
	CartDecomposition<dim,T,Memory> & operator=(CartDecomposition && cart)
	{
		static_cast<ie_loc_ghost<dim,T>*>(this)->operator=(static_cast<ie_loc_ghost<dim,T>*>(cart));
		static_cast<nn_prcs<dim,T>*>(this)->operator=(static_cast<nn_prcs<dim,T>*>(cart));
		static_cast<ie_ghost<dim,T>*>(this)->operator=(static_cast<ie_ghost<dim,T>*>(cart));

		sub_domains.swap(cart.sub_domains);
		box_nn_processor.swap(cart.box_nn_processor);
		fine_s.swap(cart.fine_s);
		gr = cart.gr;
		cd = cart.cd;
		domain = cart.domain;
		std::copy(cart.spacing,cart.spacing+3,spacing);

		//! Runtime virtual cluster
		v_cl = cart.v_cl;

		ghost = cart.ghost;

		cart.bbox = bbox;
		cart.ss_box = ss_box;

		for (size_t i = 0 ; i < dim ; i++)
			cart.bc[i] = bc[i];

		return *this;
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
	template<typename Mem, typename ofb> size_t inline processorID(encapc<1, Point<dim,T>, Mem> p)
	{
		return fine_s.get(cd.template getCell<ofb>(p));
	}

	/*! \brief Given a point return in which processor the particle should go
	 *
	 * \return processorID
	 *
	 */
	size_t inline processorID(const Point<dim,T> &p) const
	{
		return fine_s.get(cd.getCell(p));
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

	/*! \brief Given a point return in which processor the particle should go
	 *
	 * Boundary conditions are considered
	 *
	 * \return processorID
	 *
	 */
	template<typename Mem> size_t inline processorIDBC(encapc<1, Point<dim,T>, Mem> p)
	{
		Point<dim,T> pt = p;
		applyPointBC(pt);

		return fine_s.get(cd.getCell(pt));
	}

	/*! \brief Given a point return in which processor the particle should go
	 *
	 * Boundary conditions are considered
	 *
	 * \return processorID
	 *
	 */
	template<typename ofb> size_t inline processorIDBC(const Point<dim,T> &p) const
	{
		Point<dim,T> pt = p;
		applyPointBC(pt);

		return fine_s.get(cd.getCell(p));
	}

	/*! \brief Given a point return in which processor the particle should go
	 *
	 * Boundary consition are considered
	 *
	 * \return processorID
	 *
	 */
	template<typename ofb> size_t inline processorIDBC(const T (&p)[dim]) const
	{
		Point<dim,T> pt = p;
		applyPointBC(pt);

		return fine_s.get(cd.getCell(p));
	}

	/*! \brief Get the smallest subdivision of the domain on each direction
	 *
	 * \return a box p1 is set to zero
	 *
	 */
	const ::Box<dim,T> & getSmallestSubdivision()
	{
		return ss_box;
	}

	/*! \brief Get the periodicity on i dimension
	 *
	 * \param i dimension
	 *
	 * \return the periodicity in direction i
	 *
	 */
	size_t isPeriodic(size_t i)
	{
		return bc[i];
	}

	/*! \brief Set the parameter of the decomposition
	 *
	 * \param div_ storing into how many domain to decompose on each dimension
	 * \param domain_ domain to decompose
	 *
	 */
	void setParameters(const size_t (& div_)[dim], ::Box<dim,T> domain_, const size_t (& bc)[dim] ,const Ghost<dim,T> & ghost)
	{
		// set the boundary conditions
		for (size_t i = 0 ; i < dim ; i++)
			this->bc[i] = bc[i];

		// set the ghost
		this->ghost = ghost;

		// Set the decomposition parameters
		gr.setDimensions(div_);
		domain = domain_;
		cd.setDimensions(domain, div_, 0);

		// init distribution
		dist.createCartGraph(gr, domain);

	}

	/*! \brief Start decomposition
	 *
	 */
	void decompose()
	{
		computeCommunicationAndMigrationCosts(1);

		dist.decompose();

		createSubdomains(v_cl,bc);
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
	 * \return true if the re-balance has been executed, false otherwise
	 */
	bool rebalance(DLB & dlb)
	{
		// if the DLB heuristic to use is the "Unbalance Threshold" get unbalance percentage
		if (dlb.getHeurisitc() == DLB::Heuristic::UNBALANCE_THRLD)
		{
			float unbalance = dist.getUnbalance();
			dlb.setUnbalance(unbalance);
			if (v_cl.getProcessUnitID() == 0)
			{
				//std::cout << std::setprecision(3) << unbalance << "\n";
			}
		}

		if (dlb.rebalanceNeeded())
		{
			computeCommunicationAndMigrationCosts(dlb.getNTimeStepSinceDLB());
			dist.refine();
			return true;
		}
		return false;
	}

	/*! \brief Get the current un-balance value
	 *
	 * \return the un-balance percentage value
	 */
	float getUnbalance()
	{
		return dist.getUnbalance();
	}

	/*! \brief Compute the processor load counting the total weights of its vertices
	 *
	 * \return the current processor load
	 */
	size_t getProcessorLoad()
	{
		return dist.getProcessorLoad();
	}

	/*! \brief function that return the position of the cell in the space
	 *
	 * \param id vertex id
	 * \param pos vector that will contain x, y, z
	 *
	 */
	inline void getSubSubDomainPosition(size_t id, T (&pos)[dim])
	{
		dist.getSubSubDomainPosition(id, pos);
	}

	/*! \brief Get the number of sub-sub-domains in this sub-graph
	 *
	 * @return number of sub-sub-domains in this sub-graph
	 */
	size_t getNSubSubDomains()
	{
		return dist.getNSubSubDomains();
	}

	/*! \brief function that set the weight of the vertex
	 *
	 * \param id vertex id
	 *
	 */
	inline void setSubSubDomainComputationCost(size_t id, size_t weight)
	{
		dist.setComputationCost(id, weight);
	}

	/*! \brief function that set the weight of the vertex
	 *
	 * \param id vertex id
	 *
	 */
	inline size_t getSubSubDomainComputationCost(size_t id)
	{
		return dist.getComputationCost(id);
	}

	/*! \brief Operator to access the size of the sub-graph
	 *
	 * \return the size of the subgraph
	 */
	size_t subSize()
	{
		return dist.subSize();
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

	/*! \brief Return the box of the physical domain
	 *
	 * \return The physical domain box
	 *
	 */
	const ::Box<dim,T> & getDomain()
	{
		return domain;
	}

	/*! \brief Check if the particle is local
	 *
	 * \warning if the particle id outside the domain the result is unreliable
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
	 * \warning if the particle id outside the domain the result is unreliable
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

	/*! \brief Check if the particle is local considering boundary conditions
	 *
	 * \warning if the particle id outside the domain and non periodic the result
	 *          is unreliable
	 *
	 *
	 * \param p object position
	 *
	 * \return true if it is local
	 *
	 */
	template<typename Mem> bool isLocalBC(const encapc<1, Point<dim,T>, Mem> p, const size_t (& bc)[dim]) const
	{
		Point<dim,T> pt = p;

		for (size_t i = 0 ; i < dim ; i++)
		{
			if (bc[i] == PERIODIC)
				pt.get(i) = openfpm::math::periodic_l(p.template get<0>()[i],domain.getHigh(i),domain.getLow(i));
		}

		return processorID<Mem>(pt) == v_cl.getProcessUnitID();
	}

	/*! \brief Check if the particle is local considering boundary conditions
	 *
	 * \param p object position
	 *
	 * \return true if it is local
	 *
	 */
	bool isLocalBC(const T (&p)[dim], const size_t (& bc)[dim]) const
	{
		Point<dim,T> pt = p;

		for (size_t i = 0 ; i < dim ; i++)
		{
			if (bc[i] == PERIODIC)
				pt.get(i) = openfpm::math::periodic_l(p[i],domain.getHigh(i),domain.getLow(i));
		}

		return processorID(pt) == v_cl.getProcessUnitID();
	}


	/*! \brief Return the bounding box containing union of all the sub-domains for the local processor
	 *
	 * \return The bounding box
	 *
	 */
	::Box<dim, T> & getProcessorBounds()
	{
		return bbox;
	}


	/*! \brief Return the ghost
	 *
	 *
	 */
	const Ghost<dim,T> & getGhost() const
	{
		return ghost;
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
		vtk_box1.write(output + std::string("subdomains_") + std::to_string(v_cl.getProcessUnitID()) + std::string(".vtk"));

		nn_prcs<dim, T>::write(output);
		ie_ghost<dim, T>::write(output, v_cl.getProcessUnitID());
		ie_loc_ghost<dim, T>::write(output, v_cl.getProcessUnitID());

		return true;
	}

	/*! \brief Get the Virtual Cluster machine
	 *
	 * \return the Virtual cluster machine
	 *
	 */
	Vcluster & getVC() const
	{
#ifdef SE_CLASS2
		check_valid(this,8);
#endif
		return v_cl;
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
				std::cout << ie_ghost<dim, T>::getProcessorEGhostBox(p, i).toString() << "   prc=" << nn_prcs<dim, T>::IDtoProc(p) << "   id=" << ie_ghost<dim, T>::getProcessorEGhostId(p, i) << "\n";
			}
		}

		std::cout << "Internal ghost box\n";

		for (size_t p = 0; p<nn_prcs < dim, T>::getNNProcessors(); p++)
		{
			for (size_t i = 0; i<ie_ghost < dim, T>::getProcessorNIGhost(p); i++)
			{
				std::cout << ie_ghost<dim, T>::getProcessorIGhostBox(p, i).toString() << "   prc=" << nn_prcs<dim, T>::IDtoProc(p) << "   id=" << ie_ghost<dim, T>::getProcessorIGhostId(p, i) << "\n";
			}
		}
	}

	/*! \brief Check if the CartDecomposition contain the same information
	 *
	 * \param ele Element to check
	 *
	 */
	bool is_equal(CartDecomposition<dim,T,Memory> & cart)
	{
		if (static_cast<ie_loc_ghost<dim,T>*>(this)->is_equal(static_cast<ie_loc_ghost<dim,T>&>(cart)) == false)
			return false;

		if (static_cast<nn_prcs<dim,T>*>(this)->is_equal(static_cast<nn_prcs<dim,T>&>(cart)) == false)
			return false;

		if (static_cast<ie_ghost<dim,T>*>(this)->is_equal(static_cast<ie_ghost<dim,T>&>(cart)) == false)
			return false;

		if (sub_domains != cart.sub_domains)
			return false;

		if (box_nn_processor != cart.box_nn_processor)
			return false;

		if (fine_s != cart.fine_s)
			return false;

		if (gr != cart.gr)
			return false;

		if (cd != cart.cd)
			return false;

		if (domain != cart.domain)
			return false;

		if (meta_compare<T[dim]>::meta_compare_f(cart.spacing,spacing) == false)
			return false;

		if (ghost != cart.ghost)
			return false;

		return true;
	}

	/*! \brief Check if the CartDecomposition contain the same information with the exception of the ghost part
	 * It is anyway required that the ghost come from the same sub-domains decomposition
	 *
	 * \param ele Element to check
	 *
	 */
	bool is_equal_ng(CartDecomposition<dim,T,Memory> & cart)
	{
		if (static_cast<ie_loc_ghost<dim,T>*>(this)->is_equal_ng(static_cast<ie_loc_ghost<dim,T>&>(cart)) == false)
			return false;

		if (static_cast<nn_prcs<dim,T>*>(this)->is_equal(static_cast<nn_prcs<dim,T>&>(cart)) == false)
			return false;

		if (static_cast<ie_ghost<dim,T>*>(this)->is_equal_ng(static_cast<ie_ghost<dim,T>&>(cart)) == false)
			return false;

		if (sub_domains != cart.sub_domains)
			return false;

		if (box_nn_processor != cart.box_nn_processor)
			return false;

		if (fine_s != cart.fine_s)
			return false;

		if (gr != cart.gr)
			return false;

		if (cd != cart.cd)
			return false;

		if (domain != cart.domain)
			return false;

		if (meta_compare<T[dim]>::meta_compare_f(cart.spacing,spacing) == false)
			return false;

		return true;
	}

	/*! \brief Return the distribution object
	 *
	 * \return the distribution object
	 *
	 */
	Distribution & getDistribution()
	{
		return dist;
	}
};


#endif
