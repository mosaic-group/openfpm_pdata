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
#include "NN/CellList/CellDecomposer.hpp"
#include <unordered_map>
#include "NN/CellList/CellList.hpp"
#include "Space/Ghost.hpp"

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
 *  \note sub-sub-domain portion of space at finer level than the sub-domain (before optimization)
 *        (or before sub-sub-domain merging)
 *  \note sub-domain portion of space (after optimization)
 *  \note near processor sub-domain a sub-domain that live in the a near (or contiguous) processor
 *
 */

template<unsigned int dim, typename T, template<typename> class device_l=openfpm::device_cpu, typename Memory=HeapMemory, template<unsigned int, typename> class Domain=Box, template<typename, typename, typename, typename, unsigned int> class data_s = openfpm::vector>
class CartDecomposition
{
	struct N_box
	{
		// id of the processor in the nn_processor list
		size_t id;

		// Near processor sub-domains
		typename openfpm::vector<::Box<dim,T>> bx;
	};

	struct Box_proc
	{
		// Intersection between the local sub-domain enlarged by the ghost and the contiguous processor
		// sub-domains
		openfpm::vector<::Box<dim,T>> bx;

		// Intersection between the contiguous processor sub-domain enlarged by the ghost with the
		// local sub-domain
		openfpm::vector<::Box<dim,T>> nbx;


		// processor
		size_t proc;
	};

public:

	//! Type of the domain we are going to decompose
	typedef T domain_type;

	//! It simplify to access the SpaceBox element
	typedef SpaceBox<dim,T> Box;

private:

	//! This is the key type to access  data_s, for example in the case of vector
	//! acc_key is size_t
	typedef typename data_s<SpaceBox<dim,T>,device_l<SpaceBox<dim,T>>,Memory,openfpm::vector_grow_policy_default,openfpm::vect_isel<SpaceBox<dim,T>>::value >::access_key acc_key;

	//! Subspace selected
	//! access_key in case of grid is just the set of the index to access the grid
	std::vector<acc_key> id_sub;

	//! the margin of the sub-domain selected
	SpaceBox<dim,T> sub_domain;

	//! the set of all local sub-domain as vector
	openfpm::vector<SpaceBox<dim,T>> sub_domains;

	//! List of near processors
	openfpm::vector<size_t> nn_processors;

	//! for each sub-domain, contain the list of the neighborhood processors
	//! and for each processor contain the boxes calculated from the intersection
	//! of the sub-domain + ghost with the near-by processor sub-domain ()
	openfpm::vector< openfpm::vector< Box_proc > > box_nn_processor_int;

	//! for each sub-domain, contain the list of the neighborhood processors
	openfpm::vector<openfpm::vector<long unsigned int> > box_nn_processor;

	// for each near-processor store the sub-domain of the near processor
	std::unordered_map<size_t, N_box> nn_processor_subdomains;

	//! Structure that contain for each sub-domain box the processor id
	//! exist for efficient global communication
	openfpm::vector<size_t> fine_s;

	//! Structure that store the cartesian grid information
	grid_sm<dim,void> gr;

	//! Structure that decompose your structure into cell without creating them
	//! useful to convert positions to CellId or sub-domain id in this case
	CellDecomposer_sm<dim,T> cd;

	//! rectangular domain to decompose
	Domain<dim,T> domain;

	//! Ghost boxes of the processor
	//! for each Sub-domain it store the ghost boxes, or
	//! the set of boxes that enclose the the ghost space
	//! Box cannot overlap, they contain one id that is the
	//! processor the information should come from
	openfpm::vector< openfpm::vector<Domain<dim,T>> > gh_dom;

	//! Internal boxes of the processor
	//! for each Sub-domain it store the boxes enclosing the
	//! space that must be communicated when another processor
	//! require the ghost
	//! Box can overlap, they contain one id that is the
	//! processor the information should be communicated to
	openfpm::vector< openfpm::vector< Domain<dim,T>> > int_box;

	//! Box Spacing
	T spacing[dim];

	//! Runtime virtual cluster machine
	Vcluster & v_cl;

	//! Structure that store the geometrical information about intersection between the local sub-domain
	//! and the near processor sub-domains
	CellList<dim,T,FAST> geo_cell;


	/*! \brief Create internally the decomposition
	 *
     * \param v_cl Virtual cluster, used internally to handle or pipeline communication
	 *
	 */
	void CreateDecomposition(Vcluster & v_cl)
	{
		// Calculate the total number of box and and the spacing
		// on each direction
		// Get the box containing the domain
		SpaceBox<dim,T> bs = domain.getBox();

		for (unsigned int i = 0; i < dim ; i++)
		{
			// Calculate the spacing
			spacing[i] = (bs.getHigh(i) - bs.getLow(i)) / gr.size(i);
		}

		// Here we use METIS
		// Create a cartesian grid graph
		CartesianGraphFactory<dim,Graph_CSR<nm_part_v,nm_part_e>> g_factory_part;

		// Processor graph
		Graph_CSR<nm_part_v,nm_part_e> gp = g_factory_part.template construct<NO_EDGE,T,dim-1>(gr.getSize(),domain);

		// Get the number of processing units
		size_t Np = v_cl.getProcessingUnits();

		// Get the processor id
		long int p_id = v_cl.getProcessUnitID();

		// Convert the graph to metis
		Metis<Graph_CSR<nm_part_v,nm_part_e>> met(gp,Np);

		// decompose
		met.decompose<nm_part_v::id>();

		// fill the structure that store the processor id for each sub-domain
		fine_s.resize(gr.size());

		// Optimize the decomposition creating bigger spaces
		// And reducing Ghost over-stress
		dec_optimizer<dim,Graph_CSR<nm_part_v,nm_part_e>> d_o(gp,gr.getSize());

		// set of Boxes produced by the decomposition optimizer
		openfpm::vector<::Box<dim,size_t>> loc_box;

		// optimize the decomposition
		d_o.template optimize<nm_part_v::sub_id,nm_part_v::id>(gp,p_id,loc_box,box_nn_processor);

		// produce the list of the contiguous processor (nn_processors) and link nn_processor_subdomains to the
		// processor list
		for (size_t i = 0 ;  i < box_nn_processor.size() ; i++)
		{
			for (size_t j = 0 ; j < box_nn_processor.get(i).size() ; j++)
			{
				nn_processors.add(box_nn_processor.get(i).get(j));
			}
		}

		// make the list sorted and unique
	    std::sort(nn_processors.begin(), nn_processors.end());
	    auto last = std::unique(nn_processors.begin(), nn_processors.end());
	    nn_processors.erase(last, nn_processors.end());

		// produce the list of the contiguous processor (nn_processors) and link nn_processor_subdomains to the
		// processor list
		for (size_t i = 0 ;  i < box_nn_processor.size() ; i++)
		{
			for (size_t j = 0 ; j < box_nn_processor.get(i).size() ; j++)
			{
				// processor id near to this sub-domain
				size_t proc_id = box_nn_processor.get(i).get(j);

				size_t k = 0;
				// search inside near processor list
				for (k = 0 ; k < nn_processors.size() ; k++)
					if (nn_processors.get(k) == proc_id)	break;

				nn_processor_subdomains[proc_id].id = k;
			}
		}

		// Initialize ss_box and bbox
		if (loc_box.size() >= 0)
		{
			SpaceBox<dim,T> sub_d(loc_box.get(0));
			sub_d.mul(spacing);
			sub_d.expand(spacing);

			// add the sub-domain
			sub_domains.add(sub_d);

			ss_box = sub_d;
			bbox = sub_d;
		}

		// convert into sub-domain
		for (size_t s = 1 ; s < loc_box.size() ; s++)
		{
			SpaceBox<dim,T> sub_d(loc_box.get(s));

			// re-scale and add spacing (the end is the starting point of the next domain + spacing)
			sub_d.mul(spacing);
			sub_d.expand(spacing);

			// add the sub-domain
			sub_domains.add(sub_d);

			// Calculate the bound box
			bbox.enclose(sub_d);

			// Create the smallest box contained in all sub-domain
			ss_box.contained(sub_d);
		}

		//++++++++++++++++++++++++++++++++++++++++ Debug output NN boxes
		{
		VTKWriter<openfpm::vector<::SpaceBox<dim,T>>,VECTOR_BOX> vtk_box1;
		vtk_box1.add(sub_domains);
		vtk_box1.write(std::string("loc_") + std::to_string(v_cl.getProcessUnitID()) + std::string(".vtk"));
		}
		//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

		// fill fine_s structure
		// fine_s structure contain the processor id for each sub-sub-domain
		// with sub-sub-domain we mean the sub-domain decomposition before
		// running dec_optimizer (before merging sub-domains)
		auto it = gp.getVertexIterator();

		while (it.isNext())
		{
			size_t key = it.get();

			// fill with the fine decomposition
			fine_s.get(key) = gp.template vertex_p<nm_part_v::id>(key);

			++it;
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
		grid_sm<dim,void> g(div);

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

	// Heap memory receiver
	HeapMemory hp_recv;

	// vector v_proc
	openfpm::vector<size_t> v_proc;

	// Receive counter
	size_t recv_cnt;

	/*! \brief Message allocation
	 *
	 * \param message size required to receive from i
	 * \param total message size to receive from all the processors
	 * \param the total number of processor want to communicate with you
	 * \param i processor id
	 * \param ri request id (it is an id that goes from 0 to total_p, and is unique
	 *           every time message_alloc is called)
	 * \param ptr a pointer to the vector_dist structure
	 *
	 * \return the pointer where to store the message
	 *
	 */
	static void * message_alloc(size_t msg_i ,size_t total_msg, size_t total_p, size_t i, size_t ri, void * ptr)
	{
		// cast the pointer
		CartDecomposition<dim,T,device_l,Memory,Domain,data_s> * cd = static_cast< CartDecomposition<dim,T,device_l,Memory,Domain,data_s> *>(ptr);

		if (cd->v_cl.getProcessUnitID() == 0)
		{
			std::cout << "Receiving from " << i << "       msg size: " << msg_i << "\n";
		}

		// Resize the memory
		cd->nn_processor_subdomains[i].bx.resize(msg_i);

		// Return the receive pointer
		return cd->nn_processor_subdomains[i].bx.getPointer();
	}

public:

	/*! \brief Cartesian decomposition copy constructor
	 *
     * \param v_cl Virtual cluster, used internally to handle or pipeline communication
	 *
	 */
	CartDecomposition(CartDecomposition<dim,T,device_l,Memory,Domain,data_s> && cd)
	:sub_domain(cd.sub_domain),gr(cd.gr),cd(cd.cd),domain(cd.domain),v_cl(cd.v_cl)
	{
		// Reset the box to zero
		bbox.zero();

		//! Subspace selected
		//! access_key in case of grid is just the set of the index to access the grid
		id_sub.swap(cd.id_sub);

		//! the set of all local sub-domain as vector
		sub_domains.swap(cd.sub_domains);

		for (int i = 0 ; i < dim ; i++)
		{
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
	:id_sub(0),v_cl(v_cl)
	{
		// Reset the box to zero
		bbox.zero();
	}

	/*! \brief Cartesian decomposition constructor, it divide the space in boxes
	 *
	 * \param dec is a vector that store how to divide on each dimension
	 * \param domain is the domain to divide
	 * \param v_cl are information of the cluster runtime machine
	 *
	 */
	CartDecomposition(std::vector<size_t> dec, Domain<dim,T> domain, Vcluster & v_cl)
	:id_sub(0),gr(dec),cd(domain,dec,0),domain(domain),v_cl(v_cl)
	{
		// Reset the box to zero
		bbox.zero();

		// Create the decomposition
		CreateDecomposition(v_cl);
	}

	//! Cartesian decomposition destructor
	~CartDecomposition()
	{}

	openfpm::vector<size_t> ids;

	/*! \brief Given a position it return if the position belong to any neighborhood processor ghost
	 *
	 * \param p Particle position
	 *
	 * \param return the processor ids
	 *
	 */
	const openfpm::vector<size_t> ghost_processorID(Point<dim,T> & p)
	{
		ids.clear();

		// Check with geo-cell if a particle is inside one Cell caotaining boxes

		auto cell_it = geo_cell.getCellIterator(p);

		// For each element in the cell, check if the point is inside the box
		// if it is store the processor id
		while (cell_it.isNext())
		{
			size_t bid = cell_it.get();

			if (vb_int.get(bid).box.isInside(p) == true)
			{
				ids.add(vb_int.get(bid).proc);
			}
		}

		return ids;
	}

	/*! \brief Given a position it return if the position belong to any neighborhood processor ghost
	 *
	 * \param p Particle position
	 *
	 * \param return the processor ids
	 *
	 */
	template<typename Mem> inline const openfpm::vector<size_t> ghost_processorID(const encapc<1,Point<dim,T>,Mem> & p)
	{
		ids.clear();

		// Check with geo-cell if a particle is inside one Cell containing boxes

		auto cell_it = geo_cell.getCellIterator(p);

		// For each element in the cell, check if the point is inside the box
		// if it is, store the processor id
		while (cell_it.isNext())
		{
			size_t bid = cell_it.get();

			if (vb_int.get(bid).box.isInside(p) == true)
			{
				ids.add(vb_int.get(bid).proc);
			}
		}

		return ids;
	}

	// Internal boxes for this processor domain, indicated with B8_0 B9_0 ..... in the figure
	// below as a linear vector
	openfpm::vector<::Box<dim,T>> vb_int;

	/*! It calculate the ghost boxes and internal boxes
	 *
	 * Example: Processor 10 calculate
	 * B8_0 B9_0 B9_1 and B5_0
	 *
	 *
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

       and also
       G8_0 G9_0 G9_1 G5_0

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


	 *
	 *
	 *
	 * \param ghost margins for each dimensions (p1 negative part) (p2 positive part)
	 *
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

	 *
	 *
	 */
	void calculateGhostBoxes(Ghost<dim,T> & ghost)
	{
#ifdef DEBUG
		// the ghost margins are assumed to be smaller
		// than one sub-domain

		for (size_t i = 0 ; i < dim ; i++)
		{
			if (ghost.template getLow(i) >= domain.template getHigh(i) / gr.size(i) || ghost.template getHigh(i)  >= domain.template getHigh(i) / gr.size(i))
			{
				std::cerr << "Error: Ghost are bigger that one domain" << "\n";
			}
		}
#endif

		// create a buffer with the sub-domains of this processor, the informations ( the boxes )
		// of the sub-domains contiguous to the processor A are sent to the processor A and
		// the information of the contiguous sub-domains in the near processors are received
		//
		openfpm::vector< openfpm::vector< ::SpaceBox<dim,T>> > boxes(nn_processors.size());

		for (size_t b = 0 ; b < box_nn_processor.size() ; b++)
		{
			for (size_t p = 0 ; p < box_nn_processor.get(b).size() ; p++)
			{
				size_t prc = box_nn_processor.get(b).get(p);

				// id of the processor in the processor list
				// [value between 0 and the number of the near processors]
				size_t id = nn_processor_subdomains[prc].id;

				boxes.get(id).add(sub_domains.get(b));
			}
		}

		//++++++++++++++++++++++++++++++++++++++++ Debug output NN boxes
		{
		for (size_t b = 0 ; b < boxes.size() ; b++)
		{
			VTKWriter<openfpm::vector<::SpaceBox<dim,T>>,VECTOR_BOX> vtk_box1;
			vtk_box1.add(boxes.get(b));
			vtk_box1.write(std::string("Processor_") + std::to_string(v_cl.getProcessUnitID()) + "_" + std::to_string(nn_processors.get(b)) + std::string(".vtk"));
		}
		}

		//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

		// Intersect all the local sub-domains with the sub-domains of the contiguous processors

		// Get the sub-domains of the near processors
		v_cl.sendrecvMultipleMessages(nn_processors,boxes,CartDecomposition<dim,T,device_l,Memory,Domain,data_s>::message_alloc, this ,NEED_ALL_SIZE);

		// ++++++++++++++++++++++++++++++++++++++++++ Check received boxes

		{
		VTKWriter<openfpm::vector<::Box<dim,T>>,VECTOR_BOX> vtk_box1;
		for (size_t p = 0 ; p < nn_processors.size() ; p++)
		{
			size_t prc = nn_processors.get(p);

			if (v_cl.getProcessUnitID() == 0)
				std::cout << "Received from " << prc << "      n_boxes: " << nn_processor_subdomains[prc].bx.size() << "\n";

			vtk_box1.add(nn_processor_subdomains[prc].bx);
		}
		vtk_box1.write(std::string("rb_Processor_") + std::to_string(v_cl.getProcessUnitID()) + std::string(".vtk"));
		}

		// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

		box_nn_processor_int.resize(sub_domains.size());

		// For each sub-domain
		for (size_t i = 0 ; i < sub_domains.size() ; i++)
		{
			SpaceBox<dim,T> sub_with_ghost = sub_domains.get(i);

			// enlarge the sub-domain with the ghost
			sub_with_ghost.enlarge(ghost);

			// resize based on the number of contiguous processors
			box_nn_processor_int.get(i).resize(box_nn_processor.get(i).size());

			// For each processor contiguous to this sub-domain
			for (size_t j = 0 ; j < box_nn_processor.get(i).size() ; j++)
			{
				// Contiguous processor
				size_t p_id = box_nn_processor.get(i).get(j);

				// get the set of sub-domains of the contiguous processor p_id
				openfpm::vector< ::Box<dim,T> > & p_box = nn_processor_subdomains[p_id].bx;

				// near processor sub-domain intersections
				openfpm::vector< ::Box<dim,T> > & p_box_int = box_nn_processor_int.get(i).get(j).bx;

				// for each near processor sub-domain intersect with the enlarged local sub-domain and store it
				for (size_t b = 0 ; b < p_box.size() ; b++)
				{
					::Box<dim,T> bi;

					bool intersect = sub_with_ghost.Intersect(::Box<dim,T>(p_box.get(b)),bi);

					if (intersect == true)
						p_box_int.add(bi);
				}
			}

			// For each processor contiguous to this sub-domain
			for (size_t j = 0 ; j < box_nn_processor.get(i).size() ; j++)
			{
				// Contiguous processor
				size_t p_id = box_nn_processor.get(i).get(j);

				// get the set of sub-domains of the contiguous processor p_id
				openfpm::vector< ::Box<dim,T> > & nn_p_box = nn_processor_subdomains[p_id].bx;

				// near processor sub-domain intersections
				openfpm::vector< ::Box<dim,T> > & p_box_int = box_nn_processor_int.get(i).get(j).nbx;

				// For each near processor sub-domains enlarge and intersect with the local sub-domain and store the result
				for (size_t k = 0 ; k < nn_p_box.size() ; k++)
				{
					// enlarge the near-processor sub-domain
					::Box<dim,T> n_sub = nn_p_box.get(k);

					// local sub-domain
					::SpaceBox<dim,T> l_sub = sub_domains.get(i);

					// Create a margin of ghost size around the near processor sub-domain
					n_sub.enlarge(ghost);

					// Intersect with the local sub-domain

					::Box<dim,T> b_int;
					bool intersect = n_sub.Intersect(l_sub,b_int);

					// store if it intersect
					if (intersect == true)
					{
						p_box_int.add(b_int);
						vb_int.add(b_int);

						// update the geo_cell list

						// get the boxes this box span
						const grid_key_dx<dim> p1 = geo_cell.getCellGrid(b_int.getP1());
						const grid_key_dx<dim> p2 = geo_cell.getCellGrid(b_int.getP2());

						// Get the grid and the sub-iterator
						auto & gi = geo_cell.getGrid();
						grid_key_dx_iterator_sub<dim> g_sub(gi,p1,p2);

						// add the box-id to the cell list
						while (g_sub.isNext())
						{
							auto key = g_sub.get();
							geo_cell.addCell(gi.LinId(key),vb_int.size()-1);
							++g_sub;
						}
					}
				}
			}


			// ++++++++++++++++++++++++++++++++++++++++ Debug +++++++++++++++++++++++++++++

			{
			VTKWriter<openfpm::vector<::Box<dim,T>>,VECTOR_BOX> vtk_box1;
			for (size_t p = 0 ; p < box_nn_processor_int.size() ; p++)
			{
				for (size_t s = 0 ; s < box_nn_processor_int.get(p).size() ; s++)
				{
					vtk_box1.add(box_nn_processor_int.get(p).get(s).nbx);
				}
			}
			vtk_box1.write(std::string("inte_Processor_") + std::to_string(v_cl.getProcessUnitID()) + std::string(".vtk"));
			}

			// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

		}
	}

	/*! \brief processorID return in which processor the particle should go
	 *
	 * \return processorID
	 *
	 */

	template<typename Mem> size_t inline processorID(encapc<1, Point<dim,T>, Mem> p)
	{
		return fine_s.get(cd.getCell(p));
	}

	// Smallest subdivision on each direction
	::Box<dim,T> ss_box;

	/*! \brief Get the smallest subdivision of the domain on each direction
	 *
	 * \return a box p1 is set to zero
	 *
	 */
	const ::Box<dim,T> & getSmallestSubdivision()
	{
		return ss_box;
	}

	/*! \brief processorID return in which processor the particle should go
	 *
	 * \return processorID
	 *
	 */

	size_t inline processorID(T (&p)[dim])
	{
		return fine_s.get(cd.getCell(p));
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
	void setParameters(const size_t (& div_)[dim], Domain<dim,T> domain_)
	{
		// Set the decomposition parameters

		gr.setDimensions(div_);
		domain = domain_;
		cd.setDimensions(domain,div_,0);

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
		if (id >= gr.size())
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
		return gr.size();
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

	::Box<dim,T> bbox;

	/*! \brief Return the bounding box containing the processor box + smallest subdomain spacing
	 *
	 * \return The bounding box
	 *
	 */
	::Box<dim,T> & getProcessorBounds()
	{
		return bbox;
	}

	/*! \brief if the point fall into the ghost of some near processor it return the processors id's in which
	 *  it fall
	 *
	 * \param p Point
	 * \return iterator of the processors id's
	 *
	 */
	inline auto labelPoint(Point<dim,T> & p) -> decltype(geo_cell.getIterator(geo_cell.getCell(p)))
	{
		return geo_cell.getIterator(geo_cell.getCell(p));
	}

	/*! \brief if the point fall into the ghost of some near processor it return the processor number in which
	 *  it fall
	 *
	 * \param p Point
	 * \return number of processors
	 *
	 */
	inline size_t labelPointNp(Point<dim,T> & p)
	{
		return geo_cell.getNelements(geo_cell.getCell(p));
	}

	/*! \brief It return the label point cell
	 *
	 * The labeling of a point p is regulated by a Cell list, give a point it give a cell-id
	 *
	 * \param p Point
	 * \return cell-id
	 *
	 */
	inline size_t labelPointCell(Point<dim,T> & p)
	{
		return geo_cell.getCell(p);
	}

	/*! \brief Fill the ghost buffer
	 *
	 * \tparam one or more properties to get
	 *
	 */
/*	template<unsigned int ...i> void ghost_get()
	{
		// first check if a local particle must be sent to another processor
		for (size_t i = 0 ; i < ; i++)
		{

		}
	}*/

	/*! \brief Fill the ghost buffer
	 *
	 * \tparam one or more properties to get
	 *
	 */
/*	template<unsigned int ...i> void ghost_put()
	{

	}*/

};


#endif
