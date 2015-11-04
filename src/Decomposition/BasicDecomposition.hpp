/*
 * BasicDecomposition.hpp
 *
 *  Created on: Oct 07, 2015
 *      Author: Antonio Leo
 */

#ifndef BASICDECOMPOSITION_HPP
#define BASICDECOMPOSITION_HPP

#include "config.h"
#include <cmath>
#include "VCluster.hpp"
#include "Graph/CartesianGraphFactory.hpp"
#include "Decomposition.hpp"
#include "Vector/map_vector.hpp"
#include <vector>
#include <initializer_list>
#include "SubdomainGraphNodes.hpp"
#include "parmetis_util.hpp"
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

#define BASICDEC_ERROR 2000lu

// Macro that decide what to do in case of error
#ifdef STOP_ON_ERROR
#define ACTION_ON_ERROR() exit(1);
#elif defined(THROW_ON_ERROR)
#define ACTION_ON_ERROR() throw BASICDEC_ERROR;
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
 * \snippet BasicDecomposition_unit_test.hpp Create BasicDecomposition
 *
 */

template<unsigned int dim, typename T, typename Memory = HeapMemory, template<unsigned int, typename > class Domain = Box>
class BasicDecomposition: public ie_loc_ghost<dim, T>, public nn_prcs<dim, T>, public ie_ghost<dim, T> {

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

	//! Cell-list that store the geometrical information of the local internal ghost boxes
	CellList<dim, T, FAST> lgeo_cell;

	//! Convert the graph to parmetis format
	Parmetis<Graph_CSR<nm_v, nm_e>> parmetis_graph;
	
	//! Processor sub-sub-domain graph
	Graph_CSR<nm_v, nm_e> sub_g;
	
	//! Global sub-sub-domain graph
	Graph_CSR<nm_v, nm_e> gp;
	
	//! Init vtxdist needed for Parmetis
	openfpm::vector<idx_t> vtxdist;

	//! partitions
	openfpm::vector<openfpm::vector<idx_t>> partitions;
	
	//! Init data structure to keep trace of new vertices distribution in processors (needed to update main graph)
	openfpm::vector<openfpm::vector<size_t>> v_per_proc;
	
	static void * message_receive(size_t msg_i ,size_t total_msg, size_t total_p, size_t i, size_t ri, void * ptr)
	{
		openfpm::vector<openfpm::vector<idx_t>> * v = static_cast<openfpm::vector<openfpm::vector<idx_t>> *>(ptr);

		v->get(i).resize(msg_i/sizeof(idx_t));
		
		return &(v->get(i).get(0));
	}
	
	/*! \brief Constructor, it decompose and distribute the sub-domains across the processors
	 *
	 * \param v_cl Virtual cluster, used internally for communications
	 *
	 */
	void CreateDecomposition(Vcluster & v_cl) {
#ifdef SE_CLASS1
		if (&v_cl == NULL)
		{
			std::cerr << __FILE__ << ":" << __LINE__ << " error VCluster instance is null, check that you ever initialized it \n";
			ACTION_ON_ERROR()
		}
#endif
		// Calculate the total number of box and and the spacing
		// on each direction
		// Get the box containing the domain
		SpaceBox<dim, T> bs = domain.getBox();

		for (unsigned int i = 0; i < dim; i++) {
			// Calculate the spacing
			spacing[i] = (bs.getHigh(i) - bs.getLow(i)) / gr.size(i);
		}

		//! Get the processor id
		size_t p_id = v_cl.getProcessUnitID();

		// Here we use PARMETIS
		// Create a cartesian grid graph
		CartesianGraphFactory<dim, Graph_CSR<nm_v, nm_e>> g_factory_part;

		// Processor graph
		gp = g_factory_part.template construct<NO_EDGE, nm_v::id, T, dim - 1, 0, 1, 2>(gr.getSize(), domain);

		//! Add computation information to each vertex (shape can change)
		addWeights(gp, HYPERBOLOID);

		//! Init ids vector
		gp.init_map_ids();

		// Get the number of processing units
		size_t Np = v_cl.getProcessingUnits();

		// Division of vertices in Np graphs
		// Put (div+1) vertices in mod graphs
		// Put div vertices in the rest of the graphs
		size_t mod_v = gp.getNVertex() % Np;
		size_t div_v = gp.getNVertex() / Np;

		for (int i = 0; i <= Np; i++) 
		{
			if (i < mod_v)
				vtxdist.get(i) = (div_v + 1) * (i);
			else
				vtxdist.get(i) = (div_v) * (i) + mod_v;
		}

		// Just for output purpose
		if (p_id == 0) {
			VTKWriter<Graph_CSR<nm_v, nm_e>, GRAPH> gv2(gp);
			gv2.write("test_graph_0.vtk");
		}

		//TODO transform in factory


		//! Put vertices into processor graph (different for each processor)
		fillSubGraph(sub_g, gp, vtxdist, p_id, Np);

		parmetis_graph.initSubGraph(sub_g);
		
		// Decompose
		parmetis_graph.decompose<nm_v::proc_id>(vtxdist,sub_g);

		// Get result partition for this processors
		idx_t *partition = parmetis_graph.getPartition();

		// Prepare vector of arrays to contain all partitions
		
		partitions.get(p_id).resize(sub_g.getNVertex());
		std::copy(partition,partition+sub_g.getNVertex(),&partitions.get(p_id).get(0));
		
		openfpm::vector<size_t> prc;
		openfpm::vector<size_t> sz;
		openfpm::vector<void *> ptr;
		
		for (size_t i = 0 ; i < Np ; i++)
		{
			if (i != v_cl.getProcessUnitID())
			{
			  prc.add(i);
			  sz.add(sub_g.getNVertex() * sizeof(idx_t));
			  ptr.add(partitions.get(p_id).getPointer());
			}
		}
		
		v_cl.sendrecvMultipleMessagesNBX(prc.size(),&sz.get(0),&prc.get(0),&ptr.get(0),message_receive,&partitions,NONE);
		
		// Update graphs with the new distributions
		updateGraphs(partitions, gp, sub_g, v_per_proc, vtxdist, p_id, Np);
		
		if (p_id == 0) {
			VTKWriter<Graph_CSR<nm_v, nm_e>, GRAPH> gv2(gp);
			gv2.write("test_graph_2.vtk");
		}

		// Renumbering subgraph
		sub_g.reset_map_ids();
		for (size_t j = vtxdist.get(p_id), i = 0; j < vtxdist.get(p_id + 1); j++, i++) {
			sub_g.set_map_ids(j, sub_g.vertex(i).template get<nm_v::id>());
			sub_g.vertex(i).template get<nm_v::id>() = j;
		}

		gp.reset_map_ids();
		// Renumbering main graph
		for (size_t p = 0; p < Np; p++) 
		{
			for (size_t j = vtxdist.get(p), i = 0; j < vtxdist.get(p + 1); j++, i++) 
			{
				gp.set_map_ids(j, gp.vertex(v_per_proc.get(p).get(i)).template get<nm_v::id>());
				gp.vertex(v_per_proc.get(p).get(i)).template get<nm_v::id>() = j;
			}
		}

		if (p_id == 0) 
		{
			VTKWriter<Graph_CSR<nm_v, nm_e>, GRAPH> gv2(gp);
			gv2.write("test_graph_3.vtk");
		}
		
		refine();

		/*
		 // fill the structure that store the processor id for each sub-domain
		 fine_s.resize(gr.size());

		 // Optimize the decomposition creating bigger spaces
		 // And reducing Ghost over-stress
		 dec_optimizer<dim,Graph_CSR<nm_v,nm_e>> d_o(gp,gr.getSize());

		 // set of Boxes produced by the decomposition optimizer
		 openfpm::vector<::Box<dim,size_t>> loc_box;

		 // optimize the decomposition
		 d_o.template optimize<nm_v::sub_id,nm_v::id>(gp,p_id,loc_box,box_nn_processor);

		 // Initialize ss_box and bbox
		 if (loc_box.size() >= 0)
		 {
		 SpaceBox<dim,size_t> sub_dc = loc_box.get(0);
		 SpaceBox<dim,T> sub_d(sub_dc);
		 sub_d.mul(spacing);
		 sub_d.expand(spacing);

		 // Fixing sub-domains to cover all the domain

		 // Fixing sub_d
		 // if (loc_box) is a the boundary we have to ensure that the box span the full
		 // domain (avoiding rounding off error)
		 for (size_t i = 0 ; i < dim ; i++)
		 {
		 if (sub_dc.getHigh(i) == cd.getGrid().size(i) - 1)
		 {
		 sub_d.setHigh(i,domain.getHigh(i));
		 }
		 }

		 // add the sub-domain
		 sub_domains.add(sub_d);

		 ss_box = sub_d;
		 ss_box -= ss_box.getP1();
		 bbox = sub_d;
		 }

		 // convert into sub-domain
		 for (size_t s = 1 ; s < loc_box.size() ; s++)
		 {
		 SpaceBox<dim,size_t> sub_dc = loc_box.get(s);
		 SpaceBox<dim,T> sub_d(sub_dc);

		 // re-scale and add spacing (the end is the starting point of the next domain + spacing)
		 sub_d.mul(spacing);
		 sub_d.expand(spacing);

		 // Fixing sub-domains to cover all the domain

		 // Fixing sub_d
		 // if (loc_box) is a the boundary we have to ensure that the box span the full
		 // domain (avoiding rounding off error)
		 for (size_t i = 0 ; i < dim ; i++)
		 {
		 if (sub_dc.getHigh(i) == cd.getGrid().size(i) - 1)
		 {
		 sub_d.setHigh(i,domain.getHigh(i));
		 }
		 }

		 // add the sub-domain
		 sub_domains.add(sub_d);

		 // Calculate the bound box
		 bbox.enclose(sub_d);

		 // Create the smallest box contained in all sub-domain
		 ss_box.contained(sub_d);
		 }

		 nn_prcs<dim,T>::create(box_nn_processor, sub_domains);

		 // fill fine_s structure
		 // fine_s structure contain the processor id for each sub-sub-domain
		 // with sub-sub-domain we mean the sub-domain decomposition before
		 // running dec_optimizer (before merging sub-domains)
		 it = gp.getVertexIterator();

		 while (it.isNext())
		 {
		 size_t key = it.get();

		 // fill with the fine decomposition
		 fine_s.get(key) = gp.template vertex_p<nm_v::id>(key);

		 ++it;
		 }

		 // Get the smallest sub-division on each direction
		 ::Box<dim,T> unit = getSmallestSubdivision();
		 // Get the processor bounding Box
		 ::Box<dim,T> bound = getProcessorBounds();

		 // calculate the sub-divisions
		 size_t div[dim];
		 for (size_t i = 0 ; i < dim ; i++)
		 div[i] = (size_t)((bound.getHigh(i) - bound.getLow(i)) / unit.getHigh(i));

		 // Create shift
		 Point<dim,T> orig;

		 // p1 point of the Processor bound box is the shift
		 for (size_t i = 0 ; i < dim ; i++)
		 orig.get(i) = bound.getLow(i);

		 // Initialize the geo_cell structure
		 ie_ghost<dim,T>::Initialize_geo_cell(domain,div,orig);
		 lgeo_cell.Initialize(domain,div,orig);
		 */
	}

	// Save the ghost boundaries
	Ghost<dim, T> ghost;

	void refine()
	{
		size_t Np = v_cl.getProcessingUnits();
		size_t p_id = v_cl.getProcessUnitID();
	  
		// Reset parmetis graph and reconstruct it
		parmetis_graph.reset(gp,sub_g);

		// Refine
		parmetis_graph.refine<nm_v::proc_id>(vtxdist,sub_g);

		// Get result partition for this processor
		idx_t * partition = parmetis_graph.getPartition();

		partitions.get(p_id).resize(sub_g.getNVertex());
		std::copy(partition,partition+sub_g.getNVertex(),&partitions.get(p_id).get(0));

		// Reset data structure to keep trace of new vertices distribution in processors (needed to update main graph)
		for (int i = 0; i < Np; ++i) 
		{
			v_per_proc.get(i).clear();
		}

		
		openfpm::vector<size_t> prc;
		openfpm::vector<size_t> sz;
		openfpm::vector<void *> ptr;

		for (size_t i = 0 ; i < Np ; i++)
		{
			if (i != v_cl.getProcessUnitID())
			{
			  partitions.get(i).clear();
			  prc.add(i);
			  sz.add(sub_g.getNVertex() * sizeof(idx_t));

//			  std::cout << "sub_g: " << sub_g.getNVertex() * sizeof(idx_t) << "\n";
			  
			  ptr.add(partitions.get(p_id).getPointer());
			}
		}
		
		v_cl.sendrecvMultipleMessagesNBX(prc.size(),&sz.get(0),&prc.get(0),&ptr.get(0),message_receive,&partitions,NONE);
		
		// Update graphs with the new distributions
		updateGraphs(partitions, gp, sub_g, v_per_proc, vtxdist, p_id, Np);
		
		//

		if (p_id == 1) 
		{
			VTKWriter<Graph_CSR<nm_v, nm_e>, GRAPH> gv2(gp);
			gv2.write("test_graph_4.vtk");
			bool test = compare("test_graph_4.vtk","test_graph_test.vtk");
			BOOST_REQUIRE_EQUAL(test,true);
		}
	}
	
	/*! \brief Create the subspaces that decompose your domain
	 *
	 */
	void CreateSubspaces() {
		// Create a grid where each point is a space
		grid_sm<dim, void> g(div);

		// create a grid_key_dx iterator
		grid_key_dx_iterator < dim > gk_it(g);

		// Divide the space into subspaces
		while (gk_it.isNext()) {
			//! iterate through all subspaces
			grid_key_dx < dim > key = gk_it.get();

			//! Create a new subspace
			SpaceBox<dim, T> tmp;

			//! fill with the Margin of the box
			for (int i = 0; i < dim; i++) {
				tmp.setHigh(i, (key.get(i) + 1) * spacing[i]);
				tmp.setLow(i, key.get(i) * spacing[i]);
			}

			//! add the space box
			sub_domains.add(tmp);

			// add the iterator
			++gk_it;
		}
	}

	/* /brief types of weights distributions
	 *
	 */
	enum weightShape {
		UNIFORM, SPHERE, HYPERBOLOID
	};

	/* \brief add vertex weights to the main domain, follow a shape
	 *
	 ** 0 - weights are all 1 on all vertices
	 ** 1 - weights are distributed as a sphere
	 *
	 * \param i id of the shape
	 *
	 */
	void addWeights(Graph_CSR<nm_v, nm_e> & gp, int i)
	{
		float c_x = 0, c_y = 0, c_z = 0 , radius2, eq;
		float x = 0, y = 0, z = 0;

		switch (i) {
			case UNIFORM:

				// Add computation information to each vertex
				for (int i = 0; i < gp.getNVertex(); i++) {
					gp.vertex(i).template get<nm_v::computation>() = 1;
				}
				break;
			case SPHERE:

				// Fill vertices weights with a sphere (if dim=2 is a circle)
				radius2 = pow(4, 2);
				c_x = 2;
				c_y = 2;

				if(dim == 3)
					c_z = 2;

				for (int i = 0; i < gp.getNVertex(); i++) {
					x = gp.vertex(i).template get<nm_v::x>() * 10;
					y = gp.vertex(i).template get<nm_v::y>() * 10;

					if(dim == 3)
						z = gp.vertex(i).template get<nm_v::z>() * 10;

					eq = pow((x - c_x), 2) + pow((y - c_y), 2) + pow((z - c_z), 2);

					if (eq <= radius2) {
						gp.vertex(i).template get<nm_v::computation>() = 5;
					} else {
						gp.vertex(i).template get<nm_v::computation>() = 1;
					}
				}
				break;
			case HYPERBOLOID:

				// Fill vertices weights with a elliptic hyperboloid (if dim=2 is an hyperbole)
				c_x = 5;
				c_y = 5;

				if(dim == 3)
					c_z = 5;
				for (int i = 0; i < gp.getNVertex(); i++) {
					x = gp.vertex(i).template get<nm_v::x>() * 10;
					y = gp.vertex(i).template get<nm_v::y>() * 10;

					if(dim == 3)
						z = gp.vertex(i).template get<nm_v::z>() * 10;

					eq = - pow((x - c_x), 2)/3 - pow((y - c_y), 2)/3 + pow((z - c_z), 2)/2;

					if (eq >= 1) {
						gp.vertex(i).template get<nm_v::computation>() = 5;
					} else {
						gp.vertex(i).template get<nm_v::computation>() = 1;
					}
				}
				break;
		}
	}

	/* \brief fill the graph of the processor with the first decomposition (linear)
	 *
	 * Put vertices into processor graph (different for each processor)
	 *
	 * \param sub_g sub graph to fill
	 * \param gp mai graph, source for the vertices
	 * \param vtxdist array with the distribution of vertices through processors
	 * \param proc_id rank of the processor
	 * \param Np total number of processors
	 */
	void fillSubGraph(Graph_CSR<nm_v, nm_e> &sub_g, Graph_CSR<nm_v, nm_e> &gp, openfpm::vector<idx_t> &vtxdist, int proc_id, int Np)
	{

		for (size_t j = vtxdist.get(proc_id), local_j = 0; j < vtxdist.get(proc_id + 1); j++, local_j++) {

			// Add vertex
			nm_v pv = gp.vertexById(j);
			sub_g.addVertex(pv);

			// Add edges of vertex
			for (size_t s = 0; s < gp.getNChilds(j); s++) {
				sub_g.template addEdge<NoCheck>(local_j, gp.getChildByVertexId(j, s));
			}
		}

		// Just for output purpose
		if (proc_id == 0) {
			for (int i = 0; i < Np; i++) {
				for (size_t j = vtxdist.get(i); j < vtxdist.get(i + 1); j++) {
					gp.vertexById(j).template get<nm_v::proc_id>() = i;
				}
			}
			VTKWriter<Graph_CSR<nm_v, nm_e>, GRAPH> gv2(gp);
			gv2.write("test_graph_1.vtk");
		}
	}

	/* \brief exchange partitions with other processors
	 *
	 * \param partitions array to store all the partitions
	 * \param gp_nv number of vertices on main graph
	 * \param sub_g_nv number of vertices on sub graph
	 * \param requests_recv array of requests
	 * \param requests_send array of requests
	 * \param statuses array of statsu objects
	 * \param proc_id current processors rank
	 * \param Np total umber of processors
	 */
	void exchangePartitions(int** &partitions, int gp_nv, int sub_g_nv, MPI_Request* &requests_recv, MPI_Request* &requests_send,
			MPI_Status* &statuses, int proc_id, int Np)
	{

		// Receive other partitions, each partition can contain max NVertex of main graph
		for (int i = 0; i < Np; i++) {
			if (i != proc_id)
				MPI_Irecv(partitions[i], gp_nv, MPI_INT, i, 0, MPI_COMM_WORLD, &requests_recv[i]);
		}

		// Send processor partition to other processors
		for (int i = 0; i < Np; i++) {
			if (i != proc_id)
				MPI_Isend(partitions[proc_id], sub_g_nv, MPI_INT, i, 0, MPI_COMM_WORLD, &requests_send[i]);
		}

		// Wait for all partitions from other processors
		for (int i = 0; i < Np; i++) {
			if (i != proc_id)
				MPI_Wait(&requests_recv[i], &statuses[i]);
		}
	}

	/* \brief update main graph ad subgraph with the partition in partitions param
		 *
		 * \param partitions array storing all the partitions
		 * \param gp main graph
		 * \param sub_g sub graph
		 * \param v_per_proc array needed to recontruct the main graph
		 * \param vtxdist array with the distribution of vertices through processors
		 * \param statuses array of statsu objects
		 * \param proc_id current processors rank
		 * \param Np total umber of processors
		 */
	void updateGraphs(openfpm::vector<openfpm::vector<idx_t>> &partitions,Graph_CSR<nm_v, nm_e> &gp, Graph_CSR<nm_v, nm_e> &sub_g, openfpm::vector<openfpm::vector<size_t>> & v_per_proc, openfpm::vector<idx_t> & vtxdist, int proc_id, int Np) 
	{
		int local_j = 0;
		sub_g.clear();

		// Init n_vtxdist to gather informations about the new decomposition
/*		idx_t *n_vtxdist = new idx_t[Np + 1];*/

		openfpm::vector<idx_t> n_vtxdist(Np+1);
		for (int i = 0; i <= Np; i++)
			n_vtxdist.get(i) = 0;
		
		// Update main graph with other partitions made by Parmetis in other processors and the local partition
		for (int i = 0; i < Np; i++) {

			int ndata = partitions.get(i).size();

			// Update the main graph with received informations
			for (int k = 0, l = vtxdist.get(i); k < ndata && l < vtxdist.get(i + 1); k++, l++) {

				// Create new n_vtxdist (1) (just count processors vertices)
				n_vtxdist.get(partitions.get(i).get(k) + 1)++;

				// Update proc id in the vertex
				gp.vertexById(l).template get<nm_v::proc_id>() = partitions.get(i).get(k);

				// Add vertex to temporary structure of distribution (needed to update main graph)
				v_per_proc.get(partitions.get(i).get(k)).add(gp.vertexById(l).template get<nm_v::id>());

				// Add vertices belonging to this processor in sub graph
				if (partitions.get(i).get(k) == proc_id) {

					nm_v pv = gp.vertexById(l);
					sub_g.addVertex(pv);

					// Add edges of vertex
					for (size_t s = 0; s < gp.getNChildsByVertexId(l); s++) {
						sub_g.template addEdge<NoCheck>(local_j, gp.getChildByVertexId(l, s));
					}

					local_j++;
				}
			}
		}

		// Create new n_vtxdist (2) (write boundaries)
		for (int i = 2; i <= Np; i++) 
		{
			n_vtxdist.get(i) += n_vtxdist.get(i - 1);
		}

		// Copy the new decomposition in the main vtxdist
		for (int i = 0; i <= Np ; i++) 
		{
			vtxdist.get(i) = n_vtxdist.get(i);
		}
	}

	// Heap memory receiver
	HeapMemory hp_recv;

	// vector v_proc
	openfpm::vector<size_t> v_proc;

	// Receive counter
	size_t recv_cnt;

public:

	/*! \brief Basic decomposition constructor
	 *
	 * \param v_cl Virtual cluster, used internally to handle or pipeline communication
	 *
	 */
	BasicDecomposition(Vcluster & v_cl) :
	nn_prcs<dim, T>(v_cl), v_cl(v_cl),parmetis_graph(v_cl, v_cl.getProcessingUnits()),vtxdist(v_cl.getProcessingUnits() + 1),partitions(v_cl.getProcessingUnits()),v_per_proc(v_cl.getProcessingUnits())
	{
		// Reset the box to zero
		bbox.zero();
	}

	//! Basic decomposition destructor
	~BasicDecomposition() {
	}

	//	openfpm::vector<size_t> ids;

	/*! \brief class to select the returned id by ghost_processorID
	 *
	 */
	class box_id {
	public:
		/*! \brief Return the box id
		 *
		 * \param p structure containing the id informations
		 * \param b_id box_id
		 *
		 * \return box id
		 *
		 */
		inline static size_t id(p_box<dim, T> & p, size_t b_id) {
			return b_id;
		}
	};

	/*! \brief class to select the returned id by ghost_processorID
	 *
	 */
	class processor_id {
	public:
		/*! \brief Return the processor id
		 *
		 * \param p structure containing the id informations
		 * \param b_id box_id
		 *
		 * \return processor id
		 *
		 */
		inline static size_t id(p_box<dim, T> & p, size_t b_id) {
			return p.proc;
		}
	};

	/*! \brief class to select the returned id by ghost_processorID
	 *
	 */
	class lc_processor_id {
	public:
		/*! \brief Return the near processor id
		 *
		 * \param p structure containing the id informations
		 * \param b_id box_id
		 *
		 * \return local processor id
		 *
		 */
		inline static size_t id(p_box<dim, T> & p, size_t b_id) {
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
	void calculateGhostBoxes() {
#ifdef DEBUG
		// the ghost margins are assumed to be smaller
		// than one sub-domain

		for (size_t i = 0; i < dim; i++) {
			if (ghost.template getLow(i) >= domain.template getHigh(i) / gr.size(i)
					|| ghost.template getHigh(i) >= domain.template getHigh(i) / gr.size(i)) {
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
		for (size_t i = 0; i < dim; i++) {
			if (ghost.template getLow(i) >= ss_box.getHigh(i) || ghost.template getHigh(i) >= domain.template getHigh(i) / gr.size(i)) {
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
	static size_t getDefaultGrid(size_t n_sub) {
		// Calculate the number of sub-sub-domain on
		// each dimension
		return openfpm::math::round_big_2(pow(n_sub, 1.0 / dim));
	}

	/*! \brief Given a point return in which processor the particle should go
	 *
	 * \return processorID
	 *
	 */
	template<typename Mem> size_t inline processorID(encapc<1, Point<dim, T>, Mem> p) {
		return fine_s.get(cd.getCell(p));
	}

	// Smallest subdivision on each direction
	::Box<dim, T> ss_box;

	/*! \brief Get the smallest subdivision of the domain on each direction
	 *
	 * \return a box p1 is set to zero
	 *
	 */
	const ::Box<dim, T> & getSmallestSubdivision() {
		return ss_box;
	}

	/*! \brief Given a point return in which processor the particle should go
	 *
	 * \return processorID
	 *
	 */

	size_t inline processorID(const T (&p)[dim]) const {
		return fine_s.get(cd.getCell(p));
	}

	/*! \brief Set the parameter of the decomposition
	 *
	 * \param div_ storing into how many domain to decompose on each dimension
	 * \param domain_ domain to decompose
	 *
	 */
	void setParameters(const size_t (&div_)[dim], Domain<dim, T> domain_, Ghost<dim, T> ghost = Ghost<dim, T>()) {
		// set the ghost
		this->ghost = ghost;
		// Set the decomposition parameters

		gr.setDimensions(div_);
		domain = domain_;
		cd.setDimensions(domain, div_, 0);

		//! Create the decomposition

		CreateDecomposition(v_cl);
	}

	/*! \brief Get the number of local sub-domains
	 *
	 * \return the number of sub-domains
	 *
	 */
	size_t getNLocalHyperCube() {
		return sub_domains.size();
	}

	/*! \brief Get the local sub-domain
	 *
	 * \param i (each local processor can have more than one sub-domain)
	 * \return the sub-domain
	 *
	 */
	SpaceBox<dim, T> getLocalHyperCube(size_t lc) {
		// Create a space box
		SpaceBox<dim, T> sp;

		// fill the space box

		for (size_t k = 0; k < dim; k++) {
			// create the SpaceBox Low and High
			sp.setLow(k, sub_domains.template get < Box::p1 > (lc)[k]);
			sp.setHigh(k, sub_domains.template get < Box::p2 > (lc)[k]);
		}

		return sp;
	}

	/*! \brief Get the local sub-domain with ghost extension
	 *
	 * \param i (each local processor can have more than one sub-domain)
	 * \return the sub-domain
	 *
	 */
	SpaceBox<dim, T> getSubDomainWithGhost(size_t lc) {
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
	Domain<dim, T> & getDomain() {
		return domain;
	}

	/*! \brief Check if the particle is local
	 *
	 * \param p object position
	 *
	 * \return true if it is local
	 *
	 */
	template<typename Mem> bool isLocal(const encapc<1, Point<dim, T>, Mem> p) const {
		return processorID<Mem>(p) == v_cl.getProcessUnitID();
	}

	/*! \brief Check if the particle is local
	 *
	 * \param p object position
	 *
	 * \return true if it is local
	 *
	 */
	bool isLocal(const T (&pos)[dim]) const {
		return processorID(pos) == v_cl.getProcessUnitID();
	}

	::Box<dim, T> bbox;

	/*! \brief Return the bounding box containing union of all the sub-domains for the local processor
	 *
	 * \return The bounding box
	 *
	 */
	::Box<dim, T> & getProcessorBounds() {
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
	bool write(std::string output) const {
		//! subdomains_X.vtk domain for the local processor (X) as union of sub-domain
		VTKWriter<openfpm::vector<::SpaceBox<dim, T>>, VECTOR_BOX> vtk_box1;
		vtk_box1.add(sub_domains);
		vtk_box1.write(output + std::string("subdomains_") + std::to_string(v_cl.getProcessUnitID()) + std::string(".vtk"));

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
	bool check_consistency() {
		if (ie_loc_ghost<dim, T>::check_consistency(getNLocalHyperCube()) == false)
			return false;

		return true;
	}

	void debugPrint() {
		std::cout << "Subdomains\n";
		for (size_t p = 0; p < sub_domains.size(); p++) {
			std::cout << ::SpaceBox<dim, T>(sub_domains.get(p)).toString() << "\n";
		}

		std::cout << "External ghost box\n";

		for (size_t p = 0; p<nn_prcs < dim, T>::getNNProcessors(); p++) {
			for (size_t i = 0; i<ie_ghost < dim, T>::getProcessorNEGhost(p); i++) {
				std::cout << ie_ghost<dim, T>::getProcessorEGhostBox(p, i).toString() << "   prc=" << nn_prcs<dim, T>::IDtoProc(p)
						<< "   id=" << ie_ghost<dim, T>::getProcessorEGhostId(p, i) << "\n";
			}
		}

		std::cout << "Internal ghost box\n";

		for (size_t p = 0; p<nn_prcs < dim, T>::getNNProcessors(); p++) {
			for (size_t i = 0; i<ie_ghost < dim, T>::getProcessorNIGhost(p); i++) {
				std::cout << ie_ghost<dim, T>::getProcessorIGhostBox(p, i).toString() << "   prc=" << nn_prcs<dim, T>::IDtoProc(p)
						<< "   id=" << ie_ghost<dim, T>::getProcessorIGhostId(p, i) << "\n";
			}
		}
	}

};

#endif
