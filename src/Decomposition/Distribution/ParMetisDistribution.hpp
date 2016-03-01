/*
 * ParMetisDistribution.hpp
 *
 *  Created on: Nov 19, 2015
 *      Author: Antonio Leo
 */

#include "SubdomainGraphNodes.hpp"
#include "parmetis_util.hpp"
#include "Graph/dist_map_graph.hpp"
#ifndef SRC_DECOMPOSITION_PARMETISDISTRIBUTION_HPP_
#define SRC_DECOMPOSITION_PARMETISDISTRIBUTION_HPP_

#define PARMETIS_DISTRIBUTION_ERROR 100002

/*! \brief Class that distribute sub-sub-domains across processors using ParMetis Library
 *
 * Given a graph and setting Computational cost, Communication cost (on the edge) and
 * Migration cost or total Communication costs, it produce the optimal balanced distribution
 *
 * In addition to Metis it provide the functionality to refine the previously computed
 * decomposition
 *
 * ### Initialize a Cartesian graph and decompose
 * \snippet Distribution_unit_tests.hpp Initialize a ParMetis Cartesian graph and decompose
 *
 * ### Refine the decomposition
 * \snippet Distribution_unit_tests.hpp refine with parmetis the decomposition
 *
 */

template<unsigned int dim, typename T>
class ParMetisDistribution
{
	//! Vcluster
	Vcluster & v_cl;

	//! Structure that store the cartesian grid information
	grid_sm<dim, void> gr;

	//! rectangular domain to decompose
	Box<dim, T> domain;

	//! Global sub-sub-domain graph
	Graph_CSR<nm_v, nm_e> gp;

	//! Processor sub-sub-domain graph
	//
	// The subgraph has an internal map that permit to access the vertices
	// with the global-id (vtxdist) getVertexByMapId
	//
	//
	Graph_CSR<nm_v, nm_e> sub_g;

	//! Convert the graph to parmetis format
	Parmetis<Graph_CSR<nm_v, nm_e>> parmetis_graph;

	//! Init vtxdist needed for Parmetis
	//
	// vtxdist is a common array across processor, it indicate how
	// vertex are distributed across processors
	//
	// Example we have 3 processors
	//
	// processor 0 has 3 vertices
	// processor 1 has 5 vertices
	// processor 2 has 4 vertices
	//
	// vtxdist contain, 0,3,8,12
	//
	// vtx dist is the unique global-id of the vertices
	//
	openfpm::vector<idx_t> vtxdist;

	//! partitions
	openfpm::vector<openfpm::vector<idx_t>> partitions;

	//! Init data structure to keep trace of new vertices distribution in processors (needed to update main graph)
	openfpm::vector<openfpm::vector<size_t>> v_per_proc;

	//! Number of moved vertices in all iterations
	size_t g_moved = 0;

	//! Max number of moved vertices in all iterations
	size_t m_moved = 0;

	//! Flag to check if weights are used on vertices
	bool verticesGotWeights = false;

	//! Flag that indicate if we are doing a test (In general it fix the seed)
	bool testing = false;

	/*! \brief Fill the graph of the processor with the first decomposition (linear)
	 *
	 * Each processor construct his part of graph
	 *
	 */
	void fillSubGraph()
	{

		int Np = v_cl.getProcessingUnits();
		int p_id = v_cl.getProcessUnitID();

		for (size_t j = vtxdist.get(p_id), local_j = 0; j < vtxdist.get(p_id + 1); j++, local_j++)
		{
			// Add vertex

			nm_v pv = gp.vertexByMapId(j);
			sub_g.addVertex(pv, gp.vertexByMapId(j).template get<nm_v::global_id>());

			// Add edges of vertex
			for (size_t s = 0; s < gp.getNChilds(j); s++)
			{
				if (gp.vertex(gp.getChild(j, s)).template get<nm_v::proc_id>() != v_cl.getProcessUnitID())
					gp.vertex(gp.getChild(j, s)).template get<nm_v::fake_v>() = 1;
				else
					gp.vertex(gp.getChild(j, s)).template get<nm_v::fake_v>() = 0;

				// Add Edge
				nm_e pe = gp.edge(j + s);
				sub_g.template addEdge<NoCheck>(local_j, gp.getChild(j, s), pe);
			}
		}

		// Just for output purpose, we set the processor id in the big graph
		if (p_id == 0)
		{
			for (int i = 0; i < Np; i++)
			{
				for (size_t j = vtxdist.get(i); j < vtxdist.get(i + 1); j++)
				{
					// we call the vertex of the sub-graph from its global id
					gp.vertexByMapId(j).template get<nm_v::proc_id>() = i;
				}
			}
		}
	}

	/*! \brief Update main graph ad subgraph with the received data of the partitions from the other processors
	 *
	 */
	void updateGraphs()
	{

		size_t Np = v_cl.getProcessingUnits();
		size_t p_id = v_cl.getProcessUnitID();

		//stats info
		size_t moved = 0;

		// reset sub graph and local subgraph index
		size_t local_j = 0;
		sub_g.clear();

		// Init n_vtxdist to gather informations about the new decomposition
		openfpm::vector<idx_t> n_vtxdist(Np + 1);
		for (size_t i = 0; i <= Np; i++)
			n_vtxdist.get(i) = 0;

		// Update the main graph with received data from processor i
		for (size_t i = 0; i < Np; i++)
		{
			size_t ndata = partitions.get(i).size();

			// Update the main graph with the received informations
			for (size_t k = 0, l = (size_t)vtxdist.get(i); k < ndata && l < (size_t)vtxdist.get(i + 1); k++, l++)
			{
				// Create new n_vtxdist (just count processors vertices)
				n_vtxdist.get(partitions.get(i).get(k) + 1)++;

				// for statistics
				if (gp.vertexByMapId(l).template get<nm_v::proc_id>() != (size_t)partitions.get(i).get(k))
					moved++;

				// Update proc id in the vertex (using the old map)
				gp.vertexByMapId(l).template get<nm_v::proc_id>() = partitions.get(i).get(k);

				// Add vertex to temporary structure of distribution (needed to update main graph)
				v_per_proc.get(partitions.get(i).get(k)).add(gp.getVertexGlobalId(l));

				// Add vertices belonging to this processor in sub graph
				if (partitions.get(i).get(k) == p_id)
				{

					nm_v pv = gp.vertexByMapId(l);
					sub_g.addVertex(pv, pv.template get<nm_v::global_id>());

					// Add edges of vertex
					for (size_t s = 0; s < gp.getNChildsByMapId(l); s++)
					{
						if (gp.vertex(gp.getChildByVertexId(l, s)).template get<nm_v::proc_id>() != v_cl.getProcessUnitID())
							gp.vertex(gp.getChildByVertexId(l, s)).template get<nm_v::fake_v>() = 1;
						else
							gp.vertex(gp.getChildByVertexId(l, s)).template get<nm_v::fake_v>() = 0;

						nm_e pe = gp.edge(l + s);
						sub_g.template addEdge<NoCheck>(local_j, gp.getChildByVertexId(l, s), pe);
					}

					local_j++;
				}

			}
		}

		// Create new n_vtxdist (accumulate the counters)
		for (size_t i = 2; i <= Np; i++)
			n_vtxdist.get(i) += n_vtxdist.get(i - 1);

		// Copy the new decomposition in the main vtxdist
		for (size_t i = 0; i <= Np; i++)
			vtxdist.get(i) = n_vtxdist.get(i);

		// Renumbering subgraph
		sub_g.resetLocalToGlobalMap();
		for (size_t j = (size_t)vtxdist.get(p_id), i = 0; j < (size_t)vtxdist.get(p_id + 1); j++, i++)
		{
			sub_g.setMapId(j, sub_g.vertex(i).template get<nm_v::global_id>());
			sub_g.template vertex_p<nm_v::id>(i) = j;
		}

		// Renumber the main graph and re-create the map
		for (size_t p = 0; p < (size_t)Np; p++)
		{
			for (size_t j = (size_t)vtxdist.get(p), i = 0; j < (size_t)vtxdist.get(p + 1); j++, i++)
			{
				gp.setMapId(j, v_per_proc.get(p).get(i));
				gp.template vertex_p<nm_v::id>(v_per_proc.get(p).get(i)) = j;
			}
		}

		g_moved += moved;

		if (moved > m_moved)
			m_moved = moved;

	}

	void createMapsFromGlobalGraph(openfpm::vector<size_t> & vtxdist)
	{
/*		openfpm::vector<size_t> cnt_np;

		for (size_t i = 0 ; i < gp.getNVertex() ; i++)
		{
			cnt_np(gp.template vertex<nm_v::proc_id>)++;

			gp.setMapId()
		}*/
	}

	/*! \brief Callback of the sendrecv to set the size of the array received
	 *
	 * \param msg_i Index of the message
	 * \param total_msg Total numeber of messages
	 * \param total_p Total number of processors to comunicate with
	 * \param i Processor id
	 * \param ri Request id
	 * \param ptr Void pointer parameter for additional data to pass to the call-back
	 */
	static void * message_receive(size_t msg_i, size_t total_msg, size_t total_p, size_t i, size_t ri, void * ptr)
	{
		openfpm::vector < openfpm::vector < idx_t >> *v = static_cast<openfpm::vector<openfpm::vector<idx_t>> *>(ptr);

		v->get(i).resize(msg_i / sizeof(idx_t));

		return &(v->get(i).get(0));
	}

public:

	/*! Constructor for the ParMetis class
	 *
	 * \param v_cl Vcluster to use as communication object in this class
	 */
	ParMetisDistribution(Vcluster & v_cl)
	:v_cl(v_cl), parmetis_graph(v_cl, v_cl.getProcessingUnits()), vtxdist(v_cl.getProcessingUnits() + 1), partitions(v_cl.getProcessingUnits()), v_per_proc(v_cl.getProcessingUnits())
	{
	}

	/*! Copy constructor
	 *
	 * \param pm Distribution to copy
	 *
	 */
	ParMetisDistribution(const ParMetisDistribution<dim,T> & pm)
	:v_cl(pm.v_cl),parmetis_graph(v_cl, v_cl.getProcessingUnits())
	{
		this->operator=(pm);
	}

	/*! Copy constructor
	 *
	 * \param pm Distribution to copy
	 *
	 */
	ParMetisDistribution(ParMetisDistribution<dim,T> && pm)
	{
		this->operator=(pm);
	}

	/*! \brief Create the Cartesian graph
	 *
	 * \param grid info
	 * \param dom domain
	 */
	void createCartGraph(grid_sm<dim, void> & grid, Box<dim, T> dom)
	{
		size_t bc[dim];

		for (size_t i = 0 ; i < dim ; i++)
			bc[i] = NON_PERIODIC;

		// Set grid and domain
		gr = grid;
		domain = dom;

		// Create a cartesian grid graph
		CartesianGraphFactory<dim, Graph_CSR<nm_v, nm_e>> g_factory_part;
		gp = g_factory_part.template construct<NO_EDGE, nm_v::id, T, dim - 1, 0, 1, 2>(gr.getSize(), domain, bc);
		gp.initLocalToGlobalMap();

		// Create sub graphs
		DistCartesianGraphFactory<dim, Graph_CSR<nm_v, nm_e>> dist_g_factory;
		sub_g = dist_g_factory.template construct<NO_EDGE, nm_v::id, nm_v::global_id, nm_e::srcgid, nm_e::dstgid, T, dim - 1, 0, 1, 2>(gr.getSize(), domain, vtxdist);

		// Init to 0.0 axis z (to fix in graphFactory)
		if (dim < 3)
		{
			for (size_t i = 0; i < gp.getNVertex(); i++)
			{
				gp.vertex(i).template get<nm_v::x>()[2] = 0.0;
			}
		}
		for (size_t i = 0; i < gp.getNVertex(); i++)
		{
			gp.vertex(i).template get<nm_v::global_id>() = i;
		}

	}

	/*! \brief Get the current graph (main)
	 *
	 */
	Graph_CSR<nm_v, nm_e> & getGraph()
	{
		return gp;
	}

	/*! \brief Create the decomposition
	 *
	 */
	void decompose()
	{

		//! Get the processor id
		size_t p_id = v_cl.getProcessUnitID();

		//! Get the number of processing units
		size_t Np = v_cl.getProcessingUnits();

		parmetis_graph.initSubGraph(sub_g, verticesGotWeights);

		//! Decompose
		parmetis_graph.decompose<nm_v::proc_id>(vtxdist, sub_g);

		//! Get result partition for this processors
		idx_t *partition = parmetis_graph.getPartition();

		//! Prepare vector of arrays to contain all partitions
		partitions.get(p_id).resize(sub_g.getNVertex());
		std::copy(partition, partition + sub_g.getNVertex(), &partitions.get(p_id).get(0));

		// Communicate the local distribution to the other processors
		// to reconstruct individually the global graph
		openfpm::vector<size_t> prc;
		openfpm::vector<size_t> sz;
		openfpm::vector<void *> ptr;

		for (size_t i = 0; i < Np; i++)
		{
			if (i != v_cl.getProcessUnitID())
			{
				prc.add(i);
				sz.add(sub_g.getNVertex() * sizeof(idx_t));
				ptr.add(partitions.get(p_id).getPointer());
			}
		}

		v_cl.sendrecvMultipleMessagesNBX(prc.size(), &sz.get(0), &prc.get(0), &ptr.get(0), message_receive, &partitions,
		NONE);

		// Update graphs with the received data
		updateGraphs();

		// reset statistical variables, we only need it in refinement
		g_moved = 0;
		m_moved = 0;

	}

	/*! \brief Refine current decomposition
	 *
	 * It makes a refinement of the current decomposition using Parmetis function RefineKWay
	 * After that it also does the remapping of the graph
	 *
	 */
	void refine()
	{
		size_t Np = v_cl.getProcessingUnits();
		size_t p_id = v_cl.getProcessUnitID();

		// Reset parmetis graph and reconstruct it
		parmetis_graph.reset(gp, sub_g);

		// Refine
		parmetis_graph.refine<nm_v::proc_id>(vtxdist, sub_g);

		// Get result partition for this processor
		idx_t * partition = parmetis_graph.getPartition();

		partitions.get(p_id).resize(sub_g.getNVertex());
		std::copy(partition, partition + sub_g.getNVertex(), &partitions.get(p_id).get(0));

		// Reset data structure to keep trace of new vertices distribution in processors (needed to update main graph)
		for (size_t i = 0; i < Np; ++i)
		{
			v_per_proc.get(i).clear();
		}

		openfpm::vector<size_t> prc;
		openfpm::vector<size_t> sz;
		openfpm::vector<void *> ptr;

		for (size_t i = 0; i < Np; i++)
		{
			if (i != v_cl.getProcessUnitID())
			{
				partitions.get(i).clear();
				prc.add(i);
				sz.add(sub_g.getNVertex() * sizeof(idx_t));
				ptr.add(partitions.get(p_id).getPointer());
			}
		}

		// Exchange informations through processors
		v_cl.sendrecvMultipleMessagesNBX(prc.size(), &sz.get(0), &prc.get(0), &ptr.get(0), message_receive, &partitions,
		NONE);

		// Update graphs with the new distributions
		updateGraphs();
	}

	/*! \brief Compute the unbalance of the processor compared to the optimal balance
	 *
	 * \return the unbalance from the optimal one 0.01 mean 1%
	 */
	float getUnbalance()
	{
		long t_cost = 0;

		long min, max, sum;
		float unbalance;

		t_cost = getProcessorLoad();

		min = t_cost;
		max = t_cost;
		sum = t_cost;

		v_cl.min(min);
		v_cl.max(max);
		v_cl.sum(sum);
		v_cl.execute();

		unbalance = ((float) (max - min)) / (float) (sum / v_cl.getProcessingUnits());

		return unbalance * 100;
	}

	/*! \brief function that return the position of the vertex in the space
	 *
	 * \param id vertex id
	 * \param pos vector that will contain x, y, z
	 *
	 */
	void getSubSubDomainPosition(size_t id, T (&pos)[dim])
	{
		if (id >= gp.getNVertex())
			std::cerr << "Such vertex doesn't exist (id = " << id << ", " << "total size = " << gp.getNVertex() << ")\n";

		//TOACTIVATE when move to the distributed graph
		//Return the pos object only if the vertex is in this graph
		/*
		 if(sub_g.vertexIsInThisGraph(id)){
		 pos[0] = sub_g.vertex(id).template get<nm_v::x>()[0];
		 pos[1] = sub_g.vertex(id).template get<nm_v::x>()[1];
		 if (dim == 3)
		 pos[2] = sub_g.vertex(id).template get<nm_v::x>()[2];
		 }*/

		// Copy the geometrical informations inside the pos vector
		pos[0] = gp.vertex(id).template get<nm_v::x>()[0];
		pos[1] = gp.vertex(id).template get<nm_v::x>()[1];
		if (dim == 3)
			pos[2] = gp.vertex(id).template get<nm_v::x>()[2];
	}

	/*! \brief Function that set the weight of the vertex
	 *
	 * \param id vertex id
	 * \param weight to give to the vertex
	 *
	 */
	inline void setComputationCost(size_t id, size_t weight)
	{
		if (!verticesGotWeights)
			verticesGotWeights = true;

		if (id >= gp.getNVertex())
			std::cerr << "Such vertex doesn't exist (id = " << id << ", " << "total size = " << gp.getNVertex() << ")\n";

		// If the vertex is inside this processor update the value
		if (sub_g.vertexIsInThisGraph(id))
		{
			sub_g.getLocalVertexByGlobalId(id).template get<nm_v::computation>() = weight;
		}

		// Update vertex in main graph
		gp.vertex(id).template get<nm_v::computation>() = weight;
	}

	/*! \brief Checks if weights are used on the vertices
	 *
	 * \return true if weights are used in the decomposition
	 */
	bool weightsAreUsed()
	{
		return verticesGotWeights;
	}

	/*! \brief function that get the weight of the vertex
	 *
	 * \param id vertex id
	 *
	 */
	size_t getVertexWeight(size_t id)
	{
		if (id >= gp.getNVertex())
			std::cerr << "Such vertex doesn't exist (id = " << id << ", " << "total size = " << gp.getNVertex() << ")\n";

		return gp.vertex(id).template get<nm_v::computation>();
	}

	/*! \brief Compute the processor load counting the total weights of its vertices
	 *
	 * \return the computational load of the processor graph
	 */
	size_t getProcessorLoad()
	{
		size_t load = 0;

		for (size_t i = 0; i < sub_g.getNVertex(); i++)
		{
			load += sub_g.vertex(i).template get<nm_v::computation>();
		}
		//std::cout << v_cl.getProcessUnitID() << " weight " << load << " size " << sub_g.getNVertex() << "\n";
		return load;
	}

	/*! \brief return number of moved vertices in all iterations so far
	 *
	 * \param id vertex id
	 *
	 * \return vector with x, y, z
	 *
	 */
	size_t getTotalMovedV()
	{
		return g_moved;
	}

	/*! \brief return number of moved vertices in all iterations so far
	 *
	 * \param id vertex id
	 *
	 * \return vector with x, y, z
	 *
	 */
	size_t getMaxMovedV()
	{
		return m_moved;
	}

	/*! \brief Set migration cost of the vertex id
	 *
	 * \param id of the vertex to update
	 * \param migration cost of the migration
	 */
	void setMigrationCost(size_t id, size_t migration)
	{
		if (id >= gp.getNVertex())
			std::cerr << "Such vertex doesn't exist (id = " << id << ", " << "total size = " << gp.getNVertex() << ")\n";

		// If the vertex is inside this processor update the value
		if (sub_g.vertexIsInThisGraph(id))
		{
			sub_g.getLocalVertexByGlobalId(id).template get<nm_v::migration>() = migration;
		}

		gp.vertex(id).template get<nm_v::migration>() = migration;
	}

	/*! \brief Set communication cost of the edge id
	 *
	 * \param v_id Id of the source vertex of the edge
	 * \param e i child of the vertex
	 * \param communication Communication value
	 */
	void setCommunicationCost(size_t v_id, size_t e, size_t communication)
	{
		size_t e_id = v_id + e;

		if (e_id >= gp.getNEdge())
			std::cerr << "Such edge doesn't exist (id = " << e_id << ", " << "total size = " << gp.getNEdge() << ")\n";

		// If the vertex is inside this processor update the value
		if (sub_g.vertexIsInThisGraph(v_id))
		{
			// Get the local id of the vertex
			size_t local_id = sub_g.getLocalIdFromGlobalId(v_id);
			sub_g.getChildEdge(local_id, e).template get<nm_e::communication>() = communication;
		}

		gp.getChildEdge(v_id, e).template get<nm_e::communication>() = communication;
	}

	/*! \brief Returns total number of sub-sub-domains in the distribution graph
	 *
	 */
	size_t getNSubSubDomains()
	{
		return gp.getNVertex();
	}

	/*! \brief Returns total number of neighbors of the sub-sub-domain id
	 *
	 * \param i id of the sub-sub-domain
	 */
	size_t getNSubSubDomainNeighbors(size_t id)
	{
		if (id >= gp.getNVertex())
			std::cerr << "Such vertex doesn't exist (id = " << id << ", " << "total size = " << gp.getNVertex() << ")\n";

		return gp.getNChilds(id);
	}

	/*! \brief It set the Class on test mode
	 *
	 * At the moment it fix the seed to have reproducible results
	 *
	 */
	void onTest()
	{
		testing = true;
	}

	/*! \brief Print the current distribution and save it to VTK file
	 *
	 * \param file filename
	 *
	 */
	void write(const std::string & file)
	{
		if (v_cl.getProcessUnitID() == 0)
		{
			VTKWriter<Graph_CSR<nm_v, nm_e>, VTK_GRAPH> gv2(gp);
			gv2.write(file);
		}

	}

	const ParMetisDistribution<dim,T> & operator=(const ParMetisDistribution<dim,T> & dist)
	{
		v_cl = dist.v_cl;
		gr = dist.gr;
		domain = dist.domain;
		sub_g = dist.sub_g;
		gp = dist.gp;
		vtxdist = dist.vtxdist;
		partitions = dist.partitions;
		v_per_proc = dist.v_per_proc;
		verticesGotWeights = dist.verticesGotWeights;

		return *this;
	}

	const ParMetisDistribution<dim,T> & operator=(ParMetisDistribution<dim,T> && dist)
	{
		v_cl = dist.v_cl;
		gr = dist.gr;
		domain = dist.domain;
		sub_g.swap(dist.sub_g);
		gp.swap(dist.gp);
		vtxdist.swap(dist.vtxdist);
		partitions.swap(dist.partitions);
		v_per_proc.swap(dist.v_per_proc);
		verticesGotWeights = dist.verticesGotWeights;

		return *this;
	}
};

#endif /* SRC_DECOMPOSITION_PARMETISDISTRIBUTION_HPP_ */
