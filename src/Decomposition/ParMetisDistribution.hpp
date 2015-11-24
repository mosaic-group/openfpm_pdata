/*
 * ParMetisDistribution.hpp
 *
 *  Created on: Nov 19, 2015
 *      Author: Antonio Leo
 */

#include "SubdomainGraphNodes.hpp"
#include "parmetis_util.hpp"

#ifndef SRC_DECOMPOSITION_PARMETISDISTRIBUTION_HPP_
#define SRC_DECOMPOSITION_PARMETISDISTRIBUTION_HPP_

template<unsigned int dim, typename T, template<unsigned int, typename > class Domain = Box>
class ParMetisDistribution
{
	//! Vcluster
	Vcluster & v_cl;

	//! Structure that store the cartesian grid information
	grid_sm<dim, void> gr;

	//! rectangular domain to decompose
	Domain<dim, T> domain;

	//! Global sub-sub-domain graph
	Graph_CSR<nm_v, nm_e> gp;

	//! Processor sub-sub-domain graph
	Graph_CSR<nm_v, nm_e> sub_g;

	//! Convert the graph to parmetis format
	Parmetis<Graph_CSR<nm_v, nm_e>> parmetis_graph;

	//! Init vtxdist needed for Parmetis
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

	/*! \brief fill the graph of the processor with the first decomposition (linear)
	 * Put vertices into processor graph (different for each processor)
	 *
	 * \param sub_g sub graph to fill
	 * \param gp mai graph, source for the vertices
	 * \param vtxdist array with the distribution of vertices through processors
	 * \param proc_id rank of the processor
	 * \param Np total number of processors
	 */
	void fillSubGraph()
	{

		int Np = v_cl.getProcessingUnits();
		int p_id = v_cl.getProcessUnitID();

		for (size_t j = vtxdist.get(p_id), local_j = 0; j < vtxdist.get(p_id + 1); j++, local_j++)
		{
			// Add vertex
			nm_v pv = gp.vertexById(j);
			sub_g.addVertex(pv);

			// Add edges of vertex
			for (size_t s = 0; s < gp.getNChilds(j); s++)
			{
				nm_e pe = gp.edge(j + s);
				sub_g.template addEdge<NoCheck>(local_j, gp.getChild(j, s), pe);
			}
		}

		// Just for output purpose
		if (p_id == 0)
		{
			for (int i = 0; i < Np; i++)
			{
				for (size_t j = vtxdist.get(i); j < vtxdist.get(i + 1); j++)
				{
					gp.vertexById(j).template get<nm_v::proc_id>() = i;
				}
			}
		}
	}

	/*! \brief Update main graph ad subgraph with the partition in partitions param and renumber graphs
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
	void updateGraphs()
	{

		int Np = v_cl.getProcessingUnits();
		int p_id = v_cl.getProcessUnitID();

		//stats info
		size_t moved = 0;

		// reset sub graph and local subgroph index
		int local_j = 0;
		sub_g.clear();

		// Init n_vtxdist to gather informations about the new decomposition
		openfpm::vector<idx_t> n_vtxdist(Np + 1);
		for (int i = 0; i <= Np; i++)
			n_vtxdist.get(i) = 0;

		// Update main graph with other partitions made by Parmetis in other processors and the local partition
		for (int i = 0; i < Np; i++)
		{

			int ndata = partitions.get(i).size();

			// Update the main graph with received informations
			for (int k = 0, l = vtxdist.get(i); k < ndata && l < vtxdist.get(i + 1); k++, l++)
			{

				// Create new n_vtxdist (1) (just count processors vertices)
				n_vtxdist.get(partitions.get(i).get(k) + 1)++;

				if
(				gp.vertexById(l).template get<nm_v::proc_id>() != partitions.get(i).get(k))
				moved++;

				// Update proc id in the vertex
				gp.vertexById(l).template get<nm_v::proc_id>() = partitions.get(i).get(k);
				gp.vertex(l).template get<nm_v::global_id>() = l;

				// Add vertex to temporary structure of distribution (needed to update main graph)
				v_per_proc.get(partitions.get(i).get(k)).add(gp.getVertexOldId(l));

				// Add vertices belonging to this processor in sub graph
				if (partitions.get(i).get(k) == p_id)
				{

					nm_v pv = gp.vertexById(l);
					sub_g.addVertex(pv);

					// Add edges of vertex
					for (size_t s = 0; s < gp.getNChildsByVertexId(l); s++)
					{
						nm_e pe = gp.edge(l + s);
						sub_g.template addEdge<NoCheck>(local_j, gp.getChildByVertexId(l, s), pe);
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
		for (int i = 0; i <= Np; i++)
		{
			vtxdist.get(i) = n_vtxdist.get(i);
		}

		// Renumbering subgraph
		sub_g.reset_map_ids();
		for (size_t j = vtxdist.get(p_id), i = 0; j < vtxdist.get(p_id + 1); j++, i++)
		{
			sub_g.set_map_ids(j, sub_g.vertex(i).template get<nm_v::global_id>());
			sub_g.vertex(i).template get<nm_v::id>() = j;
		}

		// Renumbering main graph
		for (size_t p = 0; p < Np; p++)
		{
			for (size_t j = vtxdist.get(p), i = 0; j < vtxdist.get(p + 1); j++, i++)
			{
				gp.set_map_ids(j, v_per_proc.get(p).get(i));
				gp.vertex(v_per_proc.get(p).get(i)).template get<nm_v::id>() = j;
			}
		}

		g_moved += moved;

		if (moved > m_moved)
			m_moved = moved;

	}

	static void * message_receive(size_t msg_i, size_t total_msg, size_t total_p, size_t i, size_t ri, void * ptr)
	{
		openfpm::vector < openfpm::vector < idx_t >> *v = static_cast<openfpm::vector<openfpm::vector<idx_t>> *>(ptr);

		v->get(i).resize(msg_i / sizeof(idx_t));

		return &(v->get(i).get(0));
	}

public:

	//! constructor
	ParMetisDistribution(Vcluster & v_cl) :
			v_cl(v_cl), parmetis_graph(v_cl, v_cl.getProcessingUnits()), vtxdist(v_cl.getProcessingUnits() + 1), partitions(
					v_cl.getProcessingUnits()), v_per_proc(v_cl.getProcessingUnits())

	{
	}

	/*! \brief Initialize the distribution graph
	 *
	 */
	void init(grid_sm<dim, void> & grid, Domain<dim, T> dom)
	{
		// Set grid and domain
		gr = grid;
		domain = dom;

		// Create a cartesian grid graph
		CartesianGraphFactory<dim, Graph_CSR<nm_v, nm_e>> g_factory_part;
		gp = g_factory_part.template construct<NO_EDGE, nm_v::id, T, dim - 1, 0, 1, 2>(gr.getSize(), domain);
		gp.init_map_ids();

		// Init to 0.0 axis z (to fix in graphFactory)
		if (dim < 3)
		{
			for (size_t i = 0; i < gp.getNVertex(); i++)
			{
				gp.vertex(i).template get<nm_v::z>() = 0.0;
			}
		}
	}

	/*! \brief Get the current graph (main)
	 *
	 */
	Graph_CSR<nm_v, nm_e> & getGraph()
	{
		return gp;
	}

	/*! \brief Create first decomposition, it divides the graph in slices and give each slice to a processor
	 *
	 */
	void decompose()
	{
		//! Get the processor id
		size_t p_id = v_cl.getProcessUnitID();

		//! Get the number of processing units
		size_t Np = v_cl.getProcessingUnits();

		//! Division of vertices in Np graphs
		//! Put (div+1) vertices in mod graphs
		//! Put div vertices in the rest of the graphs
		size_t mod_v = gp.getNVertex() % Np;
		size_t div_v = gp.getNVertex() / Np;

		for (int i = 0; i <= Np; i++)
		{
			if (i < mod_v)
				vtxdist.get(i) = (div_v + 1) * (i);
			else
				vtxdist.get(i) = (div_v) * (i) + mod_v;
		}

		//TODO transform in factory

		//! Put vertices into processor graph (different for each processor)
		fillSubGraph();

		parmetis_graph.initSubGraph(sub_g);

		//! Decompose
		parmetis_graph.decompose<nm_v::proc_id>(vtxdist, sub_g);

		//! Get result partition for this processors
		idx_t *partition = parmetis_graph.getPartition();

		//! Prepare vector of arrays to contain all partitions
		partitions.get(p_id).resize(sub_g.getNVertex());
		std::copy(partition, partition + sub_g.getNVertex(), &partitions.get(p_id).get(0));

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

		// Update graphs with the new distributions
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
		for (int i = 0; i < Np; ++i)
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

	/*! \brief function that return the position of the vertex in the space
	 *
	 * \param id vertex id
	 * \param pos vector that will contain x, y, z
	 *
	 */
	void getVertexPosition(size_t id, openfpm::vector<real_t> &pos)
	{
		if (id >= gp.getNVertex())
			std::cerr << "Such vertex doesn't exist (id = " << id << ", " << "total size = " << gp.getNVertex()
					<< ")\n";

		pos.get(0) = gp.vertex(id).template get<nm_v::x>();
		pos.get(1) = gp.vertex(id).template get<nm_v::y>();

		if (dim == 3)
			pos.get(2) = gp.vertex(id).template get<nm_v::z>();
	}

	/*! \brief function that set the weight of the vertex
	 *
	 * \param id vertex id
	 * \param wieght to give to the vertex
	 *
	 */
	void setVertexWeight(size_t id, size_t weight)
	{
		if(!verticesGotWeights)
			verticesGotWeights = true;

		if (id >= gp.getNVertex())
			std::cerr << "Such vertex doesn't exist (id = " << id << ", " << "total size = " << gp.getNVertex()
					<< ")\n";

		gp.vertex(id).template get<nm_v::computation>() = weight;
	}

	/*! \brief Checks if weights are used on the vertices
	 *
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
			std::cerr << "Such vertex doesn't exist (id = " << id << ", " << "total size = " << gp.getNVertex()
					<< ")\n";

		return gp.vertex(id).template get<nm_v::computation>();
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
			std::cerr << "Such vertex doesn't exist (id = " << id << ", " << "total size = " << gp.getNVertex()
					<< ")\n";

		gp.vertex(id).template get<nm_v::migration>() = migration;
	}

	/*! \brief Set communication cost of the edge id
	 *
	 */
	void setCommunicationCost(size_t id, size_t communication)
	{
		if (id >= gp.getNEdge())
			std::cerr << "Such edge doesn't exist (id = " << id << ", " << "total size = " << gp.getNEdge() << ")\n";

		gp.edge(id).template get<nm_e::communication>() = communication;
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
			std::cerr << "Such vertex doesn't exist (id = " << id << ", " << "total size = " << gp.getNVertex()
					<< ")\n";

		return gp.getNChilds(id);
	}

	/* \brief Print current graph and save it to file with name test_graph_[id]
	 *
	 * \param id to attach to the filename
	 *
	 */
	void printCurrentDecomposition(int id)
	{
		VTKWriter<Graph_CSR<nm_v, nm_e>, GRAPH> gv2(gp);
		gv2.write("test_graph_" + std::to_string(id) + ".vtk");

	}
};

#endif /* SRC_DECOMPOSITION_PARMETISDISTRIBUTION_HPP_ */
