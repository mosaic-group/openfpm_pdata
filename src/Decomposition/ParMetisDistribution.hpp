/*
 * ParMetisDistribution.hpp
 *
 *  Created on: Nov 19, 2015
 *      Author: i-bird
 */

#ifndef SRC_DECOMPOSITION_PARMETISDISTRIBUTION_HPP_
#define SRC_DECOMPOSITION_PARMETISDISTRIBUTION_HPP_

template<unsigned int dim, typename T>
class ParMetisDistribution
{
	//! Vcluster
	Vcluster & v_cl;

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

	//! Number of moved vertices in all iterations
	size_t g_moved = 0;

	//! Max number of moved vertices in all iterations
	size_t m_moved = 0;

	//! Wn for SAR heuristic
	float w_n = 0;

	//! Computation cost for SAR heuristic
	float c_c = 5;

	//! Number of time-steps since the previous DLB
	size_t n_ts = 1;

	//! Idle time accumulated so far, needed for SAR heuristic
	openfpm::vector<float> i_times;

	// Vector to collect all timings
	openfpm::vector<long> times;

	/* \brief fill the graph of the processor with the first decomposition (linear)
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

	/* \brief Update main graph ad subgraph with the partition in partitions param and renumber graphs
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
		openfpm::vector < idx_t > n_vtxdist(Np + 1);
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
	ParMetisDistribution(Vcluster & v_cl)
	:v_cl(v_cl)
	{}

	/* ! \brief Refine current decomposition
	 *
	 * It makes a refinement of the current decomposition using Parmetis function RefineKWay
	 * After that it also does the remapping of the graph
	 *
	 */
	void refine()
	{
		size_t Np = v_cl.getProcessingUnits();
		size_t p_id = v_cl.getProcessUnitID();

		//0.01 and 1 must be given TODO
		computeCommunicationAndMigrationCosts(0.01, n_ts);

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

		openfpm::vector < size_t > prc;
		openfpm::vector < size_t > sz;
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

	/*! \brief Function that gather times informations and decides if a rebalance is needed it uses the SAR heuristic
	 *
	 * \param t
	 *
	 */
	bool balanceNeeded(long t)
	{
		float t_max = 0, t_avg = 0;

		// Exchange time informations through processors
		v_cl.allGather(t, times);
		v_cl.execute();

		t_max = *(std::max_element(std::begin(times), std::end(times)));
		//if(v_cl.getProcessUnitID())
			//std::cout << "tmax: " << t_max << "\n";

		t_avg = std::accumulate(times.begin(), times.end(), 0) / v_cl.getProcessingUnits();
		//std::cout << "tavg: " << t_avg << "\n";

		// add idle time to vector
		i_times.add(t_max - t_avg);

		// Compute Wn
		double it_sum = *(std::max_element(std::begin(i_times), std::end(i_times)));
		float nw_n = (it_sum + c_c) / n_ts;

		if(nw_n > w_n){
			i_times.clear();
			n_ts = 1;
			w_n = nw_n;
			return true;
		}else{
			++n_ts;
			w_n = nw_n;
			return false;
		}
	}

	/* ! \brief function that return the position of the vertex in the space
	 *
	 * \param id vertex id
	 * \param pos vector that will contain x, y, z
	 *
	 */
	void getVertexPosition(size_t id, openfpm::vector<real_t> &pos)
	{
		pos.get(0) = gp.vertex(id).template get<nm_v::x>();
		pos.get(1) = gp.vertex(id).template get<nm_v::y>();

		if (dim == 3)
			pos.get(2) = gp.vertex(id).template get<nm_v::z>();
	}

	/* ! \brief function that set the weight of the vertex
	 *
	 * \param id vertex id
	 *
	 * \return vector with x, y, z
	 *
	 */
	void setVertexWeight(size_t id, size_t weight)
	{
		gp.vertex(id).template get<nm_v::computation>() = weight;
	}

	/* ! \brief return number of moved vertices in all iterations so far
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

	/* ! \brief return number of moved vertices in all iterations so far
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
};

#endif /* SRC_DECOMPOSITION_PARMETISDISTRIBUTION_HPP_ */
