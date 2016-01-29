/*
 * MetisDistribution.hpp
 *
 *  Created on: Nov 19, 2015
 *      Author: Antonio Leo
 */

#ifndef SRC_DECOMPOSITION_METISDISTRIBUTION_HPP_
#define SRC_DECOMPOSITION_METISDISTRIBUTION_HPP_

#include "SubdomainGraphNodes.hpp"
#include "metis_util.hpp"

template<unsigned int dim, typename T, template<unsigned int, typename > class Domain = Box>
class MetisDistribution
{
	//! Vcluster
	Vcluster & v_cl;

	//! Structure that store the cartesian grid information
	grid_sm<dim, void> gr;

	//! rectangular domain to decompose
	Domain<dim, T> domain;

	//! Global sub-sub-domain graph
	Graph_CSR<nm_v, nm_e> gp;

	//! Flag to check if weights are used on vertices
	bool verticesGotWeights = false;

public:

	//! constructor
	MetisDistribution(Vcluster & v_cl) :
			v_cl(v_cl)
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

	/*! \brief Create first decomposition, it divides the graph in slices and give each slice to a processor
	 *
	 */
	void decompose()
	{
		Metis<Graph_CSR<nm_v, nm_e>> met(gp, v_cl.getProcessingUnits(), verticesGotWeights);

		// decompose
		met.decompose<nm_v::proc_id>();
	}

	/*! \brief Refine current decomposition (NOT AVAILABLE on Metis)
	 *
	 * It has no function
	 */
	void refine()
	{
	}

	/*! \brief function that return the position of the vertex in the space
	 *
	 * \param id vertex id
	 * \param pos vector that will contain x, y, z
	 *
	 */
	void getVertexPosition(size_t id, T (&pos)[dim])
	{
		if (id >= gp.getNVertex())
			std::cerr << "Such vertex doesn't exist (id = " << id << ", " << "total size = " << gp.getNVertex() << ")\n";

		// Copy the geometrical informations inside the pos vector
		pos[0] = gp.vertex(id).template get<nm_v::x>()[0];
		pos[1] = gp.vertex(id).template get<nm_v::x>()[1];
		if (dim == 3)
			pos[2] = gp.vertex(id).template get<nm_v::x>()[2];
	}

	/*! \brief function that set the weight of the vertex
	 *
	 * \param id vertex id
	 * \param wieght to give to the vertex
	 *
	 */
	void setVertexWeight(size_t id, size_t weight)
	{
		if (!verticesGotWeights)
			verticesGotWeights = true;

		if (id >= gp.getNVertex())
			std::cerr << "Such vertex doesn't exist (id = " << id << ", " << "total size = " << gp.getNVertex() << ")\n";

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
			std::cerr << "Such vertex doesn't exist (id = " << id << ", " << "total size = " << gp.getNVertex() << ")\n";

		return gp.vertex(id).template get<nm_v::computation>();
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

	/*! \brief Print current graph and save it to file with name test_graph_[id]
	 *
	 * \param id to attach to the filename
	 *
	 */
	void printCurrentDecomposition(int id)
	{
		VTKWriter<Graph_CSR<nm_v, nm_e>, GRAPH> gv2(gp);
		gv2.write("test_graph_" + std::to_string(id) + ".vtk");

	}

	/*! \brief Compute the unbalance value
	 *
	 * \return the unbalance value
	 */
	float getUnbalance()
	{
		long min, max, sum;
		std::vector<long> loads(v_cl.getProcessingUnits());

		for (size_t i = 0; i < loads.size(); i++)
			loads[i] = 0;

		for (size_t i = 0; i < gp.getNVertex(); i++)
			loads[gp.vertex(i).template get<nm_v::proc_id>()]++;

		max = *std::max_element(loads.begin(), loads.end());
		min = *std::min_element(loads.begin(), loads.end());
		sum = std::accumulate(loads.begin(), loads.end(), 0);

		float unbalance = ((float) (max - min)) / (float) sum;

		return unbalance;
	}

	/*! \brief Compute the processor load counting the total weights of its vertices
	 *
	 * \return the computational load of the processor graph
	 */
	size_t getProcessorLoad()
	{
		size_t load = 0;

		for (size_t i = 0; i < gp.getNVertex(); i++)
		{
			if (gp.vertex(i).template get<nm_v::proc_id>() == v_cl.getProcessUnitID())
				load += gp.vertex(i).template get<nm_v::computation>();
		}

		return load;
	}
};

#endif /* SRC_DECOMPOSITION_METISDISTRIBUTION_HPP_ */
