/*
 * P4EstDistribution.hpp
 *
 *  Created on: Mar 11, 2018
 *      Author: i-bird
 */

#ifndef SRC_DECOMPOSITION_DISTRIBUTION_P4ESTDISTRIBUTION_HPP_
#define SRC_DECOMPOSITION_DISTRIBUTION_P4ESTDISTRIBUTION_HPP_

#define P4EST_DISTRIBUTION_ERROR 100003

#include "VCluster/VCluster.hpp"
#include "Graph/dist_map_graph.hpp"
#include "SubdomainGraphNodes.hpp"

/*! \brief Class that distribute sub-sub-domains across processors using p4est library
 *
 * p4est library distribution. In general p4est produce automatically a distribution
 * using a space filling curve, so the standard p4est distribution does not do anything.
 *
 * ### Initialize a Cartesian graph and decompose
 * \snippet Distribution_unit_tests.hpp Initialize a ParMetis Cartesian graph and decompose
 *
 * ### Refine the decomposition
 * \snippet Distribution_unit_tests.hpp refine with parmetis the decomposition
 *
 */
template<unsigned int dim, typename T>
class P4estDistribution
{
	//! Is distributed
	bool is_distributed = false;

	//! Vcluster
	Vcluster & v_cl;

	//! Structure that store the initial starting grid
	grid_sm<dim, void> gr;

	//! rectangular domain to decompose
	Box<dim, T> domain;

	//! Global sub-sub-domain graph
	DistGraph_CSR<nm_v, nm_e> g;

	//! Flag to check if weights are used on vertices
	bool verticesGotWeights = false;


public:

	/*! Constructor for the ParMetis class
	 *
	 * \param v_cl Vcluster to use as communication object in this class
	 */
	P4estDistribution(Vcluster & v_cl)
	:is_distributed(false),v_cl(v_cl)
	{
	}

	/*! Copy constructor
	 *
	 * \param pm Distribution to copy
	 *
	 */
	P4estDistribution(const P4estDistribution<dim,T> & pm)
	:v_cl(pm.v_cl)
	{
		this->operator=(pm);
	}

	/*! \brief Initialize the initial distribution graph
	 *
	 * /param grid Grid
	 * /param dom Domain
	 */
	void createCartGraph(grid_sm<dim, void> & grid, Box<dim, T> dom)
	{
		//! Set grid and domain
		gr = grid;
		domain = dom;

		//! Create sub graph
		DistGraphFactory<dim, DistGraph_CSR<nm_v, nm_e>> dist_g_factory;
		g = dist_g_factory.template construct<NO_EDGE, T, dim - 1, 0>(gr.getSize(), domain);

		if (dim == 2)
			for (size_t i = 0; i < g.getNVertex(); i++)
				g.vertex(i).template get<nm_v::x>()[2] = 0;
	}

	/*! Copy constructor
	 *
	 * \param pm Distribution to copy
	 *
	 */
	P4estDistribution(P4estDistribution<dim,T> && pm)
	:v_cl(pm.v_cl)
	{
		this->operator=(pm);
	}

	/*! \brief Get the current graph (main)
	 *
	 */
	DistGraph_CSR<nm_v, nm_e> & getGraph()
	{
		return g;
	}

	/*! \brief Create the decomposition
	 *
	 */
	void decompose()
	{


		is_distributed = true;
	}

	/*! \brief Refine current decomposition
	 *
	 * It makes a refinement of the current decomposition using Parmetis function RefineKWay
	 * After that it also does the remapping of the graph
	 *
	 */
	void refine()
	{

	}

	/*! \brief Redecompose current decomposition
	 *
	 * It makes a redecomposition using Parmetis taking into consideration
	 * also migration cost
	 *
	 */
	void redecompose()
	{

	}

	/*! \brief Compute the unbalance of the processor compared to the optimal balance
	 *
	 * \return the unbalance from the optimal one 0.01 mean 1%
	 */
	float getUnbalance()
	{
		return 0.0;
	}

	/*! \brief function that return the position of the vertex in the space
	 *
	 * \param id vertex id
	 * \param pos vector that will contain x, y, z
	 *
	 */
	void getSubSubDomainPosition(size_t id, T (&pos)[dim])
	{
#ifdef SE_CLASS1
		if (id >= gp.getNVertex())
			std::cerr << __FILE__ << ":" << __LINE__ << "Such vertex doesn't exist (id = " << id << ", " << "total size = " << gp.getNVertex() << ")\n";
#endif

		// Copy the geometrical informations inside the pos vector
		pos[0] = g.vertex(id).template get<nm_v::x>()[0];
		pos[1] = g.vertex(id).template get<nm_v::x>()[1];
		if (dim == 3)
			pos[2] = g.vertex(id).template get<nm_v::x>()[2];
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

#ifdef SE_CLASS1
		if (id >= gp.getNVertex())
			std::cerr << __FILE__ << ":" << __LINE__ << "Such vertex doesn't exist (id = " << id << ", " << "total size = " << gp.getNVertex() << ")\n";
#endif

		// Update vertex in main graph
		g.vertex(id).template get<nm_v::computation>() = weight;
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
	size_t getSubSubDomainComputationCost(size_t id)
	{
#ifdef SE_CLASS1
		if (id >= gp.getNVertex())
			std::cerr << __FILE__ << ":" << __LINE__ << "Such vertex doesn't exist (id = " << id << ", " << "total size = " << gp.getNVertex() << ")\n";
#endif

		return g.vertex(id).template get<nm_v::computation>();
	}

	/*! \brief Compute the processor load counting the total weights of its vertices
	 *
	 * \return the computational load of the processor graph
	 */
	size_t getProcessorLoad()
	{
		size_t load = 0;

		return load;
	}

	/*! \brief Set migration cost of the vertex id
	 *
	 * \param id of the vertex to update
	 * \param migration cost of the migration
	 */
	void setMigrationCost(size_t id, size_t migration)
	{
#ifdef SE_CLASS1
		if (id >= gp.getNVertex())
			std::cerr << __FILE__ << ":" << __LINE__ << "Such vertex doesn't exist (id = " << id << ", " << "total size = " << gp.getNVertex() << ")\n";
#endif

		g.vertex(id).template get<nm_v::migration>() = migration;
	}

	/*! \brief Set communication cost of the edge id
	 *
	 * \param v_id Id of the source vertex of the edge
	 * \param e i child of the vertex
	 * \param communication Communication value
	 */
	void setCommunicationCost(size_t v_id, size_t e, size_t communication)
	{
#ifdef SE_CLASS1

		size_t e_id = v_id + e;

		if (e_id >= gp.getNEdge())
			std::cerr << "Such edge doesn't exist (id = " << e_id << ", " << "total size = " << gp.getNEdge() << ")\n";
#endif

		g.getChildEdge(v_id, e).template get<nm_e::communication>() = communication;
	}

	/*! \brief Returns total number of sub-sub-domains in the distribution graph
	 *
	 * \return the total number of sub-sub-domains
	 *
	 */
	size_t getNSubSubDomains() const
	{
		return g.getNVertex();
	}

	/*! \brief Return the total number of sub-sub-domains this processor own
	 *
	 * \return the total number of sub-sub-domains owned by this processor
	 *
	 */
	size_t getNOwnerSubSubDomains() const
	{
		return 0;
	}

	/*! \brief Return the global id of the owned sub-sub-domain
	 *
	 * \param id in the list of owned sub-sub-domains
	 *
	 * \return the global id
	 *
	 */
	size_t getOwnerSubSubDomain(size_t id) const
	{
		return 0;
	}

	/*! \brief Returns total number of neighbors of the sub-sub-domain id
	 *
	 * \param id id of the sub-sub-domain
	 *
	 * \return the number of neighborhood sub-sub-domains for each sub-domain
	 *
	 */
	size_t getNSubSubDomainNeighbors(size_t id)
	{
#ifdef SE_CLASS1
		if (id >= gp.getNVertex())
			std::cerr << __FILE__ << ":" << __LINE__ << "Such vertex doesn't exist (id = " << id << ", " << "total size = " << gp.getNVertex() << ")\n";
#endif

		return 0;
	}

	/*! \brief Print the current distribution and save it to VTK file
	 *
	 * \param file filename
	 *
	 */
	void write(const std::string & file)
	{
//		VTKWriter<Graph_CSR<nm_v, nm_e>, VTK_GRAPH> gv2(gp);
//		gv2.write(std::to_string(v_cl.getProcessUnitID()) + "_" + file + ".vtk");
	}

	const P4estDistribution<dim,T> & operator=(const P4estDistribution<dim,T> & dist)
	{
		is_distributed = dist.is_distributed;
		gr = dist.gr;
		domain = dist.domain;
		g = dist.g;
		verticesGotWeights = dist.verticesGotWeights;

		return *this;
	}

	const P4estDistribution<dim,T> & operator=(P4estDistribution<dim,T> && dist)
	{
		is_distributed = dist.is_distributed;
		v_cl = dist.v_cl;
		gr = dist.gr;
		domain = dist.domain;
		g.swap(dist.g);
		verticesGotWeights = dist.verticesGotWeights;

		return *this;
	}

	/*! \brief Get the decomposition counter
	 *
	 * \return the decomposition counter
	 *
	 */
	size_t get_ndec()
	{
		return 0;
	}

	/*! \brief Set the tolerance for each partition
	 *
	 * \param tol tolerance
	 *
	 */
	void setDistTol(double tol)
	{
	}

	/*! \brief Parmetis distribution is distribute sub-sub-domain on a regular grid
	 *
	 * \return true
	 *
	 */
	constexpr bool isRegularGrid()
	{
		return false;
	}

	/*! \brief Parmetis distribution is not for high scalability
	 *
	 * \return true
	 *
	 */
	constexpr bool isHighScal()
	{
		return true;
	}
};

#endif /* SRC_DECOMPOSITION_DISTRIBUTION_P4ESTDISTRIBUTION_HPP_ */
