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

#define METIS_DISTRIBUTION_ERROR 100001

/*! \brief Class that distribute sub-sub-domains across processors using Metis Library
 *
 * Given a graph and setting Computational cost, Communication cost (on the edge) and
 * Migration cost or total Communication costs, it produce the optimal distribution
 *
 * ### Initialize a Cartesian graph and decompose
 * \snippet Distribution_unit_tests.hpp Initialize a Metis Cartesian graph and decompose
 *
 * ### Set Computation Communication and Migration cost
 * \snippet Distribution_unit_tests.hpp Decomposition Metis with weights
 *
 */

template<unsigned int dim, typename T>
class MetisDistribution
{
	//! Vcluster
	Vcluster & v_cl;

	//! Structure that store the cartesian grid information
	grid_sm<dim, void> gr;

	//! rectangular domain to decompose
	Box<dim, T> domain;

	//! Global sub-sub-domain graph
	Graph_CSR<nm_v, nm_e> gp;

	//! Flag to check if weights are used on vertices
	bool useWeights = false;

	//! Flag that indicate if we are doing a test (In general it fix the seed)
	bool testing = false;

	/*! \brief Check that the sub-sub-domain id exist
	 *
	 * \param id sub-sub-domain id
	 *
	 */
	inline void check_overflow(size_t id)
	{
#ifdef SE_CLASS1
		if (id >= gp.getNVertex())
		{
			std::cerr << "Error " << __FILE__ ":" << __LINE__ << " such sub-sub-domain doesn't exist (id = " << id << ", " << "total size = " << gp.getNVertex() << ")\n";
			ACTION_ON_ERROR(METIS_DISTRIBUTION_ERROR)
		}
#endif
	}

	/*! \brief Check that the sub-sub-domain id exist
	 *
	 * \param id sub-sub-domain id
	 *
	 */
	inline void check_overflowe(size_t id, size_t e)
	{
#ifdef SE_CLASS1
		if (e >= gp.getNChilds(id))
		{
			std::cerr << "Error " << __FILE__ ":" << __LINE__ << " for the sub-sub-domain " << id << " such neighborhood doesn't exist (e = " << e << ", " << "total size = " << gp.getNChilds(id) << ")\n";
			ACTION_ON_ERROR(METIS_DISTRIBUTION_ERROR)
		}
#endif
	}

public:

	static constexpr unsigned int computation = nm_v::computation;

	//! constructor
	MetisDistribution(Vcluster & v_cl) :
			v_cl(v_cl)
	{
#ifdef SE_CLASS2
			check_new(this,8,VECTOR_EVENT,1);
#endif
	}

	/*! \brief Copy constructor
	 *
	 *
	 */
	MetisDistribution(const MetisDistribution & mt)
	:v_cl(mt.v_cl)
	{
#ifdef SE_CLASS2
			check_valid(mt);
			check_new(this,8,VECTOR_EVENT,1);
#endif
		this->operator=(mt);
	}

	/*! \brief Copy constructor
	 *
	 *
	 */
	MetisDistribution(MetisDistribution && mt)
	{
#ifdef SE_CLASS2
			check_valid(mt);
			check_new(this,8,VECTOR_EVENT,1);
#endif
		this->operator=(mt);
	}

	/*! \brief Destructor
	 *
	 *
	 */
	~MetisDistribution()
	{
#ifdef SE_CLASS2
		check_delete(this);
#endif
	}


	/*! \brief create a Cartesian distribution graph
	 *
	 * \param grid grid info (sub-sub somains on each dimension)
	 * \param dom domain (domain where the sub-sub-domains are defined)
	 *
	 */
	void createCartGraph(grid_sm<dim, void> & grid, Box<dim, T> dom)
	{
#ifdef SE_CLASS2
			check_valid(this,8);
#endif
		// NON periodic boundary conditions
		size_t bc[dim];

		for (size_t i = 0 ; i < dim ; i++)
			bc[i] = NON_PERIODIC;

		// Set grid and domain
		gr = grid;
		domain = dom;

		// Create a cartesian grid graph
		CartesianGraphFactory<dim, Graph_CSR<nm_v, nm_e>> g_factory_part;
		gp = g_factory_part.template construct<NO_EDGE, nm_v::id, T, dim - 1, 0, 1, 2>(gr.getSize(), domain, bc);

		// Init to 0.0 axis z (to fix in graphFactory)
		if (dim < 3)
		{
			for (size_t i = 0; i < gp.getNVertex(); i++)
			{
				gp.vertex(i).template get<nm_v::x>()[2] = 0.0;
			}
		}

		for (size_t i = 0; i < gp.getNVertex(); i++)
			gp.vertex(i).template get<nm_v::global_id>() = i;
	}

	/*! \brief Get the current graph (main)
	 *
	 * \return the current sub-sub domain Graph
	 *
	 */
	Graph_CSR<nm_v, nm_e> & getGraph()
	{
#ifdef SE_CLASS2
			check_valid(this,8);
#endif
		return gp;
	}

	/*! \brief Distribute the sub-sub-domains
	 *
	 */
	void decompose()
	{
#ifdef SE_CLASS2
			check_valid(this,8);
#endif
		Metis<Graph_CSR<nm_v, nm_e>> met(gp, v_cl.getProcessingUnits(), useWeights);
		met.onTest(testing);

		// decompose
		met.decompose<nm_v::proc_id>();
	}

	/*! \brief Refine current decomposition (NOT AVAILABLE on Metis)
	 *
	 * Disabled for MetisDistribution
	 *
	 */
	void refine()
	{
#ifdef SE_CLASS2
			check_valid(this,8);
#endif
		std::cerr << "Error: " << __FILE__ << ":" << __LINE__ << " MetisDistribution does not have refine functionality";
		ACTION_ON_ERROR(METIS_DISTRIBUTION_ERROR);
	}

	/*! \brief Function that return the position (point P1) of the sub-sub domain box in the space
	 *
	 * \param id vertex id
	 * \param pos vector that contain x, y, z
	 *
	 */
	void getSSDomainPos(size_t id, T (&pos)[dim])
	{
#ifdef SE_CLASS2
			check_valid(this,8);
#endif
		check_overflow(id);

		// Copy the geometrical informations inside the pos vector
		pos[0] = gp.vertex(id).template get<nm_v::x>()[0];
		pos[1] = gp.vertex(id).template get<nm_v::x>()[1];
		if (dim == 3)
			pos[2] = gp.vertex(id).template get<nm_v::x>()[2];
	}

	/*! \brief Checks if Computational/Migration/Communication Cost are used
	 *
	 * \return true if such weights are used
	 *
	 */
	bool weightsAreUsed()
	{
#ifdef SE_CLASS2
			check_valid(this,8);
#endif
		return useWeights;
	}

	/*! \brief function that get the computational cost of the sub-sub-domain
	 *
	 * \param id sub-sub-domain
	 *
	 * \return the comutational cost
	 *
	 */
	size_t getComputationalCost(size_t id)
	{
#ifdef SE_CLASS2
			check_valid(this,8);
#endif
		check_overflow(id);
		return gp.vertex(id).template get<nm_v::computation>();
	}

	/*! \brief Initialize all the weight
	 *
	 * Initialize Computation/Communication/Migration costs to 1
	 *
	 */
	void initWeights()
	{
#ifdef SE_CLASS2
			check_valid(this,8);
#endif
		for (size_t i = 0 ; i < getNSubSubDomains() ; i++)
		{
			setComputationCost(i,1);
			setMigrationCost(i,1);
			for (size_t j = 0 ; j < getNSubSubDomainNeighbors(i) ; j++)
				setCommunicationCost(i,j,1);
		}
	}

	/*! \brief Set computation cost on a sub-sub domain
	 *
	 * \param id sub-sub domain id
	 * \param cost
	 *
	 */
	void setComputationCost(size_t id, size_t cost)
	{
#ifdef SE_CLASS2
			check_valid(this,8);
#endif
		check_overflow(id);

		useWeights = true;

		gp.vertex(id).template get<nm_v::computation>() = cost;
	}

	/*! \brief Set migration cost on a sub-sub domain
	 *
	 * \param id of the sub-sub domain
	 * \param cost
	 */
	void setMigrationCost(size_t id, size_t cost)
	{
#ifdef SE_CLASS2
		check_valid(this,8);
#endif
		check_overflow(id);

		gp.vertex(id).template get<nm_v::migration>() = cost;
	}

	/*! \brief Set communication cost between neighborhood sub-sub-domains (weight on the edge)
	 *
	 * \param id sub-sub domain
	 * \param e id in the neighborhood list (id in the adjacency list)
	 * \param cost
	 */
	void setCommunicationCost(size_t id, size_t e, size_t cost)
	{
#ifdef SE_CLASS2
		check_valid(this,8);
#endif
		check_overflow(id);
		check_overflowe(id,e);

		gp.getChildEdge(id, e).template get<nm_e::communication>() = cost;
	}

	/*! \brief Returns total number of sub-sub-domains
	 *
	 * \return sub-sub domain numbers
	 *
	 */
	size_t getNSubSubDomains()
	{
#ifdef SE_CLASS2
			check_valid(this,8);
#endif
		return gp.getNVertex();
	}

	/*! \brief Returns total number of neighbors of one sub-sub-domain
	 *
	 * \param id of the sub-sub-domain
	 */
	size_t getNSubSubDomainNeighbors(size_t id)
	{
#ifdef SE_CLASS2
			check_valid(this,8);
#endif
		check_overflow(id);

		return gp.getNChilds(id);
	}

	/*! \brief Compute the unbalance of the processor compared to the optimal balance
	 *
	 * \return the unbalance from the optimal one 0.01 mean 1%
	 */
	float getUnbalance()
	{
#ifdef SE_CLASS2
			check_valid(this,8);
#endif
		long int min, max, sum;
		openfpm::vector<long int> loads(v_cl.getProcessingUnits());

		for (size_t i = 0; i < loads.size(); i++)
			loads.get(i) = 0;

		if (useWeights == false)
		{
			for (size_t i = 0; i < gp.getNVertex(); i++)
				loads.get(gp.vertex(i).template get<nm_v::proc_id>())++;
		}
		else
		{
			for (size_t i = 0; i < gp.getNVertex(); i++)
				loads.get(gp.vertex(i).template get<nm_v::proc_id>()) += (gp.vertex(i).template get<nm_v::computation>() == 0)?1:gp.vertex(i).template get<nm_v::computation>();
		}

		max = *std::max_element(loads.begin(), loads.end());
		min = *std::min_element(loads.begin(), loads.end());
		sum = std::accumulate(loads.begin(),loads.end(),0);

		float unbalance = ((float) (max - min)) / ((float) sum / v_cl.getProcessingUnits());

		return unbalance;
	}

	/*! \brief It set the Classs on test mode
	 *
	 * At the moment it fix the seed to have reproducible results
	 *
	 */
	void onTest()
	{
#ifdef SE_CLASS2
			check_valid(this,8);
#endif
		testing = true;
	}

	/*! \brief Write the distribution graph into file
	 *
	 * \param out output filename
	 *
	 */
	void write(std::string out)
	{
#ifdef SE_CLASS2
			check_valid(this,8);
#endif
		VTKWriter<Graph_CSR<nm_v, nm_e>, VTK_GRAPH> gv2(gp);
		gv2.write(out);

	}

	/*! \brief Compute the total computational cost of the processor
	 *
	 * \return the total computation cost
	 */
	size_t getProcessorLoad()
	{
#ifdef SE_CLASS2
			check_valid(this,8);
#endif
		size_t load = 0;

		for (size_t i = 0; i < gp.getNVertex(); i++)
		{
			if (gp.vertex(i).template get<nm_v::proc_id>() == v_cl.getProcessUnitID())
				load += gp.vertex(i).template get<nm_v::computation>();
		}

		return load;
	}

	/*! \brief operator=
	 *
	 *
	 */
	MetisDistribution & operator=(const MetisDistribution & mt)
	{
#ifdef SE_CLASS2
			check_valid(mt);
			check_valid(this,8);
#endif
		this->v_cl = mt.v_cl;
		this->gr = mt.gr;
		this->domain = mt.domain;
		this->gp = mt.gp;
		this->useWeights = mt.useWeights;
		return *this;
	}

	/*! \brief operator=
	 *
	 *
	 */
	MetisDistribution & operator=(MetisDistribution && mt)
	{
#ifdef SE_CLASS2
			check_valid(mt);
			check_valid(this,8);
#endif
		this->v_cl = mt.v_cl;
		this->gr = mt.gr;
		this->domain = mt.domain;
		this->gp.swap(mt.gp);
		this->useWeights = mt.useWeights;
		return *this;
	}

	/*! \brief operator==
	 *
	 * \return true if the distribution match
	 *
	 */
	inline bool operator==(const MetisDistribution & mt)
	{
#ifdef SE_CLASS2
			check_valid(mt);
			check_valid(this,8);
#endif
		bool ret = true;

		ret &= (this->gr == mt.gr);
		ret &= (this->domain == mt.domain);
		ret &= (this->gp == mt.gp);
		ret &= (this->useWeights == mt.useWeights);
		return ret;
	}
};

#endif /* SRC_DECOMPOSITION_METISDISTRIBUTION_HPP_ */
