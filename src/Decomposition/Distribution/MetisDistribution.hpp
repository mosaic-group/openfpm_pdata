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

#define METIS_DISTRIBUTION_ERROR_OBJECT std::runtime_error("Metis runtime error");

/*! \brief sub-domain list and weight
 *
 */
struct met_sub_w
{
	//! sub-domain id
	size_t id;

	//! sub-domain weight / assignment (it depend in which context is used)
	size_t w;

	static bool noPointers() {return true;}
};

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
	Vcluster<> & v_cl;

	//! Structure that store the cartesian grid information
	grid_sm<dim, void> gr;

	//! rectangular domain to decompose
	Box<dim, T> domain;

	//! Global sub-sub-domain graph
	Graph_CSR<nm_v<dim>, nm_e> & gp;

	//! Flag that indicate if we are doing a test (In general it fix the seed)
	bool testing = false;

	//! Metis decomposer utility
	Metis<Graph_CSR<nm_v<dim>, nm_e>> metis_graph;

	//! unordered map that map global sub-sub-domain to owned_cost_sub id
	std::unordered_map<size_t,size_t> owner_scs;

	//! list owned sub-sub-domains set for computation cost
	openfpm::vector<met_sub_w> owner_cost_sub;

	//! received assignment
	openfpm::vector<met_sub_w> recv_ass;

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
			ACTION_ON_ERROR(METIS_DISTRIBUTION_ERROR_OBJECT)
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
			ACTION_ON_ERROR(METIS_DISTRIBUTION_ERROR_OBJECT)
		}
#endif
	}

public:

	typedef Graph_CSR<nm_v<dim>, nm_e> graph_type; 

	static constexpr unsigned int computation = nm_v_computation;

	/*! \brief constructor
	 *
	 * \param v_cl vcluster
	 *
	 */
	MetisDistribution(Vcluster<> & v_cl, Graph_CSR<nm_v<dim>, nm_e> & gp)
	:v_cl(v_cl),gp(gp),metis_graph(gp)
	{
#ifdef SE_CLASS2
			check_new(this,8,VECTOR_EVENT,1);
#endif
	}

	/*! \brief Copy constructor
	 *
	 * \param mt distribution to copy
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
	 * \param mt distribution to copy
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
		CartesianGraphFactory<dim, Graph_CSR<nm_v<dim>, nm_e>> g_factory_part;
		gp = g_factory_part.template construct<NO_EDGE, nm_v_id, T, dim - 1, 0>(gr.getSize(), domain, bc);

		// Init to 0.0 axis z (to fix in graphFactory)
		if (dim < 3)
		{
			for (size_t i = 0; i < gp.getNVertex(); i++)
			{
				gp.vertex(i).template get<nm_v_x>()[2] = 0.0;
			}
		}

		for (size_t i = 0; i < gp.getNVertex(); i++)
			gp.vertex(i).template get<nm_v_global_id>() = i;
	}

	/*! \brief Get the current graph (main)
	 *
	 * \return the current sub-sub domain Graph
	 *
	 */
	Graph_CSR<nm_v<dim>, nm_e> & getGraph()
	{
#ifdef SE_CLASS2
			check_valid(this,8);
#endif
		return gp;
	}

	/*! \brief Distribute the sub-sub-domains
	 *
	 */
	void distribute()
	{
#ifdef SE_CLASS2
			check_valid(this,8);
#endif

		// Gather the sub-domain weight in one processor
		recv_ass.clear();
		v_cl.SGather(owner_cost_sub,recv_ass,0);

		if (v_cl.getProcessUnitID() == 0)
		{
			if (recv_ass.size() != 0)
			{
				// we fill the assignment
				for (size_t i = 0 ; i < recv_ass.size() ; i++)
					gp.template vertex_p<nm_v_computation>(recv_ass.get(i).id) = recv_ass.get(i).w;

				metis_graph.initMetisGraph(v_cl.getProcessingUnits(),true);
			}
			else
				metis_graph.initMetisGraph(v_cl.getProcessingUnits(),false);
			metis_graph.onTest(testing);

			// decompose
			metis_graph.template decompose<nm_v_proc_id>();

			if (recv_ass.size() != 0)
			{
				// we fill the assignment
				for (size_t i = 0 ; i < recv_ass.size() ; i++)
					recv_ass.get(i).w = gp.template vertex_p<nm_v_proc_id>(recv_ass.get(i).id);
			}
			else
			{
				recv_ass.resize(gp.getNVertex());

				// we fill the assignment
				for (size_t i = 0 ; i < gp.getNVertex() ; i++)
				{
					recv_ass.get(i).id = i;
					recv_ass.get(i).w = gp.template vertex_p<nm_v_proc_id>(i);
				}
			}
		}
		else
		{
			metis_graph.inc_dec();
		}

		recv_ass.resize(gp.getNVertex());

		// broad cast the result
		v_cl.Bcast(recv_ass,0);
		v_cl.execute();
		owner_scs.clear();
		owner_cost_sub.clear();

		size_t j = 0;

		// Fill the metis graph
		for (size_t i = 0 ; i < recv_ass.size() ; i++)
		{
			gp.template vertex_p<nm_v_proc_id>(recv_ass.get(i).id) = recv_ass.get(i).w;

			if (recv_ass.get(i).w == v_cl.getProcessUnitID())
			{
				owner_scs[recv_ass.get(i).id] = j;
				j++;
				owner_cost_sub.add();
				owner_cost_sub.last().id = recv_ass.get(i).id;
				owner_cost_sub.last().w = 1;
			}
		}
	}

	/*! \brief Refine current decomposition
	 *
	 * In metis case it just re-decompose
	 *
	 */
	void refine()
	{
#ifdef SE_CLASS2
			check_valid(this,8);
#endif

		distribute();
	}

	/*! \brief Redecompose current decomposition
	 *
	 */
	void redecompose()
	{
#ifdef SE_CLASS2
			check_valid(this,8);
#endif
		distribute();
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
		pos[0] = gp.vertex(id).template get<nm_v_x>()[0];
		pos[1] = gp.vertex(id).template get<nm_v_x>()[1];
		if (dim == 3)
		{pos[2] = gp.vertex(id).template get<nm_v_x>()[2];}
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
		return gp.vertex(id).template get<nm_v_computation>();
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
#ifdef SE_CLASS1
		check_overflow(id);
#endif

		auto fnd = owner_scs.find(id);
		if (fnd == owner_scs.end())
		{
			std::cerr << __FILE__ << ":" << __LINE__ << " Error you are setting a sub-sub-domain the processor does not own" << std::endl;
		}
		else
		{
			size_t id = fnd->second;
			owner_cost_sub.get(id).w = cost;
		}
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
#ifdef SE_CLASS1
		check_overflow(id);
#endif

		gp.vertex(id).template get<nm_v_migration>() = cost;
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
#ifdef SE_CLASS1
		check_overflow(id);
		check_overflowe(id,e);
#endif

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
#ifdef SE_CLASS1
		check_overflow(id);
#endif

		return gp.getNChilds(id);
	}

	/*! \brief Compute the unbalance of the processor compared to the optimal balance
	 *
	 * \warning all processor must call this function
	 *
	 * \return the unbalance from the optimal one 0.01 mean 1%
	 */
	float getUnbalance()
	{
#ifdef SE_CLASS2
			check_valid(this,8);
#endif
		size_t load_p = getProcessorLoad();

		float load_avg = load_p;
		v_cl.sum(load_avg);
		v_cl.execute();

		if (load_avg == 0)
		{
			// count the number if sub-sub-domain assigned
			load_avg = owner_cost_sub.size();

			v_cl.sum(load_avg);
			v_cl.execute();
		}

		load_avg /= v_cl.getProcessingUnits();

		return ((float)load_p - load_avg) / load_avg;
	}

	/*! \brief Return the total number of sub-sub-domains in the distribution graph
	 *
	 * \return the total number of sub-sub-domains set
	 *
	 */
	size_t getNOwnerSubSubDomains() const
	{
		return owner_cost_sub.size();
	}

	/*! \brief Return the id of the set sub-sub-domain
	 *
	 * \param id id in the list of the set sub-sub-domains
	 *
	 * \return the id
	 *
	 */
	size_t getOwnerSubSubDomain(size_t id) const
	{
		return owner_cost_sub.get(id).id;
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

		VTKWriter<Graph_CSR<nm_v<dim>, nm_e>, VTK_GRAPH> gv2(gp);
		gv2.write(std::to_string(v_cl.getProcessUnitID()) + "_" + out + ".vtk");

	}

	/*! \brief Compute the processor load
	 *
	 * \warning all processors must call this function
	 *
	 * \return the total computation cost
	 */
	size_t getProcessorLoad()
	{
#ifdef SE_CLASS2
			check_valid(this,8);
#endif
		openfpm::vector<size_t> loads(v_cl.getProcessingUnits());

		size_t load = 0;

		if (v_cl.getProcessUnitID() == 0)
		{
			for (size_t i = 0; i < gp.getNVertex(); i++)
			{loads.get(gp.template vertex_p<nm_v_proc_id>(i)) += gp.template vertex_p<nm_v_computation>(i);}

			for (size_t i = 0 ; i < v_cl.getProcessingUnits() ; i++)
			{
				v_cl.send(i,1234,&loads.get(i),sizeof(size_t));
			}
		}
		v_cl.recv(0,1234,&load,sizeof(size_t));
		v_cl.execute();

		return load;
	}

	/*! \brief operator=
	 *
	 * \param mt object to copy
	 *
	 * \return itself
	 *
	 */
	MetisDistribution & operator=(const MetisDistribution & mt)
	{
#ifdef SE_CLASS2
			check_valid(&mt,8);
			check_valid(this,8);
#endif
		this->gr = mt.gr;
		this->domain = mt.domain;
		this->gp = mt.gp;
		this->owner_cost_sub = mt.owner_cost_sub;
		this->owner_scs = mt.owner_scs;
		return *this;
	}

	/*! \brief operator=
	 *
	 * \param mt object to copy
	 *
	 * \return itself
	 *
	 */
	MetisDistribution & operator=(MetisDistribution && mt)
	{
#ifdef SE_CLASS2
			check_valid(mt);
			check_valid(this,8);
#endif
		this->gr = mt.gr;
		this->domain = mt.domain;
		this->gp.swap(mt.gp);
		this->owner_cost_sub.swap(mt.owner_cost_sub);
		this->owner_scs.swap(mt.owner_scs);
		return *this;
	}

	/*! \brief operator==
	 *
	 * \param mt Metis distribution to compare with
	 *
	 * \return true if the distribution match
	 *
	 */
	inline bool operator==(const MetisDistribution & mt)
	{
#ifdef SE_CLASS2
			check_valid(&mt,8);
			check_valid(this,8);
#endif
		bool ret = true;

		ret &= (this->gr == mt.gr);
		ret &= (this->domain == mt.domain);
		ret &= (this->gp == mt.gp);

		return ret;
	}

	/*! \brief Set the tolerance for each partition
	 *
	 * \param tol tolerance
	 *
	 */
	void setDistTol(double tol)
	{
		metis_graph.setDistTol(tol);
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

		auto fnd = owner_scs.find(id);
		if (fnd == owner_scs.end())
		{
			std::cerr << __FILE__ << ":" << __LINE__ << " Error you are setting a sub-sub-domain that the processor does not own" << std::endl;
			return 0;
		}

		size_t ids = fnd->second;
		return owner_cost_sub.get(ids).w;
	}

	/*! \brief Get the decomposition counter
	 *
	 * \return the decomposition counter
	 *
	 */
	size_t get_ndec()
	{
		return metis_graph.get_ndec();
	}
};

#endif /* SRC_DECOMPOSITION_METISDISTRIBUTION_HPP_ */
