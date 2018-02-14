/*
 * SpaceDistribution.hpp
 *
 *  Created on: Dec 3, 2016
 *      Author: i-bird
 */

#ifndef SRC_DECOMPOSITION_DISTRIBUTION_SPACEDISTRIBUTION_HPP_
#define SRC_DECOMPOSITION_DISTRIBUTION_SPACEDISTRIBUTION_HPP_

#include "util/mathutil.hpp"
#include "NN/CellList/CellDecomposer.hpp"
#include "Grid/grid_key_dx_iterator_hilbert.hpp"

/*! \brief Class that distribute sub-sub-domains across processors using an hilbert curve
 *         to divide the space
 *
 * ### Initialize a Cartesian graph and decompose
 * \snippet Distribution_unit_tests.hpp Initialize a Space Cartesian graph and decompose
 *
 *
 */
template<unsigned int dim, typename T>
class SpaceDistribution
{
	//! Vcluster
	Vcluster & v_cl;

	//! Structure that store the cartesian grid information
	grid_sm<dim, void> gr;

	//! rectangular domain to decompose
	Box<dim, T> domain;

	//! Global sub-sub-domain graph
	Graph_CSR<nm_v, nm_e> gp;


public:

	/*! Constructor
	 *
	 * \param v_cl Vcluster to use as communication object in this class
	 */
	SpaceDistribution(Vcluster & v_cl)
	:v_cl(v_cl)
	{
	}

	/*! Copy constructor
	 *
	 * \param pm Distribution to copy
	 *
	 */
	SpaceDistribution(const ParMetisDistribution<dim,T> & pm)
	:v_cl(pm.v_cl)
	{
		this->operator=(pm);
	}

	/*! Copy constructor
	 *
	 * \param pm Distribution to copy
	 *
	 */
	SpaceDistribution(SpaceDistribution<dim,T> && pm)
	:v_cl(pm.v_cl)
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
		gp = g_factory_part.template construct<NO_EDGE, nm_v::id, T, dim - 1, 0>(gr.getSize(), domain, bc);

		// Init to 0.0 axis z (to fix in graphFactory)
		if (dim < 3)
		{
			for (size_t i = 0; i < gp.getNVertex(); i++)
				gp.vertex(i).template get<nm_v::x>()[2] = 0.0;
		}
		for (size_t i = 0; i < gp.getNVertex(); i++)
			gp.vertex(i).template get<nm_v::global_id>() = i;

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
		// Get the number of processing units
		size_t Np = v_cl.getProcessingUnits();

		// Calculate the best number of sub-domains for each
		// processor
		size_t N_tot = gr.size();
		size_t N_best_each = N_tot / Np;
		size_t N_rest = N_tot % Np;

		openfpm::vector<size_t> accu(Np);
		accu.get(0) = N_best_each + ((0 < N_rest)?1:0);
		for (size_t i = 1 ; i < Np ; i++)
			accu.get(i) = accu.get(i-1) + N_best_each + ((i < N_rest)?1:0);

		// Get the maximum along dimensions and take the smallest n number
		// such that 2^n < m. n it will be order of the hilbert curve

		size_t max = 0;

		for (size_t i = 0; i < dim ; i++)
		{
			if (max < gr.size(i))
				max = gr.size(i);
		}

		// Get the order of the hilbert-curve
		size_t order = openfpm::math::log2_64(max);
		if (1ul << order < max)
			order += 1;

		size_t n = 1 << order;

		// Create the CellDecomoser

		CellDecomposer_sm<dim,T> cd_sm;
		cd_sm.setDimensions(domain, gr.getSize(), 0);

		// create the hilbert curve

		//hilbert curve iterator
		grid_key_dx_iterator_hilbert<dim> h_it(order);

		T spacing[dim];

		// Calculate the hilbert curve spacing
		for (size_t i = 0 ; i < dim ; i++)
			spacing[i] = (domain.getHigh(i) - domain.getLow(i)) / n;

		// Small grid to detect already assigned sub-sub-domains
		grid_cpu<dim,aggregate<long int>> g(gr.getSize());
		g.setMemory();

		// Reset the grid to -1
		grid_key_dx_iterator<dim> it(gr);
		while (it.isNext())
		{
			auto key = it.get();

			g.template get<0>(key) = -1;

			++it;
		}

		// Go along the hilbert-curve and divide the space

		size_t proc_d = 0;
		size_t ele_d = 0;
		while (h_it.isNext())
		{
		  auto key = h_it.get();

		  // Point p
		  Point<dim,T> p;

		  for (size_t i = 0 ; i < dim ; i++)
			  p.get(i) = key.get(i) * spacing[i] + spacing[i] / 2;

		  grid_key_dx<dim> sp = cd_sm.getCellGrid(p);

		  if (g.template get<0>(sp) == -1)
		  {
			  g.template get<0>(sp) = proc_d;
			  ele_d++;

			  if (ele_d >= accu.get(proc_d))
				  proc_d++;
		  }

		  ++h_it;
		}

		// Fill from the grid to the graph

		// Reset the grid to -1
		grid_key_dx_iterator<dim> it2(gr);
		while (it2.isNext())
		{
			auto key = it2.get();

			gp.template vertex_p<nm_v::proc_id>(gr.LinId(key)) = g.template get<0>(key);

			++it2;
		}

		return;
	}

	/*! \brief Refine current decomposition
	 *
	 * Has no effect in this case
	 *
	 */
	void refine()
	{
		std::cout << __FILE__ << ":" << __LINE__ << " You are trying to dynamicaly balance a fixed decomposition, this operation has no effect" << std::endl;
	}

	/*! \brief Compute the unbalance of the processor compared to the optimal balance
	 *
	 * \return the unbalance from the optimal one 0.01 mean 1%
	 */
	float getUnbalance()
	{
		return gr.size() % v_cl.getProcessingUnits();
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
		std::cout << __FILE__ << ":" << __LINE__ << " You are trying to set the computation cost on a fixed decomposition, this operation has no effect" << std::endl;
	}

	/*! \brief Checks if weights are used on the vertices
	 *
	 * \return true if weights are used in the decomposition
	 */
	bool weightsAreUsed()
	{
		return false;
	}

	/*! \brief function that get the weight of the vertex
	 *
	 * \param id vertex id
	 *
	 * \return the weight of the vertex
	 *
	 */
	size_t getSubSubDomainComputationCost(size_t id)
	{
		return 1.0;
	}

	/*! \brief Compute the processor load counting the total weights of its vertices
	 *
	 * \return the computational load of the processor graph
	 */
	size_t getProcessorLoad()
	{
		// Get the number of processing units
		size_t Np = v_cl.getProcessingUnits();

		// Calculate the best number of sub-domains for each
		// processor
		size_t N_tot = gr.size();
		size_t N_best_each = N_tot / Np;
		size_t N_rest = N_tot % Np;

		if (v_cl.getProcessUnitID() < N_rest)
			N_best_each += 1;

		return N_best_each;
	}

	/*! \brief Set migration cost of the vertex id
	 *
	 * \param id of the vertex to update
	 * \param migration cost of the migration
	 */
	void setMigrationCost(size_t id, size_t migration)
	{
	}

	/*! \brief Set communication cost of the edge id
	 *
	 * \param v_id Id of the source vertex of the edge
	 * \param e i child of the vertex
	 * \param communication Communication value
	 */
	void setCommunicationCost(size_t v_id, size_t e, size_t communication)
	{
	}

	/*! \brief Returns total number of sub-sub-domains in the distribution graph
	 *
	 * \return number of sub-sub-domain
	 *
	 */
	size_t getNSubSubDomains()
	{
		return gp.getNVertex();
	}

	/*! \brief Returns total number of neighbors of the sub-sub-domain id
	 *
	 * \param id id of the sub-sub-domain
	 */
	size_t getNSubSubDomainNeighbors(size_t id)
	{
		return gp.getNChilds(id);
	}

	/*! \brief Print the current distribution and save it to VTK file
	 *
	 * \param file filename
	 *
	 */
	void write(const std::string & file)
	{
		VTKWriter<Graph_CSR<nm_v, nm_e>, VTK_GRAPH> gv2(gp);
		gv2.write(std::to_string(v_cl.getProcessUnitID()) + "_" + file + ".vtk");
	}

	const SpaceDistribution<dim,T> & operator=(const SpaceDistribution<dim,T> & dist)
	{
		gr = dist.gr;
		domain = dist.domain;
		gp = dist.gp;

		return *this;
	}

	const SpaceDistribution<dim,T> & operator=(SpaceDistribution<dim,T> && dist)
	{
		v_cl = dist.v_cl;
		gr = dist.gr;
		domain = dist.domain;
		gp.swap(dist.gp);

		return *this;
	}

	/*! \brief It return the decomposition id
	 *
	 * It just return 0
	 *
	 * \return 0
	 *
	 */
	size_t get_ndec()
	{
		return 0;
	}
};


#endif /* SRC_DECOMPOSITION_DISTRIBUTION_SPACEDISTRIBUTION_HPP_ */
