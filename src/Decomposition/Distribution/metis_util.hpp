/*
 * metis_util.hpp
 *
 *  Created on: Nov 21, 2014
 *      Author: Pietro Incardona
 */

#ifndef METIS_UTIL_HPP
#define METIS_UTIL_HPP

#include <iostream>
#include "metis.h"
#include "SubdomainGraphNodes.hpp"
#include "VTKWriter/VTKWriter.hpp"

/*! \brief Metis graph structure
 *
 * Metis graph structure
 *
 */
struct Metis_graph
{
	//! The number of vertices in the graph
	idx_t * nvtxs;

	//! number of balancing constrains
	//! more practical, are the number of weights for each vertex
	//! PS even we you specify vwgt == NULL ncon must be set at leat to
	//! one
	idx_t * ncon;

	//! For each vertex it store the adjacency lost start for the vertex i
	idx_t * xadj;

	//! For each vertex it store a list of all neighborhood vertex
	idx_t * adjncy;

	//! Array that store the weight for each vertex
	idx_t * vwgt;

	//! Array of the vertex size, basically is the total communication amount
	idx_t * vsize;
	//! The weight of the edge
	idx_t * adjwgt;
	//! number of part to partition the graph
	idx_t * nparts;
	//! Desired weight for each partition (one for each constrain)
	real_t * tpwgts;
	//! For each partition load imbalance tollerated
	real_t * ubvec;
	//! Additional option for the graph partitioning
	idx_t * options;
	//! return the total comunication cost for each partition
	idx_t * objval;
	//! Is a output vector containing the partition for each vertex
	idx_t * part;
};

//! Balance communication and computation
#define BALANCE_CC 1
//! Balance communication computation and memory
#define BALANCE_CCM 2
//! Balance computation and comunication and others
#define BALANCE_CC_O(c) c+1

/*! \brief Helper class to define Metis graph
 *
 * \tparam graph structure that store the graph
 *
 */

template<typename Graph>
class Metis
{
	//! Graph in metis reppresentation
	Metis_graph Mg;

	//! Original graph
	Graph & g;

	//! indicate how many time decompose/refine/re-decompose has been called
	size_t n_dec;

	//! Distribution tolerance
	real_t dist_tol = 1.05;

	/*! \brief Construct Adjacency list
	 *
	 * \param g Graph
	 *
	 */
	void constructAdjList(Graph & g)
	{
		// create xadj and adjlist

		// create xadj, adjlist, vwgt, adjwgt and vsize
		Mg.xadj = new idx_t[g.getNVertex() + 1];
		Mg.adjncy = new idx_t[g.getNEdge()];

		//! starting point in the adjacency list
		size_t prev = 0;

		// actual position
		size_t id = 0;

		// for each vertex calculate the position of the starting point in the adjacency list
		for (size_t i = 0; i < g.getNVertex(); i++)
		{
			// Calculate the starting point in the adjacency list
			Mg.xadj[id] = prev;

			// Create the adjacency list
			for (size_t s = 0; s < g.getNChilds(i); s++)
			{
				Mg.adjncy[prev + s] = g.getChild(i, s);
			}

			// update the position for the next vertex
			prev += g.getNChilds(i);

			id++;
		}

		// Fill the last point
		Mg.xadj[id] = prev;
	}

	/*! \brief Construct Adjacency list
	 *
	 * \param g Graph
	 *
	 */
	void constructAdjListWithWeights(Graph & g)
	{
		// create xadj, adjlist, vwgt, adjwgt and vsize
		Mg.xadj = new idx_t[g.getNVertex() + 1];
		Mg.adjncy = new idx_t[g.getNEdge()];
		Mg.vwgt = new idx_t[g.getNVertex()];
		Mg.adjwgt = new idx_t[g.getNEdge()];
		Mg.vsize = new idx_t[g.getNVertex()];

		//! starting point in the adjacency list
		size_t prev = 0;

		// actual position
		size_t id = 0;

		// for each vertex calculate the position of the starting point in the adjacency list
		for (size_t i = 0; i < g.getNVertex(); i++)
		{
			// Add weight to vertex and migration cost
			Mg.vwgt[i] = g.vertex(i).template get<nm_v_computation>();
			Mg.vwgt[i] = (Mg.vwgt[i] == 0)?1:Mg.vwgt[i];
			Mg.vsize[i] = g.vertex(i).template get<nm_v_migration>();
			Mg.vsize[i] = (Mg.vsize[i] == 0)?1:Mg.vsize[i];

			// Calculate the starting point in the adjacency list
			Mg.xadj[id] = prev;

			// Create the adjacency list
			for (size_t s = 0; s < g.getNChilds(i); s++)
			{
				Mg.adjncy[prev + s] = g.getChild(i, s);

				// zero values on Metis are dangerous
				Mg.adjwgt[prev + s] = g.getChildEdge(i, s).template get<nm_e::communication>();
				Mg.adjwgt[prev + s] = (Mg.adjwgt[prev + s] == 0)?1:Mg.adjwgt[prev + s];
			}

			// update the position for the next vertex
			prev += g.getNChilds(i);

			id++;
		}

		// Fill the last point
		Mg.xadj[id] = prev;
	}

public:

	/*! \brief Constructor
	 *
	 * Construct a metis graph from Graph_CSR
	 *
	 * \param g Graph we want to convert to decompose
	 * \param nc number of partitions
	 * \param useWeights tells if weights are used or not
	 *
	 */
	Metis(Graph & g, size_t nc, bool useWeights)
	:g(g),n_dec(0)
	{
		initMetisGraph(nc,useWeights);
	}

	/*! \brief Constructor
	 *
	 * Construct a metis graph from Graph_CSR
	 *
	 * \param g Graph we want to convert to decompose
	 * \param nc number of partitions
	 *
	 */
	Metis(Graph & g, size_t nc)
	:g(g),n_dec(0)
	{
		initMetisGraph(nc,false);
	}

	/*! \brief Constructor
	 *
	 * This constructor does not initialize the internal metis graph
	 * you have to use initMetisGraph to initialize
	 *
	 * \param g Graph we want to convert to decompose
	 *
	 */
	Metis(Graph & g)
	:g(g),n_dec(0)
	{
		Mg.nvtxs = NULL;
		Mg.ncon = NULL;
		Mg.xadj = NULL;
		Mg.adjncy = NULL;
		Mg.vwgt = NULL;
		Mg.adjwgt = NULL;
		Mg.nparts = NULL;
		Mg.tpwgts = NULL;
		Mg.ubvec = NULL;
		Mg.options = NULL;
		Mg.objval = NULL;
		Mg.part = NULL;
		Mg.vsize = NULL;
	}


	/*! \brief Initialize the METIS graph
	 *
	 * \param nc number of partitions
	 * \param useWeights use the weights on the graph
	 *
	 */
	void initMetisGraph(int nc, bool useWeights)
	{

		// Get the number of vertex

		Mg.nvtxs = new idx_t[1];
		Mg.nvtxs[0] = g.getNVertex();

		// Set the number of constrains

		Mg.ncon = new idx_t[1];
		Mg.ncon[0] = 1;

		// Set to null the weight of the vertex

		Mg.vwgt = NULL;

		// Put the total communication size to NULL

		Mg.vsize = NULL;

		// Set to null the weight of the edge

		Mg.adjwgt = NULL;

		// construct the adjacency list
		if (useWeights)
			constructAdjListWithWeights(g);
		else
			constructAdjList(g);

		// Set the total number of partitions

		Mg.nparts = new idx_t[1];
		Mg.nparts[0] = nc;

		// Set to null the desired weight for each partition (one for each constrain)

		Mg.tpwgts = NULL;

		//! Set to null the partition load imbalance tolerance

		Mg.ubvec = NULL;

		//! Set tp null additional option for the graph partitioning

		Mg.options = NULL;

		//! set the objective value
		Mg.objval = new idx_t[1];

		//! Is an output vector containing the partition for each vertex
		Mg.part = new idx_t[g.getNVertex()];

		for (size_t i = 0; i < g.getNVertex(); i++)
			Mg.part[i] = 0;
	}

	/*! \brief destructor
	 *
	 * Destructor, It destroy all the memory allocated
	 *
	 */
	~Metis()
	{
		// Deallocate the Mg structure
		if (Mg.nvtxs != NULL)
		{
			delete[] Mg.nvtxs;
		}

		if (Mg.ncon != NULL)
		{
			delete[] Mg.ncon;
		}

		if (Mg.xadj != NULL)
		{
			delete[] Mg.xadj;
		}

		if (Mg.adjncy != NULL)
		{
			delete[] Mg.adjncy;
		}
		if (Mg.vwgt != NULL)
		{
			delete[] Mg.vwgt;
		}
		if (Mg.adjwgt != NULL)
		{
			delete[] Mg.adjwgt;
		}
		if (Mg.nparts != NULL)
		{
			delete[] Mg.nparts;
		}
		if (Mg.tpwgts != NULL)
		{
			delete[] Mg.tpwgts;
		}
		if (Mg.ubvec != NULL)
		{
			delete[] Mg.ubvec;
		}
		if (Mg.options != NULL)
		{
			delete[] Mg.options;
		}
		if (Mg.objval != NULL)
		{
			delete[] Mg.objval;
		}
		if (Mg.part != NULL)
		{
			delete[] Mg.part;
		}
	}

	/*! \brief Decompose the graph
	 *
	 * \tparam i which property store the decomposition
	 *
	 */

	template<unsigned int i>
	void decompose()
	{
		if (Mg.nparts[0] != 1)
		{
			// Decompose
			METIS_PartGraphRecursive(Mg.nvtxs, Mg.ncon, Mg.xadj, Mg.adjncy, Mg.vwgt, Mg.vsize, Mg.adjwgt, Mg.nparts, Mg.tpwgts, Mg.ubvec, Mg.options, Mg.objval, Mg.part);

			// vertex id

			size_t id = 0;

			// For each vertex store the processor that contain the data

			auto it = g.getVertexIterator();

			while (it.isNext())
			{
				g.vertex(it).template get<i>() = Mg.part[id];

				++id;
				++it;
			}
		}
		else
		{
			// Trivially assign all the domains to the processor 0

			auto it = g.getVertexIterator();

			while (it.isNext())
			{
				g.vertex(it).template get<i>() = 0;

				++it;
			}
		}

		n_dec++;
	}

	/*! \brief Decompose the graph
	 *
	 * \tparam i which property store the decomposition
	 *
	 */

	template<unsigned int i, typename Graph_part>
	void decompose(Graph_part & gp)
	{
		// Decompose
		METIS_PartGraphRecursive(Mg.nvtxs, Mg.ncon, Mg.xadj, Mg.adjncy, Mg.vwgt, Mg.vsize, Mg.adjwgt, Mg.nparts, Mg.tpwgts, Mg.ubvec, Mg.options, Mg.objval, Mg.part);

		// vertex id

		size_t id = 0;

		// For each vertex store the processor that contain the data

		auto it = gp.getVertexIterator();

		while (it.isNext())
		{
			gp.vertex(it).template get<i>() = Mg.part[id];

			++id;
			++it;
		}

		n_dec++;
	}

	/*! \brief It set Metis on test
	 *
	 * \param testing set to true to disable the testing
	 *
	 * At the moment disable the seed randomness to keep the result
	 * reproducible
	 *
	 */
	void onTest(bool testing)
	{
		if (testing == false)
			return;

		if (Mg.options == NULL)
		{
			// allocate
			Mg.options = new idx_t[METIS_NOPTIONS];

			// set default options
			METIS_SetDefaultOptions(Mg.options);
		}

		Mg.options[METIS_OPTION_SEED] = 0;
	}

	/*! \brief Distribution tolerance
	 *
	 * \param tol tolerance
	 *
	 */
	void setDistTol(real_t tol)
	{
		dist_tol = tol;
	}

	/*! \brief Get the decomposition counter
	 *
	 * \return the decomposition counter
	 *
	 */
	size_t get_ndec()
	{
		return n_dec;
	}

	/*! \brief Increment the decomposition counter
	 *
	 *
	 */
	void inc_dec()
	{
		n_dec++;
	}
};

#endif
