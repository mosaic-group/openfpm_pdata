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
#include "VTKWriter.hpp"

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
	// Graph in metis reppresentation
	Metis_graph Mg;

	// Original graph
	Graph & g;

	/*! \brief Construct Adjacency list
	 *
	 * \param g Graph
	 *
	 */
	void constructAdjList(Graph & g)
	{
		// create xadj and adjlist

		Mg.xadj = new idx_t[g.getNVertex()+1];
		Mg.adjncy = new idx_t[g.getNEdge()];

		//! Get a vertex iterator
		auto it = g.getVertexIterator();

		//! starting point in the adjacency list
		size_t prev = 0;

		// actual position
		size_t id = 0;

		// for each vertex calculate the position of the starting point in the adjacency list
		while (it.isNext())
		{
			// Calculate the starting point in the adjacency list
			Mg.xadj[id] = prev;

			// Create the adjacency list
			for (size_t s = 0 ; s < g.getNChilds(it) ; s++)
			{Mg.adjncy[prev+s] = g.getChild(it,s);}

			// update the position for the next vertex
			prev += g.getNChilds(it);

			++it;
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
	 *
	 */
	Metis(Graph & g, size_t nc)
	:g(g)
	{
		// Get the number of vertex

		Mg.nvtxs = new idx_t[1];
		Mg.nvtxs[0] = g.getNVertex();

		// Set the number of constrains

		Mg.ncon = new idx_t[1];
		Mg.ncon[0] = 1;

		// construct the adjacency list
		constructAdjList(g);

		// Set to null the weight of the vertex

		Mg.vwgt = NULL;

		// Put the total communication size to NULL

		Mg.vsize = NULL;

		// Set to null the weight of the edge

		Mg.adjwgt = NULL;

		// Set the total number of partitions

		Mg.nparts = new idx_t[1];
		Mg.nparts[0] = nc;

		// Set to null the desired weight for each partition (one for each constrain)

		Mg.tpwgts = NULL;

		//! Set to null the partition load imbalance tollerace

		Mg.ubvec = NULL;

		//! Set tp null additional option for the graph partitioning

		Mg.options = NULL;

		//! set the objective value
		Mg.objval = new idx_t[1];

		//! Is an output vector containing the partition for each vertex
		Mg.part = new idx_t[g.getNVertex()];
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
			delete [] Mg.nvtxs;
		}

		if (Mg.ncon != NULL)
		{
			delete [] Mg.ncon;
		}

		if (Mg.xadj != NULL)
		{
			delete [] Mg.xadj;
		}

		if (Mg.adjncy != NULL)
		{
			delete [] Mg.adjncy;
		}
		if (Mg.vwgt != NULL)
		{
			delete [] Mg.vwgt;
		}
		if (Mg.adjwgt != NULL)
		{
			delete [] Mg.adjwgt;
		}
		if (Mg.nparts != NULL)
		{
			delete [] Mg.nparts;
		}
		if (Mg.tpwgts != NULL)
		{
			delete [] Mg.tpwgts;
		}
		if (Mg.ubvec != NULL)
		{
			delete [] Mg.ubvec;
		}
		if (Mg.options != NULL)
		{
			delete [] Mg.options;
		}
		if (Mg.objval != NULL)
		{
			delete [] Mg.objval;
		}
		if (Mg.part != NULL)
		{
			delete [] Mg.part;
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
			METIS_PartGraphKway(Mg.nvtxs,Mg.ncon,Mg.xadj,Mg.adjncy,Mg.vwgt,Mg.vsize,Mg.adjwgt,
				            Mg.nparts,Mg.tpwgts,Mg.ubvec,Mg.options,Mg.objval,Mg.part);

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
		METIS_PartGraphKway(Mg.nvtxs,Mg.ncon,Mg.xadj,Mg.adjncy,Mg.vwgt,Mg.vsize,Mg.adjwgt,
				            Mg.nparts,Mg.tpwgts,Mg.ubvec,Mg.options,Mg.objval,Mg.part);

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
	}
};

#endif
