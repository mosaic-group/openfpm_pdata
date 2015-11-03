/*
 * parmetis_util.hpp
 *
 *  Created on: Oct 07, 2015
 *      Author: Antonio Leo
 */

#ifndef PARMETIS_UTIL_HPP
#define PARMETIS_UTIL_HPP

#include <iostream>
#include "parmetis.h"
#include "VTKWriter.hpp"

/*! \brief Metis graph structure
 *
 * Metis graph structure
 *
 */
struct Parmetis_graph {
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

	//! Upon successful completion, the number of edges that are cut by the partitioning is written to this parameter.
	idx_t * edgecut;

	//! This parameter describes the ratio of inter-processor communication time compared to data redistri- bution time. It should be set between 0.000001 and 1000000.0. If ITR is set high, a repartitioning with a low edge-cut will be computed. If it is set low, a repartitioning that requires little data redistri- bution will be computed. Good values for this parameter can be obtained by dividing inter-processor communication time by data redistribution time. Otherwise, a value of 1000.0 is recommended.
	real_t * itr;

	//! This is used to indicate the numbering scheme that is used for the vtxdist, xadj, adjncy, and part arrays. (0 for C-style, start from 0 index)
	idx_t * numflag;

	//! This is used to indicate if the graph is weighted. wgtflag can take one of four values:
	// 0 No weights (vwgt and adjwgt are both NULL).
	// 1 Weights on the edges only (vwgt is NULL).
	// 2 Weights on the vertices only (adjwgt is NULL).
	// 3 Weights on both the vertices and edges.
	idx_t * wgtflag;
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
class Parmetis {
	// Graph in metis reppresentation
	Parmetis_graph Mg;

	// Original graph
	Graph & g;

	// Communticator for OpenMPI
	MPI_Comm comm = NULL;

	// Process rank information
	int MPI_PROC_ID = 0;

	/*! \brief Construct Adjacency list
	 *
	 * \param g Reference graph to get informations
	 *
	 */
	void constructAdjList(Graph &refGraph) {
		// create xadj and adjlist
		Mg.vwgt = new idx_t[g.getNVertex()];
		Mg.xadj = new idx_t[g.getNVertex() + 1];
		Mg.adjncy = new idx_t[g.getNEdge()];

		//! starting point in the adjacency list
		size_t prev = 0;

		// actual position
		size_t id = 0;

		// property id
		size_t real_id;

		// boolan to check if ref is the main graph
		bool main = refGraph.getNVertex() != g.getNVertex();

		// for each vertex calculate the position of the starting point in the adjacency list
		for (size_t i = 0; i < g.getNVertex(); i++) {
			// Add weight to vertex
			Mg.vwgt[i] = g.vertex(i).template get<nm_v::computation>();

			// Calculate the starting point in the adjacency list
			Mg.xadj[id] = prev;

			if (main)
				real_id = g.get_real_id(g.vertex(i).template get<nm_v::id>());
			else
				real_id = i;

			// Create the adjacency list and the weights for edges
			for (size_t s = 0; s < refGraph.getNChilds(real_id); s++) {
				size_t child = refGraph.getChild(real_id, s);

				if (main)
					Mg.adjncy[prev + s] = refGraph.vertex(child).template get<nm_v::id>();
				else
					Mg.adjncy[prev + s] = child;
			}

			// update the position for the next vertex
			prev += refGraph.getNChilds(real_id);

			id++;
		}

		// Fill the last point
		Mg.xadj[id] = prev;

		/*
		 std::cout << MPI_PROC_ID << "\n";
		 for(int i=0; i<= g.getNVertex();i++){
		 std::cout << Mg.xadj[i] << " ";
		 }
		 std::cout << "\n\n";
		 for(int i=0; i< g.getNEdge();i++){
		 std::cout << Mg.adjncy[i] << " ";
		 }
		 std::cout << "\n\n";
		 */

	}

	/*! \brief Construct Adjacency list
	 *
	 * \param g Reference graph to get informations
	 *
	 */
	/*
	void constructAdjList(Graph &refGraph, idx_t* &old_vtxdist ) {
		// create xadj and adjlist
		Mg.vwgt = new idx_t[g.getNVertex()];
		Mg.xadj = new idx_t[g.getNVertex() + 1];
		Mg.adjncy = new idx_t[g.getNEdge()];

		//! starting point in the adjacency list
		size_t prev = 0;

		// actual position
		size_t id = 0;

		// for each vertex calculate the position of the starting point in the adjacency list
		for (size_t i = 0; i < g.getNVertex(); i++) {

			// Add weight to vertex
			Mg.vwgt[i] = g.vertex(i).template get<nm_v::computation>();

			// Calculate the starting point in the adjacency list
			Mg.xadj[id] = prev;

			// Create the adjacency list and the weights for edges
			for (size_t s = 0; s < refGraph.getNChilds(i); s++) {

				size_t child = refGraph.getChild(i, s);

				// Check if child is not in this processor
				if(child > old_vtxdist[MPI_PROC_ID+1] || child < old_vtxdist[MPI_PROC_ID])

				Mg.adjncy[prev + s] = child;
			}

			// update the position for the next vertex
			prev += refGraph.getNChilds(i);

			id++;
		}

		// Fill the last point
		Mg.xadj[id] = prev;


		 std::cout << MPI_PROC_ID << "\n";
		 for(int i=0; i<= g.getNVertex();i++){
		 std::cout << Mg.xadj[i] << " ";
		 }
		 std::cout << "\n\n";
		 for(int i=0; i< g.getNEdge();i++){
		 std::cout << Mg.adjncy[i] << " ";
		 }
		 std::cout << "\n\n";


	}
	*/
public:

	/*! \brief Constructor
	 *
	 * Construct a metis graph from Graph_CSR
	 *
	 * \param g Graph we want to convert to decompose
	 * \param nc number of partitions
	 *
	 */
	Parmetis(Graph & g, size_t nc) :
			g(g) {
		// Init OpenMPI structures

		MPI_Comm_dup(MPI_COMM_WORLD, &comm);
		MPI_Comm_rank(MPI_COMM_WORLD, &MPI_PROC_ID);

		// Get the number of vertex

		Mg.nvtxs = new idx_t[1];
		Mg.nvtxs[0] = g.getNVertex();

		// Set the number of constrains

		Mg.ncon = new idx_t[1];
		Mg.ncon[0] = 1;

		// Set to null the weight of the vertex

		Mg.vwgt = NULL;

		// construct the adjacency list

		constructAdjList(g);

		// Put the total communication size to NULL

		Mg.vsize = NULL;

		// Set to null the weight of the edge

		Mg.adjwgt = NULL;

		// Set the total number of partitions

		Mg.nparts = new idx_t[1];
		Mg.nparts[0] = nc;

		//! Set option for the graph partitioning (set as default)

		Mg.options = new idx_t[4];
		Mg.options[0] = 1;
		Mg.options[1] = 3;
		Mg.options[2] = 0;
		Mg.options[3] = 0;

		//! is an output vector containing the partition for each vertex

		Mg.part = new idx_t[g.getNVertex()];
		for (int i = 0; i < g.getNVertex(); i++)
			Mg.part[i] = MPI_PROC_ID;

		//! adaptiveRepart itr value
		Mg.itr = new real_t[1];
		Mg.itr[0] = 1000.0;

		//! init tpwgts to have balanced vertices and ubvec

		Mg.tpwgts = new real_t[Mg.nparts[0]];
		Mg.ubvec = new real_t[Mg.nparts[0]];

		for (int s = 0; s < Mg.nparts[0]; s++) {
			Mg.tpwgts[s] = 1.0 / Mg.nparts[0];
			Mg.ubvec[s] = 1.05;
		}

		Mg.edgecut = new idx_t[1];
		Mg.edgecut[0] = 0;

		//! This is used to indicate the numbering scheme that is used for the vtxdist, xadj, adjncy, and part arrays. (0 for C-style, start from 0 index)
		Mg.numflag = new idx_t[1];
		Mg.numflag[0] = 0;

		//! This is used to indicate if the graph is weighted. wgtflag can take one of four values:
		Mg.wgtflag = new idx_t[1];
		Mg.wgtflag[0] = 2;
	}

	//TODO deconstruct new variables
	/*! \brief destructor
	 *
	 * Destructor, It destroy all the memory allocated
	 *
	 */
	~Parmetis() {
		// Deallocate the Mg structure
		if (Mg.nvtxs != NULL) {
			delete[] Mg.nvtxs;
		}

		if (Mg.ncon != NULL) {
			delete[] Mg.ncon;
		}

		if (Mg.xadj != NULL) {
			delete[] Mg.xadj;
		}

		if (Mg.adjncy != NULL) {
			delete[] Mg.adjncy;
		}

		if (Mg.vwgt != NULL) {
			delete[] Mg.vwgt;
		}

		if (Mg.adjwgt != NULL) {
			delete[] Mg.adjwgt;
		}

		if (Mg.nparts != NULL) {
			delete[] Mg.nparts;
		}

		if (Mg.tpwgts != NULL) {
			delete[] Mg.tpwgts;
		}

		if (Mg.ubvec != NULL) {
			delete[] Mg.ubvec;
		}

		if (Mg.options != NULL) {
			delete[] Mg.options;
		}

		if (Mg.part != NULL) {
			delete[] Mg.part;
		}

		if (Mg.edgecut != NULL) {
			delete[] Mg.edgecut;
		}

		if (Mg.numflag != NULL) {
			delete[] Mg.numflag;
		}

		if (Mg.wgtflag != NULL) {
			delete[] Mg.wgtflag;
		}
	}

	/*! \brief Decompose the graph
	 *
	 * \tparam i which property store the decomposition
	 *
	 */

	template<unsigned int i>
	void decompose(openfpm::vector<idx_t> & vtxdist) {

		// Decompose
		ParMETIS_V3_PartKway((idx_t *)vtxdist.getPointer(), Mg.xadj, Mg.adjncy, Mg.vwgt, Mg.adjwgt, Mg.wgtflag, Mg.numflag, Mg.ncon, Mg.nparts, Mg.tpwgts,
				Mg.ubvec, Mg.options, Mg.edgecut, Mg.part, &comm);
		/*
		 ParMETIS_V3_AdaptiveRepart( vtxdist, Mg.xadj,Mg.adjncy,Mg.vwgt,Mg.vsize,Mg.adjwgt, &(Mg.wgtflag), &(Mg.numflag),
		 Mg.ncon, Mg.nparts, Mg.tpwgts, Mg.ubvec, &(Mg.itr), Mg.options, &(Mg.edgecut),
		 Mg.part, &comm );
		 */

		// For each vertex store the processor that contain the data
		for (size_t j = 0, id = 0; j < g.getNVertex(); j++, id++) {

			g.vertex(j).template get<i>() = Mg.part[id];

		}

	}

	/*! \brief Refine the graph
	 *
	 * \tparam i which property store the refined decomposition
	 *
	 */

	template<unsigned int i>
	void refine(openfpm::vector<idx_t> & vtxdist) 
	{
		// Refine

		ParMETIS_V3_RefineKway((idx_t *)vtxdist.getPointer(), Mg.xadj, Mg.adjncy, Mg.vwgt, Mg.adjwgt, Mg.wgtflag, Mg.numflag, Mg.ncon, Mg.nparts, Mg.tpwgts,
				Mg.ubvec, Mg.options, Mg.edgecut, Mg.part, &comm);

		// For each vertex store the processor that contain the data

		for (size_t j = 0, id = 0; j < g.getNVertex(); j++, id++) {

			g.vertex(j).template get<i>() = Mg.part[id];

		}
	}

	/*! \brief Get graph partition vector
	 *
	 */
	idx_t* getPartition() {
		return Mg.part;
	}

	/*! \brief Reset graph and reconstruct it
	 *
	 */
	void reset(Graph & mainGraph) {
		// Deallocate the graph structures

		if (Mg.xadj != NULL) {
			delete[] Mg.xadj;
		}

		if (Mg.adjncy != NULL) {
			delete[] Mg.adjncy;
		}

		if (Mg.vwgt != NULL) {
			delete[] Mg.vwgt;
		}

		if (Mg.adjwgt != NULL) {
			delete[] Mg.adjwgt;
		}

		if (Mg.part != NULL) {
			delete[] Mg.part;
		}

		Mg.nvtxs[0] = g.getNVertex();

		Mg.part = new idx_t[g.getNVertex()];

		for (int i = 0; i < g.getNVertex(); i++)
			Mg.part[i] = MPI_PROC_ID;

		// construct the adjacency list
		constructAdjList(mainGraph);

	}

};

#endif
