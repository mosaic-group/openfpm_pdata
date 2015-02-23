#ifndef DIST_MODEL_HPP
#define DIST_MODEL_HPP

#include "metis.h"

/*! \brief This class do graph partitioning
 *
 * This class do graph partitioning, it use METIS internaly
 *
 */

template<typename Graph, typename algorithm>
class GraphPartitioning
{
	//! Structure that store the graph
	Graph & grp;

	/*! Constructor
	 *
	 * It load the graph to partition
	 *
	 * \param g Graph to store
	 *
	 */
	GraphPartitioning(Graph & g)
	:grp(g)
	{}
};

#endif
