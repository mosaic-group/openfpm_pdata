//
// Created by aryaman on 9/2/20.
//

#ifndef SRC_GRID_DIST_ID_DEF_HPP
#define SRC_GRID_DIST_ID_DEF_HPP

#include "Decomposition/CartDecomposition.hpp"

template<unsigned int dim, typename St, typename T, typename Decomposition = CartDecomposition<dim,St>,typename Memory=HeapMemory , typename device_grid=grid_cpu<dim,T> >
class grid_dist_id;

#endif //SRC_GRID_DIST_ID_DEF_HPP
