#include <memory/ShmAllocator_manager.hpp>
#include "config.h"
#include "Vector/map_vector.hpp"

#ifndef PDATA_HPP_
#define PDATA_HPP_

constexpr int comp_host = 1;
constexpr int comp_dev = 2;

extern double tot_merge;
extern double tot_loc_merge;
extern double tot_sendrecv;
extern double tot_pack;


//! Shared memory handles for grids to be visualized
extern openfpm::vector<handle_shmem> hgrids;
//! Shared memory handle for flag indicating to visualization that grid data is to be rendered
extern handle_shmem dtype_flag;
//! Shared memory handle for storing grid extents
extern handle_shmem hgdb;
extern openfpm::vector<unsigned short *> Vis_ptrs;
extern long *Vis_header;

#endif
