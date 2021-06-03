/*
 * pdata.cpp
 *
 *  Created on: Feb 5, 2016
 *      Author: Pietro Incardona
 */
#include "pdata.hpp"
#include "SubdomainGraphNodes.hpp"
#include "memory/CudaMemory.cuh"

template<> const std::string nm_v<10>::attributes::name[] = {"x","migration","computation","global_id","id","sub_id","proc_id","id","fake_v"};
template<> const std::string nm_v<9>::attributes::name[] = {"x","migration","computation","global_id","id","sub_id","proc_id","id","fake_v"};
template<> const std::string nm_v<7>::attributes::name[] = {"x","migration","computation","global_id","id","sub_id","proc_id","id","fake_v"};
template<> const std::string nm_v<6>::attributes::name[] = {"x","migration","computation","global_id","id","sub_id","proc_id","id","fake_v"};
template<> const std::string nm_v<5>::attributes::name[] = {"x","migration","computation","global_id","id","sub_id","proc_id","id","fake_v"};
template<> const std::string nm_v<4>::attributes::name[] = {"x","migration","computation","global_id","id","sub_id","proc_id","id","fake_v"};
template<> const std::string nm_v<3>::attributes::name[] = {"x","migration","computation","global_id","id","sub_id","proc_id","id","fake_v"};
template<> const std::string nm_v<2>::attributes::name[] = {"x","migration","computation","global_id","id","sub_id","proc_id","id","fake_v"};
template<> const std::string nm_v<1>::attributes::name[] = {"x","migration","computation","global_id","id","sub_id","proc_id","id","fake_v"};


const std::string nm_e::attributes::name[] = {"communication","srcgid","dstgid"};
const std::string nm_part_v::attributes::name[] = {"id","sub_id"};
const std::string nm_part_e::attributes::name[] = {"id"};

double tot_merge = 0.0;
double tot_loc_merge = 0.0;
double tot_sendrecv = 0.0;
double tot_pack = 0.0;

/*
 * Section: Variables required for live in situ visualization
 */
openfpm::vector<handle_shmem> hgrids;
handle_shmem dtype_flag = {-1};
handle_shmem hgdb = {-1};
openfpm::vector<unsigned short*> Vis_ptrs;
long *Vis_header;
