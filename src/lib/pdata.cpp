/*
 * pdata.cpp
 *
 *  Created on: Feb 5, 2016
 *      Author: Pietro Incardona
 */

#include "SubdomainGraphNodes.hpp"

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



