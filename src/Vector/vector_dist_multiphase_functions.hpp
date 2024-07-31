/*
 * vector_dist_multiphase_functions.hpp
 *
 *  Created on: Oct 14, 2016
 *      Author: i-bird
 */

#ifndef SRC_VECTOR_VECTOR_DIST_MULTIPHASE_FUNCTIONS_HPP_
#define SRC_VECTOR_VECTOR_DIST_MULTIPHASE_FUNCTIONS_HPP_

#include "NN/CellList/multiphase/CellListM.hpp"
#include "NN/VerletList/VerletListM.hpp"

template<typename Vector, typename CL, typename T>
VerletList<Vector::dims,typename Vector::stype,VL_NON_SYMMETRIC,Mem_fast<>,shift<Vector::dims,typename Vector::stype>,typename Vector::internal_position_vector_type,CL>
createVerlet(Vector & v, Vector & v1, CL & cl, T r_cut)
{
	VerletList<Vector::dims,typename Vector::stype,VL_NON_SYMMETRIC,Mem_fast<>,shift<Vector::dims,typename Vector::stype>,typename Vector::internal_position_vector_type,CL> ver;

	// auto it = v.getPosVector().getIteratorTo(v.size_local());
	auto it = v.getDomainIterator();

	ver.Initialize(cl,r_cut,it,v1.getPosVector(),v.size_local());

	return ver;
}

template<unsigned int sh_byte, typename Vector , typename Vector1,typename CL, typename T> VerletListM<Vector::dims,typename Vector::stype,sh_byte,CL,shift<Vector::dims,typename Vector::stype>,typename Vector::internal_position_vector_type>
createVerletM(size_t pp, Vector & v, Vector1 & phases, CL & cl, T r_cut)
{
	VerletListM<Vector::dims,typename Vector::stype,sh_byte,CL,shift<Vector::dims,typename Vector::stype>,typename Vector::internal_position_vector_type> ver;

	openfpm::vector<pos_v<typename Vector::internal_position_vector_type>> v_phases;

	for (size_t i = 0 ; i < phases.size() ; i++)
	{v_phases.add(pos_v<typename Vector::internal_position_vector_type>(phases.get(i).getPosVector()));}

	ver.Initialize(cl,pp,r_cut,v.getPosVector(),v_phases,v.size_local());

	return ver;
}

template<unsigned int nbit, typename Vector, typename T>
CellListM<Vector::dims,typename Vector::stype,nbit,typename Vector::CellList_type>
createCellListM(openfpm::vector<Vector> & phases, T r_cut)
{
	size_t div[3];
	Box<Vector::dims,typename Vector::stype> box_cl;

	CellListM<Vector::dims,typename Vector::stype,nbit,typename Vector::CellList_type> NN;
	if (phases.size() == 0)
		return NN;

	box_cl = phases.get(0).getDecomposition().getProcessorBounds();
	phases.get(0).getCellListParams(r_cut,div,box_cl);

	NN.Initialize(box_cl,div);

	// for all the phases i
	for (size_t i = 0; i < phases.size() ; i++)
	{
		// iterate across all the particle of the phase i
		auto it = phases.get(i).getDomainAndGhostIterator();
		while (it.isNext())
		{
			auto key = it.get();

			Point<Vector::dims,T> xp = phases.get(i).getPos(key);

			// Add the particle of the phase i to the cell list
			NN.add(xp, key.getKey(), i);

			++it;
		}
	}

	return NN;
}


/////// Symmetric version

template<typename Vector,typename CL, typename T>
VerletList<Vector::dims,typename Vector::stype,VL_NON_SYMMETRIC,Mem_fast<>,shift<Vector::dims,typename Vector::stype>,typename Vector::internal_position_vector_type>
createVerletSym(Vector & v, Vector & v1, CL & cl, T r_cut)
{
	VerletList<Vector::dims,typename Vector::stype,VL_NON_SYMMETRIC,Mem_fast<>,shift<Vector::dims,typename Vector::stype>,typename Vector::internal_position_vector_type> ver;

	// auto it = v.getPosVector().getIteratorTo(v.size_local());
	auto it = v.getDomainIterator();

	ver.Initialize(cl,r_cut,it,v1.getPosVector(),v.size_local());

	return ver;
}

template<unsigned int sh_byte, typename Vector, typename Vector1 ,typename CL, typename T>
VerletListM<Vector::dims,typename Vector::stype,sh_byte,CL,shift<Vector::dims,typename Vector::stype>,typename Vector::internal_position_vector_type>
createVerletSymM(size_t pp, Vector & v, Vector1 & phases, CL & cl, T r_cut)
{
	VerletListM<Vector::dims,typename Vector::stype,sh_byte,CL,shift<Vector::dims,typename Vector::stype>,typename Vector::internal_position_vector_type> ver;

	openfpm::vector<pos_v<typename CL::internal_vector_pos_type>> v_phases;

	for (size_t i = 0 ; i < phases.size() ; i++)
	{v_phases.add(pos_v<typename CL::internal_vector_pos_type>(phases.get(i).getPosVector()));}

	ver.Initialize(cl,pp,r_cut,v.getPosVector(),v_phases,v.size_local(),VL_SYMMETRIC);

	return ver;
}

template<unsigned int nbit, typename Vector, typename T>
CellListM<Vector::dims,typename Vector::stype,nbit,typename Vector::CellList_type>
createCellListSymM(openfpm::vector<Vector> & phases, T r_cut)
{
	Box<Vector::dims,typename Vector::stype> box_cl;

	CellListM<Vector::dims,typename Vector::stype,nbit,typename Vector::CellList_type> NN;
	if (phases.size() == 0)
	{return NN;}

	// Calculate the Cell list division for this CellList
	CellDecomposer_sm<Vector::dims,typename Vector::stype,shift<Vector::dims,typename Vector::stype>> cd_sm;

	size_t pad = 0;
	cl_param_calculateSym(phases.get(0).getDecomposition().getDomain(),cd_sm,phases.get(0).getDecomposition().getGhost(),r_cut,pad);

	// Processor bounding box
	Box<Vector::dims, typename Vector::stype> pbox = phases.get(0).getDecomposition().getProcessorBounds();

	// Ghost padding extension
	Ghost<Vector::dims,size_t> g_ext(0);
	NN.Initialize(cd_sm,pbox,pad);

	// for all the phases i
	for (size_t i = 0; i < phases.size() ; i++)
	{
		// iterate across all the particle of the phase i
		auto it = phases.get(i).getDomainAndGhostIterator();
		while (it.isNext())
		{
			auto key = it.get();

			Point<Vector::dims,T> xp = phases.get(i).getPos(key);

			// Add the particle of the phase i to the cell list
			NN.add(xp, key.getKey(), i);

			++it;
		}
	}

	return NN;
}

#endif /* SRC_VECTOR_VECTOR_DIST_MULTIPHASE_FUNCTIONS_HPP_ */
