/*
 * vector_dist_multiphase_functions.hpp
 *
 *  Created on: Oct 14, 2016
 *      Author: i-bird
 */

#ifndef SRC_VECTOR_VECTOR_DIST_MULTIPHASE_FUNCTIONS_HPP_
#define SRC_VECTOR_VECTOR_DIST_MULTIPHASE_FUNCTIONS_HPP_

#include "NN/CellList/CellListM.hpp"
#include "NN/VerletList/VerletListM.hpp"

template<typename Vector,typename CL, typename T> VerletList<Vector::dims,typename Vector::stype,FAST,shift<Vector::dims,typename Vector::stype>> createVerlet(Vector & v, CL & cl, T r_cut)
{
	VerletList<Vector::dims,typename Vector::stype,FAST,shift<Vector::dims,typename Vector::stype>> ver;

	ver.Initialize(cl,r_cut,v.getPosVector(),v.size_local());

	return ver;
}

template<unsigned int sh_byte, typename Vector,typename CL, typename T> VerletListM<Vector::dims,typename Vector::stype,sh_byte,shift<Vector::dims,typename Vector::stype>> createVerletM(Vector & v, CL & cl, T r_cut)
{
	VerletListM<Vector::dims,typename Vector::stype,sh_byte,shift<Vector::dims,typename Vector::stype>> ver;

	ver.Initialize(cl,r_cut,v.getPosVector(),v.size_local());

	return ver;
}

template<unsigned int nbit, typename Vector, typename T> CellListM<Vector::dims,typename Vector::stype,nbit> createCellListM(openfpm::vector<Vector> & phases, T r_cut)
{
	size_t div[3];
	Box<Vector::dims,typename Vector::stype> box_cl;

	CellListM<Vector::dims,typename Vector::stype,nbit> NN;
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

			// Add the particle of the phase i to the cell list
			NN.add(phases.get(i).getPos(key), key.getKey(), i);

			++it;
		}
	}

	return NN;
}


/////// Symmetric version

template<typename Vector,typename CL, typename T> VerletList<Vector::dims,typename Vector::stype,FAST,shift<Vector::dims,typename Vector::stype>> createVerletSym(Vector & v, CL & cl, T r_cut)
{
	VerletList<Vector::dims,typename Vector::stype,FAST,shift<Vector::dims,typename Vector::stype>> ver;

	ver.Initialize(cl,r_cut,v.getPosVector(),v.size_local());

	return ver;
}

template<unsigned int sh_byte, typename Vector,typename CL, typename T> VerletListM<Vector::dims,typename Vector::stype,sh_byte,shift<Vector::dims,typename Vector::stype>> createVerletSymM(Vector & v, CL & cl, T r_cut)
{
	VerletListM<Vector::dims,typename Vector::stype,sh_byte,shift<Vector::dims,typename Vector::stype>> ver;

	ver.Initialize(cl,r_cut,v.getPosVector(),v.size_local(),VL_SYMMETRIC);

	return ver;
}

template<unsigned int nbit, typename Vector, typename T> CellListM<Vector::dims,typename Vector::stype,nbit> createCellListSymM(openfpm::vector<Vector> & phases, T r_cut)
{
	Box<Vector::dims,typename Vector::stype> box_cl;

	CellListM<Vector::dims,typename Vector::stype,nbit> NN;
	if (phases.size() == 0)
		return NN;

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

			// Add the particle of the phase i to the cell list
			NN.add(phases.get(i).getPos(key), key.getKey(), i);

			++it;
		}
	}

	return NN;
}

#endif /* SRC_VECTOR_VECTOR_DIST_MULTIPHASE_FUNCTIONS_HPP_ */
