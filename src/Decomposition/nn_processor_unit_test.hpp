/*
 * nn_processor_unit_test.hpp
 *
 *  Created on: Dec 16, 2015
 *      Author: i-bird
 */

#ifndef SRC_DECOMPOSITION_NN_PROCESSOR_UNIT_TEST_HPP_
#define SRC_DECOMPOSITION_NN_PROCESSOR_UNIT_TEST_HPP_

#include "VCluster.hpp"


BOOST_AUTO_TEST_SUITE( nn_processor_test )

BOOST_AUTO_TEST_CASE( nn_processor_np_test)
{
	constexpr unsigned int dim = 2;
	typedef float T;

	// Adjacent processor for each sub-domain
	openfpm::vector<openfpm::vector<long unsigned int> > box_nn_processor;

	// Vcluster
	Vcluster & v_cl = *global_v_cluster;

	const size_t bc[dim] = {NON_PERIODIC,NON_PERIODIC};

	SpaceBox<dim,float> domain({0.0,0.0},{1.0,1.0});

	size_t sz[dim] = {8,8};
	//! Structure that store the cartesian grid information
	grid_sm<dim,void> gr(sz);

	CellDecomposer_sm<dim,T> cd;
	cd.setDimensions(domain,sz,0);

	//! Box Spacing
	T spacing[dim];

	// Calculate the total number of box and and the spacing
	// on each direction
	// Get the box containing the domain
	SpaceBox<2,T> bs = domain.getBox();

	//! the set of all local sub-domain as vector
	openfpm::vector<SpaceBox<dim,T>> sub_domains;

	/////////// From Cart decomposition ///////////

	for (unsigned int i = 0; i < dim ; i++)
	{
		// Calculate the spacing
		spacing[i] = (bs.getHigh(i) - bs.getLow(i)) / gr.size(i);
	}

	// Here we use METIS
	// Create a cartesian grid graph
	CartesianGraphFactory<dim,Graph_CSR<nm_part_v,nm_part_e>> g_factory_part;

	// the graph has only non perdiodic boundary conditions
	size_t bc_o[dim];
	for (size_t i = 0 ; i < dim ; i++)
		bc_o[i] = NON_PERIODIC;

	// sub-sub-domain graph
	Graph_CSR<nm_part_v,nm_part_e> gp = g_factory_part.template construct<NO_EDGE,T,2-1>(gr.getSize(),domain,bc_o);

	// Get the number of processing units
	size_t Np = v_cl.getProcessingUnits();

	// Get the processor id
	long int p_id = v_cl.getProcessUnitID();

	// Convert the graph to metis
	Metis<Graph_CSR<nm_part_v,nm_part_e>> met(gp,Np);

	// decompose
	met.decompose<nm_part_v::id>();

	// Optimize the decomposition creating bigger spaces
	// And reducing Ghost over-stress
	dec_optimizer<2,Graph_CSR<nm_part_v,nm_part_e>> d_o(gp,gr.getSize());

	// set of Boxes produced by the decomposition optimizer
	openfpm::vector<::Box<2,size_t>> loc_box;

	// optimize the decomposition
	d_o.template optimize<nm_part_v::sub_id,nm_part_v::id>(gp,p_id,loc_box,box_nn_processor,bc);

	// Initialize ss_box and bbox
	if (loc_box.size() >= 0)
	{
		SpaceBox<dim,size_t> sub_dc = loc_box.get(0);
		SpaceBox<dim,T> sub_d(sub_dc);
		sub_d.mul(spacing);
		sub_d.expand(spacing);

		// Fixing sub-domains to cover all the domain

		// Fixing sub_d
		// if (loc_box) is a the boundary we have to ensure that the box span the full
		// domain (avoiding rounding off error)
		for (size_t i = 0 ; i < dim ; i++)
		{
			if (sub_dc.getHigh(i) == cd.getGrid().size(i) - 1)
			{
				sub_d.setHigh(i,domain.getHigh(i));
			}
		}

		// add the sub-domain
		sub_domains.add(sub_d);
	}

	// convert into sub-domain
	for (size_t s = 1 ; s < loc_box.size() ; s++)
	{
		SpaceBox<dim,size_t> sub_dc = loc_box.get(s);
		SpaceBox<dim,T> sub_d(sub_dc);

		// re-scale and add spacing (the end is the starting point of the next domain + spacing)
		sub_d.mul(spacing);
		sub_d.expand(spacing);

		// Fixing sub-domains to cover all the domain

		// Fixing sub_d
		// if (loc_box) is a the boundary we have to ensure that the box span the full
		// domain (avoiding rounding off error)
		for (size_t i = 0 ; i < 2 ; i++)
		{
			if (sub_dc.getHigh(i) == cd.getGrid().size(i) - 1)
			{
				sub_d.setHigh(i,domain.getHigh(i));
			}
		}

		// add the sub-domain
		sub_domains.add(sub_d);

	}

	nn_prcs<dim,T> nnp(v_cl);
	nnp.create(box_nn_processor, sub_domains);

	if (v_cl.getProcessingUnits() == 1)
	{
		BOOST_REQUIRE(nnp.getNNProcessors() == 0);
	}
	else if (v_cl.getProcessingUnits() == 2)
	{
		BOOST_REQUIRE(nnp.getNNProcessors() == 1);
	}
	else
	{
		BOOST_REQUIRE(nnp.getNNProcessors() >= 1);
	}
}

BOOST_AUTO_TEST_CASE( nn_processor_box_periodic_test)
{
	constexpr unsigned int dim = 3;
	typedef float T;

	Box<dim,T> domain({0.0,0.0,0.0},{1.0,1.0,1.0});
	Point<dim,T> middle({0.5,0.5,0.5});

	const size_t bc[dim] = {PERIODIC,PERIODIC,PERIODIC};

	// Vcluster
	Vcluster & v_cl = *global_v_cluster;

	Ghost<dim,T> ghost(0.01);

	//////////////

	nn_prcs<dim,T> nnp(v_cl);

	std::unordered_map<size_t, N_box<dim,T>> & nnp_sub = nnp.get_nn_processor_subdomains();
	openfpm::vector<size_t> & nnp_np = nnp.get_nn_processors();

	// we add the boxes

	size_t tot_n = 0;
	HyperCube<dim> hyp;

	for (long int i = dim-1 ; i >= 0 ; i--)
	{
		std::vector<comb<dim>> cmbs = hyp.getCombinations_R(i);

		for (size_t j = 0 ; j < cmbs.size() ; j++)
		{
			// Create a fake processor number
			size_t prc = i;

			Point<dim,T> p1 = (middle * toPoint<dim,T>::convert(cmbs[j]) + middle)* 1.0/1.1;
			Point<dim,T> p2 = p1 + Point<dim,T>({0.1,0.1,0.1}) * 1.0/1.1;

			Box<dim,T> bx(p1,p2);
			nnp_sub[prc+1].id = prc;
			nnp_sub[prc+1].bx.add(bx);

			tot_n++;
		}
	}

	for (size_t i = 0; i < dim-1; i++)
	{
		nnp_np.add(i+1);
	}

	// check that nn_processor contain the correct boxes

//	nnp.write("nnp_output_before");

	nnp.applyBC(domain,ghost,bc);

//	nnp.write("nnp_output_after");

	BOOST_REQUIRE_EQUAL(nnp.getNearSubdomains(nnp.IDtoProc(2)).size(),6ul);
	BOOST_REQUIRE_EQUAL(nnp.getNearSubdomains(nnp.IDtoProc(0)).size(),8ul*8ul);
	BOOST_REQUIRE_EQUAL(nnp.getNearSubdomains(nnp.IDtoProc(1)).size(),12ul*4ul);
}

BOOST_AUTO_TEST_SUITE_END()

#endif /* SRC_DECOMPOSITION_NN_PROCESSOR_UNIT_TEST_HPP_ */
