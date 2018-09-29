#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "Decomposition/CartDecomposition.hpp"
#include "util/mathutil.hpp"

BOOST_AUTO_TEST_SUITE (CartDecomposition_test)

#define SUB_UNIT_FACTOR 1024

void setComputationCosts(CartDecomposition<2, float> &dec, size_t n_v, Point<2, float> center, float radius, size_t weight_h, size_t weight_l)
{
	float radius2 = pow(radius, 2);
	float eq;

	// Position structure for the single vertex
	float pos[2];

	for (size_t i = 0; i < n_v; i++)
	{
		dec.getSubSubDomainPosition(i, pos);

		eq = pow((pos[0] - center.get(0)), 2) + pow((pos[1] - center.get(1)), 2);

		if (eq <= radius2)
			dec.setSubSubDomainComputationCost(i, weight_h);
		else
			dec.setSubSubDomainComputationCost(i, weight_l);
	}
}

void setComputationCosts3D(CartDecomposition<3, float> &dec, size_t n_v, Point<3, float> center, float radius, size_t weight_h, size_t weight_l)
{
	float radius2 = radius * radius;
	float eq;

	// Position structure for the single vertex
	float pos[3];

	for (size_t i = 0; i < n_v; i++)
	{
		dec.getSubSubDomainPosition(i, pos);

		eq = pow((pos[0] - center.get(0)), 2) + pow((pos[1] - center.get(1)), 2) + pow((pos[2] - center.get(2)), 2);

		if (eq <= radius2)
			dec.setSubSubDomainComputationCost(i, weight_h);
		else
			dec.setSubSubDomainComputationCost(i, weight_l);
	}
}



BOOST_AUTO_TEST_CASE( CartDecomposition_non_periodic_test)
{
	// Vcluster
	Vcluster<> & vcl = create_vcluster();

	CartDecomposition<3, float> dec(vcl);

	// Physical domain
	Box<3, float> box( { 0.0, 0.0, 0.0 }, { 1.0, 1.0, 1.0 });
	size_t div[3];

	// Get the number of processor and calculate the number of sub-domain
	// for each processor (SUB_UNIT_FACTOR=64)
	size_t n_proc = vcl.getProcessingUnits();
	size_t n_sub = n_proc * SUB_UNIT_FACTOR;

	// Set the number of sub-domains on each dimension (in a scalable way)
	for (int i = 0; i < 3; i++)
	{	div[i] = openfpm::math::round_big_2(pow(n_sub,1.0/3));}

	// Define ghost
	Ghost<3, float> g(0.01);

	// Boundary conditions
	size_t bc[] = { NON_PERIODIC, NON_PERIODIC, NON_PERIODIC };

	// Decompose
	dec.setParameters(div,box,bc,g);
	dec.decompose();

	// For each calculated ghost box
	for (size_t i = 0; i < dec.getNIGhostBox(); i++)
	{
		SpaceBox<3,float> b = dec.getIGhostBox(i);
		size_t proc = dec.getIGhostBoxProcessor(i);

		// sample one point inside the box
		Point<3,float> p = b.rnd();

		// Check that ghost_processorsID return that processor number
		const openfpm::vector<size_t> & pr = dec.ghost_processorID<CartDecomposition<3,float>::processor_id>(p);

		bool found = false;

		for (size_t j = 0; j < pr.size(); j++)
		{
			if (pr.get(j) == proc)
			{	found = true; break;}
		}

		if (found == false)
		{
			const openfpm::vector<size_t> pr2 = dec.ghost_processorID<CartDecomposition<3,float>::processor_id>(p);
		}

		BOOST_REQUIRE_EQUAL(found,true);
	}

	// Check the consistency

	bool val = dec.check_consistency();
	BOOST_REQUIRE_EQUAL(val,true);

	// We duplicate the decomposition
	CartDecomposition<3, float> dec2 = dec.duplicate();
	dec2.check_consistency();

	// check that dec and dec2 contain the same information
	bool ret = dec.is_equal(dec2);

	// We check if the two decomposition are equal
	BOOST_REQUIRE_EQUAL(ret,true);

	// We duplicate the decomposition redefining the ghost

	// Define ghost
	Ghost<3, float> g3(0.005);

	// We duplicate the decomposition redefining the ghost
	CartDecomposition<3, float> dec3 = dec.duplicate(g3);

	ret = dec3.check_consistency();
	BOOST_REQUIRE_EQUAL(ret,true);

	// Check that dec3 is equal to dec2 with the exception of the ghost part
	ret = dec3.is_equal_ng(dec2);
	BOOST_REQUIRE_EQUAL(ret,true);
}

BOOST_AUTO_TEST_CASE( CartDecomposition_periodic_test)
{
	// Vcluster
	Vcluster<> & vcl = create_vcluster();

	//! [Create CartDecomposition]
	CartDecomposition<3, float> dec(vcl);

	// Physical domain
	Box<3, float> box( { 0.0, 0.0, 0.0 }, { 1.0, 1.0, 1.0 });
	size_t div[3];

	// Get the number of processor and calculate the number of sub-domain
	// for each processor (SUB_UNIT_FACTOR=64)
	size_t n_proc = vcl.getProcessingUnits();
	size_t n_sub = n_proc * SUB_UNIT_FACTOR;

	// Set the number of sub-domains on each dimension (in a scalable way)
	for (int i = 0; i < 3; i++)
	{	div[i] = openfpm::math::round_big_2(pow(n_sub,1.0/3));}

	// Define ghost
	Ghost<3, float> g(0.01);

	// Boundary conditions
	size_t bc[] = { PERIODIC, PERIODIC, PERIODIC };

	// Decompose
	dec.setParameters(div,box,bc,g);
	dec.decompose();

	//! [Create CartDecomposition]

	// For each calculated ghost box
	for (size_t i = 0; i < dec.getNIGhostBox(); i++)
	{
		SpaceBox<3,float> b = dec.getIGhostBox(i);
		size_t proc = dec.getIGhostBoxProcessor(i);

		// sample one point inside the box
		Point<3,float> p = b.rnd();

		// Check that ghost_processorsID return that processor number
		const openfpm::vector<size_t> & pr = dec.ghost_processorID<CartDecomposition<3,float>::processor_id>(p);

		bool found = false;

		for (size_t j = 0; j < pr.size(); j++)
		{
			if (pr.get(j) == proc)
			{	found = true; break;}
		}

		if (found == false)
		{
			const openfpm::vector<size_t> pr2 = dec.ghost_processorID<CartDecomposition<3,float>::processor_id>(p);
		}

		BOOST_REQUIRE_EQUAL(found,true);
	}

	// Check the consistency
	bool val = dec.check_consistency();
	BOOST_REQUIRE_EQUAL(val,true);

	// We duplicate the decomposition
	CartDecomposition<3, float> dec2 = dec.duplicate();
	dec2.check_consistency();

	bool ret = dec.is_equal(dec2);

	// We check if the two decomposition are equal
	BOOST_REQUIRE_EQUAL(ret,true);

	// check that dec and dec2 contain the same information

	// We duplicate the decomposition redefining the ghost

	// Define ghost
	Ghost<3, float> g3(0.005);

	// We duplicate the decomposition refefining the ghost
	CartDecomposition<3, float> dec3 = dec.duplicate(g3);

	ret = dec3.check_consistency();
	BOOST_REQUIRE_EQUAL(ret,true);

	// Check that g3 is equal to dec2 with the exception of the ghost part
	ret = dec3.is_equal_ng(dec2);
	BOOST_REQUIRE_EQUAL(ret,true);
}


////////////////////////// CartDecomposition extended

BOOST_AUTO_TEST_CASE( CartDecomposition_ext_non_periodic_test)
{
	// Vcluster
	Vcluster<> & vcl = create_vcluster();

	CartDecomposition<3,float> dec(vcl);

	// Physical domain
	Box<3,float> box({0.0,0.0,0.0},{1.0,1.0,1.0});
	size_t div[3];

	// Get the number of processor and calculate the number of sub-domain
	// for each processor (SUB_UNIT_FACTOR=64)
	size_t n_proc = vcl.getProcessingUnits();
	size_t n_sub = n_proc * SUB_UNIT_FACTOR;

	// Set the number of sub-domains on each dimension (in a scalable way)
	for (int i = 0 ; i < 3 ; i++)
	{div[i] = openfpm::math::round_big_2(pow(n_sub,1.0/3));}

	// Define ghost
	Ghost<3,float> g(0.01);

	// Boundary conditions
	size_t bc[] = {NON_PERIODIC,NON_PERIODIC,NON_PERIODIC};

	// Decompose
	dec.setParameters(div,box,bc,g);
	dec.decompose();

	//! [Extend CartDecomposition]

	Box<3,float> box_ext({-0.1,-0.1,-0.1},{1.1,1.1,1.1});

	// Use the old decomposition to extend on a bigger domain
	CartDecomposition_ext<3,float> dec_ext(vcl);

	dec_ext.setParameters(dec,g,box_ext);

	//! [Extend CartDecomposition]

	// Check the new decomposition is fully contained in the new one, and there are
	// box not fully contained i the old box

	BOOST_REQUIRE_EQUAL(dec_ext.getNSubDomain(),dec.getNSubDomain());

	double volume = 0.0;
	for (size_t i = 0; i < dec_ext.getNSubDomain() ; i++)
	{
		volume += dec_ext.getSubDomain(i).getVolume();
		BOOST_REQUIRE_EQUAL(dec_ext.getSubDomain(i).isContained(dec.getSubDomain(i)),true);
	}

	vcl.sum(volume);
	vcl.execute();

	BOOST_REQUIRE_CLOSE(volume,1.728,0.0001);

	BOOST_REQUIRE_EQUAL(dec.getNNProcessors(),dec_ext.getNNProcessors());

	double volume_g = 0.0;
	double volume_ge = 0.0;
	for (size_t p = 0; p < dec.getNNProcessors(); p++)
	{
		BOOST_REQUIRE_EQUAL(dec.getProcessorNEGhost(p),dec_ext.getProcessorNEGhost(p));
		for (size_t i = 0; i < dec.getProcessorNEGhost(p); i++)
		{
			volume_g += dec.getProcessorEGhostBox(p,i).getVolume();
			volume_ge += dec_ext.getProcessorEGhostBox(p,i).getVolume();

			BOOST_REQUIRE_EQUAL(dec_ext.getProcessorEGhostBox(p,i).isContained(dec_ext.getProcessorEGhostBox(p,i)),true);
		}
	}

	vcl.sum(volume_g);
	vcl.sum(volume_ge);
	vcl.execute();

	if (vcl.getProcessingUnits() > 1)
	{
		BOOST_REQUIRE(volume_ge > volume_g*1.05);
	}

	volume_g = 0.0;
	volume_ge = 0.0;
	for (size_t p = 0; p < dec.getNNProcessors(); p++)
	{
		for (size_t i = 0; i< dec.getProcessorNIGhost(p); i++)
		{
			volume_g += dec.getProcessorIGhostBox(p,i).getVolume();
			volume_ge += dec_ext.getProcessorIGhostBox(p,i).getVolume();
			BOOST_REQUIRE_EQUAL(dec_ext.getProcessorIGhostBox(p,i).isContained(dec.getProcessorIGhostBox(p,i)),true);
		}
	}

	vcl.sum(volume_g);
	vcl.sum(volume_ge);
	vcl.execute();

	if (vcl.getProcessingUnits() > 1)
	{
		BOOST_REQUIRE(volume_ge > volume_g*1.05);
	}
}

BOOST_AUTO_TEST_CASE( CartDecomposition_check_cross_consistency_between_proc_idbc_and_ghost )
{
	// Vcluster
	Vcluster<> & vcl = create_vcluster();

	if (vcl.size() != 3)
	{return;}

	CartDecomposition<3, double> dec(vcl);

	size_t bc[3] = {PERIODIC,PERIODIC,PERIODIC};

	// Physical domain
	Box<3, double> box( { -0.01, -0.01, 0.0 }, { 0.01, 0.01, 0.003 });

	Ghost<3,double> g(0.0015);

	dec.setGoodParameters(box, bc, g, 512);

	dec.decompose();

	// Now we check the point

	Point<3,double> p1({-0.0067499999999999999237,-0.0012499999999999995923,0.001250000000000000026});
	Point<3,double> p2({-0.0067499999999999999237,-0.0012499999999999993755,0.001250000000000000026});

	size_t proc1 = dec.processorIDBC(p1);
	size_t proc2 = dec.processorIDBC(p2);

	const openfpm::vector<std::pair<size_t, size_t>> & vp_id1 = dec.template ghost_processorID_pair<typename CartDecomposition<3, double>::lc_processor_id, typename CartDecomposition<3, double>::shift_id>(p1, UNIQUE);
	const openfpm::vector<std::pair<size_t, size_t>> & vp_id2 = dec.template ghost_processorID_pair<typename CartDecomposition<3, double>::lc_processor_id, typename CartDecomposition<3, double>::shift_id>(p2, UNIQUE);

	if (proc1 != proc2)
	{
		if (vcl.rank() == proc2)
		{
			BOOST_REQUIRE(vp_id2.size() != 0);
			BOOST_REQUIRE(vp_id1.size() == 0);
		}

		if (vcl.rank() == proc1)
		{
			BOOST_REQUIRE(vp_id2.size() == 0 );
			BOOST_REQUIRE(vp_id1.size() != 0 );
		}
	}
}

BOOST_AUTO_TEST_CASE( CartDecomposition_check_cross_consistency_between_proc_idbc_and_ghost2 )
{
	// Vcluster
	Vcluster<> & vcl = create_vcluster();

	CartDecomposition<3, double> dec(vcl);

	size_t bc[3] = {PERIODIC,PERIODIC,PERIODIC};

	// Physical domain
	Box<3, double> box( { -0.01, -0.01, 0.0 }, { 0.01, 0.01, 0.003 });

	Ghost<3,double> g(0.0015);

	dec.setGoodParameters(box, bc, g, 512);

	dec.decompose();

	// Now we check the point

	for (size_t j = 0 ; j < 3 ; j++ )
	{
		for (size_t i = 0 ; i < dec.getNSubDomain() ; i++)
		{
			Point<3,double> p1;
			Point<3,double> p2;

			p1.get(0) = SpaceBox<3,double>(dec.getSubDomains().get(i)).getLow(0);
			p1.get(1) = SpaceBox<3,double>(dec.getSubDomains().get(i)).getLow(1);
			p1.get(2) = SpaceBox<3,double>(dec.getSubDomains().get(i)).getLow(2);

			p2 = p1;

			p2.get(j) = std::nextafter(SpaceBox<3,double>(dec.getSubDomains().get(i)).getLow(j),-1.0);

			size_t proc1 = dec.processorIDBC(p1);
			size_t proc2 = dec.processorIDBC(p2);

			BOOST_REQUIRE(proc1 < vcl.size());
			BOOST_REQUIRE(proc2 < vcl.size());

			const openfpm::vector<std::pair<size_t, size_t>> & vp_id1 = dec.template ghost_processorID_pair<typename CartDecomposition<3, double>::lc_processor_id, typename CartDecomposition<3, double>::shift_id>(p1, UNIQUE);
			const openfpm::vector<std::pair<size_t, size_t>> & vp_id2 = dec.template ghost_processorID_pair<typename CartDecomposition<3, double>::lc_processor_id, typename CartDecomposition<3, double>::shift_id>(p2, UNIQUE);

			if (proc1 != proc2)
			{
				if (vcl.rank() == proc2)
				{
					BOOST_REQUIRE(vp_id2.size() != 0);
					BOOST_REQUIRE(vp_id1.size() == 0);
				}

				if (vcl.rank() == proc1)
				{
					BOOST_REQUIRE(vp_id2.size() == 0 );
					BOOST_REQUIRE(vp_id1.size() != 0 );
				}
			}


			p1.get(0) = std::nextafter(SpaceBox<3,double>(dec.getSubDomains().get(i)).getHigh(0),SpaceBox<3,double>(dec.getSubDomains().get(i)).getLow(0));
			p1.get(1) = std::nextafter(SpaceBox<3,double>(dec.getSubDomains().get(i)).getHigh(1),SpaceBox<3,double>(dec.getSubDomains().get(i)).getLow(1));
			p1.get(2) = std::nextafter(SpaceBox<3,double>(dec.getSubDomains().get(i)).getHigh(2),SpaceBox<3,double>(dec.getSubDomains().get(i)).getLow(2));

			p2 = p1;

			p2.get(j) = std::nextafter(SpaceBox<3,double>(dec.getSubDomains().get(i)).getHigh(j),1.0);

			proc1 = dec.processorIDBC(p1);
			proc2 = dec.processorIDBC(p2);

			BOOST_REQUIRE(proc1 < vcl.size());
			BOOST_REQUIRE(proc2 < vcl.size());

			const openfpm::vector<std::pair<size_t, size_t>> & vp_id3 = dec.template ghost_processorID_pair<typename CartDecomposition<3, double>::lc_processor_id, typename CartDecomposition<3, double>::shift_id>(p1, UNIQUE);
			const openfpm::vector<std::pair<size_t, size_t>> & vp_id4 = dec.template ghost_processorID_pair<typename CartDecomposition<3, double>::lc_processor_id, typename CartDecomposition<3, double>::shift_id>(p2, UNIQUE);

			if (proc1 != proc2)
			{
				if (vcl.rank() == proc2)
				{
					BOOST_REQUIRE(vp_id4.size() != 0);
					BOOST_REQUIRE(vp_id3.size() == 0);
				}

				if (vcl.rank() == proc1)
				{
					BOOST_REQUIRE(vp_id4.size() == 0 );
					BOOST_REQUIRE(vp_id3.size() != 0 );
				}
			}

		}
	}
}



BOOST_AUTO_TEST_CASE( CartDecomposition_non_periodic_test_dist_grid)
{
	// Vcluster
	Vcluster<> & vcl = create_vcluster();

	CartDecomposition<3, float> dec(vcl);

	// Physical domain
	Box<3, float> box( { 0.0, 0.0, 0.0 }, { 1.0, 1.0, 1.0 });
	size_t div[3];
	size_t div_sub[3];

	// Get the number of processor and calculate the number of sub-domain
	// for each processor (SUB_UNIT_FACTOR=64)
	size_t n_proc = vcl.getProcessingUnits();
	size_t n_sub = n_proc * SUB_UNIT_FACTOR*4*4*4;

	// Set the number of sub-domains on each dimension (in a scalable way)
	for (int i = 0; i < 3; i++)
	{div[i] = openfpm::math::round_big_2(pow(n_sub,1.0/3));}

	// create a sub_distribution grid
	for (int i = 0; i < 3; i++)
	{div_sub[i] = div[i] / 4;}

	grid_sm<3,void> gsub(div_sub);

	// Define ghost
	Ghost<3, float> g(0.01);

	// Boundary conditions
	size_t bc[] = { NON_PERIODIC, NON_PERIODIC, NON_PERIODIC };

	// Decompose
	dec.setParameters(div,box,bc,g,gsub);
	dec.decompose();
	dec.write("Test_sub_dist2");

	// For each calculated ghost box
	for (size_t i = 0; i < dec.getNIGhostBox(); i++)
	{
		SpaceBox<3,float> b = dec.getIGhostBox(i);
		size_t proc = dec.getIGhostBoxProcessor(i);

		// sample one point inside the box
		Point<3,float> p = b.rnd();

		// Check that ghost_processorsID return that processor number
		const openfpm::vector<size_t> & pr = dec.ghost_processorID<CartDecomposition<3,float>::processor_id>(p);

		bool found = false;

		for (size_t j = 0; j < pr.size(); j++)
		{
			if (pr.get(j) == proc)
			{	found = true; break;}
		}

		if (found == false)
		{
			const openfpm::vector<size_t> pr2 = dec.ghost_processorID<CartDecomposition<3,float>::processor_id>(p);
		}

		BOOST_REQUIRE_EQUAL(found,true);
	}

	// Check the consistency

	bool val = dec.check_consistency();
	BOOST_REQUIRE_EQUAL(val,true);
}

BOOST_AUTO_TEST_CASE( CartDecomposition_nsub_algo_functions_test)
{
	size_t n_sub = 64*2;
	size_t div[3];

	nsub_to_div2<3>(div,n_sub,3);

	BOOST_REQUIRE_EQUAL(div[0],8ul);
	BOOST_REQUIRE_EQUAL(div[1],8ul);
	BOOST_REQUIRE_EQUAL(div[2],8ul);

	nsub_to_div2<3>(div,n_sub,2);

	BOOST_REQUIRE_EQUAL(div[0],16ul);
	BOOST_REQUIRE_EQUAL(div[1],16ul);
	BOOST_REQUIRE_EQUAL(div[2],1ul);

	nsub_to_div2<3>(div,n_sub,1);

	BOOST_REQUIRE_EQUAL(div[0],128ul);
	BOOST_REQUIRE_EQUAL(div[1],1ul);
	BOOST_REQUIRE_EQUAL(div[2],1ul);

	n_sub = 64*3;
	nsub_to_div<3>(div,n_sub,3);

	BOOST_REQUIRE_EQUAL(div[0],5ul);
	BOOST_REQUIRE_EQUAL(div[1],5ul);
	BOOST_REQUIRE_EQUAL(div[2],5ul);

	nsub_to_div<3>(div,n_sub,2);

	BOOST_REQUIRE_EQUAL(div[0],13ul);
	BOOST_REQUIRE_EQUAL(div[1],13ul);
	BOOST_REQUIRE_EQUAL(div[2],1ul);

	nsub_to_div<3>(div,n_sub,1);

	BOOST_REQUIRE_EQUAL(div[0],192ul);
	BOOST_REQUIRE_EQUAL(div[1],1ul);
	BOOST_REQUIRE_EQUAL(div[2],1ul);

	// Test high dimension cart decomposition subdivision

	Box<50,double> domain;
	size_t bc[50];
	Ghost<50,double> ghost(0.01);

	for(size_t i = 0 ; i < 50 ; i++)
	{
		domain.setLow(i,0.0);
		domain.setHigh(i,1.0);
		bc[i] = NON_PERIODIC;
	}

	CartDecomposition<50,double> dec(create_vcluster());

	dec.setGoodParameters(domain,bc,ghost,64);

	size_t div2[50];
	dec.getParameters(div2);

	auto & v_cl = create_vcluster();
	if (v_cl.size() == 1)
	{
		for (size_t i = 0 ; i < 50 ; i++)
		{
			if (i < 6)
			{BOOST_REQUIRE_EQUAL(div2[i],2ul);}
			else
			{BOOST_REQUIRE_EQUAL(div2[i],1ul);}
		}
	}

	if (v_cl.size() == 2)
	{
		for (size_t i = 0 ; i < 50 ; i++)
		{
			if (i < 7)
			{BOOST_REQUIRE_EQUAL(div2[i],2ul);}
			else
			{BOOST_REQUIRE_EQUAL(div2[i],1ul);}
		}
	}

	if (v_cl.size() == 3)
	{
		for (size_t i = 0 ; i < 50 ; i++)
		{
			if (i < 2)
			{BOOST_REQUIRE_EQUAL(div2[i],13ul);}
			else
			{BOOST_REQUIRE_EQUAL(div2[i],1ul);}
		}
	}

	if (v_cl.size() == 4)
	{
		for (size_t i = 0 ; i < 50 ; i++)
		{
			if (i < 8)
			{BOOST_REQUIRE_EQUAL(div2[i],2ul);}
			else
			{BOOST_REQUIRE_EQUAL(div2[i],1ul);}
		}
	}

	if (v_cl.size() == 5)
	{
		for (size_t i = 0 ; i < 50 ; i++)
		{
			if (i < 8)
			{BOOST_REQUIRE_EQUAL(div2[i],2ul);}
			else
			{BOOST_REQUIRE_EQUAL(div2[i],1ul);}
		}
	}

	if (v_cl.size() == 6)
	{
		for (size_t i = 0 ; i < 50 ; i++)
		{
			if (i < 3)
			{BOOST_REQUIRE_EQUAL(div2[i],7ul);}
			else
			{BOOST_REQUIRE_EQUAL(div2[i],1ul);}
		}
	}
}

BOOST_AUTO_TEST_SUITE_END()

