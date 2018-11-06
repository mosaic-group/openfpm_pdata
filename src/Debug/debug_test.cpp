#include "config.h"
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include "Vector/vector_dist.hpp"
#include "debug.hpp"

BOOST_AUTO_TEST_SUITE( debug_util_test )

BOOST_AUTO_TEST_CASE( debug_dist_util_test_use )
{
	Vcluster<> & v_cl = create_vcluster();

	if (v_cl.getProcessingUnits() != 1)
	{return;}

	// Boundary conditions
	size_t bc[3]={PERIODIC,PERIODIC,PERIODIC};

	// Box
	Box<3,float> box({0.0,0.0,0.0},{1.0,1.0,1.0});

	// ghost
	Ghost<3,float> ghost(0.1);

	// first phase
	vector_dist<3,float, aggregate<double,double>> vd(4,box,bc,ghost);

	vd.getPos(0)[0] = 0.1;
	vd.getPos(0)[1] = 0.2;
	vd.getPos(0)[2] = 0.3;

	vd.getPos(1)[0] = 0.2;
	vd.getPos(1)[1] = 0.3;
	vd.getPos(1)[2] = 0.4;

	vd.getPos(2)[0] = 0.3;
	vd.getPos(2)[1] = 0.4;
	vd.getPos(2)[2] = 0.5;

	vd.getPos(3)[0] = 0.4;
	vd.getPos(3)[1] = 0.5;
	vd.getPos(3)[2] = 0.6;

	Box<3,float> box1({0.0,0.0,0.0},{0.01,0.01,0.01});
	Box<3,float> box2({0.0,0.0,0.0},{0.5,0.5,0.5});

	bool test = debug_is_in_box(vd,box1);
	BOOST_REQUIRE_EQUAL(test,false);

	test = debug_is_in_box(vd,box2,debug_iterator::DOMAIN_IT,debug_run::HOST,false);
	BOOST_REQUIRE_EQUAL(test,true);
}

BOOST_AUTO_TEST_CASE( debug_util_test_use )
{
	Vcluster<> & v_cl = create_vcluster();

	if (v_cl.getProcessingUnits() != 1)
	{return;}

	// Boundary conditions
	size_t bc[3]={PERIODIC,PERIODIC,PERIODIC};

	// Box
	Box<3,float> box({0.0,0.0,0.0},{1.0,1.0,1.0});

	// ghost
	Ghost<3,float> ghost(0.1);

	// first phase
	openfpm::vector<Point<3,float>> vd(4);

	vd.get<0>(0)[0] = 0.1;
	vd.get<0>(0)[1] = 0.2;
	vd.get<0>(0)[2] = 0.3;

	vd.get<0>(1)[0] = 0.2;
	vd.get<0>(1)[1] = 0.3;
	vd.get<0>(1)[2] = 0.4;

	vd.get<0>(2)[0] = 0.3;
	vd.get<0>(2)[1] = 0.4;
	vd.get<0>(2)[2] = 0.5;

	vd.get<0>(3)[0] = 0.4;
	vd.get<0>(3)[1] = 0.5;
	vd.get<0>(3)[2] = 0.6;

	Box<3,float> box1({0.0,0.0,0.0},{0.01,0.01,0.01});
	Box<3,float> box2({0.0,0.0,0.0},{0.5,0.5,0.5});

	bool test = debug_is_in_box_single(vd,box1,vd.size());
	BOOST_REQUIRE_EQUAL(test,false);

	test = debug_is_in_box_single(vd,box2,vd.size(),"",debug_iterator::DOMAIN_IT,debug_run::HOST,false);
	BOOST_REQUIRE_EQUAL(test,true);
}

BOOST_AUTO_TEST_SUITE_END()

