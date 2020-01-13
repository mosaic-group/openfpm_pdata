/*
 * shift_vect_converter_tests.cpp
 *
 *  Created on: Feb 8, 2018
 *      Author: i-bird
 */

#include "config.h"

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include "Space/Shape/Box.hpp"

#include "Vector/map_vector.hpp"
#include "Decomposition/shift_vect_converter.hpp"

BOOST_AUTO_TEST_SUITE( shift_vect_converter_tests_suite )

BOOST_AUTO_TEST_CASE( shift_vect_converter_tests_use )
{
	{
	Box<3,double> domain({0.0,0.0,0.0},{1.0,1.0,1.0});
	shift_vect_converter<3,double,HeapMemory,memory_traits_lin> svc;
	size_t bc[3] = {PERIODIC,PERIODIC,PERIODIC};

	openfpm::vector<Point<3,double>> sv;

	svc.generateShiftVectors(domain,bc,sv);

	BOOST_REQUIRE_EQUAL(sv.size(),27ul);

	// We test that the cominations generate the correct shift vectors
	comb<3> cmb1({-1,-1,1});
	comb<3> cmb2({-1,0,1});
	comb<3> cmb3({0,0,1});

	size_t i = svc.linId(cmb1);

	BOOST_REQUIRE_EQUAL(sv.get<0>(i)[0],-1.0);
	BOOST_REQUIRE_EQUAL(sv.get<0>(i)[1],1.0);
	BOOST_REQUIRE_EQUAL(sv.get<0>(i)[2],1.0);

	i = svc.linId(cmb2);

	BOOST_REQUIRE_EQUAL(sv.get<0>(i)[0],-1.0);
	BOOST_REQUIRE_EQUAL(sv.get<0>(i)[1],0.0);
	BOOST_REQUIRE_EQUAL(sv.get<0>(i)[2],1.0);

	i = svc.linId(cmb3);

	BOOST_REQUIRE_EQUAL(sv.get<0>(i)[0],-1.0);
	BOOST_REQUIRE_EQUAL(sv.get<0>(i)[1],0.0);
	BOOST_REQUIRE_EQUAL(sv.get<0>(i)[2],0.0);

	}

	{
	openfpm::vector<Point<50,double>> sv;
	Box<50,double> domain;
	size_t bc[50];

	for (size_t i = 0 ; i < 50 ; i++)
	{
		domain.setLow(i,0.0);
		domain.setHigh(i,1.0);
		bc[i] = NON_PERIODIC;
	}

	bc[5] = PERIODIC;
	bc[17] = PERIODIC;
	bc[23] = PERIODIC;

	shift_vect_converter<50,double,HeapMemory,memory_traits_lin> svc;

	svc.generateShiftVectors(domain,bc,sv);

	BOOST_REQUIRE_EQUAL(sv.size(),27ul);

	// We test that the cominations generate the correct shift vectors
	comb<50> cmb1;
	comb<50> cmb2;
	comb<50> cmb3;

	cmb1.c[5] = 1;
	cmb1.c[17] = -1;
	cmb1.c[23] = -1;

	cmb2.c[5] = 1;
	cmb2.c[17] = 0;
	cmb2.c[23] = -1;

	cmb3.c[5] = 1;
	cmb3.c[17] = 0;
	cmb3.c[23] = 0;

	size_t i = svc.linId(cmb1);

	BOOST_REQUIRE_EQUAL(sv.get<0>(i)[5],-1.0);
	BOOST_REQUIRE_EQUAL(sv.get<0>(i)[6],0.0);
	BOOST_REQUIRE_EQUAL(sv.get<0>(i)[17],1.0);
	BOOST_REQUIRE_EQUAL(sv.get<0>(i)[23],1.0);
	BOOST_REQUIRE_EQUAL(sv.get<0>(i)[24],0.0);

	i = svc.linId(cmb2);

	BOOST_REQUIRE_EQUAL(sv.get<0>(i)[5],-1.0);
	BOOST_REQUIRE_EQUAL(sv.get<0>(i)[6],0.0);
	BOOST_REQUIRE_EQUAL(sv.get<0>(i)[17],0.0);
	BOOST_REQUIRE_EQUAL(sv.get<0>(i)[23],1.0);
	BOOST_REQUIRE_EQUAL(sv.get<0>(i)[24],0.0);

	i = svc.linId(cmb3);

	BOOST_REQUIRE_EQUAL(sv.get<0>(i)[5],-1.0);
	BOOST_REQUIRE_EQUAL(sv.get<0>(i)[6],0.0);
	BOOST_REQUIRE_EQUAL(sv.get<0>(i)[17],0.0);
	BOOST_REQUIRE_EQUAL(sv.get<0>(i)[23],0.0);
	BOOST_REQUIRE_EQUAL(sv.get<0>(i)[24],0.0);

	}
}


BOOST_AUTO_TEST_SUITE_END()

