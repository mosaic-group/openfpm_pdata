/*
 * grid_dist_id_dlb_unit_test.cpp
 *
 *  Created on: May 4, 2018
 *      Author: i-bird
 */


#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "Point_test.hpp"
#include "Grid/grid_dist_id.hpp"
#include "data_type/aggregate.hpp"
#include "grid_dist_id_util_tests.hpp"
#include "Vector/vector_dist.hpp"

BOOST_AUTO_TEST_SUITE( grid_dist_id_dlb_test )

template<typename grid, typename vector>
void test_vector_grid_dlb()
{
	// Domain
	Box<3,float> domain3({0.0,0.0,0.0},{1.0,1.0,1.0});

	Ghost<3,long int> g(1);

	size_t sz[3] = {37,37,37};

	grid gdist(sz,domain3,g,DEC_GRAN(128));


	aggregate<long int,long int,long int> bck;
	bck.template get<0>() = -57;
	bck.template get<1>() = -90;
	bck.template get<2>() = -123;

	gdist.setBackgroundValue(bck);

	vector vd(gdist.getDecomposition(),0);

	// we fill the grid with a gaussian

	auto it = gdist.getGridIterator();

	while (it.isNext())
	{
		auto p = it.get_dist();
		auto gkey = it.get();

		gdist.template insert<0>(p) = gkey.get(0);
		gdist.template insert<1>(p) = gkey.get(1);
		gdist.template insert<2>(p) = gkey.get(2);

		++it;
	}

	// fill a sphere of particles

	auto it2 = vd.getGridIterator(sz);

	while (it2.isNext())
	{
		auto gkey = it2.get();

		vd.add();

		float x = 0.2*gkey.get(0)*it2.getSpacing(0) + 0.05;
		float y = 0.2*gkey.get(1)*it2.getSpacing(1) + 0.05;
		float z = 0.2*gkey.get(2)*it2.getSpacing(2) + 0.05;

		vd.getLastPos()[0] = x;
		vd.getLastPos()[1] = y;
		vd.getLastPos()[2] = z;

		++it2;
	}

	size_t n_step = 50;
	for (size_t i = 0; i < n_step ; i++)
	{
		vd.map();
		vd.addComputationCosts();
		vd.getDecomposition().decompose();
		vd.map();

		gdist.getDecomposition() = vd.getDecomposition();
		gdist.map();

		// Check

		bool check = true;
		auto it = gdist.getDomainIterator();


		while (it.isNext())
		{
			auto p = it.get();
			auto gkey = it.getGKey(p);

			check &= gdist.template get<0>(p) == gkey.get(0);
			check &= gdist.template get<1>(p) == gkey.get(1);
			check &= gdist.template get<2>(p) == gkey.get(2);

			if (check == false)
			{
				std::cout << "Error " << gdist.template get<0>(p) << " != " << gkey.get(0) << " " << gdist.template get<1>(p) << " != " << gkey.get(1) << " " << gdist.template get<2>(p) << " != " << gkey.get(2) << std::endl;
				break;
			}

			++it;
		}

		BOOST_REQUIRE_EQUAL(check,true);

		// Calculate shift vector

		double t2 = 6.28*(double)(i+1)/n_step;
		double t = 6.28*(double)i/n_step;
		double v[3];
		v[0] = 0.7*fabs(sin(t2)) - 0.7*fabs(sin(t));
		v[1] = 0.7*fabs(sin(1.7*t2)) - 0.7*fabs(sin(1.7*t));
		v[2] = 0.7*fabs(sin(2.5*t2)) - 0.7*fabs(sin(2.5*t));

		auto it2 = vd.getDomainIterator();

		while (it2.isNext())
		{
			auto p = it2.get();

			vd.getPos(p)[0] += v[0];
			vd.getPos(p)[1] += v[1];
			vd.getPos(p)[2] += v[2];

			++it2;
		}
		vd.map();
	}
}

struct GaussianDLB
{
	double t = 0.0;

	size_t resolution(const Point<3,float> & p)
	{
		float x0 = fabs(sin(t));
		float y0 = fabs(cos(t));
		float z0 = fabs(sin(3.0*t));

		return 100 * exp( - ((p.get(0) - x0)*(p.get(0) - x0) + (p.get(1) - y0)*(p.get(1) - y0) + (p.get(2) - z0)*(p.get(2) - z0)) / 0.3);
	}

	double distributionTol()
	{
		return 1.01;
	}
};

template<typename grid, typename vector>
void test_vector_grid_dlb_resolution()
{
	// Domain
	Box<3,float> domain3({0.0,0.0,0.0},{1.0,1.0,1.0});

	Ghost<3,long int> g(1);

	size_t sz[3] = {37,37,37};

	grid gdist(sz,domain3,g,DEC_GRAN(128));


	aggregate<long int,long int,long int> bck;
	bck.template get<0>() = -57;
	bck.template get<1>() = -90;
	bck.template get<2>() = -123;

	gdist.setBackgroundValue(bck);

	// we fill the grid with a gaussian

	auto it = gdist.getGridIterator();

	while (it.isNext())
	{
		auto p = it.get_dist();
		auto gkey = it.get();

		gdist.template insert<0>(p) = gkey.get(0);
		gdist.template insert<1>(p) = gkey.get(1);
		gdist.template insert<2>(p) = gkey.get(2);

		++it;
	}

	GaussianDLB gdlb;
	gdlb.t = 0.0;

	size_t n_step = 50;
	for (size_t i = 0; i < n_step ; i++)
	{
		gdlb.t = (float)i/n_step;
		gdist.addComputationCosts(gdlb);
		gdist.getDecomposition().decompose();
		gdist.map();

		gdist.write_frame("sgrid",i);

		// Check

		bool check = true;
		auto it = gdist.getDomainIterator();


		while (it.isNext())
		{
			auto p = it.get();
			auto gkey = it.getGKey(p);

			check &= gdist.template get<0>(p) == gkey.get(0);
			check &= gdist.template get<1>(p) == gkey.get(1);
			check &= gdist.template get<2>(p) == gkey.get(2);

			++it;
		}

		BOOST_REQUIRE_EQUAL(check,true);
	}
}

BOOST_AUTO_TEST_CASE( grid_dist_dlb_test )
{
	typedef sgrid_dist_id<3,float,aggregate<long int,long int,long int>> grid_sparse;
	typedef vector_dist<3,float,aggregate<long int, long int> > particles;

	test_vector_grid_dlb<grid_sparse,particles>();
}

BOOST_AUTO_TEST_CASE( grid_dist_dlb_test_resolution )
{
	typedef sgrid_dist_id<3,float,aggregate<long int,long int,long int>> grid_sparse;
	typedef vector_dist<3,float,aggregate<long int, long int> > particles;

	test_vector_grid_dlb_resolution<grid_sparse,particles>();
}

BOOST_AUTO_TEST_SUITE_END()

