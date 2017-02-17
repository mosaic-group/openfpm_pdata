/*
 * se_class3_vector_unit_tests.hpp
 *
 *  Created on: Feb 12, 2017
 *      Author: i-bird
 */

#ifndef SRC_VECTOR_SE_CLASS3_VECTOR_UNIT_TESTS_HPP_
#define SRC_VECTOR_SE_CLASS3_VECTOR_UNIT_TESTS_HPP_

#ifdef SE_CLASS3

BOOST_AUTO_TEST_SUITE( vector_dist_class3 )

BOOST_AUTO_TEST_CASE( vector_dist_class3_check )
{
	Box<2,float> domain({0.0,0.0},{1.0,1.0});

	// Here we define the boundary conditions of our problem
    size_t bc[2]={PERIODIC,PERIODIC};

	// extended boundary around the domain, and the processor domain
	Ghost<2,float> g(0.05);

	vector_dist<2,float,aggregate<float,float[3],openfpm::vector<double>>> vd(4096,domain,bc,g);

	bool error = false;

	{
	auto it = vd.getDomainIterator();
	try
	{
		while (it.isNext())
		{
			auto p = it.get();

			// Read something not initialized ERROR
			Point<2,float> a = vd.getPosRead(p);

			// Suppress compiler error
			a.get(0) = 0;

			++it;
		}
	}
	catch (std::exception & e)
	{
		error = true;
	}
	}

	BOOST_REQUIRE_EQUAL(error,true);
	error = false;

	// Initialize

	{
	auto it = vd.getDomainIterator();

	while (it.isNext())
	{
		auto p = it.get();

		// Read something not initialized ERROR
		vd.getPosWrite(p)[0] = (double)rand()/RAND_MAX;
		vd.getPosWrite(p)[1] = (double)rand()/RAND_MAX;

		++it;
	}

	}

	{
	auto it = vd.getDomainIterator();
	try
	{
		while (it.isNext())
		{
			auto p = it.get();

			// Now has bee initialized so not error should be reported
			Point<2,float> a = vd.getPosRead(p);

			// Suppress compiler error
			a.get(0) = 0;

			++it;
		}
	}
	catch (std::exception & e)
	{
		error = true;
	}
	}

	BOOST_REQUIRE_EQUAL(error,false);
	error = false;

	// we redistribute the particle

	vd.map();

	try
	{
		// Error we are getting the cell list from ghost out of sync in position
		vd.getCellList(0.05);
	}
	catch (std::exception & e)
	{
		error = true;
	}

	BOOST_REQUIRE_EQUAL(error,true);
	error = false;

	try
	{
		// Error we are getting the cell list from ghost out of sync in position
		vd.getVerlet(0.1);
	}
	catch (std::exception & e)
	{
		error = true;
	}

	BOOST_REQUIRE_EQUAL(error,true);
	error = false;

	try
	{
		// Error we are getting the cell list from ghost out of sync in position
		vd.getCellListSym(0.05);
	}
	catch (std::exception & e)
	{
		error = true;
	}

	BOOST_REQUIRE_EQUAL(error,true);
	error = false;

	vd.ghost_get<0>();

	try
	{
		// OK
		vd.getCellList(0.05);
	}
	catch (std::exception & e)
	{
		error = true;
	}

	BOOST_REQUIRE_EQUAL(error,false);
	error = false;

	// We make dirty a ghost particle
	vd.template getPropWrite<0>(vd.size_local()) = 0.5;

	try
	{
		// Error we are destroying information
		vd.ghost_get<0>();
	}
	catch (std::exception & e)
	{
		error = true;
	}

	BOOST_REQUIRE_EQUAL(error,true);
	error = false;

	// We make dirty a ghost particle
	vd.template getPropWrite<0>(vd.size_local()) = 0.5;

	try
	{
		// Also error we are also destroying information
		vd.ghost_get<1>();
	}
	catch (std::exception & e)
	{
		error = true;
	}

	BOOST_REQUIRE_EQUAL(error,true);
	error = false;

	// We make dirty a ghost particle
	vd.template getPropWrite<0>(vd.size_local()) = 0.5;

	try
	{
		// OK we are not destroying information
		vd.ghost_get<1>(KEEP_PROPERTIES);
	}
	catch (std::exception & e)
	{
		error = true;
	}

	BOOST_REQUIRE_EQUAL(error,false);
	error = false;

	// We make dirty a ghost particle
	vd.template getPropWrite<0>(vd.size_local()) = 0.5;

	try
	{
		// Error we are destroying information
		vd.ghost_get<0>(KEEP_PROPERTIES);
	}
	catch (std::exception & e)
	{
		error = true;
	}

	BOOST_REQUIRE_EQUAL(error,true);
	error = false;

	try
	{
		// error property 0 has never been initialized
		vd.ghost_put<add_,0>();
	}
	catch (std::exception & e)
	{
		error = true;
	}

	BOOST_REQUIRE_EQUAL(error,true);
	error = false;

	{
	auto it = vd.getDomainIterator();
	try
	{
		while (it.isNext())
		{
			auto p = it.get();

			vd.getPropWrite<0>(p) = 2.0;

			++it;
		}

		// OK Property 0 has been initialized
		vd.ghost_put<add_,0>();
	}
	catch (std::exception & e)
	{
		error = true;
	}
	}

	BOOST_REQUIRE_EQUAL(error,false);
	error = false;

	try
	{
		// OK we are not destroying information
		vd.ghost_get<0>();
	}
	catch (std::exception & e)
	{
		error = true;
	}

	BOOST_REQUIRE_EQUAL(error,false);
	error = false;

	try
	{
		// Error we are getting the cell list from ghost out of sync in position
		vd.getCellList(0.05);
	}
	catch (std::exception & e)
	{
		error = true;
	}

	BOOST_REQUIRE_EQUAL(error,false);
	error = false;

	BOOST_REQUIRE_EQUAL(vd.get_se_class3().isGhostSync<1>(),NOTSYNC);

	vd.ghost_put<add_,0>();
	vd.ghost_get<1>();

	auto NN = vd.getCellList(0.05);

	BOOST_REQUIRE_EQUAL(vd.get_se_class3().isGhostSync<1>(),SYNC);


	auto it = vd.getDomainIterator();

	while (it.isNext())
	{
		auto p = it.get();

		// we set property 1
		vd.getPropWrite<1>(p)[0] = 1.0;
		vd.getPropWrite<1>(p)[1] = 1.0;
		vd.getPropWrite<1>(p)[2] = 1.0;

		++it;
	}

	BOOST_REQUIRE_EQUAL(vd.get_se_class3().isGhostSync<1>(),NOTSYNC);

	{
	auto it = vd.getDomainIterator();
	try
	{
		while (it.isNext())
		{
			auto p = it.get();

			Point<2,float> xp = vd.getPosRead(p);

			auto NNp = NN.template getNNIterator<NO_CHECK>(NN.getCell(xp));

			while (NNp.isNext())
			{
				auto q = NNp.get();

				// Error ghost is not initialized
				Point<3,float> xq = vd.template getPropRead<1>(q);

				xq.get(0) = 0.0;

				++NNp;
			}

			++it;
		}
	}
	catch (std::exception & e)
	{
		error = true;
	}
	}

	BOOST_REQUIRE_EQUAL(error,true);
	error = false;

	vd.ghost_get<1>();

	{
	auto it = vd.getDomainIterator();

	while (it.isNext())
	{
		auto p = it.get();

		// we set property 1
		vd.getPropWrite<1>(p)[0] = 1.0;
		vd.getPropWrite<1>(p)[1] = 1.0;
		vd.getPropWrite<1>(p)[2] = 1.0;

		++it;
	}
	}

	BOOST_REQUIRE_EQUAL(vd.get_se_class3().isGhostSync<1>(),NOTSYNC);

	{
	auto it = vd.getDomainIterator();
	try
	{
		while (it.isNext())
		{
			auto p = it.get();

			Point<2,float> xp = vd.getPosRead(p);

			auto NNp = NN.template getNNIterator<NO_CHECK>(NN.getCell(xp));

			while (NNp.isNext())
			{
				auto q = NNp.get();

				// Error ghost is not initialized
				Point<3,float> xq = vd.template getPropRead<1>(q);

				xq.get(0) = 0.0;

				++NNp;
			}

			++it;
		}
	}
	catch (std::exception & e)
	{
		error = true;
	}
	}

	BOOST_REQUIRE_EQUAL(error,true);
	error = false;

	vd.ghost_get<1>();

	{
	auto it = vd.getDomainIterator();
	try
	{
		while (it.isNext())
		{
			auto p = it.get();

			Point<2,float> xp = vd.getPosRead(p);

			auto NNp = NN.template getNNIterator<NO_CHECK>(NN.getCell(xp));

			while (NNp.isNext())
			{
				auto q = NNp.get();

				// Error we forgot ghost_get
				Point<3,float> xq = vd.template getPropRead<1>(q);

				xq.get(0) = 0.0;

				++NNp;
			}

			++it;
		}
	}
	catch (std::exception & e)
	{
		error = true;
	}
	}

	BOOST_REQUIRE_EQUAL(error,false);
}



///////////////////////////////////////////// Add and Remove test


BOOST_AUTO_TEST_CASE( vector_dist_class3_check_add )
{
	Box<2,float> domain({0.0,0.0},{1.0,1.0});

	// Here we define the boundary conditions of our problem
    size_t bc[2]={PERIODIC,PERIODIC};

	// extended boundary around the domain, and the processor domain
	Ghost<2,float> g(0.05);

	vector_dist<2,float,aggregate<float,float[3],openfpm::vector<double>>> vd(0,domain,bc,g);

	bool error = false;

	// Initialize

	{
	for (size_t i = 0 ; i < 1200 ; i++)
	{
		vd.add();

		// Read something not initialized ERROR
		vd.getLastPos()[0] = (double)rand()/RAND_MAX;
		vd.getLastPos()[1] = (double)rand()/RAND_MAX;
	}
	}

	{
	auto it = vd.getDomainIterator();
	try
	{
		while (it.isNext())
		{
			auto p = it.get();

			// Now has bee initialized so not error should be reported
			Point<2,float> a = vd.getPosRead(p);

			//Suppress warnings
			a.get(0) = 0.0;

			++it;
		}
	}
	catch (std::exception & e)
	{
		error = true;
	}
	}

	BOOST_REQUIRE_EQUAL(error,true);
	error = false;

	vd.map();

	{
	auto it = vd.getDomainIterator();
	try
	{
		while (it.isNext())
		{
			auto p = it.get();

			// Now has bee initialized so not error should be reported
			float a = vd.getPropRead<0>(p);

			// Suppress warnings
			a = 0;
			a++;

			++it;
		}
	}
	catch (std::exception & e)
	{
		error = true;
	}
	}

	BOOST_REQUIRE_EQUAL(error,true);
}

#if defined(CHECKFOR_POSNAN) && defined(CHECKFOR_PROPNAN) && defined(CHECKFOR_POSINF) && defined(CHECKFOR_PROPINF)


BOOST_AUTO_TEST_CASE( vector_dist_class3_check_nan_inf )
{
	Box<2,float> domain({0.0,0.0},{1.0,1.0});

	// Here we define the boundary conditions of our problem
    size_t bc[2]={PERIODIC,PERIODIC};

	// extended boundary around the domain, and the processor domain
	Ghost<2,float> g(0.05);

	vector_dist<2,float,aggregate<float,float[3],openfpm::vector<double>>> vd(0,domain,bc,g);

	bool error = false;

	// Initialize

	{
	for (size_t i = 0 ; i < 1200 ; i++)
	{
		vd.add();

		// Read something not initialized ERROR
		vd.getLastPos()[0] = 0.0/0.0;
		vd.getLastPos()[1] = 0.0/0.0;
	}
	}

	{
	try
	{
		auto it = vd.getDomainIterator();
	}
	catch (std::exception & e)
	{
		error = true;
	}
	}

	BOOST_REQUIRE_EQUAL(error,true);
	error = false;

	vector_dist<2,float,aggregate<float,float[3],openfpm::vector<double>>> vd2(0,domain,bc,g);

	{
	for (size_t i = 0 ; i < 1200 ; i++)
	{
		vd2.add();

		// Read something not initialized ERROR
		vd2.getLastPos()[0] = 5.0/0.0;
		vd2.getLastPos()[1] = 5.0/0.0;
	}
	}

	{
	try
	{
		auto it = vd2.getDomainIterator();
	}
	catch (std::exception & e)
	{
		error = true;
	}
	}

	BOOST_REQUIRE_EQUAL(error,true);
}

#endif

BOOST_AUTO_TEST_SUITE_END()

#endif

#endif /* SRC_VECTOR_SE_CLASS3_VECTOR_UNIT_TESTS_HPP_ */
