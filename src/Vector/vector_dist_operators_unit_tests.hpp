/*
 * vector_dist_operators_unit_tests.hpp
 *
 *  Created on: Jun 11, 2016
 *      Author: i-bird
 */

#ifndef SRC_VECTOR_VECTOR_DIST_OPERATORS_UNIT_TESTS_HPP_
#define SRC_VECTOR_VECTOR_DIST_OPERATORS_UNIT_TESTS_HPP_

constexpr int A = 0;
constexpr int B = 1;
constexpr int C = 2;

constexpr int PA = 3;
constexpr int PB = 4;
constexpr int PC = 5;

//////////////////// Here we define all the function to checl the operators

template <unsigned int prp, typename vector> bool check_values(const vector & v,float a)
{
	bool ret = true;
	auto it = v.getDomainIterator();

	while (it.isNext())
	{
		ret &= v.template getProp<prp>(it.get()) == a;

		++it;
	}

	BOOST_REQUIRE_EQUAL(ret,true);

	return ret;
}

template <typename vector> bool check_values_complex_expr(const vector & vd)
{
	bool ret = true;
	auto it = vd.getDomainIterator();

	while (it.isNext())
	{
		auto key = it.get();

		float base1 = vd.template getProp<B>(key) + 2.0 + vd.template getProp<B>(key) - 2.0*vd.template getProp<C>(key) / 5.0;
		float base2 = vd.template getProp<A>(key);

		ret &=  base1 == base2;

		++it;
	}

	BOOST_REQUIRE_EQUAL(ret,true);

	return ret;
}

template <typename vector> bool check_values_sum(const vector & vd, double d)
{
	bool ret = true;
	auto it = vd.getDomainIterator();

	while (it.isNext())
	{
		auto key = it.get();

		float base1 = vd.template getProp<B>(key) + d;
		float base2 = vd.template getProp<A>(key);

		ret &=  base1 == base2;

		++it;
	}

	BOOST_REQUIRE_EQUAL(ret,true);

	return ret;
}

template <typename vector> bool check_values_sum(const vector & vd1, const vector & vd2)
{
	bool ret = true;
	auto it = vd1.getDomainIterator();

	while (it.isNext())
	{
		auto key = it.get();

		float base1 = vd1.template getProp<B>(key) + vd2.template getProp<C>(key);
		float base2 = vd1.template getProp<A>(key);

		ret &=  base1 == base2;

		++it;
	}

	BOOST_REQUIRE_EQUAL(ret,true);

	return ret;
}

template <typename vector> bool check_values_sum_3(const vector & vd1)
{
	bool ret = true;
	auto it = vd1.getDomainIterator();

	while (it.isNext())
	{
		auto key = it.get();

		float base1 = vd1.template getProp<B>(key) + vd1.template getProp<C>(key) + vd1.template getProp<B>(key);
		float base2 = vd1.template getProp<A>(key);

		ret &=  base1 == base2;

		++it;
	}

	BOOST_REQUIRE_EQUAL(ret,true);

	return ret;
}

template <typename vector> bool check_values_sum_4(const vector & vd1)
{
	bool ret = true;
	auto it = vd1.getDomainIterator();

	while (it.isNext())
	{
		auto key = it.get();

		float base1 = vd1.template getProp<B>(key) + vd1.template getProp<C>(key) + vd1.template getProp<B>(key) + vd1.template getProp<C>(key);
		float base2 = vd1.template getProp<A>(key);

		ret &=  base1 == base2;

		++it;
	}

	BOOST_REQUIRE_EQUAL(ret,true);

	return ret;
}

template <typename vector> bool check_values_sub(const vector & vd, double d)
{
	bool ret = true;
	auto it = vd.getDomainIterator();

	while (it.isNext())
	{
		auto key = it.get();

		float base1 = vd.template getProp<B>(key) - d;
		float base2 = vd.template getProp<A>(key);

		ret &=  base1 == base2;

		++it;
	}

	BOOST_REQUIRE_EQUAL(ret,true);

	return ret;
}

template <typename vector> bool check_values_sub(double d, const vector & vd)
{
	bool ret = true;
	auto it = vd.getDomainIterator();

	while (it.isNext())
	{
		auto key = it.get();

		float base1 = d - vd.template getProp<B>(key);
		float base2 = vd.template getProp<A>(key);

		ret &=  base1 == base2;

		++it;
	}

	BOOST_REQUIRE_EQUAL(ret,true);

	return ret;
}

template <typename vector> bool check_values_sub(const vector & vd1, const vector & vd2)
{
	bool ret = true;
	auto it = vd1.getDomainIterator();

	while (it.isNext())
	{
		auto key = it.get();

		float base1 = vd1.template getProp<C>(key) - vd2.template getProp<B>(key);
		float base2 = vd1.template getProp<A>(key);

		ret &=  base1 == base2;

		++it;
	}

	BOOST_REQUIRE_EQUAL(ret,true);

	return ret;
}

template <typename vector> bool check_values_mul(const vector & vd, double d)
{
	bool ret = true;
	auto it = vd.getDomainIterator();

	while (it.isNext())
	{
		auto key = it.get();

		float base1 = vd.template getProp<B>(key) * d;
		float base2 = vd.template getProp<A>(key);

		ret &=  base1 == base2;

		++it;
	}

	BOOST_REQUIRE_EQUAL(ret,true);

	return ret;
}


template <typename vector> bool check_values_mul(const vector & vd1, const vector & vd2)
{
	bool ret = true;
	auto it = vd1.getDomainIterator();

	while (it.isNext())
	{
		auto key = it.get();

		float base1 = vd1.template getProp<C>(key) * vd2.template getProp<B>(key);
		float base2 = vd1.template getProp<A>(key);

		ret &=  base1 == base2;

		++it;
	}

	BOOST_REQUIRE_EQUAL(ret,true);

	return ret;
}


template <typename vector> bool check_values_div(const vector & vd, double d)
{
	bool ret = true;
	auto it = vd.getDomainIterator();

	while (it.isNext())
	{
		auto key = it.get();

		float base1 = vd.template getProp<B>(key) / d;
		float base2 = vd.template getProp<A>(key);

		ret &=  base1 == base2;

		++it;
	}

	BOOST_REQUIRE_EQUAL(ret,true);

	return ret;
}

template <typename vector> bool check_values_div(double d, const vector & vd)
{
	bool ret = true;
	auto it = vd.getDomainIterator();

	while (it.isNext())
	{
		auto key = it.get();

		float base1 = d / vd.template getProp<B>(key);
		float base2 = vd.template getProp<A>(key);

		ret &=  base1 == base2;

		++it;
	}

	BOOST_REQUIRE_EQUAL(ret,true);

	return ret;
}

template <typename vector> bool check_values_div(const vector & vd1, const vector & vd2)
{
	bool ret = true;
	auto it = vd1.getDomainIterator();

	while (it.isNext())
	{
		auto key = it.get();

		float base1 = vd1.template getProp<C>(key) / vd2.template getProp<B>(key);
		float base2 = vd1.template getProp<A>(key);

		ret &=  base1 == base2;

		++it;
	}

	BOOST_REQUIRE_EQUAL(ret,true);

	return ret;
}

template <typename vector> void fill_values(const vector & v)
{
	auto it = v.getDomainIterator();

	while (it.isNext())
	{
		auto p = it.get();

		v.template getProp<A>(p) = p.getKey()+1;
		v.template getProp<B>(p) = 2.0*p.getKey()+1;
		v.template getProp<C>(p) = 3.0*p.getKey()+1;

		++it;
	}
}

float getProp0(vector_dist<3,float,aggregate<float,float,float,VectorS<3,float>,VectorS<3,float>,VectorS<3,float>>> & v, size_t p)
{
	return v.template getProp<A>(p);
}

////////////////////////////////////////////////////////////////////////////

BOOST_AUTO_TEST_SUITE( vector_dist_test )


BOOST_AUTO_TEST_CASE( vector_dist_operators_test )
{
	if (create_vcluster().getProcessingUnits() > 3)
		return;

	Box<3,float> box({0.0,0.0,0.0},{1.0,1.0,1.0});

	// Boundary conditions
	size_t bc[3]={PERIODIC,PERIODIC,PERIODIC};

	// ghost
	Ghost<3,float> ghost(0.05);

	vector_dist<3,float,aggregate<float,float,float,VectorS<3,float>,VectorS<3,float>,VectorS<3,float>>> vd(100,box,bc,ghost);

	vd.getV<A>() = 1.0;
	vd.getV<B>() = 2.0f;
	vd.getV<C>() = 3.0;

	check_values<A>(vd,1.0);
	check_values<B>(vd,2.0);
	check_values<C>(vd,3.0);

	fill_values(vd);

	vd.getV<A>() = vd.getV<B>() + 2.0 + vd.getV<B>() - 2.0*vd.getV<C>() / 5.0;
	check_values_complex_expr(vd);

	// Various combination of 2 operator

	vd.getV<A>() = vd.getV<B>() + 2.0;
	check_values_sum(vd,2.0);
	vd.getV<A>() = 2.0 + vd.getV<B>();
	check_values_sum(vd,2.0);
	vd.getV<A>() = vd.getV<C>() + vd.getV<B>();
	check_values_sum(vd,vd);

	vd.getV<A>() = vd.getV<B>() - 2.0;
	check_values_sub(vd,2.0);
	vd.getV<A>() = 2.0 - vd.getV<B>();
	check_values_sub(2.0,vd);
	vd.getV<A>() = vd.getV<C>() - vd.getV<B>();
	check_values_sub(vd,vd);

	vd.getV<A>() = vd.getV<B>() * 2.0;
	check_values_mul(vd,2.0);
	vd.getV<A>() = 2.0 * vd.getV<B>();
	check_values_mul(vd,2.0);
	vd.getV<A>() = vd.getV<C>() * vd.getV<B>();
	check_values_mul(vd,vd);

	vd.getV<A>() = vd.getV<B>() / 2.0;
	check_values_div(vd,2.0);
	vd.getV<A>() = 2.0 / vd.getV<B>();
	check_values_div(vd,2.0);
	vd.getV<A>() = vd.getV<C>() / vd.getV<B>();
	check_values_div(vd,vd);

	// Variuos combination 3 operator

	vd.getV<A>() = vd.getV<B>() + (vd.getV<C>() + vd.getV<B>());
	check_values_sum_3(vd);
	vd.getV<A>() = (vd.getV<C>() + vd.getV<B>()) + vd.getV<B>();
	check_values_sum_3(vd);
	vd.getV<A>() = (vd.getV<C>() + vd.getV<B>()) + (vd.getV<C>() + vd.getV<B>());
	check_values_sum_4(vd);

/*	vd.getV<A>() = vd.getV<B>() - (vd.getV<C>() + vd.getV<B>());
	vd.getV<A>() = (vd.getV<C>() + vd.getV<B>()) - vd.getV<B>();
	vd.getV<A>() = (vd.getV<C>() + vd.getV<B>()) - (vd.getV<C>() + vd.getV<B>());

	vd.getV<A>() = vd.getV<B>() * (vd.getV<C>() + vd.getV<B>());
	vd.getV<A>() = (vd.getV<C>() + vd.getV<B>()) * vd.getV<B>();
	vd.getV<A>() = (vd.getV<C>() + vd.getV<B>()) * (vd.getV<C>() + vd.getV<B>());

	vd.getV<A>() = vd.getV<B>() / (vd.getV<C>() + vd.getV<B>());
	vd.getV<A>() = (vd.getV<C>() + vd.getV<B>()) / vd.getV<B>();
	vd.getV<A>() = (vd.getV<C>() + vd.getV<B>()) / (vd.getV<C>() + vd.getV<B>());*/

}

BOOST_AUTO_TEST_SUITE_END()



#endif /* SRC_VECTOR_VECTOR_DIST_OPERATORS_UNIT_TESTS_HPP_ */
