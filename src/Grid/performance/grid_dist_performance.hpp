/*
 * grid_dist_performance.hpp
 *
 *  Created on: Jun 27, 2017
 *      Author: i-bird
 */

#ifndef SRC_GRID_GRID_DIST_PERFORMANCE_HPP_
#define SRC_GRID_GRID_DIST_PERFORMANCE_HPP_

#include "../../../openfpm_numerics/src/interpolation/interpolation.hpp"
#include "Grid/grid_dist_id.hpp"
#include "Plot/GoogleChart.hpp"
#include "interpolation/mp4_kernel.hpp"

#define GRID_ITERATOR_TESTS  30
#define GRID_INTERPOLATION_TESTS  30

// Vectors to store the data for 3D
openfpm::vector<double> time_iterator_normal_mean;
openfpm::vector<double> time_iterator_normal_dev;
openfpm::vector<double> time_iterator_stencil_mean;
openfpm::vector<double> time_iterator_stencil_dev;
openfpm::vector<double> time_inte_p2m_mean;
openfpm::vector<double> time_inte_m2p_mean;
openfpm::vector<double> time_inte_p2m_dev;
openfpm::vector<double> time_inte_m2p_dev;
openfpm::vector<size_t> nk_grid_st;
openfpm::vector<size_t> nk_grid_int;

BOOST_AUTO_TEST_SUITE( grid_iterator_performance_test )

constexpr int x = 0;
constexpr int y = 1;
constexpr int z = 2;

/*! \brief
 *
 *
 */
void grid_interpolation_benchmark(openfpm::vector<size_t> & nk_grid,
		                            openfpm::vector<double> & time_interpolation_p2m_mean,
									openfpm::vector<double> & time_interpolation_m2p_mean,
									openfpm::vector<double> & time_interpolation_p2m_dev,
									openfpm::vector<double> & time_interpolation_m2p_dev)
{
	for (size_t k = 0 ; k < nk_grid.size() ; k++)
	{
		size_t np = nk_grid.get(k);
		size_t sz[] = {np,np,np};

		Ghost<3,long int> gg(3);

		Box<3,float> domain({0.0,0.0,0.0},{1.0,1.0,1.0});

		size_t tot_part = np*np*np;

		grid_dist_id<3, float, aggregate<float>> gd(sz,domain,gg);
		vector_dist<3,float,aggregate<float>> vd(gd.getDecomposition(),tot_part);

		// Fill vd with particles

		auto vd_it = vd.getDomainIterator();

		while (vd_it.isNext())
		{
			auto p = vd_it.get();


			vd.getPos(p)[0] = ((float)std::rand()) / RAND_MAX;
			vd.getPos(p)[1] = ((float)std::rand()) / RAND_MAX;
			vd.getPos(p)[2] = ((float)std::rand()) / RAND_MAX;

			vd.template getProp<0>(p) = 1.0;

			++vd_it;
		}

		vd.map();

		double mean;
		double dev;
		openfpm::vector<double> measures;
		for (size_t j = 0 ; j < GRID_INTERPOLATION_TESTS ; j++)
		{

			interpolate<decltype(vd),decltype(gd),mp4_kernel<float>> inte(vd,gd);

			timer tstl;
			tstl.start();

			inte.p2m<0,0>(vd,gd);

			tstl.stop();
			measures.add(tstl.getwct());
		}
		standard_deviation(measures,mean,dev);
		time_interpolation_p2m_mean.add(mean);
		time_interpolation_p2m_dev.add(dev);

		std::cout << "Time particles to mesh " << time_interpolation_p2m_mean.last() << std::endl;

		for (size_t j = 0 ; j < GRID_INTERPOLATION_TESTS ; j++)
		{

			interpolate<decltype(vd),decltype(gd),mp4_kernel<float>> inte(vd,gd);

			timer tstl;
			tstl.start();

			inte.m2p<0,0>(gd,vd);

			tstl.stop();
			measures.add(tstl.getwct());
		}
		standard_deviation(measures,mean,dev);
		time_interpolation_m2p_mean.add(mean);
		time_interpolation_m2p_dev.add(dev);

		std::cout << "Time mesh to particles " << time_interpolation_m2p_mean.last() << std::endl;
	}
}

/*! \brief
 *
 *
 */
double grid_iterator_benchmark_stencil(grid_dist_id<3, float, aggregate<long int>, CartDecomposition<3,float>> & g_dist, double & total)
{
	grid_key_dx<3> star_stencil_3D[7] = {{0,0,0},
	                                         {-1,0,0},
											 {1,0,0},
											 {0,-1,0},
											 {0,1,0},
											 {0,0,-1},
											 {0,0,1}};

	timer tstl;
	tstl.start();

	auto st_it = g_dist.getDomainIteratorStencil(star_stencil_3D);

	while (st_it.isNext())
	{

		total+= 6*g_dist.template get<0>(st_it.getStencil<0>()) -
					 g_dist.template get<0>(st_it.getStencil<1>()) -
					 g_dist.template get<0>(st_it.getStencil<2>()) -
					 g_dist.template get<0>(st_it.getStencil<3>()) -
					 g_dist.template get<0>(st_it.getStencil<4>()) -
					 g_dist.template get<0>(st_it.getStencil<5>()) -
					 g_dist.template get<0>(st_it.getStencil<6>());

		++st_it;
	}

	tstl.stop();
	return tstl.getwct();
}

double grid_iterator_benchmark_norm(grid_dist_id<3, float, aggregate<long int>, CartDecomposition<3,float>> & g_dist, double & total)
{
	timer tnorm;
	tnorm.start();

	auto norm_it = g_dist.getDomainIterator();

	while (norm_it.isNext())
	{
		// center point
		auto key = norm_it.get();

		total+= 6*g_dist.template get<0>(key) -
					 g_dist.template get<0>(key.move(x,-1)) -
					 g_dist.template get<0>(key.move(x,1)) -
					 g_dist.template get<0>(key.move(y,-1)) -
					 g_dist.template get<0>(key.move(y,1)) -
					 g_dist.template get<0>(key.move(z,-1)) -
					 g_dist.template get<0>(key.move(z,1));

		++norm_it;
	}

	tnorm.stop();
	return tnorm.getwct();
}

/*! \brief Function for verlet test without an Hilbert curve reordering (unordered positioning)
 *
 */
template<unsigned int dim> void grid_iterator_benchmark(openfpm::vector<size_t> & nk_grid,
														   openfpm::vector<double> & time_iterator_normal_mean,
														   openfpm::vector<double> & time_iterator_stencil_mean,
														   openfpm::vector<double> & time_iterator_normal_dev,
														   openfpm::vector<double> & time_iterator_stencil_dev)
{
	std::string str("Testing " + std::to_string(dim) + "D grid iterator stencil and normal");
	print_test_v(str);

	{
		//For different grid sizes
		for (size_t i = 0; i < nk_grid.size(); i++ )
		{
			size_t sz[dim];

			//Number of particles
			size_t k = nk_grid.get(i);


			BOOST_TEST_CHECKPOINT( "Testing " << dim << "D vector grid iterator performance k=" << k );

			Box<dim,float> box;

			for (size_t i = 0; i < dim; i++)
			{
				box.setLow(i,0.0);
				box.setHigh(i,1.0);
				sz[i] = k;
			}

			Ghost<3,long int> g(1);

			// Distributed grid with id decomposition
			grid_dist_id<3, float, aggregate<long int>, CartDecomposition<3,float>> g_dist(sz,box,g);

			// fill the grid with values

			auto it = g_dist.getDomainGhostIterator();

			while (it.isNext())
			{
				auto p = it.get();
				auto gkey = it.getGKey(p);

				g_dist.template get<0>(p) = gkey.get(0) + gkey.get(1) + gkey.get(2);

				++it;
			}

			g_dist.ghost_get<0>();

			double total = 0;

			double mean;
			double dev;
			openfpm::vector<double> measures;
			for (size_t j = 0 ; j < GRID_ITERATOR_TESTS ; j++)
			{measures.add(grid_iterator_benchmark_stencil(g_dist,total));}
			standard_deviation(measures,mean,dev);
			time_iterator_stencil_mean.add(mean);
			time_iterator_stencil_dev.add(dev);

			//// NORMAL ////

			measures.clear();
			for (size_t j = 0 ; j < GRID_ITERATOR_TESTS ;j++)
			{measures.add(grid_iterator_benchmark_norm(g_dist,total));}
			standard_deviation(measures,mean,dev);
			time_iterator_normal_mean.add(mean);
			time_iterator_normal_dev.add(dev);

			std::cout << "Size: " << nk_grid.get(i) <<  "    " << total << std::endl;
		}
	}
}


/*! \brief Function for verlet performance report
 *
 */
template<unsigned int dim>
void grid_iterator_performance_write_report(GoogleChart & cg,
											openfpm::vector<size_t> & nk_grid,
											openfpm::vector<double> & time_iterator_stencil_mean,
											openfpm::vector<double> & time_iterator_stencil_dev,
											openfpm::vector<double> & time_iterator_normal_mean,
											openfpm::vector<double> & time_iterator_normal_dev)
{
	std::string file_mean(test_dir);
	std::string file_var(test_dir);
	file_mean += std::string("/openfpm_pdata/grid_iterator_mean_" + std::to_string(dim) + std::string("_ref"));
	file_var += std::string("/openfpm_pdata/grid_iterator_dev_" + std::to_string(dim) + std::string("_ref"));

	std::string file_mean_save = std::string("grid_iterator_mean_" + std::to_string(dim) + std::to_string("_ref"));
	std::string file_var_save = std::string("grid_iterator_dev_" + std::to_string(dim) + std::to_string("_ref"));

	openfpm::vector<size_t> xp = nk_grid;

	openfpm::vector<openfpm::vector<openfpm::vector<double>>> yp_mean;
	openfpm::vector<openfpm::vector<openfpm::vector<double>>> yp_dev;

	openfpm::vector<std::string> names;
	openfpm::vector<std::string> gnames;

/*	yp_mean.add();
	yp_dev.add();
	yp_mean.last().add(time_iterator_stencil_mean);
	yp_mean.last().add(time_iterator_normal_mean);
	yp_dev.last().add(time_iterator_stencil_dev);
	yp_dev.last().add(time_iterator_normal_dev);*/

	yp_mean.resize(1);
	yp_dev.resize(1);
	for (size_t i = 0 ; i < yp_mean.size() ; i++)
	{
		yp_mean.get(i).resize(time_iterator_stencil_mean.size());
		yp_dev.get(i).resize(time_iterator_stencil_dev.size());

		for (size_t j = 0 ; j < yp_mean.get(i).size() ; j++)
		{
			yp_mean.get(i).get(j).resize(2);
			yp_dev.get(i).get(j).resize(2);

			yp_mean.get(i).get(j).get(0) = time_iterator_stencil_mean.get(j);
			yp_mean.get(i).get(j).get(1) = time_iterator_normal_mean.get(j);

			yp_dev.get(i).get(j).get(0) = time_iterator_stencil_dev.get(j);
			yp_dev.get(i).get(j).get(1) = time_iterator_normal_dev.get(j);
		}
	}

	gnames.add("Grid iterators performance for stencil");
	names.add("Stencil specialized iterator");
	names.add("Normal iterator");

	std::string y_string = std::string("Time seconds");
	std::string x_string = std::string("Number of grid poins");


	StandardPerformanceGraph(file_mean,
			                 file_var,
							 file_mean_save,
							 file_var_save,
							 cg,
							 xp,
							 yp_mean,
							 yp_dev,
							 names,
							 gnames,
							 x_string,
							 y_string,
							 true);
}

/*! \brief Function for verlet performance report
 *
 */
template<unsigned int dim>
void grid_m2p_performance_write_report(GoogleChart & cg,
											openfpm::vector<size_t> & nk_grid,
											openfpm::vector<double> & time_m2p_mean,
											openfpm::vector<double> & time_m2p_dev)
{
	std::string file_mean(test_dir);
	std::string file_var(test_dir);
	file_mean += std::string("/openfpm_pdata/grid_m2p_mean_" + std::to_string(dim) + std::string("_ref"));
	file_var += std::string("/openfpm_pdata/grid_m2p_dev_" + std::to_string(dim) + std::string("_ref"));

	std::string file_mean_save = std::string("grid_m2p_mean_" + std::to_string(dim) + std::to_string("_ref"));
	std::string file_var_save = std::string("grid_m2p_dev_" + std::to_string(dim) + std::to_string("_ref"));

	openfpm::vector<size_t> xp = nk_grid;

	openfpm::vector<openfpm::vector<openfpm::vector<double>>> yp_mean;
	openfpm::vector<openfpm::vector<openfpm::vector<double>>> yp_dev;

	openfpm::vector<std::string> names;
	openfpm::vector<std::string> gnames;

	yp_mean.resize(1);
	yp_dev.resize(1);
	for (size_t i = 0 ; i < yp_mean.size() ; i++)
	{
		yp_mean.get(i).resize(time_iterator_stencil_mean.size());
		yp_dev.get(i).resize(time_iterator_stencil_dev.size());

		for (size_t j = 0 ; j < yp_mean.get(i).size() ; j++)
		{
			yp_mean.get(i).get(j).resize(1);
			yp_dev.get(i).get(j).resize(1);

			yp_mean.get(i).get(j).get(0) = time_m2p_mean.get(j);
			yp_dev.get(i).get(j).get(0) = time_m2p_dev.get(j);
		}
	}

	gnames.add("Grid m2p performance");
	names.add("Mesh to particle performance");

	std::string y_string = std::string("Time seconds");
	std::string x_string = std::string("Number of grid poins");


	StandardPerformanceGraph(file_mean,
			                 file_var,
							 file_mean_save,
							 file_var_save,
							 cg,
							 xp,
							 yp_mean,
							 yp_dev,
							 names,
							 gnames,
							 x_string,
							 y_string,
							 true);

}



BOOST_AUTO_TEST_CASE( grid_interpolation_benchmark_test )
{
	nk_grid_int.add(96);
	nk_grid_int.add(128);
	nk_grid_int.add(192);

	//Benchmark test for 2D and 3D
	grid_interpolation_benchmark(nk_grid_int,
			                time_inte_p2m_mean,
							time_inte_m2p_mean,
							time_inte_p2m_dev,
							time_inte_m2p_dev);
}

BOOST_AUTO_TEST_CASE( grid_iterator_benchmark_test )
{
	nk_grid_st.add(96);
	nk_grid_st.add(128);
	nk_grid_st.add(192);

	//Benchmark test for 2D and 3D
	grid_iterator_benchmark<3>(nk_grid_st,
			                time_iterator_normal_mean,
							time_iterator_stencil_mean,
							time_iterator_normal_dev,
							time_iterator_stencil_dev);
}


BOOST_AUTO_TEST_CASE(grid_iterator_performance_write_report_final)
{
	GoogleChart cg;

	//Write report for 2D and 3D
	grid_iterator_performance_write_report<3>(cg,
            								  nk_grid_st,
											  time_iterator_stencil_mean,
											  time_iterator_stencil_dev,
											  time_iterator_normal_mean,
											  time_iterator_normal_dev);


	grid_m2p_performance_write_report<3>(cg,
										nk_grid_int,
										time_inte_m2p_mean,
										time_inte_m2p_dev);

	if (create_vcluster().getProcessUnitID() == 0)
	{
		addUpdtateTime(cg);

		cg.write("grid_iterator_performance.html");
	}
}


BOOST_AUTO_TEST_SUITE_END()



#endif /* SRC_GRID_GRID_DIST_PERFORMANCE_HPP_ */
