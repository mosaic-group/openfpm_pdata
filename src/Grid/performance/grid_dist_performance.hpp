/*
 * grid_dist_performance.hpp
 *
 *  Created on: Jun 27, 2017
 *      Author: i-bird
 */

#ifndef SRC_GRID_GRID_DIST_PERFORMANCE_HPP_
#define SRC_GRID_GRID_DIST_PERFORMANCE_HPP_

#define GRID_ITERATOR_TESTS  30
#define GRID_INTERPOLATION_TESTS  30

// Vectors to store the data for 3D
openfpm::vector<size_t> nk_grid_st;
openfpm::vector<size_t> nk_grid_int;

// Property tree
struct report_grid_iterator_test
{
	boost::property_tree::ptree graphs;
};

report_grid_iterator_test report_grid_iterator;

BOOST_AUTO_TEST_SUITE( grid_iterator_performance_test )

constexpr int x = 0;
constexpr int y = 1;
constexpr int z = 2;

/*! \brief
 *
 *
 */
void grid_interpolation_benchmark(openfpm::vector<size_t> & nk_grid,
									boost::property_tree::ptree & interpolation_graph)
{
	for (size_t k = 0 ; k < nk_grid.size() ; k++)
	{
		size_t np = nk_grid.get(k);
		size_t sz[] = {np,np,np};

		Ghost<3,long int> gg(3);

		Box<3,float> domain({0.0,0.0,0.0},{1.0,1.0,1.0});

		size_t tot_part = np*np*np;

		std::string base_p2m = "performance.interpolation.p2m(" + std::to_string(k) + ")";
		std::string base_m2p = "performance.interpolation.m2p(" + std::to_string(k) + ")";

		interpolation_graph.put(base_p2m + ".grid.dim",3);
		interpolation_graph.put(base_m2p + ".grid.dim",3);

		interpolation_graph.put(base_p2m + ".grid.x",np);
		interpolation_graph.put(base_p2m + ".grid.y",np);
		interpolation_graph.put(base_p2m + ".grid.z",np);

		interpolation_graph.put(base_m2p + ".grid.x",np);
		interpolation_graph.put(base_m2p + ".grid.y",np);
		interpolation_graph.put(base_m2p + ".grid.z",np);

		interpolation_graph.put(base_p2m + ".particles",tot_part);
		interpolation_graph.put(base_m2p + ".particles",tot_part);

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

		interpolation_graph.put(base_p2m + ".data.mean",mean);
		interpolation_graph.put(base_p2m + ".data.dev",dev);

		std::cout << "Time particles to mesh " << mean << std::endl;

		measures.clear();
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

		interpolation_graph.put(base_m2p + ".data.mean",mean);
		interpolation_graph.put(base_m2p + ".data.dev",dev);

		std::cout << "Time mesh to particles " << mean << std::endl;
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
														   boost::property_tree::ptree & iterator_graph)
{
	std::string str("Testing " + std::to_string(dim) + "D grid iterator stencil and normal");
	print_test_v(str,0);

	{
		//For different grid sizes
		for (size_t i = 0; i < nk_grid.size(); i++ )
		{
			size_t sz[dim];

			//Number of particles
			size_t k = nk_grid.get(i);

			std::string base = "performance.grid.iterators(" + std::to_string(i) + ")";
			iterator_graph.put(base + ".grid.dim",3);

			iterator_graph.put(base + ".grid.x",k);
			iterator_graph.put(base + ".grid.y",k);
			iterator_graph.put(base + ".grid.z",k);

			BOOST_TEST_CHECKPOINT( "Testing " << dim << "D grid iterator performance k=" << k );

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

			iterator_graph.put(base + ".stencil.data.mean",mean);
			iterator_graph.put(base + ".stencil.data.dev",dev);

			std::cout << "Size: " << nk_grid.get(i) <<  "   stencil: " << mean << std::endl;

			//// NORMAL ////

			measures.clear();
			for (size_t j = 0 ; j < GRID_ITERATOR_TESTS ;j++)
			{measures.add(grid_iterator_benchmark_norm(g_dist,total));}
			standard_deviation(measures,mean,dev);

			iterator_graph.put(base + ".normal.data.mean",mean);
			iterator_graph.put(base + ".normal.data.dev",dev);

			std::cout << "Size: " << nk_grid.get(i) <<  "   normal: " << mean << std::endl;
		}
	}
}



BOOST_AUTO_TEST_CASE( grid_interpolation_benchmark_test )
{
	nk_grid_int.add(96);
	nk_grid_int.add(128);
	nk_grid_int.add(192);

	//Benchmark test for 2D and 3D
	grid_interpolation_benchmark(nk_grid_int,
							report_grid_iterator.graphs);
}

BOOST_AUTO_TEST_CASE( grid_iterator_benchmark_test )
{
	nk_grid_st.add(96);
	nk_grid_st.add(128);
	nk_grid_st.add(192);

	//Benchmark test for 2D and 3D
	grid_iterator_benchmark<3>(nk_grid_st,
								report_grid_iterator.graphs);
}


BOOST_AUTO_TEST_CASE(grid_iterator_performance_write_report_final)
{
	report_grid_iterator.graphs.put("graphs.graph(0).type","line");
	report_grid_iterator.graphs.add("graphs.graph(0).title","Grid iterators performance for stencil");
	report_grid_iterator.graphs.add("graphs.graph(0).x.title","Number of grid points");
	report_grid_iterator.graphs.add("graphs.graph(0).y.title","Time seconds");
	report_grid_iterator.graphs.add("graphs.graph(0).y.data(0).source","performance.grid.iterators(#).normal.data.mean");
	report_grid_iterator.graphs.add("graphs.graph(0).x.data(0).source","performance.grid.iterators(#).grid.x");
	report_grid_iterator.graphs.add("graphs.graph(0).y.data(0).title","Normal iterator");
	report_grid_iterator.graphs.add("graphs.graph(0).y.data(1).source","performance.grid.iterators(#).stencil.data.mean");
	report_grid_iterator.graphs.add("graphs.graph(0).x.data(1).source","performance.grid.iterators(#).grid.x");
	report_grid_iterator.graphs.add("graphs.graph(0).y.data(1).title","Stencil specialized iterator");
	report_grid_iterator.graphs.add("graphs.graph(0).options.log_y","true");

	report_grid_iterator.graphs.put("graphs.graph(1).type","line");
	report_grid_iterator.graphs.add("graphs.graph(1).title","Grid p2m performance");
	report_grid_iterator.graphs.add("graphs.graph(1).x.title","Number of grid points");
	report_grid_iterator.graphs.add("graphs.graph(1).y.title","Time seconds");
	report_grid_iterator.graphs.add("graphs.graph(1).y.data(0).source","performance.interpolation.p2m(#).data.mean");
	report_grid_iterator.graphs.add("graphs.graph(1).y.data(0).title","Interpolation p2m");
	report_grid_iterator.graphs.add("graphs.graph(1).x.data(0).source","performance.interpolation.p2m(#).grid.x");
	report_grid_iterator.graphs.add("graphs.graph(1).options.log_y","true");

	report_grid_iterator.graphs.put("graphs.graph(2).type","line");
	report_grid_iterator.graphs.add("graphs.graph(2).title","Grid m2p performance");
	report_grid_iterator.graphs.add("graphs.graph(2).x.title","Number of grid points");
	report_grid_iterator.graphs.add("graphs.graph(2).y.title","Time seconds");
	report_grid_iterator.graphs.add("graphs.graph(2).x","performance.interpolation.m2p(#).grid.x");
	report_grid_iterator.graphs.add("graphs.graph(2).y.data(0).source","performance.interpolation.m2p(#).data.mean");
	report_grid_iterator.graphs.add("graphs.graph(2).y.data(0).title","Interpolation m2p");
	report_grid_iterator.graphs.add("graphs.graph(2).x.data(0).source","performance.interpolation.m2p(#).grid.x");
	report_grid_iterator.graphs.add("graphs.graph(2).options.log_y","true");

	if (create_vcluster().rank() == 0)
	{
		boost::property_tree::xml_writer_settings<std::string> settings(' ', 4);
		boost::property_tree::write_xml("grid_performance.xml", report_grid_iterator.graphs,std::locale(),settings);

		GoogleChart cg;

		std::string file_xml_ref(test_dir);
		file_xml_ref += std::string("/openfpm_pdata/grid_performance_ref.xml");

		StandardXMLPerformanceGraph("grid_performance.xml",file_xml_ref,cg);

		if (create_vcluster().getProcessUnitID() == 0)
		{
			addUpdtateTime(cg,create_vcluster().size(),"data","grid_performance");

			cg.write("grid_performance.html");
		}
	}
}


BOOST_AUTO_TEST_SUITE_END()



#endif /* SRC_GRID_GRID_DIST_PERFORMANCE_HPP_ */
