#ifndef VECTOR_DIST_GG_MAP_PERFORMANCE_HPP_

#include "Vector/vector_dist.hpp"
#include "data_type/aggregate.hpp"
#include "Plot/GoogleChart.hpp"
#include <functional>

// Property tree
struct report_vector_gg_map_tests
{
	boost::property_tree::ptree graphs;
};

report_vector_gg_map_tests report_ggm;

///////////////////// INPUT DATA //////////////////////

// The starting amount of particles (remember that this number is multiplied by number of processors you use for testing)
size_t k_start = 100000;
// The minimal amount of particles
size_t k_min = 15000;

///////////////////////////////////////////////////////

BOOST_AUTO_TEST_SUITE( vector_gg_map_performance_test )


// Numbers of particles vector
openfpm::vector<size_t> n_particles;

template<unsigned int dim>
double benchmark_map(size_t k_int, Vcluster<HeapMemory> & v_cl)
{
	Box<dim,float> box;

	for (size_t i = 0; i < dim; i++)
	{
		box.setLow(i,0.0);
		box.setHigh(i,1.0);
	}

	// Boundary conditions
	size_t bc[dim];

	for (size_t i = 0; i < dim; i++)
	{bc[i] = PERIODIC;}

	vector_dist<dim,float, aggregate<float[dim]> > vd(v_cl.size()*v_cl.size()*k_int,box,bc,Ghost<dim,float>(0.1));

	auto & dec = vd.getDecomposition();

	int start = 0;
	int stop = k_int;

	for (size_t i = 0 ; i < v_cl.size() ; i++)
	{
		if (i == v_cl.rank())
		{continue;}

		auto & nn_box = dec.getNearSubdomains(i);

		if (nn_box.size() != 0)
		{
			// generate all particles in the near processor

			vd_initialize_box_nomap<dim>(vd,nn_box.get(0),v_cl,start,stop);

			start += k_int;
			stop += k_int;
		}
	}

	//Timer
	timer t;
	t.start();

	// benckmark map
	vd.map();

	t.stop();

	return t.getwct();
}


/*! \brief Function for cell list test without an Hilbert curve reordering (unordered positioning)
 *
 */
template<unsigned int dim>
void vector_gg_map_benchmark(size_t cl_k_start,
		                     size_t cl_k_min,
							 openfpm::vector<size_t> & cl_n_particles)
{
	std::string str("Testing " + std::to_string(dim) + "D vector, no-order, map ghost_get");
	print_test_v(str,0);

	{

		Vcluster<> & v_cl = create_vcluster();

		if (v_cl.size() != 3)
		{return;}

		//Number of particles
		size_t k = cl_k_start * v_cl.getProcessingUnits();

		//Counter number for amounts of particles
		size_t k_count = 1 + log2(k/cl_k_min);

		int c = 0;

		//For different number of particles
		for (size_t k_int = k ; k_int >= cl_k_min ; k_int/=2 )
		{
			BOOST_TEST_CHECKPOINT( "Testing " << dim << "D vector without an Hilbert curve reordering k=" << k_int );

			report_ggm.graphs.put("performance.map_" + std::to_string(dim) + "D.npart(" + std::to_string(c) + ").n",k_int);

			if (cl_n_particles.size() < k_count)
			{cl_n_particles.add(k_int);}

			openfpm::vector<double> measures;
			double sum_map_mean = 0.0;
			double sum_map_dev = 0.0;
			for (size_t h = 0 ; h < N_STAT_TEST; h++)
			{measures.add(benchmark_map<dim>(k_int,v_cl));}
			standard_deviation(measures,sum_map_mean,sum_map_dev);

			report_ggm.graphs.put("performance.map_" + std::to_string(dim) + "D.npart(" + std::to_string(c) + ").mean",sum_map_mean);
			report_ggm.graphs.put("performance.map_" + std::to_string(dim) + "D.npart(" + std::to_string(c) + ").dev",sum_map_dev);

			c++;
		}
	}
}



BOOST_AUTO_TEST_CASE( vector_dist_map_test )
{
	//Benchmark test for 2D and 3D
	vector_gg_map_benchmark<3>(k_start,k_min,n_particles);
	vector_gg_map_benchmark<2>(k_start,k_min,n_particles);
}


BOOST_AUTO_TEST_CASE(vector_dist_gg_map_performance_write_report)
{
	// Create a graphs

	report_ggm.graphs.put("graphs.graph.type","line");
	report_ggm.graphs.add("graphs.graph.title","Map performance");
	report_ggm.graphs.add("graphs.graph.x.title","number of particles");
	report_ggm.graphs.add("graphs.graph.y.title","Time seconds");
	report_ggm.graphs.add("graphs.graph.y.data(0).source","performance.map_3D.npart(#).mean");
	report_ggm.graphs.add("graphs.graph.x.data(0).source","performance.map_3D.npart(#).n");
	report_ggm.graphs.add("graphs.graph.y.data(0).title","Map function");
	report_ggm.graphs.add("graphs.graph.options.log_y","false");

	if (create_vcluster().rank() == 0)
	{
		boost::property_tree::xml_writer_settings<std::string> settings(' ', 4);
		boost::property_tree::write_xml("particles_map_performance.xml", report_ggm.graphs,std::locale(),settings);

		std::string file_xml_ref(test_dir);
		file_xml_ref += std::string("/openfpm_pdata/particles_map_performance_ref.xml");

		GoogleChart cg;

		StandardXMLPerformanceGraph("particles_map_performance.xml",file_xml_ref,cg);

		addUpdtateTime(cg,create_vcluster().size());

		if (create_vcluster().getProcessUnitID() == 0)
		{cg.write("particles_map_performance.html");}
	}
}

BOOST_AUTO_TEST_SUITE_END()

#endif

