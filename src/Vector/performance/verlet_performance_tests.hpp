/*
 * vector_dist_verlet_performance_tests.hpp
 *
 *  Created on: Mar 9, 2016
 *      Author: Yaroslav Zaluzhnyi
 */

#ifndef SRC_VECTOR_VECTOR_DIST_VERLET_PERFORMANCE_TESTS_HPP_
#define SRC_VECTOR_VECTOR_DIST_VERLET_PERFORMANCE_TESTS_HPP_

#include "util/stat/common_statistics.hpp"

// Property tree
struct report_verlet_tests
{
	boost::property_tree::ptree graphs;
};

report_verlet_tests report_vl;

/*! \brief Print a string about the test
 *
 * \param test string to print
 * \param sz size
 *
 */
void print_test_v(std::string test, size_t sz)
{
	if (create_vcluster().getProcessUnitID() == 0)
		std::cout << test << " " << sz << "\n";
}

BOOST_AUTO_TEST_SUITE( verletlist_performance_test )

///////////////////// INPUT DATA //////////////////////

// Cut-off radiuses. Can be put different number of values
openfpm::vector<float> r_cutoff {0.004, 0.007, 0.01};
// Orders of a curve. Can be put different number of values
// The starting amount of particles (remember that this number is multiplied by number of processors you use for testing)
size_t k_start = 100000;
// The minimal amount of particles
size_t k_min = 15000;
// Ghost part of distributed vector

///////////////////////////////////////////////////////

// Numbers of particles vector
openfpm::vector<size_t> n_particles;

/*! \brief Function for verlet test without an Hilbert curve reordering (unordered positioning)
 *
 */
template<unsigned int dim> void vd_verlet_random_benchmark(size_t k_start,
		                                                   size_t k_min,
														   openfpm::vector<float> & r_cutoff,
														   openfpm::vector<size_t> & n_particles)
{
	std::string str("Testing " + std::to_string(dim) + "D vector no-order, Verlet-list");
	print_test_v(str,0);

	{
		//For different r_cut
		for (size_t r = 0; r < r_cutoff.size(); r++ )
		{
			Vcluster<> & v_cl = create_vcluster();

			//Cut-off radius
			float r_cut = r_cutoff.get(r);

			report_vl.graphs.put("performance.verletlist.getVerletList" + std::to_string(dim) + "D(" + std::to_string(r) + ").rcut",r_cut);

			//Number of particles
			size_t k = k_start * v_cl.getProcessingUnits();

			//Counter number for amounts of particles
			size_t k_count = 1 + log2(k/k_min);

			int c = 0;

			for (size_t k_int = k ; k_int >= k_min ; k_int/=2 )
			{
				BOOST_TEST_CHECKPOINT( "Testing " << dim << "D vector without an Hilbert curve reordering k=" << k_int );

				report_vl.graphs.put("performance.verletlist.getVerletList" + std::to_string(dim) + "D(" + std::to_string(r) + ").npart(" + std::to_string(c) + ").n",k_int);
				report_vl.graphs.put("performance.verletlist.calc_forces" + std::to_string(dim) + "D(" + std::to_string(r) + ").npart(" + std::to_string(c) + ").n",k_int);

				if (n_particles.size() < k_count)
					n_particles.add(k_int);

				Box<dim,float> box;

				for (size_t i = 0; i < dim; i++)
				{
					box.setLow(i,0.0);
					box.setHigh(i,1.0);
				}

				// Boundary conditions
				size_t bc[dim];
				for (size_t i = 0; i < dim; i++)
					bc[i] = PERIODIC;

				vector_dist<dim,float, aggregate<float[dim]> > vd(k_int,box,bc,Ghost<dim,float>(r_cut));

				// Initialize a dist vector
				vd_initialize<dim>(vd, v_cl);

				vd.template ghost_get<0>();

				//Get verlet list

				openfpm::vector<double> measures;
				double sum_verlet_mean = 0;
				double sum_verlet_dev = 0;
				for (size_t n = 0 ; n < N_STAT_TEST; n++)
				{measures.add(benchmark_get_verlet(vd,r_cut));}
				standard_deviation(measures,sum_verlet_mean,sum_verlet_dev);

				report_vl.graphs.put("performance.verletlist.getVerletList" + std::to_string(dim) + "D(" + std::to_string(r) + ").npart(" + std::to_string(c) + ").mean",sum_verlet_mean);
				report_vl.graphs.put("performance.verletlist.getVerletList" + std::to_string(dim) + "D(" + std::to_string(r) + ").npart(" + std::to_string(c) + ").dev",sum_verlet_dev);

				//Calculate forces

				auto NN = vd.getCellList(r_cut);
				double sum_fr_mean = 0;
				double sum_fr_dev = 0;

				measures.clear();
				for (size_t l = 0 ; l < N_STAT_TEST ; l++)
				{measures.add(benchmark_calc_forces<dim>(NN,vd,r_cut));}
				standard_deviation(measures,sum_fr_mean,sum_fr_dev);

				report_vl.graphs.put("performance.verletlist.calc_forces" + std::to_string(dim) + "D(" + std::to_string(r) + ").npart(" + std::to_string(c) + ").mean",sum_fr_mean);
				report_vl.graphs.put("performance.verletlist.calc_forces" + std::to_string(dim) + "D(" + std::to_string(r) + ").npart(" + std::to_string(c) + ").dev",sum_fr_dev);

				if (v_cl.getProcessUnitID() == 0)
				{std::cout << "Particles: " << k_int << "," << "cut-off: " << r_cut << " time to construct a Verlet list = " << sum_verlet_mean << " dev: " << sum_verlet_dev << "    calculate force = " << sum_fr_mean << " dev: " << sum_fr_dev << std::endl;}

				c++;
			}
		}
	}
}


BOOST_AUTO_TEST_CASE( vector_dist_verlet_test )
{
	//Benchmark test for 2D and 3D
	vd_verlet_random_benchmark<3>(k_start,k_min,r_cutoff,n_particles);
	vd_verlet_random_benchmark<2>(k_start,k_min,r_cutoff,n_particles);
}

BOOST_AUTO_TEST_CASE(vector_dist_verlet_performance_write_report)
{
	//For different r_cut
	for (size_t r = 0; r < r_cutoff.size(); r++ )
	{
		report_vl.graphs.put("graphs.graph(" + std::to_string(r) + ").type","line");
		report_vl.graphs.add("graphs.graph(" + std::to_string(r) + ").title","getVerletList 3D performance r_cut=" + std::to_string(r_cutoff.get(r)));
		report_vl.graphs.add("graphs.graph(" + std::to_string(r) + ").x.title","number of particles");
		report_vl.graphs.add("graphs.graph(" + std::to_string(r) + ").y.title","Time seconds");
		report_vl.graphs.add("graphs.graph(" + std::to_string(r) + ").y.data(0).source","performance.verletlist.getVerletList3D(" + std::to_string(r) + ").npart(#).mean");
		report_vl.graphs.add("graphs.graph(" + std::to_string(r) + ").x.data(0).source","performance.verletlist.getVerletList3D(" + std::to_string(r) + ").npart(#).n");
		report_vl.graphs.add("graphs.graph(" + std::to_string(r) + ").y.data(0).title","Verlet-list");
		report_vl.graphs.add("graphs.graph(" + std::to_string(r) + ").options.log_y","true");
	}

	//For different r_cut
	for (size_t r = 0; r < r_cutoff.size(); r++ )
	{
		report_vl.graphs.put("graphs.graph(" + std::to_string(r_cutoff.size() + r) + ").type","line");
		report_vl.graphs.add("graphs.graph(" + std::to_string(r_cutoff.size() + r) + ").title","calc_force 3D performance r_cut=" + std::to_string(r_cutoff.get(r)));
		report_vl.graphs.add("graphs.graph(" + std::to_string(r_cutoff.size() + r) + ").x.title","number of particles");
		report_vl.graphs.add("graphs.graph(" + std::to_string(r_cutoff.size() + r) + ").y.title","Time seconds");
		report_vl.graphs.add("graphs.graph(" + std::to_string(r_cutoff.size() + r) + ").y.data(0).source","performance.verletlist.calc_forces3D(" + std::to_string(r) + ").npart(#).mean");
		report_vl.graphs.add("graphs.graph(" + std::to_string(r_cutoff.size() + r) + ").x.data(0).source","performance.verletlist.calc_forces3D(" + std::to_string(r) + ").npart(#).n");
		report_vl.graphs.add("graphs.graph(" + std::to_string(r_cutoff.size() + r) + ").y.data(0).title","Verlet-list");
		report_vl.graphs.add("graphs.graph(" + std::to_string(r_cutoff.size() + r) + ").options.log_y","true");
	}

	//For different r_cut
	for (size_t r = 0; r < r_cutoff.size(); r++ )
	{
		report_vl.graphs.put("graphs.graph(" + std::to_string(2*r_cutoff.size() + r) + ").type","line");
		report_vl.graphs.add("graphs.graph(" + std::to_string(2*r_cutoff.size() + r) + ").title","getVerletList 2D performance r_cut=" + std::to_string(r_cutoff.get(r)));
		report_vl.graphs.add("graphs.graph(" + std::to_string(2*r_cutoff.size() + r) + ").x.title","number of particles");
		report_vl.graphs.add("graphs.graph(" + std::to_string(2*r_cutoff.size() + r) + ").y.title","Time seconds");
		report_vl.graphs.add("graphs.graph(" + std::to_string(2*r_cutoff.size() + r) + ").y.data(0).source","performance.verletlist.getVerletList2D(" + std::to_string(r) + ").npart(#).mean");
		report_vl.graphs.add("graphs.graph(" + std::to_string(2*r_cutoff.size() + r) + ").x.data(0).source","performance.verletlist.getVerletList2D(" + std::to_string(r) + ").npart(#).n");
		report_vl.graphs.add("graphs.graph(" + std::to_string(2*r_cutoff.size() + r) + ").y.data(0).title","Verlet-list");
		report_vl.graphs.add("graphs.graph(" + std::to_string(2*r_cutoff.size() + r) + ").options.log_y","true");
	}

	//For different r_cut
	for (size_t r = 0; r < r_cutoff.size(); r++ )
	{
		report_vl.graphs.put("graphs.graph(" + std::to_string(3*r_cutoff.size() + r) + ").type","line");
		report_vl.graphs.add("graphs.graph(" + std::to_string(3*r_cutoff.size() + r) + ").title","calc_force 2D performance r_cut=" + std::to_string(r_cutoff.get(r)));
		report_vl.graphs.add("graphs.graph(" + std::to_string(3*r_cutoff.size() + r) + ").x.title","number of particles");
		report_vl.graphs.add("graphs.graph(" + std::to_string(3*r_cutoff.size() + r) + ").y.title","Time seconds");
		report_vl.graphs.add("graphs.graph(" + std::to_string(3*r_cutoff.size() + r) + ").y.data(0).source","performance.verletlist.calc_forces2D(" + std::to_string(r) + ").npart(#).mean");
		report_vl.graphs.add("graphs.graph(" + std::to_string(3*r_cutoff.size() + r) + ").x.data(0).source","performance.verletlist.calc_forces2D(" + std::to_string(r) + ").npart(#).n");
		report_vl.graphs.add("graphs.graph(" + std::to_string(3*r_cutoff.size() + r) + ").y.data(0).title","Verlet-list");
		report_vl.graphs.add("graphs.graph(" + std::to_string(3*r_cutoff.size() + r) + ").options.log_y","true");
	}

	if (create_vcluster().rank() == 0)
	{
		boost::property_tree::xml_writer_settings<std::string> settings(' ', 4);
		boost::property_tree::write_xml("verletlist_performance.xml", report_vl.graphs,std::locale(),settings);

		std::string file_xml_ref(test_dir);
		file_xml_ref += std::string("/openfpm_pdata/verletlist_performance_ref.xml");

		GoogleChart cg;

		StandardXMLPerformanceGraph("verletlist_performance.xml",file_xml_ref,cg);

		addUpdateTime(cg,create_vcluster().size(),"pdata","verletlist_performance");

		if (create_vcluster().getProcessUnitID() == 0)
		{cg.write("verletlist_performance.html");}
	}
}

BOOST_AUTO_TEST_SUITE_END()

#endif /* SRC_VECTOR_VECTOR_DIST_VERLET_PERFORMANCE_TESTS_HPP_ */
