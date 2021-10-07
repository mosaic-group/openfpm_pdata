/*
 * vector_dist_cl_performance_tests.hpp
 *
 *
 *  Created on: Mar 22, 2016
 *      Author: Yaroslav Zaluzhnyi
 */

#ifndef SRC_VECTOR_VECTOR_DIST_CL_PERFORMANCE_TESTS_HPP_
#define SRC_VECTOR_VECTOR_DIST_CL_PERFORMANCE_TESTS_HPP_

#include "Vector/vector_dist.hpp"
#include "data_type/aggregate.hpp"
#include "Plot/GoogleChart.hpp"
#include <functional>

// Property tree
struct report_cell_list_preord_tests
{
	boost::property_tree::ptree graphs;
};

report_cell_list_preord_tests report_cl_preo;

BOOST_AUTO_TEST_SUITE( celllist_part_reorder_performance_test )

///////////////////// INPUT DATA //////////////////////

// Cut-off radiuses. Can be put different number of values
openfpm::vector<float> r_cutoff {0.004, 0.007, 0.01};
// Orders of a curve. Can be put different number of values
openfpm::vector<size_t> orders = {1,2,3};
// Number of steps of moving the particles
size_t n_moving = 8;
// Moving distance (step size)
double dist = 0.03;
// The starting amount of particles (remember that this number is multiplied by number of processors you use for testing)
size_t k_start = 100000;
// The minimal amount of particles
size_t k_min = 15000;

///////////////////////////////////////////////////////

// Numbers of particles vector
openfpm::vector<size_t> n_particles;


/*! \brief Function for cell list test without an Hilbert curve reordering (unordered positioning)
 *
 */
template<unsigned int dim> void cell_list_part_reorder_random_benchmark(size_t cl_k_start,
		                                                                size_t cl_k_min,
																		openfpm::vector<float> & cl_r_cutoff,
																		openfpm::vector<size_t> & cl_n_particles)
{
	std::string str("Testing " + std::to_string(dim) + "D vector, no-order, Cell-list");
	print_test_v(str,0);

	{
		//For different r_cut
		for (size_t r = 0; r < cl_r_cutoff.size(); r++ )
		{
			Vcluster<> & v_cl = create_vcluster();

			//Cut-off radius
			float r_cut = cl_r_cutoff.get(r);

			report_cl_preo.graphs.put("performance.celllist.calc_forces_reordered" + std::to_string(dim) + "D(" + std::to_string(r) + ").rcut",r_cut);

			//Number of particles
			size_t k = cl_k_start * v_cl.getProcessingUnits();

			//Counter number for amounts of particles
			size_t k_count = 1 + log2(k/cl_k_min);

			int c = 0;

			//For different number of particles
			for (size_t k_int = k ; k_int >= cl_k_min ; k_int/=2 )
			{
				BOOST_TEST_CHECKPOINT( "Testing " << dim << "D vector without an Hilbert curve reordering k=" << k_int );

				report_cl_preo.graphs.put("performance.celllist.calc_forces_reordered" + std::to_string(dim) + "D(" + std::to_string(r) + ").npart(" + std::to_string(c) + ").n",k_int);

				if (cl_n_particles.size() < k_count)
					cl_n_particles.add(k_int);

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

				//Get a cell list

				auto NN = vd.getCellList(r_cut);
				double sum_fr_mean = 0;
				double sum_fr_dev = 0;

				benchmark_get_celllist(NN,vd,r_cut);

				//Calculate forces
				size_t l = 0;

				openfpm::vector<double> measures;
				for ( ; l < N_STAT_TEST; l++)
				{measures.add(benchmark_calc_forces<dim>(NN,vd,r_cut));}
				standard_deviation(measures,sum_fr_mean,sum_fr_dev);

				report_cl_preo.graphs.put("performance.celllist.calc_forces_reordered" + std::to_string(dim) + "D(" + std::to_string(r) + ").npart(" + std::to_string(c) + ").mean",sum_fr_mean);
				report_cl_preo.graphs.put("performance.celllist.calc_forces_reordered" + std::to_string(dim) + "D(" + std::to_string(r) + ").npart(" + std::to_string(c) + ").dev",sum_fr_dev);

				if (v_cl.getProcessUnitID() == 0)
				{std::cout << "Cut-off = " << r_cut << ", Particles = " << k_int << " time to calculate forces: " << sum_fr_mean << " dev: " << sum_fr_dev << std::endl;}

				c++;
			}
		}
	}
}


/*! \brief Function for cell list test with an Hilbert curve reordering
 *
 */
template<unsigned int dim> void cell_list_part_reorder_hilbert_benchmark(size_t cl_k_start,
		                                                                 size_t cl_k_min,
																		 size_t n_moving,
																		 double dist,
																		 openfpm::vector<float> & cl_r_cutoff,
																		 openfpm::vector<size_t> & cl_n_particles,
																		 openfpm::vector<size_t> &cl_orders)
{
	{
		// Print test
		std::string str("Testing " + std::to_string(dim) + "D vector, Hilbert curve reordering, Cell-List");
		print_test_v(str,0);

		// For different r_cut
		for (size_t r = 0; r < cl_r_cutoff.size(); r++ )
		{
			Vcluster<> & v_cl = create_vcluster();

			// Cut-off radius
			float r_cut = cl_r_cutoff.get(r);

			// Number of particles
			size_t k = cl_k_start * v_cl.getProcessingUnits();

			//For different curve orders
			for ( size_t i = 0; i < cl_orders.size(); i++)
			{
				size_t m = cl_orders.get(i);
				size_t part = 0;

				int c = 0;

				for (size_t k_int = k ; k_int >= cl_k_min ; k_int/=2, part++ )
				{
					BOOST_TEST_CHECKPOINT( "Testing " << dim << "D vector with an Hilbert curve reordering k=" << k_int );

					report_cl_preo.graphs.put("performance.celllist.calc_forces_hilb(" + std::to_string(m) + ")_reordered" + std::to_string(dim) + "D(" + std::to_string(r) + ").rcut",r_cut);
					report_cl_preo.graphs.put("performance.celllist.calc_forces_hilb(" + std::to_string(m) + ")_reordered" + std::to_string(dim) + "D(" + std::to_string(r) + ").npart(" + std::to_string(c) + ").n",k_int);

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

					vector_dist<dim,float, aggregate<float[dim]> > vd(k_int,box,bc,Ghost<dim,float>(r_cut));

					// Initialize a dist vector

					vd_initialize<dim>(vd, v_cl);
					//Reorder a vector

					double sum_reorder_mean = 0;
					double sum_reorder_dev = 0;

					openfpm::vector<double> measures;
					for (size_t h = 0 ; h < N_STAT_TEST; h++)
					{measures.add(benchmark_reorder(vd,m));}
					standard_deviation(measures,sum_reorder_mean,sum_reorder_dev);

					vd.template ghost_get<0>();

					//Get cell list

					auto NN = vd.getCellList(r_cut);
					benchmark_get_celllist(NN,vd,r_cut);

					//Calculate forces

					double sum_fr_mean = 0;
					double sum_fr_dev = 0;
					measures.clear();

					for (size_t l = 0 ; l < N_STAT_TEST ; l++)
					{measures.add(benchmark_calc_forces<dim>(NN,vd,r_cut));}
					standard_deviation(measures,sum_fr_mean,sum_fr_dev);

					report_cl_preo.graphs.put("performance.celllist.calc_forces_hilb(" + std::to_string(m) + ")_reordered" + std::to_string(dim) + "D(" + std::to_string(r) + ").npart(" + std::to_string(c) + ").mean",sum_fr_mean);
					report_cl_preo.graphs.put("performance.celllist.calc_forces_hilb(" + std::to_string(m) + ")_reordered" + std::to_string(dim) + "D(" + std::to_string(r) + ").npart(" + std::to_string(c) + ").dev",sum_fr_dev);

					if (v_cl.getProcessUnitID() == 0)
					{std::cout << "Cut-off = " << r_cut << ", Particles = " << k_int << ". Time to reorder: " << sum_reorder_mean << " dev: " << sum_reorder_dev << "      time calculate forces: " << sum_fr_mean << " dev: " << sum_fr_dev << std::endl;}

					c++;
				}
			}
		}
	}
}


BOOST_AUTO_TEST_CASE( vector_dist_cl_random_test )
{
	//Benchmark test for 2D and 3D
	cell_list_part_reorder_random_benchmark<3>(k_start,k_min,r_cutoff,n_particles);
	cell_list_part_reorder_random_benchmark<2>(k_start,k_min,r_cutoff,n_particles);
}

BOOST_AUTO_TEST_CASE( vector_dist_cl_hilbert_test )
{
	//Benchmark test for 2D and 3D
	cell_list_part_reorder_hilbert_benchmark<3>(k_start,
			                                    k_min,
												n_moving,
												dist,
												r_cutoff,
												n_particles,
												orders);

	cell_list_part_reorder_hilbert_benchmark<2>(k_start,
			                                    k_min,
												n_moving,
												dist,
												r_cutoff,
												n_particles,
												orders);
}

BOOST_AUTO_TEST_CASE(vector_dist_cl_performance_write_report)
{
	// Create a graphs

	//For different r_cut
	for (size_t r = 0; r < r_cutoff.size(); r++ )
	{
		report_cl_preo.graphs.put("graphs.graph(" + std::to_string(r) + ").type","line");
		report_cl_preo.graphs.add("graphs.graph(" + std::to_string(r) + ").title","calc_force 3D with reordered particles performance r_cut=" + std::to_string(r_cutoff.get(r)));
		report_cl_preo.graphs.add("graphs.graph(" + std::to_string(r) + ").x.title","number of particles");
		report_cl_preo.graphs.add("graphs.graph(" + std::to_string(r) + ").y.title","Time seconds");
		for (size_t i  = 0 ; i < orders.size() ; i++)
		{
			size_t m = orders.get(i);
			report_cl_preo.graphs.add("graphs.graph(" + std::to_string(r) + ").y.data(" + std::to_string(i+1) + ").source","performance.celllist.calc_forces_hilb(" + std::to_string(m) + ")_reordered3D(" + std::to_string(r) + ").npart(#).mean");
			report_cl_preo.graphs.add("graphs.graph(" + std::to_string(r) + ").x.data(" + std::to_string(i+1) + ").source","performance.celllist.calc_forces_hilb(" + std::to_string(m) + ")_reordered3D(" + std::to_string(r) + ").npart(#).n");
			report_cl_preo.graphs.add("graphs.graph(" + std::to_string(r) + ").y.data(" + std::to_string(i+1) + ").title","Hilbert(" + std::to_string(m) + ") reorder");
		}
		report_cl_preo.graphs.add("graphs.graph(" + std::to_string(r) + ").y.data(0).source","performance.celllist.calc_forces_reordered3D(" + std::to_string(r) + ").npart(#).mean");
		report_cl_preo.graphs.add("graphs.graph(" + std::to_string(r) + ").x.data(0).source","performance.celllist.calc_forces_reordered3D(" + std::to_string(r) + ").npart(#).n");
		report_cl_preo.graphs.add("graphs.graph(" + std::to_string(r) + ").y.data(0).title","Random reorder");
		report_cl_preo.graphs.add("graphs.graph(" + std::to_string(r) + ").options.log_y","true");
	}

	//For different r_cut
	for (size_t r = 0; r < r_cutoff.size(); r++ )
	{
		report_cl_preo.graphs.put("graphs.graph(" + std::to_string(r + r_cutoff.size()) + ").type","line");
		report_cl_preo.graphs.add("graphs.graph(" + std::to_string(r + r_cutoff.size()) + ").title","calc_force 2D with reordered particles performance r_cut=" + std::to_string(r_cutoff.get(r)));
		report_cl_preo.graphs.add("graphs.graph(" + std::to_string(r + r_cutoff.size()) + ").x.title","number of particles");
		report_cl_preo.graphs.add("graphs.graph(" + std::to_string(r + r_cutoff.size()) + ").y.title","Time seconds");
		for (size_t i  = 0 ; i < orders.size() ; i++)
		{
			size_t m = orders.get(i);
			report_cl_preo.graphs.add("graphs.graph(" + std::to_string(r + r_cutoff.size()) + ").y.data(" + std::to_string(i+1) + ").source","performance.celllist.calc_forces_hilb(" + std::to_string(m) + ")_reordered2D(" + std::to_string(r) + ").npart(#).mean");
			report_cl_preo.graphs.add("graphs.graph(" + std::to_string(r + r_cutoff.size()) + ").x.data(" + std::to_string(i+1) + ").source","performance.celllist.calc_forces_hilb(" + std::to_string(m) + ")_reordered2D(" + std::to_string(r) + ").npart(#).n");
			report_cl_preo.graphs.add("graphs.graph(" + std::to_string(r + r_cutoff.size()) + ").y.data(" + std::to_string(i+1) + ").title","Hilbert(" + std::to_string(m) + ") reorder");
		}
		report_cl_preo.graphs.add("graphs.graph(" + std::to_string(r + r_cutoff.size()) + ").y.data(0).source","performance.celllist.calc_forces_reordered2D(" + std::to_string(r) + ").npart(#).mean");
		report_cl_preo.graphs.add("graphs.graph(" + std::to_string(r + r_cutoff.size()) + ").x.data(0).source","performance.celllist.calc_forces_reordered2D(" + std::to_string(r) + ").npart(#).n");
		report_cl_preo.graphs.add("graphs.graph(" + std::to_string(r + r_cutoff.size()) + ").y.data(0).title","Random reorder");
		report_cl_preo.graphs.add("graphs.graph(" + std::to_string(r + r_cutoff.size()) + ").options.log_y","true");
	}

	if (create_vcluster().rank() == 0)
	{
		boost::property_tree::xml_writer_settings<std::string> settings(' ', 4);
		boost::property_tree::write_xml("celllist_partreo_performance.xml", report_cl_preo.graphs,std::locale(),settings);

		std::string file_xml_ref(test_dir);
		file_xml_ref += std::string("/openfpm_pdata/celllist_partreo_performance_ref.xml");

		GoogleChart cg;

		StandardXMLPerformanceGraph("celllist_partreo_performance.xml",file_xml_ref,cg);

		addUpdateTime(cg,create_vcluster().size(),"pdata","celllist_part_ord");

		if (create_vcluster().getProcessUnitID() == 0)
		{cg.write("celllist_part_ord.html");}
	}
}

BOOST_AUTO_TEST_SUITE_END()

#endif /* SRC_VECTOR_VECTOR_DIST_CL_PERFORMANCE_TESTS_HPP_ */
