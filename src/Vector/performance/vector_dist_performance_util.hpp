/*
 * vector_dist_performance_util.hpp
 *
 *  Created on: Jun 17, 2016
 *      Author: Yaroslav Zaluzhnyi
 */

#ifndef SRC_VECTOR_VECTOR_DIST_PERFORMANCE_UTIL_HPP_
#define SRC_VECTOR_VECTOR_DIST_PERFORMANCE_UTIL_HPP_

#include "cl_comp_performance_graph.hpp"
#include "Plot/GoogleChart.hpp"
#include "cl_part_performance_graph.hpp"
#include "vector_dist_performance_common.hpp"

// Number of tests
#define N_VERLET_TEST 3
#define N_STAT_TEST 30

/*! \brief Standard deviation
 *
 * \param measures set of measures
 * \param mean the mean of the measures
 *
 * \return the standard deviation
 *
 */
static inline void standard_deviation(openfpm::vector<double> measures, double & mean, double & dev)
{
	for (size_t i = 0 ; i < measures.size() ; i++)
		mean += measures.get(i);
	mean /= measures.size();

	dev = 0;
	for (size_t i = 0 ; i < measures.size() ; i++)
		dev += (measures.get(i) - mean)*(measures.get(i) - mean);

	dev = sqrt(dev / (measures.size() - 1));
}

/*! \brief Print out only ones (no matter how many processors involved)
 *
 * \param test, sz Data to print out
 */
void print_test_v(std::string test)
{
	if (create_vcluster().getProcessUnitID() == 0)
		std::cout << test  << "\n";
}


/*! \brief Benchmark particles' forces time
 *
 * \param NN Cell list
 * \param vd Distributed vector
 * \param r_cut Cut-off radius
 *
 * \return real time
 */
template<unsigned int dim, size_t prp = 0, typename T, typename V> double benchmark_calc_forces(T & NN, V & vd, float r_cut)
{
	//Timer
	timer t;
	t.start();

	calc_forces<dim,prp>(NN,vd,r_cut);

	t.stop();

	return t.getwct();
}

/*! \brief Benchmark reordering time
 *
 * \param vd Distributed vector
 * \param m Order of an Hilbert curve
 *
 * \return real time
 */
template<typename V> double benchmark_reorder(V & vd, size_t m)
{
	//Timer
	timer t;
	t.start();

	//Reorder
	vd.reorder(m);

	t.stop();

	return t.getwct();
}

/*! \brief Benchmark celllist getting time
 *
 * \param NN Cell list
 * \param vd Distributed vector
 * \param r_cut Cut-off radius
 *
 * \return real time
 */
template<typename T, typename V> double benchmark_get_celllist(T & NN, V & vd, float r_cut)
{
	// Cell list timer
	timer t;
	t.start();

	//get cell list
	NN = vd.getCellList(r_cut);

	t.stop();

	return t.getwct();
}

/*! \brief Benchmark verlet getting time
 *
 * \param vd Distributed vector
 * \param r_cut Cut-off radius
 *
 * \return real time
 */
template<typename V> double benchmark_get_verlet(V & vd, float r_cut)
{
	openfpm::vector<openfpm::vector<size_t>> verlet;
	//Timer
	timer t;
	t.start();

	//get verlet
	auto vr = vd.getVerlet(r_cut);

	t.stop();

	return t.getwct();
}



/*! \brief Benchmark particles' forces time
 *
 * \param NN Cell list hilbert
 * \param vd Distributed vector
 * \param r_cut Cut-off radius
 *
 * \return real time
 */
template<unsigned int dim, size_t prp = 0, typename T, typename V> double benchmark_calc_forces_hilb(T & NN, V & vd, float r_cut)
{
	//Forces timer
	timer t;
	t.start();

	calc_forces_hilb<dim,prp>(NN,vd,r_cut);

	t.stop();

	return t.getwct();
}

/*! \brief Benchmark celllist hilbert getting time
 *
 * \param NN Cell list hilbert
 * \param vd Distributed vector
 * \param r_cut Cut-off radius
 *
 * \return real time
 */
template<typename T, typename V> double benchmark_get_celllist_hilb(T & NN, V & vd, float r_cut)
{
	// Cell list timer
	timer t;
	t.start();

	//get cell list
	NN = vd.getCellList_hilb(r_cut);

	t.stop();

	return t.getwct();
}

/*! \brief Move particles in random direction but the same distance
 *
 * \param vd Distributed vector
 * \param dist Distance for which the particles are moved (note that the direction is random)
 *
 */
template<unsigned int dim, typename v_dist> void move_particles(v_dist & vd, double dist)
{
	double phi;
	double theta;
	const double pi = 3.14159;

	auto it = vd.getIterator();

	while (it.isNext())
	{
		phi = rand()/double(RAND_MAX);
		theta = rand()/double(RAND_MAX);
		phi *= 2.0*pi;
		theta *= pi;

		auto key = it.get();

		if(dim == 1)
			vd.getPos(key)[0] += dist*sin(theta)*cos(phi);
		else if(dim == 2)
		{
			vd.getPos(key)[0] += dist*sin(theta)*cos(phi);
			vd.getPos(key)[1] += dist*sin(theta)*sin(phi);
		}
		else if(dim == 3)
		{
			vd.getPos(key)[0] += dist*sin(theta)*cos(phi);
			vd.getPos(key)[1] += dist*sin(theta)*sin(phi);
			vd.getPos(key)[2] += dist*cos(phi);
		}

		++it;
	}
}

////////////////////////// BENCHMARK TEST FUNCTIONS ///////////////////////////

/*! \brief Function for verlet test without an Hilbert curve reordering (unordered positioning)
 *
 */
template<unsigned int dim> void vd_verlet_random_benchmark(size_t k_start,
		                                                   size_t k_min,
														   openfpm::vector<float> & r_cutoff,
														   openfpm::vector<size_t> & n_particles,
														   openfpm::vector<openfpm::vector<double>> & time_force_mean,
														   openfpm::vector<openfpm::vector<double>> & time_create_mean,
														   openfpm::vector<openfpm::vector<double>> & time_force_dev,
														   openfpm::vector<openfpm::vector<double>> & time_create_dev)
{
	time_force_mean.resize(r_cutoff.size());
	time_create_mean.resize(r_cutoff.size());
	time_force_dev.resize(r_cutoff.size());
	time_create_dev.resize(r_cutoff.size());

	std::string str("Testing " + std::to_string(dim) + "D vector no-order, Verlet-list");
	print_test_v(str);

	{
		//For different r_cut
		for (size_t r = 0; r < r_cutoff.size(); r++ )
		{
			Vcluster & v_cl = create_vcluster();

			//Cut-off radius
			float r_cut = r_cutoff.get(r);

			//Number of particles
			size_t k = k_start * v_cl.getProcessingUnits();

			//Counter number for amounts of particles
			size_t k_count = 1 + log2(k/k_min);

			for (size_t k_int = k ; k_int >= k_min ; k_int/=2 )
			{
				BOOST_TEST_CHECKPOINT( "Testing " << dim << "D vector without an Hilbert curve reordering k=" << k_int );

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

				vector_dist<dim,float, aggregate<float[dim]>, CartDecomposition<dim,float> > vd(k_int,box,bc,Ghost<dim,float>(r_cut));

				// Initialize a dist vector
				vd_initialize<dim>(vd, v_cl, k_int);

				vd.template ghost_get<0>();

				//Get verlet list

				openfpm::vector<double> measures;
				double sum_verlet_mean = 0;
				double sum_verlet_dev = 0;
				for (size_t n = 0 ; n < N_STAT_TEST; n++)
					measures.add(benchmark_get_verlet(vd,r_cut));
				standard_deviation(measures,sum_verlet_mean,sum_verlet_dev);

				//Average total time
				time_create_mean.get(r).add(sum_verlet_mean);
				time_create_dev.get(r).add(sum_verlet_dev);

				//Calculate forces

				auto NN = vd.getCellList(r_cut);
				double sum_fr_mean = 0;
				double sum_fr_dev = 0;

				measures.clear();
				for (size_t l = 0 ; l < N_STAT_TEST ; l++)
					measures.add(benchmark_calc_forces<dim>(NN,vd,r_cut));
				standard_deviation(measures,sum_fr_mean,sum_fr_dev);
				time_force_mean.get(r).add(sum_fr_mean);
				time_force_dev.get(r).add(sum_fr_dev);

				if (v_cl.getProcessUnitID() == 0)
					std::cout << "Particles: " << k_int << "," << "cut-off: " << r_cut << " time to construct a Verlet list = " << sum_verlet_mean << " dev: " << sum_verlet_dev << "    calculate force = " << sum_fr_mean << " dev: " << sum_fr_dev << std::endl;
			}
		}
	}
}

/*! \brief Function for verlet test with an Hilbert curve reordering
 *
 */
template<unsigned int dim> void vd_verlet_hilbert_benchmark(size_t k_start, size_t k_min, double ghost_part,openfpm::vector<float> & r_cutoff, openfpm::vector<size_t> & n_particles, openfpm::vector<size_t> &orders, openfpm::vector<openfpm::vector<openfpm::vector<double>>> &time_hilb, openfpm::vector<openfpm::vector<openfpm::vector<double>>> &time_total_hilb)
{
	time_hilb.resize(r_cutoff.size());
	for (size_t r = 0; r < time_hilb.size(); r++)
	{
		time_hilb.get(r).resize(n_particles.size());
		for (size_t k = 0; k < time_hilb.get(r).size(); k++)
		{
			time_hilb.get(r).get(k).resize(orders.size());
		}
	}

	time_total_hilb.resize(r_cutoff.size());
	for (size_t r = 0; r < time_total_hilb.size(); r++)
	{
		time_total_hilb.get(r).resize(n_particles.size());
		for (size_t k = 0; k < time_total_hilb.get(r).size(); k++)
		{
			time_total_hilb.get(r).get(k).resize(orders.size());
		}
	}

	std::string str("Testing " + std::to_string(dim) + "D vector, Hilbert curve reordering, Verlet-list");
	print_test_v(str);

	// For different r_cut
	for (size_t r = 0; r < r_cutoff.size(); r++ )
	{
		Vcluster & v_cl = create_vcluster();

		//Cut-off radius
		float r_cut = r_cutoff.get(r);

		// Number of particles
		size_t k = k_start * v_cl.getProcessingUnits();

		//For different curve orders
		for ( size_t i = 0; i < orders.size(); i++)
		{
			size_t m = orders.get(i);
			size_t part = 0;

			for (size_t k_int = k ; k_int >= k_min ; k_int/=2, part++ )
			{
				BOOST_TEST_CHECKPOINT( "Testing " << dim << "D vector with an Hilbert curve reordering k=" << k_int );

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

				vector_dist<dim,float, aggregate<float[dim]>, CartDecomposition<dim,float> > vd(k_int,box,bc,Ghost<dim,float>(ghost_part));

				// Initialize a dist vector
				vd_initialize<dim>(vd, v_cl, k_int);

				vd.template ghost_get<0>();

				//Reorder a vector

				double sum_reorder = 0;
				for (size_t h = 0 ; h < N_VERLET_TEST; h++)
					sum_reorder += benchmark_reorder(vd,m);
				sum_reorder /= N_VERLET_TEST;

				//Get verlet list

				double sum_verlet = 0;

				for (size_t n = 0 ; n < N_VERLET_TEST; n++)
					sum_verlet += benchmark_get_verlet(vd,r_cut);
				sum_verlet /= N_VERLET_TEST;
				//Average total time
				time_total_hilb.get(r).get(part).get(i) = sum_verlet;

				//Calculate forces

				auto NN = vd.getCellList(r_cut);
				double sum_forces = 0;

				for (size_t l = 0 ; l < N_VERLET_TEST; l++)
					sum_forces += benchmark_calc_forces<dim>(NN,vd,r_cut);
				sum_forces /= N_VERLET_TEST;
				time_hilb.get(r).get(part).get(i) = sum_forces;

				if (v_cl.getProcessUnitID() == 0)
					std::cout << "Order = " << m << ", Cut-off = " << r_cut << ", Particles = " << k_int << ". Time to reorder: " << sum_reorder << " time to get the verlet-list: " << sum_verlet << " time to calculate forces: " << sum_forces << std::endl;
			}
		}
	}
}

/*! \brief Function for cell list test without an Hilbert curve reordering (unordered positioning)
 *
 */
template<unsigned int dim> void cell_list_part_reorder_random_benchmark(size_t cl_k_start,
		                                                                size_t cl_k_min,
																		openfpm::vector<float> & cl_r_cutoff,
																		openfpm::vector<size_t> & cl_n_particles,
																		openfpm::vector<openfpm::vector<double>> & cl_time_rand_mean,
																		openfpm::vector<openfpm::vector<double>> & cl_time_rand_dev)
{
	cl_time_rand_mean.resize(cl_r_cutoff.size());
	cl_time_rand_dev.resize(cl_r_cutoff.size());

	std::string str("Testing " + std::to_string(dim) + "D vector, no-order, Cell-list");
	print_test_v(str);

	{
		//For different r_cut
		for (size_t r = 0; r < cl_r_cutoff.size(); r++ )
		{
			Vcluster & v_cl = create_vcluster();

			//Cut-off radius
			float r_cut = cl_r_cutoff.get(r);

			//Number of particles
			size_t k = cl_k_start * v_cl.getProcessingUnits();

			//Counter number for amounts of particles
			size_t k_count = 1 + log2(k/cl_k_min);

			//For different number of particles
			for (size_t k_int = k ; k_int >= cl_k_min ; k_int/=2 )
			{
				BOOST_TEST_CHECKPOINT( "Testing " << dim << "D vector without an Hilbert curve reordering k=" << k_int );

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

				vector_dist<dim,float, aggregate<float[dim]>, CartDecomposition<dim,float> > vd(k_int,box,bc,Ghost<dim,float>(r_cut));

				// Initialize a dist vector
				vd_initialize<dim>(vd, v_cl, k_int);

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
					measures.add(benchmark_calc_forces<dim>(NN,vd,r_cut));
				standard_deviation(measures,sum_fr_mean,sum_fr_dev);

				cl_time_rand_mean.get(r).add(sum_fr_mean);
				cl_time_rand_dev.get(r).add(sum_fr_dev);

				if (v_cl.getProcessUnitID() == 0)
					std::cout << "Cut-off = " << r_cut << ", Particles = " << k_int << " time to calculate forces: " << sum_fr_mean << " dev: " << sum_fr_dev << std::endl;
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
																		 openfpm::vector<size_t> &cl_orders,
																		 openfpm::vector<openfpm::vector<openfpm::vector<double>>> &cl_time_hilb_mean,
																		 openfpm::vector<openfpm::vector<openfpm::vector<double>>> &cl_time_reorder_hilb_mean,
																		 openfpm::vector<openfpm::vector<openfpm::vector<double>>> &cl_time_hilb_dev,
																		 openfpm::vector<openfpm::vector<openfpm::vector<double>>> &cl_time_reorder_hilb_dev)
{
	{
		cl_time_hilb_mean.resize(cl_r_cutoff.size());
		cl_time_hilb_dev.resize(cl_r_cutoff.size());
		for (size_t r = 0; r < cl_time_hilb_mean.size(); r++)
		{
			cl_time_hilb_mean.get(r).resize(cl_n_particles.size());
			cl_time_hilb_dev.get(r).resize(cl_n_particles.size());
			for (size_t k = 0; k < cl_time_hilb_mean.get(r).size(); k++)
			{
				cl_time_hilb_mean.get(r).get(k).resize(cl_orders.size());
				cl_time_hilb_dev.get(r).get(k).resize(cl_orders.size());
			}
		}


		cl_time_reorder_hilb_mean.resize(cl_r_cutoff.size());
		cl_time_reorder_hilb_dev.resize(cl_r_cutoff.size());
		for (size_t r = 0; r < cl_time_reorder_hilb_mean.size(); r++)
		{
			cl_time_reorder_hilb_mean.get(r).resize(cl_n_particles.size());
			cl_time_reorder_hilb_dev.get(r).resize(cl_n_particles.size());
			for (size_t k = 0; k < cl_time_reorder_hilb_mean.get(r).size(); k++)
			{
				cl_time_reorder_hilb_mean.get(r).get(k).resize(cl_orders.size());
				cl_time_reorder_hilb_dev.get(r).get(k).resize(cl_orders.size());
			}
		}

		// Print test
		std::string str("Testing " + std::to_string(dim) + "D vector, Hilbert curve reordering, Cell-List");
		print_test_v(str);

		// For different r_cut
		for (size_t r = 0; r < cl_r_cutoff.size(); r++ )
		{
			Vcluster & v_cl = create_vcluster();

			// Cut-off radius
			float r_cut = cl_r_cutoff.get(r);

			// Number of particles
			size_t k = cl_k_start * v_cl.getProcessingUnits();

			//For different curve orders
			for ( size_t i = 0; i < cl_orders.size(); i++)
			{
				size_t m = cl_orders.get(i);
				size_t part = 0;

				for (size_t k_int = k ; k_int >= cl_k_min ; k_int/=2, part++ )
				{
					BOOST_TEST_CHECKPOINT( "Testing " << dim << "D vector with an Hilbert curve reordering k=" << k_int );

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

					vector_dist<dim,float, aggregate<float[dim]>, CartDecomposition<dim,float> > vd(k_int,box,bc,Ghost<dim,float>(r_cut));

					// Initialize a dist vector
					vd_initialize<dim>(vd, v_cl, k_int);

					//Reorder a vector

					double sum_reorder_mean = 0;
					double sum_reorder_dev = 0;

					openfpm::vector<double> measures;

					for (size_t h = 0 ; h < N_STAT_TEST; h++)
						measures.add(benchmark_reorder(vd,m));
					standard_deviation(measures,sum_reorder_mean,sum_reorder_dev);

					//Average reorder time
					cl_time_reorder_hilb_mean.get(r).get(part).get(i) = sum_reorder_mean;
					cl_time_reorder_hilb_dev.get(r).get(part).get(i) = sum_reorder_dev;

					vd.template ghost_get<0>();

					//Get cell list

					auto NN = vd.getCellList(r_cut);
					benchmark_get_celllist(NN,vd,r_cut);

					//Calculate forces

					double sum_fr_mean = 0;
					double sum_fr_dev = 0;
					measures.clear();

					for (size_t l = 0 ; l < N_STAT_TEST ; l++)
						measures.add(benchmark_calc_forces<dim>(NN,vd,r_cut));
					standard_deviation(measures,sum_fr_mean,sum_fr_dev);

					cl_time_hilb_mean.get(r).get(part).get(i) = sum_fr_mean;
					cl_time_hilb_dev.get(r).get(part).get(i) = sum_fr_dev;


					if (v_cl.getProcessUnitID() == 0)
						std::cout << "Cut-off = " << r_cut << ", Particles = " << k_int << ". Time to reorder: " << sum_reorder_mean << " dev: " << sum_reorder_dev << "      time calculate forces: " << sum_fr_mean << " dev: " << sum_fr_dev << std::endl;
				}
			}
		}
	}
}

/*! \brief Function for random cell list test
 *
 */
template<unsigned int dim> void cell_list_comp_reorder_random_benchmark(size_t cl_k_start,
		                                                                size_t cl_k_min,
																		openfpm::vector<float> & cl_r_cutoff,
																		openfpm::vector<size_t> & cl_n_particles,
																		openfpm::vector<openfpm::vector<double>> & cl_time_rand_mean,
																		openfpm::vector<openfpm::vector<double>> & cl_time_rand_dev,
																		openfpm::vector<openfpm::vector<double>> & cl_time_create_rand_mean,
																		openfpm::vector<openfpm::vector<double>> & cl_time_create_rand_dev)
{
	cl_time_rand_mean.resize(cl_r_cutoff.size());
	cl_time_create_rand_mean.resize(cl_r_cutoff.size());
	cl_time_rand_dev.resize(cl_r_cutoff.size());
	cl_time_create_rand_dev.resize(cl_r_cutoff.size());

	std::string str("Testing " + std::to_string(dim) + "D vector, no order, cell-list");
	print_test_v(str);

	{
		//For different r_cut
		for (size_t r = 0; r < cl_r_cutoff.size(); r++ )
		{
			Vcluster & v_cl = create_vcluster();

			//Cut-off radius
			float r_cut = cl_r_cutoff.get(r);

			//Number of particles
			size_t k = cl_k_start * v_cl.getProcessingUnits();

			//Counter number for amounts of particles
			size_t k_count = 1 + log2(k/cl_k_min);

			//For different number of particles
			for (size_t k_int = k ; k_int >= cl_k_min ; k_int/=2 )
			{
				BOOST_TEST_CHECKPOINT( "Testing " << dim << "D vector with a random cell list k=" << k_int );

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

				vector_dist<dim,float, aggregate<float[dim]>, CartDecomposition<dim,float> > vd(k_int,box,bc,Ghost<dim,float>(r_cut));

				// Initialize a dist vector
				vd_initialize<dim>(vd, v_cl, k_int);

				vd.template ghost_get<0>();

				//Get a cell list

				auto NN = vd.getCellList(r_cut);
				double sum_cl_mean = 0;
				double sum_cl_dev = 0;

				openfpm::vector<double> measures;
				for (size_t n = 0 ; n < N_STAT_TEST; n++)
					measures.add(benchmark_get_celllist(NN,vd,r_cut));
				standard_deviation(measures,sum_cl_mean,sum_cl_dev);
				//Average total time

				cl_time_create_rand_mean.get(r).add(sum_cl_mean);
				cl_time_create_rand_dev.get(r).add(sum_cl_dev);

				//Calculate forces

				double sum_fr_mean = 0;
				double sum_fr_dev = 0;

				measures.clear();
				for (size_t l = 0 ; l < N_STAT_TEST; l++)
					measures.add(benchmark_calc_forces<dim>(NN,vd,r_cut));
				standard_deviation(measures,sum_fr_mean,sum_fr_dev);

				cl_time_rand_mean.get(r).add(sum_fr_mean);
				cl_time_rand_dev.get(r).add(sum_fr_dev);

				if (v_cl.getProcessUnitID() == 0)
					std::cout << "Cut-off = " << r_cut << ", Particles = " << k_int << ". Time to create a cell-list: " << sum_cl_mean << " dev: " << sum_cl_dev << "    time to calculate forces: " << sum_fr_mean << " dev: " << sum_fr_dev << std::endl;
			}
		}
	}
}

/*! \brief Function for hilb cell list test
 *
 */
template<unsigned int dim> void cell_list_comp_reorder_hilbert_benchmark(size_t cl_k_start,
		                                                                 size_t cl_k_min,
																		 openfpm::vector<float> & cl_r_cutoff,
																		 openfpm::vector<size_t> & cl_n_particles,
																		 openfpm::vector<openfpm::vector<double>> & cl_time_force_mean,
																		 openfpm::vector<openfpm::vector<double>> & cl_time_force_dev,
																		 openfpm::vector<openfpm::vector<double>> & cl_time_create_mean,
																		 openfpm::vector<openfpm::vector<double>> & cl_time_create_dev)
{
	cl_time_force_mean.resize(cl_r_cutoff.size());
	cl_time_create_mean.resize(cl_r_cutoff.size());
	cl_time_force_dev.resize(cl_r_cutoff.size());
	cl_time_create_dev.resize(cl_r_cutoff.size());

	std::string str("Testing " + std::to_string(dim) + "D vector, Hilbert comp reorder, cell list");
	print_test_v(str);

	{
		//For different r_cut
		for (size_t r = 0; r < cl_r_cutoff.size(); r++ )
		{
			Vcluster & v_cl = create_vcluster();

			//Cut-off radius
			float r_cut = cl_r_cutoff.get(r);

			//Number of particles
			size_t k = cl_k_start * v_cl.getProcessingUnits();

			//Counter number for amounts of particles
			size_t k_count = 1 + log2(k/cl_k_min);

			//For different number of particles
			for (size_t k_int = k ; k_int >= cl_k_min ; k_int/=2 )
			{
				BOOST_TEST_CHECKPOINT( "Testing " << dim << "D vector with an Hilbert cell list k=" << k_int );

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

				vector_dist<dim,float, aggregate<float[dim]>, CartDecomposition<dim,float> > vd(k_int,box,bc,Ghost<dim,float>(r_cut));

				// Initialize a dist vector
				vd_initialize<dim>(vd, v_cl, k_int);

				vd.template ghost_get<0>();

				//Get a cell list hilb

				auto NN = vd.getCellList_hilb(r_cut);

				openfpm::vector<double> measures;

				double sum_cl_mean = 0;
				double sum_cl_dev = 0;
				for (size_t n = 0 ; n < N_VERLET_TEST; n++)
					measures.add(benchmark_get_celllist_hilb(NN,vd,r_cut));
				standard_deviation(measures,sum_cl_mean,sum_cl_dev);
				//Average total time
				cl_time_create_mean.get(r).add(sum_cl_mean);
				cl_time_create_dev.get(r).add(sum_cl_dev);

				//Calculate forces

				double sum_fr_mean = 0;
				double sum_fr_dev = 0;

				measures.clear();
				for (size_t l = 0 ; l < N_VERLET_TEST; l++)
					measures.add(benchmark_calc_forces_hilb<dim>(NN,vd,r_cut));
				standard_deviation(measures,sum_fr_mean,sum_fr_dev);

				cl_time_force_mean.get(r).add(sum_fr_mean);
				cl_time_force_dev.get(r).add(sum_fr_dev);

				if (v_cl.getProcessUnitID() == 0)
					std::cout << "Cut-off = " << r_cut << ", Particles = " << k_int << ". Time to create: " << sum_cl_mean << " dev: " << sum_cl_dev << " time to create a Cell-list: " << sum_fr_mean << " dev: " << sum_fr_dev << std::endl;
			}
		}
	}
}

////////////////////////// REPORT WRITING FUNCTIONS ///////////////////////////

/*! \brief Function for verlet performance report
 *
 */
template<unsigned int dim> void vd_verlet_performance_write_report(GoogleChart & cg,
																   openfpm::vector<float> & r_cutoff,
		                                                           openfpm::vector<size_t> & n_particles,
																   openfpm::vector<openfpm::vector<double>> time_force_mean,
																   openfpm::vector<openfpm::vector<double>> time_force_dev,
																   openfpm::vector<openfpm::vector<double>> time_create_mean,
																   openfpm::vector<openfpm::vector<double>> time_create_dev)
{
	// Get the test dir
	std::string file_mean(test_dir);
	std::string file_var(test_dir);
	file_mean += std::string("/openfpm_pdata/verlet_comp_create_mean_" + std::to_string(dim) + std::string("_ref"));
	file_var += std::string("/openfpm_pdata/verlet_comp_create_dev_" + std::to_string(dim) + std::string("_ref"));

	openfpm::vector<openfpm::vector<openfpm::vector<double>>> y_ref_create_mean;
	openfpm::vector<openfpm::vector<openfpm::vector<double>>> y_ref_create_dev;
	y_ref_create_mean.load(file_mean);
	y_ref_create_dev.load(file_var);

	// Get the test dir
	std::string file_mean2(test_dir);
	std::string file_var2(test_dir);
	file_mean2 += std::string("/openfpm_pdata/verlet_comp_force_mean_" + std::to_string(dim) + std::string("_ref"));
	file_var2 += std::string("/openfpm_pdata/verlet_comp_force_dev_" + std::to_string(dim) + std::string("_ref"));

	openfpm::vector<openfpm::vector<openfpm::vector<double>>> y_ref_force_mean;
	openfpm::vector<openfpm::vector<openfpm::vector<double>>> y_ref_force_dev;
	y_ref_force_mean.load(file_mean2);
	y_ref_force_dev.load(file_var2);

	// Speedup graphs data
	openfpm::vector<size_t> x;
	openfpm::vector<openfpm::vector<openfpm::vector<double>>> y;
	openfpm::vector<openfpm::vector<openfpm::vector<double>>> y_dev;
	openfpm::vector<std::string> yn;

	yn.add("Force verlet");
	for (size_t i = 0; i < n_particles.size() ; i++)
		x.add(n_particles.get(i));

	y.resize(time_force_mean.size());
	y_dev.resize(time_force_mean.size());
	for (size_t r = 0; r < time_force_mean.size(); r++)
	{
		y.get(r).resize(time_force_mean.get(r).size());
		y_dev.get(r).resize(time_force_mean.get(r).size());
		for (size_t k = 0; k < time_force_mean.get(r).size(); k++)
		{
			// Put the speedup
			y.get(r).get(k).add(time_force_mean.get(r).get(k));

			y_dev.get(r).get(k).add(time_force_dev.get(r).get(k));
		}
	}

	y.save("verlet_comp_force_mean_" + std::to_string(dim) + std::to_string("_ref"));
	y_dev.save("verlet_comp_force_dev_" + std::to_string(dim) + std::to_string("_ref"));

	if (y_ref_force_mean.size() != 0)
	{
		yn.clear();

		yn.add("Force verlet");
		yn.add("interval");
		yn.add("interval");

		y.clear();
		y.resize(time_force_mean.size());
		for (size_t r = 0; r < time_force_mean.size(); r++)
		{
			y.get(r).resize(time_force_mean.get(r).size());
			y_dev.get(r).resize(time_force_mean.get(r).size());
			for (size_t k = 0; k < time_force_mean.get(r).size(); k++)
			{
				// Put the speedup
				y.get(r).get(k).add(time_force_mean.get(r).get(k));
				y.get(r).get(k).add(y_ref_force_mean.get(r).get(k).get(0) - 3.0*y_ref_force_dev.get(r).get(k).get(0) );
				y.get(r).get(k).add(y_ref_force_mean.get(r).get(k).get(0) + 3.0*y_ref_force_dev.get(r).get(k).get(0) );
			}
		}
	}

	// Calculation time graphs data

	openfpm::vector<openfpm::vector<openfpm::vector<double>>> y2;
	openfpm::vector<openfpm::vector<openfpm::vector<double>>> y2_dev;
	openfpm::vector<std::string> yn2;

	yn2.add("Create a Verlet");

	y2.resize(time_create_mean.size());
	y2_dev.resize(time_create_mean.size());
	for (size_t r = 0; r < time_create_mean.size(); r++)
	{
		y2.get(r).resize(time_create_mean.get(r).size());
		y2_dev.get(r).resize(time_create_mean.get(r).size());
		for (size_t k = 0; k < time_create_mean.get(r).size(); k++)
		{
			// Put a random case total time
			y2.get(r).get(k).add(time_create_mean.get(r).get(k));

			y2_dev.get(r).get(k).add(time_create_dev.get(r).get(k));
		}
	}

	y2.save("verlet_comp_create_mean_" + std::to_string(dim) + std::to_string("_ref"));
	y2_dev.save("verlet_comp_create_dev_" + std::to_string(dim) + std::to_string("_ref"));

	if (y_ref_create_mean.size() != 0)
	{
		yn2.clear();
		yn2.add("Create a Verlet");
		yn2.add("interval");
		yn2.add("interval");

		y2.clear();

		y2.resize(time_create_mean.size());
		for (size_t r = 0; r < time_create_mean.size(); r++)
		{
			y2.get(r).resize(time_create_mean.get(r).size());
			for (size_t k = 0; k < time_create_mean.get(r).size(); k++)
			{
				// Put a random case total time
				y2.get(r).get(k).add(time_create_mean.get(r).get(k));

				y2.get(r).get(k).add(y_ref_create_mean.get(r).get(k).get(0) - 3.0*y_ref_create_dev.get(r).get(k).get(0) );
				y2.get(r).get(k).add(y_ref_create_mean.get(r).get(k).get(0) + 3.0*y_ref_create_dev.get(r).get(k).get(0) );
			}
		}
	}

	// Speedup graphs report

	// Google charts options
	GCoptions options;

	options.yAxis = std::string("Time (s)");
	options.xAxis = std::string("Number of particles");
	options.lineWidth = 2;
	options.more = GC_Y_LOG + "," + GC_ZOOM;

	std::string str("<h1>Verlet-list " + std::to_string(dim) + "-D performance test force calculation: </h1>");
	cg.addHTML(str);

	for (size_t i = 0; i < r_cutoff.size(); i++)
	{
		options.title = std::string("Verlet-list cut-off radius: " + std::to_string(r_cutoff.get(i)));
		cg.AddLinesGraph(x,y.get(i),yn,options);
	}

	// Calculation time graphs report

	// Google charts options
	GCoptions options2;

	options2.yAxis = std::string("Time to construct a verlet-list (s)");
	options2.xAxis = std::string("Number of particles");
	options2.lineWidth = 2;
	options2.more = GC_ZOOM;

	std::string str2("<h2>2) Total calculation time</h2>");
	cg.addHTML(str2);

	for (size_t i = 0; i < r_cutoff.size(); i++)
	{
		options2.title = std::string("Cell-list performance, cut-off radius: " + std::to_string(r_cutoff.get(i)));
		cg.AddLinesGraph(x,y2.get(i),yn2,options2);
	}
}


/*! \brief Function for cell list performance report
 *
 */
template<unsigned int dim> void cell_list_part_reorder_report(GoogleChart & cg,
		                                                      size_t n_moving,
		                                                      openfpm::vector<float> & cl_r_cutoff,
															  openfpm::vector<size_t> & cl_n_particles,
															  openfpm::vector<size_t> cl_orders,
															  openfpm::vector<openfpm::vector<openfpm::vector<double>>> cl_time_hilb_mean,
															  openfpm::vector<openfpm::vector<double>> cl_time_rand_mean,
															  openfpm::vector<openfpm::vector<openfpm::vector<double>>> cl_time_reorder_mean,
															  openfpm::vector<openfpm::vector<openfpm::vector<double>>> cl_time_hilb_dev,
															  openfpm::vector<openfpm::vector<double>> cl_time_rand_dev,
															  openfpm::vector<openfpm::vector<openfpm::vector<double>>> cl_time_reorder_dev)
{
	openfpm::vector<size_t> x;

	for (size_t i = 0; i < cl_n_particles.size() ; i++)
		x.add(cl_n_particles.get(i));

	// Speedup graphs data

	cl_part_time<dim>(cg,cl_n_particles,cl_r_cutoff,cl_orders,cl_time_hilb_mean,cl_time_rand_mean,cl_time_hilb_dev,cl_time_rand_dev);
	cl_part_reorder_time<dim>(cg,cl_n_particles,cl_r_cutoff,cl_orders,cl_time_reorder_mean,cl_time_reorder_dev);
}

/*! \brief Function for cell list hilb performance report
 *
 */
template<unsigned int dim> void cell_list_comp_reorder_report(GoogleChart & cg,
		                                                      openfpm::vector<float> & cl_r_cutoff,
															  openfpm::vector<size_t> & cl_n_particles,
															  openfpm::vector<openfpm::vector<double>> & cl_time_hilb_mean,
															  openfpm::vector<openfpm::vector<double>> & cl_time_rand_mean,
															  openfpm::vector<openfpm::vector<double>> & cl_time_hilb_dev,
															  openfpm::vector<openfpm::vector<double>> & cl_time_rand_dev,
															  openfpm::vector<openfpm::vector<double>> & cl_time_create_hilb_mean,
															  openfpm::vector<openfpm::vector<double>> & cl_time_create_rand_mean,
															  openfpm::vector<openfpm::vector<double>> & cl_time_create_hilb_dev,
															  openfpm::vector<openfpm::vector<double>> & cl_time_create_rand_dev)
{
	cl_comp_normal_vs_hilbert_force_time<dim>(cg,cl_n_particles,cl_r_cutoff,cl_time_hilb_mean,cl_time_rand_mean,cl_time_hilb_dev,cl_time_rand_dev);
	cl_comp_normal_vs_hilbert_create_time<dim>(cg,cl_n_particles,cl_r_cutoff,cl_time_create_hilb_mean,cl_time_create_rand_mean,cl_time_create_hilb_dev,cl_time_create_rand_dev);
}

#endif /* SRC_VECTOR_VECTOR_DIST_PERFORMANCE_UTIL_HPP_ */
