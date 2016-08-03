/*
 * vector_dist_performance_util.hpp
 *
 *  Created on: Jun 17, 2016
 *      Author: Yaroslav Zaluzhnyi
 */

#ifndef SRC_VECTOR_VECTOR_DIST_PERFORMANCE_UTIL_HPP_
#define SRC_VECTOR_VECTOR_DIST_PERFORMANCE_UTIL_HPP_

#include "Plot/GoogleChart.hpp"

/*! \brief Print out only ones (no matter how many processors involved)
 *
 * \param test, sz Data to print out
 */
void print_test_v(std::string test, size_t sz)
{
	if (create_vcluster().getProcessUnitID() == 0)
		std::cout << "\n" << test << " " << sz << "\n";
}

/*! \brief Initialize a distributed vector
 *
 * \param vd Distributed vector
 * \param v_cl Global vcluster
 * \param k_int Number of particles
 */
template<unsigned int dim, typename v_dist> void vd_initialize(v_dist & vd, Vcluster & v_cl, size_t k_int)
{
	// The random generator engine
	std::default_random_engine eg(v_cl.getProcessUnitID()*4313);
	std::uniform_real_distribution<float> ud(0.0f, 1.0f);

	//! [Create a vector of random elements on each processor 2D]

	auto it = vd.getIterator();

	while (it.isNext())
	{
		auto key = it.get();

		for (size_t i = 0; i < dim; i++)
			vd.getPos(key)[i] = ud(eg);

		++it;
	}

	vd.map();
}

/*! \brief Initialize 2 distributed vectors with equally positioned particles
 *
 * \param vd, vd2 Distributed vectors
 * \param v_cl Global vcluster
 * \param k_int Number of particles
 */
template<unsigned int dim, typename v_dist> void vd_initialize_double(v_dist & vd,v_dist & vd2, Vcluster & v_cl, size_t k_int)
{
	// The random generator engine
	std::default_random_engine eg(v_cl.getProcessUnitID()*4313);
	std::uniform_real_distribution<float> ud(0.0f, 1.0f);

	//! [Create a vector of random elements on each processor 2D]

	auto it = vd.getIterator();

	while (it.isNext())
	{
		auto key = it.get();

		for (size_t i = 0; i < dim; i++)
		{
			vd.getPos(key)[i] = ud(eg);
			vd2.getPos(key)[i] = vd.getPos(key)[i];
		}

		++it;
	}

	vd.map();
	vd2.map();
}

/*! \brief Calculate and put particles' forces
 *
 * \param NN Cell list
 * \param vd Distributed vector
 * \param r_cut Cut-off radius
 */
template<unsigned int dim, size_t prp = 0, typename T, typename V> void calc_forces(T & NN, V & vd, float r_cut)
{
	auto it_v = vd.getDomainIterator();

	float sum[dim];

	for (size_t i = 0; i < dim; i++)
		sum[i] = 0;

	while (it_v.isNext())
	{
		//key
		vect_dist_key_dx key = it_v.get();

    	// Get the position of the particles
		Point<dim,float> p = vd.getPos(key);

		for (size_t i = 0; i < dim; i++)
			sum[i] = 0;

    	// Get the neighborhood of the particle
    	auto cell_it = NN.template getNNIterator<NO_CHECK>(NN.getCell(p));

    	while(cell_it.isNext())
    	{
    		auto nnp = cell_it.get();

    		// p != q
    		if (nnp == key.getKey())
    		{
    			++cell_it;
    			continue;
    		}

    		Point<dim,float> q = vd.getPos(nnp);

    		if (p.distance2(q) < r_cut*r_cut)
    		{
				//Calculate the forces
    			float num[dim];
    			for (size_t i = 0; i < dim; i++)
    				num[i] = vd.getPos(key)[i] - vd.getPos(nnp)[i];

    			float denom = 0;
    			for (size_t i = 0; i < dim; i++)
					denom += num[i] * num[i];

    			float res[dim];
    			for (size_t i = 0; i < dim; i++)
					res[i] = num[i] / denom;

    			for (size_t i = 0; i < dim; i++)
					sum[i] += res[i];
    		}
			//Next particle in a cell
			++cell_it;
		}

		//Put the forces
		for (size_t i = 0; i < dim; i++)
			vd.template getProp<prp>(key)[i] += sum[i];

		//Next particle in cell list
		++it_v;
	}
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
	vd.getVerlet(verlet,r_cut);

	t.stop();

	return t.getwct();
}

/*! \brief Calculate and put particles' forces
 *
 * \param NN Cell list hilbert
 * \param vd Distributed vector
 * \param r_cut Cut-off radius
 */
template<unsigned int dim, size_t prp = 0, typename T, typename V> void calc_forces_hilb(T & NN, V & vd, float r_cut)
{
	auto it_cl = NN.getIterator();

	float sum[dim];

	for (size_t i = 0; i < dim; i++)
		sum[i] = 0;

	while (it_cl.isNext())
	{
		//key
		auto key = it_cl.get();

    	// Get the position of the particles
		Point<dim,float> p = vd.getPos(key);

		for (size_t i = 0; i < dim; i++)
			sum[i] = 0;

    	// Get the neighborhood of the particle
    	auto cell_it = NN.template getNNIterator<NO_CHECK>(NN.getCell(p));

    	while(cell_it.isNext())
    	{
    		auto nnp = cell_it.get();

    		// p != q
    		if (nnp == key)
    		{
    			++cell_it;
    			continue;
    		}

    		Point<dim,float> q = vd.getPos(nnp);

    		if (p.distance2(q) < r_cut*r_cut)
    		{
				//Calculate the forces
    			float num[dim];
    			for (size_t i = 0; i < dim; i++)
    				num[i] = vd.getPos(key)[i] - vd.getPos(nnp)[i];

    			float denom = 0;
    			for (size_t i = 0; i < dim; i++)
					denom += num[i] * num[i];

    			float res[dim];
    			for (size_t i = 0; i < dim; i++)
					res[i] = num[i] / denom;

    			for (size_t i = 0; i < dim; i++)
					sum[i] += res[i];
    		}
			//Next particle in a cell
			++cell_it;
		}

		//Put the forces
		for (size_t i = 0; i < dim; i++)
			vd.template getProp<prp>(key)[i] += sum[i];

		//Next particle in cell list
		++it_cl;
	}
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
template<unsigned int dim> void vd_verlet_random_benchmark(size_t k_start, size_t k_min, double ghost_part, openfpm::vector<float> & r_cutoff, openfpm::vector<size_t> & n_particles, openfpm::vector<openfpm::vector<double>> & time_rand, openfpm::vector<openfpm::vector<double>> & time_total_rand)
{
	time_rand.resize(r_cutoff.size());
	time_total_rand.resize(r_cutoff.size());

	{
		//For different r_cut
		for (size_t r = 0; r < r_cutoff.size(); r++ )
		{
			Vcluster & v_cl = create_vcluster();

			//Cut-off radius
			float r_cut = r_cutoff.get(r);

			//Number of particles
			size_t k = k_start * v_cl.getProcessingUnits();

			std::string str("Testing " + std::to_string(dim) + "D vector without an Hilbert curve reordering k<=");

			print_test_v(str,k);

			std::cout << std::endl << "Cut-off raidus is " << r_cut << std::endl;

			//Counter number for amounts of particles
			size_t k_count = 1 + log2(k/k_min);

			for (size_t k_int = k ; k_int >= k_min ; k_int/=2 )
			{
				BOOST_TEST_CHECKPOINT( "Testing " << dim << "D vector without an Hilbert curve reordering k=" << k_int );

				std::cout << std::endl << "Number of particles: " << k_int << std::endl;

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

				vector_dist<dim,float, aggregate<float[dim]>, CartDecomposition<dim,float> > vd(k_int,box,bc,Ghost<dim,float>(ghost_part));

				// Initialize a dist vector
				vd_initialize<dim>(vd, v_cl, k_int);

				vd.template ghost_get<0>();

				//Get verlet list

				size_t n = 0;
				double sum_verlet = 0;

				for ( ; n < 3; n++)
				{
					sum_verlet += benchmark_get_verlet(vd,r_cut);
				}
				std::cout << "Average of " << n << " calculations: verlet time = " << sum_verlet / n << std::endl;

				//Calculate forces

				auto NN = vd.getCellList(r_cut);
				double sum_forces = 0;
				size_t l = 0;

				for ( ; l < 3; l++)
				{
					sum_forces += benchmark_calc_forces<dim>(NN,vd,r_cut);
				}
				std::cout << "Average of " << l << " calculations: forces time = " << sum_forces / l << std::endl;
				time_rand.get(r).add(sum_forces / l);

				//Average total time
				time_total_rand.get(r).add(sum_forces / l + sum_verlet / n);
				std::cout << "Average total time = " << sum_forces / l + sum_verlet / n << std::endl;
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

	// For different r_cut
	for (size_t r = 0; r < r_cutoff.size(); r++ )
	{
		Vcluster & v_cl = create_vcluster();

		//Cut-off radius
		float r_cut = r_cutoff.get(r);

		// Number of particles
		size_t k = k_start * v_cl.getProcessingUnits();

		std::string str("Testing " + std::to_string(dim) + "D vector with an Hilbert curve reordering k<=");

		print_test_v(str,k);

		std::cout << std::endl << "Cut-off raidus is " << r_cut << std::endl;

		//For different curve orders
		for ( size_t i = 0; i < orders.size(); i++)
		{
			size_t m = orders.get(i);

			std::cout << std::endl << "Order of a curve: " << orders.get(i) << std::endl;

			size_t part = 0;

			for (size_t k_int = k ; k_int >= k_min ; k_int/=2, part++ )
			{
				BOOST_TEST_CHECKPOINT( "Testing " << dim << "D vector with an Hilbert curve reordering k=" << k_int );

				std::cout << std::endl << "Number of particles: " << k_int << std::endl;

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
				size_t h = 0;

				for ( ; h < 3; h++)
				{
					sum_reorder += benchmark_reorder(vd,m);
				}
				std::cout << "Average of " << h << " calculations: reordering time = " << sum_reorder / h << std::endl;

				//Get verlet list

				size_t n = 0;
				double sum_verlet = 0;

				for ( ; n < 3; n++)
				{
					sum_verlet += benchmark_get_verlet(vd,r_cut);
				}
				std::cout << "Average of " << n << " calculations: verlet time = " << sum_verlet / n << std::endl;

				//Calculate forces

				auto NN = vd.getCellList(r_cut);
				double sum_forces = 0;
				size_t l = 0;

				for ( ; l < 3; l++)
				{
					sum_forces += benchmark_calc_forces<dim>(NN,vd,r_cut);
				}
				std::cout << "Average of " << l << " calculations: forces time = " << sum_forces / l << std::endl;
				time_hilb.get(r).get(part).get(i) = sum_forces / l;

				//Average total time
				std::cout << "Average total time = " << sum_forces / l + sum_verlet / n + sum_reorder / h << std::endl;
				time_total_hilb.get(r).get(part).get(i) = sum_forces / l + sum_verlet / n + sum_reorder / h;
			}
		}
	}
}

/*! \brief Function for cell list test without an Hilbert curve reordering (unordered positioning)
 *
 */
template<unsigned int dim> void vd_cl_random_benchmark(size_t cl_k_start, size_t cl_k_min, double ghost_part, openfpm::vector<float> & cl_r_cutoff, openfpm::vector<size_t> & cl_n_particles, openfpm::vector<openfpm::vector<double>> & cl_time_rand, openfpm::vector<openfpm::vector<double>> & cl_time_total_rand)
{
	cl_time_rand.resize(cl_r_cutoff.size());
	cl_time_total_rand.resize(cl_r_cutoff.size());

	{
		//For different r_cut
		for (size_t r = 0; r < cl_r_cutoff.size(); r++ )
		{
			Vcluster & v_cl = create_vcluster();

			//Cut-off radius
			float r_cut = cl_r_cutoff.get(r);

			//Number of particles
			size_t k = cl_k_start * v_cl.getProcessingUnits();

			std::string str("Testing " + std::to_string(dim) + "D vector without an Hilbert curve reordering k<=");

			print_test_v(str,k);

			std::cout << std::endl << "Cut-off raidus is " << r_cut << std::endl;

			//Counter number for amounts of particles
			size_t k_count = 1 + log2(k/cl_k_min);

			//For different number of particles
			for (size_t k_int = k ; k_int >= cl_k_min ; k_int/=2 )
			{
				BOOST_TEST_CHECKPOINT( "Testing " << dim << "D vector without an Hilbert curve reordering k=" << k_int );

				std::cout << std::endl << "Number of particles: " << k_int << std::endl;

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

				vector_dist<dim,float, aggregate<float[dim]>, CartDecomposition<dim,float> > vd(k_int,box,bc,Ghost<dim,float>(ghost_part));

				// Initialize a dist vector
				vd_initialize<dim>(vd, v_cl, k_int);

				vd.template ghost_get<0>();

				//Get a cell list

				auto NN = vd.getCellList(r_cut);
				size_t n = 0;
				double sum_cl = 0;

				for ( ; n < 3; n++)
				{
					sum_cl += benchmark_get_celllist(NN,vd,r_cut);
				}
				std::cout << "Average of " << n << " calculations: celllist time = " << sum_cl / n << std::endl;

				//Calculate forces

				double sum_forces = 0;
				size_t l = 0;

				for ( ; l < 3; l++)
				{
					sum_forces += benchmark_calc_forces<dim>(NN,vd,r_cut);
				}
				std::cout << "Average of " << l << " calculations: forces time = " << sum_forces / l << std::endl;
				cl_time_rand.get(r).add(sum_forces / l);

				//Average total time
				cl_time_total_rand.get(r).add(sum_forces / l + sum_cl / n);
				std::cout << "Average total time = " << sum_forces / l + sum_cl / n << std::endl;
			}
		}
	}
}

/*! \brief Function for cell list test with an Hilbert curve reordering
 *
 */
template<unsigned int dim> void vd_cl_hilbert_benchmark(size_t cl_k_start, size_t cl_k_min, double ghost_part, size_t n_moving, double dist, openfpm::vector<float> & cl_r_cutoff, openfpm::vector<size_t> & cl_n_particles, openfpm::vector<size_t> &cl_orders, openfpm::vector<openfpm::vector<openfpm::vector<double>>> &cl_time_hilb, openfpm::vector<openfpm::vector<openfpm::vector<double>>> &cl_time_total_hilb, openfpm::vector<openfpm::vector<openfpm::vector<openfpm::vector<double>>>> &cl_time_hilb_moved)
{
	{
		cl_time_hilb.resize(cl_r_cutoff.size());
		for (size_t r = 0; r < cl_time_hilb.size(); r++)
		{
			cl_time_hilb.get(r).resize(cl_n_particles.size());
			for (size_t k = 0; k < cl_time_hilb.get(r).size(); k++)
			{
				cl_time_hilb.get(r).get(k).resize(cl_orders.size());
			}
		}

		cl_time_hilb_moved.resize(n_moving);
		for (size_t d = 0; d < cl_time_hilb_moved.size(); d++)
		{
			cl_time_hilb_moved.get(d).resize(cl_r_cutoff.size());
			for (size_t r = 0; r < cl_time_hilb_moved.get(d).size(); r++)
			{
				cl_time_hilb_moved.get(d).get(r).resize(cl_n_particles.size());
				for (size_t k = 0; k < cl_time_hilb_moved.get(d).get(r).size(); k++)
				{
					cl_time_hilb_moved.get(d).get(r).get(k).resize(cl_orders.size());
				}
			}
		}

		cl_time_total_hilb.resize(cl_r_cutoff.size());
		for (size_t r = 0; r < cl_time_total_hilb.size(); r++)
		{
			cl_time_total_hilb.get(r).resize(cl_n_particles.size());
			for (size_t k = 0; k < cl_time_total_hilb.get(r).size(); k++)
			{
				cl_time_total_hilb.get(r).get(k).resize(cl_orders.size());
			}
		}

		// For different r_cut
		for (size_t r = 0; r < cl_r_cutoff.size(); r++ )
		{
			Vcluster & v_cl = create_vcluster();

			// Cut-off radius
			float r_cut = cl_r_cutoff.get(r);

			// Number of particles
			size_t k = cl_k_start * v_cl.getProcessingUnits();

			std::string str("Testing " + std::to_string(dim) + "D vector with an Hilbert curve reordering k<=");

			print_test_v(str,k);

			std::cout << std::endl << "Cut-off raidus is " << r_cut << std::endl;

			//For different curve orders
			for ( size_t i = 0; i < cl_orders.size(); i++)
			{
				size_t m = cl_orders.get(i);

				std::cout << std::endl << "Order of a curve: " << cl_orders.get(i) << std::endl;

				size_t part = 0;

				for (size_t k_int = k ; k_int >= cl_k_min ; k_int/=2, part++ )
				{
					BOOST_TEST_CHECKPOINT( "Testing " << dim << "D vector with an Hilbert curve reordering k=" << k_int );

					std::cout << std::endl << "Number of particles: " << k_int << std::endl;

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

					//Reorder a vector

					double sum_reorder = 0;
					size_t h = 0;

					for ( ; h < 3; h++)
					{
						sum_reorder += benchmark_reorder(vd,m);
					}
					std::cout << "Average of " << h << " calculations: reordering time = " << sum_reorder / h << std::endl;

					vd.template ghost_get<0>();

					//Get cell list

					auto NN = vd.getCellList(r_cut);
					size_t n = 0;
					double sum_cl = 0;

					for ( ; n < 3; n++)
					{
						sum_cl += benchmark_get_celllist(NN,vd,r_cut);
					}
					std::cout << "Average of " << n << " calculations: celllist time = " << sum_cl / n << std::endl;

					//Calculate forces

					double sum_forces = 0;
					size_t l = 0;

					for ( ; l < 3; l++)
					{
						sum_forces += benchmark_calc_forces<dim>(NN,vd,r_cut);
					}
					std::cout << "Average of " << l << " calculations: forces time = " << sum_forces / l << std::endl;
					cl_time_hilb.get(r).get(part).get(i) = sum_forces / l;

					//Average total time
					std::cout << "Average total time = " << sum_forces / l + sum_cl / n + sum_reorder / h << std::endl;
					cl_time_total_hilb.get(r).get(part).get(i) = sum_forces / l + sum_cl / n + sum_reorder / h;

					//Move particles
					for ( size_t d = 0; d < n_moving; d++)
					{
						move_particles<dim>(vd,dist);

						vd.map();

						vd.template ghost_get<0>();

						auto NN = vd.getCellList(r_cut);

						//Calculate forces

						double sum_forces_moved = 0;
						size_t j = 0;

						for ( ; j < 3; j++)
						{
							sum_forces_moved += benchmark_calc_forces<dim>(NN,vd,r_cut);
						}
						std::cout << "Average of " << j << " calculations: forces time after moving = " << sum_forces_moved / j << ", iteration " << d+1 << std::endl;
						cl_time_hilb_moved.get(d).get(r).get(part).get(i) = sum_forces_moved / j;
					}
				}
			}
		}
	}
}

/*! \brief Function for random cell list test
 *
 */
template<unsigned int dim> void vd_celllist_random_benchmark(size_t cl_k_start, size_t cl_k_min, double ghost_part, openfpm::vector<float> & cl_r_cutoff, openfpm::vector<size_t> & cl_n_particles, openfpm::vector<openfpm::vector<double>> & cl_time_rand, openfpm::vector<openfpm::vector<double>> & cl_time_total_rand)
{
	cl_time_rand.resize(cl_r_cutoff.size());
	cl_time_total_rand.resize(cl_r_cutoff.size());

	{
		//For different r_cut
		for (size_t r = 0; r < cl_r_cutoff.size(); r++ )
		{
			Vcluster & v_cl = create_vcluster();

			//Cut-off radius
			float r_cut = cl_r_cutoff.get(r);

			//Number of particles
			size_t k = cl_k_start * v_cl.getProcessingUnits();

			std::string str("Testing " + std::to_string(dim) + "D vector with a random cell list k<=");

			print_test_v(str,k);

			std::cout << std::endl << "Cut-off raidus is " << r_cut << std::endl;

			//Counter number for amounts of particles
			size_t k_count = 1 + log2(k/cl_k_min);

			//For different number of particles
			for (size_t k_int = k ; k_int >= cl_k_min ; k_int/=2 )
			{
				BOOST_TEST_CHECKPOINT( "Testing " << dim << "D vector with a random cell list k=" << k_int );

				std::cout << std::endl << "Number of particles: " << k_int << std::endl;

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

				vector_dist<dim,float, aggregate<float[dim]>, CartDecomposition<dim,float> > vd(k_int,box,bc,Ghost<dim,float>(ghost_part));

				// Initialize a dist vector
				vd_initialize<dim>(vd, v_cl, k_int);

				vd.template ghost_get<0>();

				//Get a cell list

				auto NN = vd.getCellList(r_cut);
				size_t n = 0;
				double sum_cl = 0;

				for ( ; n < 3; n++)
				{
					sum_cl += benchmark_get_celllist(NN,vd,r_cut);
				}
				std::cout << "Average of " << n << " calculations: celllist time = " << sum_cl / n << std::endl;

				//Calculate forces

				double sum_forces = 0;
				size_t l = 0;

				for ( ; l < 3; l++)
				{
					sum_forces += benchmark_calc_forces<dim>(NN,vd,r_cut);
				}
				std::cout << "Average of " << l << " calculations: forces time = " << sum_forces / l << std::endl;
				cl_time_rand.get(r).add(sum_forces / l);

				//Average total time
				cl_time_total_rand.get(r).add(sum_forces / l + sum_cl / n);
				std::cout << "Average total time = " << sum_forces / l + sum_cl / n << std::endl;
			}
		}
	}
}

/*! \brief Function for hilb cell list test
 *
 */
template<unsigned int dim> void vd_celllist_hilbert_benchmark(size_t cl_k_start, size_t cl_k_min, double ghost_part, openfpm::vector<float> & cl_r_cutoff, openfpm::vector<size_t> & cl_n_particles, openfpm::vector<openfpm::vector<double>> & cl_time_hilb, openfpm::vector<openfpm::vector<double>> & cl_time_total_hilb)
{
	cl_time_hilb.resize(cl_r_cutoff.size());
	cl_time_total_hilb.resize(cl_r_cutoff.size());

	{
		//For different r_cut
		for (size_t r = 0; r < cl_r_cutoff.size(); r++ )
		{
			Vcluster & v_cl = create_vcluster();

			//Cut-off radius
			float r_cut = cl_r_cutoff.get(r);

			//Number of particles
			size_t k = cl_k_start * v_cl.getProcessingUnits();

			std::string str("Testing " + std::to_string(dim) + "D vector with an Hilbert cell list k<=");

			print_test_v(str,k);

			std::cout << std::endl << "Cut-off raidus is " << r_cut << std::endl;

			//Counter number for amounts of particles
			size_t k_count = 1 + log2(k/cl_k_min);

			//For different number of particles
			for (size_t k_int = k ; k_int >= cl_k_min ; k_int/=2 )
			{
				BOOST_TEST_CHECKPOINT( "Testing " << dim << "D vector with an Hilbert cell list k=" << k_int );

				std::cout << std::endl << "Number of particles: " << k_int << std::endl;

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

				vector_dist<dim,float, aggregate<float[dim]>, CartDecomposition<dim,float> > vd(k_int,box,bc,Ghost<dim,float>(ghost_part));

				// Initialize a dist vector
				vd_initialize<dim>(vd, v_cl, k_int);

				vd.template ghost_get<0>();

				//Get a cell list hilb

				auto NN = vd.getCellList_hilb(r_cut);

				size_t n = 0;
				double sum_cl = 0;

				for ( ; n < 3; n++)
				{
					sum_cl += benchmark_get_celllist_hilb(NN,vd,r_cut);
				}
				std::cout << "Average of " << n << " calculations: celllist time = " << sum_cl / n << std::endl;

				//Calculate forces

				double sum_forces = 0;
				size_t l = 0;

				for ( ; l < 3; l++)
				{
					sum_forces += benchmark_calc_forces_hilb<dim>(NN,vd,r_cut);
				}
				std::cout << "Average of " << l << " calculations: forces time = " << sum_forces / l << std::endl;
				cl_time_hilb.get(r).add(sum_forces / l);

				//Average total time
				cl_time_total_hilb.get(r).add(sum_forces / l + sum_cl / n);
				std::cout << "Average total time = " << sum_forces / l + sum_cl / n << std::endl;
			}
		}
	}
}

////////////////////////// REPORT WRITING FUNCTIONS ///////////////////////////

/*! \brief Function for verlet performance report
 *
 */
template<unsigned int dim> void vd_verlet_performance_write_report(openfpm::vector<float> & r_cutoff,openfpm::vector<size_t> & n_particles,openfpm::vector<size_t> orders,openfpm::vector<openfpm::vector<openfpm::vector<double>>> time_hilb,openfpm::vector<openfpm::vector<double>> time_rand,openfpm::vector<openfpm::vector<openfpm::vector<double>>> time_total_hilb,openfpm::vector<openfpm::vector<double>> time_total_rand)
{
	// Speedup graphs data
	openfpm::vector<size_t> x;
	openfpm::vector<openfpm::vector<openfpm::vector<double>>> y;
	openfpm::vector<std::string> yn;

	for (size_t i = 0; i < n_particles.size() ; i++)
		x.add(n_particles.get(i));

	for (size_t i = 0; i < orders.size(); i++)
		yn.add("Order of: " + std::to_string(orders.get(i)));

	y.resize(time_hilb.size());
	for (size_t r = 0; r < time_hilb.size(); r++)
	{
		y.get(r).resize(time_hilb.get(r).size());
		for (size_t k = 0; k < time_hilb.get(r).size(); k++)
		{
			for (size_t m = 0; m < time_hilb.get(r).get(k).size(); m++)
			{
				// Put the speedup
				y.get(r).get(k).add(time_rand.get(r).get(k)/time_hilb.get(r).get(k).get(m));
			}
		}
	}


	// Calculation time graphs data

	openfpm::vector<openfpm::vector<openfpm::vector<double>>> y2;
	openfpm::vector<std::string> yn2;

	yn2.add("Random");
	for (size_t i = 0; i < orders.size(); i++)
		yn2.add("Order of: " + std::to_string(orders.get(i)));

	y2.resize(time_total_hilb.size());
	for (size_t r = 0; r < time_total_hilb.size(); r++)
	{
		y2.get(r).resize(time_total_hilb.get(r).size());
		for (size_t k = 0; k < time_total_hilb.get(r).size(); k++)
		{
			// Put a random case total time
			y2.get(r).get(k).add(time_total_rand.get(r).get(k));
			for (size_t m = 0; m < time_total_hilb.get(r).get(k).size(); m++)
			{
				// Put an Hilbert case total time
				y2.get(r).get(k).add(time_total_hilb.get(r).get(k).get(m));
			}
		}
	}

	// Speedup graphs report

	// Google charts options
	GCoptions options;

	options.yAxis = std::string("Speedup (times)");
	options.xAxis = std::string("Number of particles");
	options.lineWidth = 2.5;

	GoogleChart cg;

	std::string str("<h1>Distributed " + std::to_string(dim) + "-D vector performance tests: </h1>");
	str += "<h2> 1) Speedup between an unordered positioning and an Hilbert curve positioning of particles</h2>";
	str += "We create a distributed vector (VD) of randomly positioned in a " + std::to_string(dim) + "D-box particles. Then we get a verlet list of VD, with a certain cut-off radius. After that we calculate the forces of each particle. Later "
			"a VD is reordered according to an Hilbert curve, then we calculate forces again and compare the forces calculation time for an unordered and an Hilbert curve cases. The speedup is calculated and shown on graphs below, depending on different numbers of particles.";
	cg.addHTML(str);

	for (size_t i = 0; i < r_cutoff.size(); i++)
	{
		options.title = std::string("Distributed vector performance speedup, cut-off radius: " + std::to_string(r_cutoff.get(i)));
		cg.AddLinesGraph(x,y.get(i),yn,options);
	}

	// Calculation time graphs report

	// Google charts options
	GCoptions options2;

	options2.yAxis = std::string("Total calculation time (s)");
	options2.xAxis = std::string("Number of particles");
	options2.lineWidth = 2.5;

	std::string str2("<h2>2) Total calculation time</h2>");
	str2 += "We count a total calculation time. In the case of unordered positioning it is: verlet list creation time + forces calculation time; in the case of an Hilbert curve positioning: reordering time + verlet list creation time + forces calculation time."
			"The total calculation time is shown on graphs below, depending on different numbers of particles.";
	cg.addHTML(str2);

	for (size_t i = 0; i < r_cutoff.size(); i++)
	{
		options2.title = std::string("Distributed vector performance, cut-off radius: " + std::to_string(r_cutoff.get(i)));
		cg.AddLinesGraph(x,y2.get(i),yn2,options2);
	}

	cg.write("Vect_dist_verlet_perf_" + std::to_string(dim) + "D.html");
}


/*! \brief Function for cell list performance report
 *
 */
template<unsigned int dim> void vd_cl_performance_write_report(size_t n_moving,openfpm::vector<float> & cl_r_cutoff,openfpm::vector<size_t> & cl_n_particles,openfpm::vector<size_t> cl_orders,openfpm::vector<openfpm::vector<openfpm::vector<double>>> cl_time_hilb,openfpm::vector<openfpm::vector<double>> cl_time_rand,openfpm::vector<openfpm::vector<openfpm::vector<double>>> cl_time_total_hilb,openfpm::vector<openfpm::vector<double>> cl_time_total_rand,openfpm::vector<openfpm::vector<openfpm::vector<openfpm::vector<double>>>> cl_time_hilb_moved)
{
	// Speedup graphs data

	openfpm::vector<size_t> x;
	openfpm::vector<openfpm::vector<openfpm::vector<double>>> y;
	openfpm::vector<std::string> yn;

	for (size_t i = 0; i < cl_n_particles.size() ; i++)
		x.add(cl_n_particles.get(i));

	for (size_t i = 0; i < cl_orders.size(); i++)
		yn.add("Order of: " + std::to_string(cl_orders.get(i)));

	y.resize(cl_time_hilb.size());
	for (size_t r = 0; r < cl_time_hilb.size(); r++)
	{
		y.get(r).resize(cl_time_hilb.get(r).size());
		for (size_t k = 0; k < cl_time_hilb.get(r).size(); k++)
		{
			for (size_t m = 0; m < cl_time_hilb.get(r).get(k).size(); m++)
			{
				// Put a speedup
				y.get(r).get(k).add(cl_time_rand.get(r).get(k) / cl_time_hilb.get(r).get(k).get(m));
			}
		}
	}


	// Calculation time graphs data

	openfpm::vector<openfpm::vector<openfpm::vector<double>>> y2;
	openfpm::vector<std::string> yn2;

	yn2.add("Random");
	for (size_t i = 0; i < cl_orders.size(); i++)
		yn2.add("Order of: " + std::to_string(cl_orders.get(i)));

	y2.resize(cl_time_total_hilb.size());
	for (size_t r = 0; r < cl_time_total_hilb.size(); r++)
	{
		y2.get(r).resize(cl_time_total_hilb.get(r).size());
		for (size_t k = 0; k < cl_time_total_hilb.get(r).size(); k++)
		{
			// Put a random case total time
			y2.get(r).get(k).add(cl_time_total_rand.get(r).get(k));
			for (size_t m = 0; m < cl_time_total_hilb.get(r).get(k).size(); m++)
			{
				// Put an Hilbert case total time
				y2.get(r).get(k).add(cl_time_total_hilb.get(r).get(k).get(m));
			}
		}
	}

	// Speedup difference graphs data

	openfpm::vector<size_t> x3;
	openfpm::vector<openfpm::vector<openfpm::vector<openfpm::vector<double>>>> y3;
	openfpm::vector<std::string> yn3;

	//'+1' for '0 moving step' case (before moving)
	for (size_t i = 0; i < n_moving + 1 ; i++)
		x3.add(i);

	yn3.add("No speedup line");
	for (size_t i = 0; i < cl_orders.size(); i++)
		yn3.add("Order of: " + std::to_string(cl_orders.get(i)));

	y3.resize(cl_n_particles.size());
	for (size_t k = 0; k < cl_n_particles.size(); k++)
	{
		y3.get(k).resize(cl_r_cutoff.size());
		for (size_t r = 0; r < cl_r_cutoff.size(); r++)
		{
			y3.get(k).get(r).resize(n_moving+1);
			for (size_t d = 0; d < n_moving+1; d++)
			{
				y3.get(k).get(r).get(d).add(1.0);
				for (size_t m = 0; m < cl_orders.size(); m++)
				{
					if (d == 0)
						// Put a "speedup before moving the particles"
						y3.get(k).get(r).get(0).add(cl_time_rand.get(r).get(k) / cl_time_hilb.get(r).get(k).get(m));
					else
						// Put a "speedup for each moving step"
						y3.get(k).get(r).get(d).add(cl_time_rand.get(r).get(k) / cl_time_hilb_moved.get(d-1).get(r).get(k).get(m));
				}
			}
		}
	}

	// Speedup graphs report

	// Google charts options
	GCoptions options;

	options.yAxis = std::string("Speedup (times)");
	options.xAxis = std::string("Number of particles");
	options.lineWidth = 2.5;
	//options.more = "hAxis: {logScale: true}";

	GoogleChart cg;

	std::string str("<h1>Distributed " + std::to_string(dim) + "-D vector performance tests: </h1>");
	str += "<h2> 1) Speedup between an unordered positioning and an Hilbert curve positioning of particles</h2>";
	str += "We create a distributed vector (VD) of randomly positioned in a " + std::to_string(dim) + "D-box particles. Then we get a cell list of VD, with a certain cut-off radius. After that we calculate the forces of each particle. Later "
			"a VD is reordered according to an Hilbert curve, then we calculate forces again and compare the forces calculation time for an unordered and an Hilbert curve cases. The speedup is calculated and shown on graphs below, depending on different numbers of particles.";
	cg.addHTML(str);

	for (size_t i = 0; i < cl_r_cutoff.size(); i++)
	{
		options.title = std::string("Distributed vector performance speedup, cut-off radius: " + std::to_string(cl_r_cutoff.get(i)));
		cg.AddLinesGraph(x,y.get(i),yn,options);
	}

	// Calculation time graphs report

	// Google charts options
	GCoptions options2;

	options2.yAxis = std::string("Total calculation time (s)");
	options2.xAxis = std::string("Number of particles");
	options2.lineWidth = 2.5;
	//options2.more = "hAxis: {logScale: true}";

	std::string str2("<h2>2) Total calculation time</h2>");
	str2 += "We count a total calculation time. In the case of unordered positioning it is: cell list creation time + forces calculation time; in the case of an Hilbert curve positioning: reordering time + cell list creation time + forces calculation time."
			"The total calculation time is shown on graphs below, depending on different numbers of particles.";
	cg.addHTML(str2);

	for (size_t i = 0; i < cl_r_cutoff.size(); i++)
	{
		options2.title = std::string("Distributed vector performance, cut-off radius: " + std::to_string(cl_r_cutoff.get(i)));
		cg.AddLinesGraph(x,y2.get(i),yn2,options2);
	}

	// Speedup difference graphs report

	// Google charts options
	GCoptions options3;

	options3.yAxis = std::string("Speedup (times)");
	options3.xAxis = std::string("A moving step number");
	options3.lineWidth = 2.5;

	std::string str3("<h2>3) Speedup degradation after moving the particles</h2>");
	str3 += "Now we move all of the particles several times in a random direction but the same distance. Then we compare the speedup we were gaining before moving (step 0) and after moving the particles. The speedup degradation is shown below as a function of a moving step number.";

	cg.addHTML(str3);

	for (size_t k = 0; k < cl_n_particles.size(); k++)
	{
		options3.title = std::string("Distributed vector performance speedup degradation, cut-off radius: " + std::to_string(cl_r_cutoff.get(cl_r_cutoff.size()-1)) + ", number of particles: " + std::to_string(cl_n_particles.get(k)));
		cg.AddLinesGraph(x3,y3.get(k).get(cl_r_cutoff.size()-1),yn3,options3);
	}

	cg.write("Vect_dist_cl_perf_" + std::to_string(dim) + "D.html");
}

/*! \brief Function for cell list hilb performance report
 *
 */
template<unsigned int dim> void vd_celllist_performance_write_report(openfpm::vector<float> & cl_r_cutoff,openfpm::vector<size_t> & cl_n_particles,openfpm::vector<openfpm::vector<double>> & cl_time_hilb,openfpm::vector<openfpm::vector<double>> & cl_time_rand,openfpm::vector<openfpm::vector<double>> & cl_time_total_hilb,openfpm::vector<openfpm::vector<double>> & cl_time_total_rand)
{
	// Speedup graphs data

	openfpm::vector<size_t> x;
	openfpm::vector<openfpm::vector<openfpm::vector<double>>> y;
	openfpm::vector<std::string> yn;

	for (size_t i = 0; i < cl_n_particles.size() ; i++)
		x.add(cl_n_particles.get(i));

	yn.add("Random cell list vs Hilbert cell list");

	y.resize(cl_time_hilb.size());
	for (size_t r = 0; r < cl_time_hilb.size(); r++)
	{
		y.get(r).resize(cl_time_hilb.get(r).size());
		for (size_t k = 0; k < cl_time_hilb.get(r).size(); k++)
		{
			// Put a speedup
			y.get(r).get(k).add(cl_time_rand.get(r).get(k) / cl_time_hilb.get(r).get(k));
		}
	}

	// Calculation time graphs data

	openfpm::vector<openfpm::vector<openfpm::vector<double>>> y2;
	openfpm::vector<std::string> yn2;

	yn2.add("Random cell list");
	yn2.add("Hilbert cell list");

	y2.resize(cl_time_total_hilb.size());
	for (size_t r = 0; r < cl_time_total_hilb.size(); r++)
	{
		y2.get(r).resize(cl_time_total_hilb.get(r).size());
		for (size_t k = 0; k < cl_time_total_hilb.get(r).size(); k++)
		{
			// Put a total time
			y2.get(r).get(k).add(cl_time_total_rand.get(r).get(k));
			y2.get(r).get(k).add(cl_time_total_hilb.get(r).get(k));
		}
	}

	// Speedup graphs report

	// Google charts options
	GCoptions options;

	options.yAxis = std::string("Speedup (times)");
	options.xAxis = std::string("Number of particles");
	options.lineWidth = 2.5;

	GoogleChart cg;

	std::string str("<h1>Distributed " + std::to_string(dim) + "-D vector performance tests: </h1>");
	str += "<h2> 1) Speedup between a random and an hilberts curve walking</h2>";
	str += "We create a distributed vector (VD) of randomly positioned in a " + std::to_string(dim) + "D-box particles. Then we get a cell list of VD, with a certain cut-off radius. After that we calculate the forces of each particle. Later "
			"we do the same for another VD, but getting an 'hilbert cell list'. Then we compare the forces calculation time for both cases. The speedup is calculated and shown on graphs below, depending on different numbers of particles.";
	cg.addHTML(str);

	for (size_t i = 0; i < cl_r_cutoff.size(); i++)
	{
		options.title = std::string("Distributed vector performance, cut-off radius: " + std::to_string(cl_r_cutoff.get(i)));
		cg.AddLinesGraph(x,y.get(i),yn,options);
	}


	// Calculation time graphs report

	// Google charts options
	GCoptions options2;

	options2.yAxis = std::string("Total calculation time (s)");
	options2.xAxis = std::string("Number of particles");
	options2.lineWidth = 2.5;

	std::string str2("<h2>2) Total calculation time</h2>");
	str2 += "We count a total calculation time. In both cases - for a 'random' and 'hilbert' cell lists - it consists of: cell list creation time + forces calculation time. "
			"The total calculation time is shown on graphs below, depending on different numbers of particles.";
	cg.addHTML(str2);

	for (size_t i = 0; i < cl_r_cutoff.size(); i++)
	{
		options2.title = std::string("Distributed vector performance, cut-off radius: " + std::to_string(cl_r_cutoff.get(i)));
		cg.AddLinesGraph(x,y2.get(i),yn2,options2);
	}

	cg.write("Vect_dist_cl_hilb_perf_" + std::to_string(dim) + "D.html");
}

#endif /* SRC_VECTOR_VECTOR_DIST_PERFORMANCE_UTIL_HPP_ */
