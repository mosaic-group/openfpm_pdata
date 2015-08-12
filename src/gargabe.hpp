/*
 * gargabe.hpp
 *
 *  Created on: Jan 13, 2015
 *      Author: i-bird
 */

#ifndef GARGABE_HPP_
#define GARGABE_HPP_



	template <unsigned int j, unsigned int i, typename Graph> void optimize(size_t start_p, Graph & graph)
	{
		// We assume that Graph is the rapresentation of a cartesian graph
		// this mean that the direction d is at the child d

		// Create an Hyper-cube

		HyperCube<dim> hyp;

		// Get the number of wavefronts

		size_t n_wf = hyp.getNumberOfElements_R(0);

		// Get the number of intersecting wavefront



		// Get the number of sub-dimensional common wavefront
		// basically are a list of all the subdomain common to two or more

		// Create n_wf wavefront queue

		openfpm::vector<wavefront> v_w;
		v.reserve(n_wf);

		// direction of expansion

		size_t domain_id = 0;
		int exp_dir = 0;
		bool can_expand = true;

		// while is possible to expand

		while (can_expand)
		{
			// for each direction of expansion expand the wavefront

			for (int d = 0 ; d < n_wf ; d++)
			{
				// get the wavefront at direction d

				openfpm::vector<size_t> & wf_d = v_w.get<wavefront::domains>(d);

				// flag to indicate if the wavefront can expand

				bool w_can_expand = true;

				// for each subdomain

				for (size_t sub = 0 ; sub < wf_d.size() ; sub++)
				{
					// check if the adjacent domain in direction d exist
					// and is of the same id

					// get the starting subdomain
					size_t sub_w = wf_d.get<0>(sub);

					// we get the processor id of the neighborhood sub-domain on direction d
					size_t exp_p = graph.getChild(sub_w,d).get<j>();

					// we check if it is the same processor id
					if (exp_p != domain_id)
					{
						w_can_expand = false;
					}
				}

				// if we can expand the wavefront expand it
				if (w_can_expand == true)
				{
					// for each subdomain
					for (size_t sub = 0 ; sub < wf_d.size() ; sub++)
					{
						// update the position of the wavefront
						wf_d.get<0>(sub) = wf_d.get<0>(sub) + gh.stride(d);
					}

					// here we add sub-domains to all the other queues
					// get the face of the hyper-cube

					SubHyperCube<dim,dim-1> sub_hyp = hyp.getSubHyperCube(d);

					std::vector<comb<dim>> q_comb = sub_hyp.getCombinations_R(dim-2);
				}
			}
		}

		// For each point in the Hyper-cube check if we can move the wave front


	}

#ifndef PARALLEL_DECOMPOSITION
//		CreateSubspaces();
#endif

#ifndef USE_METIS_GP

		// Here we do not use METIS
		// Distribute the divided domains

		// Get the number of processing units
		size_t Np = v_cl.getProcessingUnits();

		// Get the ID of this processing unit
		// and push the subspace is taking this
		// processing unit

		for (size_t p_id = v_cl.getProcessUnitID(); p_id < Np ; p_id += Np)
			id_sub.push_back(p_id);
#else


#endif



		/////////////// DEBUG /////////////////////

		// get the decomposition
		auto & dec = g_dist.getDecomposition();

		Vcluster & v_cl = *global_v_cluster;

		// check the consistency of the decomposition
		val = dec.check_consistency();
		BOOST_REQUIRE_EQUAL(val,true);

		// for each local volume
		// Get the number of local grid needed
		size_t n_grid = dec.getNLocalHyperCube();

		size_t vol = 0;

		openfpm::vector<Box<2,size_t>> v_b;

		// Allocate the grids
		for (size_t i = 0 ; i < n_grid ; i++)
		{
			// Get the local hyper-cube
			SpaceBox<2,float> sub = dec.getLocalHyperCube(i);

			Box<2,size_t> g_box = g_dist.getCellDecomposer().convertDomainSpaceIntoGridUnits(sub);
			v_b.add(g_box);

			vol += g_box.getVolumeKey();
		}

		v_cl.reduce(vol);
		v_cl.execute();

		BOOST_REQUIRE_EQUAL(vol,k*k);

		/////////////////////////////////////


		// 3D test

	//	g_dist.write("");

	/*	auto g_it = g_dist.getIteratorBulk();

		auto g_it_halo = g_dist.getHalo();

		// Let try to solve the poisson equation d2(u) = f with f = 1 and computation
		// comunication overlap (100 Jacobi iteration)

		for (int i = 0 ; i < 100 ; i++)
		{
			g_dist.ghost_get();

			// Compute the bulk

			jacobi_iteration(g_it);

			g_dist.ghost_sync();

			// Compute the halo

			jacobi_iteration(g_it_halo);
		}*/


		BOOST_AUTO_TEST_CASE( grid_dist_id_poisson_test_use)
		{
			// grid size
		/*	size_t sz[2] = {1024,1024};

			// Distributed grid with id decomposition

			grid_dist_id<2, scalar<float>, CartDecomposition<2,size_t>> g_dist(sz);

			// Create the grid on memory

			g_dist.Create();*/

		/*	auto g_it = g_dist.getIteratorBulk();

			auto g_it_halo = g_dist.getHalo();

			// Let try to solve the poisson equation d2(u) = f with f = 1 and computation
			// comunication overlap (100 Jacobi iteration)

			for (int i = 0 ; i < 100 ; i++)
			{
				g_dist.ghost_get();

				// Compute the bulk

				jacobi_iteration(g_it);

				g_dist.ghost_sync();

				// Compute the halo

				jacobi_iteration(g_it_halo);
			}*/
		}

		template<typename iterator> void jacobi_iteration(iterator g_it, grid_dist_id<2, float, scalar<float>, CartDecomposition<2,float>> & g_dist)
		{
			// scalar
			typedef scalar<float> S;

			// iterator

			while(g_it.isNext())
			{
				// Jacobi update

				auto pos = g_it.get();

				g_dist.template get<S::ele>(pos) = (g_dist.template get<S::ele>(pos.move(0,1)) +
			                             g_dist.template get<S::ele>(pos.move(0,-1)) +
			                             g_dist.template get<S::ele>(pos.move(1,1)) +
			                             g_dist.template get<S::ele>(pos.move(1,-1)) / 4.0);

				++g_it;
			}
		}

#endif /* GARGABE_HPP_ */
