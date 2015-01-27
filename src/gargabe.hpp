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

#endif /* GARGABE_HPP_ */
