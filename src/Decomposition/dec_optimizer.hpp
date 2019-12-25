#ifndef DEC_OPTIMIZER_HPP
#define DEC_OPTIMIZER_HPP

#include "Grid/iterators/grid_key_dx_iterator_sub.hpp"
#include "Grid/iterators/grid_skin_iterator.hpp"

/*! \brief this class represent a wavefront of dimension dim
 *
 * \tparam dim Dimensionality of the wavefront (dimensionality of the space
 *                                       where it live so the wavefront
 *                                       is dim-1)
 *
 * Each wavefront is identified by one starting point and one stop point.
 * More or less a wavefront is just a box defined in the integer space
 *
 */
template <unsigned int dim>
class wavefront : public Box<dim,size_t>
{
public:

	//! start point is the property with id 0 (first property)
	static const int start = 0;

	//! stop point is the property with id 1 (second property)
	static const int stop = 1;
};

///// Unfortunately it seem that nvcc it specialize incorrectly this data structure so we have to specialize for the broken cases

template<unsigned int dim>
struct is_typedef_and_data_same<true,wavefront<dim>>
{
	enum
	{
		value = 1
	};
};

/*! \brief This class take a graph representing the space decomposition and produce a
 *         simplified version
 *
 * Given a Graph_CSR and a seed point produce an alternative decomposition in boxes with
 * less sub-domain. In the following we referee with sub-domain the boxes produced by this
 * algorithm and sub-sub-domain the sub-domain before reduction
 *
 */

template <unsigned int dim, typename Graph>
class dec_optimizer
{
	//! Contain information about the grid size
	grid_sm<dim,void> gh;

private:

	/*! \brief Expand one wavefront
	 *
	 * \param v_w wavefronts
	 * \param w_comb wavefront expansion combinations
	 * \param d direction of expansion
	 *
	 */
	void expand_one_wf(openfpm::vector<wavefront<dim>> & v_w, std::vector<comb<dim>> & w_comb , size_t d)
	{
		for (size_t j = 0 ; j < dim ; j++)
		{
			v_w.template get<wavefront<dim>::stop>(d)[j] = v_w.template get<wavefront<dim>::stop>(d)[j] + w_comb[d].c[j];
			v_w.template get<wavefront<dim>::start>(d)[j] = v_w.template get<wavefront<dim>::start>(d)[j] + w_comb[d].c[j];
		}
	}


	/*! \brief Adjust the other wavefronts
	 *
	 * \param v_w array of wavefronts
	 * \param hyp Hyper cube used to adjust the wavefront
	 * \param w_comb for each wavefront indicate their position (normal to the face of the wavefront)
	 * \param d direction
	 *
	 */
	void adjust_others_wf(openfpm::vector<wavefront<dim>> & v_w,  HyperCube<dim> & hyp, std::vector<comb<dim>> & w_comb, size_t d)
	{
		// expand the intersection of the wavefronts

		std::vector<comb<dim>> q_comb = SubHyperCube<dim,dim-1>::getCombinations_R(w_comb[d],(int)dim-2);

		// Eliminate the w_comb[d] direction

		for (size_t k = 0 ; k < q_comb.size() ; k++)
		{
			for (size_t j = 0 ; j < dim ; j++)
			{
				if (w_comb[d].c[j] != 0)
				{
					q_comb[k].c[j] = 0;
				}
			}
		}

		// for all the combinations
		for (size_t j = 0 ; j < q_comb.size() ; j++)
		{
			size_t id = hyp.LinId(q_comb[j]);

			// get the combination of the direction d
			bool is_pos = hyp.isPositive(d);

			// is positive, modify the stop point or the starting point

			for (size_t s = 0 ; s < dim ; s++)
			{
				if (is_pos == true)
				{v_w.template get<wavefront<dim>::stop>(id)[s] = v_w.template get<wavefront<dim>::stop>(id)[s] + w_comb[d].c[s];}
				else
				{v_w.template get<wavefront<dim>::start>(id)[s] = v_w.template get<wavefront<dim>::start>(id)[s] + w_comb[d].c[s];}
			}
		}
	}

	/*! \brief Fill the wavefront position
	 *
	 * \tparam prp property to set
	 *
	 * \param graph we are processing
	 * \param v_w array of wavefronts
	 *
	 */
	template<unsigned int prp> void write_wavefront(Graph & graph,openfpm::vector<wavefront<dim>> & v_w)
	{
		// fill the wall domain with 0

		fill_domain<prp>(graph,gh.getBox(),0);

		// fill the wavefront

		for (int i = 0 ; i < v_w.size() ; i++)
		{
			Box<dim,size_t> box = wavefront<dim>::getBox(v_w.get(i));

			fill_domain<prp>(graph,box,1);
		}
	}

	/*! \brief Fill the domain
	 *
	 * \tparam p_sub property to set with the sub-domain id
	 *
	 * \param graph we are processing
	 * \param box Box to fill
	 * \param ids value to fill with
	 *
	 */

	template<unsigned int p_sub> void fill_domain(Graph & graph,const Box<dim,size_t> & box, long int ids)
	{
		// Create a subgrid iterator
		grid_key_dx_iterator_sub<dim,no_stencil,do_not_print_warning_on_adjustment<dim>> g_sub(gh,box.getKP1(),box.getKP2());

		// iterate through all grid points

		while (g_sub.isNext())
		{
			// get the actual key
			const grid_key_dx<dim> & gk = g_sub.get();

			// get the vertex and set the sub id

			graph.vertex(gh.LinId(gk)).template get<p_sub>() = ids;

			// next subdomain
			++g_sub;
		}
	}

	/*! \brief Add the boundary domain of id p_id to the queue
	 *
	 * \tparam p_sub property id where to store the sub-domain decomposition
	 * \tparam p_id property id where is stored the decomposition
	 *
	 * \param domains vector with sub-sub-domains still to process
	 * \param v_w array of wave-fronts
	 * \param graph we are processing
	 * \param w_comb wavefront combination, it is the normal vector to the wavefront
	 * \param pr_id processor id for which we are optimizing the decomposition
	 * \param bc boundary conditions
	 *
	 */
	template<unsigned int p_sub, unsigned int p_id> void add_to_queue(openfpm::vector<size_t> & domains, openfpm::vector<wavefront<dim>> & v_w, Graph & graph,  std::vector<comb<dim>> & w_comb, long int pr_id, const size_t(& bc)[dim])
	{
		// create a new queue
		openfpm::vector<size_t> domains_new;

		// push in the new queue, the old domains of the queue that are not assigned element

		for (size_t j = 0 ; j < domains.size() ; j++)
		{
			long int gs = graph.vertex(domains.get(j)).template get<p_sub>();
			if (gs < 0)
			{
				// not assigned push it

				domains_new.add(domains.get(j));
			}
		}

		// Create an Hyper-cube
		HyperCube<dim> hyp;

		for (size_t d = 0 ; d < v_w.size() ; d++)
		{
			expand_one_wf(v_w,w_comb,d);
			adjust_others_wf(v_w,hyp,w_comb,d);
		}

		// for each expanded wavefront create a sub-grid iterator and add the sub-domain

		for (size_t d = 0 ; d < v_w.size() ; d++)
		{
			// Create a sub-grid iterator
			grid_key_dx_iterator_sub_bc<dim,no_stencil,do_not_print_warning_on_adjustment<dim>> g_sub(gh,v_w.template get<wavefront<dim>::start>(d),v_w.template get<wavefront<dim>::stop>(d),bc);

			// iterate through all grid points

			while (g_sub.isNext())
			{
				// get the actual key
				const grid_key_dx<dim> & gk = g_sub.get();

				// get the vertex and if does not have a sub-id and is assigned ...
				long int pid = graph.vertex(gh.LinId(gk)).template get<p_sub>();

				// Get the processor id of the sub-sub-domain
				long int pp_id = graph.vertex(gh.LinId(gk)).template get<p_id>();

				// if the sub-sub-domain is not assigned
				if (pid < 0)
				{
					// ... and we are not processing the full graph
					if (pr_id != -1)
					{
						// ... and the processor id of the sub-sub-domain match the part we are processing, add to the queue

						if ( pr_id == pp_id)
							domains_new.add(gh.LinId(gk));
					}
					else
						domains_new.add(gh.LinId(gk));
				}

				++g_sub;
			}
		}

		// copy the new queue to the old one (it not copied, C++11 move semantic)
		domains.swap(domains_new);
	}

	/*! \brief Find the biggest hyper-cube
	 *
	 * starting from one initial sub-domain find the biggest hyper-cube
	 * output the box, and fill a list of neighborhood processor
	 *
	 * \tparam p_sub id of the property storing the sub-decomposition
	 * \tparam p_id id of the property containing the decomposition
	 *
	 * \param start_p initial domain
	 * \param graph representing the grid of sub-sub-domain
	 * \param box produced box
	 * \param v_w Wavefronts
	 * \param w_comb wavefronts directions (0,0,1) (0,0,-1) (0,1,0) (0,-1,0) ...
	 *
	 */
	template <unsigned int p_sub, unsigned int p_id> void expand_from_point(size_t start_p, Graph & graph, Box<dim,size_t> & box, openfpm::vector<wavefront<dim>> & v_w , std::vector<comb<dim>> & w_comb)
	{
		// We assume that Graph is the rapresentation of a cartesian graph
		// this mean that the direction d is at the child d

		// Get the number of wavefronts
		size_t n_wf = w_comb.size();

		// Create an Hyper-cube
		HyperCube<dim> hyp;

		// direction of expansion

		size_t domain_id = graph.vertex(start_p).template get<p_id>();
		bool can_expand = true;

		// while is possible to expand

		while (can_expand)
		{
			// reset can expand
			can_expand = false;

			// for each direction of expansion expand the wavefront

			for (size_t d = 0 ; d < n_wf ; d++)
			{
				// number of processed sub-domain
				size_t n_proc_sub = 0;

				// flag to indicate if the wavefront can expand
				bool w_can_expand = true;

				// Create an iterator of the expanded wavefront
				grid_key_dx<dim> start = grid_key_dx<dim>(v_w.template get<wavefront<dim>::start>(d)) + w_comb[d];
				grid_key_dx<dim> stop = grid_key_dx<dim>(v_w.template get<wavefront<dim>::stop>(d)) + w_comb[d];
				grid_key_dx_iterator_sub<dim,no_stencil,do_not_print_warning_on_adjustment<dim>> it(gh,start,stop);

				// for each sub-domain in the expanded wavefront
				while (it.isNext())
				{
					// get the wavefront sub-domain id
					size_t sub_w_e = gh.LinId(it.get());

					// we get the processor id of the neighborhood sub-domain on direction d
					// (expanded wavefront)
					size_t exp_p = graph.vertex(sub_w_e).template get<p_id>();

					// Check if already assigned
					long int ass = graph.vertex(sub_w_e).template get<p_sub>();

					// we check if it is the same processor id and is not assigned
					w_can_expand &= ((exp_p == domain_id) & (ass < 0));

					// next domain
					++it;

					// increase the number of processed sub-domain
					n_proc_sub++;
				}

				// if we did not processed sub-domain, we cannot expand
				w_can_expand &= (n_proc_sub != 0);

				// if you can expand one wavefront we did not reach the end
				can_expand |= w_can_expand;

				// if we can expand the wavefront expand it
				if (w_can_expand == true)
				{
					// expand the wavefront
					for (size_t j = 0 ; j < dim ; j++)
					{
						v_w.template get<wavefront<dim>::stop>(d)[j] = v_w.template get<wavefront<dim>::stop>(d)[j] + w_comb[d].c[j];
						v_w.template get<wavefront<dim>::start>(d)[j] = v_w.template get<wavefront<dim>::start>(d)[j] + w_comb[d].c[j];
					}

					// expand the intersection of the wavefronts

					if (dim >= 2)
					{
						std::vector<comb<dim>> q_comb = SubHyperCube<dim,dim-1>::getCombinations_R(w_comb[d],(int)dim-2);

						// Eliminate the w_comb[d] direction

						for (size_t k = 0 ; k < q_comb.size() ; k++)
						{
							for (size_t j = 0 ; j < dim ; j++)
							{
								if (w_comb[d].c[j] != 0)
								{
									q_comb[k].c[j] = 0;
								}
							}
						}

						// for all the combinations
						for (size_t j = 0 ; j < q_comb.size() ; j++)
						{
							size_t id = hyp.LinId(q_comb[j]);

							// get the combination of the direction d

							bool is_pos = hyp.isPositive(d);

							// is positive, modify the stop point or the starting point

							for (size_t s = 0 ; s < dim ; s++)
							{
								if (is_pos == true)
								{v_w.template get<wavefront<dim>::stop>(id)[s] = v_w.template get<wavefront<dim>::stop>(id)[s] + w_comb[d].c[s];}
								else
								{v_w.template get<wavefront<dim>::start>(id)[s] = v_w.template get<wavefront<dim>::start>(id)[s] + w_comb[d].c[s];}
							}
						}
					}
				}
			}
		}

		// get back the hyper-cube produced

		for (size_t i = 0 ; i < dim ; i++)
		{
			// get the index of the wavefront direction
			size_t p_f = hyp.positiveFace(i);
			size_t n_f = hyp.negativeFace(i);

			// set the box
			box.setHigh(i,v_w.template get<wavefront<dim>::stop>(p_f)[i]);
			box.setLow(i,v_w.template get<wavefront<dim>::start>(n_f)[i]);
		}
	}

	/*! \brief Initialize the wavefronts
	 *
	 * \param start_p starting point for the wavefront set
	 * \param v_w Wavefront array
	 *
	 */
	void InitializeWavefront(grid_key_dx<dim> & start_p, openfpm::vector<wavefront<dim>> & v_w)
	{
		// Wavefront to initialize

		for (size_t i = 0 ; i < v_w.size() ; i++)
		{
			for (size_t j = 0 ; j < dim ; j++)
			{
				v_w.template get<wavefront<dim>::start>(i)[j] = start_p.get(j);
				v_w.template get<wavefront<dim>::stop>(i)[j] = start_p.get(j);
			}
		}
	}

	/*! \brief Get the first seed
	 *
	 * search in the graph for one sub-domain labelled with processor id
	 * to use as seed
	 *
	 * \tparam p_id property id containing the decomposition
	 * \tparam p_sub property id that will contain the sub-domain decomposition
	 *
	 * \param graph Graph
	 * \param id processor id
	 *
	 * \return a valid seed key
	 *
	 */
	template<unsigned int p_id, unsigned int p_sub> grid_key_dx<dim> search_seed(Graph & graph, long int id)
	{
		// if no processor is selected return the first point
		if (id < -1)
		{
			grid_key_dx<dim> key;
			key.zero();

			return key;
		}

		// Create a grid iterator
		grid_key_dx_iterator<dim> g_sub(gh);

		// iterate through all grid points

		while (g_sub.isNext())
		{
			// get the actual key
			const grid_key_dx<dim> & gk = g_sub.get();

			// if the subdomain has the id we are searching stop
			if ((long int)graph.vertex(gh.LinId(gk)).template get<p_id>() == id && graph.vertex(gh.LinId(gk)).template get<p_sub>() == -1)
			{
				return gk;
			}

			++g_sub;
		}

		// If not found return an invalid key
		grid_key_dx<dim> key;
		key.invalid();

		return key;
	}


	/*! \brief optimize the graph
	 *
	 * Starting from a domain (hyper-cubic), it create wavefront at the boundary and expand
	 * the boundary until the wavefronts cannot expand any more.
	 * To the domains inside the hyper-cube one sub-id is assigned. This procedure continue until
	 * all the domain of one p_id has a sub-id
	 *
	 * \tparam p_id property containing the decomposition
	 * \tparam p_sub property to fill with the sub-domain decomposition
	 *
	 * \param start_p seed point
	 * \param graph we are processing
	 * \param pr_id Processor id (if p_id == -1 the optimization is done for all the processors)
	 * \param lb list of sub-domain boxes produced by the algorithm
	 * \param box_nn_processor for each sub-domain it list all the neighborhood processors
	 * \param ghe Ghost extension in sub-sub-domain units in each direction
	 * \param init_sub_id when true p_sub property is initially set to -1 [default true]
	 * \param sub_id starting sub_id to enumerate them [default 0]
	 * \param bc boundary conditions
	 *
	 * \return last assigned sub-id
	 *
	 */
	template <unsigned int p_sub, unsigned int p_id> size_t optimize(grid_key_dx<dim> & start_p, Graph & graph, long int pr_id, openfpm::vector<Box<dim,size_t>> & lb, openfpm::vector< openfpm::vector<size_t> > & box_nn_processor , const Ghost<dim,long int> & ghe ,const size_t (& bc)[dim], bool init_sub_id = true, size_t sub_id = 0)
	{
		// queue
		openfpm::vector<size_t> v_q;

		// box list 2
		openfpm::vector< openfpm::vector<size_t> > box_nn_processor2;

		// Create an hyper-cube
		HyperCube<dim> hyp;

		// Get the wavefront combinations
		std::vector<comb<dim>> w_comb = hyp.getCombinations_R(dim-1);

		// wavefronts
		openfpm::vector<wavefront<dim>> v_w(w_comb.size());

		// fill the sub decomposition with negative number

		if (init_sub_id == true)
			fill_domain<p_sub>(graph,gh.getBox(),-1);

		// push the first domain
		v_q.add(gh.LinId(start_p));

		while (v_q.size() != 0)
		{
			// Box
			Box<dim,size_t> box;

			// Get the grid_key position from the linearized id
			start_p = gh.InvLinId(v_q.get(0));

			// Initialize the wavefronts from the domain start_p
			InitializeWavefront(start_p,v_w);

			// Create the biggest box containing the domain
			expand_from_point<p_sub,p_id>(v_q.get(0),graph,box,v_w,w_comb);

			// Add the created box to the list of boxes
			lb.add(box);

			// fill the domain
			fill_domain<p_sub>(graph,box,sub_id);

			// add the surrounding sub-domain to the queue
			add_to_queue<p_sub,p_id>(v_q,v_w,graph,w_comb,pr_id,bc);

			// increment the sub_id
			sub_id++;
		}

		return sub_id;
	}

	/*! \brief Construct the sub-domain processor list
	 *
	 * \tparam p_id property that contain the decomposition
	 *
	 * Each entry is a sub-domain, the list of numbers indicate the neighborhood processors
	 *
	 * \param graph graph to process
	 * \param box_nn_processor for each sub-domain it list all the neighborhood processors
	 * \param subs vector of sub-domains
	 * \param ghe ghost extensions
	 * \param bc boundary conditions
	 * \param pr_id processor that we are processing
	 *
	 */
	template<unsigned int p_id> void construct_box_nn_processor(Graph & graph, openfpm::vector< openfpm::vector<size_t> > & box_nn_processor, const openfpm::vector<Box<dim,size_t>> & subs, const Ghost<dim,long int> & ghe, const size_t (& bc)[dim], long int pr_id)
	{
		std::unordered_map<size_t,size_t> map;

		for (size_t i = 0 ; i < subs.size() ; i++)
		{
			map.clear();
			Box<dim,size_t> sub = subs.get(i);
			sub.enlarge(ghe);

			grid_skin_iterator_bc<dim> gsi(gh,subs.get(i),sub,bc);

			while (gsi.isNext())
			{
				auto key = gsi.get();

				size_t pp_id = graph.vertex(gh.LinId(key)).template get<p_id>();
				if (pr_id != (long int)pp_id)
					map[pp_id] = pp_id;

				++gsi;
			}

			// Add the keys to box_nn_processors

			box_nn_processor.add();
			for ( auto it = map.begin(); it != map.end(); ++it )
			{
				box_nn_processor.last().add(it->first);
			}
		}
	}

public:

	/*! \brief Constructor
	 *
	 * \param g Graph to simplify
	 * \param sz size of the grid on each dimension
	 *
	 */

	dec_optimizer(Graph & g, const size_t (& sz)[dim])
	:gh(sz)
	{
		// The graph g is suppose to represent a cartesian grid
		// No check is performed on g
	}

	/*! \brief optimize the graph
	 *
	 * Starting from a sub-sub-domain, it create wavefronts at the boundary and expand
	 * the boundary until the wavefronts cannot expand any more, creating a sub-domain covering more sub-sub-domain.
	 * This procedure continue until all the domain is covered by a sub-domains
	 *
	 * \tparam p_id property containing the processor decomposition
	 * \tparam p_sub property to fill with the sub-domain decomposition
	 *
	 * \param start_p seed point
	 * \param graph we are processing
	 * \param ghe ghost size
	 * \param bc boundary conditions
	 *
	 */
	template <unsigned int p_sub, unsigned int p_id> void optimize(grid_key_dx<dim> & start_p, Graph & graph, const Ghost<dim,long int> & ghe , const size_t (& bc)[dim])
	{
		// temporal vector
		openfpm::vector<Box<dim,size_t>> tmp;

		// temporal vector
		openfpm::vector< openfpm::vector<size_t> > box_nn_processor;

		// optimize
		optimize<p_sub,p_id>(start_p,graph,-1,tmp, box_nn_processor,ghe,bc);
	}

	/*! \brief optimize the graph
	 *
	 * Starting from a sub-sub-domain, it create wavefronts at the boundary and expand
	 * the boundary until the wavefronts cannot expand any more, creating a sub-domain covering more sub-sub-domain.
	 * This procedure continue until all the sub-domain of the processor p_id are covered by a sub-domains
	 *
	 * \tparam p_id property containing the decomposition
	 * \tparam p_sub property to fill with the sub-domain decomposition
	 *
	 * \param graph we are processing
	 * \param pr_id Processor id (if p_id == -1 the optimization is done for all the processors)
	 * \param lb list of sub-domain boxes
	 * \param box_nn_processor for each sub-domain it list all the neighborhood processors
	 * \param ghe ghost size
	 *
	 */
	template <unsigned int p_sub, unsigned int p_id> void optimize(Graph & graph, long int pr_id, openfpm::vector<Box<dim,size_t>> & lb, openfpm::vector< openfpm::vector<size_t> > & box_nn_processor, const Ghost<dim,long int> & ghe, const size_t (& bc)[dim])
	{
		grid_key_dx<dim> key_seed;
		key_seed.zero();

		// if processor is -1 call optimize with -1 to do on all processors and exit
		if (pr_id == -1)
		{
			optimize<p_sub,p_id>(key_seed,graph,pr_id,lb,box_nn_processor,ghe,bc);

			// Construct box box_nn_processor from the constructed domain
			construct_box_nn_processor<p_id>(graph,box_nn_processor,lb,ghe,bc,pr_id);

			return;
		}

		size_t sub_id = 0;

		// fill the sub decomposition with negative number
		fill_domain<p_sub>(graph,gh.getBox(),-1);

		key_seed = search_seed<p_id,p_sub>(graph,pr_id);

		while (key_seed.isValid())
		{
			// optimize
			sub_id = optimize<p_sub,p_id>(key_seed,graph,pr_id,lb,box_nn_processor,ghe,bc,false,sub_id);

			// new seed
			key_seed = search_seed<p_id,p_sub>(graph,pr_id);
		}

		// Construct box box_nn_processor from the constructed domain
		construct_box_nn_processor<p_id>(graph,box_nn_processor,lb,ghe,bc,pr_id);
	}
};

#endif
