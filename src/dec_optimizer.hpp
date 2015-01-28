#ifndef DEC_OPTIMIZER_HPP
#define DEC_OPTIMIZER_HPP

/*! \brief this class represent a wavefront of dimension dim
 *
 * \dim Dimensionality of the wavefront (dimensionality of the space
 *                                       where it live so the wavefront
 *                                       is dim-1)
 *
 */

template <unsigned int dim>
class wavefront
{
public:

	typedef boost::fusion::vector<size_t[dim],size_t[dim]> type;

	typedef typename memory_traits_inte<type>::type memory_int;
	typedef typename memory_traits_lin<type>::type memory_lin;

	type data;

	static const int start = 0;
	static const int stop = 1;
	static const int max_prop = 2;

	/* \brief Get the key to the point 1
	 *
	 * \return the key to the point 1
	 *
	 */

	grid_key_dx<dim> getKP1()
	{
		// grid key to return
		grid_key_dx<dim> ret(boost::fusion::at_c<start>(data));

		return ret;
	}

	/* \brief Get the key to point 2
	 *
	 * \return the key to the point 2
	 *
	 */

	grid_key_dx<dim> getKP2()
	{
		// grid key to return
		grid_key_dx<dim> ret(boost::fusion::at_c<stop>(data));

		return ret;
	}

	/* \brief produce a box from an encapsulated object
	 *
	 * \param encap encapsulated object
	 *
	 */

	template<typename encap> static Box<dim,size_t> getBox(encap && enc)
	{
		Box<dim,size_t> bx;

		// Create the object from the encapsulation

		for (int i = 0 ; i < dim ; i++)
		{
			bx.setHigh(i,enc.template get<wavefront::start>()[i]);
			bx.setLow(i,enc.template get<wavefront::stop>()[i]);
		}

		return bx;
	}

	/* \brief produce a box from an encapsulated object
	 *
	 * \param encap encapsulated object
	 *
	 */

	template<typename encap> static void getBox(const encap & enc, Box<dim,size_t> & bx)
	{
		// Create the object from the encapsulation

		for (int i = 0 ; i < dim ; i++)
		{
			bx.setHigh(i,enc.template get<wavefront::start>()[i]);
			bx.setLow(i,enc.template get<wavefront::stop>()[i]);
		}
	}
};

/*! \brief This class take a graph representing the space decomposition and produce a
 *         simplified version
 *
 * Given a Graph_CSR and a seed point produce an alternative decomposition to reduce
 * the ghost over-stress
 *
 */

template <unsigned int dim, typename Graph>
class dec_optimizer
{
	// create a grid header for helping

	grid<dim,void> gh;

private:

	/* \brief Fill the wavefront position
	 *
	 * \tparam prp property to set
	 *
	 * \param graph we are processing
	 * \param Box to fill
	 * \param id value to fill with
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

	/* \brief Fill the domain
	 *
	 * \tparam p_sub property to set with the sub-domain id
	 *
	 * \param graph we are processing
	 * \param Box to fill
	 * \param ids value to fill with
	 *
	 */

	template<unsigned int p_sub> void fill_domain(Graph & graph,const Box<dim,size_t> & box, size_t ids)
	{
		// Create a subgrid iterator
		grid_key_dx_iterator_sub<dim> g_sub(gh,box.getKP1(),box.getKP2());

		// iterate through all grid points

		while (g_sub.isNext())
		{
			// get the actual key
			const grid_key_dx<dim> & gk = g_sub.get();

			// get the vertex and set the sub id

			graph.vertex(gh.LinId(gk)).template get<p_sub>() = ids;
		}
	}

	/* \brief Add the boundary domain to the queue
	 *
	 * \tparam i-property where is stored the decomposition
	 *
	 * \param domains vector with domains to process
	 * \param graph we are processing
	 * \param w_comb hyper-cube combinations
	 *
	 */
	template<unsigned int p_sub> void add_to_queue(openfpm::vector<size_t> & domains, openfpm::vector<wavefront<dim>> & v_w, Graph & graph,  std::vector<comb<dim>> & w_comb)
	{
		// create a new queue
		openfpm::vector<size_t> domains_new;

		// push in the new queue, the old domains of the queue that are not assigned element

		for (size_t j = 0 ; j < domains.size() ; j++)
		{
			if (graph.vertex(domains.get(j)).template get<p_sub>() < 0)
			{
				// not assigned push it

				domains_new.add(domains.get(j));
			}
		}

		// take the wavefront expand on direction d of one

		for (int d = 0 ; d < v_w.size() ; d++)
		{
			// expand the wavefront
			for (int j = 0 ; j < dim ; j++)
			{
				v_w.template get<wavefront<dim>::start>(d)[j] = v_w.template get<wavefront<dim>::start>(d)[j] + w_comb[d].c[j];
				v_w.template get<wavefront<dim>::stop>(d)[j] = v_w.template get<wavefront<dim>::stop>(d)[j] + w_comb[d].c[j];
			}
		}

		// for each expanded wavefront create a sub-grid iterator and add the subbomain

		for (int d = 0 ; d < v_w.size() ; d++)
		{
			// Create a sub-grid iterator
			grid_key_dx_iterator_sub<dim> g_sub(gh,v_w.template get<wavefront<dim>::start>(d),v_w.template get<wavefront<dim>::stop>(d));

			// iterate through all grid points

			while (g_sub.isNext())
			{
				// get the actual key
				const grid_key_dx<dim> & gk = g_sub.get();

				// get the vertex and check if is not assigned add to the queue

				if (graph.vertex(gh.LinId(gk)).template get<p_sub>() < 0)
				{
					// add

					domains_new.add(gh.LinId(gk));
				}
			}
		}

		// copy the new queue to the old one (it not copied, C++11 move semantic)
		domains = domains_new;
	}

	/* \brief Find the biggest hyper-cube
	 *
	 * starting from one initial sub-domain find the biggest hyper-cube
	 *
	 * \tparam j id of the property storing the sub-decomposition
	 * \tparam i id of the property containing the decomposition
	 *
	 * \param start_p initial domain
	 * \param graph representing the grid
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

		size_t domain_id = graph.vertex(gh.LinId(start_p)).template get<p_id>();
		bool can_expand = true;

		// while is possible to expand

		while (can_expand)
		{
			// reset can expand
			can_expand = false;

			// for each direction of expansion expand the wavefront

			for (int d = 0 ; d < n_wf ; d++)
			{
				// number of processed sub-domain
				size_t n_proc_sub = 0;

				// flag to indicate if the wavefront can expand
				bool w_can_expand = true;

				// Create an iterator of the wavefront
				grid_key_dx_iterator_sub<dim> it(gh,v_w.template get<wavefront<dim>::start>(d),v_w.template get<wavefront<dim>::stop>(d));

				// for each subdomain
				while (it.isNext())
				{
					// get the wavefront sub-domain id
					size_t sub_w = gh.LinId(it.get());

					// get the sub-domain id of the expanded wavefront
					grid_key_dx<dim> k_exp = it.get() + w_comb[d];
					size_t sub_w_e = gh.LinId(k_exp);

					// we get the processor id of the neighborhood sub-domain on direction d
					size_t exp_p = graph.vertex(sub_w_e).template get<p_id>();

					// we check if it is the same processor id
					w_can_expand &= exp_p == domain_id;

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
					for (int j = 0 ; j < dim ; j++)
					{
						v_w.template get<wavefront<dim>::stop>(d)[j] = v_w.template get<wavefront<dim>::stop>(d)[j] + w_comb[d].c[j];
						v_w.template get<wavefront<dim>::start>(d)[j] = v_w.template get<wavefront<dim>::start>(d)[j] + w_comb[d].c[j];
					}

					// expand the intersection of the wavefronts

					std::vector<comb<dim>> q_comb = SubHyperCube<dim,dim-1>::getCombinations_R(w_comb[d],dim-2);

					// Eliminate the w_comb[d] direction

					for (int k = 0 ; k < q_comb.size() ; k++)
					{
						for (int j = 0 ; j < dim ; j++)
						{
							if (w_comb[d].c[j] != 0)
							{
								q_comb[k].c[j] = 0;
							}
						}
					}

					// for all the combinations
					for (int j = 0 ; j < q_comb.size() ; j++)
					{
						size_t id = hyp.LinId(q_comb[j]);

						// get the combination of the direction d

						bool is_pos = hyp.isPositive(d);

						// is positive, modify the stop point or the starting point

						for (int s = 0 ; s < dim ; s++)
						{
							if (is_pos == true)
							{v_w.template get<wavefront<dim>::stop>(id)[j] = v_w.template get<wavefront<dim>::stop>(id)[j] + w_comb[id].c[j];}
							else
							{v_w.template get<wavefront<dim>::start>(id)[j] = v_w.template get<wavefront<dim>::start>(id)[j] + w_comb[id].c[j];}
						}
					}
				}
			}

			// Debug output the wavefront graph for debug

			// duplicate the graph

			Graph g_debug = graph.duplicate();
			write_wavefront<nm_part_v::id>(g_debug,v_w);

			VTKWriter<Graph> vtk(g_debug);
			vtk.template write<1>("vtk_debug.vtk");

			////////////////////////////////////
		}

		// get back the hyper-cube produced

		for (int i = 0 ; i < dim ; i++)
		{
			// get the index of the wavefront direction
			int w = 0;

			for (int j = 0 ; j < dim ; j++)
			{
				if (w_comb[i].c[j] == 1)
				{
					w = j;
					break;
				}
			}

			// set the box
			box.setHigh(i,v_w.template get<wavefront<dim>::stop>(i)[w]);
			box.setLow(i,v_w.template get<wavefront<dim>::start>(i)[w]);
		}
	}

public:

	/*! \brief Constructor
	 *
	 * \param g Graph to simplify
	 * \param sz size of the grid on each dimension
	 *
	 */

	dec_optimizer(Graph & g, std::vector<size_t> sz)
	:gh(sz)
	{
		// The graph g is suppose to represent a cartesian grid
		// No check is performed on g
	}

	/*! \brief optimize the graph
	 *
	 * Starting from a domain (hyper-cubic), it create wavefront at the boundary and expand
	 * the boundary until the wavefronts cannot expand any more.
	 * To the domains inside the hyper-cube one sub-id is assigned. This procedure continue until
	 * all the domain of one id has a sub-id
	 *
	 * \tparam j property containing the decomposition
	 * \tparam i property to fill with the sub-decomposition
	 *
	 * \param Seed point
	 * \param graph we are processing
	 *
	 */

	template <unsigned int p_sub, unsigned int p_id> void optimize(grid_key_dx<dim> & start_p, Graph & graph)
	{
		// sub-domain id
		size_t sub_id =  0;

		// queue
		openfpm::vector<size_t> v_q;

		// Create an hyper-cube
		HyperCube<dim> hyp;

		// Get the wavefront combinations
		std::vector<comb<dim>> w_comb = hyp.getCombinations_R(dim-1);

		// wavefronts
		openfpm::vector<wavefront<dim>> v_w(w_comb.size());

		// Initialize the wavefronts from the domain start_p

		for (int i = 0 ; i < v_w.size() ; i++)
		{
			for (int j = 0 ; j < dim ; j++)
			{
				v_w.template get<wavefront<dim>::start>(i)[j] = start_p.get(start_p.get(j));
				v_w.template get<wavefront<dim>::stop>(i)[j] = start_p.get(start_p.get(j));
			}
		}

		// push the first domain

		v_q.add(gh.LinId(start_p));

		while (v_q.size() != 0)
		{
			// Box
			Box<dim,size_t> box;

			// Create the biggest box containing the domain
			expand_from_point<p_sub,p_id>(gh.LinId(start_p),graph,box,v_w,w_comb);

			// fill the domain
			fill_domain<p_sub>(graph,box,sub_id);

			// create the queue
			add_to_queue<p_sub>(v_q,v_w,graph,w_comb);

			// increment the sub_id
			sub_id++;
		}
	}
};

#endif
