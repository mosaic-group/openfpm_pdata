/*
 * nn_processor.hpp
 *
 *  Created on: Aug 9, 2015
 *      Author: i-bird
 */

#ifndef SRC_DECOMPOSITION_NN_PROCESSOR_HPP_
#define SRC_DECOMPOSITION_NN_PROCESSOR_HPP_

#include "common.hpp"
#include <unordered_map>

/*! \brief This class store the adjacent processors and the adjacent sub_domains
 *
 * \tparam dim is the dimensionality of the physical domain we are going to decompose.
 * \tparam T type of the space we decompose, Real, Integer, Complex ...
 *
 * \see CartDecomposition
 *
 */
template<unsigned int dim, typename T>
class nn_prcs
{
	//! Virtual cluster
	Vcluster & v_cl;

	//! List of adjacent processors
	openfpm::vector<size_t> nn_processors;

	// for each adjacent processor store the sub-domains of the adjacent processor
	std::unordered_map<size_t, N_box<dim,T>> nn_processor_subdomains;

	// when we add new boxes, are added here
	std::unordered_map<size_t, N_box<dim,T>> nn_processor_subdomains_tmp;

	// for each near processor store the set of the sub-domains sent to the near processors,
	// the nested vector contain the id of the local sub-domain
	openfpm::vector<openfpm::vector<size_t>> proc_adj_box;

	//! contain the internal adjacent sub-domains sent to the other processors
	openfpm::vector< openfpm::vector< ::SpaceBox<dim,T>> > boxes;

	// Receive counter
	size_t recv_cnt;

	//! applyBC function is suppose to be called only one time
	bool aBC;

	/*! \brief Message allocation
	 *
	 * \param message size required to receive from i
	 * \param total message size to receive from all the processors
	 * \param the total number of processor want to communicate with you
	 * \param i processor id
	 * \param ri request id (it is an id that goes from 0 to total_p, and is unique
	 *           every time message_alloc is called)
	 * \param ptr a pointer to the vector_dist structure
	 *
	 * \return the pointer where to store the message
	 *
	 */
	static void * message_alloc(size_t msg_i ,size_t total_msg, size_t total_p, size_t i, size_t ri, void * ptr)
	{
		// cast the pointer
		nn_prcs<dim,T> * cd = static_cast< nn_prcs<dim,T> *>(ptr);

		// Resize the memory
		cd->nn_processor_subdomains[i].bx.resize(msg_i / sizeof(::Box<dim,T>) );

		// Return the receive pointer
		return cd->nn_processor_subdomains[i].bx.getPointer();
	}

	/*! Checl that the compination is valid
	 *
	 * \param cmb combination
	 * \param bc boundary conditions
	 *
	 */
	bool inline check_valid(comb<dim> cmb,const size_t (& bc)[dim])
	{
		for (size_t i = 0 ; i < dim ; i++)
		{
			if (bc[i] == NON_PERIODIC && cmb.getComb()[i] != 0)
				return false;
		}
		return true;
	}

	/*! \brief add sub-domain to processor
	 *
	 * \param i processor
	 * \param bx Box to add
	 *
	 */
	inline void add_nn_subdomain(size_t i, const Box<dim,T> & bx)
	{
		nn_processor_subdomains_tmp[i].bx.add(bx);
	}

	/*! \brief In case of periodic boundary conditions we have to add boxes
	 *         at the borders
	 *
	 * \param domain Domain box
	 * \param boundary boundary conditions
	 * \param ghost ghost part
	 *
	 */
	void add_box_periodic(const Box<dim,T> & domain, const Ghost<dim,T> & ghost, const size_t (&bc)[dim])
	{
		HyperCube<dim> hyp;

		// first we create boxes at the border of the domain used to detect the sub-domain
		// that must be adjusted, each of this boxes define a shift in case of periodic boundary condition
		for (long int i = dim-1 ; i >= 0 ; i--)
		{
			std::vector<comb<dim>> cmbs = hyp.getCombinations_R(i);

			for (size_t j = 0 ; j < cmbs.size() ; j++)
			{
				if (check_valid(cmbs[j],bc) == false)
					continue;

				Box<dim,T> bp;
				Point<dim,T> shift;

				for (size_t k = 0 ; k < dim ; k++)
				{
					switch (cmbs[j][k])
					{
					case 1:
						bp.setLow(k,domain.getHigh(k)+ghost.getLow(k));
						bp.setHigh(k,domain.getHigh(k));
						shift.get(k) = -domain.getHigh(k);
						break;
					case 0:
						bp.setLow(k,domain.getLow(k));
						bp.setHigh(k,domain.getHigh(k));
						shift.get(k) = 0;
						break;
					case -1:
						bp.setLow(k,domain.getLow(k));
						bp.setHigh(k,ghost.getHigh(k));
						shift.get(k) = domain.getHigh(k);
						break;
					}
				}

				// Detect all the sub-domain involved, shift them and add to the list
				// Detection is performed intersecting the sub-domains with the ghost
				// parts near the domain borders
				for (size_t k = 0 ; k < getNNProcessors() ; k++)
				{
					// sub-domains of the near processor
					const openfpm::vector< ::Box<dim,T> > & nn_sub = getAdjacentSubdomains(IDtoProc(k));

					for (size_t l = 0 ; l < nn_sub.size(); l++)
					{
						Box<dim,T> sub = nn_sub.get(l);
						Box<dim,T> b_int;

						if (sub.Intersect(bp,b_int) == true)
						{
							sub += shift;
							add_nn_subdomain(IDtoProc(k),sub);
						}
					}
				}
			}
		}

		flush();
	}

	/*! \brief Flush the temporal added sub-domain to the processor sub-domain
	 *
	 *
	 */
	void flush()
	{
		for ( auto it = nn_processor_subdomains_tmp.begin(); it != nn_processor_subdomains_tmp.end(); ++it )
		{
			const N_box<dim,T> & nnp_bx = it->second;

			for (size_t i = 0 ; i < nnp_bx.bx.size() ; i++)
			{
				nn_processor_subdomains[it->first].bx.add(nn_processor_subdomains_tmp[it->first].bx.get(i));
			}
		}
	}

public:

	nn_prcs(Vcluster & v_cl)
	:v_cl(v_cl),aBC(false){}

	//! Constructor from another nn_prcs
	nn_prcs(const nn_prcs<dim,T> & ilg)
	:v_cl(ilg.v_cl),aBC(false)
	{
		this->operator=(ilg);
	};

	//! Constructor from temporal ie_loc_ghost
	nn_prcs(nn_prcs<dim,T> && ilg)
	:v_cl(ilg.v_cl),aBC(false)
	{
		this->operator=(ilg);
	}

	/*! \brief Copy the object
	 *
	 * \param nnp object to copy
	 *
	 */
	nn_prcs<dim,T> & operator=(const nn_prcs<dim,T> & nnp)
	{
		v_cl = nnp.v_cl;
		nn_processors = nnp.nn_processors;
		nn_processor_subdomains = nnp.nn_processor_subdomains;
		proc_adj_box = nnp.proc_adj_box;
		boxes = nnp.boxes;

		return *this;
	}

	/*! \brief Copy the object
	 *
	 * \param nnp object to copy
	 *
	 */
	nn_prcs<dim,T> & operator=(nn_prcs<dim,T> && nnp)
	{
		v_cl = nnp.v_cl;
		nn_processors.swap(nnp.nn_processors);
		nn_processor_subdomains.swap(nnp.nn_processor_subdomains);
		proc_adj_box.swap(nnp.proc_adj_box);
		boxes = nnp.boxes;

		return *this;
	}

	/*! \brief Refine the ss_box to have the smallest size on each direction of the local subdomain and adjacent (from other processor) one
	 *
	 * \param ss_box box that store the smallest size of the sub-domain
	 *
	 */
	void refine_ss_box(Box<dim,T> & ss_box)
	{
		for (size_t p = 0 ; p < getNNProcessors() ; p++)
		{
			auto list_p_box = getExternalAdjSubdomain(IDtoProc(p));

			// Create the smallest box contained in all sub-domain
			for (size_t b = 0 ; b < list_p_box.bx.size() ; b++)
				ss_box.contained(list_p_box.bx.get(b));
		}
	}

	/*! \brief Create the list of adjacent processors and the list of adjacent sub-domains
	 *
	 * \param box_nn_processors
	 *
	 */
	void create(const openfpm::vector<openfpm::vector<long unsigned int> > & box_nn_processor, const openfpm::vector<SpaceBox<dim,T>> & sub_domains)
	{
		// produce the list of the adjacent processor (nn_processors) list
		for (size_t i = 0 ;  i < box_nn_processor.size() ; i++)
		{
			for (size_t j = 0 ; j < box_nn_processor.get(i).size() ; j++)
			{
				nn_processors.add(box_nn_processor.get(i).get(j));
			}
		}

		// make the list sorted and unique
	    std::sort(nn_processors.begin(), nn_processors.end());
	    auto last = std::unique(nn_processors.begin(), nn_processors.end());
	    nn_processors.erase(last, nn_processors.end());

        // construct nn_processor_subdomains
        for (size_t i = 0 ;  i < box_nn_processor.size() ; i++)
        {
                for (size_t j = 0 ; j < box_nn_processor.get(i).size() ; j++)
                {
                        // processor id adjacent to this sub-domain
                        size_t proc_id = box_nn_processor.get(i).get(j);

                        size_t k = 0;
                        // search inside near processor list
                        for (k = 0 ; k < nn_processors.size() ; k++)
                                if (nn_processors.get(k) == proc_id)    break;

                        nn_processor_subdomains[proc_id].id = k;
                }
        }


		// create a buffer with the sub-domains of this processor, send them to the adjacent
		// processor
		proc_adj_box.resize(getNNProcessors());
		boxes.resize(nn_processors.size());

		for (size_t b = 0 ; b < box_nn_processor.size() ; b++)
		{
			for (size_t p = 0 ; p < box_nn_processor.get(b).size() ; p++)
			{
				size_t prc = box_nn_processor.get(b).get(p);

				// id of the processor in the processor list
				// [value between 0 and the number of the near processors]
				size_t id = nn_processor_subdomains[prc].id;

				boxes.get(id).add(sub_domains.get(b));
				proc_adj_box.get(id).add(b);
			}
		}

		// Get the sub-domains of the near processors
		v_cl.sendrecvMultipleMessagesNBX(nn_processors,boxes,nn_prcs<dim,T>::message_alloc, this ,NEED_ALL_SIZE);
	}

	/*! \brief Get the number of Near processors
	 *
	 * \return the number of near processors
	 *
	 */
	inline size_t getNNProcessors() const
	{
		return nn_processors.size();
	}

	/*! \brief Return the processor id of the near processor list at place id
	 *
	 * \param id
	 *
	 * \return return the processor rank
	 *
	 */
	inline size_t IDtoProc(size_t id) const
	{
		return nn_processors.get(id);
	}

	/*! \brief Get the sub-domain pf an adjacent processor
	 *
	 * \param p_id adjacent processor rank
	 *
	 * \return the sub-domains
	 *
	 */
	inline const openfpm::vector< ::Box<dim,T> > & getAdjacentSubdomains(size_t p_id) const
	{
		auto key = nn_processor_subdomains.find(p_id);
#ifdef SE_CLASS1
		if (key == nn_processor_subdomains.end())
		{
			std::cerr << "Error " << __FILE__ << ":" << __LINE__ << " error this process rank is not adjacent to the local processor";
		}
#endif
		return key->second.bx;
	}

	/*! \brief Get the adjacent processor id
	 *
	 * \param p_id adjacent processor rank
	 *
	 * \return the processor rank
	 *
	 */
	inline size_t getAdjacentProcessor(size_t p_id) const
	{
		auto key = nn_processor_subdomains.find(p_id);
#ifdef SE_CLASS1
		if (key == nn_processor_subdomains.end())
		{
			std::cerr << "Error " << __FILE__ << ":" << __LINE__ << " error this process rank is not adjacent to the local processor";
		}
#endif
		return key->second.id;
	}


	/*! \brief Get the local sub-domains adjacent to a processor p_id
	 *
	 * \param p_id adjacent processor (id from 0 to getNNProcessors())
	 *
	 * \return the sub-domains
	 *
	 */
	inline const openfpm::vector<size_t> & getInternalAdjSubdomain(size_t p_id) const
	{
		return proc_adj_box.get(p_id);
	}

	/*! \brief Get the external sub-domain adjacent to a processor p_id
	 *
	 * \param p_id processor rank
	 * \return the set of adjacent sub-domain comming from the processor p_id
	 *
	 */
	inline const N_box<dim,T> & getExternalAdjSubdomain(size_t p_id) const
	{
		auto key = nn_processor_subdomains.find(p_id);
#ifdef SE_CLASS1
		if (key == nn_processor_subdomains.end())
		{
			std::cerr << "Error " << __FILE__ << ":" << __LINE__ << " error this process rank is not adjacent to the local processor";
		}
#endif
		return key->second;
	}

	/*! \brief Convert the processor rank to the id in the list
	 *
	 * \param p processor rank
	 *
	 * \return the id
	 *
	 */
	inline size_t ProctoID(size_t p) const
	{
		auto key = nn_processor_subdomains.find(p);
#ifdef SE_CLASS1
		if (key == nn_processor_subdomains.end())
		{
			std::cerr << "Error " << __FILE__ << ":" << __LINE__ << " error this process rank is not adjacent to the local processor";
		}
#endif

		return key->second.id;
	}

	/*! \brief Write the decomposition as VTK file
	 *
	 * The function generate several files
	 *
	 * 1) subdomains_adjacent_X.vtk sub-domains adjacent to the local processor (X)
	 *
	 * where X is the local processor rank
	 *
	 * \param output directory where to write the files
	 * \param p_id id of the local processor
	 *
	 */
	bool write(std::string output) const
	{
		//! subdomains_adjacent_X.vtk sub-domains adjacent to the local processor (X)
		VTKWriter<openfpm::vector<::Box<dim,T>>,VECTOR_BOX> vtk_box2;
		for (size_t p = 0 ; p < nn_processors.size() ; p++)
		{
			size_t prc = nn_processors.get(p);
			auto it = nn_processor_subdomains.find(prc);
			if (it != nn_processor_subdomains.end())
				vtk_box2.add(nn_processor_subdomains.at(prc).bx);
		}
		vtk_box2.write(output + std::string("subdomains_adjacent_") + std::to_string(v_cl.getProcessUnitID()) + std::string(".vtk"));

		return true;
	}

	/*! \brief Apply boundary conditions
	 *
	 * \param domain The simulation domain
	 * \param ghost ghost part
	 * \param bc Boundary conditions
	 *
	 */
	void applyBC(const Box<dim,T> & domain, const Ghost<dim,T> & ghost, const size_t (&bc)[dim])
	{
		if (aBC == true)
		{
			std::cerr << "Warning " << __FILE__ << ":" << __LINE__ << "Apply BC is suppose to be called only one time\n";
			return;
		}

		aBC=true;

		return add_box_periodic(domain,ghost,bc);
	}

	/*! \brief Check if the nn_prcs contain the same information
	 *
	 * \param ele Element to check
	 *
	 */
	bool is_equal(nn_prcs<dim,T> & np)
	{
		if (np.getNNProcessors() != getNNProcessors())
			return false;

		for (size_t p = 0 ; p < getNNProcessors() ; p++)
		{
			if (getAdjacentSubdomains(IDtoProc(p)) != np.getAdjacentSubdomains(IDtoProc(p)))
				return false;
			if (getAdjacentProcessor(IDtoProc(p)) != np.getAdjacentProcessor(IDtoProc(p)))
				return false;
			if (getInternalAdjSubdomain(p) != np.getInternalAdjSubdomain(p))
				return false;
			if (getExternalAdjSubdomain(IDtoProc(p)) != np.getExternalAdjSubdomain(IDtoProc(p)))
				return false;
		}

		return true;
	}

	//! Used for testing porpose do not use
	std::unordered_map<size_t, N_box<dim,T>> & get_nn_processor_subdomains()
	{
		return nn_processor_subdomains;
	}

	//! Used for testing porpose do not use
	openfpm::vector<size_t> & get_nn_processors()
	{
		return nn_processors;
	}
};


#endif /* SRC_DECOMPOSITION_NN_PROCESSOR_HPP_ */
