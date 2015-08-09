/*
 * nn_processor.hpp
 *
 *  Created on: Aug 9, 2015
 *      Author: i-bird
 */

#ifndef SRC_DECOMPOSITION_NN_PROCESSOR_HPP_
#define SRC_DECOMPOSITION_NN_PROCESSOR_HPP_

#include "common.hpp"

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

	// for each near-processor store the sub-domain of the near processor
	std::unordered_map<size_t, N_box<dim,T>> nn_processor_subdomains;

	// for each processor store the set of the sub-domains sent to the adjacent processors
	openfpm::vector<openfpm::vector<size_t>> proc_adj_box;

	//! contain the internal adjacent sub-domains sent to the other processors
	openfpm::vector< openfpm::vector< ::SpaceBox<dim,T>> > boxes;

	// Receive counter
	size_t recv_cnt;

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

public:

	nn_prcs(Vcluster & v_cl)
	:v_cl(v_cl){}

	/*! \brief Create the list of adjacent processors and the list of adjacent sub-domains
	 *
	 * \param box_nn_processors
	 *
	 */
	void create(const openfpm::vector<openfpm::vector<long unsigned int> > & box_nn_processor, const openfpm::vector<SpaceBox<dim,T>> & sub_domains)
	{
		// produce the list of the contiguous processor (nn_processors) and link nn_processor_subdomains to the
		// processor list
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

		// create a buffer with the sub-domains of this processor, the informations ( the boxes )
		// of the sub-domains contiguous to the processor A are sent to the processor A and
		// the information of the contiguous sub-domains in the near processors are received
		//
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

		// Intersect all the local sub-domains with the sub-domains of the contiguous processors

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
	inline size_t IDtoProc(size_t id)
	{
		return nn_processors.get(id);
	}

	/*! \brief Get the sub-domain pf an adjacent processor
	 *
	 * \param p_id adjacent processor (id from 0 to getNNProcessors())
	 *
	 * \return the sub-domains
	 *
	 */
	inline const openfpm::vector< ::Box<dim,T> > & getAdjacentSubdomain(size_t p_id)
	{
		return nn_processor_subdomains[p_id].bx;
	}

	/*! \brief Get the adjacent processor id
	 *
	 * \param p_id adjacent processor (id from 0 to getNNProcessors())
	 *
	 * \return the processor rank
	 *
	 */
	inline size_t getAdjacentProcessor(size_t p_id)
	{
		return nn_processor_subdomains[p_id].id;
	}


	/*! \brief Get the local sub-domains adjacent to a processor p_id
	 *
	 * \param p_id adjacent processor (id from 0 to getNNProcessors())
	 *
	 * \return the sub-domains
	 *
	 */
	inline const openfpm::vector<size_t> & getInternalAdjSubdomain(size_t p_id)
	{
		return proc_adj_box.get(p_id);
	}

	/*! \brief Convert the processor rank to the id in the list
	 *
	 * \param p processor rank
	 *
	 * \return the id
	 *
	 */
	inline size_t ProctoID(size_t p)
	{
		return nn_processor_subdomains[p].id;
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

};


#endif /* SRC_DECOMPOSITION_NN_PROCESSOR_HPP_ */
