/*
 * grid_dist_id_comm.hpp
 *
 *  Created on: Nov 13, 2016
 *      Author: yaroslav
 */

#ifndef SRC_GRID_GRID_DIST_ID_COMM_HPP_
#define SRC_GRID_GRID_DIST_ID_COMM_HPP_

#include "Vector/vector_dist_ofb.hpp"

/*! \brief This class is an helper for the communication of grid_dist_id
 *
 * \tparam dim Dimensionality of the grid
 * \tparam St Type of space where the grid is living
 * \tparam T object the grid is storing
 * \tparam Decomposition Class that decompose the grid for example CartDecomposition
 * \tparam Memory Is the allocator
 * \tparam device_grid of base structure is going to store the data
 *
 * \see grid_dist_id
 *
 */

template<unsigned int dim, typename St, typename T, typename Decomposition = CartDecomposition<dim,St>,typename Memory=HeapMemory , typename device_grid=grid_cpu<dim,T> >
class grid_dist_id_comm
{
	//! VCluster
	Vcluster & v_cl;

	//! Domain decomposition
	Decomposition dec;

	//! Maps the processor id with the communication request into map procedure
	openfpm::vector<size_t> p_map_req;

	//! Stores the list of processors that communicate with us (local processor)
	openfpm::vector<size_t> prc_recv_map;

	//! Stores the size of the elements added for each processor that communicate with us (local processor)
	openfpm::vector<size_t> recv_sz_map;

	//! For each near processor, outgoing Box
	//! \warning m_oBox is assumed to be an ordered list
	//! first id point Box
	//! second id is the processor id
	openfpm::vector<aggregate<Box<dim,St>,size_t>> m_oBox;

public:

	//! It process one particle
	template<typename T1, typename T2, typename T3, typename T4> inline void process_map_particle(size_t i, long int & end, long int & id_end, T1 & m_pos, T2 & m_prp, T3 & v_pos, T4 & v_prp, openfpm::vector<size_t> & cnt)
	{
/*
		long int prc_id = m_oBox.template get<1>(i);
		size_t id = m_oBox.template get<0>(i);

		if (prc_id >= 0)
		{
			size_t lbl = p_map_req.get(prc_id);

			m_pos.get(lbl).set(cnt.get(lbl), v_pos.get(id));

			cnt.get(lbl)++;

			// swap the particle
			long int id_valid = get_end_valid(end,id_end);

			if (id_valid > 0 && (long int)id < id_valid)
			{
				v_pos.set(id,v_pos.get(id_valid));
				v_prp.set(id,v_prp.get(id_valid));
			}
		}
		else
		{
			// swap the particle
			long int id_valid = get_end_valid(end,id_end);

			if (id_valid > 0 && (long int)id < id_valid)
			{
				v_pos.set(id,v_pos.get(id_valid));
				v_prp.set(id,v_prp.get(id_valid));
			}
		}
*/
	}

	/*! \brief Allocates and fills the send buffer for the map function
	 *
	 * \param v_pos vector of particle positions
	 * \param v_prp vector of particle properties
	 * \param prc_r List of processor rank involved in the send
	 * \param prc_sz_r For each processor in the list the size of the message to send
	 * \param pb send buffer
	 *
	 */
	void fill_send_map_buf(openfpm::vector<device_grid> & loc_grid, openfpm::vector<size_t> & prc_sz_r, openfpm::vector<openfpm::vector<Box<dim,St>>> & m_box)
	{
		m_box.resize(prc_sz_r.size());
		openfpm::vector<size_t> cnt(prc_sz_r.size());

		for (size_t i = 0; i < prc_sz_r.size(); i++)
		{
			// set the size and allocate, using mem warant that pos and prp is contiguous
			m_box.get(i).resize(prc_sz_r.get(i));
			cnt.get(i) = 0;
		}
/*
		// end vector point
		long int id_end = v_pos.size();

		// end opart point
		long int end = m_opart.size()-1;

		// Run through all the particles and fill the sending buffer
		for (size_t i = 0; i < m_opart.size(); i++)
		{
			process_map_particle<proc_with_prp<prp_object,prp...>>(i,end,id_end,m_pos,m_prp,v_pos,v_prp,cnt);
		}

		v_pos.resize(v_pos.size() - m_opart.size());
*/
	}

	/*! \brief Label particles for mappings
	 *
	 * \param v_pos vector of particle positions
	 * \param lbl_p Particle labeled
	 * \param prc_sz For each processor the number of particles to send
	 * \param opart id of the particles to send
	 *
	 */
	void labelIntersectionGridsProcessor(openfpm::vector<SpaceBox<dim, T>> & sub_domains_old, openfpm::vector<device_grid> & loc_grid, openfpm::vector<GBoxes<device_grid::dims>> & gdb_ext, openfpm::vector<aggregate<Box<dim,St>,size_t>> & lbl_b, openfpm::vector<size_t> & prc_sz)
	{
		// reset lbl_b
		lbl_b.clear();

		// resize the label buffer
		prc_sz.resize(v_cl.getProcessingUnits());

		// Label all the intersection boxes with the processor id where they should go
		for (size_t i = 0; i < dec.getNSubDomain(); i++)
		{
			for (size_t j = 0; j < sub_domains_old.size(); j++)
			{
				size_t p_id = 0;

				Box<dim,St> inte_box;
				bool intersect = dec.getSubDomain(i).Intersect(sub_domains_old.get(j), inte_box);

				if (intersect == true)
				{
					p_id = dec.processorID(inte_box.rnd());
					prc_sz.get(p_id)++;

					lbl_b.add();
					////////////////////
					lbl_b.last().template get<0>() = inte_box;
					////////////////////
					lbl_b.last().template get<1>() = p_id;
				}
			}
		}
	}

	/*! \brief Moves all the grids that does not belong to the local processor to the respective processor
	 *
	 * \tparam out of bound policy it specify what to do when the particles are detected out of bound
	 *
	 * In general this function is called after moving the particles to move the
	 * elements out the local processor. Or just after initialization if each processor
	 * contain non local particles
	 *
	 * \param v_pos vector of particle positions
	 * \param v_prp vector of particle properties
	 * \param g_m ghost marker
	 *
	 */
	void map_(openfpm::vector<SpaceBox<dim, T>> & sub_domains_old, openfpm::vector<device_grid> & loc_grid, openfpm::vector<GBoxes<device_grid::dims>> & gdb_ext)
	{
		// Processor communication size
		openfpm::vector<size_t> prc_sz(v_cl.getProcessingUnits());

		// Contains the processor id of each grid (basically where they have to go)
		labelIntersectionGridsProcessor(sub_domains_old,loc_grid,gdb_ext,m_oBox,prc_sz);

		// Calculate the sending buffer size for each processor, put this information in
		// a contiguous buffer
		p_map_req.resize(v_cl.getProcessingUnits());

		// Vector of number of boxes for each involved processor
		openfpm::vector<size_t> prc_sz_r;
		// Vector of ranks of involved processors
		openfpm::vector<size_t> prc_r;

		for (size_t i = 0; i < v_cl.getProcessingUnits(); i++)
		{
			if (prc_sz.get(i) != 0)
			{
				p_map_req.get(i) = prc_r.size();
				prc_r.add(i);
				prc_sz_r.add(prc_sz.get(i));
			}
		}

		//! Grids vector
		openfpm::vector<openfpm::vector<Box<dim,St>>> m_box;

		fill_send_map_buf(loc_grid, prc_sz_r, m_box);

		v_cl.SSendRecv(m_box,loc_grid,prc_r,prc_recv_map,recv_sz_map);
	}

	/*! \brief Constructor
	 *
	 * \param dec Domain decompositon
	 *
	 */
	grid_dist_id_comm(const Decomposition & dec)
	:v_cl(create_vcluster()),dec(dec)
	{

	}
};


#endif /* SRC_GRID_GRID_DIST_ID_COMM_HPP_ */
