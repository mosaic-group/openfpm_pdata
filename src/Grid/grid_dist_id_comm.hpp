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

	//! Maps the processor id with the communication request into map procedure
	openfpm::vector<size_t> p_map_req;

	//! Stores the list of processors that communicate with us (local processor)
	openfpm::vector<size_t> prc_recv_map;

	//! Stores the size of the elements added for each processor that communicate with us (local processor)
	openfpm::vector<size_t> recv_sz_map;

	//! For each near processor, outgoing intersection grid
	//! \warning m_oGrid is assumed to be an ordered list
	//! first id is grid
	//! second id is the processor id
	openfpm::vector<openfpm::vector<device_grid>> m_oGrid;

public:

	//! It processes one box
	template<typename T1, typename T2, typename T3, typename T4> inline void process_map_box(size_t i, long int & end, long int & id_end, T1 & m_pos, T2 & m_prp, T3 & v_pos, T4 & v_prp, openfpm::vector<size_t> & cnt)
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
	 *
	 */
	void fill_send_map_buf_(openfpm::vector<device_grid> & loc_grid, openfpm::vector<size_t> & prc_sz_r, openfpm::vector<openfpm::vector<SpaceBox<dim,St>>> & m_box)
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

	/*! \brief Allocates and fills the send buffer for the map function
	 *
	 *
	 */
	void fill_send_map_buf(openfpm::vector<device_grid> & loc_grid, openfpm::vector<size_t> & prc_sz_r, openfpm::vector<openfpm::vector<SpaceBox<dim,St>>> & m_box)
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

	/*! \brief Label intersection grids for mappings
	 *
	 * \param prc_sz For each processor the number of grids to send
	 */
	void labelIntersectionGridsProcessor(Box<dim,St> domain, Decomposition & dec, CellDecomposer_sm<dim,St,shift<dim,St>> & cd_sm, openfpm::vector<device_grid> & loc_grid_old, openfpm::vector<GBoxes<device_grid::dims>> & gdb_ext, openfpm::vector<GBoxes<device_grid::dims>> & gdb_ext_old, openfpm::vector<GBoxes<device_grid::dims>> & gdb_ext_global, openfpm::vector<openfpm::vector<device_grid>> & lbl_b, openfpm::vector<size_t> & prc_sz)
	{
		// resize the label buffer
		lbl_b.resize(v_cl.getProcessingUnits());

		// Label all the intersection grids with the processor id where they should go

		for (size_t i = 0; i < gdb_ext_old.size(); i++)
		{

			for (size_t j = 0; j < gdb_ext_global.size(); j++)
			{
				size_t p_id = 0;

				// Intersection box
				SpaceBox<dim,St> inte_box;

				// Local old sub-domain
				SpaceBox<dim,St> sub_dom = gdb_ext_old.get(i).Dbox;

				// Global new sub-domain
				SpaceBox<dim,St> sub_dom_new = gdb_ext_global.get(j).Dbox;

				bool intersect = sub_dom.Intersect(sub_dom_new, inte_box);

				if (intersect == true)
				{
/*
					// Grid to send size
					size_t sz1[dim];
					for (size_t l = 0; l < dim; l++)
					{
						sz1[l] = inte_box.getHigh(l) - inte_box.getLow(l);
						std::cout << "Cont. size on " << l << " dimension: " << sz1[l] << std::endl;
					}
*/
					// Get processor ID that store intersection box
					Point<dim,St> p;
					for (size_t i = 0; i < dim; i++)
						p.get(i) = (inte_box.getHigh(i) + inte_box.getLow(i))/2;
					p_id = dec.processorID(p);
					prc_sz.get(p_id)++;

					std::cout << "P_id: " << p_id << std::endl;

					// Convert intersection box from contiguous to discrete
					SpaceBox<dim,long int> inte_box_discr = cd_sm.convertDomainSpaceIntoGridUnits(inte_box,dec.periodicity());

					inte_box_discr -= gdb_ext.get(i).origin;

					// Grid corresponding for gdb_ext_old.get(i) box
					device_grid gr = loc_grid_old.get(i);

					// Grid to send size
					size_t sz[dim];
					for (size_t l = 0; l < dim; l++)
					{
						sz[l] = inte_box_discr.getHigh(l) - inte_box_discr.getLow(l);
						std::cout << " Size on " << l << " dimension: " << sz[l] << std::endl;
					}

					// Grid to send
					device_grid gr_send(sz);
					gr_send.setMemory();

					// Sub iterator across intersection box inside local grid
					grid_key_dx<dim> start = inte_box_discr.getKP1();
					grid_key_dx<dim> stop = inte_box_discr.getKP2();

					std::string start2 = start.to_string();
					std::string stop2 = stop.to_string();

					auto key_it = gr.getSubIterator(start,stop);

					// Copy selected elements into a new sub-grid
					while (key_it.isNext())
					{
						grid_key_dx<dim> key = key_it.get();
						//std::string str = key.to_string();
						//std::cout << "Key: " << str << std::endl;

						gr_send.set(key,gr,key);

						++key_it;
					}
					// Add to the labeling vector
					lbl_b.get(p_id).add(gr_send);

					std::cout << "9" << std::endl;
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
	void map_(Box<dim,St> domain, Decomposition & dec, CellDecomposer_sm<dim,St,shift<dim,St>> & cd_sm, openfpm::vector<device_grid> & loc_grid_old, openfpm::vector<GBoxes<device_grid::dims>> & gdb_ext, openfpm::vector<GBoxes<device_grid::dims>> & gdb_ext_old, openfpm::vector<GBoxes<device_grid::dims>> & gdb_ext_global)
	{
		// Processor communication size
		openfpm::vector<size_t> prc_sz(v_cl.getProcessingUnits());

		// Contains the processor id of each box (basically where they have to go)
		labelIntersectionGridsProcessor(domain,dec,cd_sm,loc_grid_old,gdb_ext,gdb_ext_old,gdb_ext_global,m_oGrid,prc_sz);

		// Calculate the sending buffer size for each processor, put this information in
		// a contiguous buffer
		p_map_req.resize(v_cl.getProcessingUnits());

		// Vector of number of sending grids for each involved processor
		openfpm::vector<size_t> prc_sz_r;
		// Vector of ranks of involved processors
		openfpm::vector<size_t> prc_r;

		for (size_t i = 0; i < v_cl.getProcessingUnits(); i++)
		{
			if (m_oGrid.get(i).size() != 0)
			{
				p_map_req.get(i) = prc_r.size();
				prc_r.add(i);
				prc_sz_r.add(m_oGrid.get(i).size());
			}
		}

		//! Boxes vector
		//openfpm::vector<openfpm::vector<device_grid>> m_box;

		//fill_send_map_buf(loc_grid_global, prc_sz_r, m_box);

		openfpm::vector<openfpm::vector<device_grid>> m_oGrid_recv;

		m_oGrid_recv.resize(m_oGrid.size());
		for (size_t i = 0; i < m_oGrid.size(); i++)
		{
			m_oGrid_recv.get(i).resize(m_oGrid.get(i).size());
		}

		v_cl.SSendRecv(m_oGrid,m_oGrid_recv,prc_r,prc_recv_map,recv_sz_map);

		std::cout << "m_oGrid_recv.size(): " << m_oGrid_recv.size() << std::endl;
	}

	/*! \brief Constructor
	 *
	 * \param dec Domain decompositon
	 *
	 */
	grid_dist_id_comm()
	:v_cl(create_vcluster())
	{

	}
};


#endif /* SRC_GRID_GRID_DIST_ID_COMM_HPP_ */
