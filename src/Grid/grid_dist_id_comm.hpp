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

	/*! \brief Reconstruct the local grids
	 *
	 * \param m_oGrid_recv Vector of labeled grids to combine into a local grid
	 */
	inline void grids_reconstruct(openfpm::vector<openfpm::vector<device_grid>> & m_oGrid_recv, openfpm::vector<device_grid> & loc_grid)
	{

	}




	/*! \brief Label intersection grids for mappings
	 *
	 * \param prc_sz For each processor the number of grids to send to
	 */
	inline void labelIntersectionGridsProcessor(Box<dim,St> domain, Decomposition & dec, CellDecomposer_sm<dim,St,shift<dim,St>> & cd_sm, openfpm::vector<device_grid> & loc_grid_old, openfpm::vector<GBoxes<device_grid::dims>> & gdb_ext, openfpm::vector<GBoxes<device_grid::dims>> & gdb_ext_old, openfpm::vector<GBoxes<device_grid::dims>> & gdb_ext_global, openfpm::vector<openfpm::vector<device_grid>> & lbl_b, openfpm::vector<size_t> & prc_sz)
	{
		// resize the label buffer
		lbl_b.resize(v_cl.getProcessingUnits());

		// Label all the intersection grids with the processor id where they should go

		for (size_t i = 0; i < gdb_ext_old.size(); i++)
		{
			// Local old sub-domain
			SpaceBox<dim,long int> sub_dom = gdb_ext_old.get(i).Dbox;
			sub_dom += gdb_ext_old.get(i).origin;

			for (size_t j = 0; j < gdb_ext_global.size(); j++)
			{
				size_t p_id = 0;

				// Intersection box
				SpaceBox<dim,long int> inte_box;

				// Global new sub-domain
				SpaceBox<dim,long int> sub_dom_new = gdb_ext_global.get(j).Dbox;
				sub_dom_new += gdb_ext_global.get(j).origin;

				bool intersect = sub_dom.Intersect(sub_dom_new, inte_box);

				if (intersect == true)
				{

					//std::cout << "Inte_box: (" << inte_box.getLow(0) << "; " << inte_box.getLow(1) << "); (" << inte_box.getHigh(0) << "; " << inte_box.getHigh(1) << ")" << std::endl;

/*
					// Grid to send size
					size_t sz1[dim];
					for (size_t l = 0; l < dim; l++)
					{
						sz1[l] = inte_box.getHigh(l) - inte_box.getLow(l);
						std::cout << "Cont. size on " << l << " dimension: " << sz1[l] << std::endl;
					}
*/

					auto inte_box_cont = cd_sm.convertCellUnitsIntoDomainSpace(inte_box);

					// Get processor ID that store intersection box
					Point<dim,St> p;
					for (size_t n = 0; n < dim; n++)
						p.get(n) = (inte_box_cont.getHigh(n) + inte_box_cont.getLow(n))/2;

					std::cout << "Point: (" << p.get(0) << "; " << p.get(1) << ")" << std::endl;

					p_id = dec.processorID(p);
					prc_sz.get(p_id)++;

					std::cout << "P_id: " << p_id << std::endl;

					// Convert intersection box from contiguous to discrete
					//SpaceBox<dim,long int> inte_box_discr = cd_sm.convertDomainSpaceIntoGridUnits(inte_box,dec.periodicity());


					//std::cout << "Beg:" << inte_box_discr.getHigh(0) << "; " << inte_box.getHigh(1) << std::endl;
					//std::cout << "End:" << inte_box_discr.getLow(0) << "; " << inte_box.getLow(1) << std::endl;

					// Transform coordinates to local
					auto inte_box_local = inte_box;

					inte_box_local -= gdb_ext_global.get(j).origin;

					//std::cout << "Inte_box_local: (" << inte_box_local.getLow(0) << "; " << inte_box_local.getLow(1) << "); (" << inte_box_local.getHigh(0) << "; " << inte_box_local.getHigh(1) << ")" << std::endl;

					// Grid corresponding for gdb_ext_old.get(i) box
					device_grid gr = loc_grid_old.get(i);

					for (size_t l = 0; l < dim; l++)
						std::cout << "GR Size on " << l << " dimension: " << gr.getGrid().size(l) << std::endl;
					// Size of the grid to send
					size_t sz[dim];
					for (size_t l = 0; l < dim; l++)
					{
						sz[l] = inte_box_local.getHigh(l) - inte_box_local.getLow(l) + 1;
						std::cout << "GR_send size on " << l << " dimension: " << sz[l] << std::endl;
					}

					// Grid to send
					device_grid gr_send(sz);
					gr_send.setMemory();

					// Sub iterator across intersection box inside local grid
					grid_key_dx<dim> start = inte_box_local.getKP1();
					grid_key_dx<dim> stop = inte_box_local.getKP2();

					//std::string start2 = start.to_string();
					//std::string stop2 = stop.to_string();

					auto it = gr.getSubIterator(start,stop);

					// Copy selected elements into a new sub-grid
					while (it.isNext())
					{
						auto key = it.get();
						std::string str = key.to_string();

						std::cout << "Key: " << str << std::endl;
						gr_send.get_o(key) = gr.get_o(key);

						//gr_send.template get<0>(key) = gr.template get<0>(key);
						//gr.template get<0>(key)

						//gr_send.set(key,gr,key);

						++it;
					}
					// Add to the labeling vector
					lbl_b.get(p_id).add(gr_send);
					//std::cout << "9" << std::endl;

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
	void map_(Box<dim,St> domain, Decomposition & dec, CellDecomposer_sm<dim,St,shift<dim,St>> & cd_sm, openfpm::vector<device_grid> & loc_grid, openfpm::vector<device_grid> & loc_grid_old, openfpm::vector<GBoxes<device_grid::dims>> & gdb_ext, openfpm::vector<GBoxes<device_grid::dims>> & gdb_ext_old, openfpm::vector<GBoxes<device_grid::dims>> & gdb_ext_global)
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

		// Vector for receiving of intersection grids
		openfpm::vector<openfpm::vector<device_grid>> m_oGrid_recv;

		m_oGrid_recv.resize(m_oGrid.size());
		for (size_t i = 0; i < m_oGrid.size(); i++)
		{
			m_oGrid_recv.get(i).resize(m_oGrid.get(i).size());
		}

		// Send and recieve intersection grids
		v_cl.SSendRecv(m_oGrid,m_oGrid_recv,prc_r,prc_recv_map,recv_sz_map);

		std::cout << "m_oGrid.size(): " << m_oGrid.size() << std::endl;

		std::cout << "m_oGrid_recv.size(): " << m_oGrid_recv.size() << std::endl;

		// Reconstruct the new local grids
		//grids_reconstruct(m_oGrid_recv,loc_grid);
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
