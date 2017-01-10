/*
 * grid_dist_id_comm.hpp
 *
 *  Created on: Nov 13, 2016
 *      Author: yaroslav
 */

#ifndef SRC_GRID_GRID_DIST_ID_COMM_HPP_
#define SRC_GRID_GRID_DIST_ID_COMM_HPP_

#include "Vector/vector_dist_ofb.hpp"
#include "data_type/scalar.hpp"


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
	openfpm::vector<openfpm::vector<aggregate<device_grid,SpaceBox<dim,long int>>>> m_oGrid;

public:

	/*! \brief Reconstruct the local grids
	 *
	 * \param m_oGrid_recv Vector of labeled grids to combine into a local grid
	 */
	inline void grids_reconstruct(openfpm::vector<openfpm::vector<aggregate<device_grid,SpaceBox<dim,long int>>>> & m_oGrid_recv, openfpm::vector<device_grid> & loc_grid, openfpm::vector<GBoxes<device_grid::dims>> & gdb_ext, CellDecomposer_sm<dim,St,shift<dim,St>> & cd_sm)
	{
		size_t count2 = 0;
		for (size_t a = 0; a < m_oGrid_recv.size(); a++)
		{
			for (size_t k = 0; k < m_oGrid_recv.get(a).size(); k++)
			{
				device_grid g = m_oGrid_recv.get(a).template get<0>(k);

				size_t count = 0;


				auto it = g.getIterator();

				while (it.isNext())
				{
					//auto key = it.get();

					//if (g.template get<0>(key) != 1)
						//std::cout << "WRONG???????" << std::endl;

					++it;
					count++;
				}

				SpaceBox<dim,long int> b = m_oGrid_recv.get(a).template get<1>(k);

				//device_grid gr_send(sz);
				//gr_send.setMemory();

				//std::cout << "B: (" << b.getLow(0) << "; " << b.getLow(1) << "); (" << b.getHigh(0) << "; " << b.getHigh(1) << "); " << "G: (" << g.getGrid().getBox().getHigh(0) << "; " << g.getGrid().getBox().getHigh(1) << ")" << std::endl;

				// Set the dimensions of the local grid
				//g.resize(l_res);

				Point<dim,St> p;
				for (size_t n = 0; n < dim; n++)
					p.get(n) = g.getGrid().getBox().getHigh(n);

				//std::cout << "G after resize: (" << g.getGrid().getBox().getLow(0) << "; " << g.getGrid().getBox().getLow(1) << "); (" << g.getGrid().getBox().getHigh(0) << "; " << g.getGrid().getBox().getHigh(1) << ")" << std::endl;

				Point<dim,St> point;
				for (size_t n = 0; n < dim; n++)
					point.get(n) = (b.getHigh(n) + b.getLow(n))/2;

				for (size_t j = 0; j < gdb_ext.size(); j++)
				{
					// Local sub-domain
					SpaceBox<dim,long int> sub = gdb_ext.get(j).Dbox;
					sub += gdb_ext.get(j).origin;

					if (sub.isInside(point) == true)
					{
						grid_key_dx<dim> start = b.getKP1() - grid_key_dx<dim>(gdb_ext.get(j).origin.asArray());
						grid_key_dx<dim> stop = b.getKP2() - grid_key_dx<dim>(gdb_ext.get(j).origin.asArray());

						std::string start2 = start.to_string();
						std::string stop2 = stop.to_string();

						auto it = loc_grid.get(j).getSubIterator(start,stop);

						// Copy selected elements into a local grid
						while (it.isNext())
						{
							auto key = it.get();
							std::string str = key.to_string();
							grid_key_dx<dim> key2 = key - start;

							//std::cout << "Key: " << str << std::endl;
							loc_grid.get(j).get_o(key) = g.get_o(key2);
							count2++;

							++it;
						}
					}
				}
			}
		}
		//std::cout << "Count after: " << count2 << std::endl;
	}


	/*! \brief Reconstruct the local grids
	 *
	 * \param m_oGrid_recv Vector of labeled grids to combine into a local grid
	 */
	inline void grids_reconstruct(openfpm::vector<openfpm::vector<aggregate<device_grid,SpaceBox<dim,long int>>>> & m_oGrid_recv, openfpm::vector<device_grid> & loc_grid, openfpm::vector<GBoxes<device_grid::dims>> & gdb_ext, CellDecomposer_sm<dim,St,shift<dim,St>> & cd_sm, openfpm::vector<size_t> & prc_r)
	{

	}

	/*! \brief Label intersection grids for mappings
	 *
	 * \param prc_sz For each processor the number of grids to send to
	 */
	inline void labelIntersectionGridsProcessor(Decomposition & dec, CellDecomposer_sm<dim,St,shift<dim,St>> & cd_sm, openfpm::vector<device_grid> & loc_grid_old, openfpm::vector<GBoxes<device_grid::dims>> & gdb_ext, openfpm::vector<GBoxes<device_grid::dims>> & gdb_ext_old, openfpm::vector<GBoxes<device_grid::dims>> & gdb_ext_global, openfpm::vector<openfpm::vector<aggregate<device_grid,SpaceBox<dim,long int>>>> & lbl_b, openfpm::vector<size_t> & prc_sz)
	{
		// resize the label buffer
		lbl_b.resize(v_cl.getProcessingUnits());

		size_t count = 0;
		size_t count2 = 0;

		// Label all the intersection grids with the processor id where they should go

		for (size_t i = 0; i < gdb_ext_old.size(); i++)
		{
			// Local old sub-domain in global coordinates
			SpaceBox<dim,long int> sub_dom = gdb_ext_old.get(i).Dbox;
			sub_dom += gdb_ext_old.get(i).origin;

			for (size_t j = 0; j < gdb_ext_global.size(); j++)
			{
				size_t p_id = 0;

				// Intersection box
				SpaceBox<dim,long int> inte_box;

				// Global new sub-domain in global coordinates
				SpaceBox<dim,long int> sub_dom_new = gdb_ext_global.get(j).Dbox;
				sub_dom_new += gdb_ext_global.get(j).origin;

				bool intersect = false;

				if (sub_dom.isValid() == true && sub_dom_new.isValid() == true)
					intersect = sub_dom.Intersect(sub_dom_new, inte_box);

				if (intersect == true)
				{
					//// DEBUG/////
					count2++;
					//////////////

					//std::cout << "Inte_box: (" << inte_box.getLow(0) << "; " << inte_box.getLow(1) << "); (" << inte_box.getHigh(0) << "; " << inte_box.getHigh(1) << ")" << std::endl;

					auto inte_box_cont = cd_sm.convertCellUnitsIntoDomainSpace(inte_box);

					// Get processor ID that store intersection box
					Point<dim,St> p;
					for (size_t n = 0; n < dim; n++)
						p.get(n) = (inte_box_cont.getHigh(n) + inte_box_cont.getLow(n))/2;

					//std::cout << "Point: (" << p.get(0) << "; " << p.get(1) << ")" << std::endl;

					p_id = dec.processorID(p);
					prc_sz.get(p_id)++;

					//std::cout << "P_id: " << p_id << std::endl;

					// Transform coordinates to local
					auto inte_box_local = inte_box;

					inte_box_local -= gdb_ext_old.get(i).origin;

					//std::cout << "gdb_ext_old.get(i): (" << sub_dom.getLow(0) << "; " << sub_dom.getLow(1) << "); (" << sub_dom.getHigh(0) << "; " << sub_dom.getHigh(1) << ")" << std::endl;

					//std::cout << "gdb_ext_global.get(j): (" << sub_dom_new.getLow(0) << "; " << sub_dom_new.getLow(1) << "); (" << sub_dom_new.getHigh(0) << "; " << sub_dom_new.getHigh(1) << ")" << std::endl;

					//std::cout << "Inte_box_local: (" << inte_box_local.getLow(0) << "; " << inte_box_local.getLow(1) << "); (" << inte_box_local.getHigh(0) << "; " << inte_box_local.getHigh(1) << ")" << std::endl;

					// Grid corresponding for gdb_ext_old.get(i) box
					device_grid & gr = loc_grid_old.get(i);

					//std::cout << "loc_grid_old.get(i): (" << gr.getGrid().getBox().getLow(0) << "; " << gr.getGrid().getBox().getLow(1) << "); (" << gr.getGrid().getBox().getHigh(0) << "; " << gr.getGrid().getBox().getHigh(1) << ")" << std::endl;

					//for (size_t l = 0; l < dim; l++)
						//std::cout << "loc_grid_old.get(i).size on " << l << " dimension: " << gr.getGrid().size(l) << std::endl;
					// Size of the grid to send
					size_t sz[dim];
					for (size_t l = 0; l < dim; l++)
					{
						sz[l] = inte_box_local.getHigh(l) - inte_box_local.getLow(l) + 1;
						//std::cout << "GR_send size on " << l << " dimension: " << sz[l] << std::endl;
					}

					// Grid to send
					device_grid gr_send(sz);
					gr_send.setMemory();

					// Sub iterator across intersection box inside local grid
					grid_key_dx<dim> start = inte_box_local.getKP1();
					grid_key_dx<dim> stop = inte_box_local.getKP2();

					Point<dim,St> p1;
					for (size_t n = 0; n < dim; n++)
						p1.get(n) = gr_send.getGrid().getBox().getLow(n);

					//std::cout << "Grid send P1: (" << p1.get(0) << "; " << p1.get(1) << ")" << std::endl;

					Point<dim,St> p2;
					for (size_t n = 0; n < dim; n++)
						p2.get(n) = gr_send.getGrid().getBox().getHigh(n);

					//std::cout << "Grid send P2: (" << p2.get(0) << "; " << p2.get(1) << ")" << std::endl;
/*
					Point<dim,St> p3;
					for (size_t n = 0; n < dim; n++)
						p3.get(n) = gr.getGrid().getBox().getLow(n);

					std::cout << "Grid local P1: (" << p3.get(0) << "; " << p3.get(1) << ")" << std::endl;

					Point<dim,St> p4;
					for (size_t n = 0; n < dim; n++)
						p4.get(n) = gr.getGrid().getBox().getHigh(n);

					std::cout << "Grid local P2: (" << p4.get(0) << "; " << p4.get(1) << ")" << std::endl;

*/
					std::string start2 = start.to_string();
					std::string stop2 = stop.to_string();

					//std::cout << "Start: " << start2 << "; Stop: " << stop2 << std::endl;

					auto it = gr.getSubIterator(start,stop);

					// Copy selected elements into a new sub-grid
					while (it.isNext())
					{
						auto key = it.get();
						grid_key_dx<dim> key2 = key - start;
						std::string str = key.to_string();

						//std::cout << "Key: " << str << std::endl;
						gr_send.get_o(key2) = gr.get_o(key);
/*
						////////// DEBUG ///////////////
						if (gr.template get<0>(key) == 1)
						{
							count++;
						}
						else if (gr_send.template get<0>(key2) != 1)
						{
							std::cout << "AHHHHHHHHHH????????" << std::endl;
						}
*/
						////////////////

						//gr_send.set(key,gr,key);

						++it;
					}

					aggregate<device_grid,SpaceBox<dim,long int>> aggr;

					aggr.template get<0>() = gr_send;
					aggr.template get<1>() = inte_box;

					// Add to the labeling vector
					lbl_b.get(p_id).add(aggr);
				}
			}
		}
		//std::cout << "Count for points: " << count << std::endl;
		//std::cout << "Count for inte_boxes: " << count2 << std::endl;
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
	void map_(Decomposition & dec, CellDecomposer_sm<dim,St,shift<dim,St>> & cd_sm, openfpm::vector<device_grid> & loc_grid, openfpm::vector<device_grid> & loc_grid_old, openfpm::vector<GBoxes<device_grid::dims>> & gdb_ext, openfpm::vector<GBoxes<device_grid::dims>> & gdb_ext_old, openfpm::vector<GBoxes<device_grid::dims>> & gdb_ext_global)
	{
		// Processor communication size
		openfpm::vector<size_t> prc_sz(v_cl.getProcessingUnits());

		// Contains the processor id of each box (basically where they have to go)
		labelIntersectionGridsProcessor(dec,cd_sm,loc_grid_old,gdb_ext,gdb_ext_old,gdb_ext_global,m_oGrid,prc_sz);
/*
		for (size_t i = 0; i < m_oGrid.size(); i++)
		{
			for (size_t k = 0; k < m_oGrid.get(i).size(); k++)
			{
				device_grid g = m_oGrid.get(i).template get<0>(k);

				auto it = g.getIterator();

				while (it.isNext())
				{
					auto key = it.get();

					if (g.template get<0>(key) != 1)
						std::cout << "WROOOOOOONG" << std::endl;

					++it;
				}
			}
		}
*/

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
/*
		for (size_t i = 0; i < m_oGrid.size(); i++)
		{
			if(m_oGrid.get(i).size() == 0)
				m_oGrid.remove(i);
		}
*/

		decltype(m_oGrid) m_oGrid_new;
		for (size_t i = 0; i < v_cl.getProcessingUnits(); i++)
		{
			if (m_oGrid.get(i).size() != 0)
				m_oGrid_new.add(m_oGrid.get(i));
		}

		// Vector for receiving of intersection grids
		openfpm::vector<openfpm::vector<aggregate<device_grid,SpaceBox<dim,long int>>>> m_oGrid_recv;

		//std::cout << "vcl.getProcessUnitID(): " << v_cl.getProcessUnitID() << "; prc_r.size(): " << prc_r.size() << std::endl;

		//std::cout << "vcl.getProcessUnitID(): " << v_cl.getProcessUnitID() << "; m_oGrid_new.size(): " << m_oGrid_new.size() << std::endl;
/*
		for (size_t i = 0; i < m_oGrid.size(); i++)
		{
			std::cout << "Processor ID:" << v_cl.getProcessUnitID() << "; I: " << i << ", Size: " << m_oGrid.get(i).size() << std::endl;
		}
*/
/*
		for (size_t i = 0; i < m_oGrid_new.size(); i++)
		{
			for (size_t k = 0; k < m_oGrid_new.get(i).size(); k++)
			{
				device_grid g = m_oGrid_new.get(i).template get<0>(k);

				auto it = g.getIterator();

				while (it.isNext())
				{
					auto key = it.get();

					if (g.template get<0>(key) != 1)
						std::cout << "WRONG BEFORE SENDRCV" << std::endl;

					++it;
				}
			}
		}
*/

		// Send and recieve intersection grids
		v_cl.SSendRecv(m_oGrid_new,m_oGrid_recv,prc_r,prc_recv_map,recv_sz_map);
/*
		for (size_t i = 0; i < m_oGrid_recv.size(); i++)
		{
			for (size_t k = 0; k < m_oGrid_recv.get(i).size(); k++)
			{
				device_grid g = m_oGrid_recv.get(i).template get<0>(k);

				auto it = g.getIterator();

				while (it.isNext())
				{
					auto key = it.get();

					if (g.template get<0>(key) != 1)
						std::cout << "WRONG AFTER SENDRCV" << std::endl;

					++it;
				}
			}
		}
*/
/*
		std::cout << "vcl.getProcessUnitID(): " << v_cl.getProcessUnitID() << "; m_oGrid_recv.size(): " << m_oGrid_recv.size() << std::endl;

		for (size_t i = 0; i < m_oGrid_recv.size(); i++)
		{
			std::cout << "Processor ID:" << v_cl.getProcessUnitID() << "; I_recv: " << i << ", Size: " << m_oGrid_recv.get(i).size() << std::endl;
		}

		for (size_t i = 0; i < prc_r.size(); i++)
			std::cout << "vcl.getProcessUnitID(): " << v_cl.getProcessUnitID() << "; prc_r: " << prc_r.get(i) << std::endl;

		for (size_t i = 0; i < prc_recv_map.size(); i++)
			std::cout << "vcl.getProcessUnitID(): " << v_cl.getProcessUnitID() << "; prc_recv_map: " << prc_recv_map.get(i) << std::endl;

		for (size_t i = 0; i < recv_sz_map.size(); i++)
			std::cout << "vcl.getProcessUnitID(): " << v_cl.getProcessUnitID() << "; recv_sz_map: " << recv_sz_map.get(i) << std::endl;
*/
		// Reconstruct the new local grids
		grids_reconstruct(m_oGrid_recv,loc_grid,gdb_ext,cd_sm);
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
