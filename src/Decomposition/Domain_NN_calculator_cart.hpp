/*
 * Domain_NN_calculator.hpp
 *
 *  Created on: Nov 26, 2016
 *      Author: i-bird
 */

#ifndef SRC_DECOMPOSITION_DOMAIN_NN_CALCULATOR_CART_HPP_
#define SRC_DECOMPOSITION_DOMAIN_NN_CALCULATOR_CART_HPP_

#include <Vector/map_vector.hpp>
#include "NN/CellList/ParticleItCRS_Cells.hpp"

/*! \brief This class calculate processor domains and neighborhood
 *  of each processor domain
 *
 * \param dim Dimensionality
 *
 */
template<unsigned int dim>
class domain_nn_calculator_cart
{
	//! True if domain and anomalous domain cells are computed
	bool are_domain_anom_computed;

    /////////////////////////////// CRS //////////////////////

	//! anomalous cell neighborhood for CRS
	openfpm::vector<subsub<dim>> anom;

	//! Set of anomalous CRS domain cells linearized
	openfpm::vector<subsub_lin<dim>> anom_lin;

	//! Set of normal domain cells for CRS
	openfpm::vector<grid_key_dx<dim>> dom;

	//! Set of normal CRS domain cells linearized
	openfpm::vector<size_t> dom_lin;

	////////////////////////////////// DOMAIN CELLS ///////////////////////

	//! Set of domain cells
	openfpm::vector<grid_key_dx<dim>> dom_cells;

	//! Set of linearized domain cells
	openfpm::vector<size_t> dom_cells_lin;

	//////////////////////////////////////////////////////////////

	//! Processor box
	Box<dim,long int> proc_box;

	//! Processor cells-grid
	grid_sm<dim,void> gs;

	//! key with all coordinates set to one
	grid_key_dx<dim> one;

	/*! \brief Calculate the subdomain that are in the skin part of the domain
	 *
       \verbatim

		+---+---+---+---+---+---+
		| 1 | 2 | 3 | 4 | 5 | 6 |
		+---+---+---+---+---+---+
		|27 |               | 7 |
		+---+               +---+
		|26 |               | 8 |
		+---+               +---+
		|25 |               | 9 |
		+---+   DOM1        +---+
		|24 |               |10 |
		+---+               +---+
		|23 |               |11 |
		+---+               +---+
		|22 |               |12 |
		+---+-----------+---+---+
		|21 |           |13 |
		+---+           +---+
		|20 |   DOM2    |14 |
		+---+---+---+---+---+
		|19 |18 |17 | 16|15 |
		+---+---+---+---+---+    <----- Domain end here
                                      |
                        ^             |
                        |_____________|


       \endverbatim
	 *
	 * In some cases like symmetric with CRS Scheme The cells indicated with numbers has a non straigh-forward
	 * neighborhood. This function calculate the list of the domain cells with normal symmetric neighborhood type
	 * and compute a list of cells with more complex neighboring cells.
	 *
	 * \param sub_keys array that contain the position of the sub-sub-domains indicated with numbers
	 *        in grid coordinates + for each its neighboring cells
	 *
	 * \param dom_cells list of all the domain cells
	 *
	 * \param dom_subsub cells with normal neighborhood
	 *
	 * \param loc_box array of local sub-sub-domain in grid coordinates
	 *
	 * \param proc_box Processor bounding box in local coordinates
	 *
	 *
	 */
	void CalculateDomAndAnomCells(openfpm::vector<subsub<dim>> & sub_keys,
			                      openfpm::vector<grid_key_dx<dim>> & dom_subsub,
								  openfpm::vector<grid_key_dx<dim>> & dom_cells,
								  const ::Box<dim,long int> & proc_box,
								  const openfpm::vector<::Box<dim, size_t>> & loc_box)
	{
		// Reset dom and dom_subsub
		sub_keys.clear();
		dom_subsub.clear();

		size_t sz[dim];

		// ----Grid size = proc_box.getHigh(j) - proc_box.getLow(j)
		// +2 is padding

		for (size_t j = 0 ; j < dim ; j++)
		{sz[j] = proc_box.getHigh(j) - proc_box.getLow(j) + 2 + 1;}

		gs.setDimensions(sz);

		// Set the grid
		grid_cpu<dim, aggregate<openfpm::vector<grid_key_dx<dim>> >> g(sz);
		g.setMemory();

		for (size_t i = 0 ; i < dim ; i++)
		{one.set_d(i,1);}

		// Calculate the csr neighborhood
		openfpm::vector<std::pair<grid_key_dx<dim>,grid_key_dx<dim>>> csr;
		NNcalc_csr(csr);

		// Draw the domain on this grid
		for (size_t i = 0 ; i < loc_box.size() ; i++)
		{
			grid_key_dx<dim> start;
			grid_key_dx<dim> stop;

			for (size_t j = 0 ; j < dim ; j++)
			{
				start.set_d(j,loc_box.template get<0>(i)[j] - proc_box.getLow(j) + 1);
				stop.set_d(j,loc_box.template get<1>(i)[j] - proc_box.getLow(j) + 1);
			}

			grid_key_dx_iterator_sub<dim> sub(g.getGrid(),start,stop);

			while (sub.isNext())
			{
				auto key = sub.get();

				for (size_t j = 0 ; j < csr.size() ; j++)
				{
					grid_key_dx<dim> src = key + csr.get(j).first;
					grid_key_dx<dim> dst = key + csr.get(j).second;
					g.template get<0>(src).add(dst);
				}

				dom_cells.add(key - one);

				++sub;
			}
		}

		// Span all the grid point and take all the sub-sub-domains that has a
		// neighborhood non-consistent to the normal symmetric NN, and the one that
		// are consistent to NN symmetric
		grid_key_dx_iterator<dim> it(g.getGrid());

		while (it.isNext())
		{
			auto key = it.get();

			// Adding in the list of the non-normal neighborhood cells
			if (g.template get<0>(key).size() == openfpm::math::pow(3,dim)/2+1)
			{
				// Add in the list of the normal neighborhood list
				dom_subsub.add(key - one);
			}
			else if (g.template get<0>(key).size() != 0)
			{
				sub_keys.add();
				sub_keys.last().subsub = key - one;
				// Adding the neighborhood of the cell

				sub_keys.last().NN_subsub.resize(g.template get<0>(key).size());

				for (size_t i = 0 ; i < g.template get<0>(key).size() ; i++)
				{sub_keys.last().NN_subsub.get(i) = g.template get<0>(key).get(i) - one;}
			}

			++it;
		}
	}

	/*! \brief Linearize the sub-sub-domains ids
	 *
	 * A subsub domain can be identified by a set of number (i,j).
	 * The linearization transform it into a number
	 *
	 * \param anom set of grid keys to linearize
	 * \param anom_lin linearized output
	 * \param shift shifting to add for the linearizaton
	 * \param gs information about the grid to linearize
	 *
	 */
	void linearize_subsub(const openfpm::vector<subsub<dim>> & anom,
			              openfpm::vector<subsub_lin<dim>> & anom_lin,
						  const grid_key_dx<dim> & shift,
						  const grid_sm<dim,void> & gs)
	{
		anom_lin.clear();
		for (size_t i = 0 ; i < anom.size() ; i++)
		{
			anom_lin.add();
			anom_lin.last().subsub = gs.LinId(anom.get(i).subsub + shift);

			long int self_cell = -1;

			for (size_t j = 0 ; j < anom.get(i).NN_subsub.size() ; j++)
			{
				anom_lin.get(i).NN_subsub.add((long int)gs.LinId(anom.get(i).NN_subsub.get(j) + shift) - anom_lin.get(i).subsub);

				// This indicate that for example in the neighborhood of one cell it-self is included in the list
				// For example the cell 100 is in the neighborhood of the cell 100
				if (anom_lin.get(i).NN_subsub.last() == 0)
					self_cell = anom_lin.get(i).NN_subsub.size() - 1;
			}

			// if exist the self interacting cell (Example cell 100 neighborhood of cell 100), this cell MUST BE ALWAYS at the beginning
			if (self_cell != -1)
			{
				// bring the self-cell into the beginning
				size_t tmp = anom_lin.get(i).NN_subsub.get(0);
				anom_lin.get(i).NN_subsub.get(0) = 0;
				anom_lin.get(i).NN_subsub.get(self_cell) = tmp;
			}
		}
	}

public:

	domain_nn_calculator_cart()
	:are_domain_anom_computed(false)
	{
	}

	/*! \brief Set parameters to calculate the cell neighborhood
	 *
	 * \param proc_box processor cells box
	 *
	 */
	void setParameters(const Box<dim,long int> & proc_box)
	{
		this->proc_box = proc_box;
	}

	/*! \brief Set parameters to calculate the cell neighborhood
	 *
	 * \param loc_box set of local sub-domains
	 * \param shift to apply in the linearization
	 * \param gs grid of cells (for the processor domain)
	 *
	 */
	void setNNParameters(openfpm::vector<::Box<dim, size_t>> & loc_box,
			             const grid_key_dx<dim> & shift,
						 const grid_sm<dim,void> & gs)
	{
		if (are_domain_anom_computed == false)
		{
			CalculateDomAndAnomCells(anom,dom,dom_cells,proc_box,loc_box);
			are_domain_anom_computed = true;

			dom_cells_lin.clear();
			for (size_t i = 0 ; i < dom_cells.size() ; i++)
				dom_cells_lin.add(gs.LinId(dom_cells.get(i) + shift));


			dom_lin.clear();
			for (size_t i = 0 ; i < dom.size() ; i++)
				dom_lin.add(gs.LinId(dom.get(i) + shift));

			linearize_subsub(anom,anom_lin,shift,gs);
		}
	}

	/*! \brief Get the domain Cells
	 *
	 *
	 * \return The set of domain cells
	 *
	 */
	openfpm::vector<size_t> & getDomainCells()
	{
		return dom_cells_lin;
	}

	/*! \brief Get the domain Cells
	 *
	 *
	 * \return The set of domain cells
	 *
	 */
	openfpm::vector<size_t> & getCRSDomainCells()
	{
		return dom_lin;
	}

	/*! \brief Get the domain anomalous cells
	 *
	 *
	 * \return The set of anomalous cells
	 *
	 */
	openfpm::vector<subsub_lin<dim>> & getCRSAnomDomainCells()
	{
		return anom_lin;
	}


	/*! \brief In case you have to recompute the indexes
	 *
	 *
	 */
	void reset()
	{
		are_domain_anom_computed = false;
	}
};

#endif /* SRC_DECOMPOSITION_DOMAIN_NN_CALCULATOR_CART_HPP_ */
