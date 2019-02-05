/*
 * CartDecomposition_ext.hpp
 *
 *  Created on: Mar 6, 2016
 *      Author: i-bird
 */

#ifndef SRC_DECOMPOSITION_CARTDECOMPOSITION_EXT_HPP_
#define SRC_DECOMPOSITION_CARTDECOMPOSITION_EXT_HPP_

#include "memory/HeapMemory.hpp"
#include "Decomposition/Distribution/ParMetisDistribution.hpp"
#include "Space/Ghost.hpp"
#include "Decomposition/nn_processor.hpp"

template<unsigned int dim, typename T, typename Memory = HeapMemory, typename Distribution = ParMetisDistribution<dim, T>>
class CartDecomposition;

/**
 * \brief This class decompose a space into sub-sub-domains and distribute them across processors
 *
 * \tparam dim is the dimensionality of the physical domain we are going to decompose.
 * \tparam T type of the space we decompose, Real, Integer, Complex ...
 * \tparam Memory Memory factory used to allocate memory
 * \tparam Distribution type of distribution, can be ParMetisDistribution or MetisDistribution
 *
 * Most of the functionality is the same as CartDecomposition so refer to that class for more information
 *
 * The additional functionality is the possibility to produce an extended decomposition, in figure is
 * show what we mean with extended
 *
 * \see CartDecomposition
 *
 *
 *
 * ### Create a Cartesian decomposition object on a Box space, distribute, calculate internal and external ghost boxes
 * \snippet CartDecomposition_unit_test.hpp Create CartDecomposition
 *
 */

template<unsigned int dim, typename T, typename Memory = HeapMemory, typename Distribution = ParMetisDistribution<dim, T>>
class CartDecomposition_ext: public CartDecomposition<dim,T,Memory,Distribution>
{
private:

	/*! \brief It copy the sub-domains into another CartesianDecomposition object extending them
	 *
	 * \see duplicate (in case of extended domain)
	 *
	 * \param dec Cartesian decomposition object
	 * \param ext_dom Extended domain
	 *
	 */
	void extend_subdomains(const CartDecomposition<dim,T,Memory,Distribution> & dec, const ::Box<dim,T> & ext_dom)
	{
		// Box
		typedef ::Box<dim,T> b;

		this->bbox.zero();

		// Extend sub-domains
		for (size_t i = 0 ; i < dec.sub_domains.size() ; i++)
		{
			::Box<dim,T> box;

			// Calculate the extended box
			for (size_t j = 0 ; j < dim ; j++)
			{
				if (dec.sub_domains.template get<b::p1>(i)[j] == dec.domain.getLow(j))
					box.setLow(j,ext_dom.getLow(j));
				else
					box.setLow(j,dec.sub_domains.template get<b::p1>(i)[j]);

				if (dec.sub_domains.template get<b::p2>(i)[j] == dec.domain.getHigh(j))
					box.setHigh(j,ext_dom.getHigh(j));
				else
					box.setHigh(j,dec.sub_domains.template get<b::p2>(i)[j]);
			}

			// add the subdomain
			this->sub_domains.add(box);

			// Calculate the bound box
			this->bbox.enclose(box);
		}
	}

	/*! \brief Extend the fines for the new Cartesian decomposition
	 *
	 * \param dec Non-extended decomposition
	 *
	 */
/*	void extend_fines(const CartDecomposition<dim,T,Memory,Distribution> & dec)
	{
		// Extension, first we calculate the extensions of the new domain compared
		// to the old one in cell units (each cell unit is a sub-sub-domain)
		::Box<dim,size_t> ext;
		// Extension of the new fines structure
		::Box<dim,size_t> n_fines_ext;
		// Extension of the old fines structure
		::Box<dim,size_t> o_fines_ext;

		size_t sz_new[dim];
		size_t sz_old[dim];

		for (size_t i = 0; i < dim ; i++)
		{
			size_t p1 = (dec.domain.getLow(i) - dec.domain.getLow(i)) / dec.cd.getCellBox().getHigh(i) + 1;
			size_t p2 = (dec.domain.getLow(i) - dec.domain.getLow(i)) / dec.cd.getCellBox().getHigh(i) + 1;

			ext.setLow(i,p1);
			ext.setHigh(i,p2);
			sz_new[i] = p1+p2+dec.cd.getGrid().size(i);
			sz_old[i] = dec.cd.getGrid().size(i);
		}

		grid_sm<dim,void> info_new(sz_new);
		grid_sm<dim,void> info_old(sz_old);

		// resize the new fines
		this->fine_s.resize(info_new.size());

		// we create an iterator that iterate across the full new fines
		grid_key_dx_iterator<dim> fines_t(info_new);

		while (fines_t.isNext())
		{
			auto key = fines_t.get();

			// new_fines is bigger than old_fines structure
			// out of bound key must be adjusted
			// The adjustment produce a natural extension
			// a representation can be seen in the figure of
			// CartDecomposition duplicate function with extended domains

			grid_key_dx<dim> key_old;
			for (size_t i = 0 ; i < dim ; i++)
			{
				key_old.set_d(i,(long int)key.get(i) - ext.getLow(i));
				if (key_old.get(i) < 0)
					key_old.set_d(i,0);
				else if(key_old.get(i) >= (long int)info_old.size(i) )
					key_old.set_d(i,info_old.size(i)-1);
			}

			this->fine_s.get(info_new.LinId(key)) = dec.fine_s.get(info_old.LinId(key_old));

			++fines_t;
		}

		this->gr.setDimensions(sz_new);

		// the new extended CellDecomposer must be consistent with the old cellDecomposer.
		this->cd.setDimensions(dec.cd,ext);
	}*/

	void reconstruct_fine_s_from_extended_domain(const ::Box<dim,T> & ext_domain)
	{
		this->initialize_fine_s(ext_domain);
		this->construct_fine_s();
	}

public:

	/*! \brief Cartesian decomposition constructor
	 *
	 * \param v_cl VCluster
	 *
	 */
	CartDecomposition_ext(Vcluster & v_cl)
	:CartDecomposition<dim,T,Memory,Distribution>(v_cl)
	{
	}

	//! The non-extended decomposition base class
	typedef CartDecomposition<dim,T,Memory,Distribution> base_type;

	/*! \brief It create another object that contain the same decomposition information but with different ghost boxes and an extended domain
	 *
	 * The domain extension is produced extending the boxes at the border like in figure
	 *
	 * \verbatim
	 *
	+--------------^--------^----------^----------+
	|              |        |          |          |
	|        A     |    E   |     F    |    N     |
	|    +-----------------------------------+---->
	|    |         |        |          |     |    |
	|  A |   A     |        |     F    |     |    |
	|    |         |        |          |     |    |
	|    |         |    E   +----------+  N  |  N |
	<--------------+        |          |     |    |
	|    |         |        |          |     |    |
	|    |         |        |     G    |     |    |
	|    |         |        |          +---------->
	|  B |   B     |        +----------+     |    |
	|    |         +--------+          |  M  |  M |
	|    |         |        |     H    |     |    |
	|    |         |        +-----+----+---------->
	<--------------+    D   |     |          |    |
	|    |         |        |  I  |     L    |  L |
	|  C |   C     |        |     |          |    |
	|    |         |        |     |          |    |
	|    +-----------------------------------+    |
	|              |        |     |               |
	|        C     |    D   |  I  |     L         |
	+--------------v--------v-----v---------------+

	 *
	 * \endverbatim
	 *
	 * \param dec Decomposition
	 * \param g ghost
	 * \param ext_domain extended domain (MUST be extended)
	 *
	 * \return a duplicated decomposition with different ghost boxes and an extended domain
	 *
	 */
	void setParameters(const CartDecomposition<dim,T,Memory,Distribution> & dec, const Ghost<dim,T> & g, const ::Box<dim,T> & ext_domain)
	{
		this->box_nn_processor = dec.box_nn_processor;

		// Calculate new sub-domains for extended domain
		extend_subdomains(dec,ext_domain);

		// Calculate fine_s structure for the extended domain
		// update the cell decomposer and gr
//		extend_fines(dec);
		reconstruct_fine_s_from_extended_domain(ext_domain);

		// Get the old sub-sub-domain grid extension

		this->domain = ext_domain;

		// spacing does not change

		for (size_t i = 0 ; i < dim ; i++)
		{this->spacing[i] = dec.spacing[i];};

		this->ghost = g;
		this->dist = dec.dist;

		for (size_t i = 0 ; i < dim ; i++)
			this->bc[i] = dec.bc[i];

		(static_cast<nn_prcs<dim,T> &>(*this)).create(this->box_nn_processor, this->sub_domains);
		(static_cast<nn_prcs<dim,T> &>(*this)).applyBC(ext_domain,g,this->bc);

		this->Initialize_geo_cell_lists();
		this->calculateGhostBoxes();
	}

};



#endif /* SRC_DECOMPOSITION_CARTDECOMPOSITION_EXT_HPP_ */
