/*
 * CartesianGraphFactory.hpp
 *
 *  Created on: Nov 28, 2014
 *      Author: i-bird
 */

#ifndef CARTESIANGRAPHFACTORY_HPP_
#define CARTESIANGRAPHFACTORY_HPP_

#include "map_vector.hpp"
#include "grid.hpp"
#include "Space/Shape/Box.hpp"

/*! \brief This class construct a cartesian graph
 *
 * This class construct a cartesian graph
 *
 * \param dim dimensionality of the cartesia ngrid
 *
 */

template<unsigned int dim, typename Graph>
class CartesianGraphFactory
{

public:

	/*!
	 *
	 * \brief Construct a cartesian graph, with V and E edge properties
	 *
	 * Construct a cartesian graph, with V and E edge properties
	 *
	 * Each vertex is a subspace (Hyper-cube) of dimension dim, each vertex is
	 * connected with an edge if two vertex (Hyper-cube) share a element of dimension grater than
	 * dim_c
	 *
	 * \param sz Vector that store the size of the grid on each dimension
	 * \param dom Box enclosing the physical domain
	 *
	 * \tparam se Indicate which properties fill with the element weight. The
	 *           element weight is the point, line , surface, d dimensional object
	 *           in contact (in common between two hyper-cube)
	 * \tparam T type of the domain like (int real complex ... )
	 * \tparam dim_c Connectivity dimension
	 * \tparam Memory class that create new memory
	 *
	 */
	template <unsigned int se,typename T, unsigned int dim_c>
	static Graph construct(std::vector<size_t> sz, Box<dim,T> dom)
	{
#ifdef DEBUG
		//! The size is wrong signal it

		if (sz.size() != dim)
		{std::cerr << "Error this factory has been specialized for catesian grid of size " << dim << "\n";}
#endif

		// Calculate the size of the hyper-cubes on each dimension

		T szd[dim];

		for (int i = 0 ; i < dim ; i++)
		{szd[i] = (dom.getHigh(i) - dom.getLow(i)) / sz[i];}

		//! Construct an hyper-cube of dimension dim

		HyperCube<dim> hc;

		// Construct a grid info

		grid<dim,void> g(sz);

		// Create a graph with the number of vertices equal to the number of
		// grid point

		//! Graph to construct

		Graph gp(g.size());

		/******************
		 *
		 * Create the edges
		 *
		 ******************/

		//! Construct a key iterator

		grid_key_dx_iterator<dim> k_it(g);

		//! Iterate through all the elements

		while (k_it.isNext())
		{
			grid_key_dx<dim> key = k_it.get();

			// Get the combinations of dimension d

			for (int d = dim-1 ; d >= dim_c ; d--)
			{
				// create the edges for that dimension

				std::vector<comb<dim>> c = hc.getCombinations_R(d);

				// for each combination calculate a safe linearization and create an edge

				for (int j = 0 ; j < c.size() ; j++)
				{
					// Calculate the element size

					T ele_sz = 0;

					// for each dimension multiply and reduce

					for (int s = 0 ; s < dim ; s++)
					{
						ele_sz += szd[s] * abs(c[j][s]);
					}

					// Calculate the end point vertex id
					// Calculate the start point id

					size_t start_v = g.LinId(key);
					size_t end_v = g.template LinId<CheckExistence>(key,c[j].getComb());

					// Add an edge and set the se edge property to the size of the face (communication weight)
					gp.template addEdge<CheckExistence>(start_v,end_v).template get<se>() = ele_sz;
				}
			}

			++k_it;
		}

		return gp;

	}

};

#endif /* CARTESIANGRAPHFACTORY_HPP_ */
