/*
 * DistGraphFactory.hpp
 *
 *  Created on: Dec 03, 2015
 *      Author: Antonio Leo
 */

#ifndef DISTGRAPHFACTORYOLD_HPP_
#define DISTGRAPHFACTORYOLD_HPP_

#include "VCluster.hpp"
#include "Vector/map_vector.hpp"
#include "Graph/map_graph.hpp"
#include "Grid/grid_sm.hpp"
#include "Space/Shape/Box.hpp"
#include "Space/Shape/HyperCube.hpp"
#include "parmetis.h"

/*! \brief This class work as a functor
 *
 * For each number in the boost::mpl::vector (for example 3 6) set the properties of the vertex at the
 * specified id (3 6) with pos[d] * spacing[d] with d running from 0 to 1, pos[d] the position id of the vertex
 * spacing the grid spacing
 *
 * Example
 *
 * if we give a grid_key of dimension 2 4x4 the expression "pos[d] * spacing[d]"
 * will assume the value
 *
 * (0.0 0.0) (0.25 0.0) ...... (1.0 0.0)
 * (0.0 0.25)................. (1.0 0.25)
 * ....................................
 * (0.0 1.0).................. (1.0 1.0)
 *
 * and the properties 3 6 will be filled with the numbers 0.0 0.0    .......  1.0 1.0
 * progressively
 *
 * \tparam dim Dimensionality of the cartesian grid
 * \tparam dT type of the domain
 * \tparam G_v vertex type object
 * \tparam v boost::mpl::vector containing all the index to fill
 * \tparam is_stub when is true, produce a trivial operator(),
 *         to use when v is an empty vector to avoid compilation error
 *
 */

template<unsigned int dim, typename dT, typename G_v, typename v, int impl>
class fill_prop_v
{
	//! Reference to an array containing the spacing
	const dT (&szd)[dim];

	//! grid_key_dx Reference containing the actual position
	grid_key_dx<dim> & gk;

	//! Vertex object to fill
	G_v & g_v;

	//! grid info
	const grid_sm<dim, void> & gs;

public:

	//! Fill the object from where to take the properties
	fill_prop_v(G_v & g_v, const dT (&szd)[dim], grid_key_dx<dim> & gk, const grid_sm<dim, void> & gs) :
			szd(szd), gk(gk), g_v(g_v), gs(gs)
	{
	}

	//! It call the function for each property we want to copy
	template<typename T>
	void operator()(T& t) const
	{
		typedef typename boost::fusion::result_of::at<v, boost::mpl::int_<T::value>>::type t_val;

		g_v.template get<t_val::value>() = gk.get(T::value) * szd[T::value];
	}
};


/*! \brief This class work as a functor
 *
 * For each number in the boost::mpl::vector (for example 3 6) set the properties of the vertex at the
 * specified id (3 6) with pos[d] * spacing[d] with d running from 0 to 1, pos[d] the position id of the vertex
 * spacing the grid spacing
 *
 * Example
 *
 * if we give a grid_key of dimension 2 4x4 the expression "pos[d] * spacing[d]"
 * will assume the value
 *
 * (0.0 0.0) (0.25 0.0) ...... (1.0 0.0)
 * (0.0 0.25)................. (1.0 0.25)
 * ....................................
 * (0.0 1.0).................. (1.0 1.0)
 *
 * and the properties 3 6 will be filled with the numbers 0.0 0.0    .......  1.0 1.0
 * progressively
 *
 * \tparam dim Dimensionality of the cartesian grid
 * \tparam dT type of the domain
 * \tparam G_v vertex type object
 * \tparam v boost::mpl::vector containing all the index to fill
 *
 */

template<unsigned int dim, typename dT, typename G_v, typename v>
class fill_prop_v<dim, dT, G_v, v, 0>
{

public:

	//! Fill the object from where to take the properties
	fill_prop_v(G_v & g_v, const dT (&szd)[dim], grid_key_dx<dim> & gk, const grid_sm<dim, void> & gs)
	{
	}

	//! It call the function for each property we want to copy
	template<typename T>
	void operator()(T& t) const
	{
	}
};

/*! \brief This class work as a functor
 *
 * For each number in the boost::mpl::vector (for example 3 6) set the properties of the vertex at the
 * specified id (3 6) with pos[d] * spacing[d] with d running from 0 to 1, pos[d] the position id of the vertex
 * spacing the grid spacing
 *
 * Example
 *
 * if we give a grid_key of dimension 2 4x4 the expression "pos[d] * spacing[d]"
 * will assume the value
 *
 * (0.0 0.0) (0.25 0.0) ...... (1.0 0.0)
 * (0.0 0.25)................. (1.0 0.25)
 * ....................................
 * (0.0 1.0).................. (1.0 1.0)
 *
 * and the properties 3 6 will be filled with the numbers 0.0 0.0    .......  1.0 1.0
 * progressively
 *
 * \tparam dim Dimensionality of the cartesian grid
 * \tparam dT type of the domain
 * \tparam G_v vertex type object
 * \tparam v boost::mpl::vector containing all the index to fill
 *
 */

template<unsigned int dim, typename dT, typename G_v, typename v>
class fill_prop_v<dim, dT, G_v, v, 2>
{

	//! Reference to an array containing the spacing
	const dT (&szd)[dim];

	//! grid_key_dx Reference containing the actual position
	grid_key_dx<dim> & gk;

	//! Vertex object to fill
	G_v & g_v;

	//! grid info
	const grid_sm<dim, void> & gs;

public:

	//! Fill the object from where to take the properties
	fill_prop_v(G_v & g_v, const dT (&szd)[dim], grid_key_dx<dim> & gk, const grid_sm<dim, void> & gs)
	:szd(szd), gk(gk), g_v(g_v), gs(gs)
	{
	}

	//! It call the function for each property we want to copy
	template<typename T>
	void operator()(T& t) const
	{
		typedef typename boost::fusion::result_of::at<v, boost::mpl::int_<0>>::type t_val;
		typedef typename boost::mpl::at<typename G_v::T_type::type,t_val>::type s_type;

		for (size_t i = 0 ; i < std::extent<s_type>::value ; i++)
			g_v.template get<t_val::value>()[i] = 0.0;

		for (size_t i = 0 ; i < dim ; i++)
			g_v.template get<t_val::value>()[i] = gk.get(i) * static_cast<float>(szd[i]);
	}
};

/*! \brief Operator for vector and scalar property
 *
 * \tparam i Size of the property
 * \tparam p Type of the property
 * \tparam Graph Graph
 * \tparam pos Array of properties
 */
template<int i, typename p, typename Graph, int ... pos>
struct fill_prop_v_by_type
{

	typedef typename boost::mpl::at<p, boost::mpl::int_<0>>::type v_element;
	typedef typename boost::mpl::at<typename Graph::V_type::type, v_element>::type pos_prop_type;

	enum
	{
		value = ((sizeof...(pos) != 0) * (std::is_array<pos_prop_type>::value + 1))
	};

};

/*! \brief Operator for vector and scalar property in the case there are no properties
 *
 * \tparam i Size of the property
 * \tparam p Type of the property
 * \tparam Graph Graph
 * \tparam pos Array of properties
 */
template<typename p, typename Graph, int ... pos>
struct fill_prop_v_by_type<0, p, Graph, pos...>
{
	enum
	{
		value = 0
	};

};

/*! \brief Graph constructor function specialization
 *
 * On C++ partial function specialization is not allowed, so we need a class to do it
 *
 * \see CartesianGraphFactory method construct
 *
 */

template<unsigned int dim, typename Graph, int se, typename T, unsigned int dim_c, int ... pos>
class DistGraph_constr_impl
{
public:
	//! Construct Cartesian graph
	static Graph construct(const size_t (&sz)[dim], Box<dim, T> dom)
	{
		Vcluster &v_cl = *global_v_cluster;

		// Calculate the size of the hyper-cubes on each dimension
		T szd[dim];

		for (size_t i = 0; i < dim; i++)
		{
			szd[i] = (dom.getHigh(i) - dom.getLow(i)) / sz[i];
		}

		//! Construct an hyper-cube of dimension dim

		HyperCube<dim> hc;

		// Construct a grid info

		grid_sm<dim, void> g(sz);

		//! Get the processor id
		size_t p_id = v_cl.getProcessUnitID();

		//! Get the number of processing units
		size_t Np = v_cl.getProcessingUnits();

		//! Division of vertices in Np graphs
		//! Put (div+1) vertices in mod graphs
		//! Put div vertices in the rest of the graphs
		size_t mod_v = g.size() % Np;
		size_t div_v = g.size() / Np;

		//! Distribution vector
		openfpm::vector<idx_t> vtxdist(v_cl.getProcessingUnits() + 1);

		for (int i = 0; i <= Np; i++)
		{
			if (i < mod_v)
				vtxdist.get(i) = (div_v + 1) * (i);
			else
				vtxdist.get(i) = (div_v) * (i) + mod_v;
		}

		//! Get size of this processor graph
		size_t gp_size = vtxdist.get(p_id + 1) - vtxdist.get(p_id);

		//! Graph to construct
		Graph gp(gp_size);

		//! Store the decomposition vector inside the distributed graph
		gp.initDistributionVector(vtxdist);

		/******************
		 *
		 * Create the edges and fill spatial
		 * information properties
		 *
		 ******************/

		//! Construct a key iterator
		grid_key_dx_iterator<dim> k_it(g);

		//! Local iterator of the graph
		size_t local_it = 0;

		//! Iterate through all the elements

		while (k_it.isNext())
		{
			size_t v_id = g.LinId(k_it.get());

			if (v_id < vtxdist.get(p_id + 1) && v_id >= vtxdist.get(p_id))
			{

				grid_key_dx<dim> key = k_it.get();

				// Vertex object

				auto obj = gp.vertex(local_it);

				typedef typename to_boost_vmpl<pos...>::type p;

				// vertex spatial properties functor

				fill_prop_v<dim, T, decltype(gp.vertex(local_it)), typename to_boost_vmpl<pos...>::type, fill_prop_v_by_type<sizeof...(pos), p, Graph, pos...>::value> flp(obj, szd, key, g);

				// fill properties

				boost::mpl::for_each<boost::mpl::range_c<int, 0, sizeof...(pos)> >(flp);

				// set map global to local in the graph, needed because vertex is already created without addVertex method

				gp.setGlobalMap(v_id, local_it, v_id);

				// Get the combinations of dimension d

				for (size_t d = dim - 1; d >= dim_c; d--)
				{
					// create the edges for that dimension

					std::vector<comb<dim>> c = hc.getCombinations_R(d);

					// for each combination calculate a safe linearization and create an edge

					for (size_t j = 0; j < c.size(); j++)
					{
						// Calculate the element size

						T ele_sz = 0;

						// for each dimension multiply and reduce

						for (size_t s = 0; s < dim; s++)
						{
							ele_sz += szd[s] * abs(c[j][s]);
						}

						// Calculate the end point vertex id
						// Calculate the start point id

						size_t start_v = local_it;
						size_t end_v = g.template LinId<CheckExistence>(key, c[j].getComb());

						// check if the end_v is valid globally
						if (end_v < g.size())
						{
							// Add an edge and set the the edge property to the size of the face (communication weight)
							gp.template addEdge<NoCheck>(start_v, end_v).template get<se>() = ele_sz;
						}
					}
				}
				++local_it;
			}
			++k_it;
		}

		return gp;
	}
};

/*! \brief Graph constructor function specialization
 *
 * On C++ partial function specialization is not allowed, so we need a class to do it
 * This specialization handle the case when we have NO_EDGE option active
 *
 * \see CartesianGraphFactory method construct
 *
 */

template<unsigned int dim, typename Graph, typename T, unsigned int dim_c, int ... pos>
class DistGraph_constr_impl<dim, Graph, NO_EDGE, T, dim_c, pos...>
{
public:
	//! Construct Cartesian graph
	static Graph construct(const size_t (&sz)[dim], Box<dim, T> dom)
	{
		Vcluster &v_cl = *global_v_cluster;

		// Calculate the size of the hyper-cubes on each dimension

		T szd[dim];

		for (size_t i = 0; i < dim; i++)
		{
			szd[i] = (dom.getHigh(i) - dom.getLow(i)) / sz[i];
		}

		//! Construct an hyper-cube of dimension dim

		HyperCube<dim> hc;

		// Construct a grid info

		grid_sm<dim, void> g(sz);

		//! Get the processor id
		size_t p_id = v_cl.getProcessUnitID();

		//! Get the number of processing units
		size_t Np = v_cl.getProcessingUnits();

		//! Division of vertices in Np graphs
		//! Put (div+1) vertices in mod graphs
		//! Put div vertices in the rest of the graphs
		size_t mod_v = g.size() % Np;
		size_t div_v = g.size() / Np;

		//! Distribution vector
		openfpm::vector<idx_t> vtxdist(v_cl.getProcessingUnits() + 1);

		for (size_t i = 0; i <= Np; i++)
		{
			if (i < mod_v)
				vtxdist.get(i) = (div_v + 1) * (i);
			else
				vtxdist.get(i) = (div_v) * (i) + mod_v;
		}

		size_t gp_size = vtxdist.get(p_id + 1) - vtxdist.get(p_id);

		//! Graph to construct
		Graph gp(gp_size);

		//! Store the decomposition vector inside the distributed graph
		gp.initDistributionVector(vtxdist);

		/******************
		 *
		 * Create the edges and fill spatial
		 * information properties
		 *
		 ******************/

		//! Construct a key iterator
		grid_key_dx_iterator<dim> k_it(g);

		//! Local iterator of the graph
		size_t local_it = 0;

		//! Iterate through all the elements

		while (k_it.isNext())
		{
			size_t v_id = g.LinId(k_it.get());

			if (v_id < (size_t)vtxdist.get(p_id + 1) && v_id >= (size_t)vtxdist.get(p_id))
			{
				grid_key_dx<dim> key = k_it.get();

				// Vertex object
				auto obj = gp.vertex(local_it);

				typedef typename to_boost_vmpl<pos...>::type p;

				// vertex spatial properties functor

				fill_prop_v<dim, T, decltype(gp.vertex(local_it)), typename to_boost_vmpl<pos...>::type, fill_prop_v_by_type<sizeof...(pos), p, Graph, pos...>::value> flp(obj, szd, key, g);

				// fill properties

				boost::mpl::for_each<boost::mpl::range_c<int, 0, sizeof...(pos)> >(flp);

				// set map global to local in the graph, needed because vertex is already created without addVertex method

				gp.setGlobalMap(v_id, local_it, v_id);

				// Get the combinations of dimension d

				for (size_t d = dim - 1; d >= dim_c; d--)
				{
					// create the edges for that dimension

					std::vector<comb<dim>> c = hc.getCombinations_R(d);

					// for each combination calculate a safe linearization and create an edge

					for (size_t j = 0; j < c.size(); j++)
					{
						// Calculate the element size

						T ele_sz = 0;

						// for each dimension multiply and reduce

						for (size_t s = 0; s < dim; s++)
						{
							ele_sz += szd[s] * abs(c[j][s]);
						}

						// Calculate the end point vertex id
						// Calculate the start point id

						size_t start_v = local_it;
						size_t end_v = g.template LinId<CheckExistence>(key, c[j].getComb());

						// check if the end_v is valid globally
						if (end_v < g.size())
						{
							// Add an edge and set the the edge property to the size of the face (communication weight)
							gp.template addEdge(start_v, end_v, v_id, end_v);
						}
					}
				}
				++local_it;
			}
			++k_it;
		}
		return gp;
	}
};

/*! \brief This class construct a cartesian graph
 *
 * This class construct a cartesian graph
 *
 * \param dim dimensionality of the cartesian grid
 *
 */

template<unsigned int dim, typename Graph>
class DistGraphFactory
{

public:

	/*! \brief Construct a cartesian graph, with V and E edge properties
	 *
	 * Each vertex is a subspace (Hyper-cube) of dimension dim, each vertex is
	 * connected with an edge if two vertex (Hyper-cube) share a element of dimension grater than
	 * dim_c. One property can be used to store the contact size or the d-dimensional
	 * surface in common between two connected hyper-cube.
	 *
	 * \param sz Vector that store the size of the grid on each dimension
	 * \param dom Box enclosing the physical domain
	 *
	 * \tparam se Indicate which properties fill with the contact size. The
	 *           contact size is the point, line , surface, d-dimensional object size
	 *           in contact (in common) between two hyper-cube. NO_EDGE indicate
	 *           no property will store this information
	 * \tparam T type of the domain like (int real complex ... )
	 * \tparam dim_c Connectivity dimension
	 * \tparam pos... (optional)one or more integer indicating the spatial properties
	 *
	 */
	template<int se, typename T, unsigned int dim_c, int ... pos>
	static Graph construct(const size_t (&sz)[dim], Box<dim, T> dom)
	{
		return DistGraph_constr_impl<dim, Graph, se, T, dim_c, pos...>::construct(sz, dom);
	}
};

#endif /* DISTGRAPHFACTORY_HPP_ */
