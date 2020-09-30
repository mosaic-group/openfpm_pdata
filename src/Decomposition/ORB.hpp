/*
 * ORB.hpp
 *
 *  Created on: Mar 13, 2015
 *      Author: Pietro Incardona
 */

#ifndef ORB_HPP_
#define ORB_HPP_

#include "util/mathutil.hpp"

/*! \brief this class is a functor for "for_each" algorithm
 *
 * This class is a functor for "for_each" algorithm. For each
 * element of a range_c the operator() is called.
 * Is mainly used to unroll the bisect cycle inside one ORB cycle
 *
 */

template<unsigned int dim, typename ORB>
struct bisect_unroll
{
	ORB & orb;

	/*! \brief constructor
	 *
	 */
	bisect_unroll(ORB & orb)
	:orb(orb)
	{};

	//! It call the copy function for each property
    template<typename T>
    void operator()(T& t)
    {
    	orb.template bisect<T::value>();
    }
};

/*! \brief This structure use SFINAE to avoid instantiation of invalid code
 *
 * birect_unroll cannot be performed when i > dim (basically you will never cut
 * on a dimension that does not exist)
 *
 * What is happening here is that at compile-time enable_if< (i < dim) >::type
 * is evaluated, if the condition is true enable_if<true>::type is a valid expression
 * the the second specialization will be used, if is false enable_if<false>::type
 * is not a valid expression, and the first specialization is chosen
 *
 */

template<unsigned int dim, unsigned int i, typename ORB, class Enable=void>
struct do_when_dim_gr_i
{
	static void bisect_loop(bisect_unroll<dim,ORB> & bu)
	{
	}
};

template<unsigned int dim, unsigned int i, typename ORB>
struct do_when_dim_gr_i<dim,i,ORB,typename boost::enable_if< boost::mpl::bool_<(i < dim)> >::type>
{
	static void bisect_loop(bisect_unroll<dim,ORB> & bu, typename boost::enable_if< boost::mpl::bool_<(i < dim)> >::type* dummy = 0)
	{
		boost::mpl::for_each< boost::mpl::range_c<int,0,i> >(bu);
	}
};


/*! \brief ORB node
 * \tparam type of space float, double, complex ...
 *
 * it store the Center of mass of the local calculated particles
 * on one direction (only one direction, mean that one scalar is enough)
 *
 */

template<typename T> class ORB_node : public aggregate<T>
{
public:

	static const unsigned int CM = 0;
};

/*! \brief This class implement orthogonal recursive bisection
 *
 * \tparam dim Dimensionality of the ORB
 * \tparam T type of the space
 * \tparam loc_wg type of structure that store the local weight of particles
 * \tparam loc_pos type of structure that store the local position of particles
 * \tparam Box type of structure that contain the domain extension
 * \tparam Tree type of structure that store the tree structure
 *
 */

template<unsigned int dim, typename T, typename loc_wg=openfpm::vector<float>, typename loc_pos=openfpm::vector<Point<dim,T>> , typename Box=Box<dim,T>, template<typename,typename> class Tree=Graph_CSR_s>
class ORB
{
	// Virtual cluster
	Vcluster<> & v_cl;

	// particle coordinate accumulator
	openfpm::vector<T> cm;
	// number of elements summed into cm
	openfpm::vector<size_t> cm_cnt;

	// Local particle vector
	loc_pos & lp;

	// Label id of the particle, the label identify where the "local" particles
	// are in the local decomposition
	openfpm::vector<size_t> lp_lbl;

	size_t dim_div;

	// Structure that store the orthogonal bisection tree
	Tree<ORB_node<T>,no_edge> grp;

	/*! \brief Calculate the local center of mass on direction dir
	 *
	 * WARNING: with CM we mean the sum of the particle coordinate over one direction
	 *
	 * \tparam dir Direction witch to calculate the center of mass
	 *
	 * \param start from where the last leafs start
	 */
	template<unsigned int dir> void local_cm(size_t start)
	{
		typedef Point<dim,T> s;

		// reset the counters and accumulators
		cm.fill(0);
		cm_cnt.fill(0);

		// Calculate the local CM
		auto it = lp.getIterator();

		while (it.isNext())
		{
			// vector key
			auto key =  it.get();

			size_t lbl = lp_lbl.get(key) - start;

			// add the particle coordinate to the CM accumulator
			cm.get(lbl) += lp.template get<s::x>(key)[dir];

			// increase the counter
			cm_cnt.get(lbl)++;

			++it;
		}
	}

	/*! \brief It label the particles
	 *
	 * It identify where the particles are
	 *
	 * \param n level
	 *
	 */

	template<unsigned int dir> inline void local_label()
	{
		typedef Point<dim,T> s;

		// direction of the previous split
		const size_t dir_cm = (dir == 0)?(dim-1):(dir-1);

		// if it is the first iteration we just have to create the labels array with zero
		if (grp.getNVertex() == 1)
		{lp_lbl.resize(lp.size());lp_lbl.fill(0); return;}

		// we check where (the node) the particles live

		auto p_it = lp.getIterator();

		while(p_it.isNext())
		{
			auto key = p_it.get();

			// Get the label
			size_t lbl = lp_lbl.get(key);

			// each node has two child's (FOR SURE)
			// No leaf are processed here
			size_t n1 = grp.getChild(lbl,0);
			size_t n2 = grp.getChild(lbl,1);

			// get the splitting center of mass
			T cm = grp.template vertex_p<ORB_node<T>::CM>(lbl);

			// if the particle n-coordinate is smaller than the CM is inside the child n1
			// otherwise is inside the child n2

			if (lp.template get<s::x>(key)[dir_cm] < cm)
			{lp_lbl.get(key) = n1;}
			else
			{lp_lbl.get(key) = n2;}

			++p_it;
		}
	}

	/*! \brief Bisect the domains along one direction
	 *
	 * \tparam dir direction where to split
	 *
	 */

	template<unsigned int dir> size_t bisect()
	{
		//
		size_t start = grp.getNVertex();

		// first label the local particles
		local_label<dir>();

		// Index from where start the first leaf
		size_t n_node = (start + 1) / 2;

		// Calculate the local cm
		local_cm<dir>(start - n_node);

		// reduce the local cm and cm_cnt
		v_cl.sum(cm);
		v_cl.sum(cm_cnt);
		v_cl.execute();

		// set the CM for the previous leaf (they are not anymore leaf)

		for (size_t i = 0 ; i < n_node ; i++)
		{
			grp.template vertex_p<ORB_node<T>::CM>(start-n_node+i) = cm.get(i) / cm_cnt.get(i);
		}

		// append on the 2^(n-1) previous end node to the 2^n leaf

		for (size_t i = start - n_node, s = 0 ; i < start ; i++, s++)
		{
			// Add 2 leaf nodes and connect them with the node
			grp.addVertex();
			grp.addVertex();
			grp.template addEdge(i,start+2*s);
			grp.template addEdge(i,start+2*s+1);
		}

		return 2*n_node;
	}

	// define the friend class bisect_unroll

	friend class bisect_unroll<dim,ORB<dim,T,loc_wg,loc_pos,Box,Tree>>;

public:

	/*! \brief constructor
	 *
	 * \param dom Box domain
	 * \param n_sub number of sub-domain to create (it is approximated to the biggest 2^n number)
	 * \param lp Local position of the particles
	 */
	ORB(Box dom, size_t n_sub, loc_pos & lp) : v_cl(create_vcluster()), lp(lp)
	{
		typedef ORB<dim,T,loc_wg,loc_pos,Box,Tree> ORB_class;

		dim_div = 0;

		n_sub = openfpm::math::round_big_2(n_sub);
		size_t nsub = log2(n_sub);

		// number of center or mass needed
		cm.resize(2 << nsub);
		cm_cnt.resize(2 << nsub);

		// Every time we divide from 0 to dim-1, calculate how many cycle from 0 to dim-1 we have
		size_t dim_cycle = nsub / dim;

		// Create the root node in the graph
		grp.addVertex();

		// unroll bisection cycle
		bisect_unroll<dim,ORB_class> bu(*this);
		for (size_t i = 0 ; i < dim_cycle ; i++)
		{
			boost::mpl::for_each< boost::mpl::range_c<int,0,dim> >(bu);
			// bu is recreated several time internaly
		}

		// calculate and execute the remaining cycles
		switch(nsub - dim_cycle * dim)
		{
		case 1:
			do_when_dim_gr_i<dim,1,ORB_class>::bisect_loop(bu);
			break;
		case 2:
			do_when_dim_gr_i<dim,2,ORB_class>::bisect_loop(bu);
			break;
		case 3:
			do_when_dim_gr_i<dim,3,ORB_class>::bisect_loop(bu);
			break;
		case 4:
			do_when_dim_gr_i<dim,4,ORB_class>::bisect_loop(bu);
			break;
		case 5:
			do_when_dim_gr_i<dim,5,ORB_class>::bisect_loop(bu);
			break;
		case 6:
			do_when_dim_gr_i<dim,6,ORB_class>::bisect_loop(bu);
			break;
		case 7:
			do_when_dim_gr_i<dim,7,ORB_class>::bisect_loop(bu);
			break;
		default:
			// To extend on higher dimension create other cases or create runtime
			// version of local_cm and bisect
			std::cerr << "Error: " << __FILE__ << ":" << __LINE__ << " ORB is not working for dimension bigger than 8";

		}
	}
};


#endif /* ORB_HPP_ */
