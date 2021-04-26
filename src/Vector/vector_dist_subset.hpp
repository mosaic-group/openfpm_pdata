//
// Created by Abhinav Singh on 24.02.20.
//

#ifndef OPENFPM_PDATA_VECTOR_DIST_SUBSET_HPP
#define OPENFPM_PDATA_VECTOR_DIST_SUBSET_HPP

#include "vector_dist.hpp"

template<unsigned int dim,
        typename St,
        typename prop,
        typename Decomposition = CartDecomposition<dim,St>,
        typename Memory = HeapMemory,
        template<typename> class layout_base = memory_traits_lin>
class vector_dist_ws : public vector_dist<dim,St,typename AggregateAppend<int,prop>::type,Decomposition,Memory,layout_base>
{
    public:

    using vector_dist<dim,St,typename AggregateAppend<int,prop>::type,Decomposition,Memory,layout_base>::vector_dist;

    typedef boost::mpl::int_<AggregateAppend<int,prop>::type::max_prop-1> flag_prop;

    void setSubset(vect_dist_key_dx key, int sub_id)
    {
        this->template getProp<flag_prop::value>(key) = sub_id;
    }

    void ghost_get_subset()
    {
        this->template ghost_get<flag_prop::value>(NO_POSITION | SKIP_LABELLING);
    }

    void getLastSubset(int sub_id)
    {
        this->template getProp<flag_prop::value>(this->size_local()-1) = sub_id;
    }

    inline bool write_frame(std::string out, size_t iteration, int opt = VTK_WRITER)
    {
        auto &prop_names=this->getPropNames();
        if(prop_names.size()<prop::max_prop+1){
            prop_names.add({"SubsetNumber"});
        }

        return vector_dist<dim,St,typename AggregateAppend<int,prop>::type,Decomposition,Memory,layout_base>::write_frame(out,iteration,opt);
    }

    inline bool write(std::string out,int opt = VTK_WRITER)
    {
        auto &prop_names=this->getPropNames();
        if(prop_names.size()<prop::max_prop+1){
            prop_names.add({"SubsetNumber"});
        }

        return vector_dist<dim,St,typename AggregateAppend<int,prop>::type,Decomposition,Memory,layout_base>::write(out,"",opt);
    }

};

template<unsigned int dim,
        typename St,
        typename prop,
        typename Decomposition = CartDecomposition<dim,St>,
        typename Memory = HeapMemory,
        template<typename> class layout_base = memory_traits_lin>
class vector_dist_subset
{
    typedef vector_dist_ws<dim,St,prop,Decomposition,Memory,layout_base> ivector_dist;

    typedef boost::mpl::int_<AggregateAppend<int,prop>::type::max_prop-1> flag_prop;

    ivector_dist & vd;

    openfpm::vector<aggregate<int>> pid;

    size_t sub_id;

#ifdef SE_CLASS1
    int subsetUpdate_ctr=0;
#endif

    void check_gm()
    {
        #ifdef SE_CLASS1
        for (size_t i = 0 ; i < pid.size() ; i++)
        {
            if (pid.template get<0>(i) >= vd.size_local())
            {
                std::cout << __FILE__ << ":" << __LINE__ << " Error you have ghost particles in your subset" << std::endl;
            }
        }
        #endif
    }

public:

    //! property object
    typedef typename ivector_dist::value_type value_type;

    typedef typename ivector_dist::Decomposition_type Decomposition_type;

    typedef typename ivector_dist::CellList_type CellList_type;

    //! space type
    typedef typename ivector_dist::stype stype;

    //! dimensions of space
    static const unsigned int dims = ivector_dist::dims;

    //!
    typedef int yes_i_am_vector_dist;

    //! Subset detection
    typedef std::integral_constant<bool,true> is_it_a_subset;

    vector_dist_subset(vector_dist_ws<dim,St,prop,Decomposition,Memory,layout_base> & vd,
                       int sub_id)
                       :vd(vd),pid(pid),sub_id(sub_id)
    {
        // construct pid vector

        auto it = vd.getDomainIterator();

        while (it.isNext())
        {
            auto p = it.get();

            if (vd.template getProp<flag_prop::value>(p) == sub_id)
            {
                pid.add();
                pid.template get<0>(pid.size()-1) = p.getKey();
            }

            ++it;
        }

    	check_gm();
    }

    void ghost_get_subset()
    {
        vd.template ghost_get<flag_prop::value>(NO_POSITION | SKIP_LABELLING);
    }

    /*! \brief Return the ids
     *
     * \return the ids of the subset
     * 
     */
    openfpm::vector<aggregate<int>> & getIds()
    {
        return pid;
    }

#ifdef SE_CLASS1
        int getUpdateCtr() const{
            return subsetUpdate_ctr;
        }
        int getMapCtr()
        {
            return vd.getMapCtr();
        }
#endif

    /*! \brief Update the subset indexes
     *
     */
    inline void update()
    {
#ifdef SE_CLASS1
        subsetUpdate_ctr=vd.getMapCtr();
#endif

        ghost_get_subset();

        pid.clear();

        auto it = vd.getDomainIterator();

        while (it.isNext())
        {
            auto p = it.get();

            if (vd.template getProp<flag_prop::value>(p) == sub_id)
            {
                pid.add();
                pid.template get<0>(pid.size()-1) = p.getKey();
            }

            ++it;
        }

        check_gm();
    }

    /*! \brief Get the decomposition
     *
     * \return
     *
     */
    inline Decomposition & getDecomposition()
    {
        return vd.getDecomposition();
    }

    /*! \brief Get the decomposition
     *
     * \return
     *
     */
    inline const Decomposition & getDecomposition() const
    {
        return vd.getDecomposition();
    }

    /*! \brief return the local size of the vector
 *
 * \return local size
 *
 */
    size_t size_local() const
    {
        return pid.size();
    }

    /*! \brief return the local size of the original vector
 *
 * \return local size
 *
 */
    size_t size_local_orig() const
    {
        return vd.size_local();
    }

#ifndef ONLY_READWRITE_GETTER

    /*! \brief Get the position of an element
     *
     * see the vector_dist iterator usage to get an element key
     *
     * \param vec_key element
     *
     * \return the position of the element in space
     *
     */
    inline auto getPos(vect_dist_key_dx vec_key) -> decltype(vd.getPos(vec_key))
    {
        return vd.getPos(vect_dist_key_dx(pid.template get<0>(vec_key.getKey())));
    }

    /*! \brief Get the position of an element
     *
     * see the vector_dist iterator usage to get an element key
     *
     * \param vec_key element
     *
     * \return the position of the element in space
     *
     */
    inline auto getPos(vect_dist_key_dx vec_key) const -> decltype(vd.getPos(vec_key))
    {
        return vd.getPos(vect_dist_key_dx(pid.template get<0>(vec_key.getKey())));
    }

    /*! \brief Get the position of an element
     *
     * see the vector_dist iterator usage to get an element key
     *
     * \param vec_key element
     *
     * \return the position of the element in space
     *
     */
    inline auto getPosOrig(vect_dist_key_dx vec_key) -> decltype(vd.getPos(vec_key))
    {
        return vd.getPos(vec_key.getKey());
    }

    /*! \brief Get the position of an element
     *
     * see the vector_dist iterator usage to get an element key
     *
     * \param vec_key element
     *
     * \return the position of the element in space
     *
     */
    inline auto getPosOrig(vect_dist_key_dx vec_key) const -> decltype(vd.getPos(vec_key))
    {
        return vd.getPos(vec_key.getKey());
    }

    /*! \brief Get the property of an element
 *
 * see the vector_dist iterator usage to get an element key
 *
 * \tparam id property id
 * \param vec_key vector element
 *
 * \return return the selected property of the vector element
 *
 */
    template<unsigned int id> inline auto getProp(vect_dist_key_dx vec_key) -> decltype(vd.template getProp<id>(vec_key))
    {
        return vd.template getProp<id>(vect_dist_key_dx(pid.template get<0>(vec_key.getKey())));
    }

    /*! \brief Get the property of an element
*
* see the vector_dist iterator usage to get an element key
*
* \tparam id property id
* \param vec_key vector element
*
* \return return the selected property of the vector element
*
*/
    template<unsigned int id> inline auto getProp(vect_dist_key_dx vec_key) const -> decltype(vd.template getProp<id>(vec_key))
    {
        return vd.template getProp<id>(vect_dist_key_dx(pid.template get<0>(vec_key.getKey())));
    }


	vect_dist_key_dx getOriginKey(vect_dist_key_dx vec_key)
	{
		return vect_dist_key_dx(pid.template get<0>(vec_key.getKey()));
	}

#endif

    /*! \brief Get an iterator that traverse the particles in the domain
     *
     * \return an iterator
     *
     */
    vector_dist_iterator_subset getDomainIterator() const
    {
#ifdef SE_CLASS3
        se3.getIterator();
#endif

        return vector_dist_iterator_subset(0,pid.size(),pid);
    }

	/*! \brief Construct a cell list starting from the stored particles
	 *
	 * \tparam CellL CellList type to construct
	 *
	 * \param r_cut interation radius, or size of each cell
	 * \param no_se3 avoid SE_CLASS3 checking
	 *
	 * \return the Cell list
	 *
	 */
	template<typename CellL = CellList_gen<dim, St, Process_keys_lin, Mem_fast<>, shift<dim, St>, typename std::remove_reference<decltype(vd.getPosVector())>::type > >
	CellL getCellList(St r_cut, bool no_se3 = false)
	{
#ifdef SE_CLASS3
		if (no_se3 == false)
		{se3.getNN();}
#endif

		// Get ghost and anlarge by 1%
		Ghost<dim,St> g = vd.getDecomposition().getGhost();
		g.magnify(1.013);

		return getCellList<CellL>(r_cut, g,no_se3);
	}

	/*! \brief Indicate that this class is not a subset
	 *
	 * \return false
	 *
	 */
	bool isSubset() const
	{
		return true;
	}

	/*! \brief Construct a cell list starting from the stored particles
	 *
	 * It differ from the get getCellList for an additional parameter, in case the
	 * domain + ghost is not big enough to contain additional padding particles, a Cell list
	 * with bigger space can be created
	 * (padding particles in general are particles added by the user out of the domains)
	 *
	 * \tparam CellL CellList type to construct
	 *
	 * \param r_cut interation radius, or size of each cell
	 * \param enlarge In case of padding particles the cell list must be enlarged, like a ghost this parameter say how much must be enlarged
	 * \param no_se3 avoid se_class3 cheking default false
	 *
	 * \return the CellList
	 *
	 */
	template<typename CellL = CellList_gen<dim, St, Process_keys_lin, Mem_fast<>, shift<dim, St> > >
	CellL getCellList(St r_cut, const Ghost<dim, St> & enlarge, bool no_se3 = false)
	{
#ifdef SE_CLASS3
		if (no_se3 == false)
		{se3.getNN();}
#endif

		CellL cell_list;

		// Division array
		size_t div[dim];

		// get the processor bounding box
		Box<dim, St> pbox = vd.getDecomposition().getProcessorBounds();

		// Processor bounding box
		cl_param_calculate(pbox, div, r_cut, enlarge);

		cell_list.Initialize(pbox, div);
		cell_list.set_gm(pid.size());
		cell_list.set_ndec(vd.getDecomposition().get_ndec());

		cell_list.clear();

		auto it = getDomainIterator();

		while (it.isNext())
		{
			auto key = it.get();

			Point<dim,St> pos = getPos(key);

			cell_list.add(pos,pid.template get<0>(key.getKey()));

			++it;
		}

        // Add also the ghost

        auto git = vd.getGhostIterator();

		while (git.isNext())
		{
			auto key = git.get();

			Point<dim,St> pos = vd.getPos(key);

            if (vd.template getProp<flag_prop::value>(key) == sub_id)
            {
			    cell_list.add(pos,key.getKey());
            }

			++git;
		}

		return cell_list;
	}

	/*! \brief Operator= for distributed vector
	 *
	 * \param v vector to copy
	 *
	 * \return itself
	 *
	 */
	vector_dist_subset<dim,St,prop,Decomposition,Memory,layout_base> &
	operator=(const vector_dist_subset<dim,St,prop,Decomposition,Memory,layout_base> & v)
	{
	    vd = v.vd;
	    pid = v.pid;

		return *this;
	}
};



#endif //OPENFPM_PDATA_VECTOR_DIST_SUBSET_HPP
