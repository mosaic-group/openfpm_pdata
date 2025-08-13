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

    typedef prop value_type;

    using vector_dist<dim,St,typename AggregateAppend<int,prop>::type,Decomposition,Memory,layout_base>::vector_dist;

    typedef boost::mpl::int_<AggregateAppend<int,prop>::type::max_prop-1> flag_prop;

    void setSubset(size_t key, int sub_id)
    {
        this->template getProp<flag_prop::value>(key) = sub_id;
    }

    int getSubset(size_t key)
    {
        return this->template getProp<flag_prop::value>(key);
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
        if(prop_names.size()==prop::max_prop){
            prop_names.add({"SubsetNumber"});
        }

        return vector_dist<dim,St,typename AggregateAppend<int,prop>::type,Decomposition,Memory,layout_base>::write_frame(out,iteration,opt);
    }
    inline bool write_frame(std::string out, size_t iteration,double time, int opt = VTK_WRITER)
    {
        auto &prop_names=this->getPropNames();
        if(prop_names.size()==prop::max_prop){
            prop_names.add({"SubsetNumber"});
        }

        return vector_dist<dim,St,typename AggregateAppend<int,prop>::type,Decomposition,Memory,layout_base>::write_frame(out,iteration,time,opt);
    }

    inline bool write(std::string out,int opt = VTK_WRITER)
    {
        auto &prop_names=this->getPropNames();
        if(prop_names.size()==prop::max_prop){
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
    typedef int yes_i_am_vector_subset;

    //! Subset detection
    typedef std::integral_constant<bool,true> is_it_a_subset;

    vector_dist_subset(vector_dist_ws<dim,St,prop,Decomposition,Memory,layout_base> & vd,
                       int sub_id)
                       :vd(vd),sub_id(sub_id)
    {
#ifdef SE_CLASS1
        subsetUpdate_ctr=vd.getMapCtr();
#endif
        // construct pid vector

        auto it = vd.getDomainIterator();

        while (it.isNext())
        {
            auto p = it.get();

            if (vd.template getProp<flag_prop::value>(p) == sub_id)
            {
                pid.add();
                pid.template get<0>(pid.size()-1) = p;
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
                pid.template get<0>(pid.size()-1) = p;
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

    /*! \brief return the local size of vector_dist_ws
    *
    * \return local size
    *
    */
    size_t size_local() const
    {
        return vd.size_local();
    }

    /*! \brief return the local size of the vector
    *
    * \return local size
    *
    */
    size_t size_local_subset() const
    {
        return pid.size();
    }

#ifndef ONLY_READWRITE_GETTER

    /*! \brief Get the position of an element
     *
     * see the vector_dist iterator usage to get an element key
     *
     * \param vector_dist_ws vec_key element
     *
     * \return the position of the element in space
     *
     */
    inline auto getPos(size_t vec_key) -> decltype(vd.getPos(vec_key))
    {
        return vd.getPos(vec_key);
    }

    /*! \brief Get the position of an element
     *
     * see the vector_dist iterator usage to get an element key
     *
     * \param vector_dist_ws vec_key element
     *
     * \return the position of the element in space
     *
     */
    inline auto getPos(size_t vec_key) const -> decltype(vd.getPos(vec_key))
    {
        return vd.getPos(vec_key);
    }

    /*! \brief Get the position of an element
     *
     * see the vector_dist iterator usage to get an element key
     *
     * \param vector_dist_subset vec_key element
     *
     * \return the position of the element in space
     *
     */
    inline auto getPosSubset(size_t vec_key) -> decltype(vd.getPos(vec_key))
    {
        return vd.getPos(pid.template get<0>(vec_key));
    }

    /*! \brief Get the position of an element
     *
     * see the vector_dist iterator usage to get an element key
     *
     * \param vector_dist_subset vec_key element
     *
     * \return the position of the element in space
     *
     */
    inline auto getPosSubset(size_t vec_key) const -> decltype(vd.getPos(vec_key))
    {
        return vd.getPos(pid.template get<0>(vec_key));
    }

    /*! \brief Move the memory from the device to host memory
     *
     * \tparam property to move use POS_PROP for position property
     *
     */
    template<unsigned int ... prp>
    inline void hostToDeviceProp()
    {
        vd.template hostToDeviceProp<prp ...>();
    }

    /*! \brief Move the memory from the device to host memory
     *
     * \tparam property to move use POS_PROP for position property
     *
     */
    inline void hostToDevicePos()
    {
        vd.template hostToDevicePos<0>();
    }

    /*! \brief Get the property of an element
    *
    * see the vector_dist iterator usage to get an element key
    *
    * \tparam id property id
    * \param vector_dist_ws vec_key vector element
    *
    * \return return the selected property of the vector element
    *
    */
    template<unsigned int id> inline auto getProp(size_t vec_key) -> decltype(vd.template getProp<id>(vec_key))
    {
        return vd.template getProp<id>(vec_key);
    }

    /*! \brief Get the property of an element
    *
    * see the vector_dist iterator usage to get an element key
    *
    * \tparam id property id
    * \param vector_dist_ws vec_key vector element
    *
    * \return return the selected property of the vector element
    *
    */
    template<unsigned int id> inline auto getProp(size_t vec_key) const -> decltype(vd.template getProp<id>(vec_key))
    {
        return vd.template getProp<id>(vec_key);
    }

    /*! \brief Get the property of an element
    *
    * see the vector_dist iterator usage to get an element key
    *
    * \tparam id property id
    * \param vector_dist_subset vec_key vector element
    *
    * \return return the selected property of the vector element
    *
    */
    template<unsigned int id> inline auto getPropSubset(size_t vec_key) -> decltype(vd.template getProp<id>(vec_key))
    {
        return vd.template getProp<id>(pid.template get<0>(vec_key));
    }

    /*! \brief Get the property of an element
    *
    * see the vector_dist iterator usage to get an element key
    *
    * \tparam id property id
    * \param vector_dist_subset vec_key vector element
    *
    * \return return the selected property of the vector element
    *
    */
    template<unsigned int id> inline auto getPropSubset(size_t vec_key) const -> decltype(vd.template getProp<id>(vec_key))
    {
        return vd.template getProp<id>(pid.template get<0>(vec_key));
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
    template<typename CellL = CellList<dim, St, Mem_fast<>, shift<dim, St>, typename std::remove_reference<decltype(vd.getPosVector())>::type > >
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

    /*! \brief Construct a Verlet List from the stored particles
     *
     * \tparam opt Verlet List option flag
     * \tparam VerletList_type Verlet List type to construct
     *
     * \param r_cut interation radius, or size of each cell
     * \param no_se3 avoid SE_CLASS3 checking
     *
     * \return the Verlet List
     *
     */
    template<
        unsigned int opt=VL_NON_SYMMETRIC,
        typename VerletList_type = VerletList<dim,St,opt,Mem_fast<>,shift<dim,St>,typename std::remove_reference<decltype(vd.getPosVector())>::type>
    >
    VerletList_type getVerlet(St r_cut, bool no_se3 = false)
    {
#ifdef SE_CLASS3
        if (no_se3 == false)
        {se3.getNN();}
#endif

        VerletList_type verletList;

        auto cellList = this->template getCellList<typename std::remove_reference<decltype(verletList.getInternalCellList())>::type>(r_cut, no_se3);

        if (verletList.getOpt() & VL_NON_SYMMETRIC)
        {
            auto it = getDomainIterator();
            verletList.Initialize(cellList, r_cut, it, vd.getPosVector(), vd.size_local());
        }

        else {
            std::cerr << __FILE__ << ":" << __LINE__ << " vector_dist_subset::getVerletList supports VL_NON_SYMMETRIC option only! " << std::endl;
        }

        return verletList;
    }


    /*! \brief Update an existing Verlet List
     *
     * \tparam opt Verlet List option flag
     * \tparam Mem_type Memory type of the Verlet List
     * \tparam vPose_type Position vector type of the Verlet List
     *
     * \param verletList Velet List to update
     * \param r_cut interation radius, or size of each cell
     *
     * \return the Verlet List
     *
     */
    template<
        unsigned int opt,
        typename Mem_type,
        typename vPos_type>
    void updateVerlet(VerletList<dim,St,opt,Mem_type,shift<dim,St>,vPos_type>& verletList, St r_cut, bool no_se3 = false)
    {
#ifdef SE_CLASS3
        if (no_se3 == false)
        {se3.getNN();}
#endif

        VerletList<dim,St,opt,Mem_type,shift<dim,St>,vPos_type> ver_tmp = getVerlet<opt, VerletList<dim,St,opt,Mem_type,shift<dim,St>,vPos_type>>(r_cut);
        verletList.swap(ver_tmp);
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
    template<typename CellL = CellList<dim, St, Mem_fast<>, shift<dim, St> > >
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
        cell_list.setGhostMarker(pid.size());
        cell_list.set_ndec(vd.getDecomposition().get_ndec());

        cell_list.clear();

        auto it = getDomainIterator();

        while (it.isNext())
        {
            size_t key = it.get();

            Point<dim,St> pos = getPos(key);

            cell_list.add(pos,key);

            ++it;
        }

        // Add also the ghost

        auto git = vd.getGhostIterator();

        while (git.isNext())
        {
            size_t key = git.get();

            Point<dim,St> pos = vd.getPos(key);

            if (vd.template getProp<flag_prop::value>(key) == sub_id)
            {
                cell_list.add(pos,key);
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


template<unsigned int dim, typename St, typename prop, typename Decomposition = CartDecomposition<dim,St,CudaMemory,memory_traits_lin>> using vector_dist_ws_gpu = vector_dist_ws<dim,St,prop,Decomposition,CudaMemory,memory_traits_lin>;
template<unsigned int dim, typename St, typename prop, typename Decomposition = CartDecomposition<dim,St,CudaMemory,memory_traits_lin>> using vector_dist_subset_gpu = vector_dist_subset<dim,St,prop,Decomposition,CudaMemory,memory_traits_lin>;

#endif //OPENFPM_PDATA_VECTOR_DIST_SUBSET_HPP
