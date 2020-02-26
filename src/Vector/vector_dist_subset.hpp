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
class vector_dist_subset
{
    typedef vector_dist<dim,St,prop,Decomposition,Memory,layout_base> ivector_dist;

    ivector_dist & vd;

    openfpm::vector<aggregate<int>> & pid;

    size_t g_m = 0;

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

    vector_dist_subset(vector_dist<dim,St,prop,Decomposition,Memory,layout_base> & vd,
                       openfpm::vector<aggregate<int>> & pid)
                       :vd(vd),pid(pid)
    {
        for (size_t i = 0 ; i < pid.size() ; i++)
        {
            g_m += (pid.template get<0>(i) < vd.size_local())?1:0;
        }
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
        return g_m;
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


#endif

    /*! \brief Get an iterator that traverse the particles in the domain
     *
     * \return an iterator
     *
     */
    vector_dist_iterator getDomainIterator() const
    {
#ifdef SE_CLASS3
        se3.getIterator();
#endif

        return vector_dist_iterator(0, g_m);
    }



};

#endif //OPENFPM_PDATA_VECTOR_DIST_SUBSET_HPP
