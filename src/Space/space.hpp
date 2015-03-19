/*
 * spaces.hpp
 *
 *  Created on: Mar 10, 2015
 *      Author: Pietro Incardona
 */

#ifndef SPACES_HPP_
#define SPACES_HPP_

template<unsigned int dim, typename T> class space
{
public:

  typedef boost::fusion::vector<T[dim]> type;
  typedef typename memory_traits_inte<type>::type memory_int;
  typedef typename memory_traits_lin<type>::type memory_lin;

  type data;

  static const unsigned int x = 0;

  static const unsigned int max_prop = 1;
};


#endif /* SPACES_HPP_ */
