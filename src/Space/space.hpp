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

  typedef boost::fusion::vector<> type;
  typedef typename memory_traits_inte<type>::type memory_int;
  typedef typename memory_traits_lin<type>::type memory_lin;

  type data;

  static const unsigned int max_prop = 0;
};


template<typename T> class space<1,T>
{
public:

  typedef boost::fusion::vector<T> type;
  typedef typename memory_traits_inte<type>::type memory_int;
  typedef typename memory_traits_lin<type>::type memory_lin;

  type data;

  static const unsigned int x = 0;

  static const unsigned int max_prop = 1;
};

template<typename T> class space<2,T>
{
public:

  typedef boost::fusion::vector<T,T> type;
  typedef typename memory_traits_inte<type>::type memory_int;
  typedef typename memory_traits_lin<type>::type memory_lin;

  type data;

  static const unsigned int x = 0;
  static const unsigned int y = 1;

  static const unsigned int max_prop = 2;
};

template<typename T> class space<3,T>
{
public:

  typedef boost::fusion::vector<T,T,T> type;
  typedef typename memory_traits_inte<type>::type memory_int;
  typedef typename memory_traits_lin<type>::type memory_lin;

  type data;

  static const unsigned int x = 0;
  static const unsigned int y = 1;
  static const unsigned int z = 2;

  static const unsigned int max_prop = 3;
};

template<typename T> class space<4,T>
{
public:

  typedef boost::fusion::vector<T,T,T,T> type;
  typedef typename memory_traits_inte<type>::type memory_int;
  typedef typename memory_traits_lin<type>::type memory_lin;

  type data;

  static const unsigned int x = 0;
  static const unsigned int y = 1;
  static const unsigned int z = 2;
  static const unsigned int t = 3;

  static const unsigned int max_prop = 4;
};

template<typename T> class space<5,T>
{
public:

  typedef boost::fusion::vector<T,T,T,T,T> type;
  typedef typename memory_traits_inte<type>::type memory_int;
  typedef typename memory_traits_lin<type>::type memory_lin;

  type data;

  static const unsigned int max_prop = 5;
};


template<typename T> class space<6,T>
{
public:

  typedef boost::fusion::vector<T,T,T,T,T,T> type;
  typedef typename memory_traits_inte<type>::type memory_int;
  typedef typename memory_traits_lin<type>::type memory_lin;

  type data;

  static const unsigned int max_prop = 6;
};


template<typename T> class space<7,T>
{
public:

  typedef boost::fusion::vector<T,T,T,T,T,T,T> type;
  typedef typename memory_traits_inte<type>::type memory_int;
  typedef typename memory_traits_lin<type>::type memory_lin;

  type data;

  static const unsigned int max_prop = 7;
};

template<typename T> class space<8,T>
{
public:

  typedef boost::fusion::vector<T,T,T,T,T,T,T,T> type;
  typedef typename memory_traits_inte<type>::type memory_int;
  typedef typename memory_traits_lin<type>::type memory_lin;

  type data;

  static const unsigned int max_prop = 8;
};

#endif /* SPACES_HPP_ */
