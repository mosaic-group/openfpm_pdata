/*
 * debug.hpp
 *
 *  Created on: Oct 7, 2016
 *      Author: i-bird
 */

#ifndef SRC_DEBUG_HPP_
#define SRC_DEBUG_HPP_

#include "Vector/map_vector.hpp"
#include "Space/Shape/Point.hpp"

Point<3,float> getPosPoint3f(openfpm::vector<Point<3,float>> & pos, size_t i)
{
	return Point<3,float>(pos.get(i));
}

Point<3,double> getPosPoint3d(openfpm::vector<Point<3,double>> & pos, size_t i)
{
	return Point<3,double>(pos.get(i));
}

Point<2,float> getPosPoint2f(openfpm::vector<Point<2,float>> & pos, size_t i)
{
	return Point<2,float>(pos.get(i));
}

Point<2,double> getPosPoint2d(openfpm::vector<Point<2,double>> & pos, size_t i)
{
	return Point<2,double>(pos.get(i));
}



#endif /* SRC_DEBUG_HPP_ */
