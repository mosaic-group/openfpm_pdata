/*
 * util.hpp
 *
 *  Created on: Nov 20, 2014
 *      Author: Pietro Incardona
 */

/*! \brief check that two array match
 *
 * check that two array match
 *
 * \param ptr1 Point to array 1
 * \param ptr2 Pointer to array 2
 * \param sz Size of the array
 *
 */

template<typename T> void boost_check_array(const T * ptr1, const T * ptr2, size_t sz)
{
	// Check the array
	for (int i = 0 ; i < sz ; i++)
	{
		BOOST_REQUIRE_EQUAL(ptr1[i],ptr2[i]);
	}
}

