/*
 * debug.hpp
 *
 *  Created on: Oct 17, 2018
 *      Author: i-bird
 */

#ifndef DEBUG_HPP_
#define DEBUG_HPP_

enum debug_run
{
	HOST,
	DEVICE
};

enum debug_iterator
{
	DOMAIN_IT,
	DOMAIN_GHOST_IT,
	GHOST_IT,
};



/*! \brief Find in a vector dist particles if "certain conditions"
 *
 * This function in in general used to debug. Many time in order to discover problem we want to control if particles/grid points has zeros or out of range values.
 * You can implement these check, by creating the loops manually or using the function debug_find with a lambda function.
 *
 *
 * @param vd distributed vector
 * @param fun_test function that define the test
 * @param fun_print function to print a message
 * @param it iterator type domain, domain + ghost, ghost
 * @param type_of_run optionally can do the scan on device
 *
 * /return This function return true if the debug function has been triggered
 *
 */
template<typename vector_dist_type, typename functor_test, typename functor_print>
bool debug_find(vector_dist_type & vd, functor_test fun_test, functor_print fun_print,
		        debug_iterator it = debug_iterator::DOMAIN_IT, debug_run type_of_run = debug_run::HOST,
		        bool print = true)
{
	vector_dist_iterator ite(0,0);

	if (it == debug_iterator::DOMAIN_IT)
	{ite = vd.getDomainIterator();}
	else if (it == debug_iterator::DOMAIN_GHOST_IT)
	{ite = vd.getDomainAndGhostIterator();}
	else
	{ite = vd.getGhostIterator();}

	bool test_tot = false;

	while (ite.isNext())
	{
		bool test = fun_test(ite.get().getKey());
		test_tot |= test;

		if (test == true && print == true)
		{std::cout << fun_print(ite.get().getKey()) << std::endl;}

		++ite;
	}

	return test_tot;
}

/*! \brief Find in a vector dist particles if "certain conditions"
 *
 * This function in in general used to debug. Many time in order to discover problem we want to control if particles/grid points has zeros or out of range values.
 * You can implement these check, by creating the loops manually or using the function debug_find with a lambda function.
 *
 *
 * @param vd distributed vector
 * @param fun_test function that define the test
 * @param fun_print function to print a message
 * @param it iterator type domain, domain + ghost, ghost
 * @param type_of_run optionally can do the scan on device
 *
 * /return This function return true if the debug function has been triggered
 *
 */
template<typename vector_type, typename functor_test, typename functor_print>
bool debug_find_single(vector_type vd, functor_test fun_test, functor_print fun_print,
		        size_t ghostMarker, debug_iterator it = debug_iterator::DOMAIN_IT, debug_run type_of_run = debug_run::HOST,
		        bool print = true)
{
	openfpm::vector_key_iterator ite(0,0);

	if (it == debug_iterator::DOMAIN_IT)
	{ite = vd.getIteratorTo(ghostMarker);}
	else if (it == debug_iterator::DOMAIN_GHOST_IT)
	{ite = vd.getIterator();}
	else
	{ite = vd.getIteratorFrom(ghostMarker);}

	bool test_tot = false;

	while (ite.isNext())
	{
		bool test = fun_test(ite.get());
		test_tot |= test;

		if (test == true && print == true)
		{std::cout << fun_print(ite.get()) << std::endl;}

		++ite;
	}

	return test_tot;
}

/*! \brief scan for a particle inside a box
 *
 *
 * \param vd vector_dist
 * \param box box
 * \param it iterator type
 * \param type_of_run check host or device memory
 * \return true if one or more is detected
 */
template<typename vector_dist_type>
bool debug_is_in_box(vector_dist_type & vd, Box<vector_dist_type::dims, typename vector_dist_type::stype> box,
		             debug_iterator it = debug_iterator::DOMAIN_IT, debug_run type_of_run = debug_run::HOST,
		             bool print = true)
{
	auto fun_test = [&vd,&box](unsigned int p){return box.isInside(vd.getPos(p));};
	auto fun_print = [&vd,&box](unsigned int p)
					 {
						std::stringstream message;
						message << "Debug info: detected particle p=" << p << " inside box: " << box.toString() << std::endl;
						return message.str();
					 };

	return debug_find(vd,fun_test,fun_print,it,type_of_run,print);
}


/*! \brief scan for a particle inside a box
 *
 *
 * \param vd vector_dist
 * \param box box
 * \param it iterator type
 * \param type_of_run check host or device memory
 * \return true if one or more is detected
 */
template<typename vector_type>
bool debug_is_in_box_single(vector_type & vd, Box<vector_type::value_type::dims, typename vector_type::value_type::coord_type> box,
		             size_t ghostMarker, std::string message = std::string() , debug_iterator it = debug_iterator::DOMAIN_IT, debug_run type_of_run = debug_run::HOST,
		             bool print = true)
{
	auto fun_test = [&vd,&box](unsigned int p){return box.isInside(vd.template get<0>(p));};
	auto fun_print = [&vd,&box,&message](unsigned int p)
					 {
						std::stringstream message_srm;
						message_srm << "Debug info: " << message << " detected particle p=" << p << " inside box: " << box.toString() << std::endl;
						return message_srm.str();
					 };

	return debug_find_single(vd,fun_test,fun_print,ghostMarker,it,type_of_run,print);
}

#endif /* DEBUG_HPP_ */
