/*
 * initialize_vcl.hpp
 *
 *  Created on: Aug 21, 2018
 *      Author: i-bird
 */

#ifndef INITIALIZE_VCL_HPP_
#define INITIALIZE_VCL_HPP_

/*! \brief If openfpm has to work on GPU we have to be sure openfpm_init is called on a file compiled with NVCC
 *
 * There are two implementation initialize.cpp and initialize.cu. In configuration stage the second implementation is chosen
 * if the test has to run on GPU
 *
 */
void openfpm_init_wrapper(int * argc, char *** argv);
void openfpm_finalize_wrapper();

#endif /* INITIALIZE_VCL_HPP_ */
