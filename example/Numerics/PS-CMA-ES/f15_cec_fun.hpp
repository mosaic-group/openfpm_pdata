/*
 * f15_cec_fun.hpp
 *
 *  Created on: Jan 14, 2018
 *      Author: i-bird
 */

#ifndef EXAMPLE_NUMERICS_PS_CMA_ES_F15_CEC_FUN_HPP_
#define EXAMPLE_NUMERICS_PS_CMA_ES_F15_CEC_FUN_HPP_

#include "f15_cec_const.hpp"
#include <limits>
#include <math.h>

template<unsigned int dim>
void Job15(int funcnr,Eigen::VectorXd & vars,double & res)
{
      // local used vars
      double sum,sum1,sum2,prod,e1,e2;
      int i,j,k;
      // weierstrass vars
      int Kmax = 20;
      const double a_c = 0.5;
      const double b_c = 3.0;


      if (funcnr < 2)
      {
        // rastrigin
         sum = 10.0 * dim;
         for (size_t i = 0 ; i < dim ; i++)
         {
            sum += vars(i)*vars(i);
            sum -= 10.0*cos(2*M_PI*vars[i]);
         }

         res = sum;
      }
      else if (funcnr < 4)
      {
        // weierstrass
         sum1 = 0.0;
         sum2 = 0.0;
         double a_k = 1.0;
         double b_k = 1.0;
         for (size_t i = 0 ; i < dim ; i++)
         {
            a_k = 1.0;
            b_k = 1.0;
            for (size_t j = 0 ; j <= Kmax ; j++, a_k *= a_c,b_k *= b_c)
            {
                sum1 = sum1 + a_k * cos((M_PI)*2.0 * b_k * (vars(i)+0.5));
            }
         }
         a_k = 1.0;
         b_k = 1.0;
         for (size_t j = 0 ; j <= Kmax ; j++, a_k *= a_c, b_k *= b_c)
         {
            sum2 = sum2 + a_k * cos((M_PI)*2.0 * b_k * (0.5));
         }
         res = sum1 - sum2*dim;
      }
      else if (funcnr < 6)
      {
         // griewank
         prod = 1;
         sum = 0.0;
         for (size_t i = 1 ; i <= dim ; i++)
         {
        	 sum= sum + (vars(i-1)*vars(i-1))/4000.0;
        	 prod=prod * cos(vars(i-1)/(sqrt(double(i))));
         }
         res = sum-prod+1;
      }
      else if (funcnr < 8)
      {
        // ackley
        e1 = 0.0;
        e2 = 0.0;
        for (size_t i = 0 ; i < dim ; i++)
        {
            e1 = e1 + vars(i)*vars(i);
            e2 = e2 + cos(2.0*M_PI*vars(i));
        }
        res = exp(1.0) + 20.0 - 20*exp(-0.2*sqrt(e1/dim));
        res = res - exp(e2/dim);
      }
      else if (funcnr <= 10)
      {
        // sphere
        sum = vars.transpose() * vars;
        res = sum;
      }
}

template<unsigned int dim>
double hybrid_composition(Eigen::VectorXd & vars)
{
	double ZBQLNOR;

	//local used vars
	double wMax,sumSqr,wSum,w1mMaxPow;
	int i,j,k;
	double sumF,t_res;
	Eigen::VectorXd job_z[10];

	for (size_t i = 0 ; i < 10 ; i++)
	{job_z[i].resize(dim);}

	double job_w[10];
	double res = 0.0;

	for (size_t i = 0 ; i < dim ; i++)
	{
		if (vars[i] < -5.0 || vars[i] > 5.0)
    	{return std::numeric_limits<double>::infinity();}
	}

	// get the raw weights
    wMax = - std::numeric_limits<double>::max();
    for (size_t i = 0; i < 10 ; i++)
    {
    	sumSqr = 0.0;
        //Shift the Input
    	job_z[i] = vars - f15_o[i];
        sumSqr += (job_z[i].transpose() * job_z[i]);

        job_w[i] = exp(-1.0 * sumSqr / (2.0 * dim));

        if (wMax < job_w[i])
        {wMax = job_w[i];}
    }

    // Modify the weights
    wSum = 0.0;

    w1mMaxPow = 1.0 - wMax*wMax*wMax*wMax*wMax*wMax*wMax*wMax*wMax*wMax;
    for (size_t i = 0; i < 10 ; i++)
    {
            if (job_w[i] != wMax)
            {job_w[i] = job_w[i]* w1mMaxPow;};

            wSum = wSum + job_w[i];
    }

    // Normalize the weights
    for (size_t i = 0; i < 10 ; i++)
    {job_w[i] /= wSum;}

    sumF = 0.0;

    for (size_t i = 0; i < 10 ; i++)
    {
    	job_z[i] = job_z[i] / job_lambda[i];

        //calling the basic functions

        Job15<dim>(i,job_z[i],t_res);

    	sumF = sumF + job_w[i] * (2000.0*t_res/f15_max[i] + bias[i]);
    }

    res = sumF + 120;

    return res;
}

template<unsigned int dim>
void prepare_f15()
{
	// load f15_o
	for (size_t j = 0 ; j < 10 ; j++)
	{
		Eigen::VectorXd fmp(dim);
		f15_o[j].resize(dim);
		for (size_t i = 0 ; i < dim ; i++)
		{
			f15_o[j](i) = f15_const[j][i];
			fmp(i) = 5.0 / job_lambda[j];
		}

		double result;
		Job15<dim>(j,fmp,result);

		f15_max[j] = fabs(result);
	}
}

#endif /* EXAMPLE_NUMERICS_PS_CMA_ES_F15_CEC_FUN_HPP_ */
