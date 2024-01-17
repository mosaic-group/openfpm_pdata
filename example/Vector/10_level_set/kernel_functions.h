#include <math.h>
// This is a collection of 1D, 2D and 3D SPH kernel functions.
//
// Created by lschulze

// Template typenames for kernel function types
struct WendlandC2 {};

template<size_t dim, typename kernel_function_type> class kernel_function
{
	public: 
		kernel_function(double smoothing_length) : H(smoothing_length)
		{
			// 1D Wendland C2
			if ((dim == 1) && (std::is_same<kernel_function_type, WendlandC2>::value))
					{
						kernel_type = 0;
						coefficient = 5.0/(8.0*H);
					}
			// 2D Wendland C2
			if ((dim == 2) && (std::is_same<kernel_function_type, WendlandC2>::value))
					{
						kernel_type = 1;
						coefficient = 7.0/(4.0*M_PI*H*H);
					}
			// 3D Wendland C2
			if ((dim == 3) && (std::is_same<kernel_function_type, WendlandC2>::value))
					{
						kernel_type = 2;
						coefficient = 21.0/(16.0*M_PI*H*H*H);
					}
		}

		double Wab(double r)
		{
    			const double q = r/H;
			// 1D Wendland C2
			if (kernel_type == 0)
			{
    				if (q <= 2.0)
				{
				        double factor = 1.0 - q/2.0;
        				factor = factor*factor*factor;
        				return(coefficient*factor*(1.5*q + 1));
    				}
    				else return(0.0);
			}
			// 2D Wendland C2
			if (kernel_type == 1)
			{
    				if (q <= 2.0)
    				{
        				double factor = 1.0 - q/2.0;
        				factor = factor*factor;
        				factor = factor*factor;
        				return(coefficient*factor*(1.0 + 2.0*q));
    				}
    				else return(0.0);
			}
			// 3D Wendland C2
			if (kernel_type == 2)
			{
    				if (q <= 2.0)
    				{
        				double factor = 1.0 - q/2.0;
        				factor = factor*factor;
        				factor = factor*factor;
        				return(coefficient*factor*(1.0 + 2.0*q));
    				}
    				else return(0.0);

			}
		}

		void DWab(Point<dim,double> & dx, Point<dim,double> & DW, double r, bool print, double & dwdrab)
		{
			const double q = r/H;
    			// 2D Wendland C2
			if(kernel_type == 1)
			{
    				if (q <= 2.0)
    				{
        				double factor = (-5.0*coefficient/(H))*q*(1.0 - q/2.0)*(1.0 - q/2.0)*(1.0 - q/2.0);

        				DW.get(0) = factor * dx.get(0)/r;
        				DW.get(1) = factor * dx.get(1)/r;

        				dwdrab = factor;
    				}
    				else
    				{
        				DW.get(0) = 0.0;
        				DW.get(1) = 0.0;

        				dwdrab = 0.0;
    				}
			}
			// 3D Wendland C2
			if (kernel_type == 2)
			{
    				if (q <= 2.0)
    				{
        				double factor = (-5.0*coefficient/(H))*q*(1.0 - q/2.0)*(1.0 - q/2.0)*(1.0 - q/2.0);

        				DW.get(0) = factor * dx.get(0)/r;
        				DW.get(1) = factor * dx.get(1)/r;
					DW.get(2) = factor * dx.get(2)/r;

        				dwdrab = factor;
    				}
    				else
    				{
        				DW.get(0) = 0.0;
        				DW.get(1) = 0.0;
					DW.get(2) = 0.0;

        				dwdrab = 0.0;
    				}
			}	
		}

	private:
		double H;
		int kernel_type;
		double coefficient;

};
