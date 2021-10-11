#include <math.h>
#include <sys/_types/_size_t.h>
#include "Vector/vector_dist.hpp"
#include "DCPSE/Dcpse.hpp"
#include "DCPSE/MonomialBasis.hpp"
#include "Draw/DrawParticles.hpp"

const double dp = 1/512.0;
const double A = 0.75;
const double B = 0.5;
const double band_width = 12.0*dp;

const int sdf = 0;
const int sdfgrad = 1;
const int curvature = 2;
const int surf_flag = 3;
const int num_neibs = 4;
const int min_sdf = 5;
const int min_sdf_x = 6;
const int interpol_coeff = 7;
const int sdf_analytical = 8;

typedef vector_dist<2, double, aggregate<double, double[2], double, int, int, double, double[2], double[16], double>> particles;
//										|		|	 	 |         |    |		 |		|	|						|
//									     sdf  sdfgrad curvature surf_flag num_neibs min sdf in support	min sdf location poly coefficients analytical sdf

typedef vector_dist<2, double, aggregate<vect_dist_key_dx>> particles_surface;

//
double get_distance(double x, double y)
{
	double a = A;
	double b = B;
	double xcp = 0.0;
	double ycp = 0.0;
	double dist = 15.0;
	double t = -b*b + b*y;
	double dt = 0.0;
	double eps = 1.0;


	int k = 0;
	std::cout<<"..................................................."<<std::endl;
	std::cout<<x<<"\t"<<y<<std::endl;

	if ((x > 0) && (y > 0))
	{
		while (eps>0.00000001)
		{
			dt = (((a*x)/(t + a*a))*((a*x)/(t + a*a)) + ((b*y)/(t + b*b))*((b*y)/(t + b*b)) - 1)/(-2*a*a*x*x/((t + a*a)*(t + a*a)*(t + a*a)) - 2*b*b*y*y/((t + b*b)*(t + b*b)*(t +b*b)));
			t = t + dt;
			eps = abs(dt);
			//std::cout<<eps<<std::endl;
			k++;
		}
	}

	if (y>0)
	{
		if (x>0)
		{
			xcp = a*a*x/(t + a*a);
			ycp = b*b*y/(t + b*b);
			dist = sqrt((xcp - x)*(xcp - x) + (ycp - y)*(ycp - y));
		}
		else
		{
			xcp = 0.0;
			ycp = b;
			dist = abs(y - b);
		}
	}
	else
	{
		if (x < (a*a - b*b)/a)
		{
			xcp = a*a*x/(a*a - b*b);
			ycp = a*sqrt(1 - (xcp/a)*(xcp/a));
			dist = a*sqrt((xcp - x)*(xcp - x) + ycp*ycp);
		}
		else
		{
			xcp = a;
			ycp = 0.0;
			dist = abs(x - a);
		}
	}
	std::cout<<dist<<std::endl;
	return dist;
}


double GetRoot ( double r0 , double z0 , double z1 , double g )
    {
		const int maxIter = 100;
        double n0 = r0*z0;
        double s0 = z1 - 1;
        double s1 = ( g < 0 ? 0 : sqrt(n0*n0+z1*z1) - 1 ) ;
        double s = 0;
        for ( int i = 0; i < maxIter; ++i ){
            s = ( s0 + s1 ) / 2 ;
            if ( s == s0 || s == s1 ) {break; }
            double ratio0 = n0 /( s + r0 );
            double ratio1 = z1 /( s + 1 );
            g = ratio0*ratio0 + ratio1*ratio1 - 1 ;
            if (g > 0) {s0 = s;} else if (g < 0) {s1 = s ;} else {break ;}
        }
        return s;
    }
double DistancePointEllipse( double e0 , double e1 , double y0 , double y1)
    {
        double distance;
        double x0 = 0.0;
        double x1 = 0.0;
        if ( y1 > 0){
            if ( y0 > 0){
                double z0 = y0 / e0;
                double z1 = y1 / e1;
                double g = z0*z0+z1*z1 - 1;
                if ( g != 0){
                    double r0 = (e0/e1)*(e0/e1);
                    double sbar = GetRoot(r0 , z0 , z1 , g);
                    x0 = r0 * y0 /( sbar + r0 );
                    x1 = y1 /( sbar + 1 );
                    distance = sqrt( (x0-y0)*(x0-y0) + (x1-y1)*(x1-y1) );
                    }else{
                        x0 = y0;
                        x1 = y1;
                        distance = 0;
                    }
                }
                else // y0 == 0
                    {x0 = 0 ; x1 = e1 ; distance = abs( y1 - e1 );}
        }else{ // y1 == 0
            double numer0 = e0*y0 , denom0 = e0*e0 - e1*e1;
            if ( numer0 < denom0 ){
                    double xde0 = numer0/denom0;
                    x0 = e0*xde0 ; x1 = e1*sqrt(1 - xde0*xde0 );
                    distance = sqrt( (x0-y0)*(x0-y0) + x1*x1 );
                }else{
                    x0 = e0;
                    x1 = 0;
                    distance = abs( y0 - e0 );
            }
        }
        return distance;
    }


int return_sign(double phi)
{
	if (phi > 0) return 1;
	if (phi < 0) return -1;
	return 0;
}

const double H = 1.3*dp;
const double a2 = 7.0/4.0/M_PI/H/H;

inline double Wab(double r)
{
	r /= H;

	if (r < 2.0) {
        double temp = 1.0 - r/2.0;
        temp = temp*temp;
        temp = temp*temp;
	    return ((temp*(1 + 2.0*r)) * a2);
    }
		else
		return 0.0;
}

const double c1 = -5.0*a2/H/H;

// Filled later
double W_dap = 0.0;

inline void DWab(Point<2,double> & dx, Point<2,double> & DW, double r, bool print)
{
	const double qq=r/H;
	double factor = 0.0;
	if (qq < 2.0) {
        double temp = 1.0 - (qq / 2.0);
        factor = c1*temp*temp*temp;
    }

    DW.get(0) = factor * dx.get(0);
    DW.get(1) = factor * dx.get(1);
}

inline double DWdr(double r)
{
	r /= H;
	if (r < 2.0) {
		double temp = 1.0 - (r / 2.0);
		return (r*H*c1*temp*temp*temp);
	}
	else return 0.0;
}


// Calculate surface normals as an SPH difference gradient (A_i - A_j). This function also tracks the number of neighbors per particle
// for the following redistancing.
template<typename CellList> inline void calc_surface_normals_sph(particles & vd, CellList & NN)
{	
	double sdf_gradient_norm = 0.0;
	int sdf_count = 0;
    auto part = vd.getDomainIterator();

    // Update the cell-list
    vd.updateCellList(NN);

    // For each particle ...
    while (part.isNext())
    {
        // ... a
        auto a = part.get();
        //reset the values
        vd.getProp<sdfgrad>(a)[0] = 0.0;
        vd.getProp<sdfgrad>(a)[1] = 0.0;
        // vd.getProp<num_neibs>(a) = 1;

        // Get the position xp of the particle
        Point<2,double> xa = vd.getPos(a);

        // For the purpose of only computing the gradient, say V_a = V_b = H^3

        // Get an iterator over the neighborhood particles of p
        auto Np = NN.template getNNIterator<NO_CHECK>(NN.getCell(vd.getPos(a)));

        // For each neighborhood particle
        while (Np.isNext() == true)
        {
            // ... q
            auto b = Np.get();

            // Get the position xp of the particle
            Point<2,double> xb = vd.getPos(b);

            // if (p == q) skip this particle
            if (a.getKey() == b)	{++Np; continue;};

            // Get the distance between p and q
            Point<2,double> dr = xa - xb;
            // take the norm of this vector
            double r2 = norm2(dr);

            // if they interact
            if (r2 < 4.0*H*H)
            {
                double r = sqrt(r2); // norm2 is norm^2

                Point<2,double> DW;
                DWab(dr,DW,r,false);

                //double cfactor = rhoa/massa*((massa/rhoa)*(massa/rhoa)+(massb/rhob)*(massb/rhob))*coloraverage(rhoa,rhob,vd.getProp<type>(a),vd.getProp<type>(b));
                double cfactor = dp*dp*(vd.getProp<sdf>(b) - vd.getProp<sdf>(a));
                vd.getProp<sdfgrad>(a)[0] += cfactor * DW.get(0);
                vd.getProp<sdfgrad>(a)[1] += cfactor * DW.get(1);

            }

            ++Np;
        }

		if (vd.template getProp<surf_flag>(a) == 1)
		{	
			Point<2, double> sdf_gradient = vd.getProp<sdfgrad>(a);
			sdf_gradient_norm += norm(sdf_gradient);
			++sdf_count;
		}
      
        ++part;
    }
    sdf_gradient_norm = sdf_gradient_norm/sdf_count;
    std::cout<<sdf_gradient_norm<<std::endl;
}

double randZeroToOne()
{
    return (rand() / (RAND_MAX + 1.));//
}

inline void perturb_pos(particles & vd)
{
	auto part = vd.getDomainIterator();

	while(part.isNext())
	{
		auto a = part.get();

		const double x = vd.getPos(a)[0];
		const double y = vd.getPos(a)[1];
		const double dist = sqrt(x*x + y*y);

		vd.getPos(a)[0] += x*randZeroToOne()*dp/3.0;
		vd.getPos(a)[1] += y*randZeroToOne()*dp/3.0;

		//vd.template getProp<d>(a) = std::sqrt(vd.getPos(a)[0]*vd.getPos(a)[0] + vd.getPos(a)[1]*vd.getPos(a)[1]);

		++part;
	}
}

template <typename CellList> inline void detect_surface_particles(particles & vd, CellList & NN, particles_surface & vd_s)
{
	auto part = vd.getDomainIterator();
	vd.updateCellList(NN);
	while (part.isNext())
	{
		auto a = part.get();
		int sgn_a = return_sign(vd. template getProp<sdf>(a));
		Point<2,double> xa = vd.getPos(a);
		int num_neibs_a = 0;
		double min_sdf = vd.getProp<sdf>(a);
		decltype(a) min_sdf_key = a;

		auto Np = NN.template getNNIterator<NO_CHECK>(NN.getCell(vd.getPos(a)));
		while (Np.isNext())
		{
			auto b = Np.get();
			int sgn_b = return_sign(vd.getProp<sdf>(b));
			Point<2, double> xb = vd.getPos(b);

            Point<2,double> dr = xa - xb;
            double r2 = norm2(dr);

            if (r2 < 4.0*H*H)
            {
            	++num_neibs_a;
            	if (sgn_a != sgn_b)
            	{
            		vd.template getProp<surf_flag>(a) = 1;
            		//break;
            	}

            	if (vd.getProp<sdf>(b) < min_sdf)
            	{
            		min_sdf = vd.getProp<sdf>(b);
            		min_sdf_key = b;
            	}
            }
            ++Np;
            
		}

		if (vd.getProp<surf_flag>(a))
		{
			vd_s.add();
			vd_s.getLastProp<0>() = min_sdf_key;
			vd_s.getLastPos()[0] = vd.getPos(min_sdf_key)[0];
			vd_s.getLastPos()[1] = vd.getPos(min_sdf_key)[1];
			vd.getProp<num_neibs>(a) = num_neibs_a;
		}
		++part;
	}

}

template <typename CellList> inline void update_sdfs(particles & vd, CellList & NN)
{
	auto part = vd.getDomainIterator();
	vd.updateCellList(NN);
	while (part.isNext())
	{
		auto a = part.get();
		Point<2,double> xa = vd.getPos(a);
		double x = xa[0];
		double y = xa[1];
		vd.getProp<sdf>(a) = 1 - sqrt((x/A)*(x/A) + (y/B)*(y/B));
		vd.getProp<sdf_analytical>(a) = return_sign(vd.getProp<sdf>(a))*DistancePointEllipse(A, B, abs(x), abs(y));

		++part;
	}
}


inline EMatrix<double, Eigen::Dynamic, Eigen::Dynamic> compute_inverse(EMatrix<double, Eigen::Dynamic, Eigen::Dynamic> H)
{
	return H;
}

inline double get_p(EMatrix<double, Eigen::Dynamic, 1> xvector, EMatrix<double, Eigen::Dynamic, 1> c)
{
	const double x = xvector[0];
	const double y = xvector[1];
	return(c[0] + c[1]*x + c[2]*x*x + c[3]*x*x*x + c[4]*y + c[5]*x*y + c[6]*x*x*y + c[7]*x*x*x*y + c[8]*y*y + c[9]*x*y*y + c[10]*x*x*y*y + c[11]*x*x*x*y*y +c[12]*y*y*y + c[13]*x*y*y*y + c[14]*x*x*y*y*y + c[15]*x*x*x*y*y*y);
}

inline double get_dpdx(EMatrix<double, Eigen::Dynamic, 1> xvector, EMatrix<double, Eigen::Dynamic, 1> c)
{	
	const double x = xvector[0];
	const double y = xvector[1];
	return(c[1] + 2*c[2]*x + 3*c[3]*x*x + c[5]*y + 2*c[6]*x*y + 3*c[7]*x*x*y + c[9]*y*y + 2*c[10]*x*y*y + 3*c[11]*x*x*y*y + c[13]*y*y*y + 2*c[14]*x*y*y*y + 3*c[15]*x*x*y*y*y);
}

inline double get_dpdy(EMatrix<double, Eigen::Dynamic, 1> xvector, EMatrix<double, Eigen::Dynamic, 1> c)
{	
	const double x = xvector[0];
	const double y = xvector[1];
	return(c[4] + c[5]*x + c[6]*x*x + c[7]*x*x*x + 2*c[8]*y + 2*c[9]*x*y + 2*c[10]*x*x*y + 2*c[11]*x*x*x*y + 3*c[12]*y*y + 3*c[13]*x*y*y + 3*c[14]*x*x*y*y + 3*c[15]*x*x*x*y*y);
}

inline double get_dpdxdx(EMatrix<double, Eigen::Dynamic, 1> xvector, EMatrix<double, Eigen::Dynamic, 1> c)
{
	const double x = xvector[0];
	const double y = xvector[1];
	return(2*c[2] + 6*c[3]*x + 2*c[6]*y + 6*c[7]*x*y + 2*c[10]*y*y + 6*c[11]*y*y*x + 2*c[14]*y*y*y + 6*c[15]*y*y*y*x);
}

inline double get_dpdydy(EMatrix<double, Eigen::Dynamic, 1> xvector, EMatrix<double, Eigen::Dynamic, 1> c)
{
	const double x = xvector[0];
	const double y = xvector[1];
	return(2*c[8] + 2*c[9]*x + 2*c[10]*x*x + 2*c[11]*x*x*x + 6*c[12]*y + 6*c[13]*x*y + 6*c[14]*x*x*y + 6*c[15]*x*x*x*y);
}

inline double get_dpdxdy(EMatrix<double, Eigen::Dynamic, 1> xvector, EMatrix<double, Eigen::Dynamic, 1> c)
{
	const double x = xvector[0];
	const double y = xvector[1];
	return(c[5] + 2*c[6]*x + 3*c[7]*x*x + 2*c[9]*y + 4*c[10]*x*y + 6*c[11]*x*x*y + 3*c[13]*y*y + 6*c[14]*x*y*y + 9*c[15]*x*x*y*y);
}

template<typename CellList> inline void calc_derivatives(particles & vd, CellList & NN)
{
	auto part = vd.getDomainIterator();
	vd.updateCellList(NN);

	while(part.isNext())
	{
		auto a = part.get();
		if (vd.getProp<surf_flag>(a) == 1)
		{
			EMatrix<double, Eigen::Dynamic, 1> c(16, 1);
			for (int k = 0; k<16; k++)
			{
				c[k] = vd.getProp<interpol_coeff>(a)[k];
			}

			EMatrix<double, Eigen::Dynamic, 1> x(2, 1);
			x[0] = vd.getPos(a)[0];
			x[1] = vd.getPos(a)[1];
			//std::cout<<"these is the position:\n"<<x<<std::endl;
			//std::cout<<"these are the coefficients:\n"<<c<<std::endl;

			vd.getProp<sdfgrad>(a)[0] = get_dpdx(x, c);
			vd.getProp<sdfgrad>(a)[1] = get_dpdy(x, c);
			const double sdfgradmag = sqrt(vd.getProp<sdfgrad>(a)[0]*vd.getProp<sdfgrad>(a)[0] + vd.getProp<sdfgrad>(a)[1]*vd.getProp<sdfgrad>(a)[1]);

			// vd.getProp<curvature>(a) = get_dpdxdx(x, c) + get_dpdydy(x, c); //could be changed.
			// vd.getProp<curvature>(a) = (get_dpdxdx(x,c) + get_dpdydy(x, c))/sdfgradmag + (get_dpdx(x,c)*(get_dpdxdx(x,c) + get_dpdxdy(x,c)) + get_dpdy(x,c)*(get_dpdxdy(x,c) + get_dpdydy(x,c)))/std::pow(sdfgradmag, 3.0);
			vd.getProp<curvature>(a) =(get_dpdxdx(x, c)*get_dpdy(x, c)*get_dpdy(x, c) - 2.0*get_dpdy(x, c)*get_dpdx(x, c)*get_dpdxdy(x, c) + get_dpdydy(x, c)*get_dpdx(x, c)*get_dpdx(x, c))/std::pow(std::sqrt(get_dpdx(x, c)*get_dpdx(x, c) + get_dpdy(x, c)*get_dpdy(x, c)), 3.0);
			//std::cout<<sdfgradmag<<std::endl;

			//std::cout<<"sdfgrad magnitude: "<<sdfgradmag<<"\ncurvature: "<<vd.getProp<curvature>(a)<<std::endl;
		}
		else
		{
			vd.getProp<sdfgrad>(a)[0] = 0.0;
			vd.getProp<sdfgrad>(a)[1] = 0.0;
			vd.getProp<curvature>(a) = 0.0;
		}

		++part;

	}
}

template <typename CellList> inline void cp_interpol(particles & vd, CellList & NN)
{
	auto part = vd.getDomainIterator();
	vd.updateCellList(NN);
	while (part.isNext())
		{
			auto a = part.get();
			if (vd.template getProp<surf_flag>(a) != 1)
			{
				++part;
				continue;
			}

			const int num_neibs_a = vd.getProp<num_neibs>(a);

			Point<2, double> xa = vd.getPos(a);

			// debug query point a
			// std::cout<<"a: "<<xa[0]<<", "<<xa[1]<<std::endl;

			double min_sdf_val = std::abs(vd.getProp<sdf>(a));
			Point<2, double> min_sdf_val_x = xa;
			int neib = 0;
			//std::vector<unsigned int> mo = {1, 0};

			// order limit 4 corresponds to bicubic basis functions
			MonomialBasis<2> m(4);
			// std::cout<<"m size"<<m.size()<<std::endl;
			VandermondeRowBuilder<2, double> vrb(m);

			EMatrix<double, Eigen::Dynamic,Eigen::Dynamic> V(num_neibs_a, m.size());
			EMatrix<double, Eigen::Dynamic, 1> phi(num_neibs_a, 1);

			auto Np = NN.template getNNIterator<NO_CHECK>(NN.getCell(vd.getPos(a)));
			while(Np.isNext())
			{
				auto b = Np.get();
				Point<2, double> xb = vd.getPos(b);
				Point<2, double> dr = xa - xb;
				double r2 = norm2(dr);

				if (r2 > 4.0*H*H)
				{
					++Np;
					continue;
				}

				// debug given data
				// std::cout<<std::setprecision(15)<<xb[0]<<"\t"<<xb[1]<<"\t"<<vd.getProp<sdf>(b)<<std::endl;


				if (std::abs(vd.getProp<sdf>(b)) < min_sdf_val)
				{
					min_sdf_val = std::abs(vd.getProp<sdf>(b));
					min_sdf_val_x = xb;
				}

				// Fill phi-vector from the right hand side
				phi[neib] = vd.getProp<sdf>(b);
				//xb[0] = 0.1;
				//xb[1] = 0.5;

				vrb.buildRow(V, neib, xb, 1.0);


				++neib;
				++Np;
			}

			// debug matrix A
			// std::cout<<"A:\n"<<V<<std::endl;

			EMatrix<double, Eigen::Dynamic, 1> c(m.size(), 1);
			c = V.completeOrthogonalDecomposition().solve(phi);

			EMatrix<double, Eigen::Dynamic, 1> xaa(2,1);
			xaa[0] = xa[0];
			xaa[1] = xa[1];
			double curvature = get_dpdxdx(xaa, c) + get_dpdydy(xaa, c);

			int print = 0;
			for (int k = 0; k<16; k++)
			{
				vd.getProp<interpol_coeff>(a)[k] = c[k];

			}

			//if ((print || (curvature > 100000.0)) && (vd.getProp<surf_flag>(a) == 1))
			if (vd.getProp<surf_flag>(a) == 10)
			{
				std::cout<<std::setprecision(16)<<"A = ["<<V<<"];"<<std::endl;
				std::cout<<"tempcond = cond(A);\ncondnumber = condnumber + tempcond;"<<std::endl;
				std::cout<<"k = k + 1;\nif(tempcond>maxcond)\nmaxcond = tempcond;\nend"<<std::endl;
				std::cout<<"if(tempcond<mincond)\nmincond = tempcond;\nend"<<std::endl;
				//std::cout<<"PHI: \n"<<phi<<std::endl;
				//std::cout<<"num neibs: "<<num_neibs_a<<std::endl;
				//std::cout<<"x: "<<xa[0]<<" "<<xa[1]<<std::endl;
				//std::cout<<"COEFF: "<<c<<std::endl;
				//std::cout<<"KAPPA: "<<curvature<<std::endl;
			}
			// vd.getProp<interpol_coeff>(a) = c;
			vd.getProp<min_sdf>(a) = min_sdf_val;
			vd.getProp<min_sdf_x>(a)[0] = min_sdf_val_x[0];
			vd.getProp<min_sdf_x>(a)[1] = min_sdf_val_x[1];

			//std::cout<<"Vmat\n"<<V<<std::endl;

			++part;
		}

}

template <typename CellList> inline void cp_optim(particles & vd, CellList & NN, particles_surface & vd_s)
{	// iterate over all particles, i.e. do closest point optimisation for all particles //
	auto NN_s = vd_s.getCellList(0.75*band_width);
	auto part = vd.getDomainIterator();
	vd.updateCellList(NN);
	int verbose = 0;
	while (part.isNext())
	{
		auto a = part.get();

		// initialise all variables
		EMatrix<double, Eigen::Dynamic, 1> xa(2,1);
		xa[0] = vd.getPos(a)[0];
		xa[1] = vd.getPos(a)[1];
		double val;
		EMatrix<double, Eigen::Dynamic, 1> x(2,1);
		EMatrix<double, Eigen::Dynamic, 1> dx(3, 1);
		dx[0] = 1.0;
		dx[1] = 1.0;
		dx[2] = 1.0;
		double lambda = 0.0;
		double norm_dx = dx.norm();

		double dpdx = 0;
		double dpdy = 0;
		double dpdxdx = 0;
		double dpdydy = 0;
		double dpdxdy = 0;
		double dpdydx = 0;

		EMatrix<double, Eigen::Dynamic, 1> c(16, 1);
		EMatrix<double, Eigen::Dynamic, 1> nabla_f(3, 1); // 3 = ndim + 1
		EMatrix<double, Eigen::Dynamic, Eigen::Dynamic> H(3, 3);

		if (vd.template getProp<surf_flag>(a) != 1)
		{
			if (return_sign(vd.getProp<sdf>(a)) < 0)
				{
					verbose = 1;
				//std::cout<<vd.getPos(a)[0]<<std::endl;
				//std::cout<<vd.getPos(a)[1]<<std::endl;}
				}

			Point<2,double> xaa = vd.getPos(a);
			//auto Np = NN_s.template getNNIterator<NO_CHECK>(NN.getCell(vd.getPos(a)));
			//auto Np = NN_s.template getNNIterator<NO_CHECK>(NN.getCell(vd.getPos(a)));
			auto Np = vd_s.getDomainIterator();
			vd_s.updateCellList(NN_s);
			double distance = 1000000.0;
			double dist_calc = 1000000000.0;
			decltype(a) b_min = a;
			//std::cout<<"huibuh2"<<std::endl;
			while (Np.isNext())
			{
				//std::cout<<"huibuh2.5"<<std::endl;
				auto b = Np.get();
				Point<2,double> xbb = vd_s.getPos(b);

				if (!vd.getProp<surf_flag>(vd_s.getProp<0>(b))) //todo: weirdly, the second cell list contains non-surface particles.
																//this needs to be investigated
				{
					++Np;
					continue;
				}
				dist_calc = xbb.distance(xaa);
				//if(verbose)std::cout<<dist_calc<<std::endl;
				if (dist_calc < distance)
				{
					distance = dist_calc;
					b_min = b;
				}
				++Np;
			}

			// set x0 to the particle closest to the surface
			//x[0] = vd.getPos(vd_s.getProp<0>(b_min))[0];
			//x[1] = vd.getPos(vd_s.getProp<0>(b_min))[1];
			val = vd.getProp<sdf>(vd_s.getProp<0>(b_min));
			x[0] = vd_s.getPos(b_min)[0];
			x[1] = vd_s.getPos(b_min)[1];
			// take the interpolation polynomial of the particle closest to the surface
			for (int k = 0 ; k < 16 ; k++)
			{
				//vd.getProp<interpol_coeff>(a)[k] = vd.getProp<interp_coeff>( vd_s.getProp<0>(b_min) )[k];
				c[k] = vd.getProp<interpol_coeff>(vd_s.getProp<0>(b_min))[k];
			}

			//++part;
			//continue;
		}
		else
		{
			x[0] = vd.getProp<min_sdf_x>(a)[0];
			x[1] = vd.getProp<min_sdf_x>(a)[1];
			val = vd.getProp<min_sdf>(a);

			for (int k = 0; k<16; k++)
				{
					c[k] = vd.getProp<interpol_coeff>(a)[k];
				}
		}
		if(verbose)
			{
			std::cout<<"VERBOSE%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%:"<<std::endl;
			std::cout<<"xa: "<<xa[0]<<", "<<xa[1]<<"\nx_0: "<<x[0]<<", "<<xa[1]<<"\nc: "<<c<<std::endl;
			}

		EMatrix<double, Eigen::Dynamic, 1> xax = x - xa;

		int k = 0;

		while(norm_dx > 0.000001)
		{
			// do optimisation //
			// gather required values //
			dpdx = get_dpdx(x, c);
			dpdy = get_dpdy(x, c);
			if (k == 0)
			{
				lambda = (-xax[0]*dpdx - xax[1]*dpdy)/(dpdx*dpdx + dpdy*dpdy);
			}
			dpdxdx = get_dpdxdx(x, c);
			dpdydy = get_dpdydy(x, c);
			dpdxdy = get_dpdxdy(x, c);
			dpdydx = dpdxdy;
			// Assemble gradient //
			nabla_f[0] = xax[0] + lambda*dpdx;
			nabla_f[1] = xax[1] + lambda*dpdy;
			nabla_f[2] = get_p(x, c);
			// Assemble Hessian matrix //
			H(0, 0) = 1 + lambda*dpdxdx;
			H(0, 1) = lambda*dpdydx;
			H(1, 0) = lambda*dpdxdy;
			H(1, 1) = 1 + lambda*dpdydy;
			H(0, 2) = dpdx;
			H(1, 2) = dpdy;
			H(2, 0) = dpdx;
			H(2, 1) = dpdy;
			H(2, 2) = 0;
			// compute and add increment // 
			//std::cout<<H.inverse()<<std::endl;
			dx = - H.inverse()*nabla_f;
			//std::cout<<dx<<std::endl;
			x[0] = x[0] + dx[0];
			x[1] = x[1] + dx[1];
			lambda = lambda + dx[2];
			// compute exit criterion and prepare xax
			xax = x - xa;
			norm_dx = dx.norm();

			//std::cout<<"H:\n"<<H<<"\nH_inv:\n"<<H_inv<<std::endl;	
			//std::cout<<"x:"<<x[0]<<"\nc:"<<std::endl;
			//std::cout<<c<<std::endl;
			//std::cout<<dpdx<<std::endl;
			++k;
			//std::cout<<k<<std::endl;
		}
//		debug optimisation
//		std::cout<<"p(x) = "<<get_p(x, c)<<std::endl;
//		std::cout<<"old sdf: "<<vd.getProp<sdf>(a)<<std::endl;
		vd.getProp<sdf>(a) = return_sign(vd.getProp<sdf>(a))*xax.norm();
		//std::cout<<"new sdf: "<<vd.getProp<sdf>(a)<<std::endl;

		//std::cout<<"c:\n"<<vd.getProp<interpol_coeff>(a)<<"\nmin_sdf: "<<vd.getProp<min_sdf>(a)<<" , min_sdf_x: "<<vd.getProp<min_sdf_x>(a)<<std::endl;

		++part;
	}

}

//################################################################# main function ##################
int main(int argc, char* argv[])
{
	openfpm_init(&argc, &argv);

	const double l = 4.0;
	Box<2, double> domain({-l/2.0, -l/2.0}, {l/2.0, l/2.0});
	size_t sz[2] = {(size_t)(l/dp + 0.5), (size_t)(l/dp + 0.5)};

	size_t bc[2] = {NON_PERIODIC, NON_PERIODIC};
	Ghost<2, double> g(0.0);

	particles vd(0, domain, bc, g, DEC_GRAN(512));
	particles_surface vd_s(vd.getDecomposition(),0);

	openfpm::vector<std::string> names({"sdf", "sdfgradient", "curvature", "surface_flag", "number_neibs", "min_sdf_neibors", "min_sdf_location", "interpol_coeff", "sdf_analytical"});
	vd.setPropNames(names);

	Box<2, double> particle_box({-l/2.0, -l/2.0}, {l/2.0, l/2.0});
	auto particle_it = DrawParticles::DrawBox(vd, sz, domain, particle_box);

	while (particle_it.isNext())
	{
		double x = particle_it.get().get(0);
		double y = particle_it.get().get(1);

		double dist = DistancePointEllipse(A, B, abs(x), abs(y));

		//std::cout<<x<<"\t"<<y<<"\t"<<dist<<std::endl;

		if (abs(dist) < band_width/2.0)
		{
			vd.add();
			vd.getLastPos()[0] = x;
			vd.getLastPos()[1] = y;
			vd.template getLastProp<sdf>() = 1 - sqrt((x/A)*(x/A) + (y/B)*(y/B));
			vd.template getLastProp<sdf_analytical>() = return_sign(1 - sqrt((x/A)*(x/A) + (y/B)*(y/B)))*dist;
			vd.template getLastProp<surf_flag>() = 0;
			vd.template getLastProp<num_neibs>() = 1;

		}

		++particle_it;

	}
//
	vd.map();

	auto NN = vd.getCellList(2*H);

	detect_surface_particles(vd, NN, vd_s);

	cp_interpol(vd, NN);

	//calc_surface_normals(vd, NN);
	std::cout<<"::::::::::::::::::::::::::::::::::::::before perturbing:"<<std::endl;
	calc_derivatives(vd, NN);
	//vd.write("init");

	perturb_pos(vd);
	vd.map();
	NN = vd.getCellList(2*H);

	update_sdfs(vd, NN);
	detect_surface_particles(vd, NN, vd_s);
	vd.write("after_perturb");
	//std::cout<<"begincode\ncondnumber = 0;\nk = 0;\nmincond = 1e99;\nmaxcond = 0;\n"<<std::endl;
	cp_interpol(vd, NN);

	cp_optim(vd, NN,vd_s);
	//detect_surface_particles(vd, NN);
	cp_interpol(vd, NN);
	std::cout<<"::::::::::::::::::::::::::::::::::::::after redistancing:"<<std::endl;
	calc_derivatives(vd, NN);

	vd.write("after_redistancing_w_sample");
	calc_derivatives(vd, NN);
	openfpm_finalize();

}
