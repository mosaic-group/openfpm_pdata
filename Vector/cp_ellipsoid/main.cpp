#define SE_CLASS1

#include <math.h>
#include <sys/_types/_size_t.h>
#include "Vector/vector_dist.hpp"
#include "DCPSE/Dcpse.hpp"
#include "DCPSE/MonomialBasis.hpp"
#include "Draw/DrawParticles.hpp"
#include "../../openfpm_numerics/src/level_set/particle_cp/particle_cp.hpp"
#include "Operators/Vector/vector_dist_operators.hpp"
#include <chrono>

const double dp = 1.0/64.0;
const double A = 0.75;//0.5;
const double B = 0.5;//1.0/3.0;
const double C = 0.5;
const double band_width = 12.0*dp;

const int sdf = 0;
const int sdfgrad = 1;
const int curvature = 2;
const int cp = 3;
const int surf_flag = 4;
const int sdf_analytical = 5;
const int cp_theoretical = 6;

typedef vector_dist<3, double, aggregate<double, double[3], double, double[3], int, double, double[3]>> particles;
//					 |	 |	    |		|       |	|	|
//				     sdf  sdfgrad  curvature 		cp surf_flag theoretical sdf	theoretical cp

//typedef vector_dist<2, double, aggregate<vect_dist_key_dx>> particles_surface;

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

double GetRoot(double r0, double r1, double z0, double z1, double z2, double g)
{	
	const int maxIter = 100;
	double n0 = r0*z0;
	double n1 = r1*z1;
	double s0 = z2 - 1;
	double s1 = (g < 0 ? 0 : sqrt(n0*n0 + n1*n1 + z2*z2) - 1) ;
	double s = s;
	for(int i = 0 ; i < maxIter ; ++i )
	{
		s = ( s0 + s1 ) / 2 ;
		if ( s == s0 || s == s1 ) {break; }
		double ratio0 = n0 / ( s + r0 ); 
		double ratio1 = n1 / ( s + r1 );
		double ratio2 = z2 / ( s + 1 );
		g = ratio0*ratio0 + ratio1*ratio1 +ratio2*ratio2 - 1;
		if ( g > 0 ) { s0 = s ;} 
		else if ( g < 0 ) { s1 = s ; }
		else {break;}
	}
	return (s);
}

double DistancePointEllipse(double e0, double e1, double y0, double y1, double& x0, double& x1)
    {
        double distance;
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

double DistancePointEllipsoid(double e0, double e1, double e2, double y0, double y1, double y2, double& x0, double& x1, double& x2)
{
	double distance;
	if( y2 > 0 )
	{
		if( y1 > 0 )
		{
			if( y0 > 0 )
			{
				double z0 = y0 / e0;
				double z1 = y1 / e1;
				double z2 = y2 / e2;
				double g = z0*z0 + z1*z1 + z2*z2 - 1 ;
				if( g != 0 )
				{
					double r0 = (e0/e2)*(e0/e2);
					double r1 = (e1/e2)*(e1/e2);
					double sbar = GetRoot ( r0 , r1 , z0 , z1 , z2 , g );
					x0 = r0 *y0 / ( sbar + r0 );
					x1 = r1 *y1 / ( sbar + r1 );
					x2 = y2 / ( sbar + 1 );
					distance = sqrt( (x0 - y0)*(x0 - y0) + (x1 - y1)*(x1 - y1) + (x2 - y2)*(x2 - y2));
				}
				else
				{
					x0 = y0;
					x1 = y1;
					x2 = y2;
					distance = 0;
				}
			}
			else // y0 == 0
			{
				x0 = 0;
				distance = DistancePointEllipse( e1 , e2 , y1 , y2, x1, x2);
			}
		}
		else // y1 == 0
		{
			if( y0 > 0 )
			{
				x1 = 0;
				distance = DistancePointEllipse( e0 , e2 , y0 , y2, x0, x2);
			}
			else // y0 == 0
			{
				x0 = 0;
				x1 = 0;
				x2 = e2;
				distance = abs(y2 - e2);
			}
		}
	}
	else // y2 == 0
	{
		double denom0 = e0*e0 - e2*e2;
		double denom1 = e1*e1 - e2*e2;
		double numer0 = e0*y0;
		double numer1 = e1*y1;
		bool computed = false;
		if((numer0 < denom0) && (numer1 < denom1))
		{
			double xde0 = numer0/denom0;
			double xde1 = numer1/denom1 ;
			double xde0sqr = xde0 *xde0;
			double xde1sqr = xde1 * xde1 ;
			double discr = 1 - xde0sqr - xde1sqr;
			if( discr > 0 )
			{
				x0 = e0*xde0;
				x1 = e1*xde1;
				x2 = e2*sqrt(discr);
				distance = sqrt((x0 - y0)*(x0 - y0) + (x1 - y1)*(x1 - y1) + x2*x2);
				computed = true;
			}
		}
		if( !computed )
		{
			x2 = 0;
			distance = DistancePointEllipse(e0 , e1 , y0 , y1, x0, x1);
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
		const double z = vd.getPos(a)[2];
		const double dist = sqrt(x*x + y*y + z*z);

		vd.getPos(a)[0] += randZeroToOne()*dp/3.0;
		vd.getPos(a)[1] += randZeroToOne()*dp/3.0;
		vd.getPos(a)[2] += randZeroToOne()*dp/3.0;
		//vd.template getProp<d>(a) = std::sqrt(vd.getPos(a)[0]*vd.getPos(a)[0] + vd.getPos(a)[1]*vd.getPos(a)[1]);

		++part;
	}
}

template <typename CellList> inline void update_sdfs(particles & vd, CellList & NN)//
{
	auto part = vd.getDomainIterator();
	vd.updateCellList(NN);
	while (part.isNext())
	{
		auto a = part.get();
		Point<3, double> xa = vd.getPos(a);
		double x = xa[0];
		double y = xa[1];
		double z = xa[2];
		vd.getProp<sdf>(a) = - 1 + sqrt((x/A)*(x/A) + (y/B)*(y/B) + (z/C)*(z/C));

		// binary distribution
		//if ((- 1 + sqrt((x/A)*(x/A) + (y/B)*(y/B))) > 0) vd.getProp<sdf>(a) = dp;
		//else vd.getProp<sdf>(a) = -dp;

		// Saye's initial configuration
		//vd.getProp<sdf>(a) = (1 - exp(-(x-0.3)*(x-0.3)-(y-0.3)*(y-0.3)))*(sqrt(4*x*x+9*y*y)-1);

		double xcp_analytical;
		double ycp_analytical;
		double zcp_analytical;
		vd.getProp<sdf_analytical>(a) = return_sign(vd.getProp<sdf>(a))*DistancePointEllipsoid(A, B, C, abs(x), abs(y), abs(z), xcp_analytical, ycp_analytical, zcp_analytical);

		vd.template getProp<cp_theoretical>(a)[0] = return_sign(x)*xcp_analytical;
		vd.template getProp<cp_theoretical>(a)[1] = return_sign(y)*ycp_analytical;
		vd.template getProp<cp_theoretical>(a)[2] = return_sign(z)*zcp_analytical;

		++part;
	}
}

inline void get_max_error(particles & vd)
{
	auto part = vd.getDomainIterator();
	double err;
	double maxerr = 0.0;
	double errkappa;
	double maxerrkappa = 0.0;
	double errnorm;
	double maxerrnorm = 0.0;
	double errcp;
	double maxerrcp = 0.0;
	double errone;
	double maxerrone = 0.0;

	while(part.isNext())
	{
		auto a = part.get();
		// L infinity norm:
		//errnorm = abs(nxa[0]);
		//if (abs(nxa[1]) > errnorm) errnorm = abs(nxa[1]);
		//if (abs(nxa[2]) > errnorm) errnorm = abs(nxa[2]);
		double xcp_analytical = vd.getProp<cp_theoretical>(a)[0];
		double ycp_analytical = vd.getProp<cp_theoretical>(a)[1];
		double zcp_analytical = vd.getProp<cp_theoretical>(a)[2];
		Point<3, double> norm_analytical;
		Point<3, double> norm_computed;
		for (int k = 0; k < 3; k++) norm_computed[k] = vd.getProp<sdfgrad>(a)[k];
		Point<3, double> cp_computed;
		Point<3, double> cp_theo;
		for (int k = 0; k < 3; k++)
		{
			cp_computed[k] = vd.getProp<cp>(a)[k];
			cp_theo[k] = vd.getProp<cp_theoretical>(a)[k];
		}

		norm_analytical[0] = 2*xcp_analytical/(A*A);
		norm_analytical[1] = 2*ycp_analytical/(B*B);
		norm_analytical[2] = 2*zcp_analytical/(C*C);
		double norm_norm_analytical = norm(norm_analytical);
		for(int k = 0; k < 3; k++) norm_analytical[k] = return_sign(vd.getProp<sdf>(a))*norm_analytical[k]/norm_norm_analytical;
		//std::cout<<norm(norm_analytical)<<std::endl;
		errnorm = norm(norm_analytical - norm_computed);

		err = abs(vd.getProp<sdf>(a) - vd.getProp<sdf_analytical>(a));
		errcp = norm(cp_computed - cp_theo);

		double kappa = (std::abs(xcp_analytical*xcp_analytical + ycp_analytical*ycp_analytical + zcp_analytical*zcp_analytical - A*A - B*B - C*C))/(2*A*A*B*B*C*C*std::pow(xcp_analytical*xcp_analytical/std::pow(A, 4) + ycp_analytical*ycp_analytical/std::pow(B, 4) + zcp_analytical*zcp_analytical/std::pow(C, 4), 1.5));
		errkappa = abs(vd.getProp<curvature>(a) - kappa);
		if (err > maxerr) maxerr = err;
		if (errkappa > maxerrkappa) maxerrkappa = errkappa;
		if (errnorm > maxerrnorm) maxerrnorm = errnorm;
		if (errcp > maxerrcp) maxerrcp = errcp;
		++part;
	}
	std::cout<<"Maximum error for sdf on processor is: "<<maxerr<<std::endl;
	std::cout<<"Maximum error for surface normal on processor is: "<<maxerrnorm<<std::endl;
	std::cout<<"Maximum error for curvature on processor is: "<<maxerrkappa<<std::endl;
	std::cout<<"Maximum error for cp(xa) is: "<<maxerrcp<<std::endl;
	std::cout<<"Maximum deviation of unit normal length from 1 is: "<<maxerrone<<std::endl;

}

//################################################################# main function ##################
int main(int argc, char* argv[])
{
	openfpm_init(&argc, &argv);

	const double l = 2.0;
	Box<3, double> domain({-l/2.0, -l/3.0, -l/3.0}, {l/2.0, l/3.0, l/3.0});
	size_t sz[3] = {(size_t)(l/dp + 0.5), (size_t)((2.0/3.0)*l/dp + 0.5), (size_t)((2.0/3.0)*l/dp + 0.5)};

	size_t bc[3] = {NON_PERIODIC, NON_PERIODIC, NON_PERIODIC};
	Ghost<3, double> g(band_width);

	particles vd(0, domain, bc, g, DEC_GRAN(512));
	//particles_surface vd_s(vd.getDecomposition(), 0);

	openfpm::vector<std::string> names({"sdf", "sdfgradient", "curvature", "cp", "surface_flag", "sdf_theoretical", "cp_theoretical"});
	vd.setPropNames(names);
	Box<3, double> particle_box({-l/2.0, -l/3.0, -l/3.0}, {l/2.0, l/3.0, l/3.0});

    int count = 0;
	auto particle_it = DrawParticles::DrawBox(vd, sz, domain, particle_box);

	while (particle_it.isNext())
	{
		double x = particle_it.get().get(0);
		double y = particle_it.get().get(1);
		double z = particle_it.get().get(2);
		double xcp_analytical;
		double ycp_analytical;
		double zcp_analytical;

		double dist = DistancePointEllipsoid(A, B, C, abs(x), abs(y), abs(z), xcp_analytical, ycp_analytical, zcp_analytical);

		//std::cout<<x<<"\t"<<y<<"\t"<<dist<<std::endl;

		if (abs(dist) < band_width/2.0)
		{
			vd.add();
			vd.getLastPos()[0] = x;
			vd.getLastPos()[1] = y;
			vd.getLastPos()[2] = z;
			vd.template getLastProp<sdf>() = 1 - sqrt((x/A)*(x/A) + (y/B)*(y/B) + (z/C)*(z/C));
			vd.template getLastProp<sdf_analytical>() = return_sign(1 - sqrt((x/A)*(x/A) + (y/B)*(y/B) + (z/C)*(z/C)))*dist;
			vd.template getLastProp<cp_theoretical>()[0] = return_sign(x)*xcp_analytical;
			vd.template getLastProp<cp_theoretical>()[1] = return_sign(y)*ycp_analytical;
			vd.template getLastProp<cp_theoretical>()[2] = return_sign(z)*zcp_analytical;
			vd.template getLastProp<surf_flag>() = 0;
			count++;
		}
		++particle_it;
	}
    std::cout<<"Total number of particles created: "<<count<<std::endl;
	vd.map();

	///

	auto NN = vd.getCellList(2*H);

	if(true)
	{
		perturb_pos(vd);
		vd.map();

		//vd.getDecomposition().decompose();
		//vd.map();
		NN = vd.getCellList(2*H);
		update_sdfs(vd, NN);

//		for debugging parallel stuff potentially (note that this set has dp = 1/32.
//		vd.load("particles_set");
//		vd.map();

		vd.write("before_redistancing_binary");
		Redist_options rdistoptions;
		//rdistoptions.max_iter = 1000;
		//rdistoptions.incremental_tolerance = 1e-7;
		rdistoptions.H = dp;
		rdistoptions.r_cutoff_factor = 2.4;
		rdistoptions.sampling_radius = 0.75*band_width;
		//rdistoptions.barrier_coefficient = 0.0;
		//rdistoptions.armijo_tau = 0.9;
		//rdistoptions.armijo_c = 0.9;
		rdistoptions.tolerance = 1e-13;//dp*dp*dp*dp*dp;
		//rdistoptions.support_prevent = 0.9;
		//rdistoptions.init_project = 1;
		rdistoptions.write_cp = 1;
		rdistoptions.compute_normals = 1;
		rdistoptions.compute_curvatures = 1;
		rdistoptions.minter_poly_degree = 4;
		rdistoptions.minter_lp_degree = 1.0;
        	static constexpr unsigned int num_coeffs = minter_lp_degree_one_num_coeffs(3,4);
	        std::cout<<"number of coefficients determined as "<<num_coeffs<<std::endl;
		//typedef vector_dist_ws<2, double, particles> particles_in_type;
		particle_cp_redistancing<particles, taylor4, num_coeffs> pcprdist(vd, rdistoptions);

		std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

		std::cout<<"starting redistancing"<<std::endl;
		pcprdist.run_redistancing();

		std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
		std::cout << "Time difference for pcp redistancing = " << std::chrono::duration_cast<std::chrono::seconds>(end - begin).count() << "[s]" << std::endl;

		vd.ghost_get<sdf,sdf_analytical>();
		get_max_error(vd);
		vd.write("after_redistancing_after_cleanup_minterdebug");


//		particle_cp_redistancing pcprdist2(vd, rdistoptions);
//
//		pcprdist2.run_redistancing();
//		rdistoptions.polynomial_degree = "taylor4";
//		rdistoptions.r_cutoff_factor = 2.5;
//
//		particle_cp_redistancing pcprdist3(vd, rdistoptions);
//
//		pcprdist3.run_redistancing();
//		particle_cp_redistancing pcprdist4(vd, rdistoptions);
//		pcprdist4.run_redistancing();
//
//		particle_cp_redistancing pcprdist5(vd, rdistoptions);
//		pcprdist5.run_redistancing();
//
//		vd.write("after_redistancing_after_cleanup_binary");

//
//		particle_cp_redistancing pcprdist5(vd, rdistoptions);
//		pcprdist5.run_redistancing();

	}

	openfpm_finalize();
}

