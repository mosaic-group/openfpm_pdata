#include "Vector/vector_dist.hpp"
#include <math.h>
#include  "Draw/DrawParticles.hpp"
#include "DCPSE/MonomialBasis.hpp"
#include "DCPSE/Dcpse.hpp"

#define PHASE_A 0
#define PHASE_B 1

const double dp = 1/256.0;
const double r = 1.0;
const double band_width = 24.0*dp;
const int np_max = 30;

const int type = 0;
const int sdf = 1;
const int prev_sdf = 2;
const int redist_sdf = 3;
const int d = 4;
const int sdfgrad = 5;
const int prev_sdfgrad = 6;
const int surf_flag = 7;
const int num_neibs = 8;
const int min_sdf = 9;
const int min_sdf_x = 10;
const int interpol_coeff = 11;

typedef vector_dist<2, double, aggregate<int, double, double, double, double, double[2], double[2], int, int, double, double[2], EMatrix<double, Eigen::Dynamic, 1>>> particles;
//										 |		|		|		|		|			|			|	   |    	|			|						|				|
//									  type    sdf  prev sdf redist_sdf distance sdfgrad prev_sdfgrad surf_flag num_neibs min sdf in support	min sdf location poly coefficients


//double randZeroToOne()
//{
//    return rand() / (RAND_MAX + 1.);
//}

// SPH definitions for the approximation of the numerical gradient of the SDF

int return_sign(double phi)
{
	if (phi > 0) return 1;
	if (phi < 0) return -1;
	return 0;
}

const double H = 2.0*dp;
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
template<typename CellList> inline void calc_surface_normals(particles & vd, CellList & NN)
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

        int num_neibs_a = 1;
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
                ++num_neibs_a;
            }

            ++Np;
        }
        vd.getProp<num_neibs>(a) = num_neibs_a;
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


inline void perturb_pos(particles & vd)
{
	auto part = vd.getDomainIterator();

	while(part.isNext())
	{
		auto a = part.get();

		const double x = vd.getPos(a)[0];
		const double y = vd.getPos(a)[1];
		const double dist = vd.template getProp<d>(a);

		vd.getPos(a)[0] += x/r*0.5*dp;
		vd.getPos(a)[1] += y/r*0.5*dp;

		vd.template getProp<d>(a) = std::sqrt(vd.getPos(a)[0]*vd.getPos(a)[0] + vd.getPos(a)[1]*vd.getPos(a)[1]);

		++part;
	}
}

template <typename CellList> inline void detect_surface_particles(particles & vd, CellList & NN)
{
	auto part = vd.getDomainIterator();
	vd.updateCellList(NN);
	while (part.isNext())
	{
		auto a = part.get();
		int sgn_a = return_sign(vd. template getProp<sdf>(a));
		Point<2,double> xa = vd.getPos(a);

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
            	if (sgn_a != sgn_b)
            	{
            		vd.template getProp<surf_flag>(a) = 1;
            		break;
            	}
            }
            ++Np;
            
		}
		++part;
	}

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

template <typename CellList> inline void cp_interpol(particles & vd, CellList & NN)
{
	auto part = vd.getDomainIterator();
	vd.updateCellList(NN);
	while (part.isNext())
		{
			auto a = part.get();
			if (vd.template getProp<surf_flag>(a) != 1)
			{	
				//++part;
				//continue;
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
			EMatrix<double, Eigen::Dynamic, 1> Atphi(num_neibs_a, 1);

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
				phi(neib, 0) = vd.getProp<sdf>(b);
				//xb[0] = 0.1;
				//xb[1] = 0.5;
				vrb.buildRow(V, neib, xb, 1.0);

				++neib;
				++Np;
			}

			EMatrix<double, Eigen::Dynamic, Eigen::Dynamic> AtA(m.size(), m.size());
			// debug matrix A
			// std::cout<<"A:\n"<<V<<std::endl;

			AtA = V.transpose()*V;
			Atphi = V.transpose()*phi;

			EMatrix<double, Eigen::Dynamic, 1> c(m.size(), 1);

			c = AtA.colPivHouseholderQr().solve(Atphi);
			// debug matrix AtA, Atphi and coefficients
			// std::cout<<"AtA:\n"<<AtA<<"\nAtphi:\n"<<Atphi<<"\nc:\n"<<c<<std::endl;

			vd.getProp<interpol_coeff>(a) = c;
			vd.getProp<min_sdf>(a) = min_sdf_val;
			vd.getProp<min_sdf_x>(a)[0] = min_sdf_val_x[0];
			vd.getProp<min_sdf_x>(a)[1] = min_sdf_val_x[1];

			++part;
		}

}

template <typename CellList> inline void cp_optim(particles & vd, CellList & NN)
{
	auto part = vd.getDomainIterator();
	vd.updateCellList(NN);
	while (part.isNext())
	{
		auto a = part.get();
		if (vd.template getProp<surf_flag>(a) != 1)
		{	
			//++part;
			//continue;
		}

		EMatrix<double, Eigen::Dynamic, 1> xa(2,1);
		xa[0] = vd.getPos(a)[0];
		xa[1] = vd.getPos(a)[1];
		double val = vd.getProp<min_sdf>(a);
		EMatrix<double, Eigen::Dynamic, 1> x(2,1);
		x[0] = vd.getProp<min_sdf_x>(a)[0];
		x[1] = vd.getProp<min_sdf_x>(a)[1];
		EMatrix<double, Eigen::Dynamic, 1> dx(3, 1);
		dx[0] = 1.0;
		dx[1] = 1.0;
		dx[2] = 1.0;
		double lambda = 0.0;
		double norm_dx = dx.norm();
		EMatrix<double, Eigen::Dynamic, 1> xax = x - xa;

		double dpdx = 0;
		double dpdy = 0;
		double dpdxdx = 0;
		double dpdydy = 0;
		double dpdxdy = 0;
		double dpdydx = 0;

		EMatrix<double, Eigen::Dynamic, 1> c = vd.getProp<interpol_coeff>(a);
		EMatrix<double, Eigen::Dynamic, 1> nabla_f(3, 1); // 3 = ndim + 1
		EMatrix<double, Eigen::Dynamic, Eigen::Dynamic> H(3, 3);

		// order limit 4 corresponds to bicubic basis functions
		//MonomialBasis<2> m(4);
		//EMatrix<double, Eigen::Dynamic, 1> nabla_f(2, 1);
		// std::cout<<x[0]<<", "<<x[1]<<std::endl;

		int k = 0;

		//while(k<1)
		while(norm_dx > 0.000000000001)
		{
			// do optimisation
			// gather required values
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

			nabla_f[0] = xax[0] + lambda*dpdx;
			nabla_f[1] = xax[1] + lambda*dpdy;
			nabla_f[2] = get_p(x, c);

			H(0, 0) = 1 + lambda*dpdxdx;
			H(0, 1) = lambda*dpdydx;
			H(1, 0) = lambda*dpdxdy;
			H(1, 1) = 1 + lambda*dpdydy;
			H(0, 2) = dpdx;
			H(1, 2) = dpdy;
			H(2, 0) = dpdx;
			H(2, 1) = dpdy;
			H(2, 2) = 0;

			dx = - H.inverse()*nabla_f;
//			std::cout<<"hi"<<std::endl;
//			std::cout<<dx<<std::endl;
			x[0] = x[0] + dx[0];
			x[1] = x[1] + dx[1];
			lambda = lambda + dx[2];
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
//		std::cout<<"new sdf: "<<vd.getProp<sdf>(a)<<std::endl;

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
	openfpm::vector<std::string> names({"part_type", "previous_sdf", "sdf", "redistanced_sdf", "distance", "sdfgradient", "previous_sdfgradient", "surface_flag", "number_neibs", "min_sdf_neibors", "min_sdf_location"});
	vd.setPropNames(names);

	Box<2, double> particle_box({-l/2.0, -l/2.0}, {l/2.0, l/2.0});
	auto particle_it = DrawParticles::DrawBox(vd, sz, domain, particle_box);


	while (particle_it.isNext())
	{
		double x = particle_it.get().get(0);
		double y = particle_it.get().get(1);
		double dist = std::sqrt(x*x + y*y);

		if ((dist < (r + band_width/2.0)) && (dist > (r - band_width/2.0)))
		{
			vd.add();
			vd.getLastPos()[0] = x;
			vd.getLastPos()[1] = y;
			if (dist > r) vd.template getLastProp<type>() = PHASE_B;
			else vd.template getLastProp<type>() = PHASE_A;
			vd.template getLastProp<prev_sdf>() = r - dist;
			vd.template getLastProp<sdf>() = r - dist;
			vd.template getLastProp<d>() = dist;
			vd.template getLastProp<surf_flag>() = 0;
			vd.template getLastProp<num_neibs>() = 1;

		}

		++particle_it;

	}

	vd.map();

	auto NN = vd.getCellList(2*H);

	detect_surface_particles(vd, NN);
	
	calc_surface_normals(vd, NN);
	std::cout<<"#########################################################"<<std::endl;

	//vd.write("init");

	detect_surface_particles(vd, NN);
	//cp_interpol(vd, NN);

	perturb_pos(vd);

	vd.map();

	NN = vd.getCellList(2*H);
	calc_surface_normals(vd, NN);
	std::cout<<"#########################################################"<<std::endl;
	//vd.write("after_perturb");

	detect_surface_particles(vd, NN);
	cp_interpol(vd, NN);
	cp_optim(vd, NN);

	calc_surface_normals(vd, NN);
	std::cout<<"#########################################################"<<std::endl;
	//vd.write("after_surf_detect");

	openfpm_finalize();

}