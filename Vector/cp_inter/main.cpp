#include "Vector/vector_dist.hpp"
#include <math.h>
#include  "Draw/DrawParticles.hpp"

#define PHASE_A 0
#define PHASE_B 1

const double dp = 1/64.0;
const double r = 1.0;
const double band_width = 12*dp;

const size_t type = 0;
const int sdf = 1;
const int prev_sdf = 2;
const int redist_sdf = 3;
const int d = 4;
const int sdfgrad = 5;
const int prev_sdfgrad = 6;

typedef vector_dist<2, double, aggregate<int, double, double, double, double, double[2], double[2]>> particles;
//										  |		|		|		|		|			|			|
//										type  prev sdf	sdf	redist_sdf	distance	sdfgrad		prev_sdfgrad


//double randZeroToOne()
//{
//    return rand() / (RAND_MAX + 1.);
//}

// SPH definitions for the approximation of the numerical gradient of the SDF

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

template<typename CellList> inline void calc_surface_normals(particles & vd, CellList & NN)
{
    auto part = vd.getDomainIterator();

    // Update the cell-list
    vd.updateCellList(NN);

    // For each particle ...
    while (part.isNext())
    {
        // ... a
        auto a = part.get();
        //reset the values
        vd.getProp<sdfgrad>(a)[0]=0.0;
        vd.getProp<sdfgrad>(a)[1]=0.0;
        // Get the position xp of the particle
        Point<2,double> xa = vd.getPos(a);

        // Take the mass of the particle dependently if it is FLUID or BOUNDARY
//        double massa = 0;
//        if (vd.getProp<type>(a) == FLUID) massa = MassFluid;
//        else if(vd.getProp<type>(a) == FLUID_B) massa = MassFluid_B;
//        else if(vd.getProp<type>(a) == BOUNDARY) massa = MassBound;
//
//        // Get the density of the of the particle a
//        double rhoa = vd.getProp<rho>(a);
        // For the purpose of only computing the gradient, say V_a = H^3, and rho_a=1000, hence
        //double rhoa = 1000.0;
        //double massa = H*H*H*rhoa;

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

//            double massb = MassFluid;
//            double rhob = vd.getProp<rho>(b);
            // same with V_b
            //double massb = massa;
            //double rhob = rhoa;

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
        ++part;
    }
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

		vd.getPos(a)[0] += x/r*dp;
		vd.getPos(a)[1] += y/r*dp;

		vd.template getProp<d>(a) = std::sqrt(vd.getPos(a)[0]*vd.getPos(a)[0] + vd.getPos(a)[1]*vd.getPos(a)[1]);

		++part;
	}
}


int main(int argc, char* argv[])
{
	openfpm_init(&argc, &argv);

	const double l = 4.0;
	Box<2, double> domain({-l/2.0, -l/2.0}, {l/2.0, l/2.0});
	size_t sz[2] = {(size_t)(l/dp + 0.5), (size_t)(l/dp + 0.5)};

	size_t bc[2] = {NON_PERIODIC, NON_PERIODIC};
	Ghost<2, double> g(0.0);

	particles vd(0, domain, bc, g, DEC_GRAN(512));
	openfpm::vector<std::string> names({"part_type", "previous_sdf", "sdf", "redistanced_sdf", "distance", "sdfgradient", "previous_sdfgradient"});
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

		}

		++particle_it;


	}

	vd.map();

	auto NN = vd.getCellList(2*H);

	calc_surface_normals(vd, NN);

	vd.write("init");

	perturb_pos(vd);

	vd.map();

	NN = vd.getCellList(2*H);
	calc_surface_normals(vd, NN);

	vd.write("init1");

	openfpm_finalize();

}
