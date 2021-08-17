#include "Vector/vector_dist.hpp"
#include <math.h>
#include  "Draw/DrawParticles.hpp"
#include "DCPSE/Dcpse.hpp"

#define PHASE_A 0
#define PHASE_B 1

const double dp = 1/64.0;
const double r = 1.0;
const double band_width = 12*dp;
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

typedef vector_dist<2, double, aggregate<int, double, double, double, double, double[2], double[2], int, int>> particles;
//										 |		|		|		|		|			|			|	   |    	|
//									  type    sdf  prev sdf redist_sdf distance sdfgrad prev_sdfgrad surf_flag num_neibs


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
			int sgn_b = return_sign(vd. template getProp<sdf>(b));
			Point<2,double> xb = vd.getPos(b);

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

template <typename CellList> inline void cp_redistance(particles & vd, CellList & NN)
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
		if (num_neibs_a > np_max)
		{
			throw std::domain_error("Particle has too many neighbours");
		}

		// Initialize phi-vector from the right hand side
		Point<np_max, double> phi_rhs = 0.0;

		Point<2, double> xa = vd.getPos(a);
		double min_sdf = std::abs(vd.getProp<sdf>(a));
		Point<2, double> min_sdf_x = xa;
		int neib = 0;

		MonomialBasis<2> m(3);

		EMatrix<double,Eigen::Dynamic,Eigen::Dynamic> V(num_neibs_a,m.size());
		EMatrix<double,Eigen::Dynamic,1> phi(num_neibs_a,1);
		EMatrix<double,Eigen::Dynamic,1> Atphi(m.size(),1);

		auto Np = NN.template getNNIterator<NO_CHECK>(NN.getCell(vd.getPos(a)));
		while(Np.isNext())
		{
			auto b = Np.get();
			Point<2,double> xb = vd.getPos(a);
			if (std::abs(vd.getProp<sdf>(b)) < min_sdf)
			{
				min_sdf = std::abs(vd.getProp<sdf>(b));
				min_sdf_x = xb;
			}

			// Fill phi-vector from the right hand side
			Atphi(neib,0) = vd.getProp<sdf>(b);

			vrb.buildRow(V,neib,xb,1.0);

			++neib;
			++Np;
		}

		EMatrix<double,Eigen::Dynamic,Eigen::Dynamic> AtA(m.size(),m.size());

		AtA = V.transpose()*V;
		Atphi = V.transpose()*phi;

		EMatrix<double,Eigen::Dynamic,1> c(m.size(),1);

		c = AtA.colPivHouseholderQr().solve(Atphi);

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
	openfpm::vector<std::string> names({"part_type", "previous_sdf", "sdf", "redistanced_sdf", "distance", "sdfgradient", "previous_sdfgradient", "surface_flag", "number_neibs"});
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

	calc_surface_normals(vd, NN);

	vd.write("init");

	perturb_pos(vd);

	vd.map();

	NN = vd.getCellList(2*H);
	calc_surface_normals(vd, NN);

	vd.write("after_perturb");

	detect_surface_particles(vd, NN);
	cp_redistance(vd, NN);

	vd.write("after_surf_detect");



	openfpm_finalize();

}
