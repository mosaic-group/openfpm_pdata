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

const int sdf = 0;
const int normal = 1;
const int curvature = 2;
const int surf_flag = 3;
const int rho = 4;
const int p = 5;
const int force = 6;
const int vel = 7;
const int vel_prev = 8;
const int pos_prev = 9;

typedef vector_dist<2, double, aggregate<double, double[2], double, int, double, double, Point<2, double>, Point<2, double>, Point<2, double>, Point<2, double>>> particles;
//											|		|	 	 |         |    |		 |		|	|						|
//										     sdf  sdfgrad curvature surf_flag 	density, previous density	pressure viscosity (dyn), force, velocity, previous velocity

//typedef vector_dist<2, double, aggregate<vect_dist_key_dx>> particles_surface;

const double A = 0.75;
const double B = 0.5;
const double dp = 1/128.0;
//const double H = std::sqrt(3.0*dp*dp);
const double H = 1.3*dp;
const double gamma_eos = 1.4;
const double c = 200.0;
const double rho_0 = 1.0;
const double rho_1 = 1.0;
const double eta_phase0 = 0.2;
const double eta_phase1 = 0.2;
const double alpha = 50.0;
const double p_background = 0.0;

const double dt = 2.5*1e-5;
double t;
int corr;
const int colorfield = 1;

const double band_width = 8.0*dp;

const double l = 2.44;
const double t_end = 1.0;
const double M_0 = A*B*M_PI*rho_0;
const double M_1 = (l*l - A*B*M_PI)*rho_1;
double np_0;
double np_1;
double m_0;
double m_1;

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

double randZeroToOne()
{
    return (rand() / (RAND_MAX + 1.));//
}

int return_sign(double phi)
{
	if (phi > 0) return 1;
	if (phi < 0) return -1;
	return 0;
}

// kernel function

const double a2 = 7.0/(4.0*M_PI*H*H);
inline double Wab(double r)
{
    const double q = r/H;
    if (q <= 2.0)
    {
        double factor = 1.0 - q/2.0;
        factor = factor*factor;
        factor = factor*factor;
        return(a2*factor*(1.0 + 2.0*q));
    }
    else return(0.0);
}

const double a3 = 5.0/(8.0*H);
inline double Wab1D(double r)
{
    const double q = r/H;
    if (q <= 2.0)
    {
        double factor = 1.0 - q/2.0;
        factor = factor*factor*factor;
        return(a3*factor*(1.5*q + 1));
    }
    else return(0.0);
}

inline void DWab(Point<2,double> & dx, Point<2,double> & DW, double r, bool print, double & dwdrab)
{
    const double qq=r/H;
    if (qq <= 2.0)
    {
        double factor = (-5.0*a2/(H))*qq*(1.0 - qq/2.0)*(1.0 - qq/2.0)*(1.0 - qq/2.0);

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

// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SPH geometric computing part %%%%%%%%%%%%%%%%%%%%%
const double eps_normal = 0.01/H;
inline double coloraverage(double rhoa, double rhob, int typea, int typeb)
{
    double cij;
    if (typea==typeb)
    {cij = 0.0;}
    else cij = 1.0;
    return (rhoa/(rhoa+rhob)*cij);
}

template<typename CellList> inline void calc_surface_normals(particles & vd, CellList & NN)
{
    vd.template ghost_get<rho, sdf>();
    auto part = vd.getDomainIterator();

    // Update the cell-list
    vd.updateCellList(NN);

    // For each particle ...
    while (part.isNext())
    {
        // ... a
        auto a = part.get();
        vd.getProp<normal>(a)[0] = 0.0;
        vd.getProp<normal>(a)[1] = 0.0;
        //if (std::abs(vd.getProp<sdf>(a))>2.0*H) {++part; continue;};
        // Get the position xp of the particle
        Point<2,double> xa = vd.getPos(a);

        double massa = (return_sign(vd.getProp<sdf>(a)) < 0)?m_0:m_1;
        double ca = return_sign(vd.getProp<sdf>(a));
        // Get the density of the particle a
        double rhoa = vd.getProp<rho>(a);

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

            double massb = (return_sign(vd.getProp<sdf>(b)) < 0)?m_0:m_1;

            double rhob = vd.getProp<rho>(b);

            // Get the distance between p and q
            Point<2,double> dr = xa - xb;
            // take the norm of this vector
            double r2 = norm2(dr);

            // if they interact
            if (r2 < 4.0*H*H)
            {
                double r = sqrt(r2); // norm2 is norm^2

                Point<2,double> DW;
                double dwdrab;
                DWab(dr,DW,r,false,dwdrab);
                double cfactor = 0.0;
                double cb = return_sign(vd.getProp<sdf>(b));
                if (colorfield) {
                    //cfactor = rhoa/massa*((massa/rhoa)*(massa/rhoa)+(massb/rhob)*(massb/rhob))*coloraverage(rhoa,rhob,return_sign(vd.getProp<sdf>(a)),return_sign(vd.getProp<sdf>(b)));
                    cfactor = massb/rhob*(cb - ca);
                }
                else cfactor = massb/rhob*(vd.getProp<sdf>(b) - vd.getProp<sdf>(a));
                vd.getProp<normal>(a)[0] += cfactor * DW.get(0);
                vd.getProp<normal>(a)[1] += cfactor * DW.get(1);
            }
            ++Np;
        }
        //std::cout<<"normal:\n"<<vd.getProp<normal>(a)[0]<<", "<<vd.getProp<normal>(a)[1]<<std::endl;
        // normalize normal to obtain unit surface normal
        if (false) {
            double colornorm = sqrt(vd.getProp<normal>(a)[0] * vd.getProp<normal>(a)[0] +
                                    vd.getProp<normal>(a)[1] * vd.getProp<normal>(a)[1]);
            if (colornorm < eps_normal) {
                vd.getProp<normal>(a)[0] = 0.0;
                vd.getProp<normal>(a)[0] = 0.0;
            } else {
                vd.getProp<normal>(a)[0] = vd.getProp<normal>(a)[0] / colornorm;
                vd.getProp<normal>(a)[1] = vd.getProp<normal>(a)[1] / colornorm;
            }
        }

        if (colorfield){
            if (return_sign(vd.getProp<sdf>(a))>0) {
                //vd.getProp<normal>(a)[0] = -1.0 * vd.getProp<normal>(a)[0];
                //vd.getProp<normal>(a)[1] = -1.0 * vd.getProp<normal>(a)[1];
            }}
        ++part;
    }
}

template<typename CellList> inline void calc_curvature(particles & vd, CellList & NN)
{
    vd.template ghost_get<rho, sdf, normal>();
    auto part = vd.getDomainIterator();

    // Update the cell-list
    vd.updateCellList(NN);

    // For each particle ...
    while (part.isNext())
    {
        // ... a
        auto a = part.get();
        //if (std::abs(vd.getProp<sdf>(a))>2.0*H) {++part; continue;};
        vd.getProp<curvature>(a) = 0.0;

        // Get the position xp of the particle
        Point<2,double> xa = vd.getPos(a);
        // Get the (unity) surface normal
        Point<2,double> na = vd.getProp<normal>(a);

        double norma = sqrt(na.get(0)*na.get(0)+na.get(1)*na.get(1));
        if (norma<eps_normal) norma = 0.0;
        if (norma<1e-10) {++part; continue;};
        na[0] = na[0]/norma;
        na[1] = na[1]/norma;

        // Take the mass of the particle dependently
        double massa = (return_sign(vd.getProp<sdf>(a)) < 0)?m_0:m_1;

        // Get the density of the of the particle a
        double rhoa = vd.getProp<rho>(a);
        double nominator = 0.0;
        double denominator = 0.0;
        double support = Wab(0.0)*massa/rhoa;

        // Get an iterator over the neighborhood particles of p
        auto Np = NN.template getNNIterator<NO_CHECK>(NN.getCell(vd.getPos(a)));

        // For each neighborhood particle
        while (Np.isNext() == true)
        {
            // ... q
            auto b = Np.get();

            // Get the position xp of the particle
            Point<2,double> xb = vd.getPos(b);
            //Get the surface normal of the particle
            Point<2,double> nb = vd.getProp<normal>(b);
            double normb = sqrt(nb.get(0)*nb.get(0)+nb.get(1)*nb.get(1));
            if (norma<eps_normal) norma = 0.0;
            if (normb<1e-10) {++Np; continue;};
            nb[0] = nb[0]/normb;
            nb[1] = nb[1]/normb;

            // if (p == q) skip this particle
            if (a.getKey() == b)	{++Np; continue;};

            double massb = (return_sign(vd.getProp<sdf>(b)) < 0)?m_0:m_1;
            double rhob = vd.getProp<rho>(b);

            // Get the distance between p and q
            Point<2,double> dr = xa - xb;
            // take the norm of this vector
            double r2 = norm2(dr);

            // if they interact
            if (r2 < 4.0*H*H)
            {
                double r = sqrt(r2); // norm2 is norm^2
                Point<2,double> nab = na - nb;

                Point<2,double> DW;
                double dwdrab;
                DWab(dr,DW,r,false,dwdrab);

                nominator += (nab.get(0)*DW.get(0) + nab.get(1)*DW.get(1))*massb/rhob;
                denominator += r*dwdrab*massb/rhob;

                vd.getProp<curvature>(a) += (-nab.get(0)*DW.get(0) + -nab.get(1)*DW.get(1))*massb/rhob;
                support += Wab(r)*massb/rhob;
            }
            ++Np;
        }
        // 2 comes from dimension
        //std::cout<<"nominator:"<<nominator<<std::endl;
        //std::cout<<"denominator:"<<denominator<<std::endl;
        //vd.getProp<curvature>(a) = 2.0*nominator/denominator;
        vd.getProp<curvature>(a) = vd.getProp<curvature>(a)/support;
        ++part;
        //std::cout<<"curvature:"<<vd.getProp<curvature>(a)<<std::endl;
        //std::cout<<"normal: "<<vd.getProp<normal>(a)[0]<<", "<<vd.getProp<normal>(a)[1]<<std::endl;
    }
}
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% End SPH geometric computing part

template <typename CellList> inline void density_summation(particles & vd, CellList & NN)
{
    vd.template ghost_get<rho>();
    auto part = vd.getDomainIterator();
    vd.updateCellList(NN);
    while (part.isNext())
    {
        auto a = part.get();
        Point<2, double> xa = vd.getPos(a);
        double rho_a = 0.0;
        double m = (return_sign(vd.getProp<sdf>(a)) < 0)?m_0:m_1;

        auto Np = NN.template getNNIterator<NO_CHECK>(NN.getCell(vd.getPos(a)));

        while (Np.isNext() == true)
        {
            auto b = Np.get();
            Point<2, double> xb = vd.getPos(b);
            Point<2, double> dr = xa - xb;
            double r2 = norm2(dr);

            if (r2 <= 4.0*H*H)
            {
                double r = std::sqrt(r2);
                rho_a += Wab(r);
            }
            ++Np;
        }
        vd.getProp<rho>(a) = m*rho_a;
        ++part;
    }
}

inline void EqOfState(particles & vd)
{
    auto it = vd.getDomainIterator();
    while (it.isNext())
    {
        auto a = it.get();
        double ref_dens = 0.0;
        if (return_sign(vd.getProp<sdf>(a)) < 0) ref_dens = rho_0;
        else ref_dens = rho_1;
        double dens =  vd.getProp<rho>(a);
        vd.getProp<p>(a) = (c*c*ref_dens/gamma_eos)*(std::pow(dens/ref_dens, gamma_eos) - 1) + p_background;

        ++it;
    }
}

template<typename CellList> inline void calc_forces(particles & vd, CellList & NN)
{
    vd.template ghost_get<rho, p, vel>();
    auto part = vd.getDomainIterator();

    // Update the cell-list
    vd.updateCellList(NN);

    double max_p = 0.0;
    double max_eta = 0.0;
    double avg_alpha = 0.0;
    double maxvel = 0.0;
    int numparticles = 0;
    int numsurfparticles = 0;

    // For each particle ...
    while (part.isNext())
    {
        auto a = part.get();

        vd.getProp<force>(a)[0] = 0.0;
        vd.getProp<force>(a)[1] = 0.0;

        Point<2, double> p_force;
        Point<2, double> eta_force;
        for(int k;k<2;k++) p_force[k] = 0.0;
        for(int k;k<2;k++) eta_force[k] = 0.0;

        // Get the position xp of the particle
        Point<2, double> xa = vd.getPos(a);

        // Take the mass of the particle dependently if it is FLUID 0 or 1
        double massa = (return_sign(vd.getProp<sdf>(a)) < 0)?m_0:m_1;
        double etaa = (return_sign(vd.getProp<sdf>(a)) < 0)?eta_phase0:eta_phase1;
        // Get the density of the of the particle a
        double rhoa = vd.getProp<rho>(a);
        double Va = massa/rhoa;

        // Get the pressure of the particle a
        double Pa = vd.getProp<p>(a);

        // Get the Velocity of the particle a
        Point<2, double> va = vd.getProp<vel>(a);

        // Get an iterator over the neighborhood particles of p
        auto Np = NN.template getNNIterator<NO_CHECK>(NN.getCell(vd.getPos(a)));

        // For each neighborhood particle
        while (Np.isNext() == true)
        {
            auto b = Np.get();

            // Get the position xp of the particle
            Point<2, double> xb = vd.getPos(b);

            // if (p == q) skip this particle
            if (a.getKey() == b)	{++Np; continue;};

            double massb = (return_sign(vd.getProp<sdf>(a)) < 0)?m_0:m_1;
            double etab = (return_sign(vd.getProp<sdf>(a)) < 0)?eta_phase0:eta_phase1;
            Point<2, double> vb = vd.getProp<vel>(b);
            double Pb = vd.getProp<p>(b);
            double rhob = vd.getProp<rho>(b);
            double Vb = massb/rhob;

            // Get the distance between p and q
            Point<2, double> dr = xa - xb;
            // take the norm of this vector
            double r2 = norm2(dr);

            // if they interact
            if (r2 < 4.0*H*H)
            {
                double r = sqrt(r2);

                Point<2, double> v_rel = va - vb;

                Point<2, double> DW;
                double dwdrab;
                DWab(dr,DW,r,false,dwdrab);

                double factor_p = - (Va*Va + Vb*Vb)*(rhoa*Pb + rhob*Pa)/(rhoa + rhob)/massa;
                //double factor_p = -(Pa/(rhoa*rhoa) + Pb/(rhob*rhob))*massb;
                //double factor_p = - rhoa*massa*(Pa/(rhoa*rhoa) + Pb/(rhob*rhob))*massb;
                double factor_visc = (2*etaa*etab/(etaa + etab))*(Va*Va + Vb*Vb)*dwdrab/r/massa;
                //if (t > 0.0002) factor_p = 0.0;

                p_force[0] += factor_p * DW.get(0);
                p_force[1] += factor_p * DW.get(1);
                eta_force[0] += factor_visc*v_rel[0];
                eta_force[1] += factor_visc*v_rel[1];

                //std::cout<<"W(0) = "<<Wab(0.0)<<" W(1.3/32) = "<<Wab(H)<<std::endl;
                //std::cout<<"Contribution by grad p: "<<factor_p*DW.get(2)<<", contribution by viscosity: "<<factor_visc*v_rel[2]<<", grad W magnitude: "<<DW.norm()<<", dwdr: "<<dwdrab<<std::endl;
            }

            ++Np;
        }

        vd.getProp<force>(a)[0] += p_force[0] + eta_force[0];
        vd.getProp<force>(a)[1] += p_force[1] + eta_force[1];

        //if (std::abs(vd.getProp<sdf>(a)) < (2.0*H))
        double colornorm;
        Point<2, double> colorgrad;
        if (colorfield)
        {
            colorgrad = vd.getProp<normal>(a);
            colornorm = sqrt(colorgrad[0]*colorgrad[0] + colorgrad[1]*colorgrad[1]);
        }
        else colornorm = (std::abs(vd.getProp<sdf>(a)) < (2.0*H))?1e10:0.0;
        if (colornorm > 0.01/H)
        {
            double stf0;
            double stf1;
            if (colorfield)
            {
                stf0 = -(alpha/rhoa)*vd.getProp<curvature>(a)*colorgrad[0];
                stf1 = -(alpha/rhoa)*vd.getProp<curvature>(a)*colorgrad[1];
            }
            else
            {
                double smoothing_factor = Wab(std::abs(vd.getProp<sdf>(a)))/(rhoa*Wab(0.0)); //latter is part of the equations, not the smoothing factor itself, but for convinience its in there
                smoothing_factor = Wab1D(std::abs(vd.getProp<sdf>(a)))/rhoa;
                double sign_corr = 1.0;
                if (return_sign(vd.getProp<sdf>(a) > 0)) sign_corr = -1.0;
                stf0 = smoothing_factor*alpha*vd.getProp<curvature>(a)*vd.getProp<normal>(a)[0]*sign_corr;
                stf1 = smoothing_factor*alpha*vd.getProp<curvature>(a)*vd.getProp<normal>(a)[1]*sign_corr;
            }

            vd.getProp<force>(a)[0] += stf0;
            vd.getProp<force>(a)[1] += stf1;

            avg_alpha += vd.getProp<curvature>(a);
            numsurfparticles++;
            //if ((a.getKey() == 6342)||(a.getKey() == 5285)) std::cout<<"Particle "<<a.getKey()<<" Position: "<<vd.getPos(a)[0]<<", "<<vd.getPos(a)[1]<<", "<<vd.getPos(a)[2]<<", surface tension force: "<<stf0<<", "<<stf1<<", "<<stf2<<std::endl;
        }
        if (va.norm() > maxvel) maxvel = va.norm();
        if (p_force.norm() > max_p) max_p = p_force.norm();
        if (eta_force.norm() > max_eta) max_eta = eta_force.norm();
        ++numparticles;
        ++part;
    }
    avg_alpha = avg_alpha/numsurfparticles;
    if ((corr) && (((int) std::round(t/dt)%5) == 0)) std::cout<<"Time step: "<<t<<", Max p force: "<<max_p<<", Max eta force: "<<max_eta<<", Max vel: "<<maxvel<<", Average curvature: "<<avg_alpha<<", number of particles: "<<numparticles<<std::endl;
    if (numparticles == 0) throw std::invalid_argument("no particles");
}

void pred_corr_int(particles & vd, double dt, int & corr)
{
    // particle iterator
    auto part = vd.getDomainIterator();

    while (part.isNext())
    {
        // ... a
        auto a = part.get();

        if (!corr)
        {
            // store the values at time n:
            vd.getProp<pos_prev>(a) = vd.getPos(a);
            vd.getProp<vel_prev>(a) = vd.getProp<vel>(a);
            // calculate intermediate values
            Point<2, double> dx = 0.5*dt*vd.getProp<vel>(a);
            vd.getPos(a)[0] += dx[0];
            vd.getPos(a)[1] += dx[1];

            vd.getProp<vel>(a) = vd.getProp<vel>(a) + 0.5*dt*vd.getProp<force>(a);

        }
        else
        {
            // correct the accelerations and velocities
            Point<2, double> x = vd.getProp<pos_prev>(a) + dt*(vd.getProp<vel_prev>(a) + 0.5*dt*vd.getProp<force>(a));
            vd.getPos(a)[0] = x[0];
            vd.getPos(a)[1] = x[1];
            vd.getProp<vel>(a) = vd.getProp<vel_prev>(a) + dt*vd.getProp<force>(a);

        }
        ++part;
    }
    corr = !corr;

}

int main(int argc, char* argv[])
{
    np_0 = 0;
    np_1 = 1;
	openfpm_init(&argc, &argv);

	Box<2, double> domain({-l/2.0, -l/2.0}, {l/2.0, l/2.0});
	size_t sz[2] = {(size_t)(l/dp),(size_t)(l/dp)};

	size_t bc[2] = {PERIODIC, PERIODIC};
	Ghost<2, double> g(3.0*band_width);

	particles vd(0, domain, bc, g, DEC_GRAN(512));
	//particles_surface vd_s(vd.getDecomposition(), 0);

	openfpm::vector<std::string> names({"sdf","normal","curvature","surf_flag","rho","p","f","vel","vel_prev"});
	vd.setPropNames(names);

	auto particle_it = vd.getGridIterator(sz);

	while (particle_it.isNext())
	{
		double x = -l/2.0 + dp/2.0 + particle_it.get().get(0)*dp;
		double y = -l/2.0 + dp/2.0 + particle_it.get().get(1)*dp;

		//x = x + randZeroToOne()*dp/3.0;
		//y = y + randZeroToOne()*dp/3.0;

		vd.add();
		vd.getLastPos()[0] = x;
		vd.getLastPos()[1] = y;

		vd.template getLastProp<surf_flag>() = 0;

		double dist = DistancePointEllipse(A, B, abs(x), abs(y));
        dist = -1.0*return_sign(1 - sqrt((x/A)*(x/A) + (y/B)*(y/B)))*dist;

		if (return_sign(dist) < 0) np_0++;
		else np_1++;

		vd.template getLastProp<sdf>() = dist;
		vd.template getLastProp<normal>()[0] = 0.0;
		vd.template getLastProp<normal>()[1] = 0.0;
		vd.template getLastProp<curvature>() = 0.0;
		vd.template getLastProp<vel>()[0] = 0.0;
		vd.template getLastProp<vel>()[1] = 0.0;
		vd.template getLastProp<vel_prev>()[0] = 0.0;
		vd.template getLastProp<vel_prev>()[1] = 0.0;
		vd.template getLastProp<force>()[0] = 0.0;
		vd.template getLastProp<force>()[1] = 0.0;
		vd.template getLastProp<rho>() = 0.0;
		vd.template getLastProp<p>() = 0.0;

		++particle_it;
	}

	m_0 = M_0/np_0;
    m_1 = M_1/np_1;
    std::cout<<"dt should be smaller than:\n"<<0.25*H/c<<"\t"<<0.25*H*H/eta_phase0<<"\t"<<0.25*sqrt(rho_0*H*H*H/(2*M_PI*alpha))<<std::endl;

    vd.map();
    vd.ghost_get<rho>();
    auto NN = vd.getCellList(2*H);
    density_summation(vd, NN);

    EqOfState(vd);

    vd.deleteGhost();
    vd.write("2Dellipse_surfacetension_init_before_redistancing");

    Redist_options rdistoptions;
    rdistoptions.H = dp;
    rdistoptions.r_cutoff_factor = 3.0;
    rdistoptions.sampling_radius = 0.75*band_width;
    rdistoptions.tolerance = 1e-14;
    rdistoptions.compute_normals = 1;
    rdistoptions.compute_curvatures = 1;
    rdistoptions.write_sdf = 0;

    particle_cp_redistancing <particles, taylor4>pcprdist_init(vd, rdistoptions);
    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

    if (colorfield)
    {
        calc_surface_normals(vd, NN);
        calc_curvature(vd, NN);
    }
    else
    {
        pcprdist_init.run_redistancing();
        //calc_surface_normals(vd, NN);
        //calc_curvature(vd, NN);
    }

    vd.write("2Dellipse_surfacetension_init");

    rdistoptions.write_sdf = 1;
    rdistoptions.tolerance = 1e-14;
    rdistoptions.r_cutoff_factor = 3.0;
    rdistoptions.compute_normals = 1;
    rdistoptions.compute_curvatures = 1;

    size_t write = 0;
    size_t write_interim = 0;
    size_t rdist = 0;
    size_t largebw = 0;
    size_t it = 0;
    size_t it_reb = 0;
    t = 0.0;
    corr = 0;

    //vd.load("pcp_sdf0_t_0_59_debug");

    while (t <= t_end)
    {
        Vcluster<> & v_cl = create_vcluster();
        timer it_time;

        vd.map();
        vd.ghost_get<rho>();
        vd.updateCellList(NN);


        // Calculate pressure from the density
        density_summation(vd, NN);
        EqOfState(vd);

        if ((rdist < t*10000000) && (t != 0.0))
        {
            if (colorfield)
            {
                calc_surface_normals(vd, NN);
                calc_curvature(vd, NN);
            }
            else
            {
                if (largebw < t*100)
                {
                    std::cout<<"Increased sampling radius for dt."<<std::endl;
                    rdistoptions.sampling_radius = 3.0*band_width;
                }
                //std::cout.setstate(std::ios_base::failbit);
                particle_cp_redistancing <particles, taylor4>pcprdist(vd, rdistoptions);
                pcprdist.run_redistancing();
                //std::cout.clear();
                if (largebw < t*100)
                {
                    rdistoptions.sampling_radius = 0.75*band_width;
                    largebw++;
                }
                //calc_surface_normals(vd, NN);
                //calc_curvature(vd, NN);
            }
            rdist++;
        }

        // Calc forces
        calc_forces(vd, NN);

        // Predictor step
        it++;
        pred_corr_int(vd, dt, corr);
        vd.map();

        density_summation(vd, NN);
        EqOfState(vd);
        calc_forces(vd, NN);

        // corrector step
        pred_corr_int(vd, dt, corr);

        t += dt;

        if (write_interim < t*10)
        {
            vd.save("2Dellipse_interim_save");
            write_interim++;
        }

        if (write < t*100)
        {
            vd.deleteGhost();

            vd.write_frame("Geometry_2Dellipse_surfacetension_v26/Geometry_2Dellipse_surfacetension_rdist_v26",write);
            write++;

            if (v_cl.getProcessUnitID() == 0)
            {std::cout << "TIME: " << t << std::endl;}
        }
/*        if (t >= 0.01)
        {
            std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
            std::cout << "Time required for simulation = " << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << "[ms]" << std::endl;
        }*/

    }
    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    std::cout << "Time required for simulation = " << std::chrono::duration_cast<std::chrono::minutes>(end - begin).count() << "[min]" << std::endl;
    openfpm_finalize();
}