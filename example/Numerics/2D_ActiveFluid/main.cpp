//
// Created by Abhinav Singh on 15.03.20.
//
//#define SE_CLASS1
#include "config.h"
#define BOOST_MPL_CFG_NO_PREPROCESSED_HEADERS
#define BOOST_MPL_LIMIT_VECTOR_SIZE 40
#include <iostream>
#include "DCPSE/DCPSE_op/DCPSE_op.hpp"
#include "DCPSE/DCPSE_op/DCPSE_Solver.hpp"
#include "Operators/Vector/vector_dist_operators.hpp"
#include "Vector/vector_dist_subset.hpp"
#include "DCPSE/DCPSE_op/EqnsStruct.hpp"
#include "OdeIntegrators/OdeIntegrators.hpp"

constexpr int x = 0;
constexpr int y = 1;
constexpr int POLARIZATION= 0,VELOCITY = 1, VORTICITY = 2, EXTFORCE = 3,PRESSURE = 4, STRAIN_RATE = 5, STRESS = 6, MOLFIELD = 7, DPOL = 8, DV = 9, VRHS = 10, F1 = 11, F2 = 12, F3 = 13, F4 = 14, F5 = 15, F6 = 16, V_T = 17, DIV = 18, DELMU = 19, HPB = 20, FE = 21, R = 22;

double eta = 1.0;
double nu = -0.5;
double gama = 0.1;
double zeta = 0.07;
double Ks = 1.0;
double Kb = 1.0;
double lambda = 0.1;

int wr_f;
int wr_at;
double V_err_eps;

void *vectorGlobal=nullptr,*vectorGlobal_bulk=nullptr,*vectorGlobal_boundary=nullptr;
const openfpm::vector<std::string>
PropNAMES={"00-Polarization","01-Velocity","02-Vorticity","03-ExternalForce","04-Pressure","05-StrainRate","06-Stress","07-MolecularField","08-DPOL","09-DV","10-VRHS","11-f1","12-f2","13-f3","14-f4","15-f5","16-f6","17-V_T","18-DIV","19-DELMU","20-HPB","21-FrankEnergy","22-R"};
typedef aggregate<VectorS<2, double>,VectorS<2, double>,double[2][2],VectorS<2, double>,double,double[2][2],double[2][2],VectorS<2, double>,VectorS<2, double>,VectorS<2, double>,VectorS<2, double>,double,double,double,double,double,double,VectorS<2, double>,double,double,double,double,double> Activegels;
typedef vector_dist_ws<2, double,Activegels> vector_type;
typedef vector_dist_subset<2, double, Activegels> vector_type2;

//Functor to Compute RHS of the time derivative of the polarity
template<typename DX,typename DY,typename DXX,typename DXY,typename DYY>
struct PolarEv
{
    DX &Dx;
    DY &Dy;
    DXX &Dxx;
    DXY &Dxy;
    DYY &Dyy;
    //Constructor
    PolarEv(DX &Dx,DY &Dy,DXX &Dxx,DXY &Dxy,DYY &Dyy):Dx(Dx),Dy(Dy),Dxx(Dxx),Dxy(Dxy),Dyy(Dyy)
    {}

    void operator()( const state_type_2d_ofp &X , state_type_2d_ofp &dxdt , const double t ) const
    {

        vector_type &Particles= *(vector_type *) vectorGlobal;
        vector_type2 &Particles_bulk= *(vector_type2 *) vectorGlobal_bulk;

        auto Pol=getV<POLARIZATION>(Particles);
        auto Pol_bulk=getV<POLARIZATION>(Particles_bulk);
        auto h = getV<MOLFIELD>(Particles);
        auto u = getV<STRAIN_RATE>(Particles);
        auto dPol = getV<DPOL>(Particles);
        auto W = getV<VORTICITY>(Particles);
        auto delmu = getV<DELMU>(Particles);
        auto H_p_b = getV<HPB>(Particles);
        auto r = getV<R>(Particles);
        auto dPol_bulk = getV<DPOL>(Particles_bulk);
        Pol[x]=X.data.get<0>();
        Pol[y]=X.data.get<1>();
        Particles.ghost_get<POLARIZATION>(SKIP_LABELLING);
        H_p_b = Pol[x] * Pol[x] + Pol[y] * Pol[y];

        h[y] = (Pol[x] * (Ks * Dyy(Pol[y]) + Kb * Dxx(Pol[y]) + (Ks - Kb) * Dxy(Pol[x])) -
                    Pol[y] * (Ks * Dxx(Pol[x]) + Kb * Dyy(Pol[x]) + (Ks - Kb) * Dxy(Pol[y])));

        h[x] = -gama * (lambda * delmu - nu * (u[x][x] * Pol[x] * Pol[x] + u[y][y] * Pol[y] * Pol[y] + 2 * u[x][y] * Pol[x] * Pol[y]) / (H_p_b));
 
        dPol_bulk[x] = ((h[x] * Pol[x] - h[y] * Pol[y]) / gama + lambda * delmu * Pol[x] -
                     nu * (u[x][x] * Pol[x] + u[x][y] * Pol[y]) + W[x][x] * Pol[x] +
                     W[x][y] * Pol[y]);
        dPol_bulk[y] = ((h[x] * Pol[y] + h[y] * Pol[x]) / gama + lambda * delmu * Pol[y] -
                     nu * (u[y][x] * Pol[x] + u[y][y] * Pol[y]) + W[y][x] * Pol[x] +
                     W[y][y] * Pol[y]);
        dxdt.data.get<0>()=dPol[x];
        dxdt.data.get<1>()=dPol[y];
    }
};

// Functor to calculate velocity and move particles with explicit euler
template<typename DX,typename DY,typename DXX,typename DXY,typename DYY>
struct CalcVelocity
{

    DX &Dx, &Bulk_Dx;
    DY &Dy, &Bulk_Dy;
    DXX &Dxx;
    DXY &Dxy;
    DYY &Dyy;

    double t_old;
    int ctr;

    //Constructor
    CalcVelocity(DX &Dx,DY &Dy,DXX &Dxx,DXY &Dxy,DYY &Dyy,DX &Bulk_Dx,DY &Bulk_Dy):Dx(Dx),Dy(Dy),Dxx(Dxx),Dxy(Dxy),Dyy(Dyy),Bulk_Dx(Bulk_Dx),Bulk_Dy(Bulk_Dy)
    {
        t_old = 0.0;
        ctr = 0;
    }

    void operator() (state_type_2d_ofp &state, double t)
    {

        double dt = t - t_old;
        vector_type &Particles= *(vector_type *) vectorGlobal;
        vector_type2 &Particles_bulk= *(vector_type2 *) vectorGlobal_bulk;
        vector_type2 &Particles_boundary= *(vector_type2 *) vectorGlobal_boundary;
        auto &v_cl = create_vcluster();

        timer tt;
        
        auto Pos = getV<PROP_POS>(Particles);
        auto Pol=getV<POLARIZATION>(Particles);
        auto V = getV<VELOCITY>(Particles);
        auto H_p_b = getV<HPB>(Particles);


        // skip in first time step
        if (dt != 0) {
            tt.start();
            //normalize polarization
            H_p_b = sqrt(Pol[x] * Pol[x] + Pol[y] * Pol[y]);
            Pol = Pol / H_p_b;
            Pos = Pos + (t-t_old)*V;
            Particles.map();
            Particles.ghost_get<POLARIZATION, EXTFORCE, DELMU>();
            Particles_bulk.update();
            Particles_boundary.update();
            tt.start();
            Dx.update(Particles);
            Dy.update(Particles);
            Dxy.update(Particles);
            Dxx.update(Particles);
            Dyy.update(Particles);

            Bulk_Dx.update(Particles_bulk);
            Bulk_Dy.update(Particles_bulk);

            state.data.get<0>()=Pol[x];
            state.data.get<1>()=Pol[y];

            tt.stop();
            if (v_cl.rank() == 0) {
                std::cout << "Updating operators took " << tt.getwct() << " seconds." << std::endl;
                std::cout << "Time step " << ctr << " : " << t << " over." << std::endl;
                std::cout << "----------------------------------------------------------" << std::endl;
            }
            ctr++;

        }
        auto Dyx = Dxy;
        t_old = t;
        tt.start();

        auto & bulk = Particles_bulk.getIds();
        auto & boundary = Particles_boundary.getIds();
        auto Pol_bulk=getV<POLARIZATION>(Particles_bulk);
        auto sigma = getV<STRESS>(Particles);
        auto r = getV<R>(Particles);
        auto h = getV<MOLFIELD>(Particles);
        auto FranckEnergyDensity = getV<FE>(Particles);
        auto f1 = getV<F1>(Particles);
        auto f2 = getV<F2>(Particles);
        auto f3 = getV<F3>(Particles);
        auto f4 = getV<F4>(Particles);
        auto f5 = getV<F5>(Particles);
        auto f6 = getV<F6>(Particles);
        auto dV = getV<DV>(Particles);
        auto delmu = getV<DELMU>(Particles); // why is delmu a property, not a constant?
        auto g = getV<EXTFORCE>(Particles);
        auto P = getV<PRESSURE>(Particles);
        auto P_bulk = getV<PRESSURE>(Particles_bulk); //Pressure only on inside
        auto RHS = getV<VRHS>(Particles);
        auto RHS_bulk = getV<VRHS>(Particles_bulk);
        auto div = getV<DIV>(Particles);
        auto V_t = getV<V_T>(Particles);
        auto u = getV<STRAIN_RATE>(Particles);
        auto W = getV<VORTICITY>(Particles);

        Pol_bulk[x]=state.data.get<0>();
        Pol_bulk[y]=state.data.get<1>();
        Particles.ghost_get<POLARIZATION>(SKIP_LABELLING);

        eq_id x_comp, y_comp;
        x_comp.setId(0);
        y_comp.setId(1);

        int n = 0,nmax = 300,errctr = 0, Vreset = 0;
        double V_err = 1,V_err_eps = 5 * 1e-3, V_err_old,sum, sum1;
        std::cout << "Calculate velocity (step t=" << t << ")" << std::endl;
        tt.start();
        petsc_solver<double> solverPetsc;
        solverPetsc.setSolver(KSPGMRES);
        solverPetsc.setPreconditioner(PCJACOBI);
        Particles.ghost_get<POLARIZATION>(SKIP_LABELLING);
        // calculate stress
        sigma[x][x] =
                -Ks * Dx(Pol[x]) * Dx(Pol[x]) - Kb * Dx(Pol[y]) * Dx(Pol[y]) + (Kb - Ks) * Dy(Pol[x]) * Dx(Pol[y]);
        sigma[x][y] =
                -Ks * Dy(Pol[y]) * Dx(Pol[y]) - Kb * Dy(Pol[x]) * Dx(Pol[x]) + (Kb - Ks) * Dx(Pol[y]) * Dx(Pol[x]);
        sigma[y][x] =
                -Ks * Dx(Pol[x]) * Dy(Pol[x]) - Kb * Dx(Pol[y]) * Dy(Pol[y]) + (Kb - Ks) * Dy(Pol[x]) * Dy(Pol[y]);
        sigma[y][y] =
                -Ks * Dy(Pol[y]) * Dy(Pol[y]) - Kb * Dy(Pol[x]) * Dy(Pol[x]) + (Kb - Ks) * Dx(Pol[y]) * Dy(Pol[x]);
        Particles.ghost_get<STRESS>(SKIP_LABELLING);

        // if R == 0 then set to 1 to avoid division by zero for defects
        r = Pol[x] * Pol[x] + Pol[y] * Pol[y];
        for (int j = 0; j < bulk.size(); j++) {
            auto p = bulk.get<0>(j);
            Particles.getProp<R>(p) = (Particles.getProp<R>(p) == 0) ? 1 : Particles.getProp<R>(p);
        }
        for (int j = 0; j < boundary.size(); j++) {
            auto p = boundary.get<0>(j);
            Particles.getProp<R>(p) = (Particles.getProp<R>(p) == 0) ? 1 : Particles.getProp<R>(p);
        }

        // calculate traversal molecular field (H_perpendicular)
        h[y] = (Pol[x] * (Ks * Dyy(Pol[y]) + Kb * Dxx(Pol[y]) + (Ks - Kb) * Dxy(Pol[x])) -
                Pol[y] * (Ks * Dxx(Pol[x]) + Kb * Dyy(Pol[x]) + (Ks - Kb) * Dxy(Pol[y])));
        Particles.ghost_get<MOLFIELD>(SKIP_LABELLING);

        // calulate FranckEnergyDensity
        FranckEnergyDensity = (Ks / 2.0) *
                              ((Dx(Pol[x]) * Dx(Pol[x])) + (Dy(Pol[x]) * Dy(Pol[x])) +
                               (Dx(Pol[y]) * Dx(Pol[y])) +
                               (Dy(Pol[y]) * Dy(Pol[y]))) +
                              ((Kb - Ks) / 2.0) * ((Dx(Pol[y]) - Dy(Pol[x])) * (Dx(Pol[y]) - Dy(Pol[x])));
        Particles.ghost_get<FE>(SKIP_LABELLING);

        // calculate preactors for LHS of Stokes Equation.
        f1 = gama * nu * Pol[x] * Pol[x] * (Pol[x] * Pol[x] - Pol[y] * Pol[y]) / (r);
        f2 = 2.0 * gama * nu * Pol[x] * Pol[y] * (Pol[x] * Pol[x] - Pol[y] * Pol[y]) / (r);
        f3 = gama * nu * Pol[y] * Pol[y] * (Pol[x] * Pol[x] - Pol[y] * Pol[y]) / (r);
        f4 = 2.0 * gama * nu * Pol[x] * Pol[x] * Pol[x] * Pol[y] / (r);
        f5 = 4.0 * gama * nu * Pol[x] * Pol[x] * Pol[y] * Pol[y] / (r);
        f6 = 2.0 * gama * nu * Pol[x] * Pol[y] * Pol[y] * Pol[y] / (r);
        Particles.ghost_get<F1, F2, F3, F4, F5, F6>(SKIP_LABELLING);
        texp_v<double> Dxf1 = Dx(f1),Dxf2 = Dx(f2),Dxf3 = Dx(f3),Dxf4 = Dx(f4),Dxf5 = Dx(f5),Dxf6 = Dx(f6),
                        Dyf1 = Dy(f1),Dyf2 = Dy(f2),Dyf3 = Dy(f3),Dyf4 = Dy(f4),Dyf5 = Dy(f5),Dyf6 = Dy(f6);

        // calculate RHS of Stokes Equation without pressure
        dV[x] = -0.5 * Dy(h[y]) + zeta * Dx(delmu * Pol[x] * Pol[x]) + zeta * Dy(delmu * Pol[x] * Pol[y]) -
                zeta * Dx(0.5 * delmu * (Pol[x] * Pol[x] + Pol[y] * Pol[y])) -
                0.5 * nu * Dx(-2.0 * h[y] * Pol[x] * Pol[y])
                - 0.5 * nu * Dy(h[y] * (Pol[x] * Pol[x] - Pol[y] * Pol[y])) - Dx(sigma[x][x]) -
                Dy(sigma[x][y]) -
                g[x]
                - 0.5 * nu * Dx(-gama * lambda * delmu * (Pol[x] * Pol[x] - Pol[y] * Pol[y]))
                - 0.5 * Dy(-2.0 * gama * lambda * delmu * (Pol[x] * Pol[y]));

        dV[y] = -0.5 * Dx(-h[y]) + zeta * Dy(delmu * Pol[y] * Pol[y]) + zeta * Dx(delmu * Pol[x] * Pol[y]) -
                zeta * Dy(0.5 * delmu * (Pol[x] * Pol[x] + Pol[y] * Pol[y])) -
                0.5 * nu * Dy(2.0 * h[y] * Pol[x] * Pol[y])
                - 0.5 * nu * Dx(h[y] * (Pol[x] * Pol[x] - Pol[y] * Pol[y])) - Dx(sigma[y][x]) -
                Dy(sigma[y][y]) -
                g[y]
                - 0.5 * nu * Dy(gama * lambda * delmu * (Pol[x] * Pol[x] - Pol[y] * Pol[y]))
                - 0.5 * Dx(-2.0 * gama * lambda * delmu * (Pol[x] * Pol[y]));
        Particles.ghost_get<DV>(SKIP_LABELLING);

        // Encode LHS of the Stokes Equations
        auto Stokes1 = eta * (Dxx(V[x]) + Dyy(V[x]))
                       + 0.5 * nu * (Dxf1 * Dx(V[x]) + f1 * Dxx(V[x]))
                       + 0.5 * nu * (Dxf2 * 0.5 * (Dx(V[y]) + Dy(V[x])) + f2 * 0.5 * (Dxx(V[y]) + Dyx(V[x])))
                       + 0.5 * nu * (Dxf3 * Dy(V[y]) + f3 * Dyx(V[y]))
                       + 0.5 * nu * (Dyf4 * Dx(V[x]) + f4 * Dxy(V[x]))
                       + 0.5 * nu * (Dyf5 * 0.5 * (Dx(V[y]) + Dy(V[x])) + f5 * 0.5 * (Dxy(V[y]) + Dyy(V[x])))
                       + 0.5 * nu * (Dyf6 * Dy(V[y]) + f6 * Dyy(V[y]));
        auto Stokes2 = eta * (Dxx(V[y]) + Dyy(V[y]))
                       - 0.5 * nu * (Dyf1 * Dx(V[x]) + f1 * Dxy(V[x]))
                       - 0.5 * nu * (Dyf2 * 0.5 * (Dx(V[y]) + Dy(V[x])) + f2 * 0.5 * (Dxy(V[y]) + Dyy(V[x])))
                       - 0.5 * nu * (Dyf3 * Dy(V[y]) + f3 * Dyy(V[y]))
                       + 0.5 * nu * (Dxf4 * Dx(V[x]) + f4 * Dxx(V[x]))
                       + 0.5 * nu * (Dxf5 * 0.5 * (Dx(V[y]) + Dy(V[x])) + f5 * 0.5 * (Dxx(V[y]) + Dyx(V[x])))
                       + 0.5 * nu * (Dxf6 * Dy(V[y]) + f6 * Dyx(V[y]));

        tt.stop();

        std::cout << "Init of Velocity took " << tt.getwct() << " seconds." << std::endl;

        tt.start();
        V_err = 1;
        n = 0;
        errctr = 0;
        if (Vreset == 1) {
            P = 0;
            Vreset = 0;
        }
        P=0;

        // integrate velocity
        Particles.ghost_get<PRESSURE>(SKIP_LABELLING);
        RHS_bulk[x] = dV[x] + Bulk_Dx(P);
        RHS_bulk[y] = dV[y] + Bulk_Dy(P);
        Particles.ghost_get<VRHS>(SKIP_LABELLING);

        // prepare solver
        DCPSE_scheme<equations2d2, vector_type> Solver(Particles);
        Solver.impose(Stokes1, bulk, RHS[0], x_comp);
        Solver.impose(Stokes2, bulk, RHS[1], y_comp);
        Solver.impose(V[x], boundary, 0, x_comp);
        Solver.impose(V[y], boundary, 0, y_comp);
        Solver.solve_with_solver(solverPetsc, V[x], V[y]);
        Particles.ghost_get<VELOCITY>(SKIP_LABELLING);
        div = -(Dx(V[x]) + Dy(V[y]));
        P_bulk = P + div;

        // approximate velocity
        while (V_err >= V_err_eps && n <= nmax) {
            Particles.ghost_get<PRESSURE>(SKIP_LABELLING);
            RHS_bulk[x] = dV[x] + Bulk_Dx(P);
            RHS_bulk[y] = dV[y] + Bulk_Dy(P);
            Particles.ghost_get<VRHS>(SKIP_LABELLING);
            Solver.reset_b();
            Solver.impose_b(bulk, RHS[0], x_comp);
            Solver.impose_b(bulk, RHS[1], y_comp);
            Solver.impose_b(boundary, 0, x_comp);
            Solver.impose_b(boundary, 0, y_comp);
            Solver.solve_with_solver(solverPetsc, V[x], V[y]);
            Particles.ghost_get<VELOCITY>(SKIP_LABELLING);
            div = -(Dx(V[x]) + Dy(V[y]));
            P_bulk = P + div;
            // calculate error
            sum = 0;
            sum1 = 0;
            for (int j = 0; j < bulk.size(); j++) {
                auto p = bulk.get<0>(j);
                sum += (Particles.getProp<V_T>(p)[0] - Particles.getProp<VELOCITY>(p)[0]) *
                       (Particles.getProp<V_T>(p)[0] - Particles.getProp<VELOCITY>(p)[0]) +
                       (Particles.getProp<V_T>(p)[1] - Particles.getProp<VELOCITY>(p)[1]) *
                       (Particles.getProp<V_T>(p)[1] - Particles.getProp<VELOCITY>(p)[1]);
                sum1 += Particles.getProp<VELOCITY>(p)[0] * Particles.getProp<VELOCITY>(p)[0] +
                        Particles.getProp<VELOCITY>(p)[1] * Particles.getProp<VELOCITY>(p)[1];
            }
            V_t = V;
            v_cl.sum(sum);
            v_cl.sum(sum1);
            v_cl.execute();
            sum = sqrt(sum);
            sum1 = sqrt(sum1);
            V_err_old = V_err;
            V_err = sum / sum1;
            if (V_err > V_err_old || abs(V_err_old - V_err) < 1e-8) {
                errctr++;
                //alpha_P -= 0.1;
            } else {
                errctr = 0;
            }
            if (n > 3) {
                if (errctr > 3) {
                    std::cout << "CONVERGENCE LOOP BROKEN DUE TO INCREASE/VERY SLOW DECREASE IN DIVERGENCE" << std::endl;
                    Vreset = 1;
                    break;
                } else {
                    Vreset = 0;
                }
            }
            n++;

        }
        tt.stop();

        Particles.ghost_get<VELOCITY>(SKIP_LABELLING);
        // calculate strain rate
        u[x][x] = Dx(V[x]);
        u[x][y] = 0.5 * (Dx(V[y]) + Dy(V[x]));
        u[y][x] = 0.5 * (Dy(V[x]) + Dx(V[y]));
        u[y][y] = Dy(V[y]);

        if (v_cl.rank() == 0) {
            std::cout << "Rel l2 cgs err in V = " << V_err << " and took " << tt.getwct() << " seconds with " << n
                      << " iterations. dt for the stepper is " << dt
                      << std::endl;
        }

        // calculate vorticity
        W[x][x] = 0;
        W[x][y] = 0.5 * (Dy(V[x]) - Dx(V[y]));
        W[y][x] = 0.5 * (Dx(V[y]) - Dy(V[x]));
        W[y][y] = 0;

        if (ctr%wr_at==0 || ctr==wr_f){
            Particles.deleteGhost();
            Particles.write_frame("Polar",ctr);
            Particles.ghost_get<POLARIZATION>();
        }
    }
};

int main(int argc, char* argv[])
{
    {   openfpm_init(&argc,&argv);
        timer tt2;
        tt2.start();
        size_t Gd = int(std::atof(argv[1]));
        double tf = std::atof(argv[2]);
        double dt = tf/std::atof(argv[3]);
        wr_f=int(std::atof(argv[3]));
        wr_at=1;
        V_err_eps = 5e-4;

        double boxsize = 10;
        const size_t sz[2] = {Gd, Gd};
        Box<2, double> box({0, 0}, {boxsize, boxsize});
        double Lx = box.getHigh(0),Ly = box.getHigh(1);
        size_t bc[2] = {NON_PERIODIC, NON_PERIODIC};
        double spacing = box.getHigh(0) / (sz[0] - 1),rCut = 3.9 * spacing;
        int ord = 2;
        Ghost<2, double> ghost(rCut);
        auto &v_cl = create_vcluster();
        vector_dist_ws<2, double,Activegels> Particles(0, box, bc, ghost);
        Particles.setPropNames(PropNAMES);

        double x0=box.getLow(0), y0=box.getLow(1), x1=box.getHigh(0), y1=box.getHigh(1);
        auto it = Particles.getGridIterator(sz);
        while (it.isNext()) {
            Particles.add();
            auto key = it.get();
            double xp = key.get(0) * it.getSpacing(0),yp = key.get(1) * it.getSpacing(1);
            Particles.getLastPos()[x] = xp;
            Particles.getLastPos()[y] = yp;
            if (xp != x0 && yp != y0 && xp != x1 && yp != y1)
                Particles.getLastSubset(0);
            else
                Particles.getLastSubset(1);
            ++it;
        }
        Particles.map();
        Particles.ghost_get<POLARIZATION>();

        //auto Pos = getV<PROP_POS>(Particles);

        auto Pol = getV<POLARIZATION>(Particles);
        auto V = getV<VELOCITY>(Particles);
        auto g = getV<EXTFORCE>(Particles);
        auto P = getV<PRESSURE>(Particles);
        auto delmu = getV<DELMU>(Particles);
        auto dPol = getV<DPOL>(Particles);

        g = 0;delmu = -1.0;P = 0;V = 0;
        auto it2 = Particles.getDomainIterator();
        while (it2.isNext()) {
            auto p = it2.get();
            Point<2, double> xp = Particles.getPos(p);
            Particles.getProp<POLARIZATION>(p)[x] = sin(2 * M_PI * (cos((2 * xp[x] - Lx) / Lx) - sin((2 * xp[y] - Ly) / Ly)));
            Particles.getProp<POLARIZATION>(p)[y] = cos(2 * M_PI * (cos((2 * xp[x] - Lx) / Lx) - sin((2 * xp[y] - Ly) / Ly)));
            ++it2;
        }
        Particles.ghost_get<POLARIZATION,EXTFORCE,DELMU>(SKIP_LABELLING);

        vector_dist_subset<2, double, Activegels> Particles_bulk(Particles,0);
        vector_dist_subset<2, double, Activegels> Particles_boundary(Particles,1);
        auto & bulk = Particles_bulk.getIds();
        auto & boundary = Particles_boundary.getIds();

        auto P_bulk = getV<PRESSURE>(Particles_bulk);//Pressure only on inside
        auto Pol_bulk = getV<POLARIZATION>(Particles_bulk);;
        auto dPol_bulk = getV<DPOL>(Particles_bulk);
        auto dV_bulk = getV<DV>(Particles_bulk);
        auto RHS_bulk = getV<VRHS>(Particles_bulk);
        auto div_bulk = getV<DIV>(Particles_bulk);

        Derivative_x Dx(Particles,ord,rCut), Bulk_Dx(Particles_bulk,ord,rCut);
        Derivative_y Dy(Particles, ord, rCut), Bulk_Dy(Particles_bulk, ord,rCut);
        Derivative_xy Dxy(Particles, ord, rCut);
        auto Dyx = Dxy;
        Derivative_xx Dxx(Particles, ord, rCut);
        Derivative_yy Dyy(Particles, ord, rCut);

        boost::numeric::odeint::runge_kutta4< state_type_2d_ofp,double,state_type_2d_ofp,double,boost::numeric::odeint::vector_space_algebra_ofp> rk4;

        vectorGlobal=(void *) &Particles;
        vectorGlobal_bulk=(void *) &Particles_bulk;
        vectorGlobal_boundary=(void *) &Particles_boundary;

        PolarEv<Derivative_x,Derivative_y,Derivative_xx,Derivative_xy,Derivative_yy> System(Dx,Dy,Dxx,Dxy,Dyy);
        CalcVelocity<Derivative_x,Derivative_y,Derivative_xx,Derivative_xy,Derivative_yy> CalcVelocityObserver(Dx,Dy,Dxx,Dxy,Dyy,Bulk_Dx,Bulk_Dy);

        state_type_2d_ofp tPol;
        tPol.data.get<0>()=Pol[x];
        tPol.data.get<1>()=Pol[y];

        eq_id vx, vy;
        vx.setId(0);
        vy.setId(1);
        timer tt;
        timer tt3;
        dPol = Pol;
        double V_err = 1, V_err_old;
        double tim=0;

        // intermediate time steps
        std::vector<double> inter_times;

        size_t steps = integrate_const(rk4 , System , tPol , tim , tf , dt, CalcVelocityObserver);


        std::cout << "Time steps: " << steps << std::endl;

        Pol_bulk[x]=tPol.data.get<0>();
        Pol_bulk[y]=tPol.data.get<1>();

        Particles.deleteGhost();
        Particles.write("Polar_Last");
        Dx.deallocate(Particles);
        Dy.deallocate(Particles);
        Dxy.deallocate(Particles);
        Dxx.deallocate(Particles);
        Dyy.deallocate(Particles);
        Bulk_Dx.deallocate(Particles_bulk);
        Bulk_Dy.deallocate(Particles_bulk);
        std::cout.precision(17);
        tt2.stop();
        if (v_cl.rank() == 0) {
            std::cout << "The simulation took " << tt2.getcputime() << "(CPU) ------ " << tt2.getwct()
                      << "(Wall) Seconds.";
        }
    }
    openfpm_finalize();
}
