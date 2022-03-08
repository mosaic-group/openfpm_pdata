#include "config.h"
#include <iostream>
#include "FiniteDifference/FD_Solver.hpp"
#include "FiniteDifference/FD_op.hpp"
#include "FiniteDifference/FD_expressions.hpp"
#include "FiniteDifference/util/EqnsStructFD.hpp"
int main(int argc, char* argv[])
{
    {   openfpm_init(&argc,&argv);
        using namespace FD;
        timer tt2;
        tt2.start();
        size_t gd_sz = 81;
        constexpr int x = 0;
        constexpr int y = 1;
        const size_t szu[2] = {gd_sz,gd_sz};
        int sz[2]={int(gd_sz),int(gd_sz)};
        Box<2, double> box({0, 0}, {1,1});
        periodicity<2> bc = {NON_PERIODIC, NON_PERIODIC};
        double spacing;
        spacing = 1.0 / (sz[0] - 1);
        Ghost<2,long int> ghost(1);
        auto &v_cl = create_vcluster();
        typedef aggregate<double, VectorS<2, double>, VectorS<2, double>,VectorS<2, double>,double,VectorS<2, double>,double,double> LidCavity;
        grid_dist_id<2, double, LidCavity> domain(szu, box,ghost,bc);
        double x0, y0, x1, y1;
        x0 = box.getLow(0);
        y0 = box.getLow(1);
        x1 = box.getHigh(0);
        y1 = box.getHigh(1);
        auto it = domain.getDomainIterator();
        while (it.isNext())
        {
            auto key = it.get();
            auto gkey = it.getGKey(key);
            double x = gkey.get(0) * domain.spacing(0);
            double y = gkey.get(1) * domain.spacing(1);
            if (y==1)
            {domain.get<1>(key) = 1.0;}
            else
            {domain.get<1>(key) = 0.0;}

            ++it;
        }
        domain.ghost_get<0>();
        Derivative_x Dx;
        Derivative_y Dy;
        Derivative_xx Dxx;
        Derivative_yy Dyy;
        auto P = getV<0>(domain);
        auto V = getV<1>(domain);
        auto RHS = getV<2>(domain);
        auto dV = getV<3>(domain);
        auto div = getV<4>(domain);
        auto V_star = getV<5>(domain);
        auto H = getV<6>(domain);
        auto dP = getV<7>(domain);        
        eq_id x_comp,y_comp;
        x_comp.setId(0);
        y_comp.setId(1);
        double sum, sum1, sum_k,V_err_old;
        auto StokesX=(V[x]*Dx(V_star[x])+V[y]*Dy(V_star[x]))-(1.0/100.0)*(Dxx(V_star[x])+Dyy(V_star[x]));
        auto StokesY=(V[x]*Dx(V_star[y])+V[y]*Dy(V_star[y]))-(1.0/100.0)*(Dxx(V_star[y])+Dyy(V_star[y]));
        RHS=0;
        dV=0;
        double V_err=1;
        int n;
        while (V_err >= 0.05 && n <= 50) {
            if (n%5==0){
                domain.ghost_get<0,1>(SKIP_LABELLING);
                domain.write_frame("LID",n,BINARY);
                domain.ghost_get<0>();
            }
            domain.ghost_get<0>(SKIP_LABELLING);
            RHS[x] = -Dx(P);
            RHS[y] = -Dy(P);
            FD_scheme<equations2d2,decltype(domain)> Solver(ghost,domain);
            Solver.impose(StokesX, {1,1},{sz[0]-2,sz[1]-2}, RHS[x], x_comp);
            Solver.impose(StokesY, {1,1},{sz[0]-2,sz[1]-2}, RHS[y], y_comp);
            Solver.impose(V_star[x], {0,sz[1]-1},{sz[0]-1,sz[1]-1}, 1.0, x_comp);
            Solver.impose(V_star[y], {0,sz[1]-1},{sz[0]-1,sz[1]-1}, 0.0, y_comp);
            Solver.impose(V_star[x], {0,1},{0,sz[1]-2}, 0.0, x_comp);
            Solver.impose(V_star[y], {0,1},{0,sz[1]-2}, 0.0, y_comp);
            Solver.impose(V_star[x], {sz[0]-1,1},{sz[0]-1,sz[1]-2}, 0.0, x_comp);
            Solver.impose(V_star[y], {sz[0]-1,1},{sz[0]-1,sz[1]-2}, 0.0, y_comp);
            Solver.impose(V_star[x], {0,0},{sz[0]-1,0}, 0.0, x_comp);
            Solver.impose(V_star[y], {0,0},{sz[0]-1,0}, 0.0, y_comp);
            Solver.solve(V_star[x], V_star[y]);

            domain.ghost_get<5>(SKIP_LABELLING);
            div = (Dx(V_star[x]) + Dy(V_star[y]));

            FD_scheme<equations2d1E,decltype(domain)> SolverH(ghost,domain);
            auto Helmholtz = Dxx(H)+Dyy(H);
            SolverH.impose(Helmholtz,{1,1},{sz[0]-2,sz[1]-2},div);
            SolverH.impose(Dy(H), {0,sz[0]-1},{sz[0]-1,sz[1]-1},0);
            SolverH.impose(Dx(H), {sz[0]-1,1},{sz[0]-1,sz[1]-2},0);
            SolverH.impose(Dx(H), {0,0},{sz[0]-1,0},0,x_comp,true);
            SolverH.impose(H, {0,0},{0,0},0);
            SolverH.impose(Dy(H), {0,1},{0,sz[1]-2},0);
            SolverH.solve(H);
            domain.ghost_get<1,4,6>(SKIP_LABELLING);
            P = P - 0.01*(div-0.5*(V[x]*Dx(H)+V[y]*Dy(H)));
            V_star[0] = V_star[0] - Dx(H);
            V_star[1] = V_star[1] - Dy(H);
            sum = 0;
            sum1 = 0;
            auto it2 = domain.getDomainIterator();
            while (it2.isNext()) {
                auto p = it2.get();
                sum += (domain.getProp<5>(p)[0] - domain.getProp<1>(p)[0]) *
                       (domain.getProp<5>(p)[0] - domain.getProp<1>(p)[0]) +
                       (domain.getProp<5>(p)[1] - domain.getProp<1>(p)[1]) *
                       (domain.getProp<5>(p)[1] - domain.getProp<1>(p)[1]);
                sum1 += domain.getProp<5>(p)[0] * domain.getProp<5>(p)[0] +
                        domain.getProp<5>(p)[1] * domain.getProp<5>(p)[1];
                ++it2;
            }
            V[x] = V_star[x];
            V[y] = V_star[y];
            v_cl.sum(sum);
            v_cl.sum(sum1);
            v_cl.execute();
            sum = sqrt(sum);
            sum1 = sqrt(sum1);
            V_err_old = V_err;
            V_err = sum / sum1;
            n++;
            if (v_cl.rank() == 0) {
                std::cout << "Rel l2 cgs err in V = " << V_err << "."<< std::endl;
            }
        }
        domain.write("LID_final");
        tt2.stop();
        if (v_cl.rank() == 0) {
            std::cout << "The entire solver  took " << tt2.getcputime() << "(CPU) ------ " << tt2.getwct()
                      << "(Wall) Seconds.";
        }
    }   
    openfpm_finalize();

}