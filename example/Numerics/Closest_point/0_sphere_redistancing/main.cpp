/* Benchmark: Redistance a sphere which is only known as color function (using OpenFPM + Algoim)
 * Date : 21 June 2021
 * Author : sachin
 */

#include <cmath>
#include <iostream>
#include "Grid/grid_dist_id.hpp"
#include "data_type/aggregate.hpp"
#include "VCluster/VCluster.hpp"
#include "Vector/Vector.hpp"
#include "FiniteDifference/util/common.hpp"
#include "FiniteDifference/eq.hpp"
#include "FiniteDifference/FD_op.hpp"
#include "level_set/closest_point/closest_point.hpp"

constexpr int SIM_DIM = 3;
constexpr int POLY_ORDER = 5;
constexpr int SIM_GRID_SIZE = 64;
constexpr double PI = 3.141592653589793;

// Fields - phi, cp
using GridDist = grid_dist_id<SIM_DIM,double,aggregate<double,double[3]>>;

// Grid size on each dimension
const long int sz[SIM_DIM] = {SIM_GRID_SIZE, SIM_GRID_SIZE, SIM_GRID_SIZE};
const size_t szu[SIM_DIM] = {(size_t) sz[0], (size_t) sz[1], (size_t) sz[2]};

// 2D physical domain
Box<SIM_DIM,double> domain({-1.5,-1.5,-1.5},{1.5,1.5,1.5});

constexpr int x = 0;
constexpr int y = 1;
constexpr int z = 2;

// Alias for properties on the grid
constexpr int phi = 0;
constexpr int cp = 1;

double nb_gamma = 0.0;


constexpr int narrow_band_half_width = 8;

typedef struct EllipseParameters{
    double origin[3];
    double radiusA;
    double radiusB;
    double radiusC;
    double eccentricity;
} EllipseParams;

// Generate the initial +1, -1 field on a sphere
template<const unsigned int phi_field, typename grid_type>
void initializeColourFunc(grid_type &gd, const EllipseParams &params)
{
    // Note: Since we use a Non-periodic boundary, ghost_get does not update ghost layers of sim box.
    // Therefore we are initializing the ghost layer with non-zero values.
    auto it = gd.getDomainGhostIterator();
    while(it.isNext())
    {
        auto key = it.get();
        Point<grid_type::dims, double> coords = gd.getPos(key);
        double posx = coords.get(x) + domain.getLow(x);
        double posy = coords.get(y) + domain.getLow(y);
        double posz = coords.get(z) + domain.getLow(z);
        
        double phi_val = 1.0 - sqrt(((posx - params.origin[0])/params.radiusA)*((posx - params.origin[0])/params.radiusA) + ((posy - params.origin[1])/params.radiusB)*((posy - params.origin[1])/params.radiusB) + ((posz - params.origin[2])/params.radiusC)*((posz - params.origin[2])/params.radiusC));
        gd.template get<phi_field>(key) = (phi_val<0)?-1.0:1.0;
        ++it;
    }
}

// Computes the error in redistancing for a unit sphere
void estimateErrorInReinit(GridDist &gd, EllipseParams &params)
{
    double max_phi_err = -1.0;


    auto it = gd.getDomainIterator();
    while(it.isNext())
    {
        auto key = it.get();
        Point<GridDist::dims, double> coords = gd.getPos(key);
        double posx = coords.get(0) + domain.getLow(x);
        double posy = coords.get(1) + domain.getLow(y);
        double posz = coords.get(2) + domain.getLow(z);

        double exact_phi_val = 1.0 - sqrt(((posx - params.origin[0])/params.radiusA)*((posx - params.origin[0])/params.radiusA) + ((posy - params.origin[1])/params.radiusB)*((posy - params.origin[1])/params.radiusB) + ((posz - params.origin[2])/params.radiusC)*((posz - params.origin[2])/params.radiusC));
        
        double phi_err = std::abs(exact_phi_val - gd.template get<phi>(key)); // / std::abs(exact_phi);

        if(phi_err > max_phi_err)
            max_phi_err = phi_err;
        ++it;
    }
    
    std::cout<<"Error in redistancing = "<<max_phi_err<<"\n";

    return;
}

int main(int argc, char* argv[])
{
	// Initialize the library
	openfpm_init(&argc, &argv);
    // Create VCluster
    Vcluster<> &v_cl = create_vcluster();

	periodicity<SIM_DIM> grid_bc = {NON_PERIODIC, NON_PERIODIC, NON_PERIODIC};
	// Ghost in grid units
	Ghost <SIM_DIM, long int> grid_ghost(2*narrow_band_half_width);
	GridDist gdist(szu, domain, grid_ghost, grid_bc);

    // Set property names
    openfpm::vector < std::string > prop_names;
    prop_names.add("ls_phi");
    prop_names.add("closest_point");
    gdist.setPropNames(prop_names);

    EllipseParams params;
    params.origin[x] = 0.0;
    params.origin[y] = 0.0;
    params.origin[z] = 0.0;
    params.radiusA = 1.0;
    params.radiusB = 1.0;
    params.radiusC = 1.0;

    initializeColourFunc<phi>(gdist, params);
    nb_gamma = narrow_band_half_width * gdist.spacing(0);

    gdist.template ghost_get<phi>();

    auto &patches = gdist.getLocalGridsInfo();

    for(int i = 0; i < patches.size();i++)
    {
        double max_x = patches.get(i).Dbox.getHigh(x) + patches.get(i).origin[x];
        double min_x = patches.get(i).Dbox.getLow(x) + patches.get(i).origin[x];
        double max_y = patches.get(i).Dbox.getHigh(y) + patches.get(i).origin[y];
        double min_y = patches.get(i).Dbox.getLow(y) + patches.get(i).origin[y];
        double max_z = patches.get(i).Dbox.getHigh(z) + patches.get(i).origin[z];
        double min_z = patches.get(i).Dbox.getLow(z) + patches.get(i).origin[z];
        std::cout<<"In processor "<<v_cl.getProcessUnitID()<<", patch "<<i<<", min = ("<<min_x<<", "<<min_y<<", "<<min_z<<"), max = ("<<max_x<<", "<<max_y<<", "<<max_z<<")"<<std::endl; 
    }
 
    std::cout<<"Grid spacing = "<<gdist.spacing(x)<<"\n";

    // gdist.write_frame("output", 0);

    estimateClosestPoint<phi, cp, POLY_ORDER>(gdist, 3.0);
    gdist.template ghost_get<cp>();

    // Redistance Levelset
    reinitializeLS<phi, cp, POLY_ORDER>(gdist, 3.0);
    gdist.template ghost_get<phi>();

    estimateErrorInReinit(gdist, params);

    // gdist.write_frame("output", 1);

	openfpm_finalize();
}


