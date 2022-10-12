#include "Grid/grid_dist_id.hpp"
#include "OpenFPM_version.hpp"
#include <chrono>
#include <fstream>


double dt = 0.05;
double t_final = 400;
double domainSize = 40.0;
const int meshSize = 64;
double meshSpacing = domainSize / meshSize;
double epsilon = meshSpacing;
double r_cut = 3 * epsilon;
double D = 0.01;
double kernel = dt * D * 15.0 * pow(meshSpacing/epsilon, 2)  / pow(epsilon * M_PI, 2);


void createStencil(std::vector<std::array<int, 2>>& stencil) {
    
    // calculate farthest neighbor distance (in mesh nodes)
    int neighborDistance = round(r_cut / meshSpacing);

    // width of nD cube
    int diameter = (2 * neighborDistance) + 1;

    // offset to shift nD cube to the center node
    int offset = floor(diameter / 2);

    // create nD cube
    for (int i = 0; i < pow(diameter, 2); ++i) {
        std::array<int, 2> neighbor;

        // compute stencil node
        for (int component = 0; component < 2; ++component) {
            neighbor[component] = ((int)round(i / pow(diameter, component)) % diameter) - offset;
        }

        // add to stencil
        stencil.push_back(neighbor);
    }

    // remove all neighbors outside of cutoff radius
    double cutoffRadius2 = r_cut * r_cut;
    double spacing = meshSpacing;
    auto it_rcut = std::remove_if(stencil.begin(), stencil.end(), [&cutoffRadius2, &spacing](std::array<int, 2> node)
    {
        // calculate squared distance
        // from center to neighbor node
        float distance2 = 0;
        for (int i = 0; i < 2; ++i) {
            distance2 += node[i] * node[i] * spacing * spacing;
        }
        return (distance2 > cutoffRadius2);
    });
    stencil.erase(it_rcut, stencil.end());

    // remove center particle
    auto it_center = std::remove_if(stencil.begin(), stencil.end(), [](std::array<int, 2> node)
    { return std::all_of(node.begin(), node.end(), [](int comp){return comp == 0;}); });
    stencil.erase(it_center, stencil.end());

    
    // make symmetric
    
    // remove half of the particles to create symmetric stencil
    auto it_symm = std::remove_if(stencil.begin(), stencil.end(), [](std::array<int, 2> node)
    {
        for (int i = 0; i < 2; ++i) {
            // remove node if component > 0
            if (node[i] > 0)
                return true;
            // keep node if component < 0
            if (node[i] < 0)
                return false;
            // if component == 0, check next dimension
        }

        // remove center node [0, 0, 0]
        return true;
    });
    stencil.erase(it_symm, stencil.end());

}

int main(int argc, char* argv[])
{

	// Initialize the library
	openfpm_init(&argc,&argv);

    auto & vcl = create_vcluster();
    size_t numCPU = vcl.getProcessingUnits();

    std::ofstream output_file;
    std::string filename = "Benchmark_OpenFPM_" + std::to_string(numCPU) + "core.csv";
    output_file.open(filename);

    output_file << numCPU << " cores,\n";
    output_file << "Time in ms,\n";

    for (int cycle = 0; cycle < 15; ++cycle) {

        int iteration = 0;

        double t = 0;

        // 3D physical domain
        Box<2,double> domain({0.0,0.0},{domainSize,domainSize});

        // Grid size on each dimension
        size_t sz[2] = {meshSize+1,meshSize+1};

        // Ghost part
        Ghost<2,double> g(r_cut);

        grid_dist_id<2, double, aggregate<double, double>> g_dist(sz,domain,g);

        std::vector<std::array<int, 2>> stencil;
        createStencil(stencil);

        // Get the iterator (No ghost)
        auto dom = g_dist.getDomainIterator();

        // Iterate over all the grid points
        while (dom.isNext()) {
            // local grid key from iterator
            auto key = dom.get();

            if (g_dist.getPos(key)[0] == domainSize / 2 && g_dist.getPos(key)[1] == domainSize / 2) {
                g_dist.get<0>(key) = 1;
            }
            else {
                g_dist.get<0>(key) = 0;

            }
            g_dist.get<1>(key) = 0;

            // next point
            ++dom;
        }

        g_dist.template ghost_get<0, 1>();

//        g_dist.write_frame("output", iteration, VTK_WRITER | FORMAT_ASCII);
        iteration++;

        auto t_start2 = std::chrono::high_resolution_clock::now();


        while (t < t_final) {


            // interact

            auto it1 = g_dist.getDomainIterator();
            while (it1.isNext()) {
                auto p = it1.get();

                // iterate through all neighbor mesh nodes
                for (std::array<int, 2> node: stencil) {
                    auto n = p;

                    // move to stencil position
                    for (int component = 0; component < 2; ++component) {
                        n = n.move(component, node[component]);
                    }

                    Point<2, double> p_pos = g_dist.getPos(p);
                    Point<2, double> n_pos = g_dist.getPos(n);
                    double distance2 = p_pos.distance2(n_pos);

                    double exchange = (g_dist.get<0>(n) - g_dist.get<0>(p))
                                      / (1 + pow(distance2 / epsilon / epsilon, 5));

                    g_dist.get<1>(p) += exchange;
                    g_dist.get<1>(n) -= exchange;


                }

                ++it1;
            }

            g_dist.template ghost_put<add_, 1>();


            // evolve

            auto it2 = g_dist.getDomainIterator();
            while (it2.isNext()) {
                auto p = it2.get();

//            std::cout << g_dist.get<1>(p) << std::endl;

                g_dist.get<0>(p) += g_dist.get<1>(p) * kernel;
                g_dist.get<1>(p) = 0;

                ++it2;
            }

            g_dist.template ghost_get<0, 1>();

//        if (iteration % 200 == 0)
//            g_dist.write_frame("output",iteration, VTK_WRITER | FORMAT_ASCII);

            iteration++;

            t += dt;
        }

        auto t_end2 = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double, std::milli> duration2 = t_end2 - t_start2;
        std::cout << std::endl << "Duration " << duration2.count() << std::endl;

        if (vcl.getProcessUnitID() == 0)
            output_file << (int)duration2.count() << ",\n";


    }


	openfpm_finalize();


}


