#include "Vector/vector_dist.hpp"
#include "OpenFPM_version.hpp"
#include <chrono>
#include <fstream>


double dt = 0.001;
double t_final = .5;
double domainSize = 9.0;
double damp = 0.9;
double r_cut = 0.08;
int number_particles = 200000;




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

        size_t bc[2]={PERIODIC,PERIODIC};

        // Ghost part
        Ghost<2,double> g(r_cut);

        vector_dist<2, double, aggregate<double[2], double[2]>> v_dist(number_particles,domain, bc,g);
        auto NN = v_dist.getCellList<CELL_MEMBAL(2,double)>(r_cut);


        std::random_device rd;  // Will be used to obtain a seed for the random number engine
        std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
        std::uniform_real_distribution<> dis_pos(0, domainSize);
        std::normal_distribution<> dis_vel(0, 1);

        // Get the iterator (No ghost)
        auto dom = v_dist.getDomainIterator();
        // Iterate over all the grid points
        while (dom.isNext()) {
            // local grid key from iterator
            auto key = dom.get();

            v_dist.getPos(key)[0] = dis_pos(gen);
            v_dist.getPos(key)[1] = dis_pos(gen);
            v_dist.getProp<0>(key)[0] = dis_vel(gen);
            v_dist.getProp<0>(key)[1] = dis_vel(gen);
            v_dist.getProp<1>(key)[0] = 0;
            v_dist.getProp<1>(key)[1] = 0;

            // next point
            ++dom;
        }

        v_dist.map();
//        v_dist.deleteGhost();
//        v_dist.write_frame("output", iteration, VTK_WRITER | FORMAT_ASCII);
        v_dist.template ghost_get<0, 1>();

        iteration++;

        auto t_start2 = std::chrono::high_resolution_clock::now();


        while (t < t_final) {

            // interact

            v_dist.updateCellList(NN);

            auto it1 = v_dist.getDomainIterator();
            while (it1.isNext()) {
                auto p = it1.get();

                // iterate through all neighbor mesh nodes

                auto Np = NN.template getNNIterator<NO_CHECK>(NN.getCell(v_dist.getPos(p)));

                while (Np.isNext()) {
                    auto n = Np.get();

                    if (n == p.getKey()) {
                        ++Np;
                        continue;
                    }

                    // Declare particle property variables
                    Point<2, double> p_pos = v_dist.getPos(p);
                    Point<2, double> n_pos = v_dist.getPos(n);
                    Point<2, double> p_vel = v_dist.getProp<0>(p);
                    Point<2, double> n_vel = v_dist.getProp<0>(n);

                    // Check cutoff radius
                    if (p_pos.distance(n_pos) > r_cut) {
                        ++Np;
                        continue;
                    }

                    // Check if particles are moving towards each other
                    Point<2, double> p_move = p_pos + (p_vel * dt);
                    Point<2, double> n_move = n_pos + (n_vel * dt);
                    float dist = p_pos.distance2(n_pos);
                    float dist_move = (p_move).distance2(n_move);
                    if (dist < dist_move) {
                        ++Np;
                        continue;
                    }

                    // Compute collision vector
                    Point<2, double> diff = n_pos - p_pos;
                    Point<2, double> diff_scaled = diff / n_pos.distance2(p_pos);
                    Point<2, double> diff_vel = n_vel - p_vel;

                    Point<2, double> diff_collision = diff * (diff_scaled[0] * diff_vel[0] + diff_scaled[1] * diff_vel[1]);
                    diff_collision = diff_collision; //* damp;

                    // Apply collision to particle acceleration
                    v_dist.getProp<1>(p)[0] += diff_collision[0];
                    v_dist.getProp<1>(p)[1] += diff_collision[1];

                    ++Np;
                }

                ++it1;
            }


            // evolve

            auto it2 = v_dist.getDomainIterator();
            while (it2.isNext()) {
                auto p = it2.get();

                v_dist.getProp<0>(p)[0] += v_dist.getProp<1>(p)[0];
                v_dist.getProp<0>(p)[1] += v_dist.getProp<1>(p)[1];

                v_dist.getProp<1>(p)[0] = 0;
                v_dist.getProp<1>(p)[1] = 0;

                v_dist.getPos(p)[0] += v_dist.getProp<0>(p)[0] * dt;
                v_dist.getPos(p)[1] += v_dist.getProp<0>(p)[1] * dt;

                ++it2;
            }


//            if (iteration % 50 == 0) {
//                v_dist.deleteGhost();
//                v_dist.write_frame("output",iteration, VTK_WRITER | FORMAT_ASCII);
//            }

            v_dist.map();
            v_dist.template ghost_get<0, 1>();

            iteration++;

            t += dt;
        }

        auto t_end2 = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double, std::milli> duration2 = t_end2 - t_start2;
        std::cout << std::endl << "Duration " << duration2.count() << std::endl;

        if (vcl.getProcessUnitID() == 0)
            output_file << (int)duration2.count() << ",\n";


    }

    output_file.close();


	openfpm_finalize();


}


