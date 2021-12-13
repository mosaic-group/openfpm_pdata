#include "Vector/map_vector.hpp"
#include "util/stat/common_statistics.hpp"


template<typename vector_type, typename vector_type2>
__attribute__((always_inline)) inline __global__ void translate_fill_prop(vector_type & vd_out, vector_type2 & vd_in)
{
	auto p = blockIdx.x * blockDim.x + threadIdx.x;

	vd_out.template get<0>(p) = vd_in.template get<0>(p)[0] + vd_in.template get<0>(p)[1];

	vd_out.template get<1>(p)[0] = vd_in.template get<0>(p)[0];
	vd_out.template get<1>(p)[1] = vd_in.template get<0>(p)[1];

	vd_out.template get<2>(p)[0][0] = vd_in.template get<0>(p)[0];
	vd_out.template get<2>(p)[0][1] = vd_in.template get<0>(p)[1];
	vd_out.template get<2>(p)[1][0] = vd_in.template get<0>(p)[0] + vd_in.template get<0>(p)[1];
	vd_out.template get<2>(p)[1][1] = vd_in.template get<0>(p)[1] - vd_in.template get<0>(p)[0];

	vd_in.template get<0>(p)[0] += 0.01f;
	vd_in.template get<0>(p)[1] += 0.01f;
}


int main(int argc, char *argv[])
{
    init_wrappers();

    openfpm::vector_gpu<aggregate<float,float[2],float[2][2]>> out;
    openfpm::vector_gpu<aggregate<float[2]>> in;

    int nele = 16777216;

    out.resize(nele);
    in.resize(nele);

    for (int i = 0 ; i < 16777216 ; i++)
    {
        in.template get<0>(i)[0] = i;
        in.template get<0>(i)[1] = i+100.0;
    }

    auto ite = out.getGPUIterator(256);

    for (int i = 0 ; i < 100 ; i++)
    {
	cudaDeviceSynchronize();
        timer t;
        t.start();

	auto vout = out.toKernel();
	auto vin = in.toKernel();

        CUDA_LAUNCH(translate_fill_prop,ite,vout,vin);

        cudaDeviceSynchronize();

        t.stop();
        std::cout << "Time: " << t.getwct() << std::endl;
        std::cout << "BW: " << nele*4*19 / t.getwct() * 1e-9 << " GB/s"  << std::endl;
    }
}
