#ifdef __NVCC__

#include "Vector/map_vector.hpp"
#include "util/stat/common_statistics.hpp"

//! Memory bandwidth with small calculations
template<typename vector_type, typename vector_type2>
inline __global__ void translate_fill_prop(vector_type vd_out, vector_type2 vd_in)
{
	auto p = blockIdx.x * blockDim.x + threadIdx.x;

	float a = vd_in.template get<0>(p)[0];
	float b = vd_in.template get<0>(p)[1];

	vd_out.template get<0>(p) = a + b;

	vd_out.template get<1>(p)[0] = a;
	vd_out.template get<1>(p)[1] = b;

	vd_out.template get<2>(p)[0][0] = a;
	vd_out.template get<2>(p)[0][1] = b;
	vd_out.template get<2>(p)[1][0] = a + b;
	vd_out.template get<2>(p)[1][1] = b - a;

	vd_in.template get<0>(p)[0] += 0.01f;
	vd_in.template get<0>(p)[1] += 0.01f;
}


int main(int argc, char *argv[])
{
    init_wrappers();

    openfpm::vector_gpu<aggregate<double,double[2],double[2][2]>> out;
    openfpm::vector_gpu<aggregate<double[2]>> in;

    int nele = 16777216;

    out.resize(nele);
    in.resize(nele);

    for (int i = 0 ; i < 16777216 ; i++)
    {
        in.template get<0>(i)[0] = i;
        in.template get<0>(i)[1] = i+100.0;
    }

    auto ite = out.getGPUIterator(256);

    openfpm::vector<double> res;
    res.resize(100);

    for (int i = 0 ; i < 101 ; i++)
    {
	cudaDeviceSynchronize();
        timer t;
        t.start();


        CUDA_LAUNCH(translate_fill_prop,ite,out.toKernel(),in.toKernel());

        cudaDeviceSynchronize();

        t.stop();

	if (i >=1)
	{res.get(i-1) = nele*8*11 / t.getwct() * 1e-9;}

        std::cout << "Time: " << t.getwct() << std::endl;
        std::cout << "BW: " << nele*8*11 / t.getwct() * 1e-9 << " GB/s"  << std::endl;
    }

    double mean = 0.0;
    double dev = 0.0;
    standard_deviation(res,mean,dev);

    std::cout << "Average: " << mean << "  deviation: " << dev << std::endl;
}

#else

int main(int argc, char *argv[])
{
}

#endif

