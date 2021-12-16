#ifdef __NVCC__

#include "Vector/map_vector.hpp"
#include "util/stat/common_statistics.hpp"

//! Memory bandwidth with small calculations
template<typename vector_type, typename vector_type2>
inline __global__ void translate_fill_prop_write(vector_type vd_out, vector_type2 vd_in)
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

	vd_in.template get<0>(p)[0] += a;
	vd_in.template get<0>(p)[1] += b;
}


template<typename vector_type, typename vector_type2>
inline __global__ void translate_fill_prop_read(vector_type vd_out, vector_type2 vd_in)
{
	auto p = blockIdx.x * blockDim.x + threadIdx.x;

	float a = vd_out.template get<0>(p);

	float b = vd_out.template get<1>(p)[0];
	float c = vd_out.template get<1>(p)[1];

	float d = vd_out.template get<2>(p)[0][0];
	float e = vd_out.template get<2>(p)[0][1];
	float f = vd_out.template get<2>(p)[1][0];
	float g = vd_out.template get<2>(p)[1][1];

	float h = vd_in.template get<0>(p)[0];
    float i = vd_in.template get<0>(p)[1];
    
	vd_in.template get<0>(p)[0] = a+b+c+d;
	vd_in.template get<0>(p)[1] = e+f+g+h+i;
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

    openfpm::vector<double> res;
    res.resize(100);

    for (int i = 0 ; i < 101 ; i++)
    {
        cudaDeviceSynchronize();
            timer t;
            t.start();


            CUDA_LAUNCH(translate_fill_prop_write,ite,out.toKernel(),in.toKernel());

            cudaDeviceSynchronize();

            t.stop();

        if (i >=1)
        {res.get(i-1) = nele*4*13 / t.getwct() * 1e-9;}

            std::cout << "Time: " << t.getwct() << std::endl;
            std::cout << "BW: " << nele*4*13 / t.getwct() * 1e-9 << " GB/s"  << std::endl;
    }

    double mean_write = 0.0;
    double dev_write = 0.0;
    standard_deviation(res,mean_write,dev_write);

    for (int i = 0 ; i < 101 ; i++)
    {
            cudaDeviceSynchronize();
            timer t;
            t.start();


            CUDA_LAUNCH(translate_fill_prop_read,ite,out.toKernel(),in.toKernel());

            cudaDeviceSynchronize();

            t.stop();

        if (i >=1)
        {res.get(i-1) = nele*4*11 / t.getwct() * 1e-9;}

            std::cout << "Time: " << t.getwct() << std::endl;
            std::cout << "BW: " << nele*4*11 / t.getwct() * 1e-9 << " GB/s"  << std::endl;
    }

    double mean_read = 0.0;
    double dev_read = 0.0;
    standard_deviation(res,mean_read,dev_read);

    std::cout << "Average READ: " << mean_read << "  deviation: " << dev_read << std::endl;
    std::cout << "Average WRITE: " << mean_write << "  deviation: " << dev_write << std::endl;
}

#else

int main(int argc, char *argv[])
{
}

#endif

