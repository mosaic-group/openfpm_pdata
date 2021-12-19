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
    
	vd_in.template get<0>(p)[0] = a+b+c+d;
	vd_in.template get<0>(p)[1] = e+f+g;
}

/////////////////////////////// Lambda based

template<typename vector_type, typename vector_type2>
inline __device__ void translate_fill_prop_write_notls(vector_type vd_out, vector_type2 vd_in, dim3 & blockIdx, dim3 & threadIdx)
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

	vd_in.template get<0>(p)[0] = a;
	vd_in.template get<0>(p)[1] = b;
}


template<typename vector_type, typename vector_type2>
inline __device__ void translate_fill_prop_read_notls(vector_type vd_out, vector_type2 vd_in, dim3 & blockIdx, dim3 & threadIdx)
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

// Arrays

inline __global__ void translate_fill_prop_write_array(float * vd_out_scal,
                                                       float * vd_out_vec,
                                                       float * vd_out_mat,
                                                       float * vd_in_vec,
                                                       int stride)
{
	auto p = blockIdx.x * blockDim.x + threadIdx.x;

	float a = vd_in_vec[p* + 0*stride];
	float b = vd_in_vec[p* + 1*stride];

	vd_out_scal[p] = a + b;

	vd_out_vec[p + 0*stride] = a;
	vd_out_vec[p + 1*stride] = b;

	vd_out_mat[p + 0*2*stride + 0*stride ] = a;
	vd_out_mat[p + 0*2*stride + 1*stride ] = b;
	vd_out_mat[p + 1*2*stride + 0*stride ] = a + b;
	vd_out_mat[p + 1*2*stride + 1*stride ] = b - a;
}


template<typename vector_type, typename vector_type2>
inline __global__ void translate_fill_prop_read_array(vector_type vd_out, vector_type2 vd_in)
{
	auto p = blockIdx.x * blockDim.x + threadIdx.x;

	float a = vd_out.template get<0>(p);

	float b = vd_out.template get<1>(p)[0];
	float c = vd_out.template get<1>(p)[1];

	float d = vd_out.template get<2>(p)[0][0];
	float e = vd_out.template get<2>(p)[0][1];
	float f = vd_out.template get<2>(p)[1][0];
	float g = vd_out.template get<2>(p)[1][1];
    
	vd_in.template get<0>(p)[0] = a+b+c+d;
	vd_in.template get<0>(p)[1] = e+f+g;
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

    // Read write test with TLS

    auto ite = out.getGPUIterator(256);

    openfpm::vector<double> res;
    res.resize(100);

    for (int i = 0 ; i < 110 ; i++)
    {
        cudaDeviceSynchronize();
        timer t;
        t.start();


        CUDA_LAUNCH(translate_fill_prop_write,ite,out.toKernel(),in.toKernel());

        cudaDeviceSynchronize();

        t.stop();

        if (i >=10)
        {res.get(i-10) = nele*4*9 / t.getwct() * 1e-9;}

        std::cout << "Time: " << t.getwct() << std::endl;
        std::cout << "BW: " << nele*4*9 / t.getwct() * 1e-9 << " GB/s"  << std::endl;
    }

    double mean_write_tls = 0.0;
    double dev_write_tls = 0.0;
    standard_deviation(res,mean_write_tls,dev_write_tls);

    for (int i = 0 ; i < 110 ; i++)
    {
        cudaDeviceSynchronize();
        timer t;
        t.start();


        CUDA_LAUNCH(translate_fill_prop_read,ite,out.toKernel(),in.toKernel());

        cudaDeviceSynchronize();

        t.stop();

        if (i >=10)
        {res.get(i-10) = nele*4*9 / t.getwct() * 1e-9;}

        std::cout << "Time: " << t.getwct() << std::endl;
        std::cout << "BW: " << nele*4*9 / t.getwct() * 1e-9 << " GB/s"  << std::endl;
    }

    double mean_read_tls = 0.0;
    double dev_read_tls = 0.0;
    standard_deviation(res,mean_read_tls,dev_read_tls);

    /////////////////////////////////////////// LAMBDA //////////////////////////////////////////


    for (int i = 0 ; i < 110 ; i++)
    {
        cudaDeviceSynchronize();
        timer t;
        t.start();

        auto vd_out = out.toKernel();
        auto vd_in = in.toKernel();

        auto lamb = [&] __device__ (dim3 & blockIdx, dim3 & threadIdx)
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
        };

        CUDA_LAUNCH_LAMBDA(ite, lamb);

        cudaDeviceSynchronize();

        t.stop();

        if (i >=10)
        {res.get(i-10) = nele*4*9 / t.getwct() * 1e-9;}

        std::cout << "Time: " << t.getwct() << std::endl;
        std::cout << "BW: " << nele*4*9 / t.getwct() * 1e-9 << " GB/s"  << std::endl;
    }

    double mean_write_lamb = 0.0;
    double dev_write_lamb = 0.0;
    standard_deviation(res,mean_write_lamb,dev_write_lamb);

    for (int i = 0 ; i < 110 ; i++)
    {
        cudaDeviceSynchronize();
        timer t;
        t.start();


        auto vd_out = out.toKernel();
        auto vd_in = in.toKernel();

        auto lamb = [vd_out,vd_in] __device__ (dim3 & blockIdx, dim3 & threadIdx)
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
                            };

        CUDA_LAUNCH_LAMBDA(ite, lamb);

        cudaDeviceSynchronize();

        t.stop();

        if (i >=10)
        {res.get(i-10) = nele*4*9 / t.getwct() * 1e-9;}

        std::cout << "Time: " << t.getwct() << std::endl;
        std::cout << "BW: " << nele*4*9 / t.getwct() * 1e-9 << " GB/s"  << std::endl;
    }

    double mean_read_lamb = 0.0;
    double dev_read_lamb = 0.0;
    standard_deviation(res,mean_read_lamb,dev_read_lamb);

    std::cout << "Average READ with TLS: " << mean_read_tls << "  deviation: " << dev_read_tls << std::endl;
    std::cout << "Average WRITE with TLS: " << mean_write_tls << "  deviation: " << dev_write_tls << std::endl;

    std::cout << "Average READ with lamb: " << mean_read_lamb << "  deviation: " << dev_read_lamb << std::endl;
    std::cout << "Average WRITE with lamb: " << mean_write_lamb << "  deviation: " << dev_write_lamb << std::endl;
}

#else

int main(int argc, char *argv[])
{
}

#endif

