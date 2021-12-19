#ifdef __NVCC__

#include "Vector/map_vector.hpp"
#include "util/stat/common_statistics.hpp"

//! Memory bandwidth with small calculations
template<typename vector_type, typename vector_type2>
__global__ void translate_fill_prop_write(vector_type vd_out, vector_type2 vd_in)
{
	auto p = blockIdx.x * blockDim.x + threadIdx.x;

	float a = vd_in.template get<0>(p)[0];

	vd_out.template get<0>(p) = a;

	vd_out.template get<1>(p)[0] = a;
	vd_out.template get<1>(p)[1] = a;

	vd_out.template get<2>(p)[0][0] = a;
	vd_out.template get<2>(p)[0][1] = a;
	vd_out.template get<2>(p)[1][0] = a;
    vd_out.template get<2>(p)[1][1] = a;
    vd_in.template get<0>(p)[1] = a;
}


template<typename vector_type, typename vector_type2>
__global__ void translate_fill_prop_read(vector_type vd_out, vector_type2 vd_in)
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
	vd_in.template get<0>(p)[1] = a+b+c+d+e+f+g+h;
}

// Arrays

__global__ void translate_fill_prop_write_array(float * vd_out_scal,
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
__global__ void translate_fill_prop_read_array(vector_type vd_out, vector_type2 vd_in)
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
	vd_in.template get<0>(p)[1] = a+b+c+d+e+f+g+h;
}

template<typename in_type, typename out_type>
void check_write(in_type & in, out_type & out)
{
    out.template deviceToHost<0,1,2>();
    in.template deviceToHost<0>();

    bool success = true;
    for (int i = 0 ; i < 16777216 ; i++)
    {
        float a = in.template get<0>(i)[0];

        success &= out.template get<0>(i) == a;

        success &= out.template get<1>(i)[0] == a;
        success &= out.template get<1>(i)[1] == a;

        success &= out.template get<2>(i)[0][0] == a;
        success &= out.template get<2>(i)[0][1] == a;
        success &= out.template get<2>(i)[1][0] == a;
        success &= out.template get<2>(i)[1][1] == a;

        success &= in.template get<0>(i)[1] == a;
    }

    if (success == false)
    {
            std::cout << "FAIL WRITE" << std::endl;
            exit(1);
    }
}

template<typename in_type, typename out_type>
void check_read(in_type & in, out_type & out)
{
    out.template deviceToHost<0,1,2>();
    in.template deviceToHost<0>();

    bool success = true;
    for (int i = 0 ; i < 16777216 ; i++)
    {
        float a = out.template get<0>(i);

        float b = out.template get<1>(i)[0];
        float c = out.template get<1>(i)[1];

        float d = out.template get<2>(i)[0][0];
        float e = out.template get<2>(i)[0][1];
        float f = out.template get<2>(i)[1][0];
        float g = out.template get<2>(i)[1][1];

        float h = in.template get<0>(i)[0];

        success &= in.template get<0>(i)[1] == (a+b+c+d+e+f+g+h);

        if (success == false)
        {
            std::cout << "FAIL READ " << i << "   " << in.template get<0>(i)[1] << " != " << a+b+c+d+e+f+g+h << std::endl;
            exit(1);
        }
    }
}

template<typename vin_type, typename vout_type>
void initialize_buf(vin_type in, vout_type out)
{
    for (int i = 0 ; i < 16777216 ; i++)
    {
        in.template get<0>(i)[0] = i;
        in.template get<0>(i)[1] = i+100.0;

        out.template get<0>(i) = i+200.0;

        out.template get<1>(i)[0] = i;
        out.template get<1>(i)[1] = i+100.0;

        out.template get<2>(i)[0][0] = i;
        out.template get<2>(i)[0][1] = i+100.0;
        out.template get<2>(i)[1][0] = i+200.0;
        out.template get<2>(i)[1][1] = i+300.0;
    }

}

int main(int argc, char *argv[])
{
    init_wrappers();

    openfpm::vector_gpu<aggregate<float,float[2],float[2][2]>> out;
    openfpm::vector_gpu<aggregate<float[2]>> in;

    int nele = 16777216;

    out.resize(nele);
    in.resize(nele);

    initialize_buf(in,out);

    // Read write test with TLS

    auto ite = out.getGPUIterator(256);

    openfpm::vector<double> res;
    res.resize(100);

    in.hostToDevice<0>();

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

    check_write(in,out);

    initialize_buf(in,out);

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

    check_read(in,out);

    //////////////

    /////////////////////////////////////////// LAMBDA //////////////////////////////////////////

    initialize_buf(in,out);

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

            float a = vd_in.template get<0>(p)[0];
            float b = vd_in.template get<0>(p)[1];

	    vd_out.template get<0>(p) = a + b;

            vd_out.template get<1>(p)[0] = a;
            vd_out.template get<1>(p)[1] = b;

            vd_out.template get<2>(p)[0][0] = a;
            vd_out.template get<2>(p)[0][1] = b;
            vd_out.template get<2>(p)[1][0] = a + b;
            vd_out.template get<2>(p)[1][1] = b - a;
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

    initialize_buf(in,out);

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
                                
                                vd_in.template get<0>(p)[0] = a+b+c+d;
                                vd_in.template get<0>(p)[1] = e+f+g;
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

    #ifdef CUDIFY_USE_CUDA

    for (int i = 0 ; i < 110 ; i++)
    {
        cudaDeviceSynchronize();
        timer t;
        t.start();

        float * a = (float *)in.getDeviceBuffer<0>();
        float * b = (float *)out.getDeviceBuffer<1>();

        cudaMemcpy(a,b,2*16777216*4,cudaMemcpyDeviceToDevice);

        cudaDeviceSynchronize();

        t.stop();

        if (i >=10)
        {res.get(i-10) = nele*4*4 / t.getwct() * 1e-9;}

        std::cout << "Time: " << t.getwct() << std::endl;
        std::cout << "BW: " << nele*4*4 / t.getwct() * 1e-9 << " GB/s"  << std::endl;
    }    

    double mean_read_mes = 0.0;
    double dev_read_mes = 0.0;
    standard_deviation(res,mean_read_mes,dev_read_mes);

    std::cout << "Average measured: " << mean_read_mes << "  deviation: " << dev_read_mes << std::endl;

    #endif

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

