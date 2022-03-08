#ifdef __NVCC__

#include "Vector/map_vector.hpp"
#include "util/stat/common_statistics.hpp"

#define NELEMENTS 67108864

//! Memory bandwidth with small calculations
template<typename vector_type, typename vector_type2>
__global__ void translate_fill_prop_write(vector_type vd_out, vector_type2 vd_in)
{
    grid_key_dx<3,int> p({blockIdx.x * blockDim.x + threadIdx.x,
                          blockIdx.y * blockDim.y + threadIdx.y,
                          blockIdx.z * blockDim.z + threadIdx.z});

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
    grid_key_dx<3,int> p({blockIdx.x * blockDim.x + threadIdx.x,
        blockIdx.y * blockDim.y + threadIdx.y,
        blockIdx.z * blockDim.z + threadIdx.z});

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
    auto it = in.getIterator();
    while (it.isNext())
    {
        auto i = it.get();

        float a = in.template get<0>(i)[0];

        if (i.get(0) == 2 && i.get(1) == 2 && i.get(2) == 2)
        {
            success &= a != 0;
        }

        success &= out.template get<0>(i) == a;

        success &= out.template get<1>(i)[0] == a;
        success &= out.template get<1>(i)[1] == a;

        success &= out.template get<2>(i)[0][0] == a;
        success &= out.template get<2>(i)[0][1] == a;
        success &= out.template get<2>(i)[1][0] == a;
        success &= out.template get<2>(i)[1][1] == a;

        success &= in.template get<0>(i)[1] == a;

        if (success == false)
        {
            std::cout << "FAIL " << a << "   " << i.to_string() << "  " << out.template get<0>(i) << " " << out.template get<1>(i)[0] 
                                                          << out.template get<1>(i)[1] << " " << out.template get<2>(i)[0][0]
                                                          << out.template get<2>(i)[0][1] << " " << out.template get<2>(i)[1][0]
                                                          << out.template get<2>(i)[1][0] << " " << out.template get<2>(i)[1][1]
                                                        << std::endl;
            break;
        }

        ++it;
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
    auto it = in.getIterator();
    while (it.isNext())
    {
        auto i = it.get();

        float a = out.template get<0>(i);

        if (i.get(0) == 2 && i.get(1) == 2 && i.get(2) == 2)
        {
            success &= a != 0;
        }

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
            std::cout << "FAIL READ " << i.to_string() << "   " << in.template get<0>(i)[1] << " != " << a+b+c+d+e+f+g+h << std::endl;
            exit(1);
        }

        ++it;
    }
}

template<typename vector_type, typename vector_type2>
__global__ void initialize_buff(vector_type vd_out, vector_type2 vd_in)
{
    grid_key_dx<3,int> i({blockIdx.x * blockDim.x + threadIdx.x,
        blockIdx.y * blockDim.y + threadIdx.y,
        blockIdx.z * blockDim.z + threadIdx.z});

    vd_in.template get<0>(i)[0] = i.get(0) + i.get(1) + i.get(2);
    vd_in.template get<0>(i)[1] = i.get(0) + i.get(1) + i.get(2)+100.0;

    vd_out.template get<0>(i) = i.get(0) + i.get(1) + i.get(2)+200.0;

    vd_out.template get<1>(i)[0] = i.get(0) + i.get(1) + i.get(2);
    vd_out.template get<1>(i)[1] = i.get(0) + i.get(1) + i.get(2)+100.0;

    vd_out.template get<2>(i)[0][0] = i.get(0) + i.get(1) + i.get(2);
    vd_out.template get<2>(i)[0][1] = i.get(0) + i.get(1) + i.get(2)+100.0;
    vd_out.template get<2>(i)[1][0] = i.get(0) + i.get(1) + i.get(2)+200.0;
    vd_out.template get<2>(i)[1][1] = i.get(0) + i.get(1) + i.get(2)+300.0;
}

template<typename vin_type, typename vout_type>
void initialize_buf(vin_type & in, vout_type & out)
{
    auto ite = out.getGPUIterator({0,0,0},{511,511,255});
    CUDA_LAUNCH(initialize_buff,ite,out.toKernel(),in.toKernel());
}

int main(int argc, char *argv[])
{
    init_wrappers();

    grid_gpu<3,aggregate<float,float[2],float[2][2]>> out;
    grid_gpu<3,aggregate<float[2]>> in;

    int nele = NELEMENTS;

    size_t sz[3] = {512,512,256};
    out.resize(sz);
    in.resize(sz);

    out.setMemory();
    in.setMemory();

    initialize_buf(in,out);

    // Read write test with TLS

    auto ite = out.getGPUIterator({0,0,0},{511,511,255});

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
        {res.get(i-10) = (double)nele*4*9 / t.getwct() * 1e-9;}

        std::cout << "Time: " << t.getwct() << std::endl;
        std::cout << "BW: " << (double)nele*4*9 / t.getwct() * 1e-9 << " GB/s"  << std::endl;
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
        {res.get(i-10) = (double)nele*4*9 / t.getwct() * 1e-9;}

        std::cout << "Time: " << t.getwct() << std::endl;
        std::cout << "BW: " << (double)nele*4*9 / t.getwct() * 1e-9 << " GB/s"  << std::endl;
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
            grid_key_dx<3,int> p({blockIdx.x * blockDim.x + threadIdx.x,
                blockIdx.y * blockDim.y + threadIdx.y,
                blockIdx.z * blockDim.z + threadIdx.z});

            float a = vd_in.template get<0>(p)[0];

            vd_out.template get<0>(p) = a;

            vd_out.template get<1>(p)[0] = a;
            vd_out.template get<1>(p)[1] = a;
        
            vd_out.template get<2>(p)[0][0] = a;
            vd_out.template get<2>(p)[0][1] = a;
            vd_out.template get<2>(p)[1][0] = a;
            vd_out.template get<2>(p)[1][1] = a;
            vd_in.template get<0>(p)[1] = a;
        };

        CUDA_LAUNCH_LAMBDA(ite, lamb);

        cudaDeviceSynchronize();

        t.stop();

        if (i >=10)
        {res.get(i-10) = (double)nele*4*9 / t.getwct() * 1e-9;}

        std::cout << "Time: " << t.getwct() << std::endl;
        std::cout << "BW: " << (double)nele*4*9 / t.getwct() * 1e-9 << " GB/s"  << std::endl;
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
                                grid_key_dx<3,int> p({blockIdx.x * blockDim.x + threadIdx.x,
                                    blockIdx.y * blockDim.y + threadIdx.y,
                                    blockIdx.z * blockDim.z + threadIdx.z});

                                float a = vd_out.template get<0>(p);

                                float b = vd_out.template get<1>(p)[0];
                                float c = vd_out.template get<1>(p)[1];
                            
                                float d = vd_out.template get<2>(p)[0][0];
                                float e = vd_out.template get<2>(p)[0][1];
                                float f = vd_out.template get<2>(p)[1][0];
                                float g = vd_out.template get<2>(p)[1][1];
                                
                                float h = vd_in.template get<0>(p)[0];
                                vd_in.template get<0>(p)[1] = a+b+c+d+e+f+g+h;
                            };

        CUDA_LAUNCH_LAMBDA(ite, lamb);

        cudaDeviceSynchronize();

        t.stop();

        if (i >=10)
        {res.get(i-10) = (double)nele*4*9 / t.getwct() * 1e-9;}

        std::cout << "Time: " << t.getwct() << std::endl;
        std::cout << "BW: " << (double)nele*4*9 / t.getwct() * 1e-9 << " GB/s"  << std::endl;
    }

    double mean_read_lamb = 0.0;
    double dev_read_lamb = 0.0;
    standard_deviation(res,mean_read_lamb,dev_read_lamb);

    // Array benchmark
    initialize_buf(in,out);

    for (int i = 0 ; i < 110 ; i++)
    {
        cudaDeviceSynchronize();
        timer t;
        t.start();

	    float * out_s = (float *)out.getDeviceBuffer<0>();
	    float * out_v = (float *)out.getDeviceBuffer<1>();
	    float * out_m = (float *)out.getDeviceBuffer<2>();
        float * in_v = (float *)in.getDeviceBuffer<0>();
        
        int sz0 = sz[0];
        int sz1 = sz[1];
        int sz2 = sz[2];
        int stride = out.size();

        auto lamb_arr_write = [out_s,out_v,out_m,in_v,sz0,sz1,sz2,stride] __device__ (dim3 & blockIdx, dim3 & threadIdx)
        {
            auto p1 = blockIdx.x * blockDim.x + threadIdx.x;
            auto p2 = blockIdx.y * blockDim.y + threadIdx.y;
            auto p3 = blockIdx.z * blockDim.z + threadIdx.z;

            float a = in_v[p1 + p2*sz0 + p3*sz0*sz1 + 0*stride];

            out_s[p1 + p2*sz0 + p3*sz0*sz1] = a;
        
            out_v[p1 + p2*sz0 + p3*sz0*sz1 + 0*stride] = a;
            out_v[p1 + p2*sz0 + p3*sz0*sz1 + 1*stride] = a;
        
            out_m[p1 + p2*sz0 + p3*sz0*sz1 + 0*2*stride + 0*stride ] = a;
            out_m[p1 + p2*sz0 + p3*sz0*sz1 + 0*2*stride + 1*stride ] = a;
            out_m[p1 + p2*sz0 + p3*sz0*sz1 + 1*2*stride + 0*stride ] = a;
            out_m[p1 + p2*sz0 + p3*sz0*sz1 + 1*2*stride + 1*stride ] = a;
            in_v[p1 + p2*sz0 + p3*sz0*sz1 + 1*stride] = a;
        };

        CUDA_LAUNCH_LAMBDA(ite,lamb_arr_write);

        cudaDeviceSynchronize();

        t.stop();

        if (i >=10)
        {res.get(i-10) = (double)nele*4*9 / t.getwct() * 1e-9;}

        std::cout << "Time ARR: " << t.getwct() << std::endl;
        std::cout << "BW ARR: " << (double)nele*4*9 / t.getwct() * 1e-9 << " GB/s"  << std::endl;
    }

    double mean_write_arr = 0.0;
    double dev_write_arr = 0.0;
    standard_deviation(res,mean_write_arr,dev_write_arr);

    check_write(in,out);

    for (int i = 0 ; i < 110 ; i++)
    {
        cudaDeviceSynchronize();
        timer t;
        t.start();

	    float * out_s = (float *)out.getDeviceBuffer<0>();
	    float * out_v = (float *)out.getDeviceBuffer<1>();
	    float * out_m = (float *)out.getDeviceBuffer<2>();
        float * in_v = (float *)in.getDeviceBuffer<0>();
        
        int sz0 = sz[0];
        int sz1 = sz[1];
        int sz2 = sz[2];
        int stride = out.size();

        auto lamb_arr_red = [out_s,out_v,out_m,in_v,sz0,sz1,sz2,stride] __device__ (dim3 & blockIdx, dim3 & threadIdx)
        {
            auto p1 = blockIdx.x * blockDim.x + threadIdx.x;
            auto p2 = blockIdx.y * blockDim.y + threadIdx.y;
            auto p3 = blockIdx.z * blockDim.z + threadIdx.z;

            float a = out_s[p1 + p2*sz0 + p3*sz0*sz1];
        
            float b = out_v[p1 + p2*sz0 + p3*sz0*sz1 + 0*stride];
            float c = out_v[p1 + p2*sz0 + p3*sz0*sz1 + 1*stride];
        
            float d = out_m[p1 + p2*sz0 + p3*sz0*sz1 + 0*2*stride + 0*stride];
            float e = out_m[p1 + p2*sz0 + p3*sz0*sz1 + 0*2*stride + 1*stride];
            float f = out_m[p1 + p2*sz0 + p3*sz0*sz1 + 1*2*stride + 0*stride];
            float g = out_m[p1 + p2*sz0 + p3*sz0*sz1 + 1*2*stride + 1*stride];
            
            float h = in_v[p1 + p2*sz0 + p3*sz0*sz1 + 0*stride];
            in_v[p1 + p2*sz0 + p3*sz0*sz1 + 1*stride] = a+b+c+d+e+f+g+h;
        };

        CUDA_LAUNCH_LAMBDA(ite,lamb_arr_red);

        cudaDeviceSynchronize();

        t.stop();

        if (i >=10)
        {res.get(i-10) = (double)nele*4*9 / t.getwct() * 1e-9;}

        std::cout << "Time ARR: " << t.getwct() << std::endl;
        std::cout << "BW ARR: " << (double)nele*4*9 / t.getwct() * 1e-9 << " GB/s"  << std::endl;
    }

    double mean_read_arr = 0.0;
    double dev_read_arr = 0.0;
    standard_deviation(res,mean_read_arr,dev_read_arr);

    check_read(in,out);

    ///////////////////

    #ifdef CUDIFY_USE_CUDA

    for (int i = 0 ; i < 110 ; i++)
    {
        cudaDeviceSynchronize();
        timer t;
        t.start();

        float * a = (float *)in.getDeviceBuffer<0>();
        float * b = (float *)out.getDeviceBuffer<1>();

        cudaMemcpy(a,b,2*NELEMENTS*4,cudaMemcpyDeviceToDevice);

        cudaDeviceSynchronize();

        t.stop();

        if (i >=10)
        {res.get(i-10) = (double)nele*4*4 / t.getwct() * 1e-9;}

        std::cout << "Time: " << t.getwct() << std::endl;
        std::cout << "BW: " << (double)nele*4*4 / t.getwct() * 1e-9 << " GB/s"  << std::endl;
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

    std::cout << "Average WRITE with array: " << mean_write_arr << "  deviation: " << dev_write_arr << std::endl;
    std::cout << "Average READ with array: " << mean_read_arr << "  deviation: " << dev_read_arr << std::endl;
}

#else

int main(int argc, char *argv[])
{
}

#endif

