# SET ENVIRONMENT VARIABLES ACCORDING TO ONE OF THE SUPPORTED MODES (cpu, gpu-nvcc, gpu-hip, gpu-emulated (CUDA_ON_CPU))

### This is a trick to avoid "Command not found if you no not have NVCC compiler". In practice the normal C++ compiler is used
### internally the example disable with the preprocessor its code if not compiled with nvcc 
CUDA_CC=
CUDA_CC_LINK=

ifdef HIP
    CUDA_CC=hipcc
    CUDA_OPTIONS=-D__NVCC__ -D__HIP__ -DCUDART_VERSION=11000 -D__CUDACC__ -D__CUDACC_VER_MAJOR__=11 -D__CUDACC_VER_MINOR__=0 -D__CUDACC_VER_BUILD__=0
    LIBS_SELECT=$(LIBS)
    CC=hipcc
    CUDA_CC_LINK=hipcc
else
    CC=mpic++
    ifdef CUDA_ON_CPU
        CUDA_CC=mpic++ -x c++ $(INCLUDE_PATH)
        INCLUDE_PATH_NVCC=
        CUDA_CC_LINK=mpic++
        CUDA_OPTIONS=-D__NVCC__ -DCUDART_VERSION=11000
        LIBS_SELECT=$(LIBS)
    else
        ifeq (, $(shell which nvcc))
            CUDA_CC=mpic++ -x c++ $(INCLUDE_PATH)
            INCLUDE_PATH_NVCC=
            CUDA_CC_LINK=mpic++
            CUDA_OPTIONS=
            LIBS_SELECT=$(LIBS)
        else
            CUDA_CC=nvcc -ccbin=mpic++
            CUDA_CC_LINK=nvcc -ccbin=mpic++
            CUDA_OPTIONS=-use_fast_math  -arch=sm_61 -lineinfo --extended-lambda --expt-relaxed-constexpr 
            LIBS_SELECT=$(LIBS_NVCC)
        endif
    endif
endif


ifeq ($(PROFILE),ON)
    CUDA_CC=scorep --nocompiler  --cuda --mpp=mpi nvcc -ccbin=mpic++
    CUDA_CC_LINK=scorep --nocompiler  --cuda --mpp=mpi nvcc -ccbin=mpic++
else
    CUDA_CC:=$(CUDA_CC)
    CUDA_CC_LINK:=$(CUDA_CC_LINK)
endif

LDIR =
OPT=

INCLUDE_PATH_NVCC:=$(INCLUDE_PATH_NVCC) -I./include 
