#ifdef __NVCC__

/**
 * BUDE CUDA kernel file
 **/

#define CUDIFY_BOOST_CONTEXT_STACK_SIZE 32768

//#define SE_CLASS1

#include <float.h>
#include <stdio.h>
#include <sys/time.h>

#include "Vector/map_vector.hpp"
#include "util/stat/common_statistics.hpp"

//#define USE_SHARED

constexpr int pos = 0;
constexpr int ind = 1;

constexpr int x = 0;
constexpr int y = 1;
constexpr int z = 2;

constexpr int   hbtype = 0;
constexpr int radius = 1;
constexpr int hphb = 2;
constexpr int elsc = 3;

#ifndef NUM_TD_PER_THREAD
// Good for CPU
//#define NUM_TD_PER_THREAD 256
// Good for GPU
#define NUM_TD_PER_THREAD 4
#endif

typedef struct
{
    float x, y, z;
    int   index;
} Atom;

typedef struct
{
    int   hbtype;
    float radius;
    float hphb;
    float elsc;
} FFParams;

typedef struct
{
    int    natlig;
    int    natpro;
    int    ntypes;
    int    nposes;
    char     * deckDir;
    int iterations;
} Params;

Params params;

typedef struct
{
  // _lin = AOS
  openfpm::vector_gpu_lin<aggregate<float[3],int>> d_protein;
  // AOS
  openfpm::vector_gpu_lin<aggregate<float[3],int>> d_ligand;
  // AOS
  openfpm::vector_gpu_lin<aggregate<int,float,float,float>> d_forcefield;
  // SOA
  openfpm::vector_gpu<aggregate<float>> d_results;
  // SOA
  openfpm::vector_gpu<aggregate<float,float,float,float,float,float>> d_poses;
  openfpm::vector<double> gflops_data;

    int deviceIndex;
    int wgsize;
    int posesPerWI;
} OpenFPM;


double getTimestamp()
{
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return tv.tv_usec + tv.tv_sec*1e6;
}

void printTimings(double start, double end, double poses_per_wi, openfpm::vector<double> & gflops_data)
{
    double ms = ((end-start)/params.iterations)*1e-3;

    // Compute FLOP/s
    double runtime   = ms*1e-3;
    double ops_per_wi = 27*poses_per_wi
        + params.natlig*(3 + 18*poses_per_wi + params.natpro*(11 + 30*poses_per_wi))
        + poses_per_wi;
    double total_ops     = ops_per_wi * (params.nposes/poses_per_wi);
    double flops      = total_ops / runtime;
    double gflops     = flops / 1e9;

    gflops_data.add(gflops);

    double interactions         =
        (double)params.nposes
        * (double)params.natlig
        * (double)params.natpro;
    double interactions_per_sec = interactions / runtime;

    // Print stats
    printf("- Total time:     %7.2lf ms\n", (end-start)*1e-3);
    printf("- Average time:   %7.2lf ms\n", ms);
    printf("- Interactions/s: %7.2lf billion\n", (interactions_per_sec / 1e9));
    printf("- GFLOP/s:        %7.2lf\n", gflops);
}

// Numeric constants
#define ZERO    0.0f
#define QUARTER 0.25f
#define HALF    0.5f
#define ONE     1.0f
#define TWO     2.0f
#define FOUR    4.0f
#define CNSTNT 45.0f

#define HBTYPE_F 70
#define HBTYPE_E 69

// The data structure for one atom - 16 bytes

typedef struct
{
  float x, y, z, w;
} Transform;

#define HARDNESS 38.0f
#define NPNPDIST  5.5f
#define NPPDIST   1.0f

__device__ inline void compute_transformation_matrix(const float transform_0,
    const float transform_1,
    const float transform_2,
    const float transform_3,
    const float transform_4,
    const float transform_5,
    Transform* transform)
{
  const float sx = sin(transform_0);
  const float cx = cos(transform_0);
  const float sy = sin(transform_1);
  const float cy = cos(transform_1);
  const float sz = sin(transform_2);
  const float cz = cos(transform_2);

  transform[0].x = cy*cz;
  transform[0].y = sx*sy*cz - cx*sz;
  transform[0].z = cx*sy*cz + sx*sz;
  transform[0].w = transform_3;
  transform[1].x = cy*sz;
  transform[1].y = sx*sy*sz + cx*cz;
  transform[1].z = cx*sy*sz - sx*cz;
  transform[1].w = transform_4;
  transform[2].x = -sy;
  transform[2].y = sx*cy;
  transform[2].z = cx*cy;
  transform[2].w = transform_5;
}



template<typename vector_atom, typename vector_ff, typename vector_tr, typename vector_out>
__global__ void fasten_main(const int natlig,
    const int natpro,
    const vector_atom protein_molecule,
    const vector_atom ligand_molecule,
    const vector_tr transforms,
    vector_out etotals,
    const vector_ff global_forcefield,
    const int num_atom_types,
    const int numTransforms)
{
  // Get index of first TD
  int ix = blockIdx.x*blockDim.x*NUM_TD_PER_THREAD + threadIdx.x;

  // Have extra threads do the last member intead of return.
  // A return would disable use of barriers, so not using return is better
  ix = ix < numTransforms ? ix : numTransforms - NUM_TD_PER_THREAD;

#ifdef USE_SHARED
  __shared__ FFParams forcefield[100];
  if(ix < num_atom_types)
  {
    forcefield[ix].hbtype = global_forcefield.template get<hbtype>(ix);
    forcefield[ix].radius = global_forcefield.template get<radius>(ix);
    forcefield[ix].hphb = global_forcefield.template get<hphb>(ix);
    forcefield[ix].elsc = global_forcefield.template get<elsc>(ix);

  }
#else
#endif

  // Compute transformation matrix to private memory
  float etot[NUM_TD_PER_THREAD];
  Transform transform[NUM_TD_PER_THREAD][3];
  const int lsz = blockDim.x;
  #pragma omp simd
  for (int i = 0; i < NUM_TD_PER_THREAD; i++)
  {
    int index = ix + i*lsz;
    compute_transformation_matrix(
        transforms.template get<0>(index),
        transforms.template get<1>(index),
        transforms.template get<2>(index),
        transforms.template get<3>(index),
        transforms.template get<4>(index),
        transforms.template get<5>(index),
        transform[i]);
    etot[i] = ZERO;
  }

#ifdef USE_SHARED
  __syncthreads();
#endif

  // Loop over ligand atoms
  int il = 0;
  do
  {
    // Load ligand atom data
    const Atom l_atom = {ligand_molecule.template get<pos>(il)[x],
    			 ligand_molecule.template get<pos>(il)[y],
    			 ligand_molecule.template get<pos>(il)[z],
    			 ligand_molecule.template get<ind>(il)}; 

    const FFParams l_params = {global_forcefield.template get<hbtype>(l_atom.index),
    			       global_forcefield.template get<radius>(l_atom.index),
    			       global_forcefield.template get<hphb>(l_atom.index),
    			       global_forcefield.template get<elsc>(l_atom.index)};
    const bool lhphb_ltz = l_params.hphb<ZERO;
    const bool lhphb_gtz = l_params.hphb>ZERO;

    float lpos_x[NUM_TD_PER_THREAD];
    float lpos_y[NUM_TD_PER_THREAD];
    float lpos_z[NUM_TD_PER_THREAD];
    const float4 linitpos = make_float4(l_atom.x,l_atom.y,l_atom.z,ONE);
    #pragma omp simd
    for (int i = 0; i < NUM_TD_PER_THREAD; i++)
    {
      // Transform ligand atom
      lpos_x[i] = transform[i][0].w + linitpos.x*transform[i][0].x + 
        linitpos.y*transform[i][0].y + linitpos.z*transform[i][0].z;
      lpos_y[i] = transform[i][1].w + linitpos.x*transform[i][1].x + 
        linitpos.y*transform[i][1].y + linitpos.z*transform[i][1].z;
      lpos_z[i] = transform[i][2].w + linitpos.x*transform[i][2].x + 
        linitpos.y*transform[i][2].y + linitpos.z*transform[i][2].z;
    }

    // Loop over protein atoms
    int ip = 0;
    do
    {
      // Load protein atom data
      const Atom p_atom = {protein_molecule.template get<pos>(ip)[x],
      			   protein_molecule.template get<pos>(ip)[y],
      			   protein_molecule.template get<pos>(ip)[z],
      			   protein_molecule.template get<ind>(ip)};

      const FFParams p_params = {global_forcefield.template get<hbtype>(p_atom.index),
	      			 global_forcefield.template get<radius>(p_atom.index),
				 global_forcefield.template get<hphb>(p_atom.index),
				 global_forcefield.template get<elsc>(p_atom.index)};

      const float radij   = p_params.radius + l_params.radius;
      const float r_radij = 1.0f/radij;

      const float elcdst  = (p_params.hbtype==HBTYPE_F && l_params.hbtype==HBTYPE_F) ? FOUR    : TWO;
      const float elcdst1 = (p_params.hbtype==HBTYPE_F && l_params.hbtype==HBTYPE_F) ? QUARTER : HALF;
      const bool type_E   = ((p_params.hbtype==HBTYPE_E || l_params.hbtype==HBTYPE_E));

      const bool phphb_ltz = p_params.hphb<ZERO;
      const bool phphb_gtz = p_params.hphb>ZERO;
      const bool phphb_nz  = p_params.hphb!=ZERO;
      const float p_hphb   = p_params.hphb * (phphb_ltz && lhphb_gtz ? -ONE : ONE);
      const float l_hphb   = l_params.hphb * (phphb_gtz && lhphb_ltz ? -ONE : ONE);
      const float distdslv = (phphb_ltz ? (lhphb_ltz ? NPNPDIST : NPPDIST) : (lhphb_ltz ? NPPDIST : -FLT_MAX) );

      float r_distdslv = 1.0f/distdslv;

      const float chrg_init = l_params.elsc * p_params.elsc;
      const float dslv_init = p_hphb + l_hphb;

      #pragma omp simd
      for (int i = 0; i < NUM_TD_PER_THREAD; i++)
      {
        // Calculate distance between atoms
        const float x      = lpos_x[i] - p_atom.x;
        const float y      = lpos_y[i] - p_atom.y;
        const float z      = lpos_z[i] - p_atom.z;
        const float distij = sqrtf(x*x + y*y + z*z);

        // Calculate the sum of the sphere radii
        const float distbb = distij - radij;
        const bool  zone1  = (distbb < ZERO);

        // Calculate steric energy
        etot[i] += (ONE - (distij*r_radij)) * (zone1 ? 2*HARDNESS : ZERO);

        // Calculate formal and dipole charge interactions
        float chrg_e = chrg_init * ((zone1 ? 1 : (ONE - distbb*elcdst1)) 
            * (distbb<elcdst ? 1 : ZERO));
        const float neg_chrg_e = -fabs(chrg_e);
        chrg_e = type_E ? neg_chrg_e : chrg_e;
        etot[i] += chrg_e*CNSTNT;

        // Calculate the two cases for Nonpolar-Polar repulsive interactions
        const float coeff  = (ONE - (distbb *r_distdslv));
        float dslv_e = dslv_init * ((distbb<distdslv && phphb_nz) ? 1 : ZERO);
        dslv_e *= (zone1 ? 1 : coeff);
        etot[i] += dslv_e;
      }
    } 
    while (++ip < natpro); // loop over protein atoms
  } 
  while (++il < natlig); // loop over ligand atoms

  // Write results
  const int td_base = blockIdx.x*blockDim.x*NUM_TD_PER_THREAD + threadIdx.x;
  if (td_base < numTransforms)
  {
    #pragma omp simd
    for (int i = 0; i < NUM_TD_PER_THREAD; i++)
    {
      etotals.template get<0>(td_base+i*blockDim.x) = etot[i]*HALF;
    }
  }
} //end of fasten_main


void runCUDA(OpenFPM & _openfpm)
{
  _openfpm.d_protein.hostToDevice<pos,ind>();
  _openfpm.d_ligand.hostToDevice<pos,ind>();
  _openfpm.d_forcefield.hostToDevice<hbtype,radius,hphb,elsc>();
  _openfpm.d_results.resize(params.nposes);
  _openfpm.d_poses.template hostToDevice<0,1,2,3,4,5>();

  size_t global = ceil(params.nposes/(double)_openfpm.posesPerWI);
  global = ceil(global/(double)_openfpm.wgsize);
  size_t local  = _openfpm.wgsize;
  size_t shared = params.ntypes * sizeof(FFParams);

  cudaDeviceSynchronize();

  double start = getTimestamp();

  for(int ii = 0; ii < params.iterations; ++ii)
  {

    CUDA_LAUNCH_DIM3(fasten_main,global, local,
        params.natlig, 
        params.natpro,
        _openfpm.d_protein.toKernel(),
        _openfpm.d_ligand.toKernel(),
        _openfpm.d_poses.toKernel(),
        _openfpm.d_results.toKernel(),
        _openfpm.d_forcefield.toKernel(),
        params.ntypes,
        params.nposes);

  }

  cudaDeviceSynchronize();

  double end = getTimestamp();

  _openfpm.d_results.deviceToHost<0>();

  printTimings(start, end, _openfpm.posesPerWI, _openfpm.gflops_data);
}

#define MAX_PLATFORMS     8
#define MAX_DEVICES      32
#define MAX_INFO_STRING 256

#define DATA_DIR          "bm1"
#define FILE_LIGAND       "/ligand.in"
#define FILE_PROTEIN      "/protein.in"
#define FILE_FORCEFIELD   "/forcefield.in"
#define FILE_POSES        "/poses.in"
#define FILE_REF_ENERGIES "/ref_energies.out"

#define REF_NPOSES 65536

// Energy evaluation parameters
#define CNSTNT   45.0f
#define HBTYPE_F 70
#define HBTYPE_E 69
#define HARDNESS 38.0f
#define NPNPDIST  5.5f
#define NPPDIST   1.0f

void printTimings(double start, double end, double poses_per_wi);
void checkError(int err, const char *op);

FILE* openFile(const char *parent, const char *child,
               const char* mode, long *length)
{
  char name[strlen(parent) + strlen(child) + 1];
  strcpy(name, parent);
  strcat(name, child);

  FILE *file = NULL;
  if (!(file = fopen(name, mode)))
  {
    fprintf(stderr, "Failed to open '%s'\n", name);
    exit(1);
  }
  if(length){
    fseek(file, 0, SEEK_END);
    *length = ftell(file);
    rewind(file);
  }
  return file;
}

int parseInt(const char *str)
{
  char *next;
  int value = strtoul(str, &next, 10);
  return strlen(next) ? -1 : value;
}

void loadParameters(int argc, char *argv[], OpenFPM & _openfpm)
{
  // Defaults
  params.deckDir        = DATA_DIR;
  params.iterations = 8;
  _openfpm.wgsize      = 256;
  _openfpm.posesPerWI  = NUM_TD_PER_THREAD;
  int nposes        = 65536;

  for (int i = 1; i < argc; i++)
  {
    if (!strcmp(argv[i], "--device") || !strcmp(argv[i], "-d"))
    {
      if (++i >= argc || (_openfpm.deviceIndex = parseInt(argv[i])) < 0)
      {
        printf("Invalid device index\n");
        exit(1);
      }
    }
    else if (!strcmp(argv[i], "--iterations") || !strcmp(argv[i], "-i"))
    {
      if (++i >= argc || (params.iterations = parseInt(argv[i])) < 0)
      {
        printf("Invalid number of iterations\n");
        exit(1);
      }
    }
    else if (!strcmp(argv[i], "--numposes") || !strcmp(argv[i], "-n"))
    {
      if (++i >= argc || (nposes = parseInt(argv[i])) < 0)
      {
        printf("Invalid number of poses\n");
        exit(1);
      }
    }
    else if (!strcmp(argv[i], "--posesperwi") || !strcmp(argv[i], "-p"))
    {
      if (++i >= argc || (_openfpm.posesPerWI = parseInt(argv[i])) < 0)
      {
        printf("Invalid poses-per-workitem value\n");
        exit(1);
      }
    }
    else if (!strcmp(argv[i], "--wgsize") || !strcmp(argv[i], "-w"))
    {
      if (++i >= argc || (_openfpm.wgsize = parseInt(argv[i])) < 0)
      {
        printf("Invalid work-group size\n");
        exit(1);
      }
    }
    else if (!strcmp(argv[i], "--deck"))
    {
      if (++i >= argc)
      {
        printf("Invalid deck\n");
        exit(1);
      }
      params.deckDir = argv[i];
    }
    else if (!strcmp(argv[i], "--help") || !strcmp(argv[i], "-h"))
    {
      printf("\n");
      printf("Usage: ./bude [OPTIONS]\n\n");
      printf("Options:\n");
      printf("  -h  --help               Print this message\n");
      printf("      --list               List available devices\n");
      printf("      --device     INDEX   Select device at INDEX\n");
      printf("  -i  --iterations I       Repeat kernel I times\n");
      printf("  -n  --numposes   N       Compute results for N poses\n");
      printf("  -p  --poserperwi PPWI    Compute PPWI poses per work-item\n");
      printf("  -w  --wgsize     WGSIZE  Run with work-group size WGSIZE\n");
      printf("      --deck       DECK    Use the DECK directory as input deck\n");
      printf("\n");
      exit(0);
    }
    else
    {
      printf("Unrecognized argument '%s' (try '--help')\n", argv[i]);
      exit(1);
    }
  }

  FILE *file = NULL;
  long length;

  file = openFile(params.deckDir, FILE_LIGAND, "rb", &length);
  params.natlig = length / sizeof(Atom);
  _openfpm.d_ligand.resize(params.natlig);

  for (int i = 0 ; i < _openfpm.d_ligand.size() ; i++)
  {
	  fread(&_openfpm.d_ligand.template get<pos>(i)[0],sizeof(float),1,file);
	  fread(&_openfpm.d_ligand.template get<pos>(i)[1],sizeof(float),1,file);
	  fread(&_openfpm.d_ligand.template get<pos>(i)[2],sizeof(float),1,file);
	  fread(&_openfpm.d_ligand.template get<ind>(i),sizeof(int),1,file);
  }

  fclose(file);

  file = openFile(params.deckDir, FILE_PROTEIN, "rb", &length);
  params.natpro = length / sizeof(Atom);

  _openfpm.d_protein.resize(params.natpro);

  for (int i = 0 ; i < _openfpm.d_protein.size() ; i++)
  {
          fread(&_openfpm.d_protein.template get<pos>(i)[0],sizeof(float),1,file);
          fread(&_openfpm.d_protein.template get<pos>(i)[1],sizeof(float),1,file);
          fread(&_openfpm.d_protein.template get<pos>(i)[2],sizeof(float),1,file);
          fread(&_openfpm.d_protein.template get<ind>(i),sizeof(int),1,file);
  }

  fclose(file);

  file = openFile(params.deckDir, FILE_FORCEFIELD, "rb", &length);
  params.ntypes = length / sizeof(FFParams);

  _openfpm.d_forcefield.resize(params.ntypes);

  for (int i = 0 ; i < _openfpm.d_forcefield.size() ; i++)
  {
          fread(&_openfpm.d_forcefield.template get<hbtype>(i),sizeof(int),1,file);
          fread(&_openfpm.d_forcefield.template get<radius>(i),sizeof(float),1,file);
          fread(&_openfpm.d_forcefield.template get<hphb>(i),sizeof(float),1,file);
          fread(&_openfpm.d_forcefield.template get<elsc>(i),sizeof(float),1,file);
  }

  fclose(file);

  file = openFile(params.deckDir, FILE_POSES, "rb", &length);
  _openfpm.d_poses.resize(nposes);

  long available = length / 6 / sizeof(float);
  params.nposes = 0;
  while (params.nposes < nposes)
  {
    long fetch = nposes - params.nposes;
    if (fetch > available)
      fetch = available;

      fseek(file, 0*available*sizeof(float), SEEK_SET);
      for (int k = 0 ; k < fetch ; k++)
      {fread(&_openfpm.d_poses.template get<0>(params.nposes+k),sizeof(float),1,file);}

      fseek(file, 1*available*sizeof(float), SEEK_SET);
      for (int k = 0 ; k < fetch ; k++)
      {fread(&_openfpm.d_poses.template get<1>(params.nposes+k),sizeof(float),1,file);}

      fseek(file, 2*available*sizeof(float), SEEK_SET);
      for (int k = 0 ; k < fetch ; k++)
      {fread(&_openfpm.d_poses.template get<2>(params.nposes+k),sizeof(float),1,file);}

      fseek(file, 3*available*sizeof(float), SEEK_SET);
      for (int k = 0 ; k < fetch ; k++)
      {fread(&_openfpm.d_poses.template get<3>(params.nposes+k),sizeof(float),1,file);}

      fseek(file, 4*available*sizeof(float), SEEK_SET);
      for (int k = 0 ; k < fetch ; k++)
      {fread(&_openfpm.d_poses.template get<4>(params.nposes+k),sizeof(float),1,file);}

      fseek(file, 5*available*sizeof(float), SEEK_SET);
      for (int k = 0 ; k < fetch ; k++)
      {fread(&_openfpm.d_poses.template get<5>(params.nposes+k),sizeof(float),1,file);}


    rewind(file);

    params.nposes += fetch;
  }
  fclose(file);
}

#ifndef __APPLE__
#include <fenv.h>
#include <xmmintrin.h>
#include <pmmintrin.h>
#endif

int main(int argc, char *argv[])
{
#ifndef __APPLE__
  _MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_ON);
  _MM_SET_DENORMALS_ZERO_MODE(_MM_DENORMALS_ZERO_ON);
#endif
  init_wrappers();

  OpenFPM _openfpm;
  loadParameters(argc, argv, _openfpm);
  printf("\n");
  printf("Poses     : %d\n", params.nposes);
  printf("Iterations: %d\n", params.iterations);
  printf("Ligands   : %d\n", params.natlig);
  printf("Proteins  : %d\n", params.natpro);
  printf("Deck      : %s\n", params.deckDir);
  float *resultsRef = (float *)malloc(params.nposes*sizeof(float));

  // We run the benchmark 30 times to get mean and variace
  for (int i = 0 ; i < 30 ; i++)
  {
    printf("Iteration %d\n",i);

    runCUDA(_openfpm);
  }

  // calculate mean and variance
  double mean;
  double dev;
  standard_deviation(_openfpm.gflops_data,mean,dev);

  printf("\n\n\nMean %f ~ %f GFlops/s \n\n\n",mean,dev);
  FILE* perf_out = openFile("./","performance_out", "w", NULL);
  char out[256];
  sprintf(out,"%f %f",mean,dev);
  fwrite(out,1,strlen(out),perf_out);
  fclose(perf_out);

  // Load reference results from file
  FILE* ref_energies = openFile(params.deckDir, FILE_REF_ENERGIES, "r", NULL);
  size_t n_ref_poses = params.nposes;
  if (params.nposes > REF_NPOSES) {
    printf("Only validating the first %d poses.\n", REF_NPOSES);
    n_ref_poses = REF_NPOSES;
  }

  for (size_t i = 0; i < n_ref_poses; i++)
    fscanf(ref_energies, "%f", &resultsRef[i]);

  fclose(ref_energies);

  float maxdiff = -100.0f;
  printf("\n Reference        CUDA   (diff)\n");
  for (int i = 0; i < n_ref_poses; i++)
  {
    if (fabs(resultsRef[i]) < 1.f && fabs(_openfpm.d_results.template get<0>(i)) < 1.f) continue;

    float diff = fabs(resultsRef[i] - _openfpm.d_results.template get<0>(i)) / _openfpm.d_results.template get<0>(i);
    if (diff > maxdiff) {
      maxdiff = diff;
      // printf ("Maxdiff: %.2f (%.3f vs %.3f)\n", maxdiff, resultsRef[i], resultsCUDA[i]);
    }

    if (i < 8)
      printf("%7.2f    vs   %7.2f  (%5.2f%%)\n", resultsRef[i], _openfpm.d_results.template get<0>(i), 100*diff);
  }
  printf("\nLargest difference was %.3f%%\n\n", maxdiff*100);

  free(resultsRef);
}

#else

int main(int argc, char *argv[])
{
}

#endif

