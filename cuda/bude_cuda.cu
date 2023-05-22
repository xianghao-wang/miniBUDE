/**
 * BUDE CUDA kernel file
 **/

#include <float.h>
#include <stdio.h>
#include "shared.h"

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

#define CudaSafeCall( err ) __cudaSafeCall( err, __FILE__, __LINE__ )

inline void __cudaSafeCall( cudaError err, const char *file, const int line )
{
#ifdef CUDA_ERROR_CHECK
    if ( cudaSuccess != err )
    {
        fprintf( stderr, "cudaSafeCall() failed at %s:%i : %s\n",
                 file, line, cudaGetErrorString( err ) );
        exit( -1 );
    }
#endif

    return;
}

// The data structure for one atom - 16 bytes

typedef struct
{
  float x, y, z, w;
} Transform;

#define HARDNESS 38.0f
#define NPNPDIST  5.5f
#define NPPDIST   1.0f

__global__ void fasten_main(const int natlig,
    const int natpro,
    const Atom* protein_molecule,
    const Atom* ligand_molecule,
    const float* transforms_0,
    const float* transforms_1,
    const float* transforms_2,
    const float* transforms_3,
    const float* transforms_4,
    const float* transforms_5,
    float* etotals,
    const FFParams* forcefield,
    const int num_atom_types,
    const int numTransforms);

  extern "C"
void runCUDA(float* results)
{
  printf("\nRunning CUDA\n");

  cudaSetDevice(0);
  CudaSafeCall(cudaMalloc((void**)&_cuda.d_protein, params.natpro*sizeof(Atom)));
  cudaDeviceSynchronize();
  CudaSafeCall(cudaMemcpy(_cuda.d_protein, params.protein, params.natpro*sizeof(Atom), cudaMemcpyHostToDevice));
  cudaDeviceSynchronize();

  CudaSafeCall(cudaMalloc((void**)&_cuda.d_ligand, params.natlig*sizeof(Atom)));
  cudaDeviceSynchronize();
  CudaSafeCall(cudaMemcpy(_cuda.d_ligand, params.ligand, params.natlig*sizeof(Atom), cudaMemcpyHostToDevice));
  cudaDeviceSynchronize();

  CudaSafeCall(cudaMalloc((void**)&_cuda.d_forcefield, params.ntypes*sizeof(FFParams)));
  cudaDeviceSynchronize();
  CudaSafeCall(cudaMemcpy(_cuda.d_forcefield, params.forcefield, params.ntypes*sizeof(FFParams), cudaMemcpyHostToDevice));
  cudaDeviceSynchronize();

  CudaSafeCall(cudaMalloc((void**)&_cuda.d_results, params.nposes*sizeof(float)));
  cudaDeviceSynchronize();

  for(int ii = 0; ii < 6; ++ii)
  {
    CudaSafeCall(cudaMalloc((void**)&_cuda.d_poses[ii], params.nposes*sizeof(float)));
    cudaDeviceSynchronize();
    CudaSafeCall(cudaMemcpy(_cuda.d_poses[ii], params.poses[ii], params.nposes*sizeof(float), cudaMemcpyHostToDevice));
    cudaDeviceSynchronize();
  }

  size_t global = ceil(params.nposes/(double)_cuda.posesPerWI);
  global = ceil(global/(double)_cuda.wgsize);
  size_t local  = _cuda.wgsize;
  size_t shared = params.ntypes * sizeof(FFParams);

  cudaDeviceSynchronize();

  double start = getTimestamp();

  for(int ii = 0; ii < params.iterations; ++ii)
  {
    fasten_main<<<global, local, shared>>>(
        params.natlig, 
        params.natpro,
        _cuda.d_protein,
        _cuda.d_ligand,
        _cuda.d_poses[0],
        _cuda.d_poses[1],
        _cuda.d_poses[2],
        _cuda.d_poses[3],
        _cuda.d_poses[4],
        _cuda.d_poses[5],
        _cuda.d_results,
        _cuda.d_forcefield,
        params.ntypes,
        params.nposes);
    CudaSafeCall(cudaPeekAtLastError());
  }

  cudaDeviceSynchronize();

  double end = getTimestamp();

  CudaSafeCall(cudaMemcpy(results, _cuda.d_results, params.nposes*sizeof(float), cudaMemcpyDeviceToHost));
  cudaDeviceSynchronize();

  printTimings(start, end, _cuda.posesPerWI);
}

__device__ void compute_transformation_matrix(const float transform_0,
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

__global__ void fasten_main(const int natlig,
    const int natpro,
    const Atom* __restrict protein_molecule,
    const Atom* __restrict ligand_molecule,
    const float* __restrict transforms_0,
    const float* __restrict transforms_1,
    const float* __restrict transforms_2,
    const float* __restrict transforms_3,
    const float* __restrict transforms_4,
    const float* __restrict transforms_5,
    float* __restrict etotals,
    const FFParams* global_forcefield,
    const int num_atom_types,
    const int numTransforms)
{
  // Get index of first TD
  int ix = blockIdx.x*blockDim.x*NUM_TD_PER_THREAD + threadIdx.x;

  // Have extra threads do the last member intead of return.
  // A return would disable use of barriers, so not using return is better
  ix = ix < numTransforms ? ix : numTransforms - NUM_TD_PER_THREAD;

#ifdef USE_SHARED
  extern __shared__ FFParams forcefield[];
  if(ix < num_atom_types)
  {
    forcefield[ix] = global_forcefield[ix];
  }
#else
  const FFParams* forcefield = global_forcefield;
#endif

  // Compute transformation matrix to private memory
  float etot[NUM_TD_PER_THREAD];
  Transform transform[NUM_TD_PER_THREAD][3];
  const int lsz = blockDim.x;
  for (int i = 0; i < NUM_TD_PER_THREAD; i++)
  {
    int index = ix + i*lsz;
    compute_transformation_matrix(
        transforms_0[index],
        transforms_1[index],
        transforms_2[index],
        transforms_3[index],
        transforms_4[index],
        transforms_5[index],
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
    const Atom l_atom = ligand_molecule[il];

    const FFParams l_params = forcefield[l_atom.index];
    const bool lhphb_ltz = l_params.hphb<ZERO;
    const bool lhphb_gtz = l_params.hphb>ZERO;

    float3 lpos[NUM_TD_PER_THREAD];
    const float4 linitpos = make_float4(l_atom.x,l_atom.y,l_atom.z,ONE);
    for (int i = 0; i < NUM_TD_PER_THREAD; i++)
    {
      // Transform ligand atom
      lpos[i].x = transform[i][0].w + linitpos.x*transform[i][0].x + 
        linitpos.y*transform[i][0].y + linitpos.z*transform[i][0].z;
      lpos[i].y = transform[i][1].w + linitpos.x*transform[i][1].x + 
        linitpos.y*transform[i][1].y + linitpos.z*transform[i][1].z;
      lpos[i].z = transform[i][2].w + linitpos.x*transform[i][2].x + 
        linitpos.y*transform[i][2].y + linitpos.z*transform[i][2].z;
    }

    // Loop over protein atoms
    int ip = 0;
    do
    {
      // Load protein atom data
      const Atom p_atom = protein_molecule[ip];

      const FFParams p_params = forcefield[p_atom.index];

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
      const float distdslv = (phphb_ltz ? (lhphb_ltz ? NPNPDIST : NPPDIST) : (lhphb_ltz ? NPPDIST : -FLT_MAX));
      const float r_distdslv = 1.0f/distdslv;

      const float chrg_init = l_params.elsc * p_params.elsc;
      const float dslv_init = p_hphb + l_hphb;

      for (int i = 0; i < NUM_TD_PER_THREAD; i++)
      {
        // Calculate distance between atoms
        const float x      = lpos[i].x - p_atom.x;
        const float y      = lpos[i].y - p_atom.y;
        const float z      = lpos[i].z - p_atom.z;
        const float distij = sqrt(x*x + y*y + z*z);

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
        const float coeff  = (ONE - (distbb*r_distdslv));
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
    for (int i = 0; i < NUM_TD_PER_THREAD; i++)
    {
      etotals[td_base+i*blockDim.x] = etot[i]*HALF;
    }
  }
} //end of fasten_main



