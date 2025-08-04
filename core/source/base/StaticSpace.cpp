module;

#include "tw_includes.h"
#include "tw_logger.h"

export module static_space;
import base;
import pic_primitives;
export import tensor;
export import tasks;
import logger;

export const tw::node5 std_packing {4,0,1,2,3};
export const tw::node4 std_coord {1,1,1,1};

/// This object is an indexing and interpolation scheme for a 5D structured grid.
/// The five dimensions are spacetime plus an internal dimension such as vector components.
/// Cells in the space are addressed as (t,x,y,z,c).
///
/// StaticSpace should only contain data that can be assumed to be static throughout a simulation.
/// Spaces with evolving properties are reserved for derivative objects.
///
/// Consider first the *topological indices*, defined as follows.
///
/// Let `ax` be an axis index, numbered from 0 to 3 (leave out 4).  Consider a one dimensional
/// strip of cells lined up along axis `ax`.  Suppose this axis uses `L` ghost cell
/// layers.  Then, for any such strip,
/// - Interior cells are labeled `[1,...,dim[ax]]`
/// - Lower ghost cells are labeled `[1-L,...,0]`
/// - Upper ghost cells are labeled `[dim[ax]+1,...,dim[ax]+L]`.
///
/// The internal dimension labeled 4 is different.  For the internal space ghost cells have no meaning.
/// So this axis has elements labeled `[0,...,dim[4]-1]` in all cases.
///
/// The topological indices are independent of all storage patterns, In other words, no matter
/// what storage pattern is used, toplogical index `(t0,x0,y0,z0,c0)` will always refer
/// to the same cell.
///
/// Storage patterns are dictated by the `packing` array, which in turn provides a cell encoding that
/// maps the topological indices to a memory location.
export struct StaticSpace
{
	protected:

	tw::vec4 spacing; ///< center-to-center cell separation along each axis (uniform grids only)
	tw::vec4 freq; ///< inverse of spacing component by component (uniform grids only)
	/// `dim[0..5]` are the number of cells along each axis.
	/// At this level we have only a static view into an empty spacetime.
	/// We will often have dim[0] = 1, i.e., the view being considered is a snapshot.
	tw::node5 dim;
	/// `num` is similar to `dim`, except ghost cells are included.
	tw::node5 num;
	/// Parameters of the cell encoding.
	/// The `encodingStride` is 0 along an ignorable axis resulting in a many-one mapping, i.e.,
	/// ghost cells and the one interior cell map to the same cell for the ignorable axis.
	tw::node5 encodingStride;
	/// The `decodingStride` differs from the `encodingStride` only in that ignorable
	/// axes have unit strides.
	tw::node5 decodingStride;
	/// Map of axes in order of large stride to small stride, it is not required, but very desirable,
	/// that the internal dimension always be the largest stride.  In that case iterators are compatible
	/// with any Field object that shares the same spatio-temporal structure.
	tw::node5 packing;
	tw::node4 lfg; ///< indices of lower far ghost cells, no internal axis
	tw::node4 ufg; ///< indices of upper far ghost cells, no internal axis
	tw::node4 lng; ///< indices of lower near ghost cells, no internal axis
	tw::node4 ung; ///< indices of upper near ghost cells, no internal axis
	tw::node4 layers; ///< number of ghost cell layers, no internal axis
	tw::node5 ignorable; ///< 1 if axis is ignorable, 0 otherwise

	public:

	/// Create an empty `StaticSpace`
	StaticSpace();
	/// Create a `StaticSpace` with purely *local* coordinates
	StaticSpace(const tw::node5& dim,const tw::vec4& size,const tw::node5& packing,const tw::node4& ghostCellLayers);
	/// Create a `StaticSpace` based on an existing one, but with new internal dimensions
	StaticSpace(const tw::Int& components,const StaticSpace& base);
	/// Resize a `StaticSpace` with the given global parameters
	void Resize(const tw::node5& domains,const tw::node5& gdim,const tw::vec4& gsize,const tw::node5& packing,const tw::node4& ghostCellLayers);
	/// Encode cell with topological indices `(n,i,j,k,0)`
	tw::Int EncodeCell(tw::Int n,tw::Int i,tw::Int j,tw::Int k) const {
		return (n-lfg[0])*encodingStride[0] + (i-lfg[1])*encodingStride[1] + (j-lfg[2])*encodingStride[2] + (k-lfg[3])*encodingStride[3];
	}
	/// Bundle version of `EncodeCell`
	void EncodeCell(tw::Int cell[tw::max_bundle_size],const tw::Int topo[4][tw::max_bundle_size]) const {
		const auto N = tw::max_bundle_size;
		const auto AB = tw::vec_align_bytes;
		#pragma omp simd aligned(cell,topo:AB)
		for (auto par=0;par<N;par++) {
			cell[par] = (topo[0][par]-lfg[0])*encodingStride[0] +
				(topo[1][par]-lfg[1])*encodingStride[1] +
				(topo[2][par]-lfg[2])*encodingStride[2] +
				(topo[3][par]-lfg[3])*encodingStride[3];
		}
	}
	/// Decode `cell` to produce topological indices `(topo[0..4])` assuming topo[4]==0
	void DecodeCell(const tw::Int& cell,tw::Int topo[4]) const
	{
		tw::Int rm = cell;
		for (auto i=0; i<5 ; i++) {
			const auto ax = packing[i];
			if (ax < 4) {
				topo[ax] = lfg[ax] + rm/decodingStride[ax];
			}
			rm %= decodingStride[ax];
		}
	}
	/// Bundle version of `DecodeCell`
	void DecodeCell(const tw::Int cell[tw::max_bundle_size],tw::Int topo[4][tw::max_bundle_size]) const
	{
		const auto N = tw::max_bundle_size;
		const auto AB = tw::vec_align_bytes;
		alignas(AB) tw::Int rm[N];
		#pragma omp simd aligned(cell,topo:AB)
		for (auto par=0;par<N;par++) {
			rm[par] = cell[par];
			for (auto i=0; i<5 ; i++) {
				const auto ax = packing[i];
				if (ax < 4) {
					topo[ax][par] = lfg[ax] + rm[par]/decodingStride[ax];
				}
				rm[par] %= decodingStride[ax];
			}
		}
	}
	/// Decode `q` to produce topological indices `(topo[0..4])`.
	void DecodeCell(const Primitive& q,tw::Int topo[4]) const { DecodeCell(q.cell,topo); }
	bool IsRefCellWithin(const Primitive& q,tw::Int inset) const;
	bool RefCellInSpatialDomain(const Primitive& q) const;
	void MinimizePrimitive(Primitive& q) const;
	void MinimizePrimitive(tw::Int cell[tw::max_bundle_size],tw::Int topo[4][tw::max_bundle_size],float x[4][tw::max_bundle_size],float domainMask[tw::max_bundle_size]) const;
	tw::Int Dim(const tw::Int& ax) const { return dim[ax]; }
	tw::Int Num(const tw::Int& ax) const { return num[ax]; }
	tw::Int Stride(const tw::Int& ax) const { return encodingStride[ax]; }
	tw::Int Ignorable(const tw::Int& ax) const { return ignorable[ax]; }
	tw::Int Layers(const tw::Int& ax) const { return layers[ax]; }
	tw::Int LNG(const tw::Int& ax) const { return lng[ax]; }
	tw::Int UNG(const tw::Int& ax) const { return ung[ax]; }
	tw::Int LFG(const tw::Int& ax) const { return lfg[ax]; }
	tw::Int UFG(const tw::Int& ax) const { return ufg[ax]; }
	tw::Int SpatialDims() { return tw::Int(3) - ignorable[1] - ignorable[2] - ignorable[3]; }
	tw::Float dx(const tw::Int& ax) const { return spacing[ax]; }
	tw::Float dk(const tw::Int& ax) const { return freq[ax]; }
	bool PowersOfTwo() {
		return (isPowerOfTwo(dim[1]) || dim[1]==1) && (isPowerOfTwo(dim[2]) || dim[2]==1) && (isPowerOfTwo(dim[3]) || dim[3]==1);
	}
	bool TransversePowersOfTwo() {
		return (isPowerOfTwo(dim[1]) || dim[1]==1) && (isPowerOfTwo(dim[2]) || dim[2]==1);
	}
	void GetWeights(weights_3D* weights,const Primitive& q);
	void GetWeights(float w[3][3][tw::max_bundle_size],float x[4][tw::max_bundle_size]);
	void GetWallWeights(float w[3][3][tw::max_bundle_size],float x[4][tw::max_bundle_size]);

	void ReadCheckpoint(std::ifstream& inFile) {
		inFile.read((char *)this,sizeof(*this));
	}
	void WriteCheckpoint(std::ofstream& outFile) {
		outFile.write((char *)this,sizeof(*this));
	}

	#ifdef USE_OPENCL
	void CellUpdateProtocol(cl_kernel k,cl_command_queue q);
	void ElementUpdateProtocol(cl_kernel k,cl_command_queue q);
	void LocalUpdateProtocol(cl_kernel k,cl_command_queue q);
	void PointUpdateProtocol(cl_kernel k,cl_command_queue q);
	#endif
};

StaticSpace::StaticSpace()
{
	ignorable[0] = 0;
	ignorable[1] = 0;
	ignorable[2] = 0;
	ignorable[3] = 0;
	ignorable[4] = 0;
}

StaticSpace::StaticSpace(const tw::node5& dim,const tw::vec4& size,const tw::node5& packing,const tw::node4& ghostCellLayers)
{
	const tw::node5 domains { 1, 1, 1, 1, 1 };
	Resize(domains,dim,size,packing,ghostCellLayers);
}

StaticSpace::StaticSpace(const tw::Int& components,const StaticSpace& base) {
	tw::node5 domains {1,1,1,1,1};
	tw::node5 dim {base.dim[0],base.dim[1],base.dim[2],base.dim[3],components};
	tw::vec4 size(base.spacing[0]*dim[0],base.spacing[1]*dim[1],base.spacing[2]*dim[2],base.spacing[3]*dim[3]);
	Resize(domains,dim,size,base.packing,base.layers);
}

void StaticSpace::Resize(const tw::node5& domains,const tw::node5& gdim,const tw::vec4& gsize,const tw::node5& packing,const tw::node4& ghostCellLayers)
{
	for (auto i=0;i<5;i++) {
		dim[i] = gdim[i]/domains[i];
		this->packing[i] = packing[i];
	}

	for (auto i=0;i<4;i++) {
		if (dim[i]==1) {
			layers[i] = 0;
			lfg[i] = 1;
			ufg[i] = 1;
			lng[i] = 1;
			ung[i] = 1;
		} else {
			layers[i] = ghostCellLayers[i];
			lfg[i] = 1 - layers[i];
			ufg[i] = dim[i] + layers[i];
			lng[i] = 1 - tw::Int(layers[i] > 0);
			ung[i] = dim[i] + tw::Int(layers[i] > 0);
		}
		num[i] = ufg[i] - lfg[i] + 1;
	}
	num[4] = dim[4];

	tw::Int stride = 1;
	for (auto i=4;i>=0;i--) {
		const tw::Int ax = packing[i];
		decodingStride[ax] = stride;
		stride *= num[ax];
	}

	for (tw::Int i=0;i<5;i++) {
		encodingStride[i] = (dim[i]==1 ? 0 : decodingStride[i]);
		ignorable[i] = (dim[i]==1 ? 1 : 0);
	}

	for (tw::Int i=0;i<4;i++) {
		spacing[i] = gsize[i]/gdim[i];
		freq[i] = 1/spacing[i];
	}
}

inline bool StaticSpace::IsRefCellWithin(const Primitive& q,tw::Int inset) const
{
	tw::Int topo[4];
	DecodeCell(q,topo);
	bool ans = true;
	for (tw::Int ax=1;ax<=3;ax++)
		ans = ans && topo[ax]>=lfg[ax]+inset*(1-ignorable[ax]) && topo[ax]<=ufg[ax]-inset*(1-ignorable[ax]);
	return ans;
}

inline bool StaticSpace::RefCellInSpatialDomain(const Primitive& q) const
{
	tw::Int topo[4];
	DecodeCell(q,topo);
	return	((topo[1]>=1 && topo[1]<=dim[1]) || dim[1]==1) &&
			((topo[2]>=1 && topo[2]<=dim[2]) || dim[2]==1) &&
			((topo[3]>=1 && topo[3]<=dim[3]) || dim[3]==1);
}

/// Bundle version of `MinimizePrimitive`.  The bundle version also sets the `domainMask` to 0
/// for particles out of the interior domain, 1 otherwise.  Also the topological indices are loaded for further use.
inline void StaticSpace::MinimizePrimitive(tw::Int cell[tw::max_bundle_size],tw::Int topo[4][tw::max_bundle_size],float x[4][tw::max_bundle_size],float domainMask[tw::max_bundle_size]) const
{
	const float eps = std::numeric_limits<float>::epsilon();
	const auto N = tw::max_bundle_size;
	const auto AB = tw::vec_align_bytes;
	alignas(AB) tw::Int pos[N],neg[N],displ[N],max_displ[N];
	DecodeCell(cell,topo);
	#pragma omp simd aligned(cell,topo:AB)
	for (auto par=0;par<N;par++) {
		domainMask[par] = 1.0f;
	}
	for (auto ax=0;ax<4;ax++)
	{
		#pragma omp simd aligned(topo,x,domainMask,pos,neg,displ,max_displ:AB)
		for (auto par=0;par<N;par++)
		{
			pos[par] = x[ax][par]>0.0f;
			neg[par] = x[ax][par]<0.0f;
			displ[par] = fabs(x[ax][par]) + 0.5f - neg[par]*eps;
			// max_displ will be 0 if axis is ignorable due to ufg=lfg=topo=1
			max_displ[par] = pos[par]*(ufg[ax] - topo[ax][par]) + neg[par]*(topo[ax][par] - lfg[ax]);
			displ[par] = (displ[par]>max_displ[par])*max_displ[par] + (displ[par]<=max_displ[par])*displ[par] + ignorable[ax]*displ[par];
			displ[par] = (pos[par]-neg[par])*displ[par];
			topo[ax][par] += displ[par] * (1-ignorable[ax]);
			x[ax][par] -= displ[par];
			domainMask[par] *= float(tw::Int(topo[ax][par]>=1 && topo[ax][par]<=dim[ax]));
		}
	}
	EncodeCell(cell,topo);
}

/// Change the reference cell to the one containing the particle centroid.
/// This leaves x[ax] within the normal range [-.5,.5) *except* when the particle
/// has left the extended domain.  In the latter case the reference cell is pegged to
/// the current extended domain, leaving x[ax] out of normal range.
inline void StaticSpace::MinimizePrimitive(Primitive& q) const
{
	tw::Int topo[4],displ,max_displ,pos,neg;
	const float eps = std::numeric_limits<float>::epsilon();
	DecodeCell(q.cell,topo);
	for (auto ax=0;ax<4;ax++)
	{
		const float x = q.x[ax];
		pos = x>0.0f;
		neg = x<0.0f;
		displ = fabs(x) + 0.5f - neg*eps;
		// max_displ will be 0 if axis is ignorable due to ufg=lfg=topo=1
		max_displ = pos*(ufg[ax] - topo[ax]) + neg*(topo[ax] - lfg[ax]);
		displ = (displ>max_displ)*max_displ + (displ<=max_displ)*displ + ignorable[ax]*displ;
		displ = (pos-neg)*displ;
		topo[ax] += displ * (1-ignorable[ax]);
		q.x[ax] -= displ;
	}
	q.cell = EncodeCell(topo[0],topo[1],topo[2],topo[3]);
}

inline void StaticSpace::GetWeights(weights_3D* weights,const Primitive& q)
{
	weights->cell = q.cell;
	weights->w[0][0] = 0.125f - 0.5f*q.x[1] + 0.5f*q.x[1]*q.x[1];
	weights->w[0][1] = 0.125f - 0.5f*q.x[2] + 0.5f*q.x[2]*q.x[2];
	weights->w[0][2] = 0.125f - 0.5f*q.x[3] + 0.5f*q.x[3]*q.x[3];
	weights->w[1][0] = 0.75f - q.x[1]*q.x[1];
	weights->w[1][1] = 0.75f - q.x[2]*q.x[2];
	weights->w[1][2] = 0.75f - q.x[3]*q.x[3];
	weights->w[2][0] = 0.125f + 0.5f*q.x[1] + 0.5f*q.x[1]*q.x[1];
	weights->w[2][1] = 0.125f + 0.5f*q.x[2] + 0.5f*q.x[2]*q.x[2];
	weights->w[2][2] = 0.125f + 0.5f*q.x[3] + 0.5f*q.x[3]*q.x[3];
}

inline void StaticSpace::GetWeights(float w[3][3][tw::max_bundle_size],float x[4][tw::max_bundle_size])
{
	tw::Int N = tw::max_bundle_size;
	#pragma omp simd aligned(w,x:tw::vec_align_bytes)
	for (tw::Int i=0;i<N;i++)
	{
		w[0][0][i] = 0.125f - 0.5f*x[1][i] + 0.5f*x[1][i]*x[1][i];
		w[0][1][i] = 0.125f - 0.5f*x[2][i] + 0.5f*x[2][i]*x[2][i];
		w[0][2][i] = 0.125f - 0.5f*x[3][i] + 0.5f*x[3][i]*x[3][i];
		w[1][0][i] = 0.75f - x[1][i]*x[1][i];
		w[1][1][i] = 0.75f - x[2][i]*x[2][i];
		w[1][2][i] = 0.75f - x[3][i]*x[3][i];
		w[2][0][i] = 0.125f + 0.5f*x[1][i] + 0.5f*x[1][i]*x[1][i];
		w[2][1][i] = 0.125f + 0.5f*x[2][i] + 0.5f*x[2][i]*x[2][i];
		w[2][2][i] = 0.125f + 0.5f*x[3][i] + 0.5f*x[3][i]*x[3][i];
	}
}

inline void StaticSpace::GetWallWeights(float w[3][3][tw::max_bundle_size],float x[4][tw::max_bundle_size])
{
	tw::Int N = tw::max_bundle_size;
	#pragma omp simd aligned(w,x:tw::vec_align_bytes)
	for (tw::Int i=0;i<N;i++)
	{
		w[0][0][i] = 0.0f;
		w[0][1][i] = 0.0f;
		w[0][2][i] = 0.0f;
		w[1][0][i] = 0.5f - x[1][i];
		w[1][1][i] = 0.5f - x[2][i];
		w[1][2][i] = 0.5f - x[3][i];
		w[2][0][i] = 0.5f + x[1][i];
		w[2][1][i] = 0.5f + x[2][i];
		w[2][2][i] = 0.5f + x[3][i];
	}
}

#ifdef USE_OPENCL

void StaticSpace::CellUpdateProtocol(cl_kernel k,cl_command_queue q)
{
	size_t cells = num[1]*num[2]*num[3];
	clEnqueueNDRangeKernel(q,k,1,NULL,&cells,NULL,0,NULL,NULL);
	clFinish(q);
}

void StaticSpace::ElementUpdateProtocol(cl_kernel k,cl_command_queue q)
{
	size_t elements = num[0]*num[1]*num[2]*num[3];
	clEnqueueNDRangeKernel(q,k,1,NULL,&elements,NULL,0,NULL,NULL);
	clFinish(q);
}

void StaticSpace::LocalUpdateProtocol(cl_kernel k,cl_command_queue q)
{
 	size_t global_offset[3] = { (size_t)layers[1],(size_t)layers[2],(size_t)layers[3] };
	size_t global_work_range[3] = {(size_t)dim[1],(size_t)dim[2],(size_t)dim[3]};
	clEnqueueNDRangeKernel(q,k,3,global_offset,global_work_range,NULL,0,NULL,NULL);
	clFinish(q);
}

void StaticSpace::PointUpdateProtocol(cl_kernel k,cl_command_queue q)
{
	size_t global_work_range[3] = { (size_t)num[1],(size_t)num[2],(size_t)num[3] };
	clEnqueueNDRangeKernel(q,k,3,NULL,global_work_range,NULL,0,NULL,NULL);
	clFinish(q);
}

#endif

