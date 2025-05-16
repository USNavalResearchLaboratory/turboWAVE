module;

#include "tw_includes.h"
#include "tw_test.h"

export module discrete_space;
import base;
import pic_primitives;
export import tensor;
export import tasks;

/// This object is an indexing and interpolation scheme for a structured grid.
/// Consider first the *topological indices*, defined as follows.
///
/// Let `ax` be an axis index, numbered from 0 to 3.  Consider a one dimensional
/// strip of cells lined up along axis `ax`.  Suppose this axis uses `L` ghost cell
/// layers.  Then, for any such strip,
/// - Interior cells are labeled `[1,...,dim[ax]]`
/// - Lower ghost cells are labeled `[1-L,...,0]`
/// - Upper ghost cells are labeled `[dim[ax]+1,...,dim[ax]+L]`.
///
/// The topological indices are independent of all storage patterns, In other words, no matter
/// what storage pattern is used, toplogical index `(t0,x0,y0,z0)` will always refer
/// to the same cell.
///
/// Storage patterns are dictated by a *cell encoding*.  This associates the
/// topological indices with a single integer that maps to ascending memory addresses.
/// The default encoding is used to sort particles.
/// Superclasses can add a secondary encoding to resolve higher
/// dimensional storage patterns. For example, `Field` adds another encoding to
/// define how its components are packed.
export struct DiscreteSpace : Testable
{
	protected:

	tw::Int stepNow,stepsToTake;
	/// There is a marker in the inertial frame, this is its velocity.
	/// Ordinary initial value problems will have (1,0,0,0), while "moving windows" will have (1,0,0,1).
	tw::vec4 solutionVelocity;
	/// There is a marker in the inertial frame, this is its position.
	/// The marker will not in general be in perfect synchronism with the window.
	tw::vec4 solutionPosition;
	tw::vec4 altSolutionPosition;
	tw::vec4 windowPosition; ///< the actual position of the computational window
	tw::vec4 altWindowPosition;
	tw::vec4 maxWindowPosition; ///< used to determine when simulation should stop
	tw::vec4 corner; ///< position where all coordinates are minimized on the local domain
	tw::vec4 size; ///< length of the local domain along each axis
	tw::vec4 globalCorner; ///< position where all coordinates are minimized on the global domain
	tw::vec4 globalSize; ///< length of the global domain along each axis
	tw::vec4 spacing; ///< center-to-center cell separation along each axis (uniform grids only)
	tw::vec4 min_spacing; ///< minimum spacing to use for adaptive grids (including time levels)
	tw::vec4 max_spacing; ///< maximum spacing to use for adaptive grids (including time levels)
	tw::vec4 critical_spacing; ///< threshold that triggers some action when we have an adaptive grid
	tw::vec4 freq; ///< inverse of spacing component by component (uniform grids only)
	/// `dim[0..=3]` are the number of cells along each axis.
	/// The notion of field components does not come in at this level.
	/// The space-time domain should be thought of as a window that can shift as a simulation advances.
	/// Often dim[0]=1 and globalCorner[0] advances by spacing[0] each step.
	tw::Int dim[4];
	/// `num` is similar to `dim`, except ghost cells are included.
	tw::Int num[4];
	/// number of interior cells along the given axis, summed over all domains on the low-side of this domain
	tw::Int lowSideCells[4];
	/// number of interior cells along the given axis, summed over all domains
	tw::Int globalCells[4];
	/// Parameters of the default cell encoding.  For the encoding see `EncodeCell`.
	/// The `encodingStride` is 0 along an ignorable axis resulting in a many-one mapping, i.e.,
	/// ghost cells and the one interior cell map to the same cell for the ignorable axis.
	tw::Int encodingStride[4];
	/// The `decodingStride` differs from the `encodingStride` only in that ignorable
	/// axes have unit strides.  See `DecodeCell` for the decoding.
	tw::Int decodingStride[4];
	tw::Int lfg[4]; ///< indices of lower far ghost cells
	tw::Int ufg[4]; ///< indices of upper far ghost cells
	tw::Int lng[4]; ///< indices of lower near ghost cells
	tw::Int ung[4]; ///< indices of upper near ghost cells
	tw::Int layers[4]; ///< number of ghost cell layers, assumed all the same some places
	tw::Int ignorable[4]; ///< 1 if axis is ignorable, 0 otherwise

	public:

	/// Create an empty `DiscreteSpace`
	DiscreteSpace();
	/// Create a `DiscreteSpace` with a single time node and purely *local* coordinates
	DiscreteSpace(tw::Int xDim,tw::Int yDim,tw::Int zDim,const tw::vec4& corner,const tw::vec4& size,tw::Int ghostCellLayers=2);
	/// Change the topology and coordinates:  `gdim`, `gcorner`, and `gsize` describe global window that exists in memory
	/// at any one time, as distinct from the union of all windows that occur during the system evolution.
	void Resize(Task *task,const tw::Int gdim[4],const tw::vec4& gcorner,const tw::vec4& gsize,tw::Int ghostCellLayers=2);
	/// Change the time step.  Use `Simulation::UpdateTimestep` to do this for all modules.
	void SetupTimeInfo(tw::Float dt0) { spacing[0] = dt0; freq[0] = 1.0/dt0; }
	void UpdateWindow(const DiscreteSpace& src) {
		solutionPosition = src.solutionPosition;
		solutionVelocity = src.solutionVelocity;
		windowPosition = src.windowPosition;
		altSolutionPosition = src.altSolutionPosition;
		altWindowPosition = src.altWindowPosition;
		maxWindowPosition = src.maxWindowPosition;
	}
	/// Encode cell with topological indices `(n,i,j,k)`
	tw::Int EncodeCell(tw::Int n,tw::Int i,tw::Int j,tw::Int k) const {
		return (n-lfg[0])*encodingStride[0] + (i-lfg[1])*encodingStride[1] + (j-lfg[2])*encodingStride[2] + (k-lfg[3])*encodingStride[3];
	}
	/// Decode `cell` to produce topological indices `(ijk[0..4])` assuming default encoding.
	void DecodeCell(const tw::Int& cell,tw::Int ijk[4]) const
	{
		ijk[0] = lfg[0] + cell / decodingStride[0];
		ijk[1] = lfg[1] + (cell % decodingStride[0]) / decodingStride[1];
		ijk[2] = lfg[2] + (cell % decodingStride[1]) / decodingStride[2];
		ijk[3] = lfg[3] + (cell % decodingStride[2]) / decodingStride[3];
	}
	/// Decode `cell` to produce topological indices `(ijk[0..4])` assuming default encoding.
	void DecodeCell(const tw::Int& cell,tw::Int *n,tw::Int *i,tw::Int *j,tw::Int *k) const
	{
		*n = lfg[0] + cell / decodingStride[0];
		*i = lfg[1] + (cell % decodingStride[0]) / decodingStride[1];
		*j = lfg[2] + (cell % decodingStride[1]) / decodingStride[2];
		*k = lfg[3] + (cell % decodingStride[2]) / decodingStride[3];
	}
	/// Decode `q` to produce topological indices `(ijk[0..4])` assuming default encoding.
	void DecodeCell(const Primitive& q,tw::Int ijk[4]) const { DecodeCell(q.cell,ijk); }
	/// Decode `q` to produce topological indices `(n,i,j,k)` assuming default encoding.
	void DecodeCell(const Primitive& q,tw::Int *n,tw::Int *i,tw::Int *j,tw::Int *k) const { DecodeCell(q.cell,n,i,j,k); }
	bool RefCellInSpatialDomain(const Primitive& q) const;
	void SetPrimitiveWithPosition(Primitive& q,const tw::vec4& pos) const;
	tw::vec4 PositionFromPrimitive(const Primitive& q) const;
	void MinimizePrimitive(Primitive& q) const;
	void MinimizePrimitive(tw::Int cell[tw::max_bundle_size],tw::Int ijk[4][tw::max_bundle_size],float x[4][tw::max_bundle_size],float domainMask[tw::max_bundle_size]) const;
	tw::Int Dim(const tw::Int& ax) const { return dim[ax]; }
	tw::Int Num(const tw::Int& ax) const { return num[ax]; }
	tw::Int GlobalDim(const tw::Int& ax) const { return globalCells[ax]; }
	tw::Int GlobalCellIndex(const tw::Int& idx,const tw::Int& ax) const { return idx + lowSideCells[ax]; }
	tw::Int LocalCellIndex(const tw::Int& idx,const tw::Int& ax) const { return idx - lowSideCells[ax]; }
	tw::Int Ignorable(const tw::Int& ax) const { return ignorable[ax]; }
	tw::Int Layers(const tw::Int& ax) const { return layers[ax]; }
	tw::Int LNG(const tw::Int& ax) const { return lng[ax]; }
	tw::Int UNG(const tw::Int& ax) const { return ung[ax]; }
	tw::Int LFG(const tw::Int& ax) const { return lfg[ax]; }
	tw::Int UFG(const tw::Int& ax) const { return ufg[ax]; }
	tw::Int StepNow() const { return stepNow; }
	tw::Int SpatialDims();
	tw::Float dx(const tw::Int& ax) const { return spacing[ax]; }
	tw::Float dk(const tw::Int& ax) const { return freq[ax]; }
	tw::vec4 PhysicalSize() { return size; }
	tw::vec4 GlobalPhysicalSize() { return globalSize; }
	tw::vec4 Corner() { return corner; }
	tw::vec4 GlobalCorner() { return globalCorner; }
	tw::Float SolutionPos(const tw::Int& ax) const { return solutionPosition[ax]; }
	tw::Float WindowPos(const tw::Int& ax) const { return windowPosition[ax]; }
	tw::Float MaxWindowPos(const tw::Int& ax) const { return maxWindowPosition[ax]; }
	tw::Float CriticalSpacing(const tw::Int& ax) const { return critical_spacing[ax]; }
	tw::Float MinSpacing(const tw::Int& ax) const { return min_spacing[ax]; }
	tw::Float MaxSpacing(const tw::Int& ax) const { return max_spacing[ax]; }
	bool IsRefCellWithin(const Primitive& q,tw::Int inset) const;
	bool IsPointWithinInterior(const tw::vec4& P);
	bool PowersOfTwo();
	bool TransversePowersOfTwo();
	void GetWeights(weights_3D* weights,const tw::vec4& P);
	void GetWeights(weights_3D* weights,const Primitive& q);
	void GetWeights(float w[3][3][tw::max_bundle_size],float x[4][tw::max_bundle_size]);
	void GetWallWeights(float w[3][3][tw::max_bundle_size],float x[4][tw::max_bundle_size]);

	void ReadCheckpoint(std::ifstream& inFile);
	void WriteCheckpoint(std::ofstream& outFile);

	#ifdef USE_OPENCL
	void CellUpdateProtocol(cl_kernel k,cl_command_queue q);
	void ElementUpdateProtocol(cl_kernel k,cl_command_queue q);
	void LocalUpdateProtocol(cl_kernel k,cl_command_queue q);
	void PointUpdateProtocol(cl_kernel k,cl_command_queue q);
	#endif
};

DiscreteSpace::DiscreteSpace()
{
	ignorable[0] = 0;
	ignorable[1] = 0;
	ignorable[2] = 0;
	ignorable[3] = 0;
	min_spacing = tw::vec4(tw::small_pos);
	max_spacing = tw::vec4(tw::big_pos);
	critical_spacing = tw::vec4(tw::small_pos);
	solutionVelocity = tw::vec4(1.0,0.0,0.0,0.0);
}

DiscreteSpace::DiscreteSpace(tw::Int xDim,tw::Int yDim,tw::Int zDim,const tw::vec4& gcorner,const tw::vec4& gsize,tw::Int ghostCellLayers)
{
	tw::Int domainIndex[4] = { 0,0,0,0 };
	tw::Int domainCount[4] = { 1,1,1,1 };
	dim[0] = 1;
	dim[1] = xDim;
	dim[2] = yDim;
	dim[3] = zDim;

	for (tw::Int i=0;i<4;i++) {
		if (dim[i]==1) {
			layers[i] = 0;
			lfg[i] = 1;
			ufg[i] = 1;
			lng[i] = 1;
			ung[i] = 1;
		} else {
			layers[i] = ghostCellLayers;
			lfg[i] = 1 - layers[i];
			ufg[i] = dim[i] + layers[i];
			lng[i] = 0;
			ung[i] = dim[i] + 1;
		}
		num[i] = ufg[i] - lfg[i] + 1;
	}

	decodingStride[0] = num[1]*num[2]*num[3];
	decodingStride[1] = num[2]*num[3];
	decodingStride[2] = num[3];
	decodingStride[3] = 1;

	for (tw::Int i=0;i<4;i++) {
		encodingStride[i] = (dim[i]==1 ? 0 : decodingStride[i]);
		ignorable[i] = (dim[i]==1 ? 1 : 0);
	}

	globalCorner = gcorner;
	globalSize = gsize;
	for (tw::Int i=0;i<4;i++) {
		globalCells[i] = dim[i];
		lowSideCells[i] = domainIndex[i]*dim[i];
		spacing[i] = globalSize[i]/dim[i];
		freq[i] = 1/spacing[i];
		size[i] = dim[i]*spacing[i];
		corner[i] = gcorner[i] + domainIndex[i]*size[i];
	}
}

void DiscreteSpace::Resize(Task *task,const tw::Int gdim[4],const tw::vec4& gcorner,const tw::vec4& gsize,tw::Int ghostCellLayers)
{
	tw::Int domainIndex[4],domainCount[4];
	task->strip[0].Get_coords(3,domainIndex);
	dim[0] = gdim[0];
	for (auto i = 1; i <= 3; i++) {
		domainCount[i] = task->strip[i].Get_size();
		dim[i] = gdim[i] / domainCount[i];
	}

	for (tw::Int i=0;i<4;i++) {
		if (dim[i]==1) {
			layers[i] = 0;
			lfg[i] = 1;
			ufg[i] = 1;
			lng[i] = 1;
			ung[i] = 1;
		} else {
			layers[i] = ghostCellLayers;
			lfg[i] = 1 - layers[i];
			ufg[i] = dim[i] + layers[i];
			lng[i] = 0;
			ung[i] = dim[i] + 1;
		}
		num[i] = ufg[i] - lfg[i] + 1;
	}

	decodingStride[0] = num[1]*num[2]*num[3];
	decodingStride[1] = num[2]*num[3];
	decodingStride[2] = num[3];
	decodingStride[3] = 1;

	for (tw::Int i=0;i<4;i++) {
		encodingStride[i] = (dim[i]==1 ? 0 : decodingStride[i]);
		ignorable[i] = (dim[i]==1 ? 1 : 0);
	}

	globalCorner = gcorner;
	globalSize = gsize;
	for (tw::Int i=0;i<4;i++) {
		globalCells[i] = gdim[i];
		lowSideCells[i] = domainIndex[i]*dim[i];
		spacing[i] = globalSize[i]/gdim[i];
		freq[i] = 1/spacing[i];
		size[i] = dim[i]*spacing[i];
		corner[i] = gcorner[i] + domainIndex[i]*size[i];
	}
}

inline bool DiscreteSpace::IsRefCellWithin(const Primitive& q,tw::Int inset) const
{
	tw::Int ijk[4];
	DecodeCell(q,ijk);
	bool ans = true;
	for (tw::Int ax=1;ax<=3;ax++)
		ans = ans && ijk[ax]>=lfg[ax]+inset*(1-ignorable[ax]) && ijk[ax]<=ufg[ax]-inset*(1-ignorable[ax]);
	return ans;
}

inline bool DiscreteSpace::IsPointWithinInterior(const tw::vec4& P)
{
	const tw::vec4 PLoc = P - corner;
	return PLoc[1]>=0.0 && PLoc[1]<size[1] && PLoc[2]>=0.0 && PLoc[2]<size[2] && PLoc[3]>=0.0 && PLoc[3]<size[3];
}

inline tw::Int DiscreteSpace::SpatialDims()
{
	return tw::Int(3) - ignorable[1] - ignorable[2] - ignorable[3];
}

inline bool DiscreteSpace::PowersOfTwo()
{
	return (isPowerOfTwo(dim[1]) || dim[1]==1) && (isPowerOfTwo(dim[2]) || dim[2]==1) && (isPowerOfTwo(dim[3]) || dim[3]==1);
}

inline bool DiscreteSpace::TransversePowersOfTwo()
{
	return (isPowerOfTwo(dim[1]) || dim[1]==1) && (isPowerOfTwo(dim[2]) || dim[2]==1);
}

inline bool DiscreteSpace::RefCellInSpatialDomain(const Primitive& q) const
{
	tw::Int ijk[4];
	DecodeCell(q,ijk);
	return	((ijk[1]>=1 && ijk[1]<=dim[1]) || dim[1]==1) &&
			((ijk[2]>=1 && ijk[2]<=dim[2]) || dim[2]==1) &&
			((ijk[3]>=1 && ijk[3]<=dim[3]) || dim[3]==1);
}

inline tw::vec4 DiscreteSpace::PositionFromPrimitive(const Primitive& q) const
{
	tw::Int ijk[4];
	tw::vec4 ans;
	DecodeCell(q,ijk);
	for (tw::Int i=0;i<4;i++)
		ans[i] = corner[i] + spacing[i] * (tw::Float(q.x[i]) + tw::Float(ijk[i]) - tw::Float(0.5));
	return ans;
}

inline void DiscreteSpace::SetPrimitiveWithPosition(Primitive& q,const tw::vec4& P) const
{
	tw::Int ijk[4];
	const tw::vec4 PLoc(P - corner);
	for (tw::Int i=0;i<4;i++)
		ijk[i] = tw::Int(floor(PLoc[i]*freq[i])) + 1;
	q.cell = EncodeCell(ijk[0],ijk[1],ijk[2],ijk[3]);
	for (tw::Int i=0;i<4;i++)
		q.x[i] = PLoc[i]*freq[i] - tw::Float(ijk[i]) + tw::Float(0.5);
}

/// Bundle version of `MinimizePrimitive`.  The bundle version also sets the `domainMask` to 0
/// for particles out of the interior domain, 1 otherwise.  Also the topological indices are loaded for further use.
inline void DiscreteSpace::MinimizePrimitive(tw::Int cell[tw::max_bundle_size],tw::Int ijk[4][tw::max_bundle_size],float x[4][tw::max_bundle_size],float domainMask[tw::max_bundle_size]) const
{
	const float eps = std::numeric_limits<float>::epsilon();
	const auto N = tw::max_bundle_size;
	const auto AB = tw::vec_align_bytes;
	alignas(AB) tw::Int pos[N],neg[N],displ[N],max_displ[N];
	#pragma omp simd aligned(cell,ijk:AB)
	for (auto par=0;par<N;par++)
	{
		DecodeCell(cell[par],&ijk[0][par],&ijk[1][par],&ijk[2][par],&ijk[3][par]);
		domainMask[par] = 1.0f;
	}
	for (auto ax=0;ax<4;ax++)
	{
		#pragma omp simd aligned(ijk,x,domainMask,pos,neg,displ,max_displ:AB)
		for (auto par=0;par<N;par++)
		{
			pos[par] = x[ax][par]>0.0f;
			neg[par] = x[ax][par]<0.0f;
			displ[par] = fabs(x[ax][par]) + 0.5f - neg[par]*eps;
			// max_displ will be 0 if axis is ignorable due to ufg=lfg=ijk=1
			max_displ[par] = pos[par]*(ufg[ax] - ijk[ax][par]) + neg[par]*(ijk[ax][par] - lfg[ax]);
			displ[par] = (displ[par]>max_displ[par])*max_displ[par] + (displ[par]<=max_displ[par])*displ[par] + ignorable[ax]*displ[par];
			displ[par] = (pos[par]-neg[par])*displ[par];
			ijk[ax][par] += displ[par] * (1-ignorable[ax]);
			x[ax][par] -= displ[par];
			domainMask[par] *= float(tw::Int(ijk[ax][par]>=1 && ijk[ax][par]<=dim[ax]));
		}
	}
	#pragma omp simd aligned(cell,ijk:AB)
	for (auto par=0;par<N;par++)
		cell[par] = EncodeCell(ijk[0][par],ijk[1][par],ijk[2][par],ijk[3][par]);
}

/// Change the reference cell to the one containing the particle centroid.
/// This leaves x[ax] within the normal range [-.5,.5) *except* when the particle
/// has left the extended domain.  In the latter case the reference cell is pegged to
/// the current extended domain, leaving x[ax] out of normal range.
inline void DiscreteSpace::MinimizePrimitive(Primitive& q) const
{
	tw::Int ijk[4],displ,max_displ,pos,neg;
	const float eps = std::numeric_limits<float>::epsilon();
	DecodeCell(q.cell,ijk);
	for (auto ax=0;ax<4;ax++)
	{
		const float x = q.x[ax];
		pos = x>0.0f;
		neg = x<0.0f;
		displ = fabs(x) + 0.5f - neg*eps;
		// max_displ will be 0 if axis is ignorable due to ufg=lfg=ijk=1
		max_displ = pos*(ufg[ax] - ijk[ax]) + neg*(ijk[ax] - lfg[ax]);
		displ = (displ>max_displ)*max_displ + (displ<=max_displ)*displ + ignorable[ax]*displ;
		displ = (pos-neg)*displ;
		ijk[ax] += displ * (1-ignorable[ax]);
		q.x[ax] -= displ;
	}
	q.cell = EncodeCell(ijk[0],ijk[1],ijk[2],ijk[3]);
}

inline void DiscreteSpace::GetWeights(weights_3D* weights,const tw::vec4& P)
{
	Primitive q;
	SetPrimitiveWithPosition(q,P);
	GetWeights(weights,q);
}

inline void DiscreteSpace::GetWeights(weights_3D* weights,const Primitive& q)
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

inline void DiscreteSpace::GetWeights(float w[3][3][tw::max_bundle_size],float x[4][tw::max_bundle_size])
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

inline void DiscreteSpace::GetWallWeights(float w[3][3][tw::max_bundle_size],float x[4][tw::max_bundle_size])
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

void DiscreteSpace::ReadCheckpoint(std::ifstream& inFile)
{
	inFile.read((char *)this,sizeof(*this));
}

void DiscreteSpace::WriteCheckpoint(std::ofstream& outFile)
{
	outFile.write((char *)this,sizeof(*this));
}

#ifdef USE_OPENCL

void DiscreteSpace::CellUpdateProtocol(cl_kernel k,cl_command_queue q)
{
	size_t cells = num[1]*num[2]*num[3];
	clEnqueueNDRangeKernel(q,k,1,NULL,&cells,NULL,0,NULL,NULL);
	clFinish(q);
}

void DiscreteSpace::ElementUpdateProtocol(cl_kernel k,cl_command_queue q)
{
	size_t elements = num[0]*num[1]*num[2]*num[3];
	clEnqueueNDRangeKernel(q,k,1,NULL,&elements,NULL,0,NULL,NULL);
	clFinish(q);
}

void DiscreteSpace::LocalUpdateProtocol(cl_kernel k,cl_command_queue q)
{
 	size_t global_offset[3] = { (size_t)layers[1],(size_t)layers[2],(size_t)layers[3] };
	size_t global_work_range[3] = {(size_t)dim[1],(size_t)dim[2],(size_t)dim[3]};
	clEnqueueNDRangeKernel(q,k,3,global_offset,global_work_range,NULL,0,NULL,NULL);
	clFinish(q);
}

void DiscreteSpace::PointUpdateProtocol(cl_kernel k,cl_command_queue q)
{
	size_t global_work_range[3] = { (size_t)num[1],(size_t)num[2],(size_t)num[3] };
	clEnqueueNDRangeKernel(q,k,3,NULL,global_work_range,NULL,0,NULL,NULL);
	clFinish(q);
}

#endif

