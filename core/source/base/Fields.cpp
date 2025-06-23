module;

#include "tw_includes.h"

export module fields;
export import discrete_space;
export import metric_space;
export import tw_iterator;
import numerics;
import fft;

export template <class T>
struct Matrix
{
	std::valarray<T> array;
	tw::Int rows, cols;
	Matrix()
	{
		rows = 3;
		cols = 3;
		array.resize(9);
	}
	Matrix(tw::Int r, tw::Int c)
	{
		rows = r;
		cols = c;
		array.resize(r * c);
	}
	void resize(tw::Int r, tw::Int c)
	{
		rows = r;
		cols = c;
		array.resize(r * c);
	}
	T& operator () (const tw::Int& r, const tw::Int& c)
	{
		return array[r * cols + c];
	}
	T operator () (const tw::Int& r, const tw::Int& c) const
	{
		return array[r * cols + c];
	}
	Matrix<T>& operator = (Matrix<T>& src)
	{
		rows = src.rows;
		cols = src.cols;
		array = src.array;
		return *this;
	}
};

export struct BoundaryCondition
{
	// Consists of a folding operation and a forcing operation.
	// The folding operation is used to clean up after source deposition operations.
	// Typically folding operates on volume integrated states.
	// Typically forcing operates on density states.
	// Folding: strip_i = sum_j(fold_ij * strip_j)
	// Forcing: strip_i = sum_j(force_ij * strip_j) + coeff_i * BC
	// Setting a particular boundary condition consists of defining the matrices.
	// To apply forcing, the caller must pass a value for BC.

	// Handling of inhomogeneous boundary conditions (only enabled for dirichlet):
	// Each boundary condition object has to be fed a value of BC for each strip.
	// Where this comes from is entirely up to the caller of ForcingOperation.
	// The primary caller in turboWAVE is Field::ApplyBoundaryCondition.
	// This can be passed an optional boolean argument ``homogeneous`` (defaults to true)
	// If ``homogeneous==true`` we have BC=0.
	// If ``homogeneous==false`` then BC is assumed to be stored in the far ghost cell.
	// In the latter case part of the forcing matrix is redundant (redundant elements should be 0)

	// WARNING: this object specialized for fields with 2 ghost cell layers

	tw::Float fold[4][4];
	tw::Float force[4][4];
	tw::Float coeff[4];
	tw::Int sgn;

	BoundaryCondition();
	void Reset();
	void Set(tw::bc::fld theBoundaryCondition, tw::grid::side whichSide);
	void FoldingOperation(tw::Float* strip, tw::Int stride)
	{
		// strip should point at the outermost ghost cell involved
		tw::Float temp[4];
		for (tw::Int i = 0; i < 4; i++)
		{
			temp[i] = strip[i * stride * sgn];
			strip[i * stride * sgn] = 0.0;
		}
		for (tw::Int i = 0; i < 4; i++)
			for (tw::Int j = 0; j < 4; j++)
				strip[i * stride * sgn] += fold[i][j] * temp[j];
	}
	void ForcingOperation(tw::Float* strip, tw::Int stride, tw::Float BC)
	{
		// strip should point at the outermost ghost cell involved
		// BC is the value at the wall (dirichletWall) or in ghost cell (dirichletCell), or the difference of adjacent cells (neumannWall)
		tw::Float temp[4];
		for (tw::Int i = 0; i < 4; i++)
		{
			temp[i] = strip[i * stride * sgn];
			strip[i * stride * sgn] = coeff[i] * BC;
		}
		for (tw::Int i = 0; i < 4; i++)
			for (tw::Int j = 0; j < 4; j++)
				strip[i * stride * sgn] += force[i][j] * temp[j];
	}
};

export struct Element
{
	// index vector components from zero
	tw::Int low, high;
	Element();
	Element(tw::Int l, tw::Int h);
	Element(tw::Int i);
	tw::Int Components() const;
	friend Element Union(const Element& e1, const Element& e2);
	Element operator () (const tw::Int& i1, const tw::Int& i2) const
	{
		return Element(low + i1, low + i2);
	}
	Element operator () (const tw::Int& i) const
	{
		return Element(low + i, low + i);
	}
};

/// The Slice contains a copy of a subset of Field data, including its location in the Field.
/// Slice objects can be sent between nodes by passing the return values of `Buffer` and `BufferSize`
/// to various Send and Receive functions  The Field object has functions for moving data
/// to and from a Slice.
export template <class T>
struct Slice
{
	Element e;
	tw::Int lb[4], ub[4]; // bounds of outer surface
private:
	tw::Int decodingStride[4];
	tw::Int encodingStride[4];
	std::valarray<T> data;
public:
	Slice() { ; }
	Slice(const Element& e, tw::Int xl, tw::Int xh, tw::Int yl, tw::Int yh, tw::Int zl, tw::Int zh,bool zero=false);
	Slice(const Element& e, tw::Int low[4], tw::Int high[4], tw::Int ignorable[4],bool zero=false);
	void Resize(const Element& e, tw::Int low[4], tw::Int high[4], tw::Int ignorable[4]);
	T* Buffer()
	{
		return &data[0];
	}
	tw::Int BufferSize()
	{
		return sizeof(T) * e.Components() * (ub[1] - lb[1] + 1) * (ub[2] - lb[2] + 1) * (ub[3] - lb[3] + 1);
	}
	void Translate(const tw::grid::axis& axis, tw::Int displ);
	void Translate(tw::Int x, tw::Int y, tw::Int z);
	T& operator () (const tw::Int& x, const tw::Int& y, const tw::Int& z, const tw::Int& c)
	{
		const tw::Int idx = (c - e.low) + (x - lb[1]) * encodingStride[1] + (y - lb[2]) * encodingStride[2] + (z - lb[3]) * encodingStride[3];
		// if (idx < 0 || idx > data.size()) {
		// 	std::cout << "bad slice access " << idx << "," << data.size() << std::endl;
		// 	std::cout << "coords " << x << "," << y << "," << z << "," << c << std::endl;
		// 	std::cout << "lbounds " << lb[1] << "," << lb[2] << "," << lb[3] << "," << e.low << std::endl;
		// 	std::cout << "ubounds " << ub[1] << "," << ub[2] << "," << ub[3] << "," << e.high << std::endl;
		// 	std::cout << "strides " << encodingStride[1] << "," << encodingStride[2] << "," << encodingStride[3] << std::endl;
		// }
		return data[idx];
	}
	T operator () (const tw::Int& x, const tw::Int& y, const tw::Int& z, const tw::Int& c) const
	{
		const tw::Int idx = (c - e.low) + (x - lb[1]) * encodingStride[1] + (y - lb[2]) * encodingStride[2] + (z - lb[3]) * encodingStride[3];
		// if (idx < 0 || idx > data.size()) {
		// 	std::cout << "bad slice access " << idx << "," << data.size() << std::endl;
		// 	std::cout << "coords " << x << "," << y << "," << z << "," << c - e.low << std::endl;
		// 	std::cout << "lbounds " << lb[1] << "," << lb[2] << "," << lb[3] << std::endl;
		// 	std::cout << "ubounds " << ub[1] << "," << ub[2] << "," << ub[3] << std::endl;
		// 	std::cout << "strides " << encodingStride[1] << "," << encodingStride[2] << "," << encodingStride[3] << std::endl;
		// }
		return data[idx];
	}
	Slice<T>& operator = (T a)
	{
		data = a;
		return *this;
	}
};

/// Field is a DiscreteSpace with data assigned to the cells, and operations on the data.
/// The data is some fixed number of floating point values per cell.
/// The storage pattern is variable, and can be specified by designating a packed axis.
/// The topology (dimensions, ghost cells) is inherited from the DiscreteSpace passed to the constructor.
/// It is possible to mix Fields with varying ghost cell layers, but it is MUCH SAFER
/// to keep the ghost cell layers the same for all Field instances in a calculation.
///
/// Note the Field does not inherit from MetricSpace.  The intention is to have only a single
/// instance of MetricSpace active at any time (mainly to avoid replicating the metric data)
export struct Field: DiscreteSpace
{
protected:

	tw::Int totalCells;
	tw::Int boundaryCells;
	tw::Int ghostCells;
	Matrix<BoundaryCondition> bc0, bc1;
	tw::Int packedAxis, bufferState;
	// Indexing in the following has 0=components, (1,2,3)=spatial axes
	// This is an encoding stride only, some may be zero
	tw::Int stride[4];

	void BoundaryDataToField();
	void FieldToBoundaryData();
	void FieldToGhostData();

public:

	std::valarray<tw::Float> array;
	std::valarray<tw::Float> boundaryData;
	std::valarray<tw::Int> boundaryDataIndexMap;
	std::valarray<tw::Float> ghostData;
	std::valarray<tw::Int> ghostDataIndexMap;
	Task* task;
#ifdef USE_OPENCL
	cl_mem computeBuffer;
	cl_mem boundaryBuffer;
	cl_mem boundaryMapBuffer;
	cl_mem ghostBuffer;
	cl_mem ghostMapBuffer;
#endif

	// SETUP

	Field();
	virtual ~Field();
	void Initialize(tw::Int components, const DiscreteSpace& ds, Task* task, const tw::grid::axis& axis = tw::grid::z);
	void SetBoundaryConditions(const Element& e, const tw::grid::axis& axis, tw::bc::fld low, tw::bc::fld high);
	friend void CopyBoundaryConditions(Field& dst, const Element& dstElement, Field& src, const Element& srcElement);

	// ACCESSORS

	// Conventional accessors using fortran style array indexing
	tw::Float& operator () (const tw::Int& x, const tw::Int& y, const tw::Int& z, const tw::Int& c)
	{
		return array[c * stride[0] + (x - lfg[1]) * stride[1] + (y - lfg[2]) * stride[2] + (z - lfg[3]) * stride[3]];
	}
	tw::Float operator () (const tw::Int& x, const tw::Int& y, const tw::Int& z, const tw::Int& c) const
	{
		return array[c * stride[0] + (x - lfg[1]) * stride[1] + (y - lfg[2]) * stride[2] + (z - lfg[3]) * stride[3]];
	}
	// Accessors for stepping through cells without concern for direction
	tw::Float& operator () (const tw::cell& cell, const tw::Int& c)
	{
		return array[cell.Index(c, stride)];
	}
	tw::Float operator () (const tw::cell& cell, const tw::Int& c) const
	{
		return array[cell.Index(c, stride)];
	}
	// Accessors for iterating across strips
	// Given a strip, returns component c in cell s as measured along the strip
	tw::Float& operator () (const tw::strip& strip, const tw::Int& s, const tw::Int& c)
	{
		return array[strip.Index(s, c, stride)];
	}
	tw::Float operator () (const tw::strip& strip, const tw::Int& s, const tw::Int& c) const
	{
		return array[strip.Index(s, c, stride)];
	}
	// Accessors for promoting compiler vectorization
	tw::Float& operator () (const tw::xstrip<1>& v, const tw::Int& s, const tw::Int& c)
	{
		return array[v.Index1(s, c, stride)];
	}
	tw::Float operator () (const tw::xstrip<1>& v, const tw::Int& s, const tw::Int& c) const
	{
		return array[v.Index1(s, c, stride)];
	}
	tw::Float& operator () (const tw::xstrip<3>& v, const tw::Int& s, const tw::Int& c)
	{
		return array[v.Index3(s, c, stride)];
	}
	tw::Float operator () (const tw::xstrip<3>& v, const tw::Int& s, const tw::Int& c) const
	{
		return array[v.Index3(s, c, stride)];
	}
	// Centered Differencing
	tw::Float operator () (const tw::cell& cell, const tw::Int& c, const tw::Int& ax) const
	{
		const tw::Int idx = cell.Index(c, stride);
		return 0.5 * freq[ax] * (array[idx + stride[ax]] - array[idx - stride[ax]]);
	}
	tw::Float operator () (const tw::xstrip<1>& v, const tw::Int& s, const tw::Int& c, const tw::Int& ax) const
	{
		const tw::Int idx = v.Index1(s, c, stride);
		return 0.5 * freq[ax] * (array[idx + stride[ax]] - array[idx - stride[ax]]);
	}
	tw::Float operator () (const tw::xstrip<3>& v, const tw::Int& s, const tw::Int& c, const tw::Int& ax) const
	{
		const tw::Int idx = v.Index3(s, c, stride);
		return 0.5 * freq[ax] * (array[idx + stride[ax]] - array[idx - stride[ax]]);
	}
	tw::Float d2(const tw::strip& strip, const tw::Int& s, const tw::Int& c, const tw::Int& ax) const
	{
		const tw::Int idx = strip.Index(s, c, stride);
		return freq[ax] * freq[ax] * (array[idx - stride[ax]] - 2.0 * array[idx] + array[idx + stride[ax]]);
	}
	tw::Float d2(const tw::xstrip<1>& v, const tw::Int& s, const tw::Int& c, const tw::Int& ax) const
	{
		const tw::Int idx = v.Index1(s, c, stride);
		return freq[ax] * freq[ax] * (array[idx - stride[ax]] - 2.0 * array[idx] + array[idx + stride[ax]]);
	}
	tw::Float d2(const tw::xstrip<3>& v, const tw::Int& s, const tw::Int& c, const tw::Int& ax) const
	{
		const tw::Int idx = v.Index3(s, c, stride);
		return freq[ax] * freq[ax] * (array[idx - stride[ax]] - 2.0 * array[idx] + array[idx + stride[ax]]);
	}
	// Unit directional offset
	tw::Float fwd(const tw::cell& cell, const tw::Int& c, const tw::Int& ax) const
	{
		return array[cell.Index(c, stride) + stride[ax]];
	}
	tw::Float bak(const tw::cell& cell, const tw::Int& c, const tw::Int& ax) const
	{
		return array[cell.Index(c, stride) - stride[ax]];
	}
	// Forward and Backward Differences and Sums
	tw::Float dfwd(const tw::strip& strip, const tw::Int& s, const tw::Int& c, const tw::Int& ax) const
	{
		tw::Int idx = strip.Index(s, c, stride);
		return freq[ax] * (array[idx + stride[ax]] - array[idx]);
	}
	tw::Float dfwd(const tw::xstrip<3>& v, const tw::Int& s, const tw::Int& c, const tw::Int& ax) const
	{
		tw::Int idx = v.Index3(s, c, stride);
		return freq[ax] * (array[idx + stride[ax]] - array[idx]);
	}
	tw::Float dbak(const tw::xstrip<3>& v, const tw::Int& s, const tw::Int& c, const tw::Int& ax) const
	{
		tw::Int idx = v.Index3(s, c, stride);
		return freq[ax] * (array[idx] - array[idx - stride[ax]]);
	}
	tw::Float sfwd(const tw::xstrip<1>& v, const tw::Int& s, const tw::Int& c, const tw::Int& ax) const
	{
		tw::Int idx = v.Index1(s, c, stride);
		return 0.5 * (array[idx + stride[ax]] + array[idx]);
	}
	tw::Float sbak(const tw::xstrip<1>& v, const tw::Int& s, const tw::Int& c, const tw::Int& ax) const
	{
		tw::Int idx = v.Index1(s, c, stride);
		return 0.5 * (array[idx] + array[idx - stride[ax]]);
	}
	tw::Float sfwd(const tw::xstrip<3>& v, const tw::Int& s, const tw::Int& c, const tw::Int& ax) const
	{
		tw::Int idx = v.Index3(s, c, stride);
		return 0.5 * (array[idx + stride[ax]] + array[idx]);
	}
	tw::Float sbak(const tw::xstrip<3>& v, const tw::Int& s, const tw::Int& c, const tw::Int& ax) const
	{
		tw::Int idx = v.Index3(s, c, stride);
		return 0.5 * (array[idx] + array[idx - stride[ax]]);
	}
	friend Element All(const Field& A) { return Element(0, A.num[0] - 1); }
	tw::Int Components() const { return num[0]; }
	tw::Int TotalCells() const { return totalCells; }
	size_t TotalBytes() const { return array.size() * sizeof(tw::Float); }
	size_t BoundaryBytes() const { return boundaryData.size() * sizeof(tw::Float); }
	size_t GhostBytes() const { return ghostData.size() * sizeof(tw::Float); }
	size_t BoundaryMapBytes() const { return boundaryDataIndexMap.size() * sizeof(tw::Int); }
	size_t GhostMapBytes() const { return ghostDataIndexMap.size() * sizeof(tw::Int); }
	tw::Int Stride(tw::Int ax) const { return stride[ax]; }
	tw::Int* Stride() { return stride; }

	// GPU Support

#ifdef USE_OPENCL
	void InitializeComputeBuffer();
	void SendToComputeBuffer();
	void ReceiveFromComputeBuffer();
	void SendBoundaryCellsToComputeBuffer();
	void ReceiveBoundaryCellsFromComputeBuffer();
	tw::Float CellValueInComputeBuffer(tw::Int i, tw::Int j, tw::Int k, tw::Int c); // an expensive call for one data point (useful for last step of reduction)
	void UpdateGhostCellsInComputeBuffer(const Element& e);
	void SendGhostCellsToComputeBuffer();
	void ZeroGhostCellsInComputeBuffer();
	friend void SwapComputeBuffers(Field& f1, Field& f2);
	friend void CopyComputeBuffer(Field& dst, Field& src);
	friend void CopyComplexComputeBufferMod2(Field& dst, Field& src);
	void FillComputeBufferVec4(const Element& e, tw::vec4& A);
	tw::Float DestructiveSumComputeBuffer();
	tw::Float DestructiveNorm1ComputeBuffer();
	void WeightComputeBufferByVolume(MetricSpace& ms, tw::Float inv);
	void DestructiveComplexMod2ComputeBuffer();
	void MADDComputeBuffer(tw::Float m, tw::Float a);
#endif

	// Transformation

	void MultiplyCellVolume(const MetricSpace& ds);
	void DivideCellVolume(const MetricSpace& ds);
	friend void Swap(Field& f1, Field& f2);
	void Swap(const Element& e1, const Element& e2);
	friend void CopyFieldData(Field& dst, const Element& e_dst, Field& src, const Element& e_src);
	friend void CopyGhostCellData(Field& dst, const Element& e_dst, Field& src, const Element& e_src);
	friend void AddFieldData(Field& dst, const Element& e_dst, Field& src, const Element& e_src);
	friend void AddMulFieldData(Field& dst, const Element& e_dst, Field& src, const Element& e_src, tw::Float mul);
	void SmoothingPass(tw::Int ax, const Element& e, const MetricSpace& ds, const tw::Float& X0, const tw::Float& X1, const tw::Float& X2);
	void Smooth(const Element& e, const MetricSpace& ds, tw::Int smoothPasses[4], tw::Int compPasses[4]);
	void Shift(const Element& e, const tw::strip& s, tw::Int cells, const tw::Float* incoming);
	void Shift(const Element& e, const tw::strip& s, tw::Int cells, const tw::Float& incoming);
	void Hankel(const Element& e, tw::Int modes, std::valarray<tw::Float>& matrix);
	void InverseHankel(const Element& e, tw::Int modes, std::valarray<tw::Float>& matrix);

	// Subsets and Slices

	void GetStrip(std::valarray<tw::Float>& cpy, const tw::strip& s, const tw::Int& c)
	{
		for (tw::Int i = 0; i <= UNG(s.Axis()); i++)
			cpy[i] = (*this)(s, i, c);
	}
	void SetStrip(std::valarray<tw::Float>& cpy, const tw::strip& s, const tw::Int& c)
	{
		for (tw::Int i = 0; i <= UNG(s.Axis()); i++)
			(*this)(s, i, c) = cpy[i];
	}
	tw::vec3 Vec3(const tw::Int& i, const tw::Int& j, const tw::Int& k, const tw::Int& c) const
	{
		return tw::vec3((*this)(i, j, k, c), (*this)(i, j, k, c + 1), (*this)(i, j, k, c + 2));
	}
	tw::vec3 Vec3(const tw::cell& cell, const tw::Int& c) const
	{
		return tw::vec3((*this)(cell, c), (*this)(cell, c + 1), (*this)(cell, c + 2));
	}
	tw::vec3 Vec3(const tw::strip& strip, const tw::Int& s, const tw::Int& c) const
	{
		return tw::vec3((*this)(strip, s, c), (*this)(strip, s, c + 1), (*this)(strip, s, c + 2));
	}
	template <class T>
	void LoadDataIntoSlice(Slice<T>* volume);
	template <class T>
	void SaveDataFromSlice(Slice<T>* volume);
	template <class T>
	void ZeroDataInField(Slice<T>* volume);
	template <class T>
	void AddDataFromSlice(Slice<T>* volume);
	template <class T>
	void AddDataFromSliceAtomic(Slice<T>* volume);

	// Boundaries, messages, smoothing

	void ApplyFoldingCondition(const Element& e);
	void ApplyBoundaryCondition(const Element& e, bool homogeneous = true);
	template <class T>
	void AdjustTridiagonalForBoundaries(const Element& e, const tw::grid::axis& axis, const tw::grid::side& side, std::valarray<T>& T1, std::valarray<T>& T2, std::valarray<T>& T3, std::valarray<T>& source, T val);
	void ZeroGhostCells(const Element& e);

	void StripCopyProtocol(tw::Int axis, tw::Int shift, Slice<tw::Float>* planeIn, Slice<tw::Float>* planeOut, bool add);
	void DownwardCopy(const tw::grid::axis& axis, const Element& e, tw::Int cells);
	void UpwardCopy(const tw::grid::axis& axis, const Element& e, tw::Int cells);
	void DownwardDeposit(const tw::grid::axis& axis, const Element& e, tw::Int cells);
	void UpwardDeposit(const tw::grid::axis& axis, const Element& e, tw::Int cells);

	void CopyFromNeighbors(const Element& e);
	void DepositFromNeighbors(const Element& e);

	Slice<tw::Float>* FormTransposeBlock(const Element& e, const tw::grid::axis& axis1, const tw::grid::axis& axis2, tw::Int start1, tw::Int end1, tw::Int start2, tw::Int end2);
	void Transpose(const Element& e, const tw::grid::axis& axis1, const tw::grid::axis& axis2, Field* target, tw::Int inversion);


	// All-element convenience functions

#ifdef USE_OPENCL
	void UpdateGhostCellsInComputeBuffer()
	{
		UpdateGhostCellsInComputeBuffer(All(*this));
	}
#endif
	void SetBoundaryConditions(const tw::grid::axis& axis, tw::bc::fld low, tw::bc::fld high)
	{
		SetBoundaryConditions(All(*this), axis, low, high);
	}
	void DownwardCopy(const tw::grid::axis& axis, tw::Int cells)
	{
		DownwardCopy(axis, All(*this), cells);
	}
	void UpwardCopy(const tw::grid::axis& axis, tw::Int cells)
	{
		UpwardCopy(axis, All(*this), cells);
	}
	void DownwardDeposit(const tw::grid::axis& axis, tw::Int cells)
	{
		DownwardDeposit(axis, All(*this), cells);
	}
	void UpwardDeposit(const tw::grid::axis& axis, tw::Int cells)
	{
		UpwardDeposit(axis, All(*this), cells);
	}
	void CopyFromNeighbors()
	{
		CopyFromNeighbors(All(*this));
	}
	void DepositFromNeighbors()
	{
		DepositFromNeighbors(All(*this));
	}
	void Transpose(const tw::grid::axis& axis1, const tw::grid::axis& axis2, Field* target, tw::Int inversion)
	{
		Transpose(All(*this), axis1, axis2, target, inversion);
	}
	void ApplyFoldingCondition()
	{
		ApplyFoldingCondition(All(*this));
	}
	void ApplyBoundaryCondition(bool homogeneous = true)
	{
		ApplyBoundaryCondition(All(*this), homogeneous);
	}
	template <class T>
	void AdjustTridiagonalForBoundaries(const tw::grid::axis& axis, const tw::grid::side& side, std::valarray<T>& T1, std::valarray<T>& T2, std::valarray<T>& T3, std::valarray<T>& source, T val)
	{
		AdjustTridiagonalForBoundaries<T>(Element(0), axis, side, T1, T2, T3, source, val);
	}
	void ZeroGhostCells()
	{
		ZeroGhostCells(All(*this));
	}
	void Smooth(const MetricSpace& ds, tw::Int smoothPasses[4], tw::Int compPasses[4])
	{
		Smooth(All(*this), ds, smoothPasses, compPasses);
	}
	void Shift(const tw::strip& s, tw::Int cells, const tw::Float* incoming)
	{
		Shift(All(*this), s, cells, incoming);
	}
	void Shift(const tw::strip& s, tw::Int cells, const tw::Float& incoming)
	{
		Shift(All(*this), s, cells, incoming);
	}
	void Hankel(tw::Int modes, std::valarray<tw::Float>& matrix)
	{
		Hankel(All(*this), modes, matrix);
	}
	void InverseHankel(tw::Int modes, std::valarray<tw::Float>& matrix)
	{
		InverseHankel(All(*this), modes, matrix);
	}

	// File operations

	void ReadCheckpoint(std::ifstream& inFile);
	void WriteCheckpoint(std::ofstream& outFile);

	// Assignment

	Field& operator = (Field& A)
	{
		array = A.array;
		return *this;
	}
	Field& operator += (Field& A)
	{
		array += A.array;
		return *this;
	}
	Field& operator -= (Field& A)
	{
		array -= A.array;
		return *this;
	}
	Field& operator *= (Field& A)
	{
		array *= A.array;
		return *this;
	}
	Field& operator /= (Field& A)
	{
		array /= A.array;
		return *this;
	}
	Field& operator = (tw::Float a)
	{
		array = a;
		return *this;
	}
	Field& operator += (tw::Float a)
	{
		array += a;
		return *this;
	}
	Field& operator -= (tw::Float a)
	{
		array -= a;
		return *this;
	}
	Field& operator *= (tw::Float a)
	{
		array *= a;
		return *this;
	}
	Field& operator /= (tw::Float a)
	{
		array /= a;
		return *this;
	}

	// Grid Interpolation

	void Interpolate(std::valarray<tw::Float>& val, const Element& e, const weights_3D& weights) const;
	void Interpolate(std::valarray<tw::Float>& val, const weights_3D& weights) const
	{
		Interpolate(val, All(*this), weights);
	}
	void InterpolateOnto(std::valarray<tw::Float>& val, const Element& e, const weights_3D& weights);
	void InterpolateOnto(std::valarray<tw::Float>& val, const weights_3D& weights)
	{
		InterpolateOnto(val, All(*this), weights);
	}

	template<tw::Int X, tw::Int Y, tw::Int Z>
	friend void add_const_vec(Field& vf, const tw::vec3& v0);
	template<tw::Int T, tw::Int X, tw::Int Y, tw::Int Z>
	friend void conserved_current_to_dens(Field& current, const MetricSpace& m);
	template <tw::Int C, tw::Int X, tw::Int Y, tw::Int Z>
	friend void assign_grad(const Field& sf, Field& vf, const MetricSpace& m, const tw::Float& scaleFactor);
	template <tw::Int C, tw::Int X, tw::Int Y, tw::Int Z>
	friend void add_grad(const Field& sf, Field& vf, const MetricSpace& m, const tw::Float& scaleFactor);
	template <tw::Int U, tw::Int V, tw::Int W, tw::Int X, tw::Int Y, tw::Int Z>
	friend void add_curlB(const Field& src, Field& dst, const MetricSpace& m, const tw::Float& scaleFactor);
	template <tw::Int U, tw::Int V, tw::Int W, tw::Int X, tw::Int Y, tw::Int Z>
	friend void add_curlE(const Field& src, Field& dst, const MetricSpace& m, const tw::Float& scaleFactor);
	template <tw::Int X, tw::Int Y, tw::Int Z>
	friend tw::Float divE(const Field& vf, tw::Int i, tw::Int j, tw::Int k, const MetricSpace& m);
	template <tw::Int X, tw::Int Y, tw::Int Z>
	friend tw::Float div(const Field& vf, tw::Int i, tw::Int j, tw::Int k, const MetricSpace& m);
	template <tw::Int C, tw::Int X, tw::Int Y, tw::Int Z>
	friend tw::Float div(const Field& coeff, const Field& vf, tw::Int i, tw::Int j, tw::Int k, const MetricSpace& m);
	template <tw::Int C, tw::Int S>
	friend tw::Float div(const Field& coeff, const Field& sf, tw::Int i, tw::Int j, tw::Int k, const MetricSpace& m);
};

template<class T>
void Field::AdjustTridiagonalForBoundaries(const Element& e, const tw::grid::axis& axis, const tw::grid::side& side, std::valarray<T>& T1, std::valarray<T>& T2, std::valarray<T>& T3, std::valarray<T>& source, T val)
{
	// Modify the tridiagonal matrix to respect the boundary condition cell_1 = force_1j * cell_j + coeff_1 * val
	// force_ij and coeff_i are defined by the boundary condition object.
	// Since val may be passed in based on values in the outermost ghost cell, we assume force_10 = 0.
	// N.b. indexing of force_ij and coeff_i starts with outermost ghost cell = 0 and works inward, regardless of tw::grid::side.
	// The element is only used to define the boundary conditions, no field data is used.
	const tw::Int ax = tw::grid::naxis(axis);
	const tw::Int N = Dim(ax) - 1;
	const tw::Int c = e.low; // component determining the BC
	assert(e.low == e.high); // asking for more than one BC is not defined
	assert(T1.size() == Dim(ax));
	if (side == tw::grid::low)
	{
		auto F = bc0(ax, c).force;
		auto C = bc0(ax, c).coeff;
		if (F[1][1] != 1.0) // if F11=1 this is presumed a node boundary
		{
			T2[0] += T1[0] * F[1][2] / (1.0 - F[1][1]);
			T3[0] += T1[0] * F[1][3] / (1.0 - F[1][1]);
			source[0] -= T1[0] * C[1] * val / (1.0 - F[1][1]);
		}
	}
	if (side == tw::grid::high)
	{
		auto F = bc1(ax, c).force;
		auto C = bc1(ax, c).coeff;
		if (F[1][1] != 1.0) // if F11=1 this is presumed a node boundary
		{
			T2[N] += T3[N] * F[1][2] / (1.0 - F[1][1]);
			T1[N] += T3[N] * F[1][3] / (1.0 - F[1][1]);
			source[N] -= T3[N] * C[1] * val / (1.0 - F[1][1]);
		}
	}
}

export template <tw::Int X, tw::Int Y, tw::Int Z>
void add_const_vec(Field& vf, const tw::vec3& v0)
{
#pragma omp parallel
	{
		for (auto v : VectorStripRange<3>(vf, true))
		{
#pragma omp simd
			for (tw::Int k = 0; k <= vf.UNG(3); k++)
				vf(v, k, X) += v0.x;
#pragma omp simd
			for (tw::Int k = 0; k <= vf.UNG(3); k++)
				vf(v, k, Y) += v0.y;
#pragma omp simd
			for (tw::Int k = 0; k <= vf.UNG(3); k++)
				vf(v, k, Z) += v0.z;
		}
	}
}

export template<tw::Int T, tw::Int X, tw::Int Y, tw::Int Z>
void conserved_current_to_dens(Field& current, const MetricSpace& m)
{
#pragma omp parallel
	{
		for (auto cell : CellRange(m, true))
		{
			current(cell, T) /= m.dS(cell, 0);
			current(cell, X) /= m.dS(cell, 1) + tw::small_pos;
			current(cell, Y) /= m.dS(cell, 2) + tw::small_pos;
			current(cell, Z) /= m.dS(cell, 3) + tw::small_pos;
		}
	}
}

export template <tw::Int C, tw::Int X, tw::Int Y, tw::Int Z>
void assign_grad(const Field& sf, Field& vf, const MetricSpace& m, const tw::Float& scaleFactor)
{
	// assign the gradient of a scalar field (sf) to a vector field (vf)
	const tw::Int xN1 = m.UNG(1), yN1 = m.UNG(2), zN1 = m.UNG(3);
#pragma omp parallel for collapse(3) schedule(static)
	for (tw::Int i = 1; i <= xN1; i++)
		for (tw::Int j = 1; j <= yN1; j++)
			for (tw::Int k = 1; k <= zN1; k++)
			{
				vf(i, j, k, X) = scaleFactor * (sf(i, j, k, C) - sf(i - 1, j, k, C)) / m.dl(i, j, k, 1);
				vf(i, j, k, Y) = scaleFactor * (sf(i, j, k, C) - sf(i, j - 1, k, C)) / m.dl(i, j, k, 2);
				vf(i, j, k, Z) = scaleFactor * (sf(i, j, k, C) - sf(i, j, k - 1, C)) / m.dl(i, j, k, 3);
			}
}

export template <tw::Int C, tw::Int X, tw::Int Y, tw::Int Z>
void add_grad(const Field& sf, Field& vf, const MetricSpace& m, const tw::Float& scaleFactor)
{
	// add the gradient of a scalar field (sf) to a vector field (vf)
	const tw::Int xN1 = m.UNG(1), yN1 = m.UNG(2), zN1 = m.UNG(3);
#pragma omp parallel for collapse(3) schedule(static)
	for (tw::Int i = 1; i <= xN1; i++)
		for (tw::Int j = 1; j <= yN1; j++)
			for (tw::Int k = 1; k <= zN1; k++)
			{
				vf(i, j, k, X) += scaleFactor * (sf(i, j, k, C) - sf(i - 1, j, k, C)) / m.dl(i, j, k, 1);
				vf(i, j, k, Y) += scaleFactor * (sf(i, j, k, C) - sf(i, j - 1, k, C)) / m.dl(i, j, k, 2);
				vf(i, j, k, Z) += scaleFactor * (sf(i, j, k, C) - sf(i, j, k - 1, C)) / m.dl(i, j, k, 3);
			}
}

export template <tw::Int U, tw::Int V, tw::Int W, tw::Int X, tw::Int Y, tw::Int Z>
void add_curlB(const Field& src, Field& dst, const MetricSpace& m, const tw::Float& scaleFactor)
{
	// curl of B-field on generalized Yee mesh
	const tw::Int xDim = m.Dim(1), yDim = m.Dim(2), zDim = m.Dim(3);
	const tw::Int xN1 = m.UNG(1), yN1 = m.UNG(2), zN1 = m.UNG(3);
#pragma omp parallel for collapse(2) schedule(static)
	for (tw::Int i = 1; i <= xN1; i++)
		for (tw::Int j = 1; j <= yDim; j++)
		{
			tw::xstrip<3> v(m, i, j, 0);
			tw::xstrip<3> vj(m, i, j + 1, 0);
#pragma omp simd
			for (tw::Int k = 1; k <= zDim; k++)
			{
				dst(v, k, X) += scaleFactor * (src(v, k, V) * m.dlh(v, k, 2) - src(v, k + 1, V) * m.dlh(v, k + 1, 2)) / m.dS(v, k, 1);
				dst(v, k, X) += scaleFactor * (src(vj, k, W) * m.dlh(vj, k, 3) - src(v, k, W) * m.dlh(v, k, 3)) / m.dS(v, k, 1);
			}
		}
#pragma omp parallel for collapse(2) schedule(static)
	for (tw::Int i = 1; i <= xDim; i++)
		for (tw::Int j = 1; j <= yN1; j++)
		{
			tw::xstrip<3> v(m, i, j, 0);
			tw::xstrip<3> vi(m, i + 1, j, 0);
#pragma omp simd
			for (tw::Int k = 1; k <= zDim; k++)
			{
				dst(v, k, Y) -= scaleFactor * (src(v, k, U) * m.dlh(v, k, 1) - src(v, k + 1, U) * m.dlh(v, k + 1, 1)) / m.dS(v, k, 2);
				dst(v, k, Y) -= scaleFactor * (src(vi, k, W) * m.dlh(vi, k, 3) - src(v, k, W) * m.dlh(v, k, 3)) / m.dS(v, k, 2);
			}
		}
#pragma omp parallel for collapse(2) schedule(static)
	for (tw::Int i = 1; i <= xDim; i++)
		for (tw::Int j = 1; j <= yDim; j++)
		{
			tw::xstrip<3> v(m, i, j, 0);
			tw::xstrip<3> vi(m, i + 1, j, 0);
			tw::xstrip<3> vj(m, i, j + 1, 0);
#pragma omp simd
			for (tw::Int k = 1; k <= zN1; k++)
			{
				dst(v, k, Z) += scaleFactor * (src(v, k, U) * m.dlh(v, k, 1) - src(vj, k, U) * m.dlh(vj, k, 1)) / m.dS(v, k, 3);
				dst(v, k, Z) += scaleFactor * (src(vi, k, V) * m.dlh(vi, k, 2) - src(v, k, V) * m.dlh(v, k, 2)) / m.dS(v, k, 3);
			}
		}
}

export template <tw::Int U, tw::Int V, tw::Int W, tw::Int X, tw::Int Y, tw::Int Z>
void add_curlE(const Field& src, Field& dst, const MetricSpace& m, const tw::Float& scaleFactor)
{
	// curl of E-field on generalized Yee mesh
	const tw::Int xN1 = m.UNG(1), yN1 = m.UNG(2), zN1 = m.UNG(3);
#pragma omp parallel for collapse(2) schedule(static)
	for (tw::Int i = 1; i <= xN1; i++)
		for (tw::Int j = 1; j <= yN1; j++)
		{
			tw::xstrip<3> v(m, i, j, 0);
			tw::xstrip<3> vi(m, i - 1, j, 0);
			tw::xstrip<3> vj(m, i, j - 1, 0);
#pragma omp simd
			for (tw::Int k = 1; k <= zN1; k++)
			{
				dst(v, k, X) += scaleFactor * (src(v, k - 1, V) * m.dl(v, k - 1, 2) - src(v, k, V) * m.dl(v, k, 2)) / m.dSh(v, k, 1);
				dst(v, k, X) += scaleFactor * (src(v, k, W) * m.dl(v, k, 3) - src(vj, k, W) * m.dl(vj, k, 3)) / m.dSh(v, k, 1);
				dst(v, k, Y) -= scaleFactor * (src(v, k - 1, U) * m.dl(v, k - 1, 1) - src(v, k, U) * m.dl(v, k, 1)) / m.dSh(v, k, 2);
				dst(v, k, Y) -= scaleFactor * (src(v, k, W) * m.dl(v, k, 3) - src(vi, k, W) * m.dl(vi, k, 3)) / m.dSh(v, k, 2);
				dst(v, k, Z) += scaleFactor * (src(vj, k, U) * m.dl(vj, k, 1) - src(v, k, U) * m.dl(v, k, 1)) / m.dSh(v, k, 3);
				dst(v, k, Z) += scaleFactor * (src(v, k, V) * m.dl(v, k, 2) - src(vi, k, V) * m.dl(vi, k, 2)) / m.dSh(v, k, 3);
			}
		}
}

export template <tw::Int X, tw::Int Y, tw::Int Z>
tw::Float divE(const Field& vf, tw::Int i, tw::Int j, tw::Int k, const MetricSpace& m)
{
	// divergence of E-field on generalized Yee mesh
	tw::Float ans, vol;
	vol = m.dS(i, j, k, 0);
	ans = vf(i + 1, j, k, X) * m.dS(i + 1, j, k, 1);
	ans -= vf(i, j, k, X) * m.dS(i, j, k, 1);
	ans += vf(i, j + 1, k, Y) * m.dS(i, j + 1, k, 2);
	ans -= vf(i, j, k, Y) * m.dS(i, j, k, 2);
	ans += vf(i, j, k + 1, Z) * m.dS(i, j, k + 1, 3);
	ans -= vf(i, j, k, Z) * m.dS(i, j, k, 3);
	return ans / vol;
}

template <tw::Int X, tw::Int Y, tw::Int Z>
inline tw::Float div(const Field& vf, tw::Int i, tw::Int j, tw::Int k, const MetricSpace& m)
{
	// divergence of a vector field
	tw::Float ans, vol;
	vol = m.dS(i, j, k, 0);
	ans = (vf(i + 1, j, k, X) + vf(i, j, k, X)) * m.dS(i + 1, j, k, 1);
	ans -= (vf(i - 1, j, k, X) + vf(i, j, k, X)) * m.dS(i, j, k, 1);
	ans += (vf(i, j + 1, k, Y) + vf(i, j, k, Y)) * m.dS(i, j + 1, k, 2);
	ans -= (vf(i, j - 1, k, Y) + vf(i, j, k, Y)) * m.dS(i, j, k, 2);
	ans += (vf(i, j, k + 1, Z) + vf(i, j, k, Z)) * m.dS(i, j, k + 1, 3);
	ans -= (vf(i, j, k - 1, Z) + vf(i, j, k, Z)) * m.dS(i, j, k, 3);
	return 0.5 * ans / vol;
}

template <tw::Int C, tw::Int X, tw::Int Y, tw::Int Z>
inline tw::Float div(const Field& coeff, const Field& vf, tw::Int i, tw::Int j, tw::Int k, const MetricSpace& m)
{
	// divergence of a scalar field times a vector field
	tw::Float ans, vol;
	vol = m.dS(i, j, k, 0);
	ans = (coeff(i + 1, j, k, C) + coeff(i, j, k, C)) * (vf(i + 1, j, k, X) + vf(i, j, k, X)) * m.dS(i + 1, j, k, 1);
	ans -= (coeff(i - 1, j, k, C) + coeff(i, j, k, C)) * (vf(i - 1, j, k, X) + vf(i, j, k, X)) * m.dS(i, j, k, 1);
	ans += (coeff(i, j + 1, k, C) + coeff(i, j, k, C)) * (vf(i, j + 1, k, Y) + vf(i, j, k, Y)) * m.dS(i, j + 1, k, 2);
	ans -= (coeff(i, j - 1, k, C) + coeff(i, j, k, C)) * (vf(i, j - 1, k, Y) + vf(i, j, k, Y)) * m.dS(i, j, k, 2);
	ans += (coeff(i, j, k + 1, C) + coeff(i, j, k, C)) * (vf(i, j, k + 1, Z) + vf(i, j, k, Z)) * m.dS(i, j, k + 1, 3);
	ans -= (coeff(i, j, k - 1, C) + coeff(i, j, k, C)) * (vf(i, j, k - 1, Z) + vf(i, j, k, Z)) * m.dS(i, j, k, 3);
	return 0.25 * ans / vol;
}

template <tw::Int C, tw::Int S>
inline tw::Float div(const Field& coeff, const Field& sf, tw::Int i, tw::Int j, tw::Int k, const MetricSpace& m)
{
	// divergence of a scalar field times the gradient of another scalar field
	tw::Float ans, vol;
	vol = m.dS(i, j, k, 0);
	ans = (coeff(i + 1, j, k, C) + coeff(i, j, k, C)) * (sf(i + 1, j, k, S) - sf(i, j, k, S)) * m.dS(i + 1, j, k, 1) / m.dl(i + 1, j, k, 1);
	ans -= (coeff(i - 1, j, k, C) + coeff(i, j, k, C)) * (sf(i, j, k, S) - sf(i - 1, j, k, S)) * m.dS(i, j, k, 1) / m.dl(i, j, k, 1);
	ans += (coeff(i, j + 1, k, C) + coeff(i, j, k, C)) * (sf(i, j + 1, k, S) - sf(i, j, k, S)) * m.dS(i, j + 1, k, 2) / m.dl(i, j + 1, k, 2);
	ans -= (coeff(i, j - 1, k, C) + coeff(i, j, k, C)) * (sf(i, j, k, S) - sf(i, j - 1, k, S)) * m.dS(i, j, k, 2) / m.dl(i, j, k, 2);
	ans += (coeff(i, j, k + 1, C) + coeff(i, j, k, C)) * (sf(i, j, k + 1, S) - sf(i, j, k, S)) * m.dS(i, j, k + 1, 3) / m.dl(i, j, k + 1, 3);
	ans -= (coeff(i, j, k - 1, C) + coeff(i, j, k, C)) * (sf(i, j, k, S) - sf(i, j, k - 1, S)) * m.dS(i, j, k, 3) / m.dl(i, j, k, 3);
	return 0.5 * ans / vol;
}

////////////////////////
//                    //
//  GATHER / SCATTER  //
//                    //
////////////////////////


inline void Field::Interpolate(std::valarray<tw::Float>& val, const Element& e, const weights_3D& weights) const
{
	tw::Int i, j, k, s, ijk[4];
	tw::Float wijk;
	DecodeCell(weights.cell, ijk);
	val = 0.0;

	for (s = e.low; s <= e.high; s++)
		for (i = 0; i < 3; i++)
			for (j = 0; j < 3; j++)
				for (k = 0; k < 3; k++)
				{
					wijk = weights.w[i][0] * weights.w[j][1] * weights.w[k][2];
					val[s] += wijk * (*this)(ijk[1] + i - 1, ijk[2] + j - 1, ijk[3] + k - 1, s);
				}
}

inline void Field::InterpolateOnto(std::valarray<tw::Float>& val, const Element& e, const weights_3D& weights)
{
	tw::Int i, j, k, s, ijk[4];
	tw::Float wijk;
	DecodeCell(weights.cell, ijk);

	for (s = e.low; s <= e.high; s++)
		for (i = 0; i < 3; i++)
			for (j = 0; j < 3; j++)
				for (k = 0; k < 3; k++)
				{
					wijk = weights.w[i][0] * weights.w[j][1] * weights.w[k][2];
#pragma omp atomic update
					(*this)(ijk[1] + i - 1, ijk[2] + j - 1, ijk[3] + k - 1, s) += wijk * val[s];
				}
}



//////////////////////////////
//                          //
//    AUTO FIELD TEMPLATE   //
//                          //
//////////////////////////////
//
// AutoField class may be used in the PARTICULAR case where the storage scheme
// is such that the components in a cell are stored contiguously in memory.
// The AutoField can then be used to treat the field as an array of type X,
// where X is some multi-component type that supports basic arithmetic operations.
// This is particularly useful for complex fields.

export template <class T>
struct AutoField : Field
{
	AutoField()
	{
		packedAxis = 0;
		num[0] = sizeof(T) / sizeof(tw::Float);
	}
	void Initialize(const DiscreteSpace& ds, Task* task)
	{
		Field::Initialize(num[0], ds, task, tw::grid::t);
	}

	// Interface to superclass accessors

	tw::Float& operator () (const tw::Int& x, const tw::Int& y, const tw::Int& z, const tw::Int& c)
	{
		return Field::operator () (x, y, z, c);
	}
	tw::Float& operator () (const tw::cell& cell, const tw::Int& c)
	{
		return Field::operator () (cell, c);
	}
	tw::Float& operator () (const tw::strip& strip, const tw::Int& x, const tw::Int& c)
	{
		return Field::operator () (strip, x, c);
	}
	tw::Float& operator () (const tw::xstrip<1>& v, const tw::Int& i, const tw::Int& c)
	{
		return Field::operator () (v, i, c);
	}
	tw::Float& operator () (const tw::xstrip<3>& v, const tw::Int& k, const tw::Int& c)
	{
		return Field::operator () (v, k, c);
	}
	tw::Float operator () (const tw::Int& x, const tw::Int& y, const tw::Int& z, const tw::Int& c) const
	{
		return Field::operator () (x, y, z, c);
	}
	tw::Float operator () (const tw::cell& cell, const tw::Int& c) const
	{
		return Field::operator () (cell, c);
	}
	tw::Float operator () (const tw::strip& strip, const tw::Int& x, const tw::Int& c) const
	{
		return Field::operator () (strip, x, c);
	}
	tw::Float operator () (const tw::xstrip<1>& v, const tw::Int& i, const tw::Int& c) const
	{
		return Field::operator () (v, i, c);
	}
	tw::Float operator () (const tw::xstrip<3>& v, const tw::Int& k, const tw::Int& c) const
	{
		return Field::operator () (v, k, c);
	}
	tw::Float operator () (const tw::cell& cell, const tw::Int& c, const tw::Int& ax) const
	{
		return Field::operator () (cell, c, ax);
	}
	tw::Float operator () (const tw::xstrip<1>& v, const tw::Int& i, const tw::Int& c, const tw::Int& ax) const
	{
		return Field::operator () (v, i, c, ax);
	}
	tw::Float operator () (const tw::xstrip<3>& v, const tw::Int& k, const tw::Int& c, const tw::Int& ax) const
	{
		return Field::operator () (v, k, c, ax);
	}

	// Auto versions are distinguished by number of arguments.
	// f(i,j,k) or f(cell) or f(strip,i) call superclass and cast result to type T

	T& operator () (const tw::Int& x, const tw::Int& y, const tw::Int& z)
	{
		return (T&)Field::operator () (x, y, z, 0);
	}
	T& operator () (const tw::cell& cell)
	{
		return (T&)Field::operator () (cell, 0);
	}
	T& operator () (const tw::strip& strip, const tw::Int& x)
	{
		return (T&)Field::operator () (strip, x, 0);
	}
	T& operator () (const tw::xstrip<1>& v, const tw::Int& i)
	{
		return (T&)Field::operator () (v, i, 0);
	}
	T& operator () (const tw::xstrip<3>& v, const tw::Int& k)
	{
		return (T&)Field::operator () (v, k, 0);
	}

	void Interpolate(T* val, const weights_3D& weights) const
	{
		std::valarray<tw::Float> temp(num[0]);
		Field::Interpolate(temp, weights);
		for (tw::Int i = 0; i < num[0]; i++)
			((tw::Float*)val)[i] = temp[i];
	}
	void InterpolateOnto(const T& val, const weights_3D& weights)
	{
		std::valarray<tw::Float> temp(num[0]);
		for (tw::Int i = 0; i < num[0]; i++)
			temp[i] = ((tw::Float*)&val)[i];
		Field::InterpolateOnto(temp, weights);
	}
	void Shift(const tw::strip& s, tw::Int cells, const T& incoming)
	{
		Field::Shift(s, cells, (tw::Float*)&incoming);
	}

	// Assignment

	AutoField<T>& operator = (T a)
	{
		for (auto cell : EntireCellRange(*this))
			(*this)(cell) = a;
		return *this;
	}
	AutoField<T>& operator += (T a)
	{
		for (auto cell : EntireCellRange(*this))
			(*this)(cell) += a;
		return *this;
	}
	AutoField<T>& operator -= (T a)
	{
		for (auto cell : EntireCellRange(*this))
			(*this)(cell) -= a;
		return *this;
	}
	AutoField<T>& operator *= (T a)
	{
		for (auto cell : EntireCellRange(*this))
			(*this)(cell) *= a;
		return *this;
	}
};



export struct ScalarField :AutoField<tw::Float>
{
	tw::Float AxialEigenvalue(tw::Int z);
	tw::Float Eigenvalue(tw::Int x, tw::Int y);
	tw::Float CyclicEigenvalue(tw::Int x, tw::Int y);
	tw::Float CyclicEigenvalue(tw::Int x, tw::Int y, tw::Int z);

	void AxialSineTransform();
	void InverseAxialSineTransform();
	void TransverseSineTransform();
	void InverseTransverseSineTransform();
	void TransverseCosineTransform();
	void InverseTransverseCosineTransform();
	void TransverseFFT();
	void InverseTransverseFFT();

	// Assignment

	ScalarField& operator = (AutoField<tw::Float>& A)
	{
		return (ScalarField&)Field::operator=(A);
	}
	ScalarField& operator += (AutoField<tw::Float>& A)
	{
		return (ScalarField&)Field::operator+=(A);
	}
	ScalarField& operator -= (AutoField<tw::Float>& A)
	{
		return (ScalarField&)Field::operator-=(A);
	}
	ScalarField& operator *= (AutoField<tw::Float>& A)
	{
		return (ScalarField&)Field::operator*=(A);
	}
	ScalarField& operator = (tw::Float a)
	{
		return (ScalarField&)Field::operator=(a);
	}
	ScalarField& operator += (tw::Float a)
	{
		return (ScalarField&)Field::operator+=(a);
	}
	ScalarField& operator -= (tw::Float a)
	{
		return (ScalarField&)Field::operator-=(a);
	}
	ScalarField& operator *= (tw::Float a)
	{
		return (ScalarField&)Field::operator*=(a);
	}
};

export struct ComplexField :AutoField<tw::Complex>
{
	tw::Float CyclicEigenvalue(tw::Int x, tw::Int y);
	tw::Float CyclicEigenvalue(tw::Int x, tw::Int y, tw::Int z);
	void FFT();
	void InverseFFT();
	void TransverseFFT();
	void InverseTransverseFFT();

	// Assignment

	ComplexField& operator = (AutoField<tw::Complex>& A)
	{
		return (ComplexField&)Field::operator=(A);
	}
	ComplexField& operator += (AutoField<tw::Complex>& A)
	{
		return (ComplexField&)Field::operator+=(A);
	}
	ComplexField& operator -= (AutoField<tw::Complex>& A)
	{
		return (ComplexField&)Field::operator-=(A);
	}
	ComplexField& operator *= (AutoField<tw::Complex>& A)
	{
		return (ComplexField&)Field::operator*=(A);
	}
	ComplexField& operator = (tw::Complex a)
	{
		return (ComplexField&)AutoField<tw::Complex>::operator=(a);
	}
	ComplexField& operator += (tw::Complex a)
	{
		return (ComplexField&)AutoField<tw::Complex>::operator+=(a);
	}
	ComplexField& operator -= (tw::Complex a)
	{
		return (ComplexField&)AutoField<tw::Complex>::operator-=(a);
	}
	ComplexField& operator *= (tw::Complex a)
	{
		return (ComplexField&)AutoField<tw::Complex>::operator*=(a);
	}
};

export struct Vec3Field :AutoField<tw::vec3>
{
	//tw::Float CyclicEigenvalue(tw::Int x,tw::Int y);
	//void TransverseFFT();
	//void InverseTransverseFFT();

	// Assignment

	Vec3Field& operator = (AutoField<tw::vec3>& A)
	{
		return (Vec3Field&)Field::operator=(A);
	}
	Vec3Field& operator += (AutoField<tw::vec3>& A)
	{
		return (Vec3Field&)Field::operator+=(A);
	}
	Vec3Field& operator -= (AutoField<tw::vec3>& A)
	{
		return (Vec3Field&)Field::operator-=(A);
	}
	Vec3Field& operator *= (AutoField<tw::vec3>& A)
	{
		return (Vec3Field&)Field::operator*=(A);
	}
	Vec3Field& operator = (tw::vec3 a)
	{
		return (Vec3Field&)AutoField<tw::vec3>::operator=(a);
	}
	Vec3Field& operator += (tw::vec3 a)
	{
		return (Vec3Field&)AutoField<tw::vec3>::operator+=(a);
	}
	Vec3Field& operator -= (tw::vec3 a)
	{
		return (Vec3Field&)AutoField<tw::vec3>::operator-=(a);
	}
	Vec3Field& operator *= (tw::vec3 a)
	{
		return (Vec3Field&)AutoField<tw::vec3>::operator*=(a);
	}
};

////////////////////
//                //
// Slice Template //
//                //
////////////////////

template <class T>
Slice<T>::Slice(const Element& e, tw::Int xl, tw::Int xh, tw::Int yl, tw::Int yh, tw::Int zl, tw::Int zh,bool zero)
{
	tw::Int low[4] = { 0 , xl , yl , zl };
	tw::Int high[4] = { 0 , xh , yh , zh };
	tw::Int ignorable[4] = { 0 , 0 , 0 , 0 };
	Resize(e, low, high, ignorable);
	if (zero)
		data = 0.0;
}

template <class T>
Slice<T>::Slice(const Element& e, tw::Int low[4], tw::Int high[4], tw::Int ignorable[4],bool zero)
{
	Resize(e, low, high, ignorable);
	if (zero)
		data = 0.0;
}

template <class T>
void Slice<T>::Resize(const Element& e, tw::Int low[4], tw::Int high[4], tw::Int ignorable[4])
{
	this->e = e;
	lb[0] = low[0];
	ub[0] = high[0];
	lb[1] = low[1];
	ub[1] = high[1];
	lb[2] = low[2];
	ub[2] = high[2];
	lb[3] = low[3];
	ub[3] = high[3];
	// Decoding strides, no zero strides
	decodingStride[3] = e.Components();
	decodingStride[2] = decodingStride[3] * (ub[3] - lb[3] + 1);
	decodingStride[1] = decodingStride[2] * (ub[2] - lb[2] + 1);
	// Encoding strides, zero strides used
	encodingStride[1] = decodingStride[1] * (1 - ignorable[1]);
	encodingStride[2] = decodingStride[2] * (1 - ignorable[2]);
	encodingStride[3] = decodingStride[3] * (1 - ignorable[3]);
	// Allocate space
	data.resize(e.Components() * (ub[1] - lb[1] + 1) * (ub[2] - lb[2] + 1) * (ub[3] - lb[3] + 1));
}

template <class T>
void Slice<T>::Translate(const tw::grid::axis& axis, tw::Int displ)
{
	tw::Int ax = tw::grid::naxis(axis);
	lb[ax] += displ;
	ub[ax] += displ;
}

template <class T>
void Slice<T>::Translate(tw::Int x, tw::Int y, tw::Int z)
{
	lb[1] += x;
	ub[1] += x;
	lb[2] += y;
	ub[2] += y;
	lb[3] += z;
	ub[3] += z;
}

template <class T>
void Field::LoadDataIntoSlice(Slice<T>* v)
{
	tw::Int i, j, k, s;
	for (s = v->e.low; s <= v->e.high; s++)
		for (i = v->lb[1]; i <= v->ub[1]; i++)
			for (j = v->lb[2]; j <= v->ub[2]; j++)
				for (k = v->lb[3]; k <= v->ub[3]; k++)
					(*v)(i, j, k, s) = (*this)(i, j, k, s);
}

template <class T>
void Field::SaveDataFromSlice(Slice<T>* v)
{
	tw::Int i, j, k, s;
	for (s = v->e.low; s <= v->e.high; s++)
		for (i = v->lb[1]; i <= v->ub[1]; i++)
			for (j = v->lb[2]; j <= v->ub[2]; j++)
				for (k = v->lb[3]; k <= v->ub[3]; k++)
					(*this)(i, j, k, s) = (*v)(i, j, k, s);
}

template <class T>
void Field::ZeroDataInField(Slice<T>* v)
{
	tw::Int i, j, k, s;
	for (s = v->e.low; s <= v->e.high; s++)
		for (i = v->lb[1]; i <= v->ub[1]; i++)
			for (j = v->lb[2]; j <= v->ub[2]; j++)
				for (k = v->lb[3]; k <= v->ub[3]; k++)
					(*this)(i, j, k, s) = 0.0;
}

template <class T>
void Field::AddDataFromSlice(Slice<T>* v)
{
	tw::Int i, j, k, s;
	for (s = v->e.low; s <= v->e.high; s++)
		for (i = v->lb[1]; i <= v->ub[1]; i++)
			for (j = v->lb[2]; j <= v->ub[2]; j++)
				for (k = v->lb[3]; k <= v->ub[3]; k++)
					(*this)(i, j, k, s) += (*v)(i, j, k, s);
}

template <class T>
void Field::AddDataFromSliceAtomic(Slice<T>* v)
{
	tw::Int i, j, k, s;
	for (s = v->e.low; s <= v->e.high; s++)
		for (i = v->lb[1]; i <= v->ub[1]; i++)
			for (j = v->lb[2]; j <= v->ub[2]; j++)
				for (k = v->lb[3]; k <= v->ub[3]; k++)
#pragma omp atomic update
				(*this)(i, j, k, s) += (*v)(i, j, k, s);
}

/////////////////////////
//                     //
// BOUNDARY CONDITIONS //
//                     //
/////////////////////////


BoundaryCondition::BoundaryCondition()
{
	// default to none (identity matrices)
	Reset();
	// default to low side
	sgn = 1;
}

void BoundaryCondition::Reset()
{
	// Set all variables to the identity operation
	for (tw::Int i=0;i<4;i++)
	{
		for (tw::Int j=0;j<4;j++)
		{
			fold[i][j] = 0.0;
			force[i][j] = 0.0;
		}
		coeff[i] = 0.0;
		fold[i][i] = 1.0;
		force[i][i] = 1.0;
	}
}

void BoundaryCondition::Set(tw::bc::fld theBoundaryCondition,tw::grid::side whichSide)
{
	// In this routine, indexing is offset from DiscreteSpace indexing convention.
	// Namely, 0 is the outer ghost cell layer, 1 is the inner, 2 is the edge cell, etc.
	// When BC is applied, negative strides are used to produce mirror image indexing on high side

	// Inhomogeneous boundary conditions are only enabled for dirichlet variants.
	// If we enabled it for neumann it would lead to unfriendly results in many cases.

	Reset();
	switch (whichSide)
	{
		case tw::grid::low:
			sgn = 1;
			break;
		case tw::grid::high:
			sgn = -1;
			break;
	}
	switch (theBoundaryCondition)
	{
		case tw::bc::fld::normalFluxFixed:
			// quantities known on cell walls that are fixed on the wall
			if (whichSide==tw::grid::low)
			{
				// If we are using 0 to store BC we can't force it
				//force[0][0] = 0.0;
				fold[3][1] = -1.0;
				force[1][1] = 0.0;
				force[2][2] = 0.0;
				coeff[2] = 0.0; // set to 1 to enable inhomogeneous
			}
			else
			{
				// If we are using 0 to store BC we can't force it
				//force[0][0] = 0.0;
				fold[2][0] = -1.0;
				force[1][1] = 0.0;
				coeff[1] = 0.0; // set to 1 to enable inhomogeneous
			}
			break;
		case tw::bc::fld::neumannWall:
			// quantities known at cell centers that are continued outside domain
			// strip_1 = strip_2 - BC <--> strip_2 - strip_1 = BC
			// If we are using 0 to store BC we can't force it
			//force[0][0] = 0.0;
			//force[0][3] = 1.0;
			//coeff[0] = -1.0;
			force[1][1] = 0.0;
			force[1][2] = 1.0;
			coeff[1] = 0.0; // set to -1 to enable inhomogeneous
			break;
		case tw::bc::fld::dirichletWall:
			// quantities known at cell centers which are fixed on cell walls
			// (e.g., a no-slip boundary condition for tangential fluid velocity)
			// strip_1 = 2*BC - strip_2 <--> (strip_1+strip_2)/2 = BC
			// If we are using 0 to store BC we can't force it
			//force[0][0] = 0.0;
			//force[0][3] = -1.0;
			//coeff[0] = 2.0;
			force[1][1] = 0.0;
			force[1][2] = -1.0;
			coeff[1] = 2.0;
			break;
		case tw::bc::fld::dirichletCell:
			// quantities known at cell centers which are fixed outside domain
			// If we are using 0 to store BC we can't force it
			//force[0][0] = 0.0;
			//coeff[0] = 1.0;
			fold[2][1] = 1.0;
			fold[3][0] = 1.0;
			force[1][1] = 0.0;
			coeff[1] = 1.0;
			break;
		case tw::bc::fld::none:
		case tw::bc::fld::periodic:
		case tw::bc::fld::natural:
			// leave as the identity
			break;
	}
}



/////////////
//         //
// ELEMENT //
//         //
/////////////


Element::Element()
{
	low = 0;
	high = 0;
}

Element::Element(tw::Int l,tw::Int h)
{
	low = l;
	high = h;
}

Element::Element(tw::Int i)
{
	low = i;
	high = i;
}

tw::Int Element::Components() const
{
	return high - low + 1;
}

Element Union(const Element& e1,const Element& e2)
{
	if (e1.low < e2.low)
		return Element(e1.low,e2.high);
	else
		return Element(e2.low,e1.high);
}


///////////////////
//               //
//  FIELD CLASS  //
//               //
///////////////////


Field::Field()
{
	for (int i=0;i<4;i++)
		num[i] = 0;
	totalCells = 0;
	boundaryCells = 0;
	ghostCells = 0;
	bufferState = 0;
	packedAxis = 3;
}

Field::~Field()
{
	#ifdef USE_OPENCL
	if (bufferState==1)
	{
		clReleaseMemObject(computeBuffer);
		clReleaseMemObject(boundaryBuffer);
		clReleaseMemObject(boundaryMapBuffer);
		clReleaseMemObject(ghostBuffer);
		clReleaseMemObject(ghostMapBuffer);
	}
	#endif
}

/// @brief initialize the field, must call this before using
/// @param components number of floats to arrange along axis 0
/// @param ds discrete space to use for this Field, data is copied
/// @param task pointer to the concurrent task with this Field's domain
/// @param axis use this as the packed axis, only affects storage pattern
void Field::Initialize(tw::Int components,const DiscreteSpace& ds,Task *task,const tw::grid::axis& axis)
{
	DiscreteSpace::operator=(ds);
	this->task = task;
	packedAxis = tw::grid::naxis(axis);
	totalCells = num[1]*num[2]*num[3];
	num[0] = components;
	bc0.resize(4,components);
	bc1.resize(4,components);
	for (tw::Int i=0;i<4;i++)
		for (tw::Int c=0;c<components;c++)
		{
			bc0(i,c).Set(tw::bc::fld::none,tw::grid::low);
			bc1(i,c).Set(tw::bc::fld::none,tw::grid::high);
		}
	if (packedAxis==0)
	{
		stride[0] = 1;
		//stride[1] = num[0];
		//stride[2] = num[0]*num[1];
		//stride[3] = num[0]*num[1]*num[2];
		stride[1] = num[0]*num[2]*num[3];
		stride[2] = num[0]*num[3];
		stride[3] = num[0];
	}
	if (packedAxis==1)
	{
		stride[0] = num[1]*num[2]*num[3];
		stride[1] = 1;
		stride[2] = num[1];
		stride[3] = num[1]*num[2];
	}
	if (packedAxis==2)
	{
		stride[0] = num[1]*num[2]*num[3];
		stride[1] = num[2];
		stride[2] = 1;
		stride[3] = num[1]*num[2];
	}
	if (packedAxis==3)
	{
		stride[0] = num[1]*num[2]*num[3];
		stride[1] = num[2]*num[3];
		stride[2] = num[3];
		stride[3] = 1;
	}
	stride[1] *= dim[1]==1 ? 0 : 1;
	stride[2] *= dim[2]==1 ? 0 : 1;
	stride[3] *= dim[3]==1 ? 0 : 1;
	// note : boundary and ghost cells have redundant information
	// (this simplifies the storage scheme, but demands operations in compute buffer be atomic)
	// boundary cells means ghost cells + as many adjacent interior cells as there are ghost cells
	tw::Int boundaryCellsVec[4],ghostCellsVec[4];
	boundaryCellsVec[1] = dim[1]>1 ? 4*layers[1]*num[2]*num[3] : 0;
	boundaryCellsVec[2] = dim[2]>1 ? 4*layers[2]*num[1]*num[3] : 0;
	boundaryCellsVec[3] = dim[3]>1 ? 4*layers[3]*num[1]*num[2] : 0;
	ghostCellsVec[1] = dim[1]>1 ? 2*layers[1]*num[2]*num[3] : 0;
	ghostCellsVec[2] = dim[2]>1 ? 2*layers[2]*num[1]*num[3] : 0;
	ghostCellsVec[3] = dim[3]>1 ? 2*layers[3]*num[1]*num[2] : 0;
	boundaryCells = boundaryCellsVec[1] + boundaryCellsVec[2] + boundaryCellsVec[3];
	ghostCells = ghostCellsVec[1] + ghostCellsVec[2] + ghostCellsVec[3];
	array.resize(totalCells*num[0]);
	boundaryData.resize(boundaryCells*num[0]);
	ghostData.resize(ghostCells*num[0]);
	boundaryDataIndexMap.resize(boundaryCells*num[0]);
	ghostDataIndexMap.resize(ghostCells*num[0]);

	// Set up the index maps going from coordinate space to memory space.
	// Ordering is not important as long as we gather all the right cells.

	tw::Int bOffset=0,gOffset=0;
	for (tw::Int c=0;c<num[0];c++)
		for (tw::Int ax=1;ax<=3;ax++)
		{
			if (dim[ax]>1)
				for (auto s : StripRange(*this,ax,strongbool::yes))
				{
					for (tw::Int i=0;i<2*layers[ax];i++)
					{
						boundaryDataIndexMap[bOffset+i] = s.Index(lfg[ax]+i,c,stride);
						boundaryDataIndexMap[bOffset+2*layers[ax]+i] = s.Index(ufg[ax]-2*layers[ax]+i+1,c,stride);
					}
					for (tw::Int i=0;i<layers[ax];i++)
					{
						ghostDataIndexMap[gOffset+i] = s.Index(lfg[ax]+i,c,stride);
						ghostDataIndexMap[gOffset+layers[ax]+i] = s.Index(ufg[ax]-layers[ax]+i+1,c,stride);
					}
					bOffset += 4*layers[ax];
					gOffset += 2*layers[ax];
				}
		}
}

void Field::MultiplyCellVolume(const MetricSpace& m)
{
	for (auto cell : EntireCellRange(*this))
		for (tw::Int c=0;c<num[0];c++)
			(*this)(cell,c) *= m.dS(cell,0);
}

void Field::DivideCellVolume(const MetricSpace& m)
{
	for (auto cell : EntireCellRange(*this))
		for (tw::Int c=0;c<num[0];c++)
			(*this)(cell,c) /= m.dS(cell,0);
}

void Field::Shift(const Element& e,const tw::strip& s,tw::Int cells,const tw::Float* incoming)
{
	// Propagate field pattern to the right if <cells> positive, to the left if negative
	// argument <incoming> is value to inject into left or right ghost cell

	tw::Int i,c,ax=s.Axis();

	if (cells>0)
	{
		for (i=UFG(ax);i>=LFG(ax)+cells;i--)
			for (c=e.low;c<=e.high;c++)
				(*this)(s,i,c) = (*this)(s,i-cells,c);
		for (i=LFG(ax);i<LFG(ax)+cells;i++)
			for (c=e.low;c<=e.high;c++)
				(*this)(s,i,c) = incoming[c-e.low];
	}
	if (cells<0)
	{
		for (i=LFG(ax);i<=UFG(ax)+cells;i++)
			for (c=e.low;c<=e.high;c++)
				(*this)(s,i,c) = (*this)(s,i-cells,c);
		for (i=UFG(ax)+1+cells;i<=UFG(ax);i++)
			for (c=e.low;c<=e.high;c++)
				(*this)(s,i,c) = incoming[c-e.low];
	}
}

void Field::Shift(const Element& e,const tw::strip& s,tw::Int cells,const tw::Float& incoming)
{
	// Propagate field pattern to the right if <cells> positive, to the left if negative
	// argument <incoming> is value to inject into left or right ghost cell

	tw::Int i,c,ax=s.Axis();

	if (cells>0)
	{
		for (i=UFG(ax);i>=LFG(ax)+cells;i--)
			for (c=e.low;c<=e.high;c++)
				(*this)(s,i,c) = (*this)(s,i-cells,c);
		for (i=LFG(ax);i<LFG(ax)+cells;i++)
			for (c=e.low;c<=e.high;c++)
				(*this)(s,i,c) = incoming;
	}
	if (cells<0)
	{
		for (i=LFG(ax);i<=UFG(ax)+cells;i++)
			for (c=e.low;c<=e.high;c++)
				(*this)(s,i,c) = (*this)(s,i-cells,c);
		for (i=UFG(ax)+1+cells;i<=UFG(ax);i++)
			for (c=e.low;c<=e.high;c++)
				(*this)(s,i,c) = incoming;
	}
}

//////////////////////////////////
//                              //
//  Global Boundary Conditions  //
//                              //
//////////////////////////////////


void Field::SetBoundaryConditions(const Element& e,const tw::grid::axis& axis,tw::bc::fld low,tw::bc::fld high)
{
	tw::Int i;
	tw::Int ax = tw::grid::naxis(axis);

	for (i=e.low;i<=e.high;i++)
	{
		if (task->n0[ax]==MPI_PROC_NULL && dim[ax]>1)
			bc0(ax,i).Set(low,tw::grid::low);
		if (task->n1[ax]==MPI_PROC_NULL && dim[ax]>1)
			bc1(ax,i).Set(high,tw::grid::high);
	}
}

void Field::ZeroGhostCells(const Element& e)
{
	tw::Int ax,c,s;
	for (c=e.low;c<=e.high;c++)
		for (ax=1;ax<=3;ax++)
			for (auto strip : StripRange(*this,ax,strongbool::yes))
				for (s=0;s<layers[ax];s++)
				{
					(*this)(strip,lfg[ax]+s,c) = 0.0;
					(*this)(strip,dim[ax]+s+1,c) = 0.0;
				}
}

void Field::ApplyFoldingCondition(const Element& e)
{
	tw::Int ax,c;
	for (c=e.low;c<=e.high;c++)
		for (ax=1;ax<=3;ax++)
			if (num[ax]>1)
				for (auto strip : StripRange(*this,ax,strongbool::yes))
				{
					bc0(ax,c).FoldingOperation(&(*this)(strip,lfg[ax],c),stride[ax]);
					bc1(ax,c).FoldingOperation(&(*this)(strip,ufg[ax],c),stride[ax]);
				}
}

void Field::ApplyBoundaryCondition(const Element& e,bool homogeneous)
{
	tw::Int ax,c;
	tw::Float mult = homogeneous ? 0.0 : 1.0;
	for (c=e.low;c<=e.high;c++)
		for (ax=1;ax<=3;ax++)
			if (num[ax]>1)
				for (auto strip : StripRange(*this,ax,strongbool::yes))
				{
					bc0(ax,c).ForcingOperation(&(*this)(strip,lfg[ax],c),stride[ax],(*this)(strip,lfg[ax],c) * mult);
					bc1(ax,c).ForcingOperation(&(*this)(strip,ufg[ax],c),stride[ax],(*this)(strip,ufg[ax],c) * mult);
				}
}

void Field::BoundaryDataToField()
{
	for (tw::Int i=0;i<boundaryData.size();i++)
		array[boundaryDataIndexMap[i]] = boundaryData[i];
}

void Field::FieldToBoundaryData()
{
	for (tw::Int i=0;i<boundaryData.size();i++)
		boundaryData[i] = array[boundaryDataIndexMap[i]];
}

void Field::FieldToGhostData()
{
	for (tw::Int i=0;i<ghostData.size();i++)
		ghostData[i] = array[ghostDataIndexMap[i]];
}

/////////////////////////////
//                         //
//  GPU Related Routines   //
//                         //
/////////////////////////////

#ifdef USE_OPENCL

void Field::InitializeComputeBuffer()
{
	cl_int err;
	if (bufferState==1)
	{
		clReleaseMemObject(computeBuffer);
		clReleaseMemObject(boundaryBuffer);
		clReleaseMemObject(boundaryMapBuffer);
		clReleaseMemObject(ghostBuffer);
		clReleaseMemObject(ghostMapBuffer);
	}
	computeBuffer = clCreateBuffer(task->context,CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR,TotalBytes(),&array[0],&err);
	boundaryBuffer = clCreateBuffer(task->context,CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR,BoundaryBytes(),&boundaryData[0],&err);
	boundaryMapBuffer = clCreateBuffer(task->context,CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,BoundaryMapBytes(),&boundaryDataIndexMap[0],&err);
	ghostBuffer = clCreateBuffer(task->context,CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR,GhostBytes(),&ghostData[0],&err);
	ghostMapBuffer = clCreateBuffer(task->context,CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,GhostMapBytes(),&ghostDataIndexMap[0],&err);
	bufferState = 1;

	if (err)
		throw tw::FatalError("Device buffer failure in Field.");
}

void Field::SendToComputeBuffer()
{
	clEnqueueWriteBuffer(task->commandQueue,computeBuffer,CL_TRUE,0,TotalBytes(),&array[0],  0,NULL,NULL);
	clFinish(task->commandQueue);
}

void Field::ReceiveFromComputeBuffer()
{
	clEnqueueReadBuffer(task->commandQueue,computeBuffer,CL_TRUE,0,TotalBytes(),&array[0],  0,NULL,NULL);
	clFinish(task->commandQueue);
}

void Field::SendBoundaryCellsToComputeBuffer()
{
	FieldToBoundaryData();
	size_t elements = boundaryData.size();
	clEnqueueWriteBuffer(task->commandQueue,boundaryBuffer,CL_TRUE,0,BoundaryBytes(),&boundaryData[0]  ,0,NULL,NULL);
	clFinish(task->commandQueue);
	clSetKernelArg(task->k_boundaryToField,0,sizeof(computeBuffer),&computeBuffer);
	clSetKernelArg(task->k_boundaryToField,1,sizeof(boundaryBuffer),&boundaryBuffer);
	clSetKernelArg(task->k_boundaryToField,2,sizeof(boundaryMapBuffer),&boundaryMapBuffer);
	clEnqueueNDRangeKernel(task->commandQueue,task->k_boundaryToField,1,NULL,&elements,NULL,  0,NULL,NULL);
	clFinish(task->commandQueue);
}

void Field::ReceiveBoundaryCellsFromComputeBuffer()
{
	size_t elements = boundaryData.size();
	clSetKernelArg(task->k_fieldToBoundary,0,sizeof(computeBuffer),&computeBuffer);
	clSetKernelArg(task->k_fieldToBoundary,1,sizeof(boundaryBuffer),&boundaryBuffer);
	clSetKernelArg(task->k_fieldToBoundary,2,sizeof(boundaryMapBuffer),&boundaryMapBuffer);
	clEnqueueNDRangeKernel(task->commandQueue,task->k_fieldToBoundary,1,NULL,&elements,NULL,  0,NULL,NULL);
	clFinish(task->commandQueue);
	clEnqueueReadBuffer(task->commandQueue,boundaryBuffer,CL_TRUE,0,BoundaryBytes(),&boundaryData[0]  ,0,NULL,NULL);
	clFinish(task->commandQueue);
	BoundaryDataToField();
}

tw::Float Field::CellValueInComputeBuffer(tw::Int i,tw::Int j,tw::Int k,tw::Int c)
{
	tw::Float ans;
	tw::Int idx = c*stride[0] + (i-lfg[1])*stride[1] + (j-lfg[2])*stride[2] + (k-lfg[3])*stride[3];
	clEnqueueReadBuffer(task->commandQueue,computeBuffer,CL_TRUE,sizeof(tw::Float)*idx,sizeof(tw::Float),&ans,  0,NULL,NULL);
	clFinish(task->commandQueue);
	return ans;
}

void Field::UpdateGhostCellsInComputeBuffer(const Element& e)
{
	ReceiveBoundaryCellsFromComputeBuffer();
	CopyFromNeighbors(e);
	ApplyBoundaryCondition(e);
	SendBoundaryCellsToComputeBuffer();
}

void Field::SendGhostCellsToComputeBuffer()
{
	FieldToGhostData();
	size_t elements = ghostData.size();
	clEnqueueWriteBuffer(task->commandQueue,ghostBuffer,CL_TRUE,0,GhostBytes(),&ghostData[0]  ,0,NULL,NULL);
	clFinish(task->commandQueue);
	clSetKernelArg(task->k_ghostToField,0,sizeof(computeBuffer),&computeBuffer);
	clSetKernelArg(task->k_ghostToField,1,sizeof(ghostBuffer),&ghostBuffer);
	clSetKernelArg(task->k_ghostToField,2,sizeof(ghostMapBuffer),&ghostMapBuffer);
	clEnqueueNDRangeKernel(task->commandQueue,task->k_ghostToField,1,NULL,&elements,NULL,  0,NULL,NULL);
	clFinish(task->commandQueue);
}

void Field::ZeroGhostCellsInComputeBuffer()
{
	size_t cells = ghostCells * num[0];
	clSetKernelArg(task->k_zeroGhostCells,0,sizeof(computeBuffer),&computeBuffer);
	clSetKernelArg(task->k_zeroGhostCells,1,sizeof(ghostMapBuffer),&ghostMapBuffer);
	clEnqueueNDRangeKernel(task->commandQueue,task->k_zeroGhostCells,1,NULL,&cells,NULL,  0,NULL,NULL);
	clFinish(task->commandQueue);
}

void SwapComputeBuffers(Field& f1,Field& f2)
{
	clSetKernelArg(f1.task->k_swapBuffers,0,sizeof(cl_mem),&f1.computeBuffer);
	clSetKernelArg(f1.task->k_swapBuffers,1,sizeof(cl_mem),&f2.computeBuffer);
	f1.ElementUpdateProtocol(f1.task->k_swapBuffers,f1.task->commandQueue);
}

void CopyComputeBuffer(Field& dst,Field& src)
{
	clEnqueueCopyBuffer(dst.task->commandQueue,src.computeBuffer,dst.computeBuffer,0,0,dst.TotalBytes(),  0,NULL,NULL);
	clFinish(dst.task->commandQueue);
}

void CopyComplexComputeBufferMod2(Field& dst,Field& src)
{
	// dst is a real scalar buffer, src is a c-packed complex scalar buffer.
	// ElementUpdateProtocol must be called on the real buffer, not the complex buffer.
	clSetKernelArg(dst.task->k_complexMod2,0,sizeof(cl_mem),&dst.computeBuffer);
	clSetKernelArg(dst.task->k_complexMod2,1,sizeof(cl_mem),&src.computeBuffer);
	dst.ElementUpdateProtocol(dst.task->k_complexMod2,dst.task->commandQueue);
}

void Field::DestructiveComplexMod2ComputeBuffer()
{
	// The field is assumed to be a c-packed complex scalar buffer.
	clSetKernelArg(task->k_destructiveComplexMod2,0,sizeof(cl_mem),&computeBuffer);
	CellUpdateProtocol(task->k_destructiveComplexMod2,task->commandQueue);
}

void Field::FillComputeBufferVec4(const Element& e,tw::vec4& A)
{
	// Overwrite elements in the set e
	// Preserve elements not in the set e
	// e must be in the range of a 4-vector
	tw::Int i;
	// Complex arguments are a device to encode which elements to preserve
	tw::Complex Ac[4];
	for (i=0;i<4;i++)
		Ac[i] = tw::Complex(0.0,1.0);
	for (i=e.low;i<=e.high;i++)
		Ac[i] = tw::Complex(A[i],0.0);
	clSetKernelArg(task->k_fillVec4Field,0,sizeof(cl_mem),&computeBuffer);
	clSetKernelArg(task->k_fillVec4Field,1,sizeof(tw::Complex),&Ac[0]);
	clSetKernelArg(task->k_fillVec4Field,2,sizeof(tw::Complex),&Ac[1]);
	clSetKernelArg(task->k_fillVec4Field,3,sizeof(tw::Complex),&Ac[2]);
	clSetKernelArg(task->k_fillVec4Field,4,sizeof(tw::Complex),&Ac[3]);
	CellUpdateProtocol(task->k_fillVec4Field,task->commandQueue);
}

tw::Float Field::DestructiveSumComputeBuffer()
{
	// add up all the components in all the cells (ghost cells included)
	size_t elements = totalCells*num[0];
	clSetKernelArg(task->k_destructiveSum,0,sizeof(cl_mem),&computeBuffer);
	while (elements>1)
	{
		clEnqueueNDRangeKernel(task->commandQueue,task->k_destructiveSum,1,NULL,&elements,NULL,0,NULL,NULL);
		clFinish(task->commandQueue);
		elements /= 2;
	}
	return CellValueInComputeBuffer(0,0,0,0);
}

tw::Float Field::DestructiveNorm1ComputeBuffer()
{
	// add up absolute value of all the components in all the cells (ghost cells included)
	size_t elements = totalCells*num[0];
	clSetKernelArg(task->k_destructiveNorm1,0,sizeof(cl_mem),&computeBuffer);
	while (elements>1)
	{
		clEnqueueNDRangeKernel(task->commandQueue,task->k_destructiveNorm1,1,NULL,&elements,NULL,0,NULL,NULL);
		clFinish(task->commandQueue);
		elements /= 2;
	}
	return CellValueInComputeBuffer(0,0,0,0);
}

void Field::WeightComputeBufferByVolume(MetricSpace& ms,tw::Float inv)
{
	clSetKernelArg(task->k_weightByVolume,0,sizeof(cl_mem),&computeBuffer);
	clSetKernelArg(task->k_weightByVolume,1,sizeof(cl_mem),&ms.metricsBuffer);
	clSetKernelArg(task->k_weightByVolume,3,sizeof(tw::Float),&inv);
	for (tw::Int c=0;c<num[0];c++)
	{
		clSetKernelArg(task->k_weightByVolume,2,sizeof(tw::Int),&c);
		LocalUpdateProtocol(task->k_weightByVolume,task->commandQueue);
	}
}

void Field::MADDComputeBuffer(tw::Float m,tw::Float a)
{
	clSetKernelArg(task->k_MADD,0,sizeof(cl_mem),&computeBuffer);
	clSetKernelArg(task->k_MADD,1,sizeof(tw::Float),&m);
	clSetKernelArg(task->k_MADD,2,sizeof(tw::Float),&a);
	ElementUpdateProtocol(task->k_MADD,task->commandQueue);
}

#endif

///////////////////////////////////////
//                                   //
//  Boundaries with Message Passing  //
//                                   //
///////////////////////////////////////


void Field::StripCopyProtocol(tw::Int axis,tw::Int shift,Slice<tw::Float> *planeIn,Slice<tw::Float> *planeOut,bool add)
{
	tw::Int src,dst;
	bool odd;

	task->strip[axis].Shift(1,shift,&src,&dst);
	odd = task->strip[axis].Get_rank() % 2;

	if (dst!=MPI_PROC_NULL)
		LoadDataIntoSlice<tw::Float>(planeOut);

	if (odd)
	{
		task->strip[axis].Recv(planeIn->Buffer(),planeIn->BufferSize(),src);
		task->strip[axis].Send(planeOut->Buffer(),planeOut->BufferSize(),dst);
	}
	else
	{
		task->strip[axis].Send(planeOut->Buffer(),planeOut->BufferSize(),dst);
		task->strip[axis].Recv(planeIn->Buffer(),planeIn->BufferSize(),src);
	}

	if (src!=MPI_PROC_NULL)
	{
		if (add)
			AddDataFromSlice<tw::Float>(planeIn);
		else
			SaveDataFromSlice<tw::Float>(planeIn);
	}
}

void Field::DownwardCopy(const tw::grid::axis& axis,const Element& e,tw::Int cells)
{
	Slice<tw::Float> *planeIn,*planeOut;
	const tw::Int ax = tw::grid::naxis(axis);
	tw::Int b[6];

	if (dim[ax]==1) return;
	for (tw::Int i=1;i<=3;i++)
	{
		b[(i-1)*2] = ax==i ? ung[i] : lfg[i];
		b[(i-1)*2+1] = ax==i ? ung[i] + cells - 1 : ufg[i];
	}
	planeIn = new Slice<tw::Float>(e,b[0],b[1],b[2],b[3],b[4],b[5]);
	planeOut = new Slice<tw::Float>(e,b[0],b[1],b[2],b[3],b[4],b[5]);
	planeOut->Translate(axis,lng[ax] + 1 - ung[ax]);

	StripCopyProtocol(tw::grid::naxis(axis),-1,planeIn,planeOut,false);

	delete planeOut;
	delete planeIn;
}

void Field::UpwardCopy(const tw::grid::axis& axis,const Element& e,tw::Int cells)
{
	Slice<tw::Float> *planeIn,*planeOut;
	const tw::Int ax = tw::grid::naxis(axis);
	tw::Int b[6];

	if (dim[ax]==1) return;
	for (tw::Int i=1;i<=3;i++)
	{
		b[(i-1)*2] = ax==i ? lng[i] - cells + 1 : lfg[i];
		b[(i-1)*2+1] = ax==i ? lng[i] : ufg[i];
	}
	planeIn = new Slice<tw::Float>(e,b[0],b[1],b[2],b[3],b[4],b[5]);
	planeOut = new Slice<tw::Float>(e,b[0],b[1],b[2],b[3],b[4],b[5]);
	planeOut->Translate(axis,ung[ax] - lng[ax] - 1);

	StripCopyProtocol(tw::grid::naxis(axis),1,planeIn,planeOut,false);

	delete planeOut;
	delete planeIn;
}

void Field::DownwardDeposit(const tw::grid::axis& axis,const Element& e,tw::Int cells)
{
	// deposits work by sending both ghost cells and edge cells one way and adding.
	// Before accepting the return message these cells are zeroed.
	// The return message then effectively overwrites them with the correct data.
	Slice<tw::Float> *planeIn,*planeOut;
	const tw::Int ax = tw::grid::naxis(axis);
	tw::Int b[6];

	if (dim[ax]==1) return;
	for (tw::Int i=1;i<=3;i++)
	{
		b[(i-1)*2] = ax==i ? ung[i] - cells : lfg[i];
		b[(i-1)*2+1] = ax==i ? ung[i] + cells - 1 : ufg[i];
	}
	planeIn = new Slice<tw::Float>(e,b[0],b[1],b[2],b[3],b[4],b[5]);
	planeOut = new Slice<tw::Float>(e,b[0],b[1],b[2],b[3],b[4],b[5]);
	planeOut->Translate(axis,lng[ax] + 1 - ung[ax]);

	StripCopyProtocol(tw::grid::naxis(axis),-1,planeIn,planeOut,true);

	delete planeOut;
	delete planeIn;
}

void Field::UpwardDeposit(const tw::grid::axis& axis,const Element& e,tw::Int cells)
{
	// deposits work by sending both ghost cells and edge cells one way and adding.
	// Before accepting the return message these cells are zeroed.
	// The return message then effectively overwrites them with the correct data.
	Slice<tw::Float> *planeIn,*planeOut;
	const tw::Int ax = tw::grid::naxis(axis);
	tw::Int b[6];

	if (dim[ax]==1) return;
	for (tw::Int i=1;i<=3;i++)
	{
		b[(i-1)*2] = ax==i ? lng[i] + 1 - cells : lfg[i];
		b[(i-1)*2+1] = ax==i ? lng[i] + cells : ufg[i];
	}
	planeIn = new Slice<tw::Float>(e,b[0],b[1],b[2],b[3],b[4],b[5]);
	planeOut = new Slice<tw::Float>(e,b[0],b[1],b[2],b[3],b[4],b[5]);
	planeOut->Translate(axis,ung[ax] - lng[ax] - 1);

	StripCopyProtocol(tw::grid::naxis(axis),1,planeIn,planeOut,true);

	delete planeOut;
	delete planeIn;
}

void Field::CopyFromNeighbors(const Element& e)
{
	DownwardCopy(tw::grid::x,e,1);
	UpwardCopy(tw::grid::x,e,1);

	DownwardCopy(tw::grid::y,e,1);
	UpwardCopy(tw::grid::y,e,1);

	DownwardCopy(tw::grid::z,e,1);
	UpwardCopy(tw::grid::z,e,1);
}

void Field::DepositFromNeighbors(const Element& e)
{
	for (tw::Int ax=1;ax<=3;ax++)
		if (dim[ax]>1)
		{
			DownwardDeposit(tw::grid::enumaxis(ax),e,layers[ax]);
			if (task->n0[ax]!=MPI_PROC_NULL)
			{
				tw::Int ub = lfg[ax] + 2*layers[ax] - 1;
				tw::Int x1 = ax==1 ? ub : ufg[1];
				tw::Int y1 = ax==2 ? ub : ufg[2];
				tw::Int z1 = ax==3 ? ub : ufg[3];
				Slice<tw::Float>* plane = new Slice<tw::Float>(e,lfg[1],x1,lfg[2],y1,lfg[3],z1);
				ZeroDataInField<tw::Float>(plane);
				delete plane;
			}
			UpwardDeposit(tw::grid::enumaxis(ax),e,layers[ax]);
		}
}

Slice<tw::Float>* Field::FormTransposeBlock(const Element& e,const tw::grid::axis& axis1,const tw::grid::axis& axis2,tw::Int start1,tw::Int end1,tw::Int start2,tw::Int end2)
{
	Slice<tw::Float> *ans;
	if (axis1==tw::grid::x)
	{
		if (axis2==tw::grid::y)
			ans = new Slice<tw::Float>(e,start1,end1,start2,end2,lfg[3],ufg[3]);
		else
			ans = new Slice<tw::Float>(e,start1,end1,lfg[2],ufg[2],start2,end2);
	}
	if (axis1==tw::grid::y)
	{
		if (axis2==tw::grid::x)
			ans = new Slice<tw::Float>(e,start2,end2,start1,end1,lfg[3],ufg[3]);
		else
			ans = new Slice<tw::Float>(e,lfg[1],ufg[1],start1,end1,start2,end2);
	}
	if (axis1==tw::grid::z)
	{
		if (axis2==tw::grid::x)
			ans = new Slice<tw::Float>(e,start2,end2,lfg[2],ufg[2],start1,end1);
		else
			ans = new Slice<tw::Float>(e,lfg[1],ufg[1],start2,end2,start1,end1);
	}
	return ans;
}


/// We have in mind a matrix whose rows consist of NODES arranged along axis1,
/// and whose columns consist of BLOCKS arranged along axis2.
/// In this picture, the global image of the data does not change.
/// Instead, the nodes are lined up along a new axis.
/// We choose the blocks so that #blocks <= #nodes.
/// 'target' is an externally owned field that will receive the transposed data
/// This routine resizes target when the forward transpose is invoked.
/// Upon reverse transposing, the same target should be passed in.
/// 'inversion' is 1 if forward transpose, -1 if reverse transpose
///
/// Ghost cell policy:
/// We send the ghost cells in the directions perp. to 'axis1'
/// The ghost cells are not sent in the direction parallel to 'axis1'
/// However, the ghost cells are sent in all directions if inversion=-1
void Field::Transpose(const Element& e,const tw::grid::axis& axis1,const tw::grid::axis& axis2,Field *target,tw::Int inversion)
{
	tw::Int dN0,dN1;
	Slice<tw::Float>* block;

	const tw::Int ax1 = tw::grid::naxis(axis1);
	const tw::Int ax2 = tw::grid::naxis(axis2);
	const tw::Int nodes = task->domains[ax1];
	const tw::Int thisNode = task->strip[ax1].Get_rank(); // assumes rank=coord
	const tw::Int cellsPerBlock = num[ax2]/nodes + 1;
	const tw::Int interiorCellsPerBlock = cellsPerBlock - 2*layers[ax2];
	const tw::Int fullBlocks = num[ax2]/cellsPerBlock;
	const tw::Int cellsRemaining = num[ax2] - cellsPerBlock*fullBlocks;
	std::vector<tw::Int> blockSize(nodes),start(nodes),end(nodes),offset(nodes);

	if (inversion==1)
	{
		dN0 = 1;
		dN1 = dim[ax1];
	}
	else
	{
		dN0 = 0;
		dN1 = dim[ax1] + 1;
	}

	for (tw::Int i=0;i<nodes;i++)
	{
		start[i] = i*cellsPerBlock + 1 - layers[ax2];
		if (i<fullBlocks)
			blockSize[i] = cellsPerBlock;
		if (i==fullBlocks)
			blockSize[i] = cellsRemaining;
		if (i>fullBlocks)
			blockSize[i] = 0;
		end[i] = start[i] + blockSize[i] - 1;
		offset[i] = start[i] + layers[ax2] - 1;
	}

	if (inversion==1)
	{
		tw::Int ax3,dim_T[4];
		ax3 = 1;
		while (ax3==ax1 || ax3==ax2) ax3++;
		dim_T[ax1] = dim[ax1]*nodes;
		dim_T[ax2] = interiorCellsPerBlock;
		dim_T[ax3] = dim[ax3];
		target->Initialize(e.Components(),DiscreteSpace(dim_T[1],dim_T[2],dim_T[3],corner,size),task);
	}

  // The message passing pattern is to have simultaneous exchanges between pairs.
  // Each pair is made up of a sub-diagonal and corresponding super-diagonal element

	// First do the diagonal, which involves no message passing

	if (blockSize[thisNode]>0)
	{
		const tw::Int j = thisNode;
		block = FormTransposeBlock(e,axis1,axis2,dN0,dN1,start[j],end[j]);
		if (inversion==-1)
		{
	    block->Translate(axis1,j*dim[ax1]);
	    block->Translate(axis2,-offset[j]);
	    target->LoadDataIntoSlice<tw::Float>(block);
	    block->Translate(axis1,-j*dim[ax1]);
	    block->Translate(axis2,offset[j]);
	    SaveDataFromSlice<tw::Float>(block);
		}
		else
		{
	    LoadDataIntoSlice<tw::Float>(block);
	    block->Translate(axis1,j*dim[ax1]);
	    block->Translate(axis2,-offset[j]);
	    target->SaveDataFromSlice<tw::Float>(block);
		}
		delete block;
	}

  // Now exchange data between sub/super diagonal pairs

	for (tw::Int i=1;i<nodes;i++) // Loop over sub-diagonals
		for (tw::Int k=0;k<2;k++) // each sub-diagonal has 2 independent groups
	    for (tw::Int l=k*i;l<nodes-i;l+=2*i) // Loop over super-blocks of this independent group
	      for (tw::Int s=0;s<i && s<nodes-l-i;s++) // Loop over blocks in each super-block
				{
					// sub-diagonal : (l+s , l+s+i) = (old node , old block) = (new block , new node)
					// super-diagonal : (l+s+i , l+s) = (old node , old block) = (new block , new node)
					if (thisNode==l+s+i)
					{
						// My role as super-diagonal old node
						if (blockSize[l+s]>0)
						{
							block = FormTransposeBlock(e,axis1,axis2,dN0,dN1,start[l+s],end[l+s]);
							if (inversion==-1)
							{
								task->strip[ax1].Recv(block->Buffer(),block->BufferSize(),l+s);
								SaveDataFromSlice<tw::Float>(block);
							}
							else
							{
								LoadDataIntoSlice<tw::Float>(block);
								task->strip[ax1].Send(block->Buffer(),block->BufferSize(),l+s);
							}
							delete block;
						}

						// My role as sub-diagonal new node
						if (blockSize[thisNode]>0)
						{
							block = FormTransposeBlock(e,axis1,axis2,dN0,dN1,start[thisNode],end[thisNode]);
							block->Translate(axis1,(l+s)*dim[ax1]);
							block->Translate(axis2,-offset[thisNode]);
							if (inversion==-1)
							{
								target->LoadDataIntoSlice<tw::Float>(block);
								task->strip[ax1].Send(block->Buffer(),block->BufferSize(),l+s);
							}
							else
							{
								task->strip[ax1].Recv(block->Buffer(),block->BufferSize(),l+s);
								target->SaveDataFromSlice<tw::Float>(block);
							}
							delete block;
						}
					}
					if (thisNode==l+s)
					{
						// My role as super-diagonal new node
						if (blockSize[thisNode]>0)
						{
							block = FormTransposeBlock(e,axis1,axis2,dN0,dN1,start[thisNode],end[thisNode]);
							block->Translate(axis1,(l+s+i)*dim[ax1]);
							block->Translate(axis2,-offset[thisNode]);
							if (inversion==-1)
							{
								target->LoadDataIntoSlice<tw::Float>(block);
								task->strip[ax1].Send(block->Buffer(),block->BufferSize(),l+s+i);
							}
							else
							{
								task->strip[ax1].Recv(block->Buffer(),block->BufferSize(),l+s+i);
								target->SaveDataFromSlice<tw::Float>(block);
							}
							delete block;
						}

						// My role as sub-diagonal old node
						if (blockSize[l+s+i]>0)
						{
							block = FormTransposeBlock(e,axis1,axis2,dN0,dN1,start[l+s+i],end[l+s+i]);
							if (inversion==-1)
							{
								task->strip[ax1].Recv(block->Buffer(),block->BufferSize(),l+s+i);
								SaveDataFromSlice<tw::Float>(block);
							}
							else
							{
								LoadDataIntoSlice<tw::Float>(block);
								task->strip[ax1].Send(block->Buffer(),block->BufferSize(),l+s+i);
							}
							delete block;
						}
					}
				}
}

/////////////////
//             //
//  Smoothing  //
//             //
/////////////////


void Field::SmoothingPass(tw::Int ax,const Element& e,const MetricSpace& ds,const tw::Float& X0,const tw::Float& X1,const tw::Float& X2)
{
	tw::Int i,c;
	tw::Float ansNow,temp;
	tw::grid::axis axs[4] = { tw::grid::t, tw::grid::x, tw::grid::y, tw::grid::z };

	for (c=e.low;c<=e.high;c++)
		if (dim[ax]>1)
			for (auto s : StripRange(*this,ax,strongbool::yes))
			{
				temp = (*this)(s,0,c);
				for (i=1;i<=dim[ax];i++)
				{
					ansNow = X0*temp + X1*(*this)(s,i,c) + X2*(*this)(s,i+1,c);
					temp = (*this)(s,i,c);
					(*this)(s,i,c) = ansNow;
				}
				bc0(ax,c).ForcingOperation(&(*this)(s,lfg[ax],c),stride[ax],0.0);
				bc1(ax,c).ForcingOperation(&(*this)(s,ufg[ax],c),stride[ax],0.0);
			}

	UpwardCopy(axs[ax],e,1);
	DownwardCopy(axs[ax],e,1);
}

void Field::Smooth(const Element& e,const MetricSpace& ds,tw::Int smoothPasses[4],tw::Int compPasses[4])
{
	for (tw::Int ax=1;ax<=3;ax++)
	{
		for (tw::Int ipass=0;ipass<smoothPasses[ax];ipass++)
			SmoothingPass(ax,e,ds,0.25,0.5,0.25);
		for (tw::Int ipass=0;ipass<compPasses[ax];ipass++)
			SmoothingPass(ax,e,ds,-1.25,3.5,-1.25);
	}
}

void Field::ReadCheckpoint(std::ifstream& inFile)
{
	DiscreteSpace::ReadCheckpoint(inFile);
	inFile.read((char*)&packedAxis,sizeof(tw::Int));
	inFile.read((char*)&bc0(0,0),sizeof(BoundaryCondition)*num[0]*4);
	inFile.read((char*)&bc1(0,0),sizeof(BoundaryCondition)*num[0]*4);
	inFile.read((char *)&array[0],sizeof(tw::Float)*totalCells*num[0]);
}

void Field::WriteCheckpoint(std::ofstream& outFile)
{
	DiscreteSpace::WriteCheckpoint(outFile);
	outFile.write((char *)&packedAxis,sizeof(tw::Int));
	outFile.write((char*)&bc0(0,0),sizeof(BoundaryCondition)*num[0]*4);
	outFile.write((char*)&bc1(0,0),sizeof(BoundaryCondition)*num[0]*4);
	outFile.write((char *)&array[0],sizeof(tw::Float)*totalCells*num[0]);
}

export void CopyBoundaryConditions(Field& dst,const Element& dstElement,Field& src,const Element& srcElement)
{
	tw::Int n;
	n = dstElement.Components();
	if (n!=srcElement.Components())
		throw tw::FatalError("Incompatible copy operation (boundary conditions)");
	dst.bc0 = src.bc0;
	dst.bc1 = src.bc1;
}

export void Swap(Field& f1,Field& f2)
{
	tw::Int i;
	tw::Float temp;
	for (i=0;i<f1.array.size();i++)
	{
		temp = f1.array[i];
		f1.array[i] = f2.array[i];
		f2.array[i] = temp;
	}
}

void Field::Swap(const Element& e1,const Element& e2)
{
	tw::Int c;
	tw::Float temp;
	for (c=0;c<e1.Components();c++)
	{

		for(auto cell : EntireCellRange(*this))
		{
			temp = (*this)(cell,e1.low+c);
			(*this)(cell,e1.low+c) = (*this)(cell,e2.low+c);
			(*this)(cell,e2.low+c) = temp;
		}
	}
}

export void CopyFieldData(Field& dst,const Element& e_dst,Field& src,const Element& e_src)
{
	for (tw::Int c=0;c<e_dst.Components();c++)
		for (auto cell : CellRange(dst,true))
			dst(cell,e_dst.low+c) = src(cell,e_src.low+c);
}

export void CopyGhostCellData(Field& dst,const Element& e_dst,Field& src,const Element& e_src)
{
	for (tw::Int c=0;c<e_dst.Components();c++)
		for (tw::Int ax=1;ax<=3;ax++)
			for (auto strip : StripRange(dst,ax,strongbool::yes))
				for (tw::Int i=0;i<dst.layers[ax];i++)
				{
					dst(strip,dst.lfg[ax]+i,e_dst.low+c) = src(strip,src.lfg[ax]+i,e_src.low+c);
					dst(strip,dst.ufg[ax]-i,e_dst.low+c) = src(strip,src.ufg[ax]-i,e_src.low+c);
				}
}

export void AddFieldData(Field& dst,const Element& e_dst,Field& src,const Element& e_src)
{
	tw::Int c;
	for (c=0;c<e_dst.Components();c++)
		for (auto cell : CellRange(dst,true))
			dst(cell,e_dst.low+c) += src(cell,e_src.low+c);
}

export void AddMulFieldData(Field& dst,const Element& e_dst,Field& src,const Element& e_src,tw::Float mul)
{
	tw::Int c;
	for (c=0;c<e_dst.Components();c++)
		for (auto cell : CellRange(dst,true))
			dst(cell,e_dst.low+c) += mul*src(cell,e_src.low+c);
}


//////////////////////////
//                      //
// SCALAR FIELD METHODS //
//                      //
//////////////////////////


tw::Float ScalarField::AxialEigenvalue(tw::Int z)
{
	// eigenvalues for sine and cosine transforms are the same
	return eigenvalue_FST(GlobalCellIndex(z,3)-1,GlobalDim(3),freq[3]);
}

tw::Float ScalarField::Eigenvalue(tw::Int x,tw::Int y)
{
	// eigenvalues for sine and cosine transforms are the same
	return eigenvalue_FST(GlobalCellIndex(x,1)-1,GlobalDim(1),freq[1]) + eigenvalue_FST(GlobalCellIndex(y,2)-1,GlobalDim(2),freq[2]);
}

tw::Float ScalarField::CyclicEigenvalue(tw::Int x,tw::Int y)
{
	tw::Float ans;
	x = GlobalCellIndex(x,1);
	y = GlobalCellIndex(y,2);

	ans = eigenvalue_RFFT(x-1,GlobalDim(1),freq[1]);

	if (x==1 || x==2)
		ans += eigenvalue_RFFT(y-1,GlobalDim(2),freq[2]);
	else
		ans += eigenvalue_CFFT(y-1,GlobalDim(2),freq[2]);

	return ans;
}

tw::Float ScalarField::CyclicEigenvalue(tw::Int x,tw::Int y,tw::Int z)
{
	tw::Float ans;
	x = GlobalCellIndex(x,1);
	y = GlobalCellIndex(y,2);
	z = GlobalCellIndex(z,3);

	ans = eigenvalue_RFFT(x-1,GlobalDim(1),freq[1]);

	if (x==1 || x==2)
		ans += eigenvalue_RFFT(y-1,GlobalDim(2),freq[2]);
	else
		ans += eigenvalue_CFFT(y-1,GlobalDim(2),freq[2]);

	if ((x==1 || x==2) && (y==1 || y==2))
		ans += eigenvalue_RFFT(z-1,GlobalDim(3),freq[3]);
	else
		ans += eigenvalue_CFFT(z-1,GlobalDim(3),freq[3]);

	return ans;
}

void ScalarField::AxialSineTransform()
{
	Field T;
	if (GlobalDim(3)>1)
	{
		Transpose(tw::grid::z,tw::grid::x,&T,1);
		#pragma omp parallel
		{
			for (auto strip : StripRange(T,3,strongbool::yes))
				SineTransform( &T(strip,1,0), T.Dim(3), T.Stride(3), 1 );
		}
		Transpose(tw::grid::z,tw::grid::x,&T,-1);
	}
}

void ScalarField::InverseAxialSineTransform()
{
	Field T;
	if (GlobalDim(3)>1)
	{
		Transpose(tw::grid::z,tw::grid::x,&T,1);
		#pragma omp parallel
		{
			for (auto strip : StripRange(T,3,strongbool::yes))
			{
				SineTransform( &T(strip,1,0), T.Dim(3), T.Stride(3), -1 );
				T(strip,T.LNG(3),0) = -T(strip,2,0);
				T(strip,T.UNG(3),0) = 0.0;
			}
		}
		Transpose(tw::grid::z,tw::grid::x,&T,-1);
	}
}

void ScalarField::TransverseCosineTransform()
{
	Field T;
	tw::grid::axis axis1,axis2;
	for (tw::Int ax=1;ax<=2;ax++)
	{
		axis1 = tw::grid::enumaxis(ax);
		axis2 = tw::grid::z;
		if (GlobalDim(ax)>1)
		{
			Transpose(axis1,axis2,&T,1);
			#pragma omp parallel
			{
				for (auto strip : StripRange(T,ax,strongbool::yes))
					CosineTransform( &T(strip,1,0), T.Dim(ax), T.Stride(ax), 1 );
			}
			Transpose(axis1,axis2,&T,-1);
		}
	}
}

void ScalarField::InverseTransverseCosineTransform()
{
	Field T;
	tw::grid::axis axis1,axis2;
	for (tw::Int ax=2;ax>=1;ax--)
	{
		axis1 = tw::grid::enumaxis(ax);
		axis2 = tw::grid::z;
		if (GlobalDim(ax)>1)
		{
			Transpose(axis1,axis2,&T,1);
			#pragma omp parallel
			{
				for (auto strip : StripRange(T,ax,strongbool::yes))
				{
					CosineTransform( &T(strip,1,0), T.Dim(ax), T.Stride(ax), -1 );
					T(strip,T.LNG(ax),0) = T(strip,1,0);
					T(strip,T.UNG(ax),0) = T(strip,T.Dim(ax),0);
				}
			}
			Transpose(axis1,axis2,&T,-1);
		}
	}
}

void ScalarField::TransverseSineTransform()
{
	Field T;
	tw::grid::axis axis1,axis2;
	for (tw::Int ax=1;ax<=2;ax++)
	{
		axis1 = tw::grid::enumaxis(ax);
		axis2 = tw::grid::z;
		if (GlobalDim(ax)>1)
		{
			Transpose(axis1,axis2,&T,1);
			#pragma omp parallel
			{
				for (auto strip : StripRange(T,ax,strongbool::yes))
					SineTransform( &T(strip,1,0), T.Dim(ax), T.Stride(ax), 1 );
			}
			Transpose(axis1,axis2,&T,-1);
		}
	}
}

void ScalarField::InverseTransverseSineTransform()
{
	Field T;
	tw::grid::axis axis1,axis2;
	for (tw::Int ax=2;ax>=1;ax--)
	{
		axis1 = tw::grid::enumaxis(ax);
		axis2 = tw::grid::z;
		if (GlobalDim(ax)>1)
		{
			Transpose(axis1,axis2,&T,1);
			#pragma omp parallel
			{
				for (auto strip : StripRange(T,ax,strongbool::yes))
				{
					SineTransform( &T(strip,1,0), T.Dim(ax), T.Stride(ax), -1 );
					T(strip,T.LNG(ax),0) = -T(strip,2,0);
					T(strip,T.UNG(ax),0) = 0.0;
				}
			}
			Transpose(axis1,axis2,&T,-1);
		}
	}
}

void ScalarField::TransverseFFT()
{
	Field T;

	if (GlobalDim(1)>1)
	{
		Transpose(tw::grid::x,tw::grid::z,&T,1);
		#pragma omp parallel
		{
			for (auto strip : StripRange(T,1,strongbool::yes))
				RealFFT( &T(strip,1,0), T.Dim(1), T.Stride(1), 1);
		}
		Transpose(tw::grid::x,tw::grid::z,&T,-1);
	}

	if (GlobalDim(2)>1)
	{
		Transpose(tw::grid::y,tw::grid::z,&T,1);
		const tw::Int xDim=T.Dim(1),zN0=T.LFG(3),zN1=T.UFG(3);
		#pragma omp parallel for collapse(2) schedule(static)
		for (tw::Int i=1;i<=xDim;i+=2) // can't include ghost cells due to complex numbers; instead do copy ops below
			for (tw::Int k=zN0;k<=zN1;k++)
			{
				if (GlobalCellIndex(i,1)==1)
				{
					RealFFT( &T(i,1,k,0), T.Dim(2), T.Stride(2), 1);
					RealFFT( &T(i+1,1,k,0), T.Dim(2), T.Stride(2), 1);
				}
				else
					ComplexFFT( &T(i,1,k,0), &T(i+1,1,k,0), T.Dim(2), T.Stride(2), 1.0);
			}
		T.UpwardCopy(tw::grid::x,1);
		T.DownwardCopy(tw::grid::x,1);
		Transpose(tw::grid::y,tw::grid::z,&T,-1);
	}
}

void ScalarField::InverseTransverseFFT()
{
	Field T;

	if (GlobalDim(2)>1)
	{
		Transpose(tw::grid::y,tw::grid::z,&T,1);
		const tw::Int xDim=T.Dim(1),zN0=T.LFG(3),zN1=T.UFG(3);
		#pragma omp parallel for collapse(2) schedule(static)
		for (tw::Int i=1;i<=xDim;i+=2) // can't include ghost cells due to complex numbers; instead do copy ops below
			for (tw::Int k=zN0;k<=zN1;k++)
			{
				if (GlobalCellIndex(i,1)==1)
				{
					RealFFT( &T(i,1,k,0), T.Dim(2), T.Stride(2), -1);
					T(i,T.LNG(2),k,0) = T(i,T.Dim(2),k,0);
					T(i,T.UNG(2),k,0) = T(i,1,k,0);
					RealFFT( &T(i+1,1,k,0), T.Dim(2), T.Stride(2), -1);
					T(i+1,T.LNG(2),k,0) = T(i+1,T.Dim(2),k,0);
					T(i+1,T.UNG(2),k,0) = T(i+1,1,k,0);
				}
				else
				{
					ComplexFFT( &T(i,1,k,0), &T(i+1,1,k,0), T.Dim(2), T.Stride(2), -1.0);
					T(i,T.LNG(2),k,0) = T(i,T.Dim(2),k,0);
					T(i,T.UNG(2),k,0) = T(i,1,k,0);
					T(i+1,T.LNG(2),k,0) = T(i+1,T.Dim(2),k,0);
					T(i+1,T.UNG(2),k,0) = T(i+1,1,k,0);
				}
			}
		T.UpwardCopy(tw::grid::x,1);
		T.DownwardCopy(tw::grid::x,1);
		Transpose(tw::grid::y,tw::grid::z,&T,-1);
	}

	if (GlobalDim(1)>1)
	{
		Transpose(tw::grid::x,tw::grid::z,&T,1);
		#pragma omp parallel
		{
			for (auto strip : StripRange(T,1,strongbool::yes))
			{
				RealFFT( &T(strip,1,0), T.Dim(1), T.Stride(1), -1);
				T(strip,T.LNG(1),0) = T(strip,T.Dim(1),0);
				T(strip,T.UNG(1),0) = T(strip,1,0);
			}
		}
		Transpose(tw::grid::x,tw::grid::z,&T,-1);
	}
}


///////////////////////////
//                       //
// COMPLEX FIELD METHODS //
//                       //
///////////////////////////


tw::Float ComplexField::CyclicEigenvalue(tw::Int x,tw::Int y)
{
	x = GlobalCellIndex(x,1);
	y = GlobalCellIndex(y,2);
	return eigenvalue_CFFT(x-1,GlobalDim(1),freq[1]) + eigenvalue_CFFT(y-1,GlobalDim(2),freq[2]);
}

tw::Float ComplexField::CyclicEigenvalue(tw::Int x,tw::Int y,tw::Int z)
{
	x = GlobalCellIndex(x,1);
	y = GlobalCellIndex(y,2);
	z = GlobalCellIndex(z,3);
	return eigenvalue_CFFT(x-1,GlobalDim(1),freq[1]) + eigenvalue_CFFT(y-1,GlobalDim(2),freq[2]) + eigenvalue_CFFT(z-1,GlobalDim(3),freq[3]);
}

void ComplexField::TransverseFFT()
{
	Field T;
	tw::grid::axis axis1,axis2;
	for (tw::Int ax=1;ax<=2;ax++)
	{
		axis1 = tw::grid::enumaxis(ax);
		axis2 = tw::grid::z;
		if (GlobalDim(ax)>1)
		{
			Transpose(axis1,axis2,&T,1);
			#pragma omp parallel
			{
				for (auto strip : StripRange(T,ax,strongbool::yes))
					ComplexFFT( &T(strip,1,0), &T(strip,1,1), T.Dim(ax), T.Stride(ax), 1.0 );
			}
			Transpose(axis1,axis2,&T,-1);
		}
	}
}

void ComplexField::InverseTransverseFFT()
{
	Field T;
	tw::grid::axis axis1,axis2;
	for (tw::Int ax=2;ax>=1;ax--)
	{
		axis1 = tw::grid::enumaxis(ax);
		axis2 = tw::grid::z;
		if (GlobalDim(ax)>1)
		{
			Transpose(axis1,axis2,&T,1);
			#pragma omp parallel
			{
				for (auto strip : StripRange(T,ax,strongbool::yes))
				{
					ComplexFFT( &T(strip,1,0), &T(strip,1,1), T.Dim(ax), T.Stride(ax), -1.0 );
					T(strip,T.LNG(ax),0) = T(strip,T.Dim(ax),0);
					T(strip,T.LNG(ax),1) = T(strip,T.Dim(ax),1);
					T(strip,T.UNG(ax),0) = T(strip,1,0);
					T(strip,T.UNG(ax),1) = T(strip,1,1);
				}
			}
			Transpose(axis1,axis2,&T,-1);
		}
	}
}

void ComplexField::FFT()
{
	Field T;
	tw::grid::axis axis1,axis2;
	for (tw::Int ax=1;ax<=3;ax++)
	{
		axis1 = tw::grid::enumaxis(ax);
		axis2 = tw::grid::enumaxis(ax==3 ? 1 : ax+1);
		if (GlobalDim(ax)>1)
		{
			Transpose(axis1,axis2,&T,1);
			#pragma omp parallel
			{
				for (auto strip : StripRange(T,ax,strongbool::yes))
					ComplexFFT( &T(strip,1,0), &T(strip,1,1), T.Dim(ax), T.Stride(ax), 1.0 );
			}
			Transpose(axis1,axis2,&T,-1);
		}
	}
}

void ComplexField::InverseFFT()
{
	Field T;
	tw::grid::axis axis1,axis2;
	for (tw::Int ax=3;ax>=1;ax--)
	{
		axis1 = tw::grid::enumaxis(ax);
		axis2 = tw::grid::enumaxis(ax==3 ? 1 : ax+1);
		if (GlobalDim(ax)>1)
		{
			Transpose(axis1,axis2,&T,1);
			#pragma omp parallel
			{
				for (auto strip : StripRange(T,ax,strongbool::yes))
				{
					ComplexFFT( &T(strip,1,0), &T(strip,1,1), T.Dim(ax), T.Stride(ax), -1.0 );
					T(strip,T.LNG(ax),0) = T(strip,T.Dim(ax),0);
					T(strip,T.LNG(ax),1) = T(strip,T.Dim(ax),1);
					T(strip,T.UNG(ax),0) = T(strip,1,0);
					T(strip,T.UNG(ax),1) = T(strip,1,1);
				}
			}
			Transpose(axis1,axis2,&T,-1);
		}
	}
}

void Field::Hankel(const Element& e,tw::Int modes,std::valarray<tw::Float>& matrix)
{
	Field T;
	Transpose(e,tw::grid::x,tw::grid::z,&T,1);
	#pragma omp parallel
	{
		for (auto strip : StripRange(T,1,strongbool::yes))
			for (tw::Int c=e.low;c<=e.high;c++)
				Transform(&T(strip,1,c),T.Dim(1),modes,T.Stride(1),matrix);
	}
	Transpose(e,tw::grid::x,tw::grid::z,&T,-1);
}

void Field::InverseHankel(const Element& e,tw::Int modes,std::valarray<tw::Float>& matrix)
{
	Field T;
	Transpose(e,tw::grid::x,tw::grid::z,&T,1);
	#pragma omp parallel
	{
		for (auto strip : StripRange(T,1,strongbool::yes))
			for (tw::Int c=e.low;c<=e.high;c++)
				ReverseTransform(&T(strip,1,c),T.Dim(1),modes,T.Stride(1),matrix);
	}
	Transpose(e,tw::grid::x,tw::grid::z,&T,-1);
}
