module;

#include "tw_includes.h"
#include "tw_logger.h"

export module fields:base;
import :primitives;

import static_space;
import dyn_space;
import metric_space;
import tw_iterator;
import numerics;
import logger;

/// Field is a StaticSpace with data assigned to the cells, and operations on the data.
/// The data is some fixed number of floating point values per cell.
/// The storage pattern is variable, and can be specified by designating a packed axis.
/// The topology (dimensions, ghost cells) is inherited from the StaticSpace passed to the constructor.
/// It is possible to mix Fields with varying ghost cell layers, but it is MUCH SAFER
/// to keep the ghost cell layers the same for all Field instances in a calculation.
///
/// Note the Field does not inherit from DynSpace or MetricSpace, because these are
/// too variable and too heavyweight to carry around with every Field.
export struct Field: StaticSpace
{
protected:

	tw::Int totalCells;
	tw::Int boundaryCells;
	tw::Int ghostCells;
	Matrix<BoundaryCondition> bc0, bc1;
	tw::Int bufferState;

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
	void Initialize(const StaticSpace& ss, Task* task);
	void Initialize(const tw::Int& components,const StaticSpace& ss, Task* task);
	void SetBoundaryConditions(const Rng& rng, const tw::grid::axis& axis, tw::bc::fld low, tw::bc::fld high);
	friend void CopyBoundaryConditions(Field& dst, const Rng& r_dst, Field& src, const Rng& r_src);

	// ACCESSORS

	// Conventional accessors using fortran style array indexing
	tw::Float& operator () (const tw::Int& n, const tw::Int& i, const tw::Int& j, const tw::Int& k, const tw::Int& c)
	{
		return array[(n - lfg[0]) * encodingStride[0] + (i - lfg[1]) * encodingStride[1] + (j - lfg[2]) * encodingStride[2] + (k - lfg[3]) * encodingStride[3] + c * encodingStride[4]];
	}
	tw::Float operator () (const tw::Int& n, const tw::Int& i, const tw::Int& j, const tw::Int& k, const tw::Int& c) const
	{
		return array[(n - lfg[0]) * encodingStride[0] + (i - lfg[1]) * encodingStride[1] + (j - lfg[2]) * encodingStride[2] + (k - lfg[3]) * encodingStride[3] + c * encodingStride[4]];
	}
	// Accessors for stepping through spatial cells without concern for direction
	tw::Float& operator () (const tw::cell& cell, const tw::Int& c)
	{
		BOUNDS(cell.Index(c));
		return array[cell.Index(c)];
	}
	tw::Float operator () (const tw::cell& cell, const tw::Int& c) const
	{
		BOUNDS(cell.Index(c));
		return array[cell.Index(c)];
	}
	// Accessors for iterating across strips
	// Given a strip, returns component c in cell s as measured along the strip, t is the secondary visible axis, often time.
	tw::Float& operator () (const tw::strip& strip, const tw::Int& s, const tw::Int& c)
	{
		BOUNDS(strip.Index(s,c));
		return array[strip.Index(s,c)];
	}
	tw::Float operator () (const tw::strip& strip, const tw::Int& s, const tw::Int& c) const
	{
		BOUNDS(strip.Index(s,c));
		return array[strip.Index(s,c)];
	}
	// Accessors for promoting compiler vectorization
	tw::Float& operator () (const tw::xstrip<1>& v, const tw::Int& s, const tw::Int& c)
	{
		BOUNDS(v.Index1(s,c));
		return array[v.Index1(s,c)];
	}
	tw::Float operator () (const tw::xstrip<1>& v, const tw::Int& s, const tw::Int& c) const
	{
		BOUNDS(v.Index1(s,c));
		return array[v.Index1(s,c)];
	}
	tw::Float& operator () (const tw::xstrip<3>& v, const tw::Int& s, const tw::Int& c)
	{
		BOUNDS(v.Index1(s,c));
		return array[v.Index1(s,c)];
	}
	tw::Float operator () (const tw::xstrip<3>& v, const tw::Int& s, const tw::Int& c) const
	{
		BOUNDS(v.Index1(s,c));
		return array[v.Index1(s,c)];
	}
	// Centered Differencing
	tw::Float d1(const tw::cell& cell, const tw::Int& c, const tw::Int& ax) const
	{
		const tw::Int idx = cell.Index(c);
		return 0.5 * freq[ax] * (array[idx + encodingStride[ax]] - array[idx - encodingStride[ax]]);
	}
	tw::Float d1(const tw::xstrip<1>& v, const tw::Int& s, const tw::Int& c, const tw::Int& ax) const
	{
		const tw::Int idx = v.Index1(s,c);
		return 0.5 * freq[ax] * (array[idx + encodingStride[ax]] - array[idx - encodingStride[ax]]);
	}
	tw::Float d1(const tw::xstrip<3>& v, const tw::Int& s, const tw::Int& c, const tw::Int& ax) const
	{
		const tw::Int idx = v.Index1(s,c);
		return 0.5 * freq[ax] * (array[idx + encodingStride[ax]] - array[idx - encodingStride[ax]]);
	}
	tw::Float d2(const tw::strip& strip, const tw::Int& s, const tw::Int& c, const tw::Int& ax) const
	{
		const tw::Int idx = strip.Index(s,c);
		return freq[ax] * freq[ax] * (array[idx - encodingStride[ax]] - 2.0 * array[idx] + array[idx + encodingStride[ax]]);
	}
	tw::Float d2(const tw::xstrip<1>& v, const tw::Int& s, const tw::Int& c, const tw::Int& ax) const
	{
		const tw::Int idx = v.Index1(s,c);
		return freq[ax] * freq[ax] * (array[idx - encodingStride[ax]] - 2.0 * array[idx] + array[idx + encodingStride[ax]]);
	}
	tw::Float d2(const tw::xstrip<3>& v, const tw::Int& s, const tw::Int& c, const tw::Int& ax) const
	{
		const tw::Int idx = v.Index1(s,c);
		return freq[ax] * freq[ax] * (array[idx - encodingStride[ax]] - 2.0 * array[idx] + array[idx + encodingStride[ax]]);
	}
	// Unit directional offset
	tw::Float fwd(const tw::cell& cell, const tw::Int& c, const tw::Int& ax) const
	{
		return array[cell.Index(c) + encodingStride[ax]];
	}
	tw::Float bak(const tw::cell& cell, const tw::Int& c, const tw::Int& ax) const
	{
		return array[cell.Index(c) - encodingStride[ax]];
	}
	// Forward and Backward Differences and Sums
	tw::Float dfwd(const tw::strip& strip, const tw::Int& s, const tw::Int& c, const tw::Int& ax) const
	{
		tw::Int idx = strip.Index(s,c);
		return freq[ax] * (array[idx + encodingStride[ax]] - array[idx]);
	}
	tw::Float dfwd(const tw::xstrip<3>& v, const tw::Int& s, const tw::Int& c, const tw::Int& ax) const
	{
		tw::Int idx = v.Index1(s,c);
		return freq[ax] * (array[idx + encodingStride[ax]] - array[idx]);
	}
	tw::Float dbak(const tw::xstrip<3>& v, const tw::Int& s, const tw::Int& c, const tw::Int& ax) const
	{
		tw::Int idx = v.Index1(s,c);
		return freq[ax] * (array[idx] - array[idx - encodingStride[ax]]);
	}
	tw::Float sfwd(const tw::xstrip<1>& v, const tw::Int& s, const tw::Int& c, const tw::Int& ax) const
	{
		tw::Int idx = v.Index1(s,c);
		return 0.5 * (array[idx + encodingStride[ax]] + array[idx]);
	}
	tw::Float sbak(const tw::xstrip<1>& v, const tw::Int& s, const tw::Int& c, const tw::Int& ax) const
	{
		tw::Int idx = v.Index1(s,c);
		return 0.5 * (array[idx] + array[idx - encodingStride[ax]]);
	}
	tw::Float sfwd(const tw::xstrip<3>& v, const tw::Int& s, const tw::Int& c, const tw::Int& ax) const
	{
		tw::Int idx = v.Index1(s,c);
		return 0.5 * (array[idx + encodingStride[ax]] + array[idx]);
	}
	tw::Float sbak(const tw::xstrip<3>& v, const tw::Int& s, const tw::Int& c, const tw::Int& ax) const
	{
		tw::Int idx = v.Index1(s,c);
		return 0.5 * (array[idx] + array[idx - encodingStride[ax]]);
	}
	friend Rng All(const Field& f) { return Rng(0,f.num[4]); }
	tw::Int Components() const { return num[4]; }
	tw::Int TotalCells() const { return totalCells; }
	size_t TotalBytes() const { return array.size() * sizeof(tw::Float); }
	size_t BoundaryBytes() const { return boundaryData.size() * sizeof(tw::Float); }
	size_t GhostBytes() const { return ghostData.size() * sizeof(tw::Float); }
	size_t BoundaryMapBytes() const { return boundaryDataIndexMap.size() * sizeof(tw::Int); }
	size_t GhostMapBytes() const { return ghostDataIndexMap.size() * sizeof(tw::Int); }

	// GPU Support

#ifdef USE_OPENCL
	void InitializeComputeBuffer();
	void SendToComputeBuffer();
	void ReceiveFromComputeBuffer();
	void SendBoundaryCellsToComputeBuffer();
	void ReceiveBoundaryCellsFromComputeBuffer();
	tw::Float CellValueInComputeBuffer(tw::Int i, tw::Int j, tw::Int k, tw::Int c); // an expensive call for one data point (useful for last step of reduction)
	void UpdateGhostCellsInComputeBuffer(const Rng04& r);
	void SendGhostCellsToComputeBuffer();
	void ZeroGhostCellsInComputeBuffer();
	friend void SwapComputeBuffers(Field& f1, Field& f2);
	friend void CopyComputeBuffer(Field& dst, Field& src);
	friend void CopyComplexComputeBufferMod2(Field& dst, Field& src);
	void FillComputeBufferVec4(const Rng04& r, tw::vec4& A);
	tw::Float DestructiveSumComputeBuffer();
	tw::Float DestructiveNorm1ComputeBuffer();
	void WeightComputeBufferByVolume(MetricSpace& ms, tw::Float inv);
	void DestructiveComplexMod2ComputeBuffer();
	void MADDComputeBuffer(tw::Float m, tw::Float a);
#endif

	// Transformation

	void MultiplyCellVolume(const Rng04& r,const MetricSpace& ms);
	void DivideCellVolume(const Rng04& r,const MetricSpace& ms);
	friend void Swap(Field& f1, Field& f2);
	void Swap(const Rng04& r1, const Rng04& r2);
	friend void CopyFieldData(Field& dst, const Rng04& r_dst, Field& src, const Rng04& r_src);
	friend void CopyGhostCellData(Field& dst, const Rng04& r_dst, Field& src, const Rng04& r_src);
	friend void AddFieldData(Field& dst, const Rng04& r_dst, Field& src, const Rng04& r_src);
	friend void AddMulFieldData(Field& dst, const Rng04& r_dst, Field& src, const Rng04& r_src, tw::Float mul);
	void SmoothingPass(const Rng04& r, tw::Int ax, const MetricSpace& ms, const tw::Float& X0, const tw::Float& X1, const tw::Float& X2);
	void Smooth(const Rng04& r, const MetricSpace& ms, tw::Int smoothPasses[4], tw::Int compPasses[4]);
	void Shift(const Rng& r, const tw::strip& s, tw::Int cells, const tw::Float* incoming);
	void Shift(const Rng& r, const tw::strip& s, tw::Int cells, const tw::Float& incoming);

	tw::Float RealAxialEigenvalue(tw::Int z,const DynSpace& ds);
	tw::Float RealEigenvalue(tw::Int x, tw::Int y,const DynSpace& ds);
	tw::Float RealCyclicEigenvalue(tw::Int x, tw::Int y,const DynSpace& ds);
	tw::Float RealCyclicEigenvalue(tw::Int x, tw::Int y, tw::Int z,const DynSpace& ds);
	void RealAxialSineTransform(const Rng04& r,const DynSpace& ds);
	void RealInverseAxialSineTransform(const Rng04& r,const DynSpace& ds);
	void RealTransverseSineTransform(const Rng04& r,const DynSpace& ds);
	void RealInverseTransverseSineTransform(const Rng04& r,const DynSpace& ds);
	void RealTransverseCosineTransform(const Rng04& r,const DynSpace& ds);
	void RealInverseTransverseCosineTransform(const Rng04& r,const DynSpace& ds);
	void RealTransverseFFT(const Rng04& r,const DynSpace& ds);
	void RealInverseTransverseFFT(const Rng04& r,const DynSpace& ds);

	tw::Float ComplexCyclicEigenvalue(tw::Int x, tw::Int y,const DynSpace& ds);
	tw::Float ComplexCyclicEigenvalue(tw::Int x, tw::Int y, tw::Int z,const DynSpace& ds);
	void ComplexFFT(const Rng04& r,const DynSpace& ds);
	void ComplexInverseFFT(const Rng04& r,const DynSpace& ds);
	void ComplexTransverseFFT(const Rng04& r,const DynSpace& ds);
	void ComplexInverseTransverseFFT(const Rng04& r,const DynSpace& ds);

	void Hankel(const Rng04& r, tw::Int modes, std::valarray<tw::Float>& matrix);
	void InverseHankel(const Rng04& r, tw::Int modes, std::valarray<tw::Float>& matrix);

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
	const tw::vec3 Vec3(const tw::Int& n, const tw::Int& i, const tw::Int& j, const tw::Int& k, const tw::Int& c) const
	{
		return tw::vec3((*this)(n, i, j, k, c), (*this)(n, i, j, k, c + 1), (*this)(n, i, j, k, c + 2));
	}
	const tw::vec3 Vec3(const tw::cell& cell, const tw::Int& c) const
	{
		return tw::vec3((*this)(cell, c), (*this)(cell, c + 1), (*this)(cell, c + 2));
	}
	const tw::vec3 Vec3(const tw::strip& strip, const tw::Int& s, const tw::Int& c) const
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

	// Boundaries and messages

	void ApplyFoldingCondition(const Rng04& r);
	void ApplyBoundaryCondition(const Rng04& r, bool homogeneous = true);
	template <class T>
	void AdjustTridiagonalForBoundaries(const Rng& r, const tw::grid::axis& axis, const tw::grid::side& side, std::valarray<T>& T1, std::valarray<T>& T2, std::valarray<T>& T3, std::valarray<T>& source, T val);
	void ZeroGhostCells(const Rng04& r);

	void StripCopyProtocol(tw::Int axis, tw::Int shift, Slice<tw::Float>* planeIn, Slice<tw::Float>* planeOut, bool add);
	void DownwardCopy(const Rng04& r, const tw::grid::axis& axis, tw::Int cells);
	void UpwardCopy(const Rng04& r, const tw::grid::axis& axis, tw::Int cells);
	void DownwardDeposit(const Rng04& r, const tw::grid::axis& axis, tw::Int cells);
	void UpwardDeposit(const Rng04& r, const tw::grid::axis& axis, tw::Int cells);

	void CopyFromNeighbors(const Rng04& r);
	void DepositFromNeighbors(const Rng04& r);

	Slice<tw::Float>* FormTransposeBlock(const Rng04& r, const tw::grid::axis& axis1, const tw::grid::axis& axis2, tw::Int start1, tw::Int stop1, tw::Int start2, tw::Int stop2);
	void Transpose(const Rng04& r, const tw::grid::axis& axis1, const tw::grid::axis& axis2, Field* target, tw::Int inversion);

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

	void Interpolate(const Rng04& r, std::valarray<tw::Float>& val, const weights_3D& weights) const;
	void InterpolateOnto(const Rng04& r, std::valarray<tw::Float>& val, const weights_3D& weights);
};

Field::Field()
{
	totalCells = 0;
	boundaryCells = 0;
	ghostCells = 0;
	bufferState = 0;
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

/// @brief initialize the field with a StaticSpace, but with a new internal dimension
/// @param components update the internal dimension relative to `ss`
/// @param ss static space to use for this Field, data is copied
/// @param task pointer to the concurrent task with this Field's domain
void Field::Initialize(const tw::Int& components,const StaticSpace& ss,Task *task) {
	auto new_ss = StaticSpace(components,ss);
	Initialize(new_ss,task);
}

/// @brief initialize the field with a StaticSpace
/// @param ss static space to use for this Field, data is copied
/// @param task pointer to the concurrent task with this Field's domain
void Field::Initialize(const StaticSpace& ss,Task *task) {
	StaticSpace::operator = (ss);
	//logger::TRACE(std::format("field init dim {}",dim));
	//logger::TRACE(std::format("field init estride {}",encodingStride));
	this->task = task;
	totalCells = num[0]*num[1]*num[2]*num[3];
	bc0.resize(4,num[4]);
	bc1.resize(4,num[4]);
	for (tw::Int i=0;i<4;i++)
		for (tw::Int c=0;c<num[4];c++)
		{
			bc0(i,c).Set(tw::bc::fld::none,tw::grid::low);
			bc1(i,c).Set(tw::bc::fld::none,tw::grid::high);
		}
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
	array.resize(totalCells*num[4]);
	boundaryData.resize(boundaryCells*num[4]);
	ghostData.resize(ghostCells*num[4]);
	boundaryDataIndexMap.resize(boundaryCells*num[4]);
	ghostDataIndexMap.resize(ghostCells*num[4]);

	// Set up the index maps going from coordinate space to memory space.
	// Ordering is not important as long as we gather all the right cells.

	tw::Int bOffset=0,gOffset=0;
	for (tw::Int c=0;c<num[4];c++)
		for (tw::Int ax=1;ax<=3;ax++)
		{
			if (dim[ax]>1)
				for (auto s : StripRange(*this,ax,0,0,strongbool::yes))
				{
					for (tw::Int i=0;i<2*layers[ax];i++)
					{
						boundaryDataIndexMap[bOffset+i] = s.Index(lfg[ax]+i,c);
						boundaryDataIndexMap[bOffset+2*layers[ax]+i] = s.Index(ufg[ax]-2*layers[ax]+i+1,c);
					}
					for (tw::Int i=0;i<layers[ax];i++)
					{
						ghostDataIndexMap[gOffset+i] = s.Index(lfg[ax]+i,c);
						ghostDataIndexMap[gOffset+layers[ax]+i] = s.Index(ufg[ax]-layers[ax]+i+1,c);
					}
					bOffset += 4*layers[ax];
					gOffset += 2*layers[ax];
				}
		}
}

void Field::ReadCheckpoint(std::ifstream& inFile)
{
	StaticSpace::ReadCheckpoint(inFile);
	inFile.read((char*)&bc0(0,0),sizeof(BoundaryCondition)*num[0]*4);
	inFile.read((char*)&bc1(0,0),sizeof(BoundaryCondition)*num[0]*4);
	inFile.read((char *)&array[0],sizeof(tw::Float)*totalCells*num[0]);
}

void Field::WriteCheckpoint(std::ofstream& outFile)
{
	StaticSpace::WriteCheckpoint(outFile);
	outFile.write((char*)&bc0(0,0),sizeof(BoundaryCondition)*num[0]*4);
	outFile.write((char*)&bc1(0,0),sizeof(BoundaryCondition)*num[0]*4);
	outFile.write((char *)&array[0],sizeof(tw::Float)*totalCells*num[0]);
}

void Field::MultiplyCellVolume(const Rng04& r, const MetricSpace& m)
{
	for (auto n=r.b0; n<r.e0; n++)
		for (auto cell : EntireCellRange(*this,n))
			for (auto c=r.b4; c<r.e4; c++)
				(*this)(cell,c) *= m.dS(cell,0);
}

void Field::DivideCellVolume(const Rng04& r, const MetricSpace& m)
{
	for (auto n=r.b0; n<r.e0; n++)
		for (auto cell : EntireCellRange(*this,n))
			for (auto c=r.b4; c<r.e4; c++)
				(*this)(cell,c) /= m.dS(cell,0);
}

void Field::Shift(const Rng& r,const tw::strip& s,tw::Int cells,const tw::Float* incoming)
{
	// Propagate field pattern to the right if <cells> positive, to the left if negative
	// argument <incoming> is value to inject into left or right ghost cell

	tw::Int i,c,ax=s.Axis();

	if (cells>0)
	{
		for (i=UFG(ax);i>=LFG(ax)+cells;i--)
			for (c=r.beg; c<r.end; c++)
				(*this)(s,i,c) = (*this)(s,i-cells,c);
		for (i=LFG(ax);i<LFG(ax)+cells;i++)
			for (c=r.beg; c<r.end; c++)
				(*this)(s,i,c) = incoming[c-r.beg];
	}
	if (cells<0)
	{
		for (i=LFG(ax);i<=UFG(ax)+cells;i++)
			for (c=r.beg; c<r.end; c++)
				(*this)(s,i,c) = (*this)(s,i-cells,c);
		for (i=UFG(ax)+1+cells;i<=UFG(ax);i++)
			for (c=r.beg; c<r.end; c++)
				(*this)(s,i,c) = incoming[c-r.beg];
	}
}

void Field::Shift(const Rng& r,const tw::strip& s,tw::Int cells,const tw::Float& incoming)
{
	// Propagate field pattern to the right if <cells> positive, to the left if negative
	// argument <incoming> is value to inject into left or right ghost cell

	tw::Int i,c,ax=s.Axis();

	if (cells>0)
	{
		for (i=UFG(ax);i>=LFG(ax)+cells;i--)
			for (c=r.beg; c<r.end; c++)
				(*this)(s,i,c) = (*this)(s,i-cells,c);
		for (i=LFG(ax);i<LFG(ax)+cells;i++)
			for (c=r.beg; c<r.end; c++)
				(*this)(s,i,c) = incoming;
	}
	if (cells<0)
	{
		for (i=LFG(ax);i<=UFG(ax)+cells;i++)
			for (c=r.beg; c<r.end; c++)
				(*this)(s,i,c) = (*this)(s,i-cells,c);
		for (i=UFG(ax)+1+cells;i<=UFG(ax);i++)
			for (c=r.beg; c<r.end; c++)
				(*this)(s,i,c) = incoming;
	}
}

//////////////////////////////////
//                              //
//  Global Boundary Conditions  //
//                              //
//////////////////////////////////

template<class T>
void Field::AdjustTridiagonalForBoundaries(const Rng& r, const tw::grid::axis& axis, const tw::grid::side& side, std::valarray<T>& T1, std::valarray<T>& T2, std::valarray<T>& T3, std::valarray<T>& source, T val)
{
	// Modify the tridiagonal matrix to respect the boundary condition cell_1 = force_1j * cell_j + coeff_1 * val
	// force_ij and coeff_i are defined by the boundary condition object.
	// Since val may be passed in based on values in the outermost ghost cell, we assume force_10 = 0.
	// N.b. indexing of force_ij and coeff_i starts with outermost ghost cell = 0 and works inward, regardless of tw::grid::side.
	// The element is only used to define the boundary conditions, no field data is used.
	const tw::Int ax = tw::grid::naxis(axis);
	const tw::Int N = Dim(ax) - 1;
	const tw::Int c = r.beg; // component determining the BC
	assert(r.beg == r.end+1); // asking for more than one BC is not defined
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

void Field::SetBoundaryConditions(const Rng& r,const tw::grid::axis& axis,tw::bc::fld low,tw::bc::fld high)
{
	tw::Int ax = tw::grid::naxis(axis);

	for (auto i=r.beg; i<r.end; i++)
	{
		if (task->n0[ax]==MPI_PROC_NULL && dim[ax]>1)
			bc0(ax,i).Set(low,tw::grid::low);
		if (task->n1[ax]==MPI_PROC_NULL && dim[ax]>1)
			bc1(ax,i).Set(high,tw::grid::high);
	}
}

export void CopyBoundaryConditions(Field& dst,const Rng& r_dst,Field& src,const Rng& r_src)
{
	if (r_src.Components() != r_dst.Components()) {
		throw tw::FatalError("bad copy operation (boundary conditions)");
	}
	for (auto c=0; c<r_src.Components(); c++) {
		for (auto ax=0;ax<4;ax++) {
			dst.bc0(ax,r_dst.beg + c) = src.bc0(ax,r_src.beg + c);
			dst.bc1(ax,r_dst.beg + c) = src.bc1(ax,r_src.beg + c);
		}
	}
}

void Field::ZeroGhostCells(const Rng04& r)
{
	for (auto n=r.b0; n<r.e0; n++)
		for (auto c=r.b4; c<r.e4; c++)
			for (auto ax=1; ax<4; ax++)
				for (auto strip : StripRange(*this,ax,0,n,strongbool::yes))
					for (auto s=0; s<layers[ax]; s++)
					{
						(*this)(strip,lfg[ax]+s,c) = 0.0;
						(*this)(strip,dim[ax]+s+1,c) = 0.0;
					}
}

void Field::ApplyFoldingCondition(const Rng04& r)
{
	for (auto n=r.b0; n<r.e0; n++)
		for (auto c=r.b4; c<r.e4; c++)
			for (auto ax=1; ax<4; ax++)
				if (num[ax]>1)
					for (auto strip : StripRange(*this,ax,0,n,strongbool::yes))
					{
						bc0(ax,c).FoldingOperation(&(*this)(strip,lfg[ax],c),encodingStride[ax]);
						bc1(ax,c).FoldingOperation(&(*this)(strip,ufg[ax],c),encodingStride[ax]);
					}
}

void Field::ApplyBoundaryCondition(const Rng04& r,bool homogeneous)
{
	tw::Float mult = homogeneous ? 0.0 : 1.0;
	for (auto n=r.b0; n<r.e0; n++)
		for (auto c=r.b4; c<r.e4; c++)
			for (auto ax=1; ax<4; ax++)
				if (num[ax]>1)
					for (auto strip : StripRange(*this,ax,0,n,strongbool::yes))
					{
						bc0(ax,c).ForcingOperation(&(*this)(strip,lfg[ax],c),encodingStride[ax],(*this)(strip,lfg[ax],c) * mult);
						bc1(ax,c).ForcingOperation(&(*this)(strip,ufg[ax],c),encodingStride[ax],(*this)(strip,ufg[ax],c) * mult);
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

void Field::UpdateGhostCellsInComputeBuffer(const Rng04& r)
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

void Field::FillComputeBufferVec4(const Rng04& r,tw::vec4& A)
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

void Field::DownwardCopy(const Rng04& r,const tw::grid::axis& axis,tw::Int cells)
{
	Slice<tw::Float> *planeIn,*planeOut;
	const tw::Int ax = tw::grid::naxis(axis);
	tw::node5 beg,end;

	if (dim[ax]==1) return;
	beg[0] = r.b0; end[0] = r.e0;
	beg[4] = r.b4; end[4] = r.e4;
	for (tw::Int i=1;i<=3;i++) {
		beg[i] = ax==i ? ung[i] : lfg[i];
		end[i] = ax==i ? ung[i] + cells : ufg[i] + 1;
	}
	planeIn = new Slice<tw::Float>(beg,end);
	planeOut = new Slice<tw::Float>(beg,end);
	planeOut->Translate(axis,lng[ax] + 1 - ung[ax]);

	StripCopyProtocol(tw::grid::naxis(axis),-1,planeIn,planeOut,false);

	delete planeOut;
	delete planeIn;
}

void Field::UpwardCopy(const Rng04& r,const tw::grid::axis& axis,tw::Int cells)
{
	Slice<tw::Float> *planeIn,*planeOut;
	const tw::Int ax = tw::grid::naxis(axis);
	tw::node5 beg,end;

	if (dim[ax]==1) return;
	beg[0] = r.b0; end[0] = r.e0;
	beg[4] = r.b4; end[4] = r.e4;
	for (tw::Int i=1;i<=3;i++)
	{
		beg[i] = ax==i ? lng[i] - cells + 1 : lfg[i];
		end[i] = ax==i ? lng[i] + 1 : ufg[i] + 1;
	}
	planeIn = new Slice<tw::Float>(beg,end);
	planeOut = new Slice<tw::Float>(beg,end);
	planeOut->Translate(axis,ung[ax] - lng[ax] - 1);

	StripCopyProtocol(tw::grid::naxis(axis),1,planeIn,planeOut,false);

	delete planeOut;
	delete planeIn;
}

void Field::DownwardDeposit(const Rng04& r,const tw::grid::axis& axis,tw::Int cells)
{
	// deposits work by sending both ghost cells and edge cells one way and adding.
	// Before accepting the return message these cells are zeroed.
	// The return message then effectively overwrites them with the correct data.
	Slice<tw::Float> *planeIn,*planeOut;
	const tw::Int ax = tw::grid::naxis(axis);
	tw::node5 beg,end;

	if (dim[ax]==1) return;
	beg[0] = r.b0; end[0] = r.e0;
	beg[4] = r.b4; end[4] = r.e4;
	for (tw::Int i=1;i<=3;i++)
	{
		beg[i] = ax==i ? ung[i] - cells : lfg[i];
		end[i] = ax==i ? ung[i] + cells : ufg[i] + 1;
	}
	planeIn = new Slice<tw::Float>(beg,end);
	planeOut = new Slice<tw::Float>(beg,end);
	planeOut->Translate(axis,lng[ax] + 1 - ung[ax]);

	StripCopyProtocol(tw::grid::naxis(axis),-1,planeIn,planeOut,true);

	delete planeOut;
	delete planeIn;
}

void Field::UpwardDeposit(const Rng04& r,const tw::grid::axis& axis,tw::Int cells)
{
	// deposits work by sending both ghost cells and edge cells one way and adding.
	// Before accepting the return message these cells are zeroed.
	// The return message then effectively overwrites them with the correct data.
	Slice<tw::Float> *planeIn,*planeOut;
	const tw::Int ax = tw::grid::naxis(axis);
	tw::node5 beg,end;

	if (dim[ax]==1) return;
	beg[0] = r.b0; end[0] = r.e0;
	beg[4] = r.b4; end[4] = r.e4;
	for (tw::Int i=1;i<=3;i++)
	{
		beg[i] = ax==i ? lng[i] + 1 - cells : lfg[i];
		end[i] = ax==i ? lng[i] + cells + 1 : ufg[i] + 1;
	}
	planeIn = new Slice<tw::Float>(beg,end);
	planeOut = new Slice<tw::Float>(beg,end);
	planeOut->Translate(axis,ung[ax] - lng[ax] - 1);

	StripCopyProtocol(tw::grid::naxis(axis),1,planeIn,planeOut,true);

	delete planeOut;
	delete planeIn;
}

void Field::CopyFromNeighbors(const Rng04& r)
{
	DownwardCopy(r,tw::grid::x,1);
	UpwardCopy(r,tw::grid::x,1);

	DownwardCopy(r,tw::grid::y,1);
	UpwardCopy(r,tw::grid::y,1);

	DownwardCopy(r,tw::grid::z,1);
	UpwardCopy(r,tw::grid::z,1);
}

void Field::DepositFromNeighbors(const Rng04& r)
{
	tw::node5 beg,end;
	beg[0] = r.b0; end[0] = r.e0;
	beg[4] = r.b4; end[4] = r.e4;
	for (tw::Int ax=1;ax<=3;ax++)
		if (dim[ax]>1)
		{
			DownwardDeposit(r,tw::grid::enumaxis(ax),layers[ax]);
			if (task->n0[ax]!=MPI_PROC_NULL)
			{
				for (auto i=1; i<=3; i++) {
					beg[i] = lfg[i];
					end[i] = ax==i ? lfg[ax] + 2*layers[ax] : ufg[i] + 1;
				}
				Slice<tw::Float>* plane = new Slice<tw::Float>(beg,end);
				ZeroDataInField<tw::Float>(plane);
				delete plane;
			}
			UpwardDeposit(r,tw::grid::enumaxis(ax),layers[ax]);
		}
}


///////////////////////////
//                       //
// SIMPLE TRANSFORMATION //
//                       //
///////////////////////////


void Field::SmoothingPass(const Rng04& r,tw::Int ax,const MetricSpace& ms,const tw::Float& X0,const tw::Float& X1,const tw::Float& X2)
{
	tw::Int i,c;
	tw::Float ansNow,temp;
	tw::grid::axis axs[4] = { tw::grid::t, tw::grid::x, tw::grid::y, tw::grid::z };

	for (auto n=r.b0; n<r.e0; n++)
		for (auto c=r.b4; c<=r.e4; c++)
			if (dim[ax]>1)
				for (auto s : StripRange(*this,ax,0,n,strongbool::yes))
				{
					temp = (*this)(s,0,c);
					for (i=1;i<=dim[ax];i++)
					{
						ansNow = X0*temp + X1*(*this)(s,i,c) + X2*(*this)(s,i+1,c);
						temp = (*this)(s,i,c);
						(*this)(s,i,c) = ansNow;
					}
					bc0(ax,c).ForcingOperation(&(*this)(s,lfg[ax],c),encodingStride[ax],0.0);
					bc1(ax,c).ForcingOperation(&(*this)(s,ufg[ax],c),encodingStride[ax],0.0);
				}

	UpwardCopy(r,axs[ax],1);
	DownwardCopy(r,axs[ax],1);
}

void Field::Smooth(const Rng04& r,const MetricSpace& ms,tw::Int smoothPasses[4],tw::Int compPasses[4])
{
	for (tw::Int ax=1;ax<=3;ax++)
	{
		for (tw::Int ipass=0;ipass<smoothPasses[ax];ipass++)
			SmoothingPass(r,ax,ms,0.25,0.5,0.25);
		for (tw::Int ipass=0;ipass<compPasses[ax];ipass++)
			SmoothingPass(r,ax,ms,-1.25,3.5,-1.25);
	}
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

void Field::Swap(const Rng04& r1,const Rng04& r2)
{
	if (r1.Components() != r2.Components() || r1.Times() != r2.Times() || r1.b0 != r2.b0) {
		throw tw::FatalError("bad swap operation");
	}
	for (auto c=0; c<r1.Components(); c++) {
		for (auto n = 0; n<r1.Times(); n++) {
			for(auto cell : EntireCellRange(*this,n))
			{
				tw::Float temp = (*this)(cell,r1.b4+c);
				(*this)(cell,r1.b4+c) = (*this)(cell,r2.b4+c);
				(*this)(cell,r2.b4+c) = temp;
			}
		}
	}
}

export void CopyFieldData(Field& dst,const Rng04& r_dst,Field& src,const Rng04& r_src)
{
	if (r_src.Components() != r_dst.Components() || r_src.Times() != r_dst.Times() || r_src.b0 != r_dst.b0) {
		throw tw::FatalError("bad copy operation");
	}
	for (auto n=0;n<r_dst.Times();n++)
		for (auto c=0;c<r_dst.Components();c++)
			for (auto cell : EntireCellRange(dst,n))
				dst(cell,r_dst.b4+c) = src(cell,r_src.b4+c);
}

export void CopyGhostCellData(Field& dst,const Rng04& r_dst,Field& src,const Rng04& r_src)
{
	if (r_src.Components() != r_dst.Components() || r_src.Times() != r_dst.Times() || r_src.b0 != r_dst.b0) {
		throw tw::FatalError("bad copy operation");
	}
	for (auto n=0;n<r_dst.Times();n++)
		for (auto c=0;c<r_dst.Components();c++)
			for (auto ax=1;ax<=3;ax++)
				for (auto strip : StripRange(dst,ax,0,n,strongbool::yes))
					for (tw::Int i=0;i<dst.layers[ax];i++)
					{
						dst(strip,dst.lfg[ax]+i,r_dst.b4+c) = src(strip,src.lfg[ax]+i,r_src.b4+c);
						dst(strip,dst.ufg[ax]-i,r_dst.b4+c) = src(strip,src.ufg[ax]-i,r_src.b4+c);
					}
}

export void AddFieldData(Field& dst,const Rng04& r_dst,Field& src,const Rng04& r_src)
{
	if (r_src.Components() != r_dst.Components() || r_src.Times() != r_dst.Times() || r_src.b0 != r_dst.b0) {
		throw tw::FatalError("bad add operation");
	}
	for (auto n=0;n<r_dst.Times();n++)
		for (auto c=0;c<r_dst.Components();c++)
			for (auto cell : EntireCellRange(dst,n))
				dst(cell,r_dst.b4+c) += src(cell,r_src.b4+c);
}

export void AddMulFieldData(Field& dst,const Rng04& r_dst,Field& src,const Rng04& r_src,tw::Float mul)
{
	if (r_src.Components() != r_dst.Components() || r_src.Times() != r_dst.Times() || r_src.b0 != r_dst.b0) {
		throw tw::FatalError("bad add operation");
	}
	for (auto n=0;n<r_dst.Times();n++)
		for (auto c=0;c<r_dst.Components();c++)
			for (auto cell : EntireCellRange(dst,n))
				dst(cell,r_dst.b4+c) += mul*src(cell,r_src.b4+c);
}


////////////////////////
//                    //
//  GATHER / SCATTER  //
//                    //
////////////////////////


inline void Field::Interpolate(const Rng04& r, std::valarray<tw::Float>& val, const weights_3D& weights) const
{
	tw::Int n = r.b0;
	tw::Int topo[4];
	tw::Float wijk;
	DecodeCell(weights.cell, topo);
	val = 0.0;

	for (auto c = r.b4; c < r.e4; c++)
		for (auto i = 0; i < 3; i++)
			for (auto j = 0; j < 3; j++)
				for (auto k = 0; k < 3; k++)
				{
					wijk = weights.w[i][0] * weights.w[j][1] * weights.w[k][2];
					val[c] += wijk * (*this)(n, topo[1] + i - 1, topo[2] + j - 1, topo[3] + k - 1, c);
				}
}

inline void Field::InterpolateOnto(const Rng04& r, std::valarray<tw::Float>& val, const weights_3D& weights)
{
	tw::Int n = r.b0;
	tw::Int topo[4];
	tw::Float wijk;
	DecodeCell(weights.cell, topo);

	for (auto c = r.b4; c < r.e4; c++)
		for (auto i = 0; i < 3; i++)
			for (auto j = 0; j < 3; j++)
				for (auto k = 0; k < 3; k++)
				{
					wijk = weights.w[i][0] * weights.w[j][1] * weights.w[k][2];
#pragma omp atomic update
					(*this)(n, topo[1] + i - 1, topo[2] + j - 1, topo[3] + k - 1, c) += wijk * val[c];
				}
}


template <class T>
void Field::LoadDataIntoSlice(Slice<T>* v)
{
	for (auto c = v->beg[4]; c < v->end[4]; c++)
		for (auto n = v->beg[0]; n < v->end[0]; n++)
			for (auto i = v->beg[1]; i < v->end[1]; i++)
				for (auto j = v->beg[2]; j < v->end[2]; j++)
					for (auto k = v->beg[3]; k < v->end[3]; k++)
						(*v)(n, i, j, k, c) = (*this)(n, i, j, k, c);
}

template <class T>
void Field::SaveDataFromSlice(Slice<T>* v)
{
	for (auto c = v->beg[4]; c < v->end[4]; c++)
		for (auto n = v->beg[0]; n < v->end[0]; n++)
			for (auto i = v->beg[1]; i < v->end[1]; i++)
				for (auto j = v->beg[2]; j < v->end[2]; j++)
					for (auto k = v->beg[3]; k < v->end[3]; k++)
						(*this)(n, i, j, k, c) = (*v)(n, i, j, k, c);
}

template <class T>
void Field::ZeroDataInField(Slice<T>* v)
{
	for (auto c = v->beg[4]; c < v->end[4]; c++)
		for (auto n = v->beg[0]; n < v->end[0]; n++)
			for (auto i = v->beg[1]; i < v->end[1]; i++)
				for (auto j = v->beg[2]; j < v->end[2]; j++)
					for (auto k = v->beg[3]; k < v->end[3]; k++)
						(*this)(n, i, j, k, c) = 0.0;
}

template <class T>
void Field::AddDataFromSlice(Slice<T>* v)
{
	for (auto c = v->beg[4]; c < v->end[4]; c++)
		for (auto n = v->beg[0]; n < v->end[0]; n++)
			for (auto i = v->beg[1]; i < v->end[1]; i++)
				for (auto j = v->beg[2]; j < v->end[2]; j++)
					for (auto k = v->beg[3]; k < v->end[3]; k++)
						(*this)(n, i, j, k, c) += (*v)(n, i, j, k, c);
}

template <class T>
void Field::AddDataFromSliceAtomic(Slice<T>* v)
{
	for (auto c = v->beg[4]; c < v->end[4]; c++)
		for (auto n = v->beg[0]; n < v->end[0]; n++)
			for (auto i = v->beg[1]; i < v->end[1]; i++)
				for (auto j = v->beg[2]; j < v->end[2]; j++)
					for (auto k = v->beg[3]; k < v->end[3]; k++)
#pragma omp atomic update
						(*this)(n, i, j, k, c) += (*v)(n, i, j, k, c);
}
