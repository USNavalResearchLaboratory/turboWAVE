module;

#include "tw_includes.h"
#include "tw_logger.h"

export module fields:primitives;
import base;
import pic_primitives;
import static_space;
import logger;

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

/// Range for constraining temporal and internal dimensions while operating on space.
/// N.b. the bounds are of the [begin,end) form rather than the [low,high] used in the old `Element`.
export struct Rng
{
	tw::Int beg,end;
	Rng() {
		beg = 1;
		end = 2;
	}
	Rng(tw::Int beg,tw::Int end) {
		this->beg = beg;
		this->end = end;
	}
	Rng(tw::Int beg) {
		this->beg = beg;
		this->end = beg+1;
	}
	tw::Int Components() const {
		return end - beg;
	}
};


/// Range for constraining temporal and internal dimensions while operating on space
/// N.b. the bounds are of the [begin,end) form rather than the [low,high] used in the old `Element`.
export struct Rng04
{
	// index vector components from zero
	tw::Int b0,e0,b4,e4;
	Rng04() {
		b0 = 1;
		e0 = 2;
		b4 = 0;
		e4 = 1;
	}
	Rng04(tw::Int b0,tw::Int e0,tw::Int b4,tw::Int e4) {
		this->b0 = b0;
		this->e0 = e0;
		this->b4 = b4;
		this->e4 = e4;
	}
	Rng04(tw::Int b4) {
		this->b0 = 1;
		this->e0 = 2;
		this->b4 = b4;
		this->e4 = b4+1;
	}
	Rng04(Rng r) {
		b0 = 1;
		e0 = 2;
		b4 = r.beg;
		e4 = r.end;
	}
	tw::Int Times() const {
		return e0 - b0;
	}
	tw::Int Components() const {
		return e4 - b4;
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

	BoundaryCondition() {
		Reset(); // identity matrices
		sgn = 1; // default to low side
	}
	void Reset() {
		// Set all variables to the identity operation
		for (auto i=0;i<4;i++)
		{
			for (auto j=0;j<4;j++)
			{
				fold[i][j] = 0.0;
				force[i][j] = 0.0;
			}
			coeff[i] = 0.0;
			fold[i][i] = 1.0;
			force[i][i] = 1.0;
		}
	}
	void Set(tw::bc::fld theBoundaryCondition, tw::grid::side whichSide) {
		// In this routine, indexing is offset from StaticSpace indexing convention.
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


/// The Slice contains a copy of a subset of Field data, including its location in the Field.
/// Slice objects can be sent between nodes by passing the return values of `Buffer` and `BufferSize`
/// to various Send and Receive functions  The Field object has functions for moving data
/// to and from a Slice.
export template <class T>
struct Slice
{
	tw::Int beg[5], end[5]; // bounds of outer surface
private:
	tw::Int decodingStride[5];
	tw::Int encodingStride[5];
	std::valarray<T> data;
public:
	Slice() { ; }
	Slice(const tw::node5& beg, const tw::node5& end,bool zero=false) {
        Resize(beg, end);
        if (zero)
            data = 0.0;
    }
	void Resize(const tw::node5& beg, const tw::node5& end) {
		for (auto i=0; i<5; i++) {
			this->beg[i] = beg[i];
			this->end[i] = end[i];
		}
        // Decoding strides, no zero strides
        decodingStride[3] = 1;
        decodingStride[2] = decodingStride[3] * (end[3] - beg[3]);
        decodingStride[1] = decodingStride[2] * (end[2] - beg[2]);
		decodingStride[0] = decodingStride[1] * (end[1] - beg[1]);
		decodingStride[4] = decodingStride[0] * (end[0] - beg[0]);
        // Encoding strides, zero strides used
		for (auto i=0; i<5; i++) {
			encodingStride[i] = decodingStride[i] * tw::Int(beg[i]+1!=end[i]);
		}
        // Allocate space
        data.resize((end[0] - beg[0]) * (end[1] - beg[1]) * (end[2] - beg[2]) * (end[3] - beg[3]) * (end[4] - beg[4]));
    }
	T* Buffer()
	{
		return &data[0];
	}
	tw::Int BufferSize()
	{
		return sizeof(T) * (end[0] - beg[0]) * (end[1] - beg[1]) * (end[2] - beg[2]) * (end[3] - beg[3]) * (end[4] - beg[4]);
	}
	void Translate(const tw::grid::axis& axis, tw::Int displ) {
        tw::Int ax = tw::grid::naxis(axis);
        beg[ax] += displ;
        end[ax] += displ;
    }
	void Translate(tw::Int x, tw::Int y, tw::Int z) {
        beg[1] += x;
        end[1] += x;
        beg[2] += y;
        end[2] += y;
        beg[3] += z;
        end[3] += z;
    }
	T& operator () (const tw::Int& n, const tw::Int& i, const tw::Int& j, const tw::Int& k, const tw::Int& c)
	{
		const tw::Int idx = 
			(n - beg[0]) * encodingStride[0] +
			(i - beg[1]) * encodingStride[1] +
			(j - beg[2]) * encodingStride[2] +
			(k - beg[3]) * encodingStride[3] +
			(c - beg[4]) * encodingStride[4];
		#ifndef NDEBUG
		if (idx < 0 || idx > data.size()) {
			logger::ERROR(std::format("bad slice access {} >= {}",idx,data.size()));
			logger::DEBUG(std::format("coords {},{},{},{},{}",n,i,j,k,c));
			//logger::DEBUG(std::format("beg {} end {} strides {}",beg,end,encodingStride));
		}
		#endif
		return data[idx];
	}
	T operator () (const tw::Int& n, const tw::Int& i, const tw::Int& j, const tw::Int& k, const tw::Int& c) const
	{
		const tw::Int idx = 
			(n - beg[0]) * encodingStride[0] +
			(i - beg[1]) * encodingStride[1] +
			(j - beg[2]) * encodingStride[2] +
			(k - beg[3]) * encodingStride[3] +
			(c - beg[4]) * encodingStride[4];
		#ifndef NDEBUG
		if (idx < 0 || idx > data.size()) {
			logger::ERROR(std::format("bad slice access {} >= {}",idx,data.size()));
			logger::DEBUG(std::format("coords {},{},{},{},{}",n,i,j,k,c));
			//logger::DEBUG(std::format("beg {} end {} strides {}",beg,end,encodingStride));
		}
		#endif
		return data[idx];
	}
	Slice<T>& operator = (T a)
	{
		data = a;
		return *this;
	}
};
