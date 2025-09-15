module;

#include "tw_includes.h"

export module tw_iterator;
import base;
import static_space;

// TurboWAVE iterators serve the specific purpose of accessing the Field class.
// They are not intended for generalized use.  They imitate the standard library, but there are differences.
// They allow for range based loops, automatic parallelism, and compact strip processing.

// There are 3 types of classes involved:
// 1. Range classes : Representation of cells or strips to iterate over
// 2. tw::iterator : used under the hood to move through the range, and support range based loops
// 3. Reference classes : used as arguments in Field accessors, much like an array index

// If you are used to standard containers, compare as follows:
// 1. Range classes behave like a container class, e.g., have begin() and end() methods.  However,
// they have no storage of their own, they are a proxy for Field instances.
// 2. tw::iterator plays essentially the same role as std::iterator, except that dereferencing it
// does not give an element of Field, but rather a reference that can be passed to any Field object.
// 3. Reference classes essentially point to some part of the data represented by Field objects.  They are
// what results from dereferencing tw::iterator.

// Range classes are aware of OpenMP state and automatically split the work accross threads.
// Range classes are created from a StaticSpace, and therefore may only be used in Field objects created from
// the same StaticSpace.

export namespace tw
{
	template <class RNG,class REF>
	class iterator
	{
		RNG *range;
		// The RNG class is responsible for decoding the count.
		// The decoding produces a tuple that Field can use to retrieve floating point data.
		// The decoded indices are encapsulated in the REF class.
		tw::Int curr,first,last;
		// curr is the encoded index of the current object, called the count.
		// The global encoding is such that count ranges from (0...N-1), where N is the number of objects.
		// curr contains the GLOBAL count, which the user does not interact with.
		// count() returns the LOCAL count, which ranges from (0...M-1), where M is the number of objects to be processed by the current thread.
		// first and last are defined such that first <= curr <= last.

	public:
		iterator(RNG *r,tw::Int c,tw::Int f,tw::Int l)
		{
			range = r;
			curr = c;
			first = f;
			last = l;
		}
		tw::Int global_count() const { return curr; }
		tw::Int count() const { return curr-first; }
		REF operator * ()
		{
			return range->GetReference(curr);
		}
		iterator& operator ++ ()
		{
			curr++;
			return *this;
		}
		friend bool operator < (const iterator& i1,const iterator& i2)
		{
			return i1.curr < i2.curr;
		}
		friend bool operator != (const iterator& i1,const iterator& i2)
		{
			return i1.curr != i2.curr;
		}
	};

	class cell
	{
		tw::Int stride[5];
		tw::Int n0,i0,j0,k0; // zero offset indexing
		tw::Int n,i,j,k; // standard indexing

	public:
		cell(const StaticSpace& ss,const tw::Int& n,const tw::Int& i,const tw::Int& j,const tw::Int& k)
		{
			for (auto ax=0;ax<5;ax++) {
				stride[ax] = ss.Stride(ax);
			}
			this->n = n;
			this->i = i;
			this->j = j;
			this->k = k;
			n0 = n - ss.LFG(0);
			i0 = i - ss.LFG(1);
			j0 = j - ss.LFG(2);
			k0 = k - ss.LFG(3);
		}
		void Decode(tw::Int *x,tw::Int *y,tw::Int *z) const
		{
			*x = i;
			*y = j;
			*z = k;
		}
		tw::Int dcd1() const { return i; }
		tw::Int dcd2() const { return j; }
		tw::Int dcd3() const { return k; }
		tw::Int Index(const tw::Int& c) const
		{
			return stride[0]*n0 + stride[1]*i0 + stride[2]*j0 + stride[3]*k0 + stride[4]*c;
		}
	};

	class strip
	{
		tw::Int strip_stride,comp_stride;
		tw::Int strip_lfg,off;
		tw::Int lfg[4],x[4],iter_ax[4],strip_ax[4],fixed_ax[4];

	public:
		strip(const StaticSpace& ss,const tw::Int& strip_ax,const tw::Int& fixed_ax,const node4& coord)
		{
			strip_stride = ss.Stride(strip_ax);
			strip_lfg = ss.LFG(strip_ax);
			comp_stride = ss.Stride(4);
			off = 0;
			for (auto i=0; i<4; i++) {
				x[i] = coord[i];
				this->strip_ax[i] = i==strip_ax;
				this->fixed_ax[i] = i==fixed_ax;
				this->iter_ax[i] = i!=strip_ax && i!=fixed_ax;
				off += (i!=strip_ax) * (coord[i] - ss.LFG(i)) * ss.Stride(i);
			}
		}
		tw::Int Axis() const
		{
			for (auto i=0; i<4; i++) {
				if (strip_ax[i]) {
					return i;
				}
			}
			return 0;
		}
		/// get standard cell indices
		void Decode(tw::Int s, tw::Int coord[4]) const
		{
			for (auto i=0; i<4; i++) {
				coord[i] = x[i] * iter_ax[i] + x[i] * fixed_ax[i] + s * strip_ax[i];
			}
		}
		/// decode the standard cell index for axis 1
		tw::Int dcd1(tw::Int s) const { return x[1]*iter_ax[1] + x[1]*fixed_ax[1] + s*strip_ax[1]; }
		/// decode the standard cell index for axis 2
		tw::Int dcd2(tw::Int s) const { return x[2]*iter_ax[2] + x[2]*fixed_ax[2] + s*strip_ax[2]; }
		/// decode the standard cell index for axis 3
		tw::Int dcd3(tw::Int s) const { return x[3]*iter_ax[3] + x[3]*fixed_ax[3] + s*strip_ax[3]; }
		/// Get the index into the field data
		tw::Int Index(const tw::Int& s,const tw::Int& c) const
		{
			return off + (s-strip_lfg)*strip_stride + c*comp_stride;
		}
		/// Get the index into the field data assuming the strip axis is packed.
		// The compiler should be able to vectorize efficiently with respect to s.
		tw::Int Index1(const tw::Int& s,const tw::Int& c) const
		{
			return off + s - strip_lfg + c*comp_stride;
		}
	};

	/// Strip reference with explicit packed axis, for promoting compiler vectorization.
	/// Performance assumption is that AX is the packed axis.
	/// This allows Field accessors to drop the stride multiplier so the compiler can figure
	/// out that the data is packed.
	/// This specialization happens in the Field class.
	/// That is, this class appearing as an argument triggers special treatment.
	template <tw::Int AX>
	class xstrip : public strip
	{
	public:
		/// general constructor
		xstrip(const StaticSpace& ss,const tw::Int& fixed_ax,const tw::node4& coord) : strip(ss,AX,fixed_ax,coord) {
			// no need to do anything
		}
		/// special constructor for spatial strips at time 1, (i,j) map to the coordinates in the normal
		/// plane in cyclic fashion, i.e. (y,z), (z,x), or (x,y).
		xstrip(const StaticSpace& ss,const tw::Int& i,const tw::Int& j) : strip(ss,AX,0,get_coord(i,j)) {
			// no need to do anything
		}
		constexpr node4 get_coord(const tw::Int& i,const tw::Int& j) {
			if (AX==1) {
				return node4{1,1,i,j};
			} else if (AX==2) {
				return node4{1,j,1,i};
			} else {
				return node4{1,i,j,1};
			}
		}
	};
}

/// abstract base class providing automatic parallel threads
class TWRange
{
protected:
	const StaticSpace *ss;
	tw::Int first,last;
public:
	void SetCountRange(tw::Int global_size)
	{
		// The whole interaction with OpenMP occurs here.
		// We merely have to compute first and last for the current thread.
		tw::GetOMPTaskLoopRange(tw::GetOMPThreadNum(),global_size,tw::GetOMPNumThreads(),&first,&last);
	}
	tw::Int size() const { return last-first+1; }
};

/// Steps through spatial cells without concern for order or geometry.
/// The time position is fixed upon construction.
export class CellRange:public TWRange
{
	// io,jo,ko,is,js,ks are decoding offsets.
	// iCells,jCells,kCells measure offsets in index space as the count increases.

	tw::Int n;
	tw::Int io,jo,ko,is,js,ks,iCells,jCells,kCells;

public:
	CellRange(const StaticSpace& space,const tw::Int& n,strongbool include_ghost_cells)
	{
		ss = &space;
		this->n = n;
		io = include_ghost_cells==strongbool::yes ? 0 : space.Layers(1);
		jo = include_ghost_cells==strongbool::yes ? 0 : space.Layers(2);
		ko = include_ghost_cells==strongbool::yes ? 0 : space.Layers(3);
		is = space.LFG(1);
		js = space.LFG(2);
		ks = space.LFG(3);
		iCells = include_ghost_cells==strongbool::yes ? space.Num(1) : space.Dim(1);
		jCells = include_ghost_cells==strongbool::yes ? space.Num(2) : space.Dim(2);
		kCells = include_ghost_cells==strongbool::yes ? space.Num(3) : space.Dim(3);
		SetCountRange(iCells*jCells*kCells);
	}
	tw::iterator<CellRange,tw::cell> begin() { return tw::iterator<CellRange,tw::cell>(this,first,first,last); }
	tw::iterator<CellRange,tw::cell> end() { return tw::iterator<CellRange,tw::cell>(this,last+1,first,last); }
	tw::cell GetReference(tw::Int global_count)
	{
		// the access order is considered unimportant
		const tw::Int i0 = io + global_count/(jCells*kCells);
		const tw::Int j0 = jo + (global_count%(jCells*kCells)) / kCells;
		const tw::Int k0 = ko + (global_count%(jCells*kCells)) % kCells;
		return tw::cell(*ss,n,i0+is,j0+js,k0+ks);
	}
};

export class EntireCellRange:public CellRange
{
public:
	EntireCellRange(const StaticSpace& space,const tw::Int& n) : CellRange(space,n,strongbool::yes)
	{
	}
};

export class InteriorCellRange:public CellRange
{
public:
	InteriorCellRange(const StaticSpace& space,const tw::Int& n) : CellRange(space,n,strongbool::no)
	{
	}
};

/// Take a 3D subpace and iterate over parallel strips.
/// Strips of ghost cells are included if includeGhostCells=true in constructor.
/// Ordinary integer loop is intended to be nested within for stepping through cells along the strip.
/// A cell is indexed by counting cells along the strip in the usual way (interior cells start at 1).
export class StripRange:public TWRange
{
protected:
	tw::Int strip_ax,fixed_ax;
	tw::Int ref[4],D[4],M[4];
public:
	StripRange(const StaticSpace& space,tw::Int strip_ax,tw::Int fixed_ax,tw::Int fixed_pos,strongbool include_ghost_cells)
	{
		tw::Int N[4];
		// This sets us up so we can do coord[i] = global_count/D[i]%M[i] for any axis.
		// (it will come to 0 for strip_ax or fixed_ax)
		this->ss = &space;
		this->fixed_ax = fixed_ax;
		this->strip_ax = strip_ax;
		tw::Int max_count = 1, last_hidden_dim = 1;
		for (auto i=0; i<4; i++) {
			N[i] = include_ghost_cells==strongbool::yes ? space.Num(i) : space.Dim(i);
			if (i==strip_ax) {
				ref[i] = space.LFG(i);
				D[i] = 1;
				M[i] = 1;
			} else if (i==fixed_ax) {
				ref[i] = fixed_pos;
				D[i] = 1;
				M[i] = 1;
			} else {
				ref[i] = include_ghost_cells==strongbool::yes ? space.LFG(i) : 1;
				D[i] = last_hidden_dim == 1 ? 1 : last_hidden_dim;
				M[i] = last_hidden_dim == 1 ? N[i] : 0xffffffff; // second alt is anything > max(dim)
				max_count *= N[i];
				last_hidden_dim = N[i];
			}
		}
		SetCountRange(max_count);
	}
	tw::iterator<StripRange,tw::strip> begin() { return tw::iterator<StripRange,tw::strip>(this,first,first,last); }
	tw::iterator<StripRange,tw::strip> end() { return tw::iterator<StripRange,tw::strip>(this,last+1,first,last); }
	tw::strip GetReference(tw::Int global_count)
	{
		// The access order across strips is considered unimportant.
		// We only care about access order along strips, which the caller will control.
		tw::node4 coord = {
			ref[0] + (global_count/D[0]) % M[0],
			ref[1] + (global_count/D[1]) % M[1],
			ref[2] + (global_count/D[2]) % M[2],
			ref[3] + (global_count/D[3]) % M[3]
		};
		return tw::strip(*ss,strip_ax,fixed_ax,coord);
	}
};

export template <tw::Int AX>
class VectorStripRange:public StripRange
{
public:
	VectorStripRange(const StaticSpace& space,tw::Int fixed_ax,tw::Int fixed_pos,bool include_ghost_cells):
		StripRange(space,AX,fixed_ax,fixed_pos,include_ghost_cells==true ? strongbool::yes : strongbool::no)
	{
		// no need to do anything
	}
	tw::iterator<VectorStripRange,tw::xstrip<AX>> begin()
	{
		return tw::iterator<VectorStripRange,tw::xstrip<AX>>(this,first,first,last);
	}
	tw::iterator<VectorStripRange,tw::xstrip<AX>> end()
	{
		return tw::iterator<VectorStripRange,tw::xstrip<AX>>(this,last+1,first,last);
	}
	tw::xstrip<AX> GetReference(tw::Int global_count)
	{
		// The access order across strips is considered unimportant.
		// We only care about access order along strips, which the caller will control.
		tw::node4 coord = {
			ref[0] + (global_count/D[0]) % M[0],
			ref[1] + (global_count/D[1]) % M[1],
			ref[2] + (global_count/D[2]) % M[2],
			ref[3] + (global_count/D[3]) % M[3]
		};
		return tw::xstrip<AX>(*ss,fixed_ax,coord);
	}
};
