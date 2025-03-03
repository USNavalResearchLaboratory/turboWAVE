module;

#include "tw_includes.h"

export module tw_iterator;
import base;
import discrete_space;

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
// Range classes are created from a DiscreteSpace, and therefore may only be used in Field objects created from
// the same DiscreteSpace.

export namespace tw
{
	template <class RNG,class REF>
	class iterator
	{
		RNG *range;
		// The RNG class is responsible for decoding the count.
		// The decoding produces a triple that Field can use to retrieve floating point data.
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
		tw::Int i0,j0,k0; // zero offset indexing
		tw::Int i,j,k; // standard indexing

	public:
		cell(const tw::Int& i0,const tw::Int& j0,const tw::Int& k0,const tw::Int& i,const tw::Int& j,const tw::Int& k)
		{
			this->i0 = i0;
			this->j0 = j0;
			this->k0 = k0;
			this->i = i;
			this->j = j;
			this->k = k;
		}
		cell(const DiscreteSpace& ds,const tw::Int& i,const tw::Int& j,const tw::Int& k)
		{
			this->i = i;
			this->j = j;
			this->k = k;
			i0 = i - ds.LFG(1);
			j0 = j - ds.LFG(2);
			k0 = k - ds.LFG(3);
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
		tw::Int Index(const tw::Int& c,const tw::Int stride[4]) const
		{
			return stride[0]*c + stride[1]*i0 + stride[2]*j0 + stride[3]*k0;
		}
	};

	class strip
	{
		tw::Int ax;
		tw::Int i0,j0,k0; // zero offset indexing
		tw::Int i,j,k; // standard indexing
		tw::Int di,dj,dk; // id the strip direction

	public:
		strip(const tw::Int& ax,const tw::Int& i0,const tw::Int& j0,const tw::Int& k0,const tw::Int& i,const tw::Int& j,const tw::Int& k)
		{
			this->ax = ax;
			this->i0 = i0;
			this->j0 = j0;
			this->k0 = k0;
			this->i = i;
			this->j = j;
			this->k = k;
			di = tw::Int(ax==1);
			dj = tw::Int(ax==2);
			dk = tw::Int(ax==3);
		}
		strip(const tw::Int& ax,const DiscreteSpace& ds,const tw::Int& i,const tw::Int& j,const tw::Int& k)
		{
			this->ax = ax;
			this->i = i;
			this->j = j;
			this->k = k;
			di = tw::Int(ax==1);
			dj = tw::Int(ax==2);
			dk = tw::Int(ax==3);
			i0 = i - ds.LFG(1);
			j0 = j - ds.LFG(2);
			k0 = k - ds.LFG(3);
		}
		tw::Int Axis() const { return ax; } // axis parallel to the strip
		void Decode(tw::Int s,tw::Int *x,tw::Int *y,tw::Int *z) const
		{
			// Get the standard cell indices
			*x = (1-di)*i + di*s;
			*y = (1-dj)*j + dj*s;
			*z = (1-dk)*k + dk*s;
		}
		// The following decode the standard cell indices one at a time
		tw::Int dcd1(tw::Int s) const { return (1-di)*i + di*s; }
		tw::Int dcd2(tw::Int s) const { return (1-dj)*j + dj*s; }
		tw::Int dcd3(tw::Int s) const { return (1-dk)*k + dk*s; }
		tw::Int Index(const tw::Int& s,const tw::Int& c,const tw::Int stride[4]) const
		{
			return stride[0]*c + stride[1]*(s*di+i0) + stride[2]*(s*dj+j0) + stride[3]*(s*dk+k0);
		}
		// Following 3 functions assume packing along the axis indicated
		// The compiler should be able to vectorize efficiently with respect to s
		tw::Int Index1(const tw::Int& s,const tw::Int& c,const tw::Int stride[4]) const
		{
			return s + i0 + stride[0]*c + stride[2]*j0 + stride[3]*k0;
		}
		tw::Int Index2(const tw::Int& s,const tw::Int& c,const tw::Int stride[4]) const
		{
			return s + j0 + stride[0]*c + stride[1]*i0 + stride[3]*k0;
		}
		tw::Int Index3(const tw::Int& s,const tw::Int& c,const tw::Int stride[4]) const
		{
			return s + k0 + stride[0]*c + stride[1]*i0 + stride[2]*j0;
		}
	};

	template <tw::Int AX>
	class xstrip : public strip
	{
		// Strip reference with explicit packed axis, for promoting compiler vectorization.
		// Performance assumption is that AX is the packed axis.
		// This allows Field accessors to drop the stride multiplier so the compiler can figure
		// out that the data is packed.
		// This specialization happens in the Field class.
		// That is, this class appearing as an argument triggers special treatment.

	public:
		xstrip(const tw::Int& i0,const tw::Int& j0,const tw::Int& k0,const tw::Int& i,const tw::Int& j,const tw::Int& k) : strip(AX,i0,j0,k0,i,j,k)
		{
			// no need to do anything
		}
		xstrip(const DiscreteSpace& ds,const tw::Int& i,const tw::Int& j,const tw::Int& k) : strip(AX,ds,i,j,k)
		{
			// no need to do anything
		}
	};
}


class TWRange
{
	// TWRange is an abstract base class
	// The abstact class defines no encoding; the encoding is particular to derived classes.

protected:
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

export class CellRange:public TWRange
{
	// Step through cells without concern for geometry.
	// Ghost cells are included if includeGhostCells=true in constructor.

	// INTERNALS:

	// io,jo,ko,is,js,ks are decoding offsets.
	// iCells,jCells,kCells measure offsets in index space as the count increases.

	tw::Int io,jo,ko,is,js,ks,iCells,jCells,kCells;

public:
	CellRange(const DiscreteSpace& space,bool includeGhostCells)
	{
		io = includeGhostCells ? 0 : space.Layers(1);
		jo = includeGhostCells ? 0 : space.Layers(2);
		ko = includeGhostCells ? 0 : space.Layers(3);
		is = space.LFG(1);
		js = space.LFG(2);
		ks = space.LFG(3);
		iCells = includeGhostCells ? space.Num(1) : space.Dim(1);
		jCells = includeGhostCells ? space.Num(2) : space.Dim(2);
		kCells = includeGhostCells ? space.Num(3) : space.Dim(3);
		SetCountRange(iCells*jCells*kCells);
	}
	tw::iterator<CellRange,tw::cell> begin() { return tw::iterator<CellRange,tw::cell>(this,first,first,last); }
	tw::iterator<CellRange,tw::cell> end() { return tw::iterator<CellRange,tw::cell>(this,last+1,first,last); }
	tw::cell GetReference(tw::Int global_count)
	{
		const tw::Int i0 = io + global_count/(jCells*kCells);
		const tw::Int j0 = jo + (global_count%(jCells*kCells)) / kCells;
		const tw::Int k0 = ko + (global_count%(jCells*kCells)) % kCells;
		return tw::cell(i0,j0,k0,i0+is,j0+js,k0+ks);
	}
};

export class EntireCellRange:public CellRange
{
public:
	EntireCellRange(const DiscreteSpace& space) : CellRange(space,true)
	{
	}
};

export class InteriorCellRange:public CellRange
{
public:
	InteriorCellRange(const DiscreteSpace& space) : CellRange(space,false)
	{
	}
};

export class StripRange:public TWRange
{
	// Step across and along strips in a given direction.
	// Strips of ghost cells are included if includeGhostCells=true in constructor.
	// Ordinary integer loop is intended to be nested within for stepping through cells along the strip.
	// A cell is indexed by counting cells along the strip in the usual way (interior cells start at 1).

	// INTERNALS:

	// io,jo,ko,is,js,ks are decoding offsets.
	// di,dj,dk are 0 except for the element in the strip direction which is 1

protected:
	tw::Int ax,di,dj,dk,io,jo,ko,is,js,ks,iCells,jCells,kCells;

public:
	StripRange(const DiscreteSpace& space,tw::Int axis,strongbool includeGhostCells)
	{
		ax = axis;
		di = tw::Int(ax==1);
		dj = tw::Int(ax==2);
		dk = tw::Int(ax==3);
		is = space.LFG(1);
		js = space.LFG(2);
		ks = space.LFG(3);
		io = includeGhostCells==strongbool::yes ? 0 : space.Layers(1);
		jo = includeGhostCells==strongbool::yes ? 0 : space.Layers(2);
		ko = includeGhostCells==strongbool::yes ? 0 : space.Layers(3);
		// the reference cell index in the strip direction is indirectly visible to caller because it interacts with
		// the caller's cell index in the strip direction.  The caller expects the cell index to be consistent
		// with the DiscreteSpace standard indexing (i.e., 1...dim label the interior cells).
		// Therefore the reference index in the strip direction has to be treated specially.
		io = (1-di)*io - di*is;
		jo = (1-dj)*jo - dj*js;
		ko = (1-dk)*ko - dk*ks;
		iCells = includeGhostCells==strongbool::yes ? space.Num(1) : space.Dim(1);
		jCells = includeGhostCells==strongbool::yes ? space.Num(2) : space.Dim(2);
		kCells = includeGhostCells==strongbool::yes ? space.Num(3) : space.Dim(3);
		SetCountRange(di*jCells*kCells + dj*iCells*kCells + dk*iCells*jCells);
	}
	tw::iterator<StripRange,tw::strip> begin() { return tw::iterator<StripRange,tw::strip>(this,first,first,last); }
	tw::iterator<StripRange,tw::strip> end() { return tw::iterator<StripRange,tw::strip>(this,last+1,first,last); }
	tw::strip GetReference(tw::Int global_count)
	{
		const tw::Int i0 = di*io + dj*(io + global_count/kCells) + dk*(io + global_count/jCells);
		const tw::Int j0 = di*(jo + global_count%jCells) + dj*jo + dk*(jo + global_count%jCells);
		const tw::Int k0 = di*(ko + global_count/jCells) + dj*(ko + global_count%kCells) + dk*ko;
		return tw::strip(ax,i0,j0,k0,i0+is,j0+js,k0+ks);
	}
};

export template <tw::Int AX>
class VectorStripRange:public StripRange
{
public:
	VectorStripRange(const DiscreteSpace& space,bool includeGhostCells) : StripRange(space,AX,includeGhostCells==true ? strongbool::yes : strongbool::no)
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
		const tw::Int i0 = di*io + dj*(io + global_count/kCells) + dk*(io + global_count/jCells);
		const tw::Int j0 = di*(jo + global_count%jCells) + dj*jo + dk*(jo + global_count%jCells);
		const tw::Int k0 = di*(ko + global_count/jCells) + dj*(ko + global_count%kCells) + dk*ko;
		return tw::xstrip<AX>(i0,j0,k0,i0+is,j0+js,k0+ks);
	}
};
