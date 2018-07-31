// TurboWAVE iterators work differently from standard library, i.e., are non-conforming.
// Rather than dereferencing the iterator to get an object, the iterator is passed
// as an argument to the Field class to get a floating point number.
// This seems more suitable since we often want to access different Field instances
// using the same iterator (i.e., grabbing a number from different fields at the same space-point).
// This also allows us to create special iterators designed to encourage various optimizations.

// Iterators are aware of OpenMP state and automatically split the work accross threads.
// Iterators are compatible only with Fields derived from the same DiscreteSpace

class TWIterator
{
	// TWIterator is an abstract base class

	// curr is the encoded index of the current object, called the count.
	// The global encoding is such that count ranges from (0...N-1), where N is the number of objects.
	// curr contains the GLOBAL count, which the user does not interact with.
	// count() returns the LOCAL count, which ranges from (0...M-1), where M is the number of objects to be processed by the current thread.
	// first and last bound the global count range for the thread (i.e., curr = first + local_count <= last)

	// The abstact class defines no encoding; the encoding is particular to derived classes.
	// This leads to derived classes defining their own assignment and increment operators.
	// One could use virtual functions to make this more elegant, but efficiency might suffer.

	protected:
	tw::Int curr;
	tw::Int first,last;

	public:
	void SetCountRange(tw::Int global_size)
	{
		// The whole interaction with OpenMP occurs here.
		// We merely have to compute first and last for the current thread.
		tw::GetOMPTaskLoopRange(tw::GetOMPThreadNum(),global_size,tw::GetOMPNumThreads(),&first,&last);
		curr = first;
	}
	tw::Int global_count() const { return curr; }
	tw::Int count() const { return curr-first; }
	tw::Int size() const { return last-first+1; }
	tw::Int begin() const { return 0; }
	tw::Int end() const { return last-first+1; }
	friend bool operator < (const TWIterator& cell,tw::Int n)
	{
		return cell.curr < cell.first + n;
	}
};

class CellIterator:public TWIterator
{
	// Step through cells without concern for geometry.
	// Ghost cells are included if includeGhostCells=true in constructor.

	// INTERNALS:

	// i0,j0,k0 labels the current cell; io,jo,ko,is,js,ks are decoding offsets.
	// i0,j0,k0 are indexed with zero-offset regardless of ghost cell layers (first cell of any kind is labelled 0)
	// This should not confuse users since i0,j0,k0 are not visible to the caller

	tw::Int i0,j0,k0,io,jo,ko,is,js,ks,iCells,jCells,kCells;

	public:
	CellIterator(const DiscreteSpace& space,bool includeGhostCells)
	{
		io = includeGhostCells ? 0 : space.Layers(1);
		jo = includeGhostCells ? 0 : space.Layers(2);
		ko = includeGhostCells ? 0 : space.Layers(3);
		is = 1 - space.Layers(1);
		js = 1 - space.Layers(2);
		ks = 1 - space.Layers(3);
		iCells = includeGhostCells ? space.Num(1) : space.Dim(1);
		jCells = includeGhostCells ? space.Num(2) : space.Dim(2);
		kCells = includeGhostCells ? space.Num(3) : space.Dim(3);
		SetCountRange(iCells*jCells*kCells);
		SetDecodingRefs();
	}
	CellIterator& operator = (tw::Int local_count)
	{
		curr = first + local_count;
		SetDecodingRefs();
		return *this;
	}
	CellIterator& operator ++ ()
	{
		curr++;
		SetDecodingRefs();
		return *this;
	}
	void SetDecodingRefs()
	{
		i0 = io + curr/(jCells*kCells);
		j0 = jo + (curr%(jCells*kCells)) / kCells;
		k0 = ko + (curr%(jCells*kCells)) % kCells;
	}
	void Decode(tw::Int *x,tw::Int *y,tw::Int *z) const
	{
		// Get the standard cell indices
		*x = i0 + is;
		*y = j0 + js;
		*z = k0 + ks;
	}
	// The following decode the standard cell indices one at a time
	tw::Int dcd1() const { return i0 + is; }
	tw::Int dcd2() const { return j0 + js; }
	tw::Int dcd3() const { return k0 + ks; }
	void SetCell(tw::Int i,tw::Int j,tw::Int k)
	{
		// Avoid using.
		// Iterating after an explicit SetCell is not supported.
		i0 = i - is;
		j0 = j - js;
		k0 = k - ks;
	}
	tw::Int Index(const tw::Int& c,const tw::Int stride[4]) const
	{
		return stride[0]*c + stride[1]*i0 + stride[2]*j0 + stride[3]*k0;
	}
};

class StripIterator:public TWIterator
{
	// Step across and along strips in a given direction.
	// Strips of ghost cells are included if includeGhostCells=true in constructor.
	// Ordinary integer loop is intended to be nested within for stepping through cells along the strip.
	// A cell is indexed by counting cells along the strip in the usual way (interior cells start at 1).

	// INTERNALS:

	// i0,j0,k0 labels the current reference cell; io,jo,ko,is,js,ks are decoding offsets.
	// i0,j0,k0 are indexed with zero-offset regardless of ghost cell layers (first cell of any kind is labelled 0)
	// This should not confuse users since i0,j0,k0 are not visible to the caller
	// di,dj,dk are 0 except for the element in the strip direction which is 1

	// the reference cell index in the strip direction is indirectly visible to callers because it interacts with
	// the caller's cell index in the strip direction.  The caller expects the cell index to be consistent
	// with the DiscreteSpace standard indexing (i.e., 1...dim label the interior cells).
	// Therefore the reference index in the strip direction has to be treated specially.

	private:
	tw::Int ax,di,dj,dk,io,jo,ko,is,js,ks,iCells,jCells,kCells;
	tw::Int dim,num,n0,n1;
	protected:
	tw::Int i0,j0,k0;

	public:
	StripIterator(const DiscreteSpace& space,tw::Int axis,strongbool includeGhostCells)
	{
		ax = axis;
		di = tw::Int(ax==1);
		dj = tw::Int(ax==2);
		dk = tw::Int(ax==3);
		is = 1 - space.Layers(1);
		js = 1 - space.Layers(2);
		ks = 1 - space.Layers(3);
		io = includeGhostCells==strongbool::yes ? 0 : space.Layers(1);
		jo = includeGhostCells==strongbool::yes ? 0 : space.Layers(2);
		ko = includeGhostCells==strongbool::yes ? 0 : space.Layers(3);
		io = (1-di)*io - di*is;
		jo = (1-dj)*jo - dj*js;
		ko = (1-dk)*ko - dk*ks;
		iCells = includeGhostCells==strongbool::yes ? space.Num(1) : space.Dim(1);
		jCells = includeGhostCells==strongbool::yes ? space.Num(2) : space.Dim(2);
		kCells = includeGhostCells==strongbool::yes ? space.Num(3) : space.Dim(3);
		dim = space.Dim(ax);
		num = space.Num(ax);
		n0 = space.N0(ax);
		n1 = space.N1(ax);
		SetCountRange(di*jCells*kCells + dj*iCells*kCells + dk*iCells*jCells);
		SetDecodingRefs();
	}
	StripIterator& operator = (tw::Int local_count)
	{
		curr = first + local_count;
		SetDecodingRefs();
		return *this;
	}
	StripIterator& operator ++ ()
	{
		curr++;
		SetDecodingRefs();
		return *this;
	}
	void SetDecodingRefs()
	{
		i0 = di*io;
		//j0 = di*(jo + curr/kCells);
		//k0 = di*(ko + curr%kCells);
		// alternate encoding (choice may affect performance)
		j0 = di*(jo + curr%jCells);
		k0 = di*(ko + curr/jCells);

		i0 += dj*(io + curr/kCells);
		j0 += dj*jo;
		k0 += dj*(ko + curr%kCells);

		i0 += dk*(io + curr/jCells);
		j0 += dk*(jo + curr%jCells);
		k0 += dk*ko;
	}
	void Decode(tw::Int s,tw::Int *x,tw::Int *y,tw::Int *z) const
	{
		// Get the standard cell indices
		*x = (1-di)*(i0 + is) + di*s;
		*y = (1-dj)*(j0 + js) + dj*s;
		*z = (1-dk)*(k0 + ks) + dk*s;
	}
	// The following decode the standard cell indices one at a time
	tw::Int dcd1(tw::Int s) const { return (1-di)*(i0 + is) + di*s; }
	tw::Int dcd2(tw::Int s) const { return (1-dj)*(j0 + js) + dj*s; }
	tw::Int dcd3(tw::Int s) const { return (1-dk)*(k0 + ks) + dk*s; }
	void SetStrip(tw::Int i,tw::Int j,tw::Int k)
	{
		// Avoid using.
		// Iterating after an explicit SetStrip is not supported.
		i0 = (1-di)*(i - is) + di*io;
		j0 = (1-dj)*(j - js) + dj*jo;
		k0 = (1-dk)*(k - ks) + dk*ko;
	}
	tw::Int Axis() const { return ax; } // axis parallel to the strip
	tw::Int Dim() const { return dim; } // interior cells parallel to the strip
	tw::Int Num() const { return num; } // total cells parallel to the strip
	tw::Int N0() const { return n0; } // index of first cell along strip (may be ghost cell)
	tw::Int N1() const { return n1; } // index of last cell along strip (may be ghost cell)
	tw::Int Index(const tw::Int& s,const tw::Int& c,const tw::Int stride[4]) const
	{
		return stride[0]*c + stride[1]*(s*di+i0) + stride[2]*(s*dj+j0) + stride[3]*(s*dk+k0);
	}
};

template <tw::Int AX>
class VectorizingIterator : public StripIterator
{
	// Strip iterator for promoting compiler vectorization.
	// Performance assumption is that AX is the packed axis.
	// This allows Field accessors to drop the stride multiplier so the compiler can figure
	// out that the data is packed.
	// This specialization happens in the Field class.
	// That is, this class appearing as an argument triggers special treatment.

	public:
	VectorizingIterator(const DiscreteSpace& space,bool includeGhostCells) : StripIterator(space,AX,includeGhostCells?strongbool::yes:strongbool::no)
	{
		// no need to do anything
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
