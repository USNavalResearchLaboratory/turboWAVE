// tw_iterator.h is included mid-file
enum tw_geometry {cartesian,cylindrical,spherical};
enum tw_boundary_spec {cyclic,reflecting,absorbing,emitting,axisymmetric};

struct Primitive
{
	float x[3];	// x[0],x[1],x[2] = relative position in cell using interval [-0.5,0.5); linearly mapped to curvilinear coordinates in cell
	tw::Int cell; // encoded index of the reference cell
	Primitive()
	{
		cell = 0;
		x[0] = x[1] = x[2] = 0.0;
	}
	Primitive(tw::Int c,float x0,float x1,float x2)
	{
		cell = c;
		x[0] = x0;
		x[1] = x1;
		x[2] = x2;
	}
};

struct weights_3D
{
	tw::Float w[3][3];
	tw::Int cell;
};

struct DiscreteSpace
{
	// This object is an indexing and interpolation scheme for a structured grid.
	// It defines a cell encoding, which labels cells with a single integer in unbroken sequence.
	// Decoded cell indices are triple integers defined as follows:
	// The interior cells are labeled [1,...,dim[ax]], where ax is an axis numbered from 1 to 3
	// Ghost cell layers are labeled [1-L,...,0] and [dim[ax]+1,...,dim[ax]+L]
	// Note that the cell encoding is not required to be related to storage of field data in the cells.

	protected:

	tw::Float dt; // timestep, only for constant timestep
	tw::vec3 corner,size,globalCorner,globalSize,spacing,freq; // spacing and freq are only for uniform grids
	// Indexing in the following has 0=components, (1,2,3)=spatial axes
	// num includes ghost cells, dim does not
	// Element 0 only becomes meaningful when a Field class is derived
	tw::Int dim[4],num[4];
	// Cell strides help encode a cell; they have nothing to do with storage patterns
	// decodingStride does not use 0 strides, encodingStride does
	// The Field class governs packing of quantities in cells
	// matching the encoding and storage patterns may boost performance
	tw::Int encodingStride[4],decodingStride[4];
	// The following are the indices of the extreme cells (which may be ghost cells)
	tw::Int lb[4],ub[4];
	// The following is the number of layers of ghost cells, layers[0] holds the maximum layers
	tw::Int layers[4];
	// Identify ignorable dimensions, 1 if axis is ignorable, 0 otherwise
	tw::Int ignorable[4];

	public:

	DiscreteSpace();
	DiscreteSpace(tw::Float dt,tw::Int x,tw::Int y,tw::Int z,const tw::vec3& corner,const tw::vec3& size,tw::Int ghostCellLayers);
	void Resize(tw::Int x,tw::Int y,tw::Int z,const tw::vec3& corner,const tw::vec3& size,tw::Int ghostCellLayers);
	void Resize(tw::Int x,tw::Int y,tw::Int z,const tw::vec3& corner,const tw::vec3& size);
	tw::Int EncodeCell(tw::Int i,tw::Int j,tw::Int k) const { return (i-lb[1])*encodingStride[1] + (j-lb[2])*encodingStride[2] + (k-lb[3])*encodingStride[3]; }
	void DecodeCell(const tw::Int& cell,tw::Int ijk[4]) const
	{
		// Assumes z-packed encoding pattern, i.e. decodingStride[3]=1
		ijk[1] = lb[1] + cell / decodingStride[1];
		ijk[2] = lb[2] + (cell % decodingStride[1]) / decodingStride[2];
		ijk[3] = lb[3] + (cell % decodingStride[1]) % decodingStride[2];
	}
	void DecodeCell(const tw::Int& cell,tw::Int *i,tw::Int *j,tw::Int *k) const
	{
		*i = lb[1] + cell / decodingStride[1];
		*j = lb[2] + (cell % decodingStride[1]) / decodingStride[2];
		*k = lb[3] + (cell % decodingStride[1]) % decodingStride[2];
	}
	void DecodeCell(const Primitive& q,tw::Int ijk[4]) const { DecodeCell(q.cell,ijk); }
	void DecodeCell(const Primitive& q,tw::Int *i,tw::Int *j,tw::Int *k) const { DecodeCell(q.cell,i,j,k); }
	bool PrimitiveInCell(const Primitive& q) const { return q.x[0]>=-0.5 && q.x[0]<0.5 && q.x[1]>=-0.5 && q.x[1]<0.5 && q.x[2]>=-0.5 && q.x[2]<0.5; }
	bool PrimitiveInDomain(const Primitive& q) const;
	void SetPrimitiveWithPosition(Primitive& q,const tw::vec3& pos) const;
	void UpdatePrimitiveWithPosition(Primitive& q,const tw::vec3& pos) const;
	tw::vec3 PositionFromPrimitive(const Primitive& q) const;
	void MinimizePrimitive(Primitive& q) const;
	void MinimizePrimitive(tw::Int *cell,tw::Int ijk[3][tw::max_bundle_size],float x[3][tw::max_bundle_size],float domainMask[tw::max_bundle_size]) const;
	//void MinimizePrimitiveBohm(tw::Int *cell,tw::Int ijk[3][tw::max_bundle_size],float x[3][tw::max_bundle_size],float domainMask[tw::max_bundle_size]) const;
	tw::Int Dim(const tw::Int& ax) const { return dim[ax]; }
	tw::Int Num(const tw::Int& ax) const { return num[ax]; }
	tw::Int Layers(const tw::Int& ax) const { return layers[ax]; }
	tw::Int N0(const tw::Int& ax) const { return lb[ax]; }
	tw::Int N1(const tw::Int& ax) const { return ub[ax]; }
	friend tw::Float dt(const DiscreteSpace& A) { return A.dt; }
	friend tw::Float dx(const DiscreteSpace& A)  { return A.spacing.x; }
	friend tw::Float dy(const DiscreteSpace& A)  { return A.spacing.y; }
	friend tw::Float dz(const DiscreteSpace& A)  { return A.spacing.z; }
	friend tw::Float dxi(const DiscreteSpace& A)  { return A.freq.x; }
	friend tw::Float dyi(const DiscreteSpace& A)  { return A.freq.y; }
	friend tw::Float dzi(const DiscreteSpace& A)  { return A.freq.z; }
	friend tw::vec3 PhysicalSize(const DiscreteSpace& A)  { return A.size; }
	friend tw::vec3 GlobalPhysicalSize(const DiscreteSpace& A)  { return A.globalSize; }
	friend tw::vec3 Corner(const DiscreteSpace& A) { return A.corner; }
	friend tw::vec3 GlobalCorner(const DiscreteSpace& A)  { return A.globalCorner; }
	bool IsWithinInset(const Primitive& q,tw::Int inset) const;
	bool IsPointValid(const tw::vec3& P);
	tw::Int Dimensionality();
	bool PowersOfTwo();
	bool TransversePowersOfTwo();
	void GetWeights(weights_3D* weights,const tw::vec3& P);
	void GetWeights(weights_3D* weights,const Primitive& q);
	void GetWeights(float w[3][3][tw::max_bundle_size],float x[3][tw::max_bundle_size]);
	void GetWallWeights(float w[3][3][tw::max_bundle_size],float x[3][tw::max_bundle_size]);

	void ReadData(std::ifstream& inFile);
	void WriteData(std::ofstream& outFile);

	#ifdef USE_OPENCL
	void CellUpdateProtocol(cl_kernel k,cl_command_queue q);
	void ElementUpdateProtocol(cl_kernel k,cl_command_queue q);
	void LocalUpdateProtocol(cl_kernel k,cl_command_queue q);
	void PointUpdateProtocol(cl_kernel k,cl_command_queue q);
	#endif
};

inline bool DiscreteSpace::IsWithinInset(const Primitive& q,tw::Int inset) const
{
	tw::Int i,j,k;
	DecodeCell(q,&i,&j,&k);
	return (dim[1]==1 || (i>=lb[1]+inset && i<=ub[1]-inset))
		&& (dim[2]==1 || (j>=lb[2]+inset && j<=ub[2]-inset))
		&& (dim[3]==1 || (k>=lb[3]+inset && k<=ub[3]-inset));
}

inline bool DiscreteSpace::IsPointValid(const tw::vec3& P)
{
	const tw::vec3 PLoc = P - corner;
	return (PLoc.x>=0.0 && PLoc.x<size.x && PLoc.y>=0.0 && PLoc.y<size.y && PLoc.z>=0.0 && PLoc.z<size.z);
}

inline tw::Int DiscreteSpace::Dimensionality()
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

inline bool DiscreteSpace::PrimitiveInDomain(const Primitive& q) const
{
	const tw::vec3 r = PositionFromPrimitive(q) - corner;
	return	((r.x>=0.0 && r.x<size.x) || dim[1]==1) &&
			((r.y>=0.0 && r.y<size.y) || dim[2]==1) &&
			((r.z>=0.0 && r.z<size.z) || dim[3]==1);
}

inline tw::vec3 DiscreteSpace::PositionFromPrimitive(const Primitive& q) const
{
	tw::Int ijk[4];
	tw::vec3 ans;
	DecodeCell(q,ijk);
	ans[0] = corner[0] + spacing[0] * (tw::Float(q.x[0]) + tw::Float(ijk[1]) - tw::Float(0.5));
	ans[1] = corner[1] + spacing[1] * (tw::Float(q.x[1]) + tw::Float(ijk[2]) - tw::Float(0.5));
	ans[2] = corner[2] + spacing[2] * (tw::Float(q.x[2]) + tw::Float(ijk[3]) - tw::Float(0.5));
	return ans;
}

inline void DiscreteSpace::SetPrimitiveWithPosition(Primitive& q,const tw::vec3& P) const
{
	const tw::vec3 PLoc = P - corner;
	const tw::Int i = tw::Int(floor(PLoc.x*freq.x)) + 1;
	const tw::Int j = tw::Int(floor(PLoc.y*freq.y)) + 1;
	const tw::Int k = tw::Int(floor(PLoc.z*freq.z)) + 1;
	//const tw::Int i = MyFloor(PLoc.x*freq.x) + 1;
	//const tw::Int j = MyFloor(PLoc.y*freq.y) + 1;
	//const tw::Int k = MyFloor(PLoc.z*freq.z) + 1;
	q.cell = EncodeCell(i,j,k);
	q.x[0] = PLoc.x*freq.x - tw::Float(i) + tw::Float(0.5);
	q.x[1] = PLoc.y*freq.y - tw::Float(j) + tw::Float(0.5);
	q.x[2] = PLoc.z*freq.z - tw::Float(k) + tw::Float(0.5);
}

inline void DiscreteSpace::UpdatePrimitiveWithPosition(Primitive& q,const tw::vec3& P) const
{
	// This routine does NOT change the reference cell, i.e., primitive is not minimized
	tw::Int ijk[4];
	const tw::vec3 PLoc = P - corner;
	DecodeCell(q,ijk);
	q.x[0] = PLoc[0]*freq[0] - tw::Float(ijk[1]) + tw::Float(0.5);
	q.x[1] = PLoc[1]*freq[1] - tw::Float(ijk[2]) + tw::Float(0.5);
	q.x[2] = PLoc[2]*freq[2] - tw::Float(ijk[3]) + tw::Float(0.5);
}

inline void DiscreteSpace::MinimizePrimitive(tw::Int cell[tw::max_bundle_size],tw::Int ijk[3][tw::max_bundle_size],float x[3][tw::max_bundle_size],float domainMask[tw::max_bundle_size]) const
{
	// Bundle version of below.
	// Note non-standard indexing of decoded cell ijk[] is tolerated for sake of efficiency.
	// This also sets the domain mask to 0 for particles out of interior domain, 1 otherwise.
	const tw::Int N = tw::max_bundle_size;
	const tw::Int AB = tw::vec_align_bytes;
	alignas(AB) tw::Int itest[N];
	alignas(AB) float ftest[N];
	tw::Int par,ax,repeat;
	#pragma omp simd aligned(cell,ijk:AB)
	for (par=0;par<N;par++)
		DecodeCell(cell[par],&ijk[0][par],&ijk[1][par],&ijk[2][par]);
	#pragma omp simd aligned(domainMask:AB)
	for (par=0;par<N;par++)
		domainMask[par] = 1.0; // in domain
	for (ax=0;ax<3;ax++)
	{
		do
		{
			for (par=0;par<N;par++)
			{
				itest[par] = tw::Int(x[ax][par]<-0.5 && ijk[ax][par]>lb[ax+1]) - tw::Int(x[ax][par]>=0.5 && ijk[ax][par]<ub[ax+1]);
				ftest[par] = float(itest[par]);
			}
			#pragma omp simd aligned(itest,ijk:AB)
			for (par=0;par<N;par++)
				ijk[ax][par] -= itest[par];
			#pragma omp simd aligned(ftest,x:AB)
			for (par=0;par<N;par++)
				x[ax][par] += ftest[par];
			repeat = 0;
			for (par=0;par<N;par++)
				repeat += tw::Int(x[ax][par]<-0.5 && ijk[ax][par]>lb[ax+1]) + tw::Int(x[ax][par]>=0.5 && ijk[ax][par]<ub[ax+1]);
		} while (repeat);
		for (par=0;par<N;par++)
			ftest[par] = float(tw::Int(ijk[ax][par]>=1 && ijk[ax][par]<=dim[ax+1]));
		#pragma omp simd aligned(ftest,domainMask:AB)
		for (par=0;par<N;par++)
			domainMask[par] *= ftest[par];
	}
	#pragma omp simd aligned(cell,ijk:AB)
	for (par=0;par<N;par++)
		cell[par] = EncodeCell(ijk[0][par],ijk[1][par],ijk[2][par]);
}

inline void DiscreteSpace::MinimizePrimitive(Primitive& q) const
{
	// Change the reference cell to the one containing the particle centroid
	// Absolute position remains unchanged
	// The reference cell is pegged such that the cell is in the extended domain
	// This means x stays out of range in the event the reference cell is pegged
	tw::Int ax,ijk[4],test;
	DecodeCell(q.cell,ijk);
	for (ax=1;ax<=3;ax++)
	{
		do
		{
			test = tw::Int(q.x[ax-1]<-0.5 && ijk[ax]>lb[ax]);
			test -= tw::Int(q.x[ax-1]>=0.5 && ijk[ax]<ub[ax]);
			ijk[ax] -= test;
			q.x[ax-1] += float(test);
		} while (test);
	}
	q.cell = EncodeCell(ijk[1],ijk[2],ijk[3]);
}

inline void DiscreteSpace::GetWeights(weights_3D* weights,const tw::vec3& P)
{
	Primitive q;
	SetPrimitiveWithPosition(q,P);
	GetWeights(weights,q);
}

inline void DiscreteSpace::GetWeights(weights_3D* weights,const Primitive& q)
{
	weights->cell = q.cell;
	weights->w[0][0] = 0.125f - 0.5f*q.x[0] + 0.5f*q.x[0]*q.x[0];
	weights->w[0][1] = 0.125f - 0.5f*q.x[1] + 0.5f*q.x[1]*q.x[1];
	weights->w[0][2] = 0.125f - 0.5f*q.x[2] + 0.5f*q.x[2]*q.x[2];
	weights->w[1][0] = 0.75f - q.x[0]*q.x[0];
	weights->w[1][1] = 0.75f - q.x[1]*q.x[1];
	weights->w[1][2] = 0.75f - q.x[2]*q.x[2];
	weights->w[2][0] = 0.125f + 0.5f*q.x[0] + 0.5f*q.x[0]*q.x[0];
	weights->w[2][1] = 0.125f + 0.5f*q.x[1] + 0.5f*q.x[1]*q.x[1];
	weights->w[2][2] = 0.125f + 0.5f*q.x[2] + 0.5f*q.x[2]*q.x[2];
}

inline void DiscreteSpace::GetWeights(float w[3][3][tw::max_bundle_size],float x[3][tw::max_bundle_size])
{
	tw::Int N = tw::max_bundle_size;
	#pragma omp simd aligned(w,x:tw::vec_align_bytes)
	for (tw::Int i=0;i<N;i++)
	{
		w[0][0][i] = 0.125f - 0.5f*x[0][i] + 0.5f*x[0][i]*x[0][i];
		w[0][1][i] = 0.125f - 0.5f*x[1][i] + 0.5f*x[1][i]*x[1][i];
		w[0][2][i] = 0.125f - 0.5f*x[2][i] + 0.5f*x[2][i]*x[2][i];
		w[1][0][i] = 0.75f - x[0][i]*x[0][i];
		w[1][1][i] = 0.75f - x[1][i]*x[1][i];
		w[1][2][i] = 0.75f - x[2][i]*x[2][i];
		w[2][0][i] = 0.125f + 0.5f*x[0][i] + 0.5f*x[0][i]*x[0][i];
		w[2][1][i] = 0.125f + 0.5f*x[1][i] + 0.5f*x[1][i]*x[1][i];
		w[2][2][i] = 0.125f + 0.5f*x[2][i] + 0.5f*x[2][i]*x[2][i];
	}
}

inline void DiscreteSpace::GetWallWeights(float w[3][3][tw::max_bundle_size],float x[3][tw::max_bundle_size])
{
	tw::Int N = tw::max_bundle_size;
	#pragma omp simd aligned(w,x:tw::vec_align_bytes)
	for (tw::Int i=0;i<N;i++)
	{
		w[0][0][i] = 0.0f;
		w[0][1][i] = 0.0f;
		w[0][2][i] = 0.0f;
		w[1][0][i] = 0.5f - x[0][i];
		w[1][1][i] = 0.5f - x[1][i];
		w[1][2][i] = 0.5f - x[2][i];
		w[2][0][i] = 0.5f + x[0][i];
		w[2][1][i] = 0.5f + x[1][i];
		w[2][2][i] = 0.5f + x[2][i];
	}
}

#include "tw_iterator.h"

struct MetricSpace:DiscreteSpace
{
	tw::Float car,cyl,sph; // set variable corresponding to coordinate system to 1.0, all others to 0.0
	tw::Int mnum[4]; // num for metric arrays (see DiscreteSpace)
	tw::Int mlb[4],mub[4]; // lb and ub for metric arrays (see DiscreteSpace)
	tw::Int I3x3[3][3];

	std::valarray<tw::Float> gpos; // global positions in parameter space
	std::valarray<tw::Float> width; // cell sizes in parameter space
	// Cell metrics are packed assuming separable functional forms
	// It is assumed there is no dependence of metrics on y (true for the typical 3 systems)
	std::valarray<tw::Float> cell_area_x;
	std::valarray<tw::Float> cell_area_z;
	// Elements 0,1,2,3 are the volume and lower wall areas of the cell respectively
	// Elements 4,5,6,7 are the volume and upper wall areas of a cell shifted back by 1/2
	// External access of areas is through dS and dSh
	std::valarray<tw::Float> cell_arc_x;
	std::valarray<tw::Float> cell_arc_z;
	// Elements 0,1,2 are from cell center to cell center
	// Elements 3,4,5 are offset by 1/2 cell forward in arc direction, back in other 2
	// External access of arcs is through dl and dlh, and uses spatial indexing 1,2,3

	#ifdef USE_OPENCL
	cl_mem metricsBuffer;
	cl_mem stripBuffer[4];
	#endif

	MetricSpace();
	~MetricSpace();

	void Resize(tw::Int x,tw::Int y,tw::Int z,const tw::vec3& corner,const tw::vec3& size,tw::Int ghostCellLayers);
	void Resize(tw::Int x,tw::Int y,tw::Int z,const tw::vec3& corner,const tw::vec3& size);
	void ReadData(std::ifstream& inFile);
	void WriteData(std::ofstream& outFile);

	void SetupPositionArrays();
	void SetCartesianGeometry();
	void SetCylindricalGeometry();
	void SetSphericalGeometry();

	#ifdef USE_OPENCL
	void InitializeMetricsBuffer(cl_context ctx,tw::Float dt);
	void StripUpdateProtocol(cl_kernel k,cl_command_queue q,tw::Int axis,tw::Int stripArgument);
	#endif

	tw::vec3 ScaleFactor(const tw::vec3& r) const;
	tw::Float ScaleFactor(const tw::Int& a,const tw::vec3& r) const
	{
		return ScaleFactor(r)[a-1];
	}
	tw::Float CylindricalRadius(const tw::vec3& r) const;
	tw::Float SphericalRadius(const tw::vec3& r) const;

	tw::Int kdelta(const tw::Int& ax1,const tw::Int& ax2) const
	{
		// 3x3 Kronecker delta function
		return I3x3[ax1-1][ax2-1];
	}
	tw::Float X(const tw::Int& i,const tw::Int& ax) const
	{
		// position in parameter space
		return gpos[mnum[0]*ax-mlb[ax]+i];
	}
	tw::Float dX(const tw::Int& i,const tw::Int& ax) const
	{
		// cell size in parameter space (not an arc length)
		return width[mnum[0]*ax-mlb[ax]+i];
	}
	tw::Float& X(const tw::Int& i,const tw::Int& ax)
	{
		// position in parameter space
		return gpos[mnum[0]*ax-mlb[ax]+i];
	}
	tw::Float& dX(const tw::Int& i,const tw::Int& ax)
	{
		// cell size in parameter space (not an arc length)
		return width[mnum[0]*ax-mlb[ax]+i];
	}
	tw::vec3 Pos(const tw::Int& x,const tw::Int& y,const tw::Int& z) const
	{
		return tw::vec3(gpos[mnum[0]*1-mlb[1]+x],gpos[mnum[0]*2-mlb[2]+y],gpos[mnum[0]*3-mlb[3]+z]);
	}
	tw::vec3 Pos(const tw::cell& cell) const
	{
		return Pos(cell.dcd1(),cell.dcd2(),cell.dcd3());
	}
	tw::vec3 Pos(const tw::strip& s,const tw::Int& i) const
	{
		return Pos(s.dcd1(i),s.dcd2(i),s.dcd3(i));
	}
	tw::vec3 Pos(const tw::xstrip<3>& v,const tw::Int& k) const
	{
		return Pos(v.dcd1(k),v.dcd2(k),k);
	}
	tw::vec3 dPos(const tw::Int& x,const tw::Int& y,const tw::Int& z) const
	{
		return tw::vec3(width[mnum[0]*1-mlb[1]+x],width[mnum[0]*2-mlb[2]+y],width[mnum[0]*3-mlb[3]+z]);
	}
	tw::vec3 dPos(const tw::cell& cell) const
	{
		return dPos(cell.dcd1(),cell.dcd2(),cell.dcd3());
	}
	tw::vec3 dPos(const tw::strip& s,const tw::Int& i) const
	{
		return dPos(s.dcd1(i),s.dcd2(i),s.dcd3(i));
	}
	tw::vec3 dPos(const tw::xstrip<3>& v,const tw::Int& k) const
	{
		return dPos(v.dcd1(k),v.dcd2(k),k);
	}
	tw::Float dS(const tw::Int& x,const tw::Int& y,const tw::Int& z,const tw::Int& ax) const
	{
		// returns wall area for ax = axis normal to wall.  ax = 0 returns cell volume.
		return cell_area_x[ax*mnum[1] - mlb[1] + x] * cell_area_z[ax*mnum[3] - mlb[3] + z];
	}
	tw::Float dS(const tw::cell& cell,const tw::Int& ax) const
	{
		return dS(cell.dcd1(),cell.dcd2(),cell.dcd3(),ax);
	}
	tw::Float dS(const tw::strip& s,const tw::Int& i,const tw::Int& ax) const
	{
		return dS(s.dcd1(i),s.dcd2(i),s.dcd3(i),ax);
	}
	tw::Float dS(const tw::xstrip<1>& v,const tw::Int& i,const tw::Int& ax) const
	{
		return dS(i,v.dcd2(i),v.dcd3(i),ax);
	}
	tw::Float dS(const tw::xstrip<3>& v,const tw::Int& k,const tw::Int& ax) const
	{
		return dS(v.dcd1(k),v.dcd2(k),k,ax);
	}
	tw::Float dl(const tw::Int& x,const tw::Int& y,const tw::Int& z,const tw::Int& ax) const
	{
		// returns arc length from cell center to cell center, along axis=ax, from low side
		return cell_arc_x[(ax-1)*mnum[1] - mlb[1] + x] * cell_arc_z[(ax-1)*mnum[3] - mlb[3] + z];
	}
	tw::Float dl(const tw::cell& cell,const tw::Int& ax) const
	{
		return dl(cell.dcd1(),cell.dcd2(),cell.dcd3(),ax);
	}
	tw::Float dl(const tw::strip& s,const tw::Int& i,const tw::Int& ax) const
	{
		return dl(s.dcd1(i),s.dcd2(i),s.dcd3(i),ax);
	}
	tw::Float dl(const tw::xstrip<3>& v,const tw::Int& k,const tw::Int& ax) const
	{
		return dl(v.dcd1(k),v.dcd2(k),k,ax);
	}
	tw::Float dSh(const tw::Int& x,const tw::Int& y,const tw::Int& z,const tw::Int& ax) const
	{
		// returns wall area for ax = axis normal to wall.  ax = 0 returns cell volume.
		return dS(x,y,z,ax+4);
	}
	tw::Float dSh(const tw::cell& cell,const tw::Int& ax) const
	{
		return dS(cell,ax+4);
	}
	tw::Float dSh(const tw::xstrip<3>& v,const tw::Int& k,const tw::Int& ax) const
	{
		return dS(v,k,ax+4);
	}
	tw::Float dlh(const tw::Int& x,const tw::Int& y,const tw::Int& z,const tw::Int& ax) const
	{
		// returns arc length from cell wall to cell wall, along axis=ax, along low side edge
		return dl(x,y,z,ax+3);
	}
	tw::Float dlh(const tw::cell& cell,const tw::Int& ax) const
	{
		return dl(cell,ax+3);
	}
	tw::Float dlh(const tw::xstrip<3>& v,const tw::Int& k,const tw::Int& ax) const
	{
		return dl(v,k,ax+3);
	}
	tw::Float dL(const tw::Int& x,const tw::Int& y,const tw::Int& z,const tw::Int& ax) const
	{
		// returns arc length between 2 cell centers adjacent to this cell center, along axis=ax
		return dl(x,y,z,ax) + dl(x+kdelta(1,ax),y,z+kdelta(3,ax),ax);
	}
	tw::Float dL(const tw::cell& cell,const tw::Int& ax) const
	{
		return dL(cell.dcd1(),cell.dcd2(),cell.dcd3(),ax);
	}
	tw::Float dL(const tw::xstrip<3>& v,const tw::Int& k,const tw::Int& ax) const
	{
		return dL(v.dcd1(k),v.dcd2(k),k,ax);
	}
	void GetCellMetrics(const tw::cell& cell,const tw::Int& ax,tw::Float *dV,tw::Float *dS0,tw::Float *dS1,tw::Float *dl0,tw::Float *dl1) const
	{
		const tw::Int i1 = cell.dcd1() - mlb[1];
		const tw::Int i3 = cell.dcd3() - mlb[3];
		const tw::Int d1 = kdelta(1,ax);
		const tw::Int d3 = kdelta(3,ax);
		*dV = cell_area_x[i1] * cell_area_z[i3];
		*dS0 = cell_area_x[ax*mnum[1] + i1] * cell_area_z[ax*mnum[3] + i3];
		*dS1 = cell_area_x[ax*mnum[1] + i1 + d1] * cell_area_z[ax*mnum[3] + i3 + d3];
		*dl0 = cell_arc_x[(ax-1)*mnum[1] + i1] * cell_arc_z[(ax-1)*mnum[3] + i3];
		*dl1 = cell_arc_x[(ax-1)*mnum[1] + i1 + d1] * cell_arc_z[(ax-1)*mnum[3] + i3 + d3];
	}
	void CurvilinearToCartesian(tw::vec3 *r) const;
	void CartesianToCurvilinear(tw::vec3 *r) const;
	void CurvilinearToSpherical(tw::vec3 *r) const;
	void CurvilinearToCylindrical(tw::vec3 *r) const;
	void GetTangentVectorBasis(tw::basis *xfrm,const tw::vec3& r) const; // r expressed in curvilinear coordinates
	void TangentVectorToCartesian(tw::vec3 *v,const tw::vec3& r) const; // r expressed in curvilinear coordinates
	void TangentVectorToCurvilinear(tw::vec3 *v,const tw::vec3& r) const; // r expressed in curvilinear coordinates
	void LaplacianParameters(const tw::Int& a,const tw::Int& x,const tw::Int& y,const tw::Int& z,tw::Float *D1,tw::Float *D2,tw::Float *l1,tw::Float *l2) const;
};

inline tw::vec3 MetricSpace::ScaleFactor(const tw::vec3& r) const
{
	return tw::vec3(
				1.0,
				car + cyl*r.x + sph*r.x*sin(r.z),
				car + cyl + sph*r.x
				);
}

inline tw::Float MetricSpace::CylindricalRadius(const tw::vec3& r) const
{
	return sqrt((car+cyl)*sqr(r.x) + car*sqr(r.y)) + sph*r.x*sin(r.z);
}

inline tw::Float MetricSpace::SphericalRadius(const tw::vec3& r) const
{
	return sqrt(sqr(r.x) + car*sqr(r.y) + (car+cyl)*sqr(r.z));
}

inline void MetricSpace::CurvilinearToCartesian(tw::vec3 *r) const
{
	tw::vec3 temp = *r;
	const tw::Float cy = cos(temp.y); const tw::Float sy = sin(temp.y);
	const tw::Float cz = cos(temp.z); const tw::Float sz = sin(temp.z);
	r->x = car*temp.x + cyl*temp.x*cy + sph*temp.x*sz*cy;
	r->y = car*temp.y + cyl*temp.x*sy + sph*temp.x*sz*sy;
	r->z = car*temp.z + cyl*temp.z + sph*temp.x*cz;
}

inline void MetricSpace::CartesianToCurvilinear(tw::vec3 *r) const
{
	tw::vec3 temp = *r;
	const tw::Float rho = sqrt(sqr(temp.x) + sqr(temp.y));
	const tw::Float R = sqrt(rho*rho + sqr(temp.z));
	const tw::Float phi = atan2(temp.y,temp.x);
	const tw::Float theta = atan2(rho,temp.z);
	r->x = car*temp.x + cyl*rho + sph*R;
	r->y = car*temp.y + cyl*phi + sph*phi;
	r->z = car*temp.z + cyl*temp.z + sph*theta;
}

inline void MetricSpace::CurvilinearToCylindrical(tw::vec3 *r) const
{
	tw::vec3 temp = *r;
	const tw::Float rho = sqrt(sqr(temp.x) + sqr(temp.y));
	const tw::Float phi = atan2(temp.y,temp.x);
	r->x = car*rho + cyl*temp.x + sph*temp.x*sin(temp.z);
	r->y = car*phi + cyl*temp.y + sph*temp.y;
	r->z = car*temp.z + cyl*temp.z + sph*temp.x*cos(temp.z);
}

inline void MetricSpace::CurvilinearToSpherical(tw::vec3 *r) const
{
	tw::vec3 temp = *r;
	const tw::Float rho = car*sqrt(sqr(temp.x) + sqr(temp.y)) + cyl*temp.x;
	const tw::Float R = sqrt(rho*rho + temp.z*temp.z);
	r->x = (car+cyl)*R + sph*temp.x;
	r->y = car*atan2(temp.y,temp.x) + cyl*temp.y + sph*temp.y;
	r->z = (car+cyl)*atan2(rho,temp.z) + sph*temp.z;
}

inline void MetricSpace::GetTangentVectorBasis(tw::basis *b,const tw::vec3& r) const
{
	const tw::Float cy = cos(r.y); const tw::Float sy = sin(r.y);
	const tw::Float cz = cos(r.z); const tw::Float sz = sin(r.z);
	b->u.x = car + cyl*cy + sph*sz*cy;
	b->u.y = cyl*sy + sph*sz*sy;
	b->u.z = sph*cz;
	b->v.x = -cyl*sy - sph*sy;
	b->v.y = car + cyl*cy + sph*cy;
	b->v.z = 0.0;
	b->w.x = sph*cz*cy;
	b->w.y = sph*cz*sy;
	b->w.z = car + cyl - sph*sz;
}

inline void MetricSpace::TangentVectorToCartesian(tw::vec3 *v,const tw::vec3& r) const
{
	tw::basis b;
	GetTangentVectorBasis(&b,r);
	b.ExpressInStdBasis(v);
}

inline void MetricSpace::TangentVectorToCurvilinear(tw::vec3 *v,const tw::vec3& r) const
{
	tw::basis b;
	GetTangentVectorBasis(&b,r);
	b.ExpressInBasis(v);
}
