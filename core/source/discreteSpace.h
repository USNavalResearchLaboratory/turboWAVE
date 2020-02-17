enum tw_geometry {cartesian,cylindrical,spherical};
enum tw_boundary_spec {cyclic,reflecting,absorbing,emitting,axisymmetric};

struct Primitive
{
	float x[3];	// x[0],x[1],x[2] = relative position in cell using interval [-0.5,0.5); linearly mapped to curvilinear coordinates in cell
	tw::Int cell; // encoded index of the reference cell
	Primitive() noexcept
	{
		cell = 0;
		x[0] = x[1] = x[2] = 0.0;
	}
	Primitive(tw::Int c,float x0,float x1,float x2) noexcept
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

	tw::Float dt,dth,dti; // timestep, half timestep, inverse timestep
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
	// The following are the indices of the lower and upper far ghost cells
	tw::Int lfg[4],ufg[4];
	// The following are the indices of the lower and upper near ghost cells
	tw::Int lng[4],ung[4];
	// The following is the number of layers of ghost cells, layers[0] holds the maximum layers
	tw::Int layers[4];
	// Identify ignorable dimensions, 1 if axis is ignorable, 0 otherwise
	tw::Int ignorable[4];

	public:

	DiscreteSpace();
	DiscreteSpace(tw::Float dt,tw::Int x,tw::Int y,tw::Int z,const tw::vec3& corner,const tw::vec3& size,tw::Int ghostCellLayers);
	void Resize(tw::Int x,tw::Int y,tw::Int z,const tw::vec3& corner,const tw::vec3& size,tw::Int ghostCellLayers);
	void Resize(tw::Int x,tw::Int y,tw::Int z,const tw::vec3& corner,const tw::vec3& size);
	void SetupTimeInfo(tw::Float dt0) { dt = dt0; dth = 0.5*dt0; dti = 1.0/dt0; }
	tw::Int EncodeCell(tw::Int i,tw::Int j,tw::Int k) const { return (i-lfg[1])*encodingStride[1] + (j-lfg[2])*encodingStride[2] + (k-lfg[3])*encodingStride[3]; }
	void DecodeCell(const tw::Int& cell,tw::Int ijk[4]) const
	{
		// Assumes z-packed encoding pattern, i.e. decodingStride[3]=1
		ijk[1] = lfg[1] + cell / decodingStride[1];
		ijk[2] = lfg[2] + (cell % decodingStride[1]) / decodingStride[2];
		ijk[3] = lfg[3] + (cell % decodingStride[1]) % decodingStride[2];
	}
	void DecodeCell(const tw::Int& cell,tw::Int *i,tw::Int *j,tw::Int *k) const
	{
		*i = lfg[1] + cell / decodingStride[1];
		*j = lfg[2] + (cell % decodingStride[1]) / decodingStride[2];
		*k = lfg[3] + (cell % decodingStride[1]) % decodingStride[2];
	}
	void DecodeCell(const Primitive& q,tw::Int ijk[4]) const { DecodeCell(q.cell,ijk); }
	void DecodeCell(const Primitive& q,tw::Int *i,tw::Int *j,tw::Int *k) const { DecodeCell(q.cell,i,j,k); }
	bool RefCellInDomain(const Primitive& q) const;
	void SetPrimitiveWithPosition(Primitive& q,const tw::vec3& pos) const;
	void UpdatePrimitiveWithPosition(Primitive& q,const tw::vec3& pos) const;
	tw::vec3 PositionFromPrimitive(const Primitive& q) const;
	void MinimizePrimitive(Primitive& q) const;
	void MinimizePrimitive(tw::Int *cell,tw::Int ijk[3][tw::max_bundle_size],float x[3][tw::max_bundle_size],float domainMask[tw::max_bundle_size]) const;
	//void MinimizePrimitiveBohm(tw::Int *cell,tw::Int ijk[3][tw::max_bundle_size],float x[3][tw::max_bundle_size],float domainMask[tw::max_bundle_size]) const;
	tw::Int Dim(const tw::Int& ax) const { return dim[ax]; }
	tw::Int Num(const tw::Int& ax) const { return num[ax]; }
	tw::Int Layers(const tw::Int& ax) const { return layers[ax]; }
	tw::Int LNG(const tw::Int& ax) const { return lng[ax]; }
	tw::Int UNG(const tw::Int& ax) const { return ung[ax]; }
	tw::Int LFG(const tw::Int& ax) const { return lfg[ax]; }
	tw::Int UFG(const tw::Int& ax) const { return ufg[ax]; }
	friend tw::Float timestep(const DiscreteSpace& A) { return A.dt; }
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
	return (dim[1]==1 || (i>=lfg[1]+inset && i<=ufg[1]-inset))
		&& (dim[2]==1 || (j>=lfg[2]+inset && j<=ufg[2]-inset))
		&& (dim[3]==1 || (k>=lfg[3]+inset && k<=ufg[3]-inset));
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

inline bool DiscreteSpace::RefCellInDomain(const Primitive& q) const
{
	tw::Int ijk[4];
	DecodeCell(q,ijk);
	return	((ijk[1]>=1 && ijk[1]<=dim[1]) || dim[1]==1) &&
			((ijk[2]>=1 && ijk[2]<=dim[2]) || dim[2]==1) &&
			((ijk[3]>=1 && ijk[3]<=dim[3]) || dim[3]==1);
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
				itest[par] = tw::Int(x[ax][par]<-0.5 && ijk[ax][par]>lfg[ax+1]) - tw::Int(x[ax][par]>=0.5 && ijk[ax][par]<ufg[ax+1]);
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
				repeat += tw::Int(x[ax][par]<-0.5 && ijk[ax][par]>lfg[ax+1]) + tw::Int(x[ax][par]>=0.5 && ijk[ax][par]<ufg[ax+1]);
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
			test = tw::Int(q.x[ax-1]<-0.5 && ijk[ax]>lfg[ax]);
			test -= tw::Int(q.x[ax-1]>=0.5 && ijk[ax]<ufg[ax]);
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
