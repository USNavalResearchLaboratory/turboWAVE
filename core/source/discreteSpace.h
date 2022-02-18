namespace tw
{
	/// boundary conditions
	namespace bc
	{
		/// boundary conditions for particles
		enum class par {none,periodic,reflecting,absorbing,emitting,axisymmetric};
		/// Maps input file identifiers to particle B.C. enumeration identifiers
		inline std::map<std::string,par> par_map()
		{
			return {{"periodic",par::periodic},{"reflecting",par::reflecting},{"absorbing",par::absorbing},{"open",par::absorbing},{"emitting",par::emitting},{"axisymmetric",par::axisymmetric}};
		}
		/// boundary conditions for fields
		enum class fld	{none,periodic,normalFluxFixed,dirichletWall,neumannWall,dirichletCell,natural};
		/// Maps input file identifiers to field B.C. enumeration identifiers
		inline std::map<std::string,fld> fld_map()
		{
			return {{"periodic",fld::periodic},{"neumann",fld::neumannWall},{"dirichlet",fld::dirichletCell},{"open",fld::natural}};
		}
	}
	/// Types and maps for working with grid axes
	namespace grid
	{
		enum geometry {cartesian,cylindrical,spherical};
		enum axis { t,x,y,z,mass,px,py,pz,g,gbx,gby,gbz }; // (x,y,z) map to (r,phi,z) or (r,phi,theta) in cases of curvilinear geometry
		enum side { low , high };
		inline std::map<std::string,axis> axis_map()
		{
			return {{"t",t},{"x",x},{"y",y},{"z",z},{"mass",mass},{"px",px},{"py",py},{"pz",pz},{"g",g},{"gbx",gbx},{"gby",gby},{"gbz",gbz}};
		}
		inline std::string pretty_axis_label(const tw::grid::axis& axis)
		{
			std::map<tw::grid::axis,std::string> m = {{t,"$t$"},{x,"$x$"},{y,"$y$"},{z,"$z$"},{mass,"$\\gamma m$"},{px,"$p_x$"},{py,"$p_y$"},{pz,"$p_z$"},{g,"$\\gamma$"},{gbx,"$\\gamma\\beta_x$"},{gby,"$\\gamma\\beta_y$"},{gbz,"$\\gamma\\beta_z$"}};
			return m[axis];
		}
		inline tw::Int naxis(const tw::grid::axis& axis)
		{
			std::map<tw::grid::axis,tw::Int> M = {{t,0},{x,1},{y,2},{z,3},{mass,4},{px,5},{py,6},{pz,7},{g,8},{gbx,9},{gby,10},{gbz,11}};
			return M[axis];
		}

		inline axis enumaxis(tw::Int ax)
		{
			std::map<tw::Int,axis> M = {{0,y},{1,x},{2,y},{3,z},{4,mass},{5,px},{6,py},{7,pz},{8,g},{9,gbx},{10,gby},{11,gbz}};
			return M[ax];
		}
	}
}

/// Abstraction for location on a grid
struct Primitive
{
	/// This is the relative position in a cell referenced to the interval [-0.5,0.5).
	/// The components can represent arbitrary coordinates.
	float x[3];
	tw::Int cell; ///< encoded index of the reference cell
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
	friend std::ostream& operator << (std::ostream& os, const Primitive& q)
	{
		os << '<' << q.cell << ": (" << q.x[0] << ',' << q.x[1] << ',' << q.x[2] << ")>";
		return os;
	}
};

/// Data describing any kind of particle
struct Particle
{
	Primitive q; ///< abstraction for the spatial coordinate
	tw::vec3 p; ///< momentum , always known in Cartesian coordinates
	float number; ///< number = particles per macroparticle divided by \f$n_0(c/wp)^3\f$
	float aux1; ///< auxiliary, typically node of origin
	float aux2; ///< auxiliary, typically particle tag

	/// Constructor, parameters shadow the member variables
	Particle(const tw::vec3& p,const Primitive& q,const float number,const float aux1,const float aux2) noexcept;
	/// Assertions verifying that the primitive is minimized and in bounds
	void Assert(tw::Int beg,tw::Int end)
	{
		const float max_displ = 0.5f + std::numeric_limits<float>::epsilon();
		assert(q.cell>=beg && q.cell<end);
		assert(fabs(q.x[0]) < max_displ);
		assert(fabs(q.x[1]) < max_displ);
		assert(fabs(q.x[2]) < max_displ);
	}

	void ReadCheckpoint(std::ifstream& inFile);
	void WriteCheckpoint(std::ofstream& outFile);

	/// Used to define the ordering of particles for `std::sort`
	friend bool operator < (const Particle& p1,const Particle& p2)
	{
		return p1.q.cell < p2.q.cell;
	}
};

/// Used to pack data for message passing of particles
/// To avoid inconsistencies arising from FP comparisons on different nodes
/// the destination information is computed on the source domain and packaged with the particle
/// The position is kept as a double precision global coordinate until final call to `AddParticle`
struct TransferParticle
{
	/// `dst[0]` is rank of starting domain upon construction; gets set to destination domain later.
	/// `dst[1..3]` are -1, 0, or 1, giving direction of movement or no movement.
	tw::Int dst[4];
	tw::vec4 x,p; // can use x[0] or p[0] to pack extra info
	float number,aux1,aux2;
};

/// Used to create sorting map within a thread for subsets of particle lists
struct ParticleRef
{
	tw::Int idx,cell;
	ParticleRef() noexcept
	{
		idx = 0;
		cell = 0;
	}
	ParticleRef(tw::Int list_index,const Particle& par) noexcept
	{
		idx = list_index;
		cell = par.q.cell;
	}
	friend bool operator < (const ParticleRef& r1,const ParticleRef& r2)
	{
		return r1.cell < r2.cell;
	}
};

/// This holds the coefficients used to spread the particle cloud across
/// 3x3x3 grid cells.  Due to assumptions of separability
/// this only requires storing 3x3 coefficients.
struct weights_3D
{
	/// The weights are packed in a 3x3 matrix.  The first index (row) selects
	/// a cell from a 3-cell strip along a given axis, the second index (column)
	/// selects the axis, and the value is the weight factor in the cell.
	tw::Float w[3][3];
	tw::Int cell; ///< Encoded representation of the cell in which the particle center resides.
};

/// This object is an indexing and interpolation scheme for a structured grid.
/// Consider first the *topological indices*, defined as follows.
///
/// Let `ax` be an axis index, numbered from 1 to 3.  Consider a one dimensional
/// strip of cells lined up along axis `ax`.  Suppose this axis uses `L` ghost cell
/// layers.  Then, for any such strip,
/// - Interior cells are labeled `[1,...,dim[ax]]`
/// - Lower ghost cells are labeled `[1-L,...,0]`
/// - Upper ghost cells are labeled `[dim[ax]+1,...,dim[ax]+L]`.
///
/// The topological indices are independent of all storage patterns, In other words, no matter
/// what storage pattern is used, toplogical index `(i0,j0,k0)` will always refer
/// to the same cell.
///
/// Storage patterns are dictated by a *cell encoding*.  This associates the
/// topological indices with a single integer that maps to ascending memory addresses.
/// The default encoding is used to sort particles.
/// Superclasses can add a secondary encoding to resolve higher
/// dimensional storage patterns. For example, `Field` adds another encoding to
/// define how its components are packed.
struct DiscreteSpace
{
	protected:

	tw::Float dt; ///< timestep
	tw::Float dth; ///< half timestep
	tw::Float dti; ///< inverse timestep
	tw::vec3 corner; ///< position where all coordinates are minimized on the local domain
	tw::vec3 size; ///< length of the local domain along each axis
	tw::vec3 globalCorner; ///< position where all coordinates are minimized on the global domain
	tw::vec3 globalSize; ///< length of the global domain along each axis
	tw::vec3 spacing; ///< center-to-center cell separation along each axis (uniform grids only)
	tw::vec3 freq; ///< inverse of spacing component by component (uniform grids only)
	/// `dim[1..3]` are the number of cells along each axis.
	/// When a `Field` is derived, `dim[0]` becomes the number of components.
	tw::Int dim[4];
	/// `num` is similar to `dim`, except ghost cells are included.
	tw::Int num[4];
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
	tw::Int layers[4]; ///< number of ghost cell layers, `layers[0]` holds maximum layers
	tw::Int ignorable[4]; ///< 1 if axis is ignorable, 0 otherwise

	public:

	/// Create an empty `DiscreteSpace`
	DiscreteSpace();
	/// Create a `DiscreteSpace` with purely *local* coordinates.
	DiscreteSpace(tw::Int xDim,tw::Int yDim,tw::Int zDim,const tw::vec3& corner,const tw::vec3& size,tw::Int ghostCellLayers=2);
	/// Change the topology and coordinates.
	void Resize(const tw::Int dim[4],const tw::Int gdim[4],const tw::Int dom[4],const tw::vec3& gcorner,const tw::vec3& gsize,tw::Int ghostCellLayers=2);
	/// Change the coordinates, inheriting the `Task` topology.
	void Resize(Task& task,const tw::vec3& gcorner,const tw::vec3& gsize,tw::Int ghostCellLayers=2);
	/// Change the time step.  Use `Simulation::UpdateTimestep` to do this for all modules.
	void SetupTimeInfo(tw::Float dt0) { dt = dt0; dth = 0.5*dt0; dti = 1.0/dt0; }
	/// Encode the cell with topological indices `(i,j,k)`
	tw::Int EncodeCell(tw::Int i,tw::Int j,tw::Int k) const { return (i-lfg[1])*encodingStride[1] + (j-lfg[2])*encodingStride[2] + (k-lfg[3])*encodingStride[3]; }
	/// Decode `cell` to produce topological indices `(ijk[1],ijk[2],ijk[3])`. This assumes z-packed encoding, i.e., `decodingStride[3]=1`.
	void DecodeCell(const tw::Int& cell,tw::Int ijk[4]) const
	{
		ijk[1] = lfg[1] + cell / decodingStride[1];
		ijk[2] = lfg[2] + (cell % decodingStride[1]) / decodingStride[2];
		ijk[3] = lfg[3] + (cell % decodingStride[1]) % decodingStride[2];
	}
	/// Decode `cell` to produce topological indices `(i,j,k)`. This assumes z-packed encoding, i.e., `decodingStride[3]=1`.
	void DecodeCell(const tw::Int& cell,tw::Int *i,tw::Int *j,tw::Int *k) const
	{
		*i = lfg[1] + cell / decodingStride[1];
		*j = lfg[2] + (cell % decodingStride[1]) / decodingStride[2];
		*k = lfg[3] + (cell % decodingStride[1]) % decodingStride[2];
	}
	/// Decode `q` to produce topological indices `(ijk[1],ijk[2],ijk[3])`. This assumes z-packed encoding, i.e., `decodingStride[3]=1`.
	void DecodeCell(const Primitive& q,tw::Int ijk[4]) const { DecodeCell(q.cell,ijk); }
	/// Decode `q` to produce topological indices `(i,j,k)`. This assumes z-packed encoding, i.e., `decodingStride[3]=1`.
	void DecodeCell(const Primitive& q,tw::Int *i,tw::Int *j,tw::Int *k) const { DecodeCell(q.cell,i,j,k); }
	bool RefCellInDomain(const Primitive& q) const;
	void SetPrimitiveWithPosition(Primitive& q,const tw::vec3& pos) const;
	void UpdatePrimitiveWithPosition(Primitive& q,const tw::vec3& pos) const;
	tw::vec3 PositionFromPrimitive(const Primitive& q) const;
	void MinimizePrimitive(Primitive& q) const;
	void MinimizePrimitive(tw::Int *cell,tw::Int ijk[4][tw::max_bundle_size],float x[4][tw::max_bundle_size],float domainMask[tw::max_bundle_size]) const;
	tw::Int Dim(const tw::Int& ax) const { return dim[ax]; }
	tw::Int Num(const tw::Int& ax) const { return num[ax]; }
	tw::Int Ignorable(const tw::Int& ax) const { return ignorable[ax]; }
	tw::Int Layers(const tw::Int& ax) const { return layers[ax]; }
	tw::Int LNG(const tw::Int& ax) const { return lng[ax]; }
	tw::Int UNG(const tw::Int& ax) const { return ung[ax]; }
	tw::Int LFG(const tw::Int& ax) const { return lfg[ax]; }
	tw::Int UFG(const tw::Int& ax) const { return ufg[ax]; }
	tw::Float dx0(const tw::Int& ax) const { return spacing[ax-1]; }
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

inline bool DiscreteSpace::IsWithinInset(const Primitive& q,tw::Int inset) const
{
	tw::Int ijk[4];
	DecodeCell(q,ijk);
	bool ans = true;
	for (tw::Int ax=1;ax<=3;ax++)
		ans = ans && (dim[ax]==1 || (ijk[ax]>=lfg[ax]+inset && ijk[ax]<=ufg[ax]-inset));
	return ans;
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

/// Bundle version of `MinimizePrimitive`.  The bundle version also sets the `domainMask` to 0
/// for particles out of the interior domain, 1 otherwise.  Also the topological indices are loaded for further use.
inline void DiscreteSpace::MinimizePrimitive(tw::Int cell[tw::max_bundle_size],tw::Int ijk[4][tw::max_bundle_size],float x[4][tw::max_bundle_size],float domainMask[tw::max_bundle_size]) const
{
	const float eps = std::numeric_limits<float>::epsilon();
	const tw::Int N = tw::max_bundle_size;
	const tw::Int AB = tw::vec_align_bytes;
	alignas(AB) tw::Int pos[N],neg[N],displ[N],max_displ[N];
	tw::Int par,ax;
	#pragma omp simd aligned(cell,ijk:AB)
	for (par=0;par<N;par++)
		DecodeCell(cell[par],&ijk[1][par],&ijk[2][par],&ijk[3][par]);
	for (ax=1;ax<=3;ax++)
	{
		#pragma omp simd aligned(ijk,x,domainMask,pos,neg,displ,max_displ:AB)
		for (par=0;par<N;par++)
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
			domainMask[par] = float(tw::Int(ijk[ax][par]>=1 && ijk[ax][par]<=dim[ax]));
		}
	}
	#pragma omp simd aligned(cell,ijk:AB)
	for (par=0;par<N;par++)
		cell[par] = EncodeCell(ijk[1][par],ijk[2][par],ijk[3][par]);
}

/// Change the reference cell to the one containing the particle centroid.
/// This leaves x[ax] within the normal range [-.5,.5) *except* when the particle
/// has left the domain.  In the latter case the reference cell is pegged to
/// the current domain, leaving x[ax] out of normal range.
inline void DiscreteSpace::MinimizePrimitive(Primitive& q) const
{
	tw::Int ax,ijk[4],displ,max_displ,pos,neg;
	const float eps = std::numeric_limits<float>::epsilon();
	DecodeCell(q.cell,ijk);
	for (ax=1;ax<=3;ax++)
	{
		const float x = q.x[ax-1];
		pos = x>0.0f;
		neg = x<0.0f;
		displ = fabs(x) + 0.5f - neg*eps;
		// max_displ will be 0 if axis is ignorable due to ufg=lfg=ijk=1
		max_displ = pos*(ufg[ax] - ijk[ax]) + neg*(ijk[ax] - lfg[ax]);
		displ = (displ>max_displ)*max_displ + (displ<=max_displ)*displ + ignorable[ax]*displ;
		displ = (pos-neg)*displ;
		ijk[ax] += displ * (1-ignorable[ax]);
		q.x[ax-1] -= displ;
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
