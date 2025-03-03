module;

#include "tw_includes.h"
#include "tw_test.h"

export module discrete_space;
import base;
export import tensor;
export import tasks;

export namespace tw
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
		enum axis { t,x,y,z,mass,px,py,pz,g,gbx,gby,gbz,Qp }; // (x,y,z) map to (r,phi,z) or (r,phi,theta) in cases of curvilinear geometry
		enum side { low , high };
		inline std::map<std::string,axis> axis_map()
		{
			return {{"t",t},{"x",x},{"y",y},{"z",z},{"mass",mass},{"px",px},{"py",py},{"pz",pz},{"g",g},{"gbx",gbx},{"gby",gby},{"gbz",gbz},{"Qp",Qp}};
		}
		inline std::string pretty_axis_label(const tw::grid::axis& axis)
		{
			std::map<tw::grid::axis,std::string> m = {{t,"$t$"},{x,"$x$"},{y,"$y$"},{z,"$z$"},{mass,"$\\gamma m$"},{px,"$p_x$"},{py,"$p_y$"},{pz,"$p_z$"},{g,"$\\gamma$"},{gbx,"$\\gamma\\beta_x$"},{gby,"$\\gamma\\beta_y$"},{gbz,"$\\gamma\\beta_z$"},{Qp,"$\\chi$"}};
			return m[axis];
		}
		inline tw::Int naxis(const tw::grid::axis& axis)
		{
			std::map<tw::grid::axis,tw::Int> M = {{t,0},{x,1},{y,2},{z,3},{mass,4},{px,5},{py,6},{pz,7},{g,8},{gbx,9},{gby,10},{gbz,11},{Qp,12}};
			return M[axis];
		}

		inline axis enumaxis(tw::Int ax)
		{
			std::map<tw::Int,axis> M = {{0,y},{1,x},{2,y},{3,z},{4,mass},{5,px},{6,py},{7,pz},{8,g},{9,gbx},{10,gby},{11,gbz},{12,Qp}};
			return M[ax];
		}
	}
}

/// Abstraction for location on a grid
export struct Primitive
{
	/// The reference cell encoded as a single integer
	tw::Int cell;
	/// This is the relative position in a space-time cell, referenced to the interval [-0.5,0.5).
	/// The components can represent arbitrary coordinates.
	/// The time component x[0] is not involved with cell encoding.
	float x[4];
	Primitive() noexcept
	{
		cell = 0;
		x[0] = x[1] = x[2] = x[3] = 0.0;
	}
	Primitive(tw::Int c,float x0,float x1,float x2,float x3) noexcept
	{
		cell = c;
		x[0] = x0;
		x[1] = x1;
		x[2] = x2;
		x[3] = x3;
	}
	friend std::ostream& operator << (std::ostream& os, const Primitive& q)
	{
		os << '<' << q.cell << ": (" << q.x[0] << ',' << q.x[1] << ',' << q.x[2] << ',' << q.x[3] << ")>";
		return os;
	}
};

/// Data describing any kind of particle
export struct Particle
{
	float number; ///< particles per macroparticle divided by \f$n_0(c/wp)^3\f$
	Primitive q; ///< abstraction for the spatial coordinate
	tw::vec4 p; ///< momentum , always known in Cartesian coordinates
	tw::vec4 s; ///< polarization, always known in Cartesian coordinates
	uint64_t tag; ///< unique identifier, low 32 bits is the node of origin
	tw::Float Qparam; //< quantum parameter

	/// Constructor, parameters shadow the member variables
	Particle(const float number,const Primitive& q,const tw::vec4& p,const tw::vec4& s,const uint64_t tag,const tw::Float& Qparam) noexcept;
	void ReadCheckpoint(std::ifstream& inFile);
	void WriteCheckpoint(std::ofstream& outFile);

	/// Used to define the ordering of particles for `std::sort`
	friend bool operator < (const Particle& p1,const Particle& p2)
	{
		return p1.q.cell < p2.q.cell;
	}
};

/// Used to pack data for message passing of particles
/// The destination information is computed on the source node and packaged with
/// the particle.  No floating point operations other than copying should be needed.
export struct TransferParticle
{
	/// `dst[0]` is rank of starting domain upon construction; gets set to destination domain later.
	/// `dst[1..3]` are -1, 0, or 1, giving direction of movement or no movement.
	tw::Int dst[4];
	float number; ///< particles per macroparticle divided by \f$n_0(c/wp)^3\f$
	tw::Int ijk[4]; ///< topological indices referenced on the source node
	float x[4]; ///< for transfers, the relative cell position can be kept without change
	tw::vec4 p; ///< for tansfers, momentum can be kept unchanged
	tw::vec4 s; ///< for transfers, polarization can be kept unchanged
	uint64_t tag; ///< for transfers, tag can be kept unchanged
	tw::Float Qparam; //< for transfers, quantum parameter can be kept unchanged
};

/// Used to create sorting map within a thread for subsets of particle lists
export struct ParticleRef
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
export struct weights_3D
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
export struct DiscreteSpace : Testable
{
	protected:

	tw::vec4 corner; ///< position where all coordinates are minimized on the local domain
	tw::vec4 size; ///< length of the local domain along each axis
	tw::vec4 globalCorner; ///< position where all coordinates are minimized on the global domain
	tw::vec4 globalSize; ///< length of the global domain along each axis
	tw::vec4 spacing; ///< center-to-center cell separation along each axis (uniform grids only)
	tw::vec4 freq; ///< inverse of spacing component by component (uniform grids only)
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
	tw::Int layers[4]; ///< number of ghost cell layers, assumed all the same some places
	tw::Int ignorable[4]; ///< 1 if axis is ignorable, 0 otherwise

	public:

	/// Create an empty `DiscreteSpace`
	DiscreteSpace();
	/// Create a `DiscreteSpace` with purely *local* coordinates.
	DiscreteSpace(tw::Int xDim,tw::Int yDim,tw::Int zDim,const tw::vec4& corner,const tw::vec4& size,tw::Int ghostCellLayers=2);
	/// Change the topology and coordinates.
	void Resize(const tw::Int dim[4],const tw::Int gdim[4],const tw::Int dom[4],const tw::vec4& gcorner,const tw::vec4& gsize,tw::Int ghostCellLayers=2);
	/// Change the coordinates, inheriting the `Task` topology.
	void Resize(Task& task,const tw::vec4& gcorner,const tw::vec4& gsize,tw::Int ghostCellLayers=2);
	/// Change the time step.  Use `Simulation::UpdateTimestep` to do this for all modules.
	void SetupTimeInfo(tw::Float dt0) { spacing[0] = dt0; freq[0] = 1.0/dt0; }
	/// Encode the cell with topological indices `(n,i,j,k)`
	tw::Int EncodeCell(tw::Int n,tw::Int i,tw::Int j,tw::Int k) const {
		return (n-lfg[0])*encodingStride[0] + (i-lfg[1])*encodingStride[1] + (j-lfg[2])*encodingStride[2] + (k-lfg[3])*encodingStride[3];
	}
	/// Decode `cell` to produce topological indices `(ijk[0..4])`. This assumes z-packed encoding, i.e., `decodingStride[3]=1`.
	void DecodeCell(const tw::Int& cell,tw::Int ijk[4]) const
	{
		ijk[0] = lfg[0] + cell / decodingStride[0];
		ijk[1] = lfg[1] + (cell % decodingStride[0]) / decodingStride[1];
		ijk[2] = lfg[2] + (cell % decodingStride[1]) / decodingStride[2];
		ijk[3] = lfg[3] + (cell % decodingStride[2]) / decodingStride[3];
	}
	/// Decode `cell` to produce topological indices `(n,i,j,k)`. This assumes z-packed encoding, i.e., `decodingStride[3]=1`.
	void DecodeCell(const tw::Int& cell,tw::Int *n,tw::Int *i,tw::Int *j,tw::Int *k) const
	{
		*n = lfg[0] + cell / decodingStride[0];
		*i = lfg[1] + (cell % decodingStride[0]) / decodingStride[1];
		*j = lfg[2] + (cell % decodingStride[1]) / decodingStride[2];
		*k = lfg[3] + (cell % decodingStride[2]) / decodingStride[3];
	}
	/// Decode `q` to produce topological indices `(ijk[0..4])`. This assumes z-packed encoding, i.e., `decodingStride[3]=1`.
	void DecodeCell(const Primitive& q,tw::Int ijk[4]) const { DecodeCell(q.cell,ijk); }
	/// Decode `q` to produce topological indices `(n,i,j,k)`. This assumes z-packed encoding, i.e., `decodingStride[3]=1`.
	void DecodeCell(const Primitive& q,tw::Int *n,tw::Int *i,tw::Int *j,tw::Int *k) const { DecodeCell(q.cell,n,i,j,k); }
	bool RefCellInDomain(const Primitive& q) const;
	void SetPrimitiveWithPosition(Primitive& q,const tw::vec4& pos) const;
	tw::vec4 PositionFromPrimitive(const Primitive& q) const;
	void MinimizePrimitive(Primitive& q) const;
	void MinimizePrimitive(tw::Int cell[tw::max_bundle_size],tw::Int ijk[4][tw::max_bundle_size],float x[4][tw::max_bundle_size],float domainMask[tw::max_bundle_size]) const;
	tw::Int Dim(const tw::Int& ax) const { return dim[ax]; }
	tw::Int Num(const tw::Int& ax) const { return num[ax]; }
	tw::Int Ignorable(const tw::Int& ax) const { return ignorable[ax]; }
	tw::Int Layers(const tw::Int& ax) const { return layers[ax]; }
	tw::Int LNG(const tw::Int& ax) const { return lng[ax]; }
	tw::Int UNG(const tw::Int& ax) const { return ung[ax]; }
	tw::Int LFG(const tw::Int& ax) const { return lfg[ax]; }
	tw::Int UFG(const tw::Int& ax) const { return ufg[ax]; }
	tw::Float dx(const tw::Int& ax) const { return spacing[ax]; }
	tw::Float dk(const tw::Int& ax) const { return freq[ax]; }
	tw::vec4 PhysicalSize() { return size; }
	tw::vec4 GlobalPhysicalSize() { return globalSize; }
	tw::vec4 Corner() { return corner; }
	tw::vec4 GlobalCorner() { return globalCorner; }
	bool IsRefCellWithin(const Primitive& q,tw::Int inset) const;
	bool IsPointWithinInterior(const tw::vec4& P);
	tw::Int Dimensionality();
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
	const tw::Int N = tw::max_bundle_size;
	const tw::Int AB = tw::vec_align_bytes;
	alignas(AB) tw::Int pos[N],neg[N],displ[N],max_displ[N];
	tw::Int par,ax;
	#pragma omp simd aligned(cell,ijk:AB)
	for (par=0;par<N;par++)
	{
		DecodeCell(cell[par],&ijk[0][par],&ijk[1][par],&ijk[2][par],&ijk[3][par]);
		domainMask[par] = 1.0f;
	}
	for (ax=0;ax<4;ax++)
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
			domainMask[par] *= float(tw::Int(ijk[ax][par]>=1 && ijk[ax][par]<=dim[ax]));
		}
	}
	#pragma omp simd aligned(cell,ijk:AB)
	for (par=0;par<N;par++)
		cell[par] = EncodeCell(ijk[0][par],ijk[1][par],ijk[2][par],ijk[3][par]);
}

/// Change the reference cell to the one containing the particle centroid.
/// This leaves x[ax] within the normal range [-.5,.5) *except* when the particle
/// has left the extended domain.  In the latter case the reference cell is pegged to
/// the current extended domain, leaving x[ax] out of normal range.
inline void DiscreteSpace::MinimizePrimitive(Primitive& q) const
{
	tw::Int ax,ijk[4],displ,max_displ,pos,neg;
	const float eps = std::numeric_limits<float>::epsilon();
	DecodeCell(q.cell,ijk);
	for (ax=0;ax<4;ax++)
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

DiscreteSpace::DiscreteSpace()
{
	ignorable[0] = 0;
	ignorable[1] = 0;
	ignorable[2] = 0;
	ignorable[3] = 0;
}

DiscreteSpace::DiscreteSpace(tw::Int xDim,tw::Int yDim,tw::Int zDim,const tw::vec4& corner,const tw::vec4& size,tw::Int ghostCellLayers)
{
	const tw::Int ldim[4] = { 1, xDim, yDim, zDim };
	const tw::Int gdim[4] = { 1, xDim, yDim, zDim };
	const tw::Int dom[4] = { 0, 0, 0, 0 };
	Resize(ldim,gdim,dom,corner,size,ghostCellLayers);
}

// DiscreteSpace::DiscreteSpace(tw::Float dt0,Task& task,const tw::vec3& gcorner,const tw::vec3& gsize,tw::Int ghostCellLayers)
// {
// 	SetupTimeInfo(dt0);
// 	Resize(task,gcorner,gsize,ghostCellLayers);
// }

void DiscreteSpace::Resize(Task& task,const tw::vec4& gcorner,const tw::vec4& gsize,tw::Int ghostCellLayers)
{
	Resize(task.localCells,task.globalCells,task.domainIndex,gcorner,gsize,ghostCellLayers);
}

void DiscreteSpace::Resize(const tw::Int dim[4],const tw::Int gdim[4],const tw::Int dom[4],const tw::vec4& gcorner,const tw::vec4& gsize,tw::Int ghostCellLayers)
{
	this->dim[0] = dim[0];
	this->dim[1] = dim[1];
	this->dim[2] = dim[2];
	this->dim[3] = dim[3];

	for (tw::Int i=0;i<4;i++)
	{
		if (dim[i]==1)
		{
			layers[i] = 0;
			lfg[i] = 1;
			ufg[i] = 1;
			lng[i] = 1;
			ung[i] = 1;
		}
		else
		{
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

	for (tw::Int i=0;i<4;i++)
	{
		encodingStride[i] = (dim[i]==1 ? 0 : decodingStride[i]);
		ignorable[i] = (dim[i]==1 ? 1 : 0);
	}

	globalCorner = gcorner;
	globalSize = gsize;
	for (tw::Int i=0;i<4;i++)
	{
		spacing[i] = globalSize[i]/gdim[i];
		freq[i] = 1/spacing[i];
		size[i] = dim[i]*spacing[i];
		corner[i] = gcorner[i] + dom[i]*size[i];
	}
}

void DiscreteSpace::ReadCheckpoint(std::ifstream& inFile)
{
	inFile.read((char*)&corner,sizeof(tw::vec4));
	inFile.read((char*)&size,sizeof(tw::vec4));
	inFile.read((char*)&globalCorner,sizeof(tw::vec4));
	inFile.read((char*)&globalSize,sizeof(tw::vec4));
	inFile.read((char*)&spacing,sizeof(tw::vec4));
	inFile.read((char*)&freq,sizeof(tw::vec4));
	inFile.read((char*)num,sizeof(num));
	inFile.read((char*)dim,sizeof(dim));
	inFile.read((char*)lfg,sizeof(lfg));
	inFile.read((char*)ufg,sizeof(ufg));
	inFile.read((char*)lng,sizeof(lng));
	inFile.read((char*)ung,sizeof(ung));
	inFile.read((char*)ignorable,sizeof(ignorable));
	inFile.read((char*)encodingStride,sizeof(encodingStride));
	inFile.read((char*)decodingStride,sizeof(decodingStride));
	inFile.read((char*)layers,sizeof(layers));
}

void DiscreteSpace::WriteCheckpoint(std::ofstream& outFile)
{
	outFile.write((char*)&corner,sizeof(tw::vec4));
	outFile.write((char*)&size,sizeof(tw::vec4));
	outFile.write((char*)&globalCorner,sizeof(tw::vec4));
	outFile.write((char*)&globalSize,sizeof(tw::vec4));
	outFile.write((char*)&spacing,sizeof(tw::vec4));
	outFile.write((char*)&freq,sizeof(tw::vec4));
	outFile.write((char*)num,sizeof(num));
	outFile.write((char*)dim,sizeof(dim));
	outFile.write((char*)lfg,sizeof(lfg));
	outFile.write((char*)ufg,sizeof(ufg));
	outFile.write((char*)lng,sizeof(lng));
	outFile.write((char*)ung,sizeof(ung));
	outFile.write((char*)ignorable,sizeof(ignorable));
	outFile.write((char*)encodingStride,sizeof(encodingStride));
	outFile.write((char*)decodingStride,sizeof(decodingStride));
	outFile.write((char*)layers,sizeof(layers));
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

////////////////////
// PARTICLE CLASS //
////////////////////


Particle::Particle(const float number,const Primitive& q,const tw::vec4& p,const tw::vec4& s,const uint64_t tag,const tw::Float& Qparam) noexcept
{
	this->number = number;
	this->q = q;
	this->p = p;
	this->s = s;
	this->tag = tag;
	this->Qparam = Qparam;
}

void Particle::ReadCheckpoint(std::ifstream& inFile)
{
	inFile.read((char *)this,sizeof(Particle));
}

void Particle::WriteCheckpoint(std::ofstream& outFile)
{
	outFile.write((char *)this,sizeof(Particle));
}
