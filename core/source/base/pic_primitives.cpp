module;

#include "tw_includes.h"

export module pic_primitives;
import base;
import tensor;

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
	Particle(const float number,const Primitive& q,const tw::vec4& p,const tw::vec4& s,const uint64_t tag,const tw::Float& Qparam) noexcept {
        this->number = number;
        this->q = q;
        this->p = p;
        this->s = s;
        this->tag = tag;
        this->Qparam = Qparam;
    }
	void ReadCheckpoint(std::ifstream& inFile) {
	    inFile.read((char *)this,sizeof(Particle));
    }
	void WriteCheckpoint(std::ofstream& outFile) {
        outFile.write((char *)this,sizeof(Particle));
    }

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
