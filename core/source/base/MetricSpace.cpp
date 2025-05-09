module;

#include "tw_includes.h"

export module metric_space;
export import base;
export import tensor;
export import units;
export import tasks;
import discrete_space;
import tw_iterator;

// This abstract class is needed to anticipate ComputeTool::Warp.
export struct warp_base
{
	tw::grid::axis ax;
	virtual tw::Float AddedCellWidth(tw::Int globalCell) = 0;
};

export struct MetricSpace:DiscreteSpace
{
	tw::grid::geometry gridGeometry;
	tw::Float car,cyl,sph; // set variable corresponding to coordinate system to 1.0, all others to 0.0
	tw::Int mnum[4]; // num for metric arrays (see DiscreteSpace)
	tw::Int mlb[4],mub[4]; // lfg and ufg for metric arrays (see DiscreteSpace)
	tw::Int I3x3[3][3];
	bool adaptiveTimestep,adaptiveGrid;
	
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

	tw::UnitConverter units;
	std::vector<warp_base*> warps;

	#ifdef USE_OPENCL
	cl_mem metricsBuffer;
	cl_mem stripBuffer[4];
	#endif

	MetricSpace();
	~MetricSpace();

private:
	void SetTopology(Task *task,
		const tw::Int gdim[4],
		const tw::vec4& gcorner,
		const tw::vec4& gsize,
		tw::Int ghostCellLayers);
	void Allocate();
	void SetSpacings();
	void UpdateHulls(Task *task);
	void SetupPositionArrays();
	void SetCartesianGeometry();
	void SetCylindricalGeometry();
	void SetSphericalGeometry();
public:
	void Resize(Task *task,
		const tw::Int gdim[4],
		const tw::vec4& gcorner,
		const tw::vec4& gsize,
		tw::Int ghostCellLayers=2,
		tw::grid::geometry geo=tw::grid::cartesian);
	void ReadCheckpoint(std::ifstream& inFile);
	void WriteCheckpoint(std::ofstream& outFile);
	void AttachUnits(tw::units sys,tw::Float unitDensityCGS);
	void AttachWarp(warp_base *w);
	bool Test(tw::Int& id);

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

	/// 3x3 Kronecker delta function
	tw::Int kdelta(const tw::Int& ax1,const tw::Int& ax2) const
	{
		return I3x3[ax1-1][ax2-1];
	}
	/// position in parameter space
	tw::Float X(const tw::Int& i,const tw::Int& ax) const
	{
		return gpos[mnum[0]*ax-mlb[ax]+i];
	}
	/// cell size in parameter space (not an arc length)
	tw::Float dX(const tw::Int& i,const tw::Int& ax) const
	{
		return width[mnum[0]*ax-mlb[ax]+i];
	}
	/// position in parameter space
	tw::Float& X(const tw::Int& i,const tw::Int& ax)
	{
		return gpos[mnum[0]*ax-mlb[ax]+i];
	}
	/// cell size in parameter space (not an arc length)
	tw::Float& dX(const tw::Int& i,const tw::Int& ax)
	{
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
	/// returns wall area for ax = axis normal to wall.  ax = 0 returns cell volume.
	tw::Float dS(const tw::Int& x,const tw::Int& y,const tw::Int& z,const tw::Int& ax) const
	{
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
	/// returns arc length from cell center to cell center, along axis=ax, from low side
	tw::Float dl(const tw::Int& x,const tw::Int& y,const tw::Int& z,const tw::Int& ax) const
	{
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
	/// returns midplane area for ax = axis normal to area.  ax = 0 returns cell volume.
	tw::Float dSh(const tw::Int& x,const tw::Int& y,const tw::Int& z,const tw::Int& ax) const
	{
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
	/// returns arc length from cell wall to cell wall, along axis=ax, along low side edge
	tw::Float dlh(const tw::Int& x,const tw::Int& y,const tw::Int& z,const tw::Int& ax) const
	{
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
	/// returns arc length between 2 cell centers adjacent to this cell center, along axis=ax
	tw::Float dL(const tw::Int& x,const tw::Int& y,const tw::Int& z,const tw::Int& ax) const
	{
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

	tw::Float ToLab(const Evolution& evo,tw::Float zeta,tw::Float relativeTime);
	tw::Float ToLight(const Evolution& evo,tw::Float z,tw::Float relativeTime);
	template <class T,class U>
	U ValueOnLabGrid(const Evolution& evo,T& A,tw::strip s,tw::Int k,tw::Float relativeTime);
	template <class T,class U>
	U ValueOnLightGrid(const Evolution& evo,T& A,tw::strip s,tw::Int k,tw::Float relativeTime);
};

MetricSpace::MetricSpace()
{
	gridGeometry = tw::grid::cartesian;
	car = 1.0;
	cyl = sph = 0.0;
	adaptiveTimestep = false;
	adaptiveGrid = false;
	I3x3[0][0] = 1;
	I3x3[0][1] = 0;
	I3x3[0][2] = 0;
	I3x3[1][0] = 0;
	I3x3[1][1] = 1;
	I3x3[1][2] = 0;
	I3x3[2][0] = 0;
	I3x3[2][1] = 0;
	I3x3[2][2] = 1;
	units = tw::UnitConverter(tw::units::plasma,1e19);
	#ifdef USE_OPENCL
	metricsBuffer = NULL;
	#endif
}

MetricSpace::~MetricSpace()
{
	#ifdef USE_OPENCL
	if (metricsBuffer!=NULL)
	{
		clReleaseMemObject(metricsBuffer);
		clReleaseMemObject(stripBuffer[0]);
		clReleaseMemObject(stripBuffer[1]);
		clReleaseMemObject(stripBuffer[2]);
		clReleaseMemObject(stripBuffer[3]);
	}
	#endif
}

/// Fully initialize object, may involve message passing.
/// If there are warps they should be attached before calling.
/// The `task` passed as the first argument must itself be initialized.
/// The `gsize` should be the size of the hull assuming uniform spacing.
/// It will be adjusted automatically to account for warps.
void MetricSpace::Resize(Task *task,
	const tw::Int gdim[4],
	const tw::vec4& gcorner,
	const tw::vec4& gsize,
	tw::Int ghostCellLayers,
	tw::grid::geometry geo)
{
	SetTopology(task,gdim,gcorner,gsize,ghostCellLayers);
	Allocate();
	SetSpacings();
	UpdateHulls(task);
	switch (geo)
	{
		case tw::grid::cartesian:
			SetCartesianGeometry();
			break;
		case tw::grid::cylindrical:
			SetCylindricalGeometry();
			break;
		case tw::grid::spherical:
			SetSphericalGeometry();
			break;
	}
	#ifdef USE_OPENCL
	InitializeMetricsBuffer(task->context,dt);
	#endif
}

void MetricSpace::AttachUnits(tw::units sys,tw::Float unitDensityCGS)
{
	units = tw::UnitConverter(sys,unitDensityCGS);
}

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

/// Argument will typically point to a ComputeTool.
/// Do not update refCount since ultimately the Simulation
/// object, which is derived from this object, owns the tool.
void MetricSpace::AttachWarp(warp_base *w)
{
	warps.push_back(w);
}

#ifdef USE_OPENCL

void MetricSpace::InitializeMetricsBuffer(cl_context ctx,tw::Float dt)
{
	cl_int err;
	tw::Float metrics[12] = { dt , spacing.x , spacing.y , spacing.z , 0.0 , corner.x , corner.y , corner.z , car , cyl , sph , 0.0 };
	tw::Int tStrip[7] = { 0,0,0,0,0,0 }; // not generally used; can signal a special operation
	tw::Int xStrip[7] = { 1,1,0,0,encodingStride[1],dim[1] }; // which-axis , x-axis? , y-axis? , z-axis? , stride , dim (there is deliberate redundancy)
	tw::Int yStrip[7] = { 2,0,1,0,encodingStride[2],dim[2] };
	tw::Int zStrip[7] = { 3,0,0,1,encodingStride[3],dim[3] };

	metricsBuffer = clCreateBuffer(ctx,CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,sizeof(tw::Float)*12,metrics,&err);
	stripBuffer[0] = clCreateBuffer(ctx,CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,sizeof(tw::Int)*6,tStrip,&err);
	stripBuffer[1] = clCreateBuffer(ctx,CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,sizeof(tw::Int)*6,xStrip,&err);
	stripBuffer[2] = clCreateBuffer(ctx,CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,sizeof(tw::Int)*6,yStrip,&err);
	stripBuffer[3] = clCreateBuffer(ctx,CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,sizeof(tw::Int)*6,zStrip,&err);

	if (err)
		throw tw::FatalError("Device buffer failure in MetricSpace.");
}

void MetricSpace::StripUpdateProtocol(cl_kernel k,cl_command_queue q,tw::Int axis,tw::Int stripArgument)
{
	size_t global_work_range[3] = { (size_t)num[1],(size_t)num[2],(size_t)num[3] };
	global_work_range[axis-1] = 1;
	clSetKernelArg(k,stripArgument,sizeof(stripBuffer[axis]),&stripBuffer[axis]);
	clEnqueueNDRangeKernel(q,k,3,NULL,global_work_range,NULL,0,NULL,NULL);
	clFinish(q);
}

#endif

/// Setup topological information only (step 1)
void MetricSpace::SetTopology(Task *task,
	const tw::Int gdim[4],
	const tw::vec4& gcorner,
	const tw::vec4& gsize,
	tw::Int ghostCellLayers)
{
	DiscreteSpace::Resize(task,gdim,gcorner,gsize,ghostCellLayers);

	// Metric arrays have ghost cell layers even when dim=1.
	// Hence the topology of the metric data differs from that of other data

	tw::Int maxLayers = 0;
	for (auto i=0;i<4;i++)
		maxLayers = layers[i] > maxLayers ? layers[i] : maxLayers;

	mlb[1] = 1 - maxLayers;
	mlb[2] = 1 - maxLayers;
	mlb[3] = 1 - maxLayers;

	mub[1] = dim[1] + maxLayers;
	mub[2] = dim[2] + maxLayers;
	mub[3] = dim[3] + maxLayers;

	mnum[1] = dim[1]+2*maxLayers;
	mnum[2] = dim[2]+2*maxLayers;
	mnum[3] = dim[3]+2*maxLayers;

	mnum[0] = mnum[1];
	if (mnum[2]>mnum[0])
		mnum[0] = mnum[2];
	if (mnum[3]>mnum[0])
		mnum[0] = mnum[3];
}

/// Allocate space for the metric data (step 2).
/// When checkpointing this has to be called during the read.
void MetricSpace::Allocate()
{
	gpos.resize(4*mnum[0]); // node coordinates
	width.resize(4*mnum[0]); // parametric cell size (not an arc length in general)

	cell_area_x.resize(mnum[1]*8);
	cell_arc_x.resize(mnum[1]*6);
	cell_area_z.resize(mnum[3]*8);
	cell_arc_z.resize(mnum[3]*6);
}

/// Setup uniform spacings and add warps (step 3).
/// N.b. spacings are coordinates, not arc lengths.
void MetricSpace::SetSpacings()
{
	for (tw::Int ax=1;ax<=3;ax++)
		for (tw::Int i=mlb[ax];i<=mub[ax];i++)
			dX(i,ax) = spacing[ax];
	for (auto warp : warps)
	{
		tw::Int ax = tw::grid::naxis(warp->ax);
		for (tw::Int i=lfg[ax];i<=ufg[ax];i++)
			dX(i,ax) += warp->AddedCellWidth(lowSideCells[ax]+i);
	}
}

/// Update the parameters describing global and local hulls (step 4).
/// This is only needed if non-uniform spacings are imposed.
/// This deals only with parameters, not metrics.
/// Requires message passing.
void MetricSpace::UpdateHulls(Task *task)
{
	// Get the local hull size using only local data

	size = 0.0;
	for (tw::Int ax=1;ax<=3;ax++)
		for (tw::Int i=1;i<=dim[ax];i++)
			size[ax] += dX(i,ax);

	// Perform message passing to determine the effect of non-uniform cell widths
	// on the global domain size and the coordinates of the local domain corner.
	// Assumes local domain sizes are already calculated.

	tw::Int axis,src,dst;
	tw::Float inData,outData;
	tw::Float lsize[4] = { size[0] , size[1] , size[2] , size[3] };
	tw::Float gcorn[4] = { globalCorner[0] , globalCorner[1] , globalCorner[2] , globalCorner[3] };
	tw::Float lcorn[4];
	tw::Float gsize[4];

	for (axis=1;axis<=3;axis++)
	{
		task->finiteStrip[axis].Shift(1,1,&src,&dst);
		if (task->domainIndex[axis]==0)
		{
			task->finiteStrip[axis].Send(&lsize[axis],sizeof(tw::Float),dst);
			lcorn[axis] = gcorn[axis];
			outData = lsize[axis]; // in case domains=1
		}
		else
		{
			task->finiteStrip[axis].Recv(&inData,sizeof(tw::Float),src);
			lcorn[axis] = gcorn[axis] + inData;
			outData = inData + lsize[axis];
			if (task->domainIndex[axis]!=task->domains[axis]-1)
				task->finiteStrip[axis].Send(&outData,sizeof(tw::Float),dst);
		}
		task->finiteStrip[axis].Shift(1,-1,&src,&dst);
		if (task->domainIndex[axis]==task->domains[axis]-1)
		{
			gsize[axis] = outData;
			task->finiteStrip[axis].Send(&gsize[axis],sizeof(tw::Float),dst);
		}
		else
		{
			task->finiteStrip[axis].Recv(&inData,sizeof(tw::Float),src);
			gsize[axis] = inData;
			if (task->domainIndex[axis]!=0)
				task->finiteStrip[axis].Send(&gsize[axis],sizeof(tw::Float),dst);
		}
		corner[axis] = lcorn[axis];
		globalSize[axis] = gsize[axis];
	}
}

/// Used by Set*Geometry variants to build positions from spacings
void MetricSpace::SetupPositionArrays()
{
	for (tw::Int ax=1;ax<=3;ax++)
	{
		X(mlb[ax],ax) = corner[ax] - dX(mlb[ax]+1,ax) - 0.5*dX(mlb[ax],ax);
		for (tw::Int i=mlb[ax]+1;i<=mub[ax];i++)
			X(i,ax) = X(i-1,ax) + 0.5*dX(i-1,ax) + 0.5*dX(i,ax);
	}
}

void MetricSpace::SetCartesianGeometry()
{
	SetupPositionArrays();
	gridGeometry = tw::grid::cartesian;
	car = 1.0; cyl = 0.0; sph = 0.0;
	tw::Int i,sx=mnum[1],sz=mnum[3];
	tw::Float *xpos = &gpos[mnum[0]*1];
	//tw::Float *ypos = &gpos[mnum[0]*2];
	tw::Float *zpos = &gpos[mnum[0]*3];
	tw::Float *xwidth = &width[mnum[0]*1];
	tw::Float *ywidth = &width[mnum[0]*2];
	tw::Float *zwidth = &width[mnum[0]*3];
	for (i=0;i<sx;i++)
	{
		cell_area_x[i + 0*sx] = xwidth[i]*ywidth[2];
		cell_area_x[i + 1*sx] = ywidth[2];
		cell_area_x[i + 2*sx] = xwidth[i];
		cell_area_x[i + 3*sx] = xwidth[i]*ywidth[2];
		cell_arc_x[i + 1*sx] = ywidth[2];
		cell_arc_x[i + 2*sx] = 1.0;
	}
	for (i=0;i<sz;i++)
	{
		cell_area_z[i + 0*sz] = zwidth[i];
		cell_area_z[i + 1*sz] = zwidth[i];
		cell_area_z[i + 2*sz] = zwidth[i];
		cell_area_z[i + 3*sz] = 1.0;
		cell_arc_z[i + 0*sz] = 1.0;
		cell_arc_z[i + 1*sz] = 1.0;
	}
	for (i=1;i<sx;i++)
	{
		cell_arc_x[i + 0*sx] = xpos[i] - xpos[i-1];
	}
	for (i=1;i<sz;i++)
	{
		cell_arc_z[i + 2*sz] = zpos[i] - zpos[i-1];
	}
	for (i=1;i<sx;i++)
	{
		cell_area_x[i + 4*sx] = (xpos[i]-xpos[i-1])*ywidth[2];
		cell_area_x[i + 6*sx] = (xpos[i]-xpos[i-1]);
		cell_area_x[i + 7*sx] = (xpos[i]-xpos[i-1])*ywidth[2];
	}
	for (i=1;i<sz;i++)
	{
		cell_area_z[i + 4*sz] = zpos[i]-zpos[i-1];
		cell_area_z[i + 5*sz] = zpos[i]-zpos[i-1];
		cell_area_z[i + 6*sz] = zpos[i]-zpos[i-1];
	}
	for (i=0;i<sx;i++)
	{
		cell_area_x[i + 5*sx] = ywidth[2];
		cell_arc_x[i + 3*sx] = xwidth[i];
		cell_arc_x[i + 4*sx] = ywidth[2];
		cell_arc_x[i + 5*sx] = 1.0;
	}
	for (i=0;i<sz;i++)
	{
		cell_area_z[i + 7*sz] = 1.0;
		cell_arc_z[i + 3*sz] = 1.0;
		cell_arc_z[i + 4*sz] = 1.0;
		cell_arc_z[i + 5*sz] = zwidth[i];
	}
}

void MetricSpace::SetCylindricalGeometry()
{
	// we are using x = r, y = phi, z = z
	SetupPositionArrays();
	gridGeometry = tw::grid::cylindrical;
	car = 0.0; cyl = 1.0; sph = 0.0;
	tw::Int i,sx=mnum[1],sz=mnum[3];
	tw::Float *xpos = &gpos[mnum[0]*1];
	//tw::Float *ypos = &gpos[mnum[0]*2];
	tw::Float *zpos = &gpos[mnum[0]*3];
	tw::Float *xwidth = &width[mnum[0]*1];
	tw::Float *ywidth = &width[mnum[0]*2];
	tw::Float *zwidth = &width[mnum[0]*3];
	for (i=0;i<sx;i++)
	{
		const tw::Float x0 = xpos[i]; const tw::Float x1 = xpos[i] - 0.5*xwidth[i];
		cell_area_x[i + 0*sx] = x0*xwidth[i]*ywidth[2];
		cell_area_x[i + 1*sx] = x1*ywidth[2];
		cell_area_x[i + 2*sx] = xwidth[i];
		cell_area_x[i + 3*sx] = x0*xwidth[i]*ywidth[2];
		cell_arc_x[i + 1*sx] = x0*ywidth[2];
		cell_arc_x[i + 2*sx] = 1.0;
	}
	for (i=0;i<sz;i++)
	{
		cell_area_z[i + 0*sz] = zwidth[i];
		cell_area_z[i + 1*sz] = zwidth[i];
		cell_area_z[i + 2*sz] = zwidth[i];
		cell_area_z[i + 3*sz] = 1.0;
		cell_arc_z[i + 0*sz] = 1.0;
		cell_arc_z[i + 1*sz] = 1.0;
	}
	for (i=1;i<sx;i++)
	{
		cell_arc_x[i + 0*sx] = xpos[i] - xpos[i-1];
	}
	for (i=1;i<sz;i++)
	{
		cell_arc_z[i + 2*sz] = zpos[i] - zpos[i-1];
	}
	for (i=1;i<sx;i++)
	{
		const tw::Float x0 = xpos[i] - 0.5*xwidth[i];
		cell_area_x[i + 4*sx] = x0*(xpos[i]-xpos[i-1])*ywidth[2];
		cell_area_x[i + 6*sx] = (xpos[i]-xpos[i-1]);
		cell_area_x[i + 7*sx] = x0*(xpos[i]-xpos[i-1])*ywidth[2];
	}
	for (i=1;i<sz;i++)
	{
		cell_area_z[i + 4*sz] = zpos[i]-zpos[i-1];
		cell_area_z[i + 5*sz] = zpos[i]-zpos[i-1];
		cell_area_z[i + 6*sz] = zpos[i]-zpos[i-1];
	}
	for (i=0;i<sx;i++)
	{
		const tw::Float x0 = xpos[i] - 0.5*xwidth[i]; const tw::Float x1 = xpos[i];
		cell_area_x[i + 5*sx] = x1*ywidth[2];
		cell_arc_x[i + 3*sx] = xwidth[i];
		cell_arc_x[i + 4*sx] = x0*ywidth[2];
		cell_arc_x[i + 5*sx] = 1.0;
	}
	for (i=0;i<sz;i++)
	{
		cell_area_z[i + 7*sz] = 1.0;
		cell_arc_z[i + 3*sz] = 1.0;
		cell_arc_z[i + 4*sz] = 1.0;
		cell_arc_z[i + 5*sz] = zwidth[i];
	}

	// Don't allow any wall area to strictly collapse
	for (i=0;i<sx;i++)
		for (tw::Int c=0;c<8;c++)
			if (cell_area_x[i + c*sx]==0.0)
				cell_area_x[i + c*sx] = tw::small_pos;
}

void MetricSpace::SetSphericalGeometry()
{
	// we are using q1 = rho, q2 = phi (atan(y/x)), q3 = theta (acos(z/rho))
	SetupPositionArrays();
	gridGeometry = tw::grid::spherical;
	car = 0.0; cyl = 0.0; sph = 1.0;
	tw::Int i,sx=mnum[1],sz=mnum[3];
	tw::Float *xpos = &gpos[mnum[0]*1];
	//tw::Float *ypos = &gpos[mnum[0]*2];
	tw::Float *zpos = &gpos[mnum[0]*3];
	tw::Float *xwidth = &width[mnum[0]*1];
	tw::Float *ywidth = &width[mnum[0]*2];
	tw::Float *zwidth = &width[mnum[0]*3];
	for (i=0;i<sx;i++)
	{
		const tw::Float x0 = xpos[i];
		const tw::Float dx = xwidth[i]; const tw::Float dy = ywidth[2];
		const tw::Float x1 = x0 - 0.5*dx;
		cell_area_x[i + 0*sx] = dx*dy*(dx*dx/6.0 + 2.0*x0*x0);
		cell_area_x[i + 1*sx] = 2.0*x1*x1*dy;
		cell_area_x[i + 2*sx] = dx*x0;
		cell_area_x[i + 3*sx] = dx*dy*x0;
		cell_arc_x[i + 1*sx] = x0*ywidth[2];
		cell_arc_x[i + 2*sx] = x0;
	}
	for (i=0;i<sz;i++)
	{
		const tw::Float z0 = zpos[i];
		const tw::Float dz = zwidth[i];
		const tw::Float z1 = z0 - 0.5*dz;
		cell_area_z[i + 0*sz] = sin(0.5*dz)*sin(z0);
		cell_area_z[i + 1*sz] = sin(0.5*dz)*sin(z0);
		cell_area_z[i + 2*sz] = dz;
		cell_area_z[i + 3*sz] = sin(z1);
		cell_arc_z[i + 0*sz] = 1.0;
		cell_arc_z[i + 1*sz] = sin(z0);
	}
	for (i=1;i<sx;i++)
	{
		cell_arc_x[i + 0*sx] = xpos[i] - xpos[i-1];
	}
	for (i=1;i<sz;i++)
	{
		cell_arc_z[i + 2*sz] = zpos[i] - zpos[i-1];
	}
	for (i=1;i<sx;i++)
	{
		const tw::Float x0 = xpos[i] - 0.5*xwidth[i];
		const tw::Float dx = xpos[i] - xpos[i-1]; const tw::Float dy = ywidth[2];
		cell_area_x[i + 4*sx] = dx*dy*(dx*dx/6.0 + 2.0*x0*x0);
		cell_area_x[i + 6*sx] = dx*x0;
		cell_area_x[i + 7*sx] = dx*dy*x0;
	}
	for (i=1;i<sz;i++)
	{
		const tw::Float z0 = zpos[i] - 0.5*zwidth[i];
		const tw::Float dz = zpos[i] - zpos[i-1];
		cell_area_z[i + 4*sz] = sin(0.5*dz)*sin(z0);
		cell_area_z[i + 5*sz] = sin(0.5*dz)*sin(z0);
		cell_area_z[i + 6*sz] = dz;
	}
	for (i=0;i<sx;i++)
	{
		const tw::Float x0 = xpos[i] - 0.5*xwidth[i];
 		const tw::Float dy = ywidth[2];
		const tw::Float x1 = xpos[i];
		cell_area_x[i + 5*sx] = 2.0*x1*x1*dy;
		cell_arc_x[i + 3*sx] = xwidth[i];
		cell_arc_x[i + 4*sx] = x0*dy;
		cell_arc_x[i + 5*sx] = x0;
	}
	for (i=0;i<sz;i++)
	{
		const tw::Float z0 = zpos[i] - 0.5*zwidth[i];
		const tw::Float z1 = zpos[i];
		cell_area_z[i + 7*sz] = sin(z1);
		cell_arc_z[i + 3*sz] = 1.0;
		cell_arc_z[i + 4*sz] = sin(z0);
		cell_arc_z[i + 5*sz] = zwidth[i];
	}

	// Don't allow any wall area to strictly collapse
	// N.b. however tw::small_pos^2 = underflow
	for (i=0;i<sx;i++)
		for (tw::Int c=0;c<8;c++)
			if (cell_area_x[i + c*sx]==0.0)
				cell_area_x[i + c*sx] = tw::small_pos;
	for (i=0;i<sz;i++)
		for (tw::Int c=0;c<8;c++)
			if (cell_area_z[i + c*sz]==0.0)
				cell_area_z[i + c*sz] = tw::small_pos;
}

tw::Float MetricSpace::ToLab(const Evolution& evo,tw::Float zeta,tw::Float relativeTime)
{
	return zeta + evo.antiWindowPosition - evo.windowPosition + evo.signalSpeed*(evo.elapsedTime + relativeTime - 0.5*spacing[0]);
}

tw::Float MetricSpace::ToLight(const Evolution& evo,tw::Float z,tw::Float relativeTime)
{
	return z + evo.windowPosition - evo.antiWindowPosition - evo.signalSpeed*(evo.elapsedTime + relativeTime - 0.5*spacing[0]);
}

template <class T,class U>
U MetricSpace::ValueOnLabGrid(const Evolution& evo,T& A,tw::strip s,tw::Int k,tw::Float relativeTime)
{
	// Take a quantity known on the light grid and get its value in a cell of the lab grid
	tw::Int klight;
	tw::Float z,zeta,w;
	z = Pos(s,k).z - corner[3];
	zeta = ToLight(evo,z,relativeTime);
	klight = MyFloor(zeta*freq[3] + 0.5001);
	w = 0.5 - zeta*freq[3] + tw::Float(klight);
	return w*A(s,klight) + (one - w)*A(s,klight+1);
}

template <class T,class U>
U MetricSpace::ValueOnLightGrid(const Evolution& evo,T& A,tw::strip s,tw::Int k,tw::Float relativeTime)
{
	// Take a quantity known on the lab grid and get its value in a cell of the light grid
	tw::Int klab;
	tw::Float z,zeta,w;
	zeta = Pos(s,k).z - corner[3];
	z = ToLab(evo,zeta,relativeTime);
	klab = MyFloor(z*freq[3] + 0.4999);
	w = 0.5 - z*freq[3] + tw::Float(klab);
	return w*A(s,klab) + (one - w)*A(s,klab+1);
}

void MetricSpace::ReadCheckpoint(std::ifstream& inFile)
{
	DiscreteSpace::ReadCheckpoint(inFile);
	inFile.read((char *)mnum,sizeof(mnum));
	inFile.read((char *)mlb,sizeof(mlb));
	inFile.read((char *)mub,sizeof(mub));
	inFile.read((char *)&car,sizeof(tw::Float));
	inFile.read((char *)&cyl,sizeof(tw::Float));
	inFile.read((char *)&sph,sizeof(tw::Float));
	Allocate();
	inFile.read((char *)&gpos[0],sizeof(tw::Float)*gpos.size());
	inFile.read((char *)&width[0],sizeof(tw::Float)*width.size());
	inFile.read((char *)&cell_area_x[0],sizeof(tw::Float)*cell_area_x.size());
	inFile.read((char *)&cell_arc_x[0],sizeof(tw::Float)*cell_arc_x.size());
	inFile.read((char *)&cell_area_z[0],sizeof(tw::Float)*cell_area_z.size());
	inFile.read((char *)&cell_arc_z[0],sizeof(tw::Float)*cell_arc_z.size());
}

void MetricSpace::WriteCheckpoint(std::ofstream& outFile)
{
	DiscreteSpace::WriteCheckpoint(outFile);
	outFile.write((char *)mnum,sizeof(mnum));
	outFile.write((char *)mlb,sizeof(mlb));
	outFile.write((char *)mub,sizeof(mub));
	outFile.write((char *)&car,sizeof(tw::Float));
	outFile.write((char *)&cyl,sizeof(tw::Float));
	outFile.write((char *)&sph,sizeof(tw::Float));
	outFile.write((char *)&gpos[0],sizeof(tw::Float)*gpos.size());
	outFile.write((char *)&width[0],sizeof(tw::Float)*width.size());
	outFile.write((char *)&cell_area_x[0],sizeof(tw::Float)*cell_area_x.size());
	outFile.write((char *)&cell_arc_x[0],sizeof(tw::Float)*cell_arc_x.size());
	outFile.write((char *)&cell_area_z[0],sizeof(tw::Float)*cell_area_z.size());
	outFile.write((char *)&cell_arc_z[0],sizeof(tw::Float)*cell_arc_z.size());
}
