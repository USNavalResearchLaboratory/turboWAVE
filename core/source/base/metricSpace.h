// This abstract class is needed to anticipate ComputeTool::Warp.
struct warp_base
{
	tw::grid::axis ax;
	virtual tw::Float AddedCellWidth(tw::Int globalCell) = 0;
};

struct MetricSpace:DiscreteSpace
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
	void SetTopology(Task& task,const tw::vec3& gcorner,const tw::vec3& gsize,tw::Int ghostCellLayers);
	void Allocate();
	void SetSpacings(Task& task);
	void UpdateHulls(Task& task);
	void SetupPositionArrays();
	void SetCartesianGeometry();
	void SetCylindricalGeometry();
	void SetSphericalGeometry();
public:
	void Resize(Task& task,
		const tw::vec3& gcorner,
		const tw::vec3& gsize,
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
