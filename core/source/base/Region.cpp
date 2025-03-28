module;

#include <tree_sitter/api.h>
#include "tw_includes.h"

export module region;
import input;
import metric_space;
import tensor;

enum regionSpec { baseRegion , entireRegion , rectRegion , prismRegion , circRegion , cylinderRegion , cylindricalShellRegion ,
	roundedCylinderRegion , ellipsoidRegion , trueSphereRegion , boxArrayRegion , torusRegion , coneRegion , tangentOgiveRegion };

export struct Region
{
	regionSpec rgnType;
	std::string name;
	tw::vec3 center;
	tw::vec3 rbox;
	tw::basis orientation;
	tw::Int rawBounds[6];
	tw::Int localBounds[6];
	tw::Int globalBounds[6];
	bool intersectsDomain;
	bool complement,intersection,moveWithWindow;
	std::vector<Region*> composite;
	std::vector<Region*>& masterList;

	// input processing aids
	tw::input::DirectiveReader directives;
	tw::vec3 temp_vec3;
	tw::Float temp_Float;
	tw::Float temp_bounds[6];

	Region(std::vector<Region*>& ml);
	virtual ~Region() {}
	virtual bool Inside(const tw::vec3& pos,const MetricSpace& ds) const
	{
		bool ans = intersection;
		tw::vec3 p = pos - center;
		orientation.ExpressInBasis(&p);
		for (tw::Int i=0;i<composite.size();i++)
			if (intersection)
				ans = ans && composite[i]->Inside(p,ds);
			else
				ans = ans || composite[i]->Inside(p,ds);
		return complement ^ ans;
	}
	virtual void Initialize(const MetricSpace& ds,Task *tsk);
	virtual void Translate(const tw::vec3& dr)
	{
		center += dr;
	}
	void GetBoxLim(tw::Float *low,tw::Float *high,tw::Int ax)
	{
		*low = center[ax-1] - rbox[ax-1];
		*high = center[ax-1] + rbox[ax-1];
	}
	void GetGlobalCellBounds(tw::Int g[6]) const
	{
		for (tw::Int i=0;i<6;i++)
			g[i] = globalBounds[i];
	}
	void GetLocalCellBounds(tw::Int *x0,tw::Int *x1,tw::Int *y0,tw::Int *y1,tw::Int *z0,tw::Int *z1) const
	{
		*x0 = localBounds[0]; *x1 = localBounds[1];
		*y0 = localBounds[2]; *y1 = localBounds[3];
		*z0 = localBounds[4]; *z1 = localBounds[5];
	}
	void GetRawCellBounds(tw::Int *x0,tw::Int *x1,tw::Int *y0,tw::Int *y1,tw::Int *z0,tw::Int *z1) const
	{
		*x0 = rawBounds[0]; *x1 = rawBounds[1];
		*y0 = rawBounds[2]; *y1 = rawBounds[3];
		*z0 = rawBounds[4]; *z1 = rawBounds[5];
	}

	static Region* CreateObjectFromString(std::vector<Region*>& ml,const std::string& str);
	static Region* FindRegion(std::vector<Region*>& ml,const std::string& name);

	virtual void ReadInputFileBlock(TSTreeCursor *curs,const std::string& src);
	virtual void ReadInputFileDirective(const TSTreeCursor *curs,const std::string& src);
	virtual void ReadCheckpoint(std::ifstream& inFile);
	virtual void WriteCheckpoint(std::ofstream& outFile);
};

export struct EntireRegion:Region
{
	EntireRegion(std::vector<Region*>& ml) : Region(ml)
	{
		rgnType = entireRegion;
	}
	virtual bool Inside(const tw::vec3& pos,const MetricSpace& ds) const
	{
		return complement ^ true;
	}
};

export struct RectRegion:Region
{
	RectRegion(std::vector<Region*>& ml) : Region(ml)
	{
		rgnType = rectRegion;
	}
	virtual bool Inside(const tw::vec3& pos,const MetricSpace& ds) const
	{
		tw::vec3 p = pos - center;
		orientation.ExpressInBasis(&p);
		return complement ^ (fabs(p.x)<rbox.x && fabs(p.y)<rbox.y && fabs(p.z)<rbox.z);
	}
};

export struct PrismRegion:Region
{
	PrismRegion(std::vector<Region*>& ml) : Region(ml)
	{
		rgnType = prismRegion;
	}
	virtual bool Inside(const tw::vec3& pos,const MetricSpace& ds) const
	{
		tw::vec3 p = pos - center;
		orientation.ExpressInBasis(&p);
		return complement ^ (fabs(p.x)<rbox.x && fabs(p.y)<rbox.y && fabs(p.z)<rbox.z*0.5*(rbox.x-p.x)/rbox.x);
	}
};

export struct CircRegion:Region
{
	CircRegion(std::vector<Region*>& ml) : Region(ml)
	{
		rgnType = circRegion;
	}
	virtual bool Inside(const tw::vec3& pos,const MetricSpace& ds) const
	{
		return complement ^ (Norm(pos-center)<sqr(rbox.x));
	}
};

export struct CylinderRegion:Region
{
	CylinderRegion(std::vector<Region*>& ml) : Region(ml)
	{
		rgnType = cylinderRegion;
	}
	virtual bool Inside(const tw::vec3& pos,const MetricSpace& ds) const
	{
		tw::vec3 p = pos - center;
		orientation.ExpressInBasis(&p);
		return complement ^ (sqr(p.x) + sqr(p.y) < sqr(rbox.x) && fabs(p.z) < rbox.z);
	}
};

export struct CylindricalShellRegion:Region
{
	tw::Float innerRadius,outerRadius;
	CylindricalShellRegion(std::vector<Region*>& ml);
	virtual bool Inside(const tw::vec3& pos,const MetricSpace& ds) const
	{
		tw::Float rho;
		tw::vec3 p = pos - center;
		orientation.ExpressInBasis(&p);
		rho = sqrt(p.x*p.x + p.y*p.y);
		return complement ^ ( rho < outerRadius && rho > innerRadius && sqr(p.z) < sqr(rbox.z));
	}
	virtual void ReadInputFileDirective(const TSTreeCursor *curs,const std::string& src);
	virtual void ReadCheckpoint(std::ifstream& inFile);
	virtual void WriteCheckpoint(std::ofstream& outFile);
};

export struct RoundedCylinderRegion:Region
{
	RoundedCylinderRegion(std::vector<Region*>& ml) : Region(ml)
	{
		rgnType = roundedCylinderRegion;
	}
	virtual bool Inside(const tw::vec3& pos,const MetricSpace& ds) const
	{
		bool ans;
		tw::vec3 p = pos - center;
		orientation.ExpressInBasis(&p);
		ans = sqr(p.x) + sqr(p.y) < sqr(rbox.x) && fabs(p.z) < rbox.z;
		ans = ans || Norm(p - tw::vec3(0,0,rbox.z)) < sqr(rbox.x);
		ans = ans || Norm(p + tw::vec3(0,0,rbox.z)) < sqr(rbox.x);
		return complement ^ ans;
	}
};

export struct EllipsoidRegion:Region
{
	EllipsoidRegion(std::vector<Region*>& ml) : Region(ml)
	{
		rgnType = ellipsoidRegion;
	}
	virtual bool Inside(const tw::vec3& pos,const MetricSpace& ds) const
	{
		tw::vec3 p = pos - center;
		orientation.ExpressInBasis(&p);
		return complement ^ (sqr(p.x/rbox.x) + sqr(p.y/rbox.y) + sqr(p.z/rbox.z) < 1.0);
	}
};

export struct TrueSphere:CircRegion
{
	TrueSphere(std::vector<Region*>& ml) : CircRegion(ml)
	{
		rgnType = trueSphereRegion;
	}
	virtual bool Inside(const tw::vec3& pos,const MetricSpace& ds) const;
};

export struct BoxArrayRegion:Region
{
	tw::vec3 size,spacing;

	BoxArrayRegion(std::vector<Region*>& ml);
	virtual bool Inside(const tw::vec3& pos,const MetricSpace& ds) const
	{
		tw::vec3 p = pos - center;
		orientation.ExpressInBasis(&p);
		p += 0.5*size;
		tw::Int i = MyFloor(p.x/spacing.x);
		tw::Int j = MyFloor(p.y/spacing.y);
		tw::Int k = MyFloor(p.z/spacing.z);
		bool ans = p.x - tw::Float(i)*spacing.x < size.x && p.y - tw::Float(j)*spacing.y < size.y && p.z - tw::Float(k)*spacing.z < size.z;
		return complement ^ ans;
	}
	virtual void ReadInputFileDirective(const TSTreeCursor *curs,const std::string& src);
	virtual void ReadCheckpoint(std::ifstream& inFile);
	virtual void WriteCheckpoint(std::ofstream& outFile);
};

export struct TorusRegion:Region
{
	tw::Float majorRadius,minorRadius;
	TorusRegion(std::vector<Region*>& ml);
	virtual bool Inside(const tw::vec3& pos,const MetricSpace& ds) const
	{
		tw::Float rho;
		tw::vec3 p = pos - center;
		orientation.ExpressInBasis(&p);
		rho = sqrt(p.x*p.x + p.y*p.y);
		return complement ^ ( rho > majorRadius-minorRadius && rho < majorRadius+minorRadius &&
			sqr(p.z) < sqr(minorRadius)-sqr(rho-majorRadius));
	}
	virtual void ReadInputFileDirective(const TSTreeCursor *curs,const std::string& src);
	virtual void ReadCheckpoint(std::ifstream& inFile);
	virtual void WriteCheckpoint(std::ofstream& outFile);
};

export struct ConeRegion:Region
{
	tw::Float majorRadius,minorRadius;
	ConeRegion(std::vector<Region*>& ml);
	virtual bool Inside(const tw::vec3& pos,const MetricSpace& ds) const
	{
		tw::Float rho,rOfz;
		tw::vec3 p = pos - center;
		orientation.ExpressInBasis(&p);
		rho = sqrt(p.x*p.x + p.y*p.y);
		rOfz = minorRadius - (p.z - rbox.z)*(majorRadius-minorRadius)/(2.0*rbox.z);
		return complement ^ ( rho < rOfz && sqr(p.z) < sqr(rbox.z));
	}
	virtual void ReadInputFileDirective(const TSTreeCursor *curs,const std::string& src);
	virtual void ReadCheckpoint(std::ifstream& inFile);
	virtual void WriteCheckpoint(std::ofstream& outFile);
};

export struct TangentOgiveRegion:Region
{
	tw::Float tipRadius,bodyRadius;
	TangentOgiveRegion(std::vector<Region*>& ml);
	virtual bool Inside(const tw::vec3& pos,const MetricSpace& ds) const
	{
		tw::Float rho,rOfz,ogiveRadius,x0,xt,yt;
		tw::vec3 p = pos - center;
		orientation.ExpressInBasis(&p);

		x0 = 2.0*rbox.z - tipRadius;
		ogiveRadius = (x0*x0+sqr(bodyRadius)-sqr(tipRadius))/(2.0*(bodyRadius-tipRadius));
		yt = tipRadius*(ogiveRadius-bodyRadius)/(ogiveRadius-tipRadius);
		xt = x0 + sqrt(sqr(tipRadius)-yt*yt);

		rho = sqrt(p.x*p.x + p.y*p.y);
		p.z += rbox.z;
		rOfz = -1.0;
		if (p.z>0.0 && p.z<xt)
			rOfz = sqrt(sqr(ogiveRadius)-p.z*p.z) + bodyRadius - ogiveRadius;
		if (p.z>=xt && p.z<2.0*rbox.z)
			rOfz = sqrt(sqr(tipRadius)-sqr(p.z-x0));
		return complement ^ (rho < rOfz);
	}
	virtual void ReadInputFileDirective(const TSTreeCursor *curs,const std::string& src);
	virtual void ReadCheckpoint(std::ifstream& inFile);
	virtual void WriteCheckpoint(std::ofstream& outFile);
};

Region::Region(std::vector<Region*>& ml) : masterList(ml)
{
	name = "entire";
	rgnType = baseRegion;
	center = tw::vec3(0.0);
	rbox = tw::vec3(tw::big_pos);
	intersectsDomain = true;
	complement = false;
	intersection = false;
	moveWithWindow = true;
	orientation.u = tw::vec3(1,0,0);
	orientation.v = tw::vec3(0,1,0);
	orientation.w = tw::vec3(0,0,1);
	// Regions use directives somewhat differently.
	// Dependent geometric variables are updated as the input file is read in.
	// This happens in ReadInputFileDirective (normally used only for custom inputs).
	// It might be better to make this more conventional by updating dependent variables in the Initialize function instead.
	directives.Add("bounds",new tw::input::Numbers<tw::Float>(&temp_bounds[0],6),false);
	directives.Add("radius",new tw::input::Float(&temp_Float),false);
	directives.Add("length",new tw::input::Float(&temp_Float),false);
	directives.Add("elements",new tw::input::Custom,false);
	directives.Add("translation",new tw::input::Vec3(&temp_vec3),false);
	directives.Add("rotation about x",new tw::input::Float(&temp_Float),false);
	directives.Add("rotation about y",new tw::input::Float(&temp_Float),false);
	directives.Add("rotation about z",new tw::input::Float(&temp_Float),false);
	directives.Add("move with window",new tw::input::Bool(&moveWithWindow),false);
	directives.Add("complement",new tw::input::Bool(&complement),false);
}

void Region::Initialize(const MetricSpace& ds,Task *task)
{
	// Initialization is needed to set up index limits
	// The index limits define a bounding rectangle drawn around the region

	tw::Int dims[3];
	tw::Float lims[6];
	tw::Int i,d,low,high;
	std::vector<tw::Float> temp;

	// Set up the index limits such that the region is in cells [xl,xh] * [yl,yh] * [zl,zh] (inclusive)
	// std::lower_bound returns first element >= to the test data
	// if r0.x < xpos[0], xl = 0   ,	if r1.x < xpos[0], xh = -1
	// if r0.x > xpos[xN1], xl = xN1+1   ,   if r1.x > xpos[xN1], xh = xN1

	for (d=0;d<3;d++)
	{
		dims[d] = ds.Dim(d+1);
		lims[2*d] = center[d] - rbox[d];
		lims[2*d+1] = center[d] + rbox[d];
	}

	for (d=0;d<3;d++)
	{
		temp.resize(dims[d]+2);
		for (i=0;i<=dims[d]+1;i++)
			temp[i] = ds.X(i,d+1);
		rawBounds[2*d] = tw::Int(std::lower_bound(temp.begin(),temp.end(),lims[2*d]) - temp.begin());
		rawBounds[2*d+1] = tw::Int(std::lower_bound(temp.begin(),temp.end(),lims[2*d+1]) - temp.begin())-1;
	}

	// Bounds checking
	// localBounds will imply at least one null loop (low>high) in case entire region is outside domain
	// low values fall in the range [1,dim+2] and high values fall in the range [-1,dim]

	for (d=0;d<3;d++)
	{
		localBounds[2*d] = rawBounds[2*d]<1 ? 1 : rawBounds[2*d];
		localBounds[2*d+1] = rawBounds[2*d+1]>dims[d] ? dims[d] : rawBounds[2*d+1];

		intersectsDomain = intersectsDomain && (localBounds[2*d]<=localBounds[2*d+1]);
	}

	// Global bounds

	for (d=0;d<3;d++)
	{
		if (localBounds[2*d]>=1 && localBounds[2*d]<=dims[d])
			low = ds.GlobalCellIndex(localBounds[2*d],d+1);
		else
			low = ds.GlobalDim(d+1);
		globalBounds[2*d] = task->strip[d+1].GetMin(low);
	}

	for (d=0;d<3;d++)
	{
		if (localBounds[2*d+1]>=1 && localBounds[2*d+1]<=dims[d])
			high = ds.GlobalCellIndex(localBounds[2*d+1],d+1);
		else
			high = 1;
		globalBounds[2*d+1] = task->strip[d+1].GetMax(high);
	}
}

Region* Region::CreateObjectFromString(std::vector<Region*>& ml,const std::string& word)
{
	if (word=="entire")
		return new EntireRegion(ml);
	if (word=="rect")
		return new RectRegion(ml);
	if (word=="prism")
		return new PrismRegion(ml);
	if (word=="ellipsoid")
		return new EllipsoidRegion(ml);
	if (word=="circ")
		return new CircRegion(ml);
	if (word=="cylinder")
		return new CylinderRegion(ml);
	if (word=="rounded_cylinder")
		return  new RoundedCylinderRegion(ml);
	if (word=="true_sphere")
		return new TrueSphere(ml);
	if (word=="box_array")
		return new BoxArrayRegion(ml);
	if (word=="torus")
		return new TorusRegion(ml);
	if (word=="cone")
		return new ConeRegion(ml);
	if (word=="tangent_ogive")
		return new TangentOgiveRegion(ml);
	if (word=="cylindrical_shell")
		return new CylindricalShellRegion(ml);

	if (word=="union")
	{
		Region *ans = new Region(ml);
		ans->intersection = false;
		return ans;
	}
	if (word=="intersection")
	{
		Region *ans = new Region(ml);
		ans->intersection = true;
		return ans;
	}

	return new EntireRegion(ml);
}

Region* Region::FindRegion(std::vector<Region*>& ml,const std::string& name)
{
	tw::Int i;

	for (i=0;i<ml.size();i++)
		if (ml[i]->name==name)
			return ml[i];
	throw tw::FatalError("Could not find region <"+name+">.");
	return NULL;
}

/// @brief read all directives in the block
/// @param curs can be on block or on first child of block
/// @param src source document
void Region::ReadInputFileBlock(TSTreeCursor *curs,const std::string& src)
{
	if (tw::input::node_kind(curs)=="block") {
		ts_tree_cursor_goto_first_child(curs);
	}
	do
	{
		directives.ReadNext(curs,src);
		// following is unconditional, see comments up top
		TSTreeCursor curs1 = ts_tree_cursor_copy(curs);
		ReadInputFileDirective(&curs1,src);
	} while (ts_tree_cursor_goto_next_sibling(curs));
	directives.ThrowErrorIfMissingKeys(name);
}

void Region::ReadInputFileDirective(const TSTreeCursor *curs0,const std::string& src)
{
	TSTreeCursor curs = ts_tree_cursor_copy(curs0);
	ts_tree_cursor_goto_first_child(&curs);
	std::string com = tw::input::node_text(&curs,src);

	if (com=="bounds")
	{
		tw::vec3 r0,r1;
		r0.x = temp_bounds[0]; r1.x = temp_bounds[1];
		r0.y = temp_bounds[2]; r1.y = temp_bounds[3];
		r0.z = temp_bounds[4]; r1.z = temp_bounds[5];
		if (r0.x>r1.x) std::swap(r0.x,r1.x);
		if (r0.y>r1.y) std::swap(r0.y,r1.y);
		if (r0.z>r1.z) std::swap(r0.z,r1.z);
		center = 0.5*(r0 + r1);
		rbox = 0.5*(r1 - r0);
	}
	if (com=="radius")
	{
		rbox.x = temp_Float;
		if (rbox.y > 0.99*tw::big_pos)
			rbox.y = rbox.x;
		if (rbox.z > 0.99*tw::big_pos)
			rbox.z = rbox.x;
	}
	if (com=="length")
	{
		rbox.z = 0.5*temp_Float;
	}
	if (com=="elements")
	{
		Region *curr;
		tw::input::next_named_node(&curs,false);
		if (tw::input::node_kind(&curs) != "list") {
			throw tw::FatalError("expected `elements` to be a list while processing region");
		}
		ts_tree_cursor_goto_first_child(&curs);
		while (tw::input::next_named_node(&curs,false)) {
			std::string elName = tw::input::node_text(&curs,src);
			tw::input::StripQuotes(elName);
			curr = FindRegion(masterList,elName);
			composite.push_back(curr);
		}
	}
	if (com=="translation") // eg, translation = 1 0 0
	{
		center += temp_vec3;
	}
	if (com=="rotation about x")
	{
		orientation.u.RotateX(temp_Float);
		orientation.v.RotateX(temp_Float);
		orientation.w.RotateX(temp_Float);
		center.RotateX(temp_Float);
	}
	if (com=="rotation about y")
	{
		orientation.u.RotateY(temp_Float);
		orientation.v.RotateY(temp_Float);
		orientation.w.RotateY(temp_Float);
		center.RotateY(temp_Float);
	}
	if (com=="rotation about z")
	{
		orientation.u.RotateZ(temp_Float);
		orientation.v.RotateZ(temp_Float);
		orientation.w.RotateZ(temp_Float);
		center.RotateZ(temp_Float);
	}
}

void Region::ReadCheckpoint(std::ifstream& inFile)
{
	inFile.read((char *)&center,sizeof(center));
	inFile.read((char *)&rbox,sizeof(rbox));
	inFile.read((char *)&orientation,sizeof(orientation));
	inFile.read((char *)rawBounds,sizeof(rawBounds));
	inFile.read((char *)localBounds,sizeof(localBounds));
	inFile.read((char *)globalBounds,sizeof(globalBounds));
	inFile.read((char *)&intersectsDomain,sizeof(bool));
	inFile.read((char *)&complement,sizeof(bool));
	inFile.read((char *)&intersection,sizeof(bool));
	inFile.read((char *)&moveWithWindow,sizeof(bool));
}

void Region::WriteCheckpoint(std::ofstream& outFile)
{
	outFile << name << " ";
	outFile.write((char *)&center,sizeof(center));
	outFile.write((char *)&rbox,sizeof(rbox));
	outFile.write((char *)&orientation,sizeof(orientation));
	outFile.write((char *)rawBounds,sizeof(rawBounds));
	outFile.write((char *)localBounds,sizeof(localBounds));
	outFile.write((char *)globalBounds,sizeof(globalBounds));
	outFile.write((char *)&intersectsDomain,sizeof(bool));
	outFile.write((char *)&complement,sizeof(bool));
	outFile.write((char *)&intersection,sizeof(bool));
	outFile.write((char *)&moveWithWindow,sizeof(bool));
}

bool TrueSphere::Inside(const tw::vec3& pos,const MetricSpace& ds) const
{
	tw::vec3 c_cart = center;
	tw::vec3 p_cart = pos;
	ds.CurvilinearToCartesian(&c_cart);
	ds.CurvilinearToCartesian(&p_cart);
	return complement ^ (Norm(p_cart-c_cart)<sqr(rbox.x));
}

BoxArrayRegion::BoxArrayRegion(std::vector<Region*>& ml) : Region(ml)
{
	rgnType = boxArrayRegion;
	size = tw::vec3(1,1,1);
	spacing = tw::vec3(2,2,2);
	directives.Add("size",new tw::input::Vec3(&size));
	directives.Add("spacing",new tw::input::Vec3(&spacing));
}

void BoxArrayRegion::ReadInputFileDirective(const TSTreeCursor *curs,const std::string& src)
{
	Region::ReadInputFileDirective(curs,src);
	rbox = 0.5*size;
}

void BoxArrayRegion::ReadCheckpoint(std::ifstream& inFile)
{
	Region::ReadCheckpoint(inFile);
	inFile.read((char *)&size,sizeof(tw::vec3));
	inFile.read((char *)&spacing,sizeof(tw::vec3));
}

void BoxArrayRegion::WriteCheckpoint(std::ofstream& outFile)
{
	Region::WriteCheckpoint(outFile);
	outFile.write((char *)&size,sizeof(tw::vec3));
	outFile.write((char *)&spacing,sizeof(tw::vec3));
}

TorusRegion::TorusRegion(std::vector<Region*>& ml) : Region(ml)
{
	rgnType = torusRegion;
	majorRadius = 1.0;
	minorRadius = 0.1;
	rbox = tw::vec3(majorRadius+minorRadius,majorRadius+minorRadius,minorRadius);
	directives.Add("minor radius",new tw::input::Float(&minorRadius));
	directives.Add("major radius",new tw::input::Float(&majorRadius));
}

void TorusRegion::ReadInputFileDirective(const TSTreeCursor *curs,const std::string& src)
{
	Region::ReadInputFileDirective(curs,src);
	rbox = tw::vec3(majorRadius+minorRadius,majorRadius+minorRadius,minorRadius);
}

void TorusRegion::ReadCheckpoint(std::ifstream& inFile)
{
	Region::ReadCheckpoint(inFile);
	inFile.read((char *)&minorRadius,sizeof(tw::Float));
	inFile.read((char *)&majorRadius,sizeof(tw::Float));
}

void TorusRegion::WriteCheckpoint(std::ofstream& outFile)
{
	Region::WriteCheckpoint(outFile);
	outFile.write((char *)&minorRadius,sizeof(tw::Float));
	outFile.write((char *)&majorRadius,sizeof(tw::Float));
}

ConeRegion::ConeRegion(std::vector<Region*>& ml) : Region(ml)
{
	rgnType = coneRegion;
	majorRadius = 1.0;
	minorRadius = 0.1;
	rbox = tw::vec3(majorRadius,majorRadius,tw::big_pos);
	directives.Add("tip radius",new tw::input::Float(&minorRadius));
	directives.Add("base radius",new tw::input::Float(&majorRadius));
}

void ConeRegion::ReadInputFileDirective(const TSTreeCursor *curs,const std::string& src)
{
	Region::ReadInputFileDirective(curs,src);
	rbox.x = rbox.y = majorRadius;
}

void ConeRegion::ReadCheckpoint(std::ifstream& inFile)
{
	Region::ReadCheckpoint(inFile);
	inFile.read((char *)&minorRadius,sizeof(tw::Float));
	inFile.read((char *)&majorRadius,sizeof(tw::Float));
}

void ConeRegion::WriteCheckpoint(std::ofstream& outFile)
{
	Region::WriteCheckpoint(outFile);
	outFile.write((char *)&minorRadius,sizeof(tw::Float));
	outFile.write((char *)&majorRadius,sizeof(tw::Float));
}

TangentOgiveRegion::TangentOgiveRegion(std::vector<Region*>& ml) : Region(ml)
{
	rgnType = tangentOgiveRegion;
	tipRadius = 1.0;
	bodyRadius = 5.0;
	rbox = tw::vec3(bodyRadius,bodyRadius,tw::big_pos);
	directives.Add("tip radius",new tw::input::Float(&tipRadius));
	directives.Add("body radius",new tw::input::Float(&bodyRadius));
}

void TangentOgiveRegion::ReadInputFileDirective(const TSTreeCursor *curs,const std::string& src)
{
	Region::ReadInputFileDirective(curs,src);
	rbox.x = rbox.y = bodyRadius;
}

void TangentOgiveRegion::ReadCheckpoint(std::ifstream& inFile)
{
	Region::ReadCheckpoint(inFile);
	inFile.read((char *)&tipRadius,sizeof(tw::Float));
	inFile.read((char *)&bodyRadius,sizeof(tw::Float));
}

void TangentOgiveRegion::WriteCheckpoint(std::ofstream& outFile)
{
	Region::WriteCheckpoint(outFile);
	outFile.write((char *)&tipRadius,sizeof(tw::Float));
	outFile.write((char *)&bodyRadius,sizeof(tw::Float));
}

CylindricalShellRegion::CylindricalShellRegion(std::vector<Region*>& ml) : Region(ml)
{
	rgnType = cylindricalShellRegion;
	innerRadius = 1.0;
	outerRadius = 2.0;
	rbox = tw::vec3(outerRadius,outerRadius,tw::big_pos);
	directives.Add("inner radius",new tw::input::Float(&innerRadius));
	directives.Add("outer radius",new tw::input::Float(&outerRadius));
}

void CylindricalShellRegion::ReadInputFileDirective(const TSTreeCursor *curs,const std::string& src)
{
	Region::ReadInputFileDirective(curs,src);
	rbox.x = rbox.y = outerRadius;
}

void CylindricalShellRegion::ReadCheckpoint(std::ifstream& inFile)
{
	Region::ReadCheckpoint(inFile);
	inFile.read((char *)&innerRadius,sizeof(tw::Float));
	inFile.read((char *)&outerRadius,sizeof(tw::Float));
}

void CylindricalShellRegion::WriteCheckpoint(std::ofstream& outFile)
{
	Region::WriteCheckpoint(outFile);
	outFile.write((char *)&innerRadius,sizeof(tw::Float));
	outFile.write((char *)&outerRadius,sizeof(tw::Float));
}
