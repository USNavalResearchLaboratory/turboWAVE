enum regionSpec { baseRegion , entireRegion , rectRegion , prismRegion , circRegion , cylinderRegion , cylindricalShellRegion ,
	roundedCylinderRegion , ellipsoidRegion , trueSphereRegion , boxArrayRegion , torusRegion , coneRegion , tangentOgiveRegion };

struct Region
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

	Region(std::vector<Region*>& ml) : masterList(ml)
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
	}
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
	void GetGlobalCellBounds(tw::Int *x0,tw::Int *x1,tw::Int *y0,tw::Int *y1,tw::Int *z0,tw::Int *z1) const
	{
		*x0 = globalBounds[0]; *x1 = globalBounds[1];
		*y0 = globalBounds[2]; *y1 = globalBounds[3];
		*z0 = globalBounds[4]; *z1 = globalBounds[5];
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
	static Region* CreateObjectFromFile(std::vector<Region*>& ml,std::ifstream& inFile);
	static Region* ReadRegion(std::vector<Region*>& ml,Region *curr,std::stringstream& source,const std::string& word);
	static Region* FindRegion(std::vector<Region*>& ml,const std::string& name);

	virtual void ReadInputFileBlock(std::stringstream& inputString);
	virtual void ReadInputFileDirective(std::stringstream& inputString,const std::string& command);
	virtual void ReadData(std::ifstream& inFile);
	virtual void WriteData(std::ofstream& outFile);
};

struct EntireRegion:Region
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

struct RectRegion:Region
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

struct PrismRegion:Region
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

struct CircRegion:Region
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

struct CylinderRegion:Region
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

struct CylindricalShellRegion:Region
{
	tw::Float innerRadius,outerRadius;
	CylindricalShellRegion(std::vector<Region*>& ml) : Region(ml)
	{
		rgnType = cylindricalShellRegion;
		innerRadius = 1.0;
		outerRadius = 2.0;
		rbox = tw::vec3(outerRadius,outerRadius,tw::big_pos);
	}
	virtual bool Inside(const tw::vec3& pos,const MetricSpace& ds) const
	{
		tw::Float rho;
		tw::vec3 p = pos - center;
		orientation.ExpressInBasis(&p);
		rho = sqrt(p.x*p.x + p.y*p.y);
		return complement ^ ( rho < outerRadius && rho > innerRadius && sqr(p.z) < sqr(rbox.z));
	}
	virtual void ReadInputFileDirective(std::stringstream& inputString,const std::string& command);
	virtual void ReadData(std::ifstream& inFile);
	virtual void WriteData(std::ofstream& outFile);
};

struct RoundedCylinderRegion:Region
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

struct EllipsoidRegion:Region
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

struct TrueSphere:CircRegion
{
	TrueSphere(std::vector<Region*>& ml) : CircRegion(ml)
	{
		rgnType = trueSphereRegion;
	}
	virtual bool Inside(const tw::vec3& pos,const MetricSpace& ds) const;
};

struct BoxArrayRegion:Region
{
	tw::vec3 size,spacing;

	BoxArrayRegion(std::vector<Region*>& ml) : Region(ml)
	{
		rgnType = boxArrayRegion;
		size = tw::vec3(1,1,1);
		spacing = tw::vec3(2,2,2);
	}
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
	virtual void ReadInputFileDirective(std::stringstream& inputString,const std::string& command);
	virtual void ReadData(std::ifstream& inFile);
	virtual void WriteData(std::ofstream& outFile);
};

struct TorusRegion:Region
{
	tw::Float majorRadius,minorRadius;
	TorusRegion(std::vector<Region*>& ml) : Region(ml)
	{
		rgnType = torusRegion;
		majorRadius = 1.0;
		minorRadius = 0.1;
		rbox = tw::vec3(majorRadius+minorRadius,majorRadius+minorRadius,minorRadius);
	}
	virtual bool Inside(const tw::vec3& pos,const MetricSpace& ds) const
	{
		tw::Float rho;
		tw::vec3 p = pos - center;
		orientation.ExpressInBasis(&p);
		rho = sqrt(p.x*p.x + p.y*p.y);
		return complement ^ ( rho > majorRadius-minorRadius && rho < majorRadius+minorRadius &&
			sqr(p.z) < sqr(minorRadius)-sqr(rho-majorRadius));
	}
	virtual void ReadInputFileDirective(std::stringstream& inputString,const std::string& command);
	virtual void ReadData(std::ifstream& inFile);
	virtual void WriteData(std::ofstream& outFile);
};

struct ConeRegion:Region
{
	tw::Float majorRadius,minorRadius;
	ConeRegion(std::vector<Region*>& ml) : Region(ml)
	{
		rgnType = coneRegion;
		majorRadius = 1.0;
		minorRadius = 0.1;
		rbox = tw::vec3(majorRadius,majorRadius,tw::big_pos);
	}
	virtual bool Inside(const tw::vec3& pos,const MetricSpace& ds) const
	{
		tw::Float rho,rOfz;
		tw::vec3 p = pos - center;
		orientation.ExpressInBasis(&p);
		rho = sqrt(p.x*p.x + p.y*p.y);
		rOfz = minorRadius - (p.z - rbox.z)*(majorRadius-minorRadius)/(2.0*rbox.z);
		return complement ^ ( rho < rOfz && sqr(p.z) < sqr(rbox.z));
	}
	virtual void ReadInputFileDirective(std::stringstream& inputString,const std::string& command);
	virtual void ReadData(std::ifstream& inFile);
	virtual void WriteData(std::ofstream& outFile);
};

struct TangentOgiveRegion:Region
{
	tw::Float tipRadius,bodyRadius;
	TangentOgiveRegion(std::vector<Region*>& ml) : Region(ml)
	{
		rgnType = tangentOgiveRegion;
		tipRadius = 1.0;
		bodyRadius = 5.0;
		rbox = tw::vec3(bodyRadius,bodyRadius,tw::big_pos);
}
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
	virtual void ReadInputFileDirective(std::stringstream& inputString,const std::string& command);
	virtual void ReadData(std::ifstream& inFile);
	virtual void WriteData(std::ofstream& outFile);
};
