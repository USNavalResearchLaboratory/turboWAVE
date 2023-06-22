#include "meta_base.h"

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
			low = task->GlobalCellIndex(localBounds[2*d],d+1);
		else
			low = task->globalCells[d+1];
		globalBounds[2*d] = task->strip[d+1].GetMin(low);
	}

	for (d=0;d<3;d++)
	{
		if (localBounds[2*d+1]>=1 && localBounds[2*d+1]<=dims[d])
			high = task->GlobalCellIndex(localBounds[2*d+1],d+1);
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

void Region::ReadInputFileBlock(std::stringstream& inputString)
{
	std::string com;
	do
	{
		com = directives.ReadNext(inputString);
		if (com=="tw::EOF")
			throw tw::FatalError("Encountered EOF while processing <"+name+">.");
		ReadInputFileDirective(inputString,com);
	} while (com!="}");
	directives.ThrowErrorIfMissingKeys(name);
}

void Region::ReadInputFileDirective(std::stringstream& inputString,const std::string& com)
{
	std::string word;

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
		inputString >> word;
		if (word!="=")
			throw tw::FatalError("Expected <=> after <elements>.");
		inputString >> word;
		if (word!="{")
			throw tw::FatalError("Expected <{> at start of list.");
		do {
			inputString >> word;
			if (inputString.eof())
				throw tw::FatalError("Encountered EOF while processing <"+name+">.");
			if (word!="}")
			{
				tw::input::StripQuotes(word);
				curr = FindRegion(masterList,word);
				composite.push_back(curr);
			}
		} while (word!="}");
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

void BoxArrayRegion::ReadInputFileDirective(std::stringstream& inputString,const std::string& com)
{
	Region::ReadInputFileDirective(inputString,com);
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

void TorusRegion::ReadInputFileDirective(std::stringstream& inputString,const std::string& com)
{
	Region::ReadInputFileDirective(inputString,com);
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

void ConeRegion::ReadInputFileDirective(std::stringstream& inputString,const std::string& com)
{
	Region::ReadInputFileDirective(inputString,com);
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

void TangentOgiveRegion::ReadInputFileDirective(std::stringstream& inputString,const std::string& com)
{
	Region::ReadInputFileDirective(inputString,com);
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

void CylindricalShellRegion::ReadInputFileDirective(std::stringstream& inputString,const std::string& com)
{
	Region::ReadInputFileDirective(inputString,com);
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
