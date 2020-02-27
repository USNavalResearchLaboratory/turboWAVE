#include "meta_base.h"

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

Region* Region::CreateObjectFromFile(std::vector<Region*>& ml,std::ifstream& inFile)
{
	Region *ans;
	regionSpec rgnType;
	inFile.read((char *)&rgnType,sizeof(regionSpec));
	switch (rgnType)
	{
		case baseRegion:
			ans = new Region(ml);
			break;
		case entireRegion:
			ans = new EntireRegion(ml);
			break;
		case rectRegion:
			ans = new RectRegion(ml);
			break;
		case prismRegion:
			ans = new PrismRegion(ml);
			break;
		case circRegion:
			ans = new CircRegion(ml);
			break;
		case cylinderRegion:
			ans = new CylinderRegion(ml);
			break;
		case roundedCylinderRegion:
			ans = new RoundedCylinderRegion(ml);
			break;
		case ellipsoidRegion:
			ans = new EllipsoidRegion(ml);
			break;
		case trueSphereRegion:
			ans = new TrueSphere(ml);
			break;
		case boxArrayRegion:
			ans = new BoxArrayRegion(ml);
			break;
		case torusRegion:
			ans = new TorusRegion(ml);
			break;
		case coneRegion:
			ans = new ConeRegion(ml);
			break;
		case tangentOgiveRegion:
			ans = new TangentOgiveRegion(ml);
			break;
		case cylindricalShellRegion:
			ans = new CylindricalShellRegion(ml);
			break;
	}
	ans->ReadData(inFile);
	return ans;
}

Region* Region::ReadRegion(std::vector<Region*>& ml,Region *curr,std::stringstream& source,const std::string& com)
{
	// This should go in favor of standard directive reader
	std::string word;
	if (com=="clipping") // eg, clipping region = region1
	{
		source >> word >> word >> word;
		std::cout << word << std::endl;
		tw::input::StripQuotes(word);
		curr = FindRegion(ml,word);
		if (curr==NULL)
			throw tw::FatalError("Could not read region <"+word+">.");
	}
	return curr;
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
		inputString >> com;
		ReadInputFileDirective(inputString,com);
	} while (com!="}");
}

void Region::ReadInputFileDirective(std::stringstream& inputString,const std::string& com)
{
	std::string word;
	tw::Float theta;

	if (com=="bounds")
	{
		tw::vec3 r0,r1;
		inputString >> word >> r0.x >> r1.x >> r0.y >> r1.y >> r0.z >> r1.z;
		if (r0.x>r1.x) std::swap(r0.x,r1.x);
		if (r0.y>r1.y) std::swap(r0.y,r1.y);
		if (r0.z>r1.z) std::swap(r0.z,r1.z);
		center = 0.5*(r0 + r1);
		rbox = 0.5*(r1 - r0);
	}
	if (com=="center")
	{
		throw tw::FatalError("Direct assignment to center is illegal.  Use translation instead.");
	}
	if (com=="radius")
	{
		inputString >> word >> rbox.x;
		if (rbox.y > 0.99*tw::big_pos)
			rbox.y = rbox.x;
		if (rbox.z > 0.99*tw::big_pos)
			rbox.z = rbox.x;
	}
	if (com=="length")
	{
		inputString >> word >> rbox.z;
		rbox.z *= 0.5;
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
			if (word!="}")
			{
				tw::input::StripQuotes(word);
				curr = FindRegion(masterList,word);
				composite.push_back(curr);
			}
		} while (word!="}");
	}

	if (com=="move") // eg, move with window = true
	{
		inputString >> word >> word >> word >> word;
		moveWithWindow = (word=="on" || word=="true" || word=="yes");
	}
	if (com=="translation") // eg, translation = 1 0 0
	{
		tw::vec3 dr;
		inputString >> word >> dr.x >> dr.y >> dr.z;
		center += dr;
	}
	if (com=="rotation")
	{
		inputString >> word >> word; // eg, rotation about x = 45
		if (word=="x")
		{
			inputString >> word >> theta;
			orientation.u.RotateX(theta);
			orientation.v.RotateX(theta);
			orientation.w.RotateX(theta);
			center.RotateX(theta);
		}
		if (word=="y")
		{
			inputString >> word >> theta;
			orientation.u.RotateY(theta);
			orientation.v.RotateY(theta);
			orientation.w.RotateY(theta);
			center.RotateY(theta);
		}
		if (word=="z")
		{
			inputString >> word >> theta;
			orientation.u.RotateZ(theta);
			orientation.v.RotateZ(theta);
			orientation.w.RotateZ(theta);
			center.RotateZ(theta);
		}
	}

	if (com=="complement")
	{
		inputString >> word >> word; // eg, complement = true
		complement = (word=="true" || word=="yes" || word=="on");
	}
}

void Region::ReadData(std::ifstream& inFile)
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
	tw::Int num,idx;
	inFile.read((char *)&num,sizeof(tw::Int));

	for (tw::Int i=0;i<num;i++)
	{
		inFile.read((char *)&idx,sizeof(tw::Int));
		composite.push_back(masterList[idx]);
	}
	inFile >> name;
	inFile.ignore();
}

void Region::WriteData(std::ofstream& outFile)
{
	outFile.write((char *)&rgnType,sizeof(regionSpec));
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
	tw::Int num,idx;
	num = composite.size();
	outFile.write((char *)&num,sizeof(tw::Int));
	for (tw::Int i=0;i<num;i++)
	{
		idx = std::find(masterList.begin(),masterList.end(),composite[i]) - masterList.begin();
		outFile.write((char *)&idx,sizeof(tw::Int));
	}
	outFile << name << " ";
}

bool TrueSphere::Inside(const tw::vec3& pos,const MetricSpace& ds) const
{
	tw::vec3 c_cart = center;
	tw::vec3 p_cart = pos;
	ds.CurvilinearToCartesian(&c_cart);
	ds.CurvilinearToCartesian(&p_cart);
	return complement ^ (Norm(p_cart-c_cart)<sqr(rbox.x));
}

void BoxArrayRegion::ReadInputFileDirective(std::stringstream& inputString,const std::string& com)
{
	Region::ReadInputFileDirective(inputString,com);
	std::string word;
	if (com=="size")
		inputString >> word >> size.x >> size.y >> size.z;
	if (com=="spacing")
		inputString >> word >> spacing.x >> spacing.y >> spacing.z;
}

void BoxArrayRegion::ReadData(std::ifstream& inFile)
{
	Region::ReadData(inFile);
	inFile.read((char *)&size,sizeof(tw::vec3));
	inFile.read((char *)&spacing,sizeof(tw::vec3));
}

void BoxArrayRegion::WriteData(std::ofstream& outFile)
{
	Region::WriteData(outFile);
	outFile.write((char *)&size,sizeof(tw::vec3));
	outFile.write((char *)&spacing,sizeof(tw::vec3));
}

void TorusRegion::ReadInputFileDirective(std::stringstream& inputString,const std::string& com)
{
	Region::ReadInputFileDirective(inputString,com);
	std::string word;
	if (com=="minor")
		inputString >> word >> word >> minorRadius;
	if (com=="major")
		inputString >> word >> word >> majorRadius;
	rbox = tw::vec3(majorRadius+minorRadius,majorRadius+minorRadius,minorRadius);
}

void TorusRegion::ReadData(std::ifstream& inFile)
{
	Region::ReadData(inFile);
	inFile.read((char *)&minorRadius,sizeof(tw::Float));
	inFile.read((char *)&majorRadius,sizeof(tw::Float));
}

void TorusRegion::WriteData(std::ofstream& outFile)
{
	Region::WriteData(outFile);
	outFile.write((char *)&minorRadius,sizeof(tw::Float));
	outFile.write((char *)&majorRadius,sizeof(tw::Float));
}

void ConeRegion::ReadInputFileDirective(std::stringstream& inputString,const std::string& com)
{
	Region::ReadInputFileDirective(inputString,com);
	std::string word;
	if (com=="tip" || com=="minor")
		inputString >> word >> word >> minorRadius;
	if (com=="base" || com=="major")
	{
		inputString >> word >> word >> majorRadius;
		rbox.x = rbox.y = majorRadius;
	}
}

void ConeRegion::ReadData(std::ifstream& inFile)
{
	Region::ReadData(inFile);
	inFile.read((char *)&minorRadius,sizeof(tw::Float));
	inFile.read((char *)&majorRadius,sizeof(tw::Float));
}

void ConeRegion::WriteData(std::ofstream& outFile)
{
	Region::WriteData(outFile);
	outFile.write((char *)&minorRadius,sizeof(tw::Float));
	outFile.write((char *)&majorRadius,sizeof(tw::Float));
}

void TangentOgiveRegion::ReadInputFileDirective(std::stringstream& inputString,const std::string& com)
{
	Region::ReadInputFileDirective(inputString,com);
	std::string word;
	if (com=="tip" || com=="minor")
		inputString >> word >> word >> tipRadius;
	if (com=="base" || com=="major" || com=="body")
	{
		inputString >> word >> word >> bodyRadius;
		rbox.x = rbox.y = bodyRadius;
	}
}

void TangentOgiveRegion::ReadData(std::ifstream& inFile)
{
	Region::ReadData(inFile);
	inFile.read((char *)&tipRadius,sizeof(tw::Float));
	inFile.read((char *)&bodyRadius,sizeof(tw::Float));
}

void TangentOgiveRegion::WriteData(std::ofstream& outFile)
{
	Region::WriteData(outFile);
	outFile.write((char *)&tipRadius,sizeof(tw::Float));
	outFile.write((char *)&bodyRadius,sizeof(tw::Float));
}

void CylindricalShellRegion::ReadInputFileDirective(std::stringstream& inputString,const std::string& com)
{
	Region::ReadInputFileDirective(inputString,com);
	std::string word;
	if (com=="inner" || com=="minor")
		inputString >> word >> word >> innerRadius;
	if (com=="outer" || com=="major")
	{
		inputString >> word >> word >> outerRadius;
		rbox.x = rbox.y = outerRadius;
	}
}

void CylindricalShellRegion::ReadData(std::ifstream& inFile)
{
	Region::ReadData(inFile);
	inFile.read((char *)&innerRadius,sizeof(tw::Float));
	inFile.read((char *)&outerRadius,sizeof(tw::Float));
}

void CylindricalShellRegion::WriteData(std::ofstream& outFile)
{
	Region::WriteData(outFile);
	outFile.write((char *)&innerRadius,sizeof(tw::Float));
	outFile.write((char *)&outerRadius,sizeof(tw::Float));
}
