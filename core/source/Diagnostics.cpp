#include "sim.h"
//#define ALT_BOX_DIAG

void WriteDVHeader(std::ofstream& outFile,tw::Int version,tw::Int xDim,tw::Int yDim,tw::Int zDim,float x0,float x1,float y0,float y1,float z0,float z1)
{
	int32_t xDim32 = xDim;
	int32_t yDim32 = yDim;
	int32_t zDim32 = zDim;

	switch (version)
	{
		case 1:
			WriteBigEndian((char *)&xDim32,sizeof(int32_t),0,outFile);
			WriteBigEndian((char *)&yDim32,sizeof(int32_t),0,outFile);
			WriteBigEndian((char *)&zDim32,sizeof(int32_t),0,outFile);
			break;
		case 2:
			const char buffer[17] = "DataViewer 2.0.0";
			outFile.write(buffer,16);
			WriteBigEndian((char *)&xDim32,sizeof(int32_t),0,outFile);
			WriteBigEndian((char *)&yDim32,sizeof(int32_t),0,outFile);
			WriteBigEndian((char *)&zDim32,sizeof(int32_t),0,outFile);
			WriteBigEndian((char*)&x0,sizeof(float),0,outFile);
			WriteBigEndian((char*)&x1,sizeof(float),0,outFile);
			WriteBigEndian((char*)&y0,sizeof(float),0,outFile);
			WriteBigEndian((char*)&y1,sizeof(float),0,outFile);
			WriteBigEndian((char*)&z0,sizeof(float),0,outFile);
			WriteBigEndian((char*)&z1,sizeof(float),0,outFile);
			break;
	}
}

////////////////////////////
// DIAGNOSTIC DESCRIPTORS //
////////////////////////////


DiagnosticDescriptor::DiagnosticDescriptor(std::vector<Region*>& rlist) : rgnList(rlist)
{
	period = 0;
	t0 = 0.0;
	t1 = tw::big_pos;
	tRef = tw::big_neg;
	timePeriod = 0.0;
	filename = "diagnostic";
	theRgn = rgnList[0];
	moveWithWindow = true;
}

bool DiagnosticDescriptor::WriteThisStep(tw::Float elapsedTime,tw::Float dt,tw::Int stepNow)
{
	if (timePeriod!=0.0)
	{
		if (tRef==elapsedTime)
			return true;
		if (elapsedTime >= t0 && elapsedTime <= t1)
		{
			tRef = elapsedTime;
			t0 += timePeriod;
			return true;
		}
		return false;
	}

	if (period!=0)
	{
		tw::Int startStep = tw::Int((t0+dt)/dt);
		if ( (elapsedTime>=t0) && (elapsedTime<=t1) && ((stepNow - startStep + 1) % period == 0) )
			return true;
		else
			return false;
	}

	return false;
}

void DiagnosticDescriptor::ReadData(std::ifstream& inFile)
{
	tw::Int rgnIndex;
	inFile.read((char *)&rgnIndex,sizeof(tw::Int));
	theRgn = rgnList[rgnIndex];
	inFile.read((char *)&period,sizeof(tw::Int));
	inFile.read((char *)&timePeriod,sizeof(tw::Float));
	inFile.read((char *)&tRef,sizeof(tw::Float));
	inFile.read((char *)&t0,sizeof(tw::Float));
	inFile.read((char *)&t1,sizeof(tw::Float));
	inFile >> filename;
	inFile.ignore();
}

void DiagnosticDescriptor::WriteData(std::ofstream& outFile)
{
	tw::Int rgnIndex = std::find(rgnList.begin(),rgnList.end(),theRgn) - rgnList.begin();
	outFile.write((char *)&rgnIndex,sizeof(tw::Int));
	outFile.write((char *)&period,sizeof(tw::Int));
	outFile.write((char *)&timePeriod,sizeof(tw::Float));
	outFile.write((char *)&tRef,sizeof(tw::Float));
	outFile.write((char *)&t0,sizeof(tw::Float));
	outFile.write((char *)&t1,sizeof(tw::Float));
	outFile << filename << " ";
}

EnergySeriesDescriptor::EnergySeriesDescriptor(std::vector<Region*>& rlist) : DiagnosticDescriptor(rlist)
{
	numSigFigs = 6;
	filename = "energy";
}

void EnergySeriesDescriptor::ReadInputFile(std::stringstream& inputString)
{
	std::string word;
	do
	{
		inputString >> word;
		theRgn = Region::ReadRegion(rgnList,theRgn,inputString,word);
		if (word=="period")
			inputString >> word >> period;
		if (word=="time") // eg, time period = 100.0
			inputString >> word >> word >> timePeriod;
		if (word=="filename")
			inputString >> word >> filename;
		if (word=="precision")
			inputString >> word >> numSigFigs;
		if (word=="t0")
			inputString >> word >> t0;
		if (word=="t1")
			inputString >> word >> t1;
	} while (word!="}");
}

void EnergySeriesDescriptor::ReadData(std::ifstream& inFile)
{
	DiagnosticDescriptor::ReadData(inFile);
	inFile.read((char *)&numSigFigs,sizeof(tw::Int));
}

void EnergySeriesDescriptor::WriteData(std::ofstream& outFile)
{
	DiagnosticDescriptor::WriteData(outFile);
	outFile.write((char *)&numSigFigs,sizeof(tw::Int));
}

PointSeriesDescriptor::PointSeriesDescriptor(std::vector<Region*>& rlist) : DiagnosticDescriptor(rlist)
{
	thePoint = tw::vec3(1.0);
	filename = "pt1";
}

void PointSeriesDescriptor::ReadInputFile(std::stringstream& inputString)
{
	std::string word;
	do
	{
		inputString >> word;
		if (word=="move") // eg, move with window = true
		{
			inputString >> word >> word >> word >> word;
			moveWithWindow = (word=="on" || word=="true" || word=="yes");
		}
		if (word=="point")
			inputString >> word >> thePoint.x >> thePoint.y >> thePoint.z;
		if (word=="period")
			inputString >> word >> period;
		if (word=="time") // eg, time period = 100.0
			inputString >> word >> word >> timePeriod;
		if (word=="filename")
			inputString >> word >> filename;
		if (word=="t0")
			inputString >> word >> t0;
		if (word=="t1")
			inputString >> word >> t1;
	} while (word!="}");
}

void PointSeriesDescriptor::ReadData(std::ifstream& inFile)
{
	DiagnosticDescriptor::ReadData(inFile);
	inFile.read((char *)&thePoint,sizeof(tw::vec3));
}

void PointSeriesDescriptor::WriteData(std::ofstream& outFile)
{
	DiagnosticDescriptor::WriteData(outFile);
	outFile.write((char *)&thePoint,sizeof(tw::vec3));
}

GridDataDescriptor::GridDataDescriptor(std::vector<Region*>& rlist) : DiagnosticDescriptor(rlist)
{
	average = false;
	xSkip = ySkip = zSkip = 1;
	filename = "full";
}

void GridDataDescriptor::ReadInputFile(std::stringstream& inputString)
{
	std::string word;
	do
	{
		inputString >> word;
		theRgn = Region::ReadRegion(rgnList,theRgn,inputString,word);
		if (word=="average")
		{
			inputString >> word >> word;
			average = (word=="yes" || word=="true" || word=="on");
		}
		if (word=="skip")
			inputString >> word >> xSkip >> ySkip >> zSkip;
		if (word=="period")
			inputString >> word >> period;
		if (word=="time") // eg, time period = 100.0
			inputString >> word >> word >> timePeriod;
		if (word=="filename")
			inputString >> word >> filename;
		if (word=="t0")
			inputString >> word >> t0;
		if (word=="t1")
			inputString >> word >> t1;
	} while (word!="}");
}

void GridDataDescriptor::ReadData(std::ifstream& inFile)
{
	DiagnosticDescriptor::ReadData(inFile);
	inFile.read((char *)&average,sizeof(bool));
	inFile.read((char *)&xSkip,sizeof(tw::Int));
	inFile.read((char *)&ySkip,sizeof(tw::Int));
	inFile.read((char *)&zSkip,sizeof(tw::Int));
}

void GridDataDescriptor::WriteData(std::ofstream& outFile)
{
	DiagnosticDescriptor::WriteData(outFile);
	outFile.write((char *)&average,sizeof(bool));
	outFile.write((char *)&xSkip,sizeof(tw::Int));
	outFile.write((char *)&ySkip,sizeof(tw::Int));
	outFile.write((char *)&zSkip,sizeof(tw::Int));
}

ParticleDetectorDescriptor::ParticleDetectorDescriptor(std::vector<Region*>& rlist) : DiagnosticDescriptor(rlist)
{
	minGamma = 1.0;
	filename = "det";
}

void ParticleDetectorDescriptor::ReadInputFile(std::stringstream& inputString)
{
	std::string word;
	do
	{
		inputString >> word;
		theRgn = Region::ReadRegion(rgnList,theRgn,inputString,word);
		if (word=="minimum")
			inputString >> word >> word >> minGamma;
		if (word=="filename")
			inputString >> word >> filename;
		if (word=="position")
			inputString >> word >> position.x >> position.y >> position.z;
		if (word=="t0")
			inputString >> word >> t0;
		if (word=="t1")
			inputString >> word >> t1;
	} while (word!="}");
}

void ParticleDetectorDescriptor::ReadData(std::ifstream& inFile)
{
	DiagnosticDescriptor::ReadData(inFile);
	inFile.read((char *)&minGamma,sizeof(tw::Float));
	inFile.read((char *)&position,sizeof(tw::vec3));
}

void ParticleDetectorDescriptor::WriteData(std::ofstream& outFile)
{
	DiagnosticDescriptor::WriteData(outFile);
	outFile.write((char *)&minGamma,sizeof(tw::Float));
	outFile.write((char *)&position,sizeof(tw::vec3));
}

void ParticleDetectorDescriptor::WriteRecords(std::valarray<float>& data)
{
	tw::Int i,recs;
	recs = data.size()/7;
	for (i=0;i<recs;i++)
		outFile << std::scientific << data[i*7+0] << " " << data[i*7+1] << " " << data[i*7+2] << " " <<
			data[i*7+3] << " " << data[i*7+4] << " " << data[i*7+5] << " " << data[i*7+6] << std::endl;

}

ParticleOrbitDescriptor::ParticleOrbitDescriptor(std::vector<Region*>& rlist) : DiagnosticDescriptor(rlist)
{
	filename = "par";
	minGamma = 1.0;
}

void ParticleOrbitDescriptor::ReadInputFile(std::stringstream& inputString)
{
	std::string word;
	do
	{
		inputString >> word;
		theRgn = Region::ReadRegion(rgnList,theRgn,inputString,word);
		if (word=="minimum")
			inputString >> word >> word >> minGamma;
		if (word=="period")
			inputString >> word >> period;
		if (word=="time") // eg, time period = 100.0
			inputString >> word >> word >> timePeriod;
		if (word=="filename")
			inputString >> word >> filename;
		if (word=="t0")
			inputString >> word >> t0;
		if (word=="t1")
			inputString >> word >> t1;
	} while (word!="}");
}

void ParticleOrbitDescriptor::ReadData(std::ifstream& inFile)
{
	DiagnosticDescriptor::ReadData(inFile);
	inFile.read((char *)&minGamma,sizeof(tw::Float));
}

void ParticleOrbitDescriptor::WriteData(std::ofstream& outFile)
{
	DiagnosticDescriptor::WriteData(outFile);
	outFile.write((char *)&minGamma,sizeof(tw::Float));
}

PhaseSpaceDescriptor::PhaseSpaceDescriptor(std::vector<Region*>& rlist) : DiagnosticDescriptor(rlist)
{
	hAxis = "x";
	vAxis = "px";
	filename = "xpx";
	min = tw::vec3(0.0,0.0,0.0);
	max = tw::vec3(1.0,1.0,1.0);
	hDim = vDim = 100;
}

void PhaseSpaceDescriptor::SetupGeometry(tw_geometry theGeometry)
{
	// This is called in Species::Initialize.  The Volume factors don't need to be saved in restart file.
	tw::Float dh,dv;
	horVolFactor.resize(hDim);
	verVolFactor.resize(vDim);
	dh = (max.x - min.x) / tw::Float(hDim);
	dv = (max.y - min.y) / tw::Float(vDim);
	horVolFactor = dh;
	verVolFactor = dv;
	// Currently all phase space data is cartesian.
	// particle coordinates are transformed to cartesian in Speces::CustomDiagnose
	// particle momentum is always in cartesian, so no transformation is necessary
	/*tw::Int i,j;
	if (theGeometry==cylindrical)
	{
		if (hAxis=="x")
			for (i=0;i<hDim;i++)
				horVolFactor[i] = pi*(sqr(tw::Float(i+1)*dh) - sqr(tw::Float(i)*dh)); // assumes axisymmetry
		if (vAxis=="x")
			for (j=0;j<vDim;j++)
				verVolFactor[j] = pi*(sqr(tw::Float(j+1)*dv) - sqr(tw::Float(j)*dv)); // assumes axisymmetry
	}*/
}

void PhaseSpaceDescriptor::ReadInputFile(std::stringstream& inputString)
{
	std::string word;
	do
	{
		inputString >> word;
		theRgn = Region::ReadRegion(rgnList,theRgn,inputString,word);
		if (word=="period")
			inputString >> word >> period;
		if (word=="time") // eg, time period = 100.0
			inputString >> word >> word >> timePeriod;
		if (word=="abcissa")
			inputString >> word >> hAxis;
		if (word=="ordinate")
			inputString >> word >> vAxis;
		if (word=="minimum")
			inputString >> word >> min.x >> min.y;
		if (word=="maximum")
			inputString >> word >> max.x >> max.y;
		if (word=="dimensions")
			inputString >> word >> hDim >> vDim;
		if (word=="filename")
			inputString >> word >> filename;
		if (word=="t0")
			inputString >> word >> t0;
		if (word=="t1")
			inputString >> word >> t1;
	} while (word!="}");
}

void PhaseSpaceDescriptor::ReadData(std::ifstream& inFile)
{
	DiagnosticDescriptor::ReadData(inFile);
	inFile.read((char *)&min,sizeof(tw::vec3));
	inFile.read((char *)&max,sizeof(tw::vec3));
	inFile.read((char *)&hDim,sizeof(tw::Int));
	inFile.read((char *)&vDim,sizeof(tw::Int));
	inFile >> hAxis >> vAxis;
	inFile.ignore();
}

void PhaseSpaceDescriptor::WriteData(std::ofstream& outFile)
{
	DiagnosticDescriptor::WriteData(outFile);
	outFile.write((char *)&min,sizeof(tw::vec3));
	outFile.write((char *)&max,sizeof(tw::vec3));
	outFile.write((char *)&hDim,sizeof(tw::Int));
	outFile.write((char *)&vDim,sizeof(tw::Int));
	outFile << hAxis << " " << vAxis << " ";
}


FarFieldDetectorDescriptor::FarFieldDetectorDescriptor(Grid *theGrid) : DiagnosticDescriptor(theGrid->clippingRegion)
{
	filename = "far-field";
	radius = 100;
	theta0 = 0;
	theta1 = 3.14;
	thetaPts = 10;
	phi0 = -3.14;
	phi1 = 3.14;
	phiPts = 1;
	t0 = 0;
	t1 = 100;
	timePts = 100;
	this->theGrid = theGrid;
}

void FarFieldDetectorDescriptor::ReadInputFile(std::stringstream& inputString)
{
	std::string word;
	do
	{
		inputString >> word;
		theRgn = Region::ReadRegion(rgnList,theRgn,inputString,word);
		if (word=="period")
			inputString >> word >> period;
		if (word=="filename")
			inputString >> word >> filename;
		if (word=="radius")
			inputString >> word >> radius;
		if (word=="theta") // eg, theta , start 0 , step 0.01 , points 100
		{
			inputString >> word >> theta0 >> word >> theta1 >> word >> thetaPts;
			theta0 -= 0.5*theta1;
			theta1 = theta0 + theta1*tw::Float(thetaPts);
		}
		if (word=="phi") // eg, phi , start 0 , step 0.01 , points 100
		{
			inputString >> word >> phi0 >> word >> phi1 >> word >> phiPts;
			phi0 -= 0.5*phi1;
			phi1 = phi0 + phi1*tw::Float(phiPts);
		}
		if (word=="time")
		{
			inputString >> word;
			if (word=="start") // eg, time , start 0 , step 0.01 , points 100
			{
				inputString >> t0 >> word >> t1 >> word >> timePts;
				t0 -= 0.5*t1;
				t1 = t0 + t1*tw::Float(timePts);
			}
			if (word=="period") // eg, time period = 100.0
				inputString >> word >> timePeriod;
		}
	} while (word!="}");

	DiscreteSpace farFieldLayout;
	farFieldLayout.Resize(timePts,thetaPts,phiPts,tw::vec3(0.0,0.0,0.0),tw::vec3(t1-t0,theta1-theta0,phi1-phi0));
	A.Initialize(farFieldLayout,theGrid);
}

void FarFieldDetectorDescriptor::ReadData(std::ifstream& inFile)
{
	DiagnosticDescriptor::ReadData(inFile);
	inFile.read((char *)&radius,sizeof(tw::Float));
	inFile.read((char *)&theta0,sizeof(tw::Float));
	inFile.read((char *)&theta1,sizeof(tw::Float));
	inFile.read((char *)&thetaPts,sizeof(tw::Int));
	inFile.read((char *)&phi0,sizeof(tw::Float));
	inFile.read((char *)&phi1,sizeof(tw::Float));
	inFile.read((char *)&phiPts,sizeof(tw::Int));
	inFile.read((char *)&timePts,sizeof(tw::Int));
	A.ReadData(inFile);
}

void FarFieldDetectorDescriptor::WriteData(std::ofstream& outFile)
{
	DiagnosticDescriptor::WriteData(outFile);
	outFile.write((char *)&radius,sizeof(tw::Float));
	outFile.write((char *)&theta0,sizeof(tw::Float));
	outFile.write((char *)&theta1,sizeof(tw::Float));
	outFile.write((char *)&thetaPts,sizeof(tw::Int));
	outFile.write((char *)&phi0,sizeof(tw::Float));
	outFile.write((char *)&phi1,sizeof(tw::Float));
	outFile.write((char *)&phiPts,sizeof(tw::Int));
	outFile.write((char *)&timePts,sizeof(tw::Int));
	A.WriteData(outFile);
}

void FarFieldDetectorDescriptor::AccumulateField(Field& J4)
{
	tw::Int i,j,k,ip,jp;
	tw::vec3 rVec,rp;
	std::valarray<tw::Float> j4(4);
	tw::Float dS;
	weights_3D w;
	tw::Float thetaNow,phiNow,tNow;
	tw::Float dt = (t1 - t0)/tw::Float(timePts);
	tw::Float dtheta = (theta1 - theta0)/tw::Float(thetaPts);
	tw::Float dphi = (phi1 - phi0)/tw::Float(phiPts);

	tw::Float zmin = theGrid->X(1,3) - 0.5*theGrid->dX(1,3);
	tw::Float zmax = theGrid->X(theGrid->Dim(3),3) + 0.5*theGrid->dX(theGrid->Dim(3),3);
	tw::Float tp = theGrid->elapsedTime;
	tw::Float dtau = timePeriod==0.0 ? theGrid->dt * tw::Float(period) : timePeriod;

	phiNow = phi0 + 0.5*dphi;
	for (k=1;k<=phiPts;k++)
	{
		thetaNow = theta0 + 0.5*dtheta;
		for (j=1;j<=thetaPts;j++)
		{
			rVec.x = sin(thetaNow)*cos(phiNow);
			rVec.y = sin(thetaNow)*sin(phiNow);
			rVec.z = cos(thetaNow);
			tNow = t0 + 0.5*dt;
			for (i=1;i<=timePts;i++)
			{
				for (jp=1;jp<=theGrid->Dim(2);jp++)
					for (ip=1;ip<=theGrid->Dim(1);ip++)
					{
						rp.x = theGrid->X(ip,1);
						rp.y = theGrid->X(jp,2);;
						dS = theGrid->dX(ip,1)*theGrid->dX(jp,2)/rVec.z;
						rp.z = ((tp - tNow) - rp.x*rVec.x - rp.y*rVec.y)/rVec.z;
						if (rp.z > zmin && rp.z < zmax)
						{
							theGrid->GetWeights(&w,rp);
							J4.Interpolate(j4,w);
							A(i,j,k) += tw::vec3(j4[1],j4[2],j4[3]) * dS * dtau / radius;
						}
					}
				tNow += dt;
			}
			thetaNow += dtheta;
		}
		phiNow += dphi;
	}
}

//////////////////////
// GRID DIAGNOSTICS //
//////////////////////



#ifndef ALT_BOX_DIAG

void Grid::GetGlobalBoxDataIndexing(GridDataDescriptor* theBox,tw::Int pts[4],tw::Int glb[6],tw::Int skip[4])
{
	// On output, pts has the number of cells, glb has the global index bounds, skip has the skip corrected for ignorable dimensions
	tw::Int s[4] = { 0 , theBox->xSkip , theBox->ySkip , theBox->zSkip };

	theBox->theRgn->GetGlobalCellBounds(&glb[0],&glb[1],&glb[2],&glb[3],&glb[4],&glb[5]);

	for (tw::Int ax=1;ax<=3;ax++)
	{
		const tw::Int lb = 2*ax-2;
		const tw::Int ub = 2*ax-1;
		skip[ax] = dim[ax]==1 ? 1 : s[ax];
		// Force global bounds into skipping sequence
		// The upper bound is adjusted to respect the lower bound
		pts[ax] = (glb[ub] - glb[lb] + 1)/skip[ax];
		glb[ub] = glb[lb] + skip[ax]*(pts[ax]-1);
	}
}

void Grid::GetLocalBoxDataIndexing(GridDataDescriptor* theBox,const tw::Int pts[4],const tw::Int glb[6],const tw::Int skip[4],tw::Int loc[6],const tw::Int coords[4])
{
	// On output, pts has the number of cells, glb has the global index bounds, loc has the local index bounds
	// coords is an input with the cartesian MPI domain indices (not necessarily the domain of execution)
	tw::Int i0;
	for (tw::Int ax=1;ax<=3;ax++)
	{
		const tw::Int cell0 = coords[ax]*localCells[ax];
		const tw::Int lb = 2*ax-2;
		const tw::Int ub = 2*ax-1;
		// Get the lower local cell bounds and force into global skipping sequence
		i0 = 1 + cell0;
		i0 = i0 < glb[lb] ? glb[lb] : i0;
		i0 -= cell0 + ((i0 - glb[lb]) % skip[ax]);
		if (i0<1) i0 += skip[ax];
		loc[lb] = i0;
		// Get the upper local cell bounds and force into global skipping sequence
		i0 = localCells[ax] + cell0;
		i0 = i0 > glb[ub] ? glb[ub] : i0;
		i0 -= cell0 + ((i0 - glb[lb]) % skip[ax]);
		loc[ub] = i0;
	}
}

void Grid::WriteBoxDataHeader(const std::string& quantity,GridDataDescriptor* theBox)
{
	if (appendMode && restarted)
		return;

	// Only the task at index (0,0,0) will write data

	if (strip[0].Get_rank()!=0)
		return;

	std::ofstream outFile;
	std::stringstream fileName;
	tw::Int pts[4],glb[6],skip[4];

	fileName.str("");
	if (theBox->filename=="full")
		fileName << quantity;
	else
		fileName << theBox->filename << "_" << quantity;
	fileName << ".dvdat";

	GetGlobalBoxDataIndexing(theBox,pts,glb,skip);

	outFile.open(fileName.str().c_str(),std::ios::binary);
	WriteDVHeader(outFile,2,pts[1],pts[2],pts[3],
		globalCorner.x + spacing.x*(glb[0]-0.5-0.5*skip[1]),
		globalCorner.x + spacing.x*(glb[1]-0.5+0.5*skip[1]),
		globalCorner.y + spacing.y*(glb[2]-0.5-0.5*skip[2]),
		globalCorner.y + spacing.y*(glb[3]-0.5+0.5*skip[2]),
		globalCorner.z + spacing.z*(glb[4]-0.5-0.5*skip[3]),
		globalCorner.z + spacing.z*(glb[5]-0.5+0.5*skip[3]));
	outFile.close();
}

void Grid::WriteBoxData(const std::string& quantity,GridDataDescriptor* theBox,tw::Float* theData,const tw::Int *stride)
{
	tw::Int i,j,k,ax,idom,jdom,kdom,buffSize,ready,cell0,i0,i1;
	std::valarray<float> buffer,gData;
	std::ofstream outFile;
	std::stringstream fileName;
	tw::Int s[4],coords[4],pts[4],loc[6],glb[6];
	tw::Int thisNode,curr,master;

	thisNode = strip[0].Get_rank();
	master = 0;

	fileName.str("");
	if (theBox->filename=="full")
		fileName << quantity;
	else
		fileName << theBox->filename << "_" << quantity;
	fileName << ".dvdat";

	GetGlobalBoxDataIndexing(theBox,pts,glb,s);

	if (thisNode==master)
		gData.resize(pts[1]*pts[2]*pts[3]);

	for (kdom=(glb[4]-1)/dim[3];kdom<=(glb[5]-1)/dim[3];kdom++)
		for (jdom=(glb[2]-1)/dim[2];jdom<=(glb[3]-1)/dim[2];jdom++)
			for (idom=(glb[0]-1)/dim[1];idom<=(glb[1]-1)/dim[1];idom++)
			{
				curr = strip[0].Cart_rank(idom,jdom,kdom); // domain being written out
				coords[1] = idom; coords[2] = jdom; coords[3] = kdom;
				GetLocalBoxDataIndexing(theBox,pts,glb,s,loc,coords);

				if (loc[0]<=dim[1] && loc[1]>=1 && loc[2]<=dim[2] && loc[3]>=1 && loc[4]<=dim[3] && loc[5]>=1)
				{
					if (curr==thisNode)
					{
						if (curr==master)
						{
							for (k=loc[4];k<=loc[5];k+=s[3])
								for (j=loc[2];j<=loc[3];j+=s[2])
									for (i=loc[0];i<=loc[1];i+=s[1])
									{
										i0 = (GlobalCellIndex(i,1) - glb[0])/s[1];
										i0 += pts[1]*(GlobalCellIndex(j,2) - glb[2])/s[2];
										i0 += pts[1]*pts[2]*(GlobalCellIndex(k,3) - glb[4])/s[3];
										gData[i0] = theData[i*stride[1] + j*stride[2] + k*stride[3]];
									}
						}
						else
						{
							strip[0].Recv(&ready,sizeof(tw::Int),master);
							buffSize = (loc[1]-loc[0]+s[1])*(loc[3]-loc[2]+s[2])*(loc[5]-loc[4]+s[3])/(s[1]*s[2]*s[3]);
							buffer.resize(buffSize);
							for (k=loc[4];k<=loc[5];k+=s[3])
								for (j=loc[2];j<=loc[3];j+=s[2])
									for (i=loc[0];i<=loc[1];i+=s[1])
									{
										i0 = (i - loc[0])/s[1];
										i0 += (loc[1]-loc[0]+s[1])*(j - loc[2])/(s[1]*s[2]);
										i0 += (loc[1]-loc[0]+s[1])*(loc[3]-loc[2]+s[2])*(k - loc[4])/(s[1]*s[2]*s[3]);
										buffer[i0] = theData[i*stride[1] + j*stride[2] + k*stride[3]];
									}
							strip[0].Send(&buffer[0],sizeof(float)*buffSize,master);
						}
					}
					if (thisNode==master && curr!=master)
					{
						strip[0].Send(&ready,sizeof(tw::Int),curr);
						buffSize = (loc[1]-loc[0]+s[1])*(loc[3]-loc[2]+s[2])*(loc[5]-loc[4]+s[3])/(s[1]*s[2]*s[3]);
						buffer.resize(buffSize);
						strip[0].Recv(&buffer[0],sizeof(float)*buffSize,curr);
						for (k=loc[4];k<=loc[5];k+=s[3])
							for (j=loc[2];j<=loc[3];j+=s[2])
								for (i=loc[0];i<=loc[1];i+=s[1])
								{
									i0 = (i - loc[0])/s[1];
									i0 += (loc[1]-loc[0]+s[1])*(j - loc[2])/(s[1]*s[2]);
									i0 += (loc[1]-loc[0]+s[1])*(loc[3]-loc[2]+s[2])*(k - loc[4])/(s[1]*s[2]*s[3]);
									i1 = (coords[1]*localCells[1]+i - glb[0])/s[1];
									i1 += pts[1]*(coords[2]*localCells[2]+j - glb[2])/s[2];
									i1 += pts[1]*pts[2]*(coords[3]*localCells[3]+k - glb[4])/s[3];
									gData[i1] = buffer[i0];
								}
					}
				}
			}

	if (thisNode==master)
	{
		outFile.open(fileName.str().c_str(),std::ios::app | std::ios::binary);
		WriteBigEndian((char*)&gData[0],sizeof(float)*pts[1]*pts[2]*pts[3],sizeof(float),outFile);
		outFile.close();
	}
}

#endif

void Grid::WriteMomentumDataHeader(const std::string& quantity,GridDataDescriptor* theBox)
{
	if (appendMode && restarted)
		return;

	tw::Int x0,x1,y0,y1,z0,z1;
	theBox->theRgn->GetGlobalCellBounds(&x0,&x1,&y0,&y1,&z0,&z1);

	// Only the task at index (0,0,0) will write data

	if (strip[0].Get_rank()!=0)
		return;

	std::ofstream outFile;
	std::stringstream fileName;
	tw::Int xPoints,yPoints,zPoints;

	fileName.str("");
	if (theBox->filename=="full")
		fileName << quantity;
	else
		fileName << theBox->filename << "_" << quantity;
	fileName << ".dvdat";

	if (dim[1]==1) theBox->xSkip = 1;
	if (dim[2]==1) theBox->ySkip = 1;
	if (dim[3]==1) theBox->zSkip = 1;

	xPoints = (x1 - x0 + 1)/theBox->xSkip;
	yPoints = (y1 - y0 + 1)/theBox->ySkip;
	zPoints = (z1 - z0 + 1)/theBox->zSkip;

	outFile.open(fileName.str().c_str(),std::ios::binary);
	WriteDVHeader(outFile,2,xPoints,yPoints,zPoints,
		-pi/spacing.x-pi/globalSize.x + 2.0*pi*(x0-1)/globalSize.x,
		-pi/spacing.x-pi/globalSize.x + 2.0*pi*x1/globalSize.x,
		-pi/spacing.y-pi/globalSize.y + 2.0*pi*(y0-1)/globalSize.y,
		-pi/spacing.y-pi/globalSize.y + 2.0*pi*y1/globalSize.y,
		-pi/spacing.z-pi/globalSize.z + 2.0*pi*(z0-1)/globalSize.z,
		-pi/spacing.z-pi/globalSize.z + 2.0*pi*z1/globalSize.z);
	outFile.close();
}

void Grid::WriteMomentumData(const std::string& quantity,GridDataDescriptor* theBox,tw::Float* theData,const tw::Int *stride)
{
	WriteBoxData(quantity,theBox,theData,stride);
}


void Grid::WriteCellDataHeader(GridDataDescriptor* theBox)
{
/*	if (appendMode && restarted)
		return;

	if (global->Index(local)!=0)
		return;

	std::ofstream gridStream;
	std::stringstream fileName;

	fileName.str("");
	if (theBox->filename=="full")
		fileName << "gridx";
	else
		fileName << theBox->filename << "_gridx";
	fileName << ".dvdat";

	gridStream.open(fileName.str().c_str(),std::ios::binary);
	WriteDVHeader(gridStream,2,global->xDim,1,1,0.5,float(global->xDim)+0.5,0.0,1.0,0.0,1.0);
	gridStream.close();

	fileName.str("");
	if (theBox->filename=="full")
		fileName << "gridy";
	else
		fileName << theBox->filename << "_gridy";
	fileName << ".dvdat";

	gridStream.open(fileName.str().c_str(),std::ios::binary);
	WriteDVHeader(gridStream,2,global->yDim,1,1,0.5,float(global->yDim)+0.5,0.0,1.0,0.0,1.0);
	gridStream.close();

	fileName.str("");
	if (theBox->filename=="full")
		fileName << "gridz";
	else
		fileName << theBox->filename << "_gridz";
	fileName << ".dvdat";

	gridStream.open(fileName.str().c_str(),std::ios::binary);
	WriteDVHeader(gridStream,2,global->zDim,1,1,0.5,float(global->zDim)+0.5,0.0,1.0,0.0,1.0);
	gridStream.close();*/
}

void Grid::WriteCellData(GridDataDescriptor* theBox)
{
/*	LocalDomain *master;
	master = (*global)(0,0,0);

	tw::Int i;
	std::ofstream gridStream;
	std::stringstream fileName;
	std::valarray<float> buffer;

	// X GRID
	// ------

	buffer.resize(global->xDim);
	fileName.str("");
	if (theBox->filename=="full")
		fileName << "gridx";
	else
		fileName << theBox->filename << "_gridx";
	fileName << ".dvdat";

	if (local==master)
	{
		gridStream.open(fileName.str().c_str(),std::ios::app | std::ios::binary);
		for (i=0;i<xDim;i++)
			buffer[i] = xpos[i+1];
		for (i=1;i<global->xDomains;i++)
			ReceiveAndWait(&buffer[i*xDim],sizeof(float)*xDim,(*global)(i,0,0)->thread);
		WriteBigEndian((char *)&buffer[0],sizeof(float)*global->xDim,sizeof(float),gridStream);
		gridStream.close();
	}
	else
	{
		if (local->y==0 && local->z==0)
		{
			for (i=0;i<xDim;i++)
				buffer[i] = xpos[i+1];
			SendAndWait(&buffer[0],sizeof(float)*xDim,master->thread);
		}
	}

	// Y GRID
	// ------

	buffer.resize(global->yDim);
	fileName.str("");
	if (theBox->filename=="full")
		fileName << "gridy";
	else
		fileName << theBox->filename << "_gridy";
	fileName << ".dvdat";

	if (local==master)
	{
		gridStream.open(fileName.str().c_str(),std::ios::app | std::ios::binary);
		for (i=0;i<yDim;i++)
			buffer[i] = ypos[i+1];
		for (i=1;i<global->yDomains;i++)
			ReceiveAndWait(&buffer[i*yDim],sizeof(float)*yDim,(*global)(0,i,0)->thread);
		WriteBigEndian((char *)&buffer[0],sizeof(float)*global->yDim,sizeof(float),gridStream);
		gridStream.close();
	}
	else
	{
		if (local->x==0 && local->z==0)
		{
			for (i=0;i<yDim;i++)
				buffer[i] = ypos[i+1];
			SendAndWait(&buffer[0],sizeof(float)*yDim,master->thread);
		}
	}

	// Z GRID
	// ------

	buffer.resize(global->zDim);
	fileName.str("");
	if (theBox->filename=="full")
		fileName << "gridz";
	else
		fileName << theBox->filename << "_gridz";
	fileName << ".dvdat";

	if (local==master)
	{
		gridStream.open(fileName.str().c_str(),std::ios::app | std::ios::binary);
		for (i=0;i<zDim;i++)
			buffer[i] = zpos[i+1];
		for (i=1;i<global->zDomains;i++)
			ReceiveAndWait(&buffer[i*zDim],sizeof(float)*zDim,(*global)(0,0,i)->thread);
		WriteBigEndian((char *)&buffer[0],sizeof(float)*global->zDim,sizeof(float),gridStream);
		gridStream.close();
	}
	else
	{
		if (local->x==0 && local->y==0)
		{
			for (i=0;i<zDim;i++)
				buffer[i] = zpos[i+1];
			SendAndWait(&buffer[0],sizeof(float)*zDim,master->thread);
		}
	}*/
}

void Grid::Diagnose()
{
	tw::Int i,s;
	std::stringstream fileName;
	std::ofstream restartFile,energyFile,pointFile;
	std::vector<tw::Float> cols;
	std::vector<bool> avg;
	std::valarray<tw::Float> energyBuffer;
	tw::Int master = 0;
	tw::Int curr = strip[0].Get_rank();

	// DIAGNOSTIC PREP
	// need something to handle custom

	bool doing_restart=false,doing_energy=false,doing_box=false,doing_point=false;

	doing_restart = dumpPeriod>0 && stepNow%dumpPeriod==0;
	for (s=0;s<energyDiagnostic.size();s++)
		doing_energy = doing_energy || energyDiagnostic[s]->WriteThisStep(elapsedTime,dt,stepNow);
	for (s=0;s<pointDiagnostic.size();s++)
		doing_point = doing_point || pointDiagnostic[s]->WriteThisStep(elapsedTime,dt,stepNow);
	for (s=0;s<boxDiagnostic.size();s++)
		doing_box = doing_box || boxDiagnostic[s]->WriteThisStep(elapsedTime,dt,stepNow);
	if (doing_restart || doing_energy || doing_point || doing_box)
		for (i=0;i<module.size();i++)
			module[i]->StartDiagnostics();

	// RESTART MECHANISM

	Lock();
	if (doing_restart)
	{
		fileName.str("");
		fileName << strip[0].Get_rank() << "_dump" << stepNow/dumpPeriod;
		restartFile.open(fileName.str().c_str(),std::ios::binary);
		WriteData(restartFile);
		restartFile.close();
	}
	Unlock();

	// ENERGY DIAGNOSTIC

	for (s=0;s<energyDiagnostic.size();s++)
	{
		fileName.str("");
		fileName << energyDiagnostic[s]->filename << ".txt";
		if (IsFirstStep() && !(appendMode&&restarted) && curr==master)
		{
			energyFile.open(fileName.str().c_str());
			energyFile << "Time dt ";
			for (i=0;i<module.size();i++)
				module[i]->EnergyHeadings(energyFile);
			energyFile << std::endl;
			energyFile.close();
		}
		if (energyDiagnostic[s]->WriteThisStep(elapsedTime,dt,stepNow))
		{
			if (curr==master)
			{
				energyFile.open(fileName.str().c_str(),std::ios::app);
				energyFile << elapsedTime << " " << dt << " ";
			}
			cols.clear();
			avg.clear();
			for (i=0;i<module.size();i++)
				module[i]->EnergyColumns(cols,avg,*energyDiagnostic[s]->theRgn);
			if (cols.size())
			{
				energyBuffer.resize(cols.size());
				for (i=0;i<cols.size();i++)
					energyBuffer[i] = cols[i];
				strip[0].Sum(&energyBuffer[0],&energyBuffer[0],cols.size()*sizeof(tw::Float),master);
				if (curr==master)
				{
					energyFile.precision(energyDiagnostic[s]->numSigFigs);
					for (i=0;i<energyBuffer.size();i++)
						energyFile << energyBuffer[i]/(avg[i]?strip[0].Get_size():1.0) << " ";
				}
			}
			if (curr==master)
			{
				energyFile << std::endl;
				energyFile.close();
			}
		}
	}

	// POINT DIAGNOSTIC

	Lock();
	tw::vec3 r;
	weights_3D w;
	for (s=0;s<pointDiagnostic.size();s++)
	{
		fileName.str("");
		fileName << pointDiagnostic[s]->filename << ".txt";
		if (IsFirstStep() && !(appendMode&&restarted) && curr==master)
		{
			pointFile.open(fileName.str().c_str());
			pointFile << "Time ";
			for (i=0;i<module.size();i++)
				module[i]->PointDiagnosticHeader(pointFile);
			pointFile << std::endl;
			pointFile.close();
		}
		r = pointDiagnostic[s]->thePoint;
		if (pointDiagnostic[s]->moveWithWindow && movingWindow)
			r.z += signalPosition;
		if (IsPointValid(r)) // assumes uniform grid
		{
			GetWeights(&w,r);
			if (pointDiagnostic[s]->WriteThisStep(elapsedTime,dt,stepNow))
			{
				pointFile.open(fileName.str().c_str(),std::ios::app);
				pointFile << elapsedTime << " ";
				for (i=0;i<module.size();i++)
					module[i]->PointDiagnose(pointFile,w);
				pointFile << std::endl;
				pointFile.close();
			}
		}
	}
	Unlock();

	// BOX DIAGNOSTIC

	for (s=0;s<boxDiagnostic.size();s++)
	{
		if (IsFirstStep())
		{
			WriteCellDataHeader(boxDiagnostic[s]);
			for (i=0;i<module.size();i++)
				module[i]->BoxDiagnosticHeader(boxDiagnostic[s]);
		}
		if (boxDiagnostic[s]->WriteThisStep(elapsedTime,dt,stepNow))
		{
			WriteCellData(boxDiagnostic[s]);
			for (i=0;i<module.size();i++)
				module[i]->BoxDiagnose(boxDiagnostic[s]);
		}
	}

	// CUSTOM DIAGNOSTICS

	for (i=0;i<module.size();i++)
		module[i]->CustomDiagnose();
}

// Alternate box diagnostic routines to write out all data with ghost cells

#ifdef ALT_BOX_DIAG

void Grid::WriteBoxDataHeader(const std::string& quantity,GridDataDescriptor* theBox)
{
	if (appendMode && restarted)
		return;

	Lock();
	std::ofstream outFile;
	std::stringstream fileName;
	tw::Int xPoints,yPoints,zPoints;

	fileName.str("");
	if (theBox->filename=="full")
		fileName << strip[0].Get_rank() << "_" << quantity;
	else
		fileName << strip[0].Get_rank() << "_" << theBox->filename << "_" << quantity;
	fileName << ".dvdat";

	xPoints = ub[1] - lb[1] + 1;
	yPoints = ub[2] - lb[2] + 1;
	zPoints = ub[3] - lb[3] + 1;

	outFile.open(fileName.str().c_str(),std::ios::binary);
	WriteDVHeader(outFile,2,xPoints,yPoints,zPoints,
		corner.x + spacing.x*(lb[1]-1),
		corner.x + spacing.x*ub[1],
		corner.y + spacing.y*(lb[2]-1),
		corner.y + spacing.y*ub[2],
		corner.z + spacing.z*(lb[3]-1),
		corner.z + spacing.z*ub[3]);
	outFile.close();
	Unlock();
}

void Grid::WriteBoxData(const std::string& quantity,GridDataDescriptor* theBox,tw::Float* theData,const tw::Int *stride)
{
	Lock();

	tw::Int i,j,k;
	std::ofstream outFile;
	std::stringstream fileName;
	float data;

	fileName.str("");
	if (theBox->filename=="full")
		fileName << strip[0].Get_rank() << "_" << quantity;
	else
		fileName << strip[0].Get_rank() << "_" << theBox->filename << "_" << quantity;
	fileName << ".dvdat";

	outFile.open(fileName.str().c_str(),std::ios::app | std::ios::binary);
	for (k=lb[3];k<=ub[3];k++)
		for (j=lb[2];j<=ub[2];j++)
			for (i=lb[1];i<=ub[1];i++)
			{
				data = theData[i*stride[1] + j*stride[2] + k*stride[3]];
				WriteBigEndian((char*)&data,sizeof(float),sizeof(float),outFile);
			}
	outFile.close();

	Unlock();
}

#endif

void Grid::EmergencyDump()
{
	tw::Int i,s;

	/*** Unconditionally write out all box diagnostics ***/

	for (s=0;s<boxDiagnostic.size();s++)
	{
		WriteCellData(boxDiagnostic[s]);
		for (i=0;i<module.size();i++)
			module[i]->BoxDiagnose(boxDiagnostic[s]);
	}
}
