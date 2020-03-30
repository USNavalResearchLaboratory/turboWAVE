#include "simulation.h"

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

//////////////////////
// DIAGNOSTIC TOOLS //
//////////////////////


Diagnostic::Diagnostic(const std::string& name,MetricSpace *ms,Task *tsk) : ComputeTool(name,ms,tsk)
{
	skip[0] = 0;
	skip[1] = skip[2] = skip[3] = 1;
	t = 0.0;
	t0 = 0.0;
	t1 = tw::big_pos;
	tRef = tw::big_neg;
	timePeriod = 0.0;
	filename = "diagnostic";
	vGalileo = 0.0;
	headerWritten = false;
	directives.Add("period",new tw::input::Int(&skip[0]),false);
	directives.Add("time period",new tw::input::Float(&timePeriod),false);
	directives.Add("skip",new tw::input::Numbers<tw::Int>(&skip[1],3),false);
	directives.Add("filename",new tw::input::String(&filename),false);
	directives.Add("t0",new tw::input::Float(&t0),false);
	directives.Add("t1",new tw::input::Float(&t1),false);
	directives.Add("galilean velocity",new tw::input::Vec3(&vGalileo),false);
}

bool Diagnostic::WriteThisStep(tw::Float elapsedTime,tw::Float dt,tw::Int stepNow)
{
	t = elapsedTime; // save for use in reports

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

	if (skip[0]!=0)
	{
		tw::Int startStep = tw::Int(t0/dt);
		if ( (elapsedTime>=t0) && (elapsedTime<=t1) && ((stepNow - startStep) % skip[0] == 0) )
			return true;
		else
			return false;
	}

	return false;
}

void Diagnostic::StartGridFile(std::ofstream& grid)
{
	std::string xname;
	if (filename=="tw::none")
		xname = "grid_warp.txt";
	else
		xname = filename + "_grid_warp.txt";
	if (headerWritten)
		grid.open(xname.c_str(),std::ios::app);
	else
		grid.open(xname.c_str());
	grid << "t = " << t << std::endl;
}

void Diagnostic::Start()
{
}

void Diagnostic::Finish()
{
}

void Diagnostic::Float(const std::string& label,tw::Float val,bool average)
{
}

tw::Float Diagnostic::VolumeIntegral(const std::string& fieldName,const struct Field& F,const tw::Int c)
{
	tw::Float ans = 0.0;
	tw::Int x0,x1,y0,y1,z0,z1;
	theRgn->GetLocalCellBounds(&x0,&x1,&y0,&y1,&z0,&z1);
	for (tw::Int k=z0;k<=z1;k++)
		for (tw::Int j=y0;j<=y1;j++)
			for (tw::Int i=x0;i<=x1;i++)
				if (theRgn->Inside(space->Pos(i,j,k),*space))
					ans += F(i,j,k,c) * space->dS(i,j,k,0);
	return ans;
}

tw::Float Diagnostic::FirstMoment(const std::string& fieldName,const struct Field& F,const tw::Int c,const tw::vec3& r0,const tw::grid::axis axis)
{
	tw::Float ans = 0.0;
	const tw::Int ax = tw::grid::naxis(axis);
	tw::Int x0,x1,y0,y1,z0,z1;
	theRgn->GetLocalCellBounds(&x0,&x1,&y0,&y1,&z0,&z1);
	for (tw::Int k=z0;k<=z1;k++)
		for (tw::Int j=y0;j<=y1;j++)
			for (tw::Int i=x0;i<=x1;i++)
			{
				const tw::vec3 pos = space->Pos(i,j,k);
				if (theRgn->Inside(pos,*space))
				{
					tw::vec3 r1 = r0;
					tw::vec3 r2 = pos;
					space->CurvilinearToCartesian(&r1);
					space->CurvilinearToCartesian(&r2);
					ans += F(i,j,k,c) * (r2[ax-1]-r1[ax-1]) * space->dS(i,j,k,0);
				}
			}
	return ans;
}

void Diagnostic::Field(const std::string& fieldName,const struct Field& F,const tw::Int c)
{
}

void Diagnostic::Particle(const struct Particle& par,tw::Float m0,tw::Float tp)
{
}

void Diagnostic::ReadCheckpoint(std::ifstream& inFile)
{
	ComputeTool::ReadCheckpoint(inFile);
	inFile.read((char *)&timePeriod,sizeof(tw::Float));
	inFile.read((char *)&t,sizeof(tw::Float));
	inFile.read((char *)&tRef,sizeof(tw::Float));
	inFile.read((char *)&t0,sizeof(tw::Float));
	inFile.read((char *)&t1,sizeof(tw::Float));
	inFile.read((char *)&headerWritten,sizeof(headerWritten));
}

void Diagnostic::WriteCheckpoint(std::ofstream& outFile)
{
	ComputeTool::WriteCheckpoint(outFile);
	outFile.write((char *)&timePeriod,sizeof(tw::Float));
	outFile.write((char *)&t,sizeof(tw::Float));
	outFile.write((char *)&tRef,sizeof(tw::Float));
	outFile.write((char *)&t0,sizeof(tw::Float));
	outFile.write((char *)&t1,sizeof(tw::Float));
	outFile.write((char *)&headerWritten,sizeof(headerWritten));
}

TextTableBase::TextTableBase(const std::string& name,MetricSpace *ms,Task *tsk) : Diagnostic(name,ms,tsk)
{
	numSigFigs = 6;
	filename = "table";
	directives.Add("precision",new tw::input::Int(&numSigFigs),false);
}

void TextTableBase::Start()
{
	Diagnostic::Start();
	labels.clear();
	values.clear();
	avg.clear();
}

void TextTableBase::Finish()
{
	Diagnostic::Finish();
	if (values.size())
	{
		const tw::Int master = 0;
		const tw::Int curr = task->strip[0].Get_rank();
		const tw::Int numRanks = task->strip[0].Get_size();
		std::valarray<tw::Float> buff;
		buff.resize(values.size());
		for (tw::Int i=0;i<values.size();i++)
			buff[i] = values[i]/(avg[i]?numRanks:1.0);
		task->strip[0].Sum(&buff[0],&buff[0],values.size()*sizeof(tw::Float),master);
		for (tw::Int i=0;i<values.size();i++)
			values[i] = buff[i];
		if (curr==master)
		{
			std::string fileName(filename+".txt");
			std::ofstream file;
			if (!headerWritten)
			{
				file.open(fileName.c_str());
				for (auto s : labels)
					file << s << " ";
				file << std::endl;
				headerWritten = true;
			}
			else
				file.open(fileName.c_str(),std::ios::app);
			file.precision(numSigFigs);
			for (auto v : values)
				file << v << " ";
			file << std::endl;
			file.close();
		}
	}
}

void TextTableBase::Float(const std::string& label,tw::Float val,bool average)
{
	Diagnostic::Float(label,val,average);
	labels.push_back(label);
	values.push_back(val);
	avg.push_back(average);
}

VolumeDiagnostic::VolumeDiagnostic(const std::string& name,MetricSpace *ms,Task *tsk) : TextTableBase(name,ms,tsk)
{
	filename = "energy";
}

tw::Float VolumeDiagnostic::VolumeIntegral(const std::string& fieldName,const struct Field& F,const tw::Int c)
{
	const tw::Float ans = Diagnostic::VolumeIntegral(fieldName,F,c);
	labels.push_back(fieldName);
	values.push_back(ans);
	avg.push_back(false);
	return ans;
}

tw::Float VolumeDiagnostic::FirstMoment(const std::string& fieldName,const struct Field& F,const tw::Int c,const tw::vec3& r0,const tw::grid::axis axis)
{
	const tw::Float ans = Diagnostic::FirstMoment(fieldName,F,c,r0,axis);
	labels.push_back(fieldName);
	values.push_back(ans);
	avg.push_back(false);
	return ans;
}

PointDiagnostic::PointDiagnostic(const std::string& name,MetricSpace *ms,Task *tsk) : TextTableBase(name,ms,tsk)
{
	thePoint = tw::vec3(1.0);
	filename = "pt1";
	directives.Add("point",new tw::input::Vec3(&thePoint));
}

void PointDiagnostic::Field(const std::string& fieldName,const struct Field& F,const tw::Int c)
{
	tw::vec3 r = thePoint + vGalileo*t;
	if (space->IsPointValid(r)) // assumes uniform grid
	{
		std::valarray<tw::Float> ans(1);
		weights_3D w;
		space->GetWeights(&w,r);
		F.Interpolate(ans,Element(c),w);
		labels.push_back(fieldName);
		values.push_back(ans[0]);
		avg.push_back(false);
	}
	else
	{
		labels.push_back(fieldName);
		values.push_back(0.0);
		avg.push_back(false);
	}
}

void PointDiagnostic::ReadCheckpoint(std::ifstream& inFile)
{
	TextTableBase::ReadCheckpoint(inFile);
	inFile.read((char *)&thePoint,sizeof(tw::vec3));
}

void PointDiagnostic::WriteCheckpoint(std::ofstream& outFile)
{
	TextTableBase::WriteCheckpoint(outFile);
	outFile.write((char *)&thePoint,sizeof(tw::vec3));
}

BoxDiagnostic::BoxDiagnostic(const std::string& name,MetricSpace *ms,Task *tsk) : Diagnostic(name,ms,tsk)
{
	average = false;
	filename = "tw::none";
	directives.Add("average",new tw::input::Bool(&average),false);
	directives.Add("reports",new tw::input::List<std::vector<std::string>>(&reports),false);
}

void BoxDiagnostic::GetGlobalIndexing(tw::Int pts[4],tw::Int glb[6])
{
	// On output, pts has the number of cells, glb has the global index bounds
	// skip is corrected for ignorable dimensions
	theRgn->GetGlobalCellBounds(glb);

	for (tw::Int ax=1;ax<=3;ax++)
	{
		const tw::Int lb = 2*ax-2;
		const tw::Int ub = 2*ax-1;
		skip[ax] = space->Dim(ax)==1 ? 1 : skip[ax];
		// Force global bounds into skipping sequence
		// The upper bound is adjusted to respect the lower bound
		pts[ax] = (glb[ub] - glb[lb] + 1)/skip[ax];
		glb[ub] = glb[lb] + skip[ax]*(pts[ax]-1);
	}
}

void BoxDiagnostic::GetLocalIndexing(const tw::Int pts[4],const tw::Int glb[6],tw::Int loc[6],const tw::Int coords[4])
{
	// On output, pts has the number of cells, glb has the global index bounds, loc has the local index bounds
	// coords is an input with the cartesian MPI domain indices (not necessarily the domain of execution)
	tw::Int i0;
	for (tw::Int ax=1;ax<=3;ax++)
	{
		const tw::Int cell0 = coords[ax]*task->localCells[ax];
		const tw::Int lb = 2*ax-2;
		const tw::Int ub = 2*ax-1;
		// Get the lower local cell bounds and force into global skipping sequence
		i0 = 1 + cell0;
		i0 = i0 < glb[lb] ? glb[lb] : i0;
		i0 -= cell0 + ((i0 - glb[lb]) % skip[ax]);
		if (i0<1) i0 += skip[ax];
		loc[lb] = i0;
		// Get the upper local cell bounds and force into global skipping sequence
		i0 = task->localCells[ax] + cell0;
		i0 = i0 > glb[ub] ? glb[ub] : i0;
		i0 -= cell0 + ((i0 - glb[lb]) % skip[ax]);
		loc[ub] = i0;
	}
}

void BoxDiagnostic::Field(const std::string& fieldName,const struct Field& F,const tw::Int c)
{
	if (reports.size()>0)
		if (std::find(reports.begin(),reports.end(),fieldName)==reports.end())
			return;

	tw::Int buffSize,ready,i0,i1;
	std::valarray<float> buffer,gData;
	std::ofstream outFile;
	std::string xname;
	tw::Int coords[4],pts[4],loc[6],glb[6],dim[4],s[4];
	tw::Int thisNode,curr,master;

	thisNode = task->strip[0].Get_rank();
	master = 0;

	if (filename=="tw::none")
		xname = fieldName + ".dvdat";
	else
		xname = filename + "_" + fieldName + ".dvdat";

	GetGlobalIndexing(pts,glb);
	for (tw::Int i=1;i<=3;i++) dim[i] = space->Dim(i);
	for (tw::Int i=1;i<=3;i++) s[i] = skip[i];

	if (thisNode==master)
		gData.resize(pts[1]*pts[2]*pts[3]);

	for (tw::Int kdom=(glb[4]-1)/dim[3];kdom<=(glb[5]-1)/dim[3];kdom++)
		for (tw::Int jdom=(glb[2]-1)/dim[2];jdom<=(glb[3]-1)/dim[2];jdom++)
			for (tw::Int idom=(glb[0]-1)/dim[1];idom<=(glb[1]-1)/dim[1];idom++)
			{
				curr = task->strip[0].Cart_rank(idom,jdom,kdom); // domain being written out
				coords[1] = idom; coords[2] = jdom; coords[3] = kdom;
				GetLocalIndexing(pts,glb,loc,coords);

				if (loc[0]<=dim[1] && loc[1]>=1 && loc[2]<=dim[2] && loc[3]>=1 && loc[4]<=dim[3] && loc[5]>=1)
				{
					if (curr==thisNode)
					{
						if (curr==master)
						{
							for (tw::Int k=loc[4];k<=loc[5];k+=s[3])
								for (tw::Int j=loc[2];j<=loc[3];j+=s[2])
									for (tw::Int i=loc[0];i<=loc[1];i+=s[1])
									{
										i0 = (task->GlobalCellIndex(i,1) - glb[0])/s[1];
										i0 += pts[1]*(task->GlobalCellIndex(j,2) - glb[2])/s[2];
										i0 += pts[1]*pts[2]*(task->GlobalCellIndex(k,3) - glb[4])/s[3];
										gData[i0] = F(i,j,k,c);
									}
						}
						else
						{
							task->strip[0].Recv(&ready,sizeof(tw::Int),master);
							buffSize = (loc[1]-loc[0]+s[1])*(loc[3]-loc[2]+s[2])*(loc[5]-loc[4]+s[3])/(s[1]*s[2]*s[3]);
							buffer.resize(buffSize);
							for (tw::Int k=loc[4];k<=loc[5];k+=s[3])
								for (tw::Int j=loc[2];j<=loc[3];j+=s[2])
									for (tw::Int i=loc[0];i<=loc[1];i+=s[1])
									{
										i0 = (i - loc[0])/s[1];
										i0 += (loc[1]-loc[0]+s[1])*(j - loc[2])/(s[1]*s[2]);
										i0 += (loc[1]-loc[0]+s[1])*(loc[3]-loc[2]+s[2])*(k - loc[4])/(s[1]*s[2]*s[3]);
										buffer[i0] = F(i,j,k,c);
									}
							task->strip[0].Send(&buffer[0],sizeof(float)*buffSize,master);
						}
					}
					if (thisNode==master && curr!=master)
					{
						task->strip[0].Send(&ready,sizeof(tw::Int),curr);
						buffSize = (loc[1]-loc[0]+s[1])*(loc[3]-loc[2]+s[2])*(loc[5]-loc[4]+s[3])/(s[1]*s[2]*s[3]);
						buffer.resize(buffSize);
						task->strip[0].Recv(&buffer[0],sizeof(float)*buffSize,curr);
						for (tw::Int k=loc[4];k<=loc[5];k+=s[3])
							for (tw::Int j=loc[2];j<=loc[3];j+=s[2])
								for (tw::Int i=loc[0];i<=loc[1];i+=s[1])
								{
									i0 = (i - loc[0])/s[1];
									i0 += (loc[1]-loc[0]+s[1])*(j - loc[2])/(s[1]*s[2]);
									i0 += (loc[1]-loc[0]+s[1])*(loc[3]-loc[2]+s[2])*(k - loc[4])/(s[1]*s[2]*s[3]);
									i1 = (coords[1]*task->localCells[1]+i - glb[0])/s[1];
									i1 += pts[1]*(coords[2]*task->localCells[2]+j - glb[2])/s[2];
									i1 += pts[1]*pts[2]*(coords[3]*task->localCells[3]+k - glb[4])/s[3];
									gData[i1] = buffer[i0];
								}
					}
				}
			}

	if (thisNode==master)
	{
		if (!headerWritten)
		{
			outFile.open(xname.c_str(),std::ios::binary);
			WriteDVHeader(outFile,2,pts[1],pts[2],pts[3],
			GlobalCorner(*space).x + dx(*space)*(glb[0]-0.5-0.5*skip[1]),
			GlobalCorner(*space).x + dx(*space)*(glb[1]-0.5+0.5*skip[1]),
			GlobalCorner(*space).y + dy(*space)*(glb[2]-0.5-0.5*skip[2]),
			GlobalCorner(*space).y + dy(*space)*(glb[3]-0.5+0.5*skip[2]),
			GlobalCorner(*space).z + dz(*space)*(glb[4]-0.5-0.5*skip[3]),
			GlobalCorner(*space).z + dz(*space)*(glb[5]-0.5+0.5*skip[3]));
			outFile.close();
		}
		outFile.open(xname.c_str(),std::ios::app | std::ios::binary);
		WriteBigEndian((char*)&gData[0],sizeof(float)*pts[1]*pts[2]*pts[3],sizeof(float),outFile);
		outFile.close();
	}
}

void BoxDiagnostic::Finish()
{
	// Write out the grid data

	std::ofstream gridFile;
	const tw::Int master = 0;
	const tw::Int curr_global = task->strip[0].Get_rank();
	tw::Int pts[4],glb[6];
	GetGlobalIndexing(pts,glb);
	if (curr_global==master)
		StartGridFile(gridFile);
	if (!headerWritten) // assuming static grid no need to keep writing spatial points
	{
		for (tw::Int ax=1;ax<=3;ax++)
		{
			// Message passing is needed to get the spatial points
			std::valarray<tw::Float> X(task->globalCells[ax]);
			const tw::Int offset = space->Dim(ax)*task->strip[ax].Get_rank();
			for (tw::Int i=1;i<=space->Dim(ax);i++)
				X[i-1+offset] = space->X(i,ax);
			task->strip[ax].Gather(&X[offset],&X[offset],task->localCells[ax]*sizeof(tw::Float),master);
			if (curr_global==master)
			{
				const tw::Int lb = glb[2*ax-2];
				const tw::Int ub = glb[2*ax-1];
				gridFile << "axis" << ax << " = ";
				for (tw::Int i=lb;i<ub;i+=skip[ax])
					gridFile << X[i-1] << " ";
				gridFile << X[ub-1] << std::endl;
			}
		}
	}
	if (curr_global==master)
		gridFile.close();
	headerWritten = true;
}

ParticleOrbits::ParticleOrbits(const std::string& name,MetricSpace *ms,Task *tsk) : Diagnostic(name,ms,tsk)
{
	filename = "par";
	minGamma = 1.0;
	directives.Add("minimum gamma",new tw::input::Float(&minGamma));
}

void ParticleOrbits::Start()
{
	Diagnostic::Start();
	parData.clear();
}

void ParticleOrbits::Finish()
{
	const tw::Int master = 0;
	std::valarray<float> parBuffer;
	parBuffer.resize(parData.size());
	for (tw::Int i=0;i<parData.size();i++)
		parBuffer[i] = parData[i];
	tw::Int pts = parData.size();
	if (task->strip[0].Get_rank()!=master)
	{
		task->strip[0].Send(&pts,sizeof(tw::Int),master);
		task->strip[0].Send(&parBuffer[0],sizeof(float)*pts,master);
	}
	if (task->strip[0].Get_rank()==master)
	{
		std::string fileName = filename + ".dvpar";
		std::ofstream outFile;
		if (!headerWritten)
		{
			outFile.open(fileName.c_str(),std::ios::binary);
			headerWritten = true;
		}
		else
			outFile.open(fileName.c_str(),std::ios::binary | std::ios::app);
		WriteBigEndian((char*)&parBuffer[0],sizeof(float)*pts,sizeof(float),outFile);
		for (tw::Int i=0;i<task->strip[0].Get_size();i++)
		{
			if (i!=master)
			{
				task->strip[0].Recv(&pts,sizeof(tw::Int),i);
				parBuffer.resize(pts);
				task->strip[0].Recv(&parBuffer[0],sizeof(float)*pts,i);
				WriteBigEndian((char*)&parBuffer[0],sizeof(float)*pts,sizeof(float),outFile);
			}
		}
		float zero = 0.0;
		for (tw::Int i=0;i<8;i++)
			WriteBigEndian((char*)&zero,sizeof(float),sizeof(float),outFile);
		outFile.close();
	}
}

void ParticleOrbits::Particle(const struct Particle& par,tw::Float m0,tw::Float tp)
{
	tw::vec3 position = space->PositionFromPrimitive(par.q);
	if (theRgn->Inside(position,*space))
		if (sqrt(1.0 + Norm(par.p/m0)) >= minGamma)
		{
			space->CurvilinearToCartesian(&position);
			parData.push_back(position.x - vGalileo.x*tp);
			parData.push_back(par.p.x);
			parData.push_back(position.y - vGalileo.y*tp);
			parData.push_back(par.p.y);
			parData.push_back(position.z - vGalileo.z*tp);
			parData.push_back(par.p.z);
			parData.push_back(par.aux1);
			parData.push_back(par.aux2);
		}
}

void ParticleOrbits::ReadCheckpoint(std::ifstream& inFile)
{
	Diagnostic::ReadCheckpoint(inFile);
	inFile.read((char *)&minGamma,sizeof(tw::Float));
}

void ParticleOrbits::WriteCheckpoint(std::ofstream& outFile)
{
	Diagnostic::WriteCheckpoint(outFile);
	outFile.write((char *)&minGamma,sizeof(tw::Float));
}

PhaseSpaceDiagnostic::PhaseSpaceDiagnostic(const std::string& name,MetricSpace *ms,Task *tsk) : Diagnostic(name,ms,tsk)
{
	filename = "xpx";
	ax[0] = tw::grid::t;
	ax[1] = tw::grid::x;
	ax[2] = tw::grid::px;
	ax[3] = tw::grid::py;
	bounds[0] = bounds[2] = bounds[4] = 0.0;
	bounds[1] = bounds[3] = bounds[5] = 1.0;
	dims[0] = 1;
	dims[1] = 100;
	dims[2] = 100;
	dims[3] = 1;
	directives.Add("axes",new tw::input::Enums<tw::grid::axis>(tw::grid::axis_map(),&ax[1],&ax[2],&ax[3]));
	directives.Add("dimensions",new tw::input::Numbers<tw::Int>(&dims[1],3));
	directives.Add("bounds",new tw::input::Numbers<tw::Float>(&bounds[0],6));
}

void PhaseSpaceDiagnostic::Start()
{
	Diagnostic::Start();

	if (!headerWritten && task->strip[0].Get_rank()==0)
	{
		std::string fileName = filename + ".dvdat";
		std::ofstream outFile;
		outFile.open(fileName.c_str(),std::ios::binary);
		WriteDVHeader(outFile,2,dims[1],dims[2],dims[3],bounds[0],bounds[1],bounds[2],bounds[3],bounds[4],bounds[5]);
		outFile.close();
		// don't set headerWritten until end of Finish() due to grid file
	}

	tw::vec3 phaseSpaceMin(bounds[0],bounds[2],bounds[4]);
	tw::vec3 phaseSpaceSize(bounds[1]-bounds[0],bounds[3]-bounds[2],bounds[5]-bounds[4]);
	plot_layout.Resize(dims[1],dims[2],dims[3],phaseSpaceMin,phaseSpaceSize,1);
	fxp.Initialize(plot_layout,task);
	fxp.SetBoundaryConditions(tw::grid::x,tw::bc::fld::dirichletCell,tw::bc::fld::dirichletCell);
	fxp.SetBoundaryConditions(tw::grid::y,tw::bc::fld::dirichletCell,tw::bc::fld::dirichletCell);
	fxp.SetBoundaryConditions(tw::grid::z,tw::bc::fld::dirichletCell,tw::bc::fld::dirichletCell);
}

void PhaseSpaceDiagnostic::Finish()
{
	tw::Int master = 0;
	task->strip[0].Sum(&fxp(0,0,0),&fxp(0,0,0),fxp.TotalBytes(),master);
	if (task->strip[0].Get_rank()==master)
	{
		fxp.ApplyBoundaryCondition();
		std::string fileName = filename + ".dvdat";
		std::ofstream outFile;
		outFile.open(fileName.c_str(),std::ios::binary | std::ios::app);
		// Convention is to write in FORTRAN order, so need an explicit nested loop
		for (tw::Int k=1;k<=dims[3];k++)
			for (tw::Int j=1;j<=dims[2];j++)
				for (tw::Int i=1;i<=dims[1];i++)
				{
					float data = fxp(i,j,k);
					WriteBigEndian((char *)&data,sizeof(float),0,outFile);
		}
		outFile.close();
	}

	// Write out the phase space grid data

	std::ofstream gridFile;
	if (task->strip[0].Get_rank()==0)
	{
		StartGridFile(gridFile);
		if (!headerWritten) // assuming static grid no need to keep writing spatial points
		{
			for (tw::Int ax=1;ax<=3;ax++)
			{
				// Unlike box diagnostic the grid is always uniform, no message passing needed
				gridFile << "axis" << ax << " = ";
				const tw::Float ds = (bounds[2*ax-1] - bounds[2*ax-2])/tw::Float(dims[ax]);
				const tw::Float s1 = bounds[2*ax-2] + 0.5*ds;
				const tw::Float sN = bounds[2*ax-1] - 0.5*ds;
				for (tw::Int i=0;i<dims[ax]-1;i++)
					gridFile << s1 + ds*tw::Float(i) << " ";
				gridFile << sN << std::endl;
			}
		}
		gridFile.close();
	}
	headerWritten = true;
}

void PhaseSpaceDiagnostic::Particle(const struct Particle& par,tw::Float m0,tw::Float tp)
{
	weights_3D weights;
	tw::vec3 r = space->PositionFromPrimitive(par.q);
	const tw::vec3 p = par.p;
	const tw::Float dV = (bounds[1]-bounds[0])*(bounds[3]-bounds[2])*(bounds[5]-bounds[4])/(dims[1]*dims[2]*dims[3]);
	if (theRgn->Inside(r,*space))
	{
		space->CurvilinearToCartesian(&r);
		r -= vGalileo*tp;
		std::map<tw::grid::axis,tw::Float> m = {{tw::grid::t,tp},{tw::grid::x,r.x},{tw::grid::y,r.y},{tw::grid::z,r.z},
			{tw::grid::mass,sqrt(m0*m0+Norm(p))},{tw::grid::px,p.x},{tw::grid::py,p.y},{tw::grid::pz,p.z},
			{tw::grid::g,sqrt(1+Norm(p)/(m0*m0))},{tw::grid::gbx,p.x/m0},{tw::grid::gby,p.y/m0},{tw::grid::gbz,p.z/m0}};
		tw::vec3 q(m[ax[1]],m[ax[2]],m[ax[3]]);
		if (dims[1]==1) q.x = 0.0;
		if (dims[2]==1) q.y = 0.0;
		if (dims[3]==1) q.z = 0.0;
		if (q.x>bounds[0] && q.x<bounds[1] && q.y>bounds[2] && q.y<bounds[3] && q.z>bounds[4] && q.z<bounds[5])
		{
			fxp.GetWeights(&weights,q);
			fxp.InterpolateOnto( par.number/dV, weights );
		}
	}
}
