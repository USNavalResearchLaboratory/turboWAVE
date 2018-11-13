#include "sim.h"
#include "particles.h"
#include "fieldSolve.h"
#include "electrostatic.h"
#include "laserSolve.h"
#include "fluid.h"
#include "quantum.h"
#include "solidState.h"
//#include <unistd.h>

//////////////////////
//  TEXT PROCESSOR  //
//////////////////////


void ReduceInputFile(std::ifstream& inFile,std::stringstream& out)
{
	// Strips comments
	// Strips decorative punctuation tokens (commas,parenthesis)
	// Inserts spaces before and after equals signs and curly braces
	bool ignoreUntilMark,ignoreUntilEOL;
	char charNow,nextChar;

	ignoreUntilMark = false;
	ignoreUntilEOL = false;
	while (inFile.get(charNow))
	{
		if (charNow=='/' || charNow=='*' || charNow=='\n' || charNow=='\r')
		{
			if (charNow=='/')
			{
				inFile.get(nextChar);
				if (nextChar=='/')
					ignoreUntilEOL = true;
				if (nextChar=='*' && !ignoreUntilEOL)
					ignoreUntilMark = true;
				if (nextChar!='/' && nextChar!='*')
				{
					if (!ignoreUntilMark && !ignoreUntilEOL)
						out << charNow;
					inFile.putback(nextChar);
				}
			}
			if (charNow=='*')
			{
				inFile.get(nextChar);
				if (nextChar=='/')
					ignoreUntilMark = false;
				else
				{
					if (!ignoreUntilMark && !ignoreUntilEOL)
						out << charNow;
					inFile.putback(nextChar);
				}
			}
			if (charNow=='\n' || charNow=='\r')
			{
				if (!ignoreUntilMark)
					out << charNow;
				ignoreUntilEOL = false;
			}
		}
		else
		{
			if (!ignoreUntilMark && !ignoreUntilEOL)
			{
				if (charNow==',' || charNow=='(' || charNow==')')
					charNow = ' ';
				if (charNow=='{' || charNow=='}' || charNow=='=')
					out << ' ';
				out << charNow;
				if (charNow=='{' || charNow=='}' || charNow=='=')
					out << ' ';
			}
		}
	}
}

void PreprocessInputFile(std::ifstream& inFile,std::stringstream& out)
{
	// This routine erroneously adds closing braces to the end of "out".
	// This does not corrupt subsequent parsing, however.
	std::ifstream *includedFile;
	std::stringstream stage1,stage2,stage3;
	std::string word;
	tw::Float unitDensityCGS;

	// Handle included files, strip comments, clean whitespace
	ReduceInputFile(inFile,stage1);
	stage1.seekg(0);
	while (!stage1.eof())
	{
		stage1 >> word;
		if (word=="include")
		{
			stage1 >> word;
			includedFile = new std::ifstream(word.c_str());
			if (includedFile->rdstate() & std::ios::failbit)
				throw tw::FatalError("couldn't open " + word);
			ReduceInputFile(*includedFile,stage2);
			delete includedFile;
		}
		else
		{
			stage2 << word << " ";
		}
	}

	// Set up the unit converter
	stage2.seekg(0);
	while (!stage2.eof())
	{
		stage2 >> word;
		stage3 << word << " ";
		if (word=="unit")
		{
			stage2 >> word;
			stage3 << word << " ";
			if (word=="density")
			{
				stage2 >> word >> unitDensityCGS;
				stage3 << " = " << unitDensityCGS << " ";
			}
		}
	}
	UnitConverter uc(unitDensityCGS);

	// Unit conversion macro substitution
	stage3.seekg(0);
	while (!stage3.eof())
	{
		stage3 >> word;
		if (word[0]=='%')
			tw::input::NormalizeInput(uc,word);
		out << word << " ";
	}
}




///////////////////////////
//  NON-UNIFORM REGIONS  //
///////////////////////////


NonUniformRegion::NonUniformRegion(tw::Int first,tw::Int last,tw::Float length,tw::Float dz0)
{
	i1 = first;
	i2 = last;
	ih = (i1+i2)/2;
	L = length;
	dz = dz0;
	N = i2 - i1 + 1;

	tw::Int i;
	gridSum = 0.0;
	for (i=1;i<=N;i++)
		gridSum += QuinticPulse(tw::Float(i-1)/tw::Float(N-1));
}

tw::Float NonUniformRegion::AddedCellWidth(tw::Int globalCell)
{
	tw::Float A = (1.0/gridSum)*(L/dz - N);
	if (globalCell>=i1 && globalCell<=i2)
		return dz*A*QuinticPulse(tw::Float(globalCell-i1)/tw::Float(N-1));
	else
		return 0.0;
}

tw::Float NonUniformRegion::ACoefficient(tw::Float length)
{
	return (1.0/gridSum)*(length/dz - N);
}

void Grid::SetCellWidthsAndLocalSize()
{
	tw::Int i,j;

	// Overwrite uniform cell widths with requested variable widths.
	// Correct the local domain size to reflect the variable cell widths.

	if (radialProgressionFactor!=1.0)
	{
		for (i=cornerCell[1]-1;i<=cornerCell[1]+dim[1];i++)
		{
			if (i > globalCells[1]/3)
				dX(i-cornerCell[1]+1,1) = spacing.x*pow(radialProgressionFactor,tw::Float(i-globalCells[1]/3));
		}
	}

	if (dim[3]>1 && region.size())
	{
		for (i=lb[3];i<=ub[3];i++)
		{
			dX(i,3) = spacing.z;
			for (j=0;j<region.size();j++)
				dX(i,3) += region[j]->AddedCellWidth(cornerCell[3] - 1 + i);
		}
	}

	size = 0.0;
	for (i=1;i<=dim[1];i++)
		size.x += dX(i,1);
	for (i=1;i<=dim[2];i++)
		size.y += dX(i,2);
	for (i=1;i<=dim[3];i++)
		size.z += dX(i,3);
}

void Grid::SetGlobalSizeAndLocalCorner()
{
	// Perform message passing to determine the effect of non-uniform cell widths
	// on the global domain size and the coordinates of the local domain corner.
	// Assumes local domain sizes are already calculated.

	tw::Int axis,src,dst;
	tw::Float inData,outData;
	tw::Float lsize[4] = { 0.0 , size.x , size.y , size.z };
	tw::Float gcorn[4] = { 0.0 , globalCorner.x , globalCorner.y , globalCorner.z };
	tw::Float lcorn[4];
	tw::Float gsize[4];

	for (axis=1;axis<=3;axis++)
	{
		finiteStrip[axis].Shift(1,1,&src,&dst);
		if (domainIndex[axis]==0)
		{
			finiteStrip[axis].Send(&lsize[axis],sizeof(tw::Float),dst);
			lcorn[axis] = gcorn[axis];
			outData = lsize[axis]; // in case domains=1
		}
		else
		{
			finiteStrip[axis].Recv(&inData,sizeof(tw::Float),src);
			lcorn[axis] = gcorn[axis] + inData;
			outData = inData + lsize[axis];
			if (domainIndex[axis]!=domains[axis]-1)
				finiteStrip[axis].Send(&outData,sizeof(tw::Float),dst);
		}
		finiteStrip[axis].Shift(1,-1,&src,&dst);
		if (domainIndex[axis]==domains[axis]-1)
		{
			gsize[axis] = outData;
			finiteStrip[axis].Send(&gsize[axis],sizeof(tw::Float),dst);
		}
		else
		{
			finiteStrip[axis].Recv(&inData,sizeof(tw::Float),src);
			gsize[axis] = inData;
			if (domainIndex[axis]!=0)
				finiteStrip[axis].Send(&gsize[axis],sizeof(tw::Float),dst);
		}
		((tw::Float*)&corner)[axis-1] = lcorn[axis];
		((tw::Float*)&globalSize)[axis-1] = gsize[axis];
	}
}



////////////
//  GRID  //
////////////


Grid::Grid()
{
	clippingRegion.push_back(new EntireRegion(clippingRegion));

	gridGeometry = cartesian;

	dt0 = 0.1;
	dt = 0.1;
	dth = 0.05;
	dtMin = tw::small_pos;
	dtMax = tw::big_pos;
	elapsedTime = 0.0;
	elapsedTimeMax = tw::big_pos;
	signalPosition = 0.0;
	windowPosition = 0.0;
	signalSpeed = 1.0;
	antiSignalPosition = 0.0;
	antiWindowPosition = 0.0;

	radialProgressionFactor = 1.0;

	appendMode = true;
	fullOutput = false;
	neutralize = true;
	smoothing = 0;
	compensation = 0;
	movingWindow = false;
	restarted = false;
	completed = false;
	adaptiveTimestep = false;
	adaptiveGrid = false;

	stepNow = 1;
	stepsToTake = 32;
	lastTime = 0;
	binaryFormat = 3;

	dumpPeriod = 0;

	bc0[1] = cyclic;
	bc1[1] = cyclic;
	bc0[2] = cyclic;
	bc1[2] = cyclic;
	bc0[3] = absorbing;
	bc1[3] = absorbing;

	#ifdef USE_OPENCL
	waveBuffer = NULL;
	#endif
}

Grid::~Grid()
{
	tw::Int i;

	for (i=0;i<wave.size();i++)
		delete wave[i];
	for (i=0;i<pulse.size();i++)
		delete pulse[i];
	for (i=0;i<energyDiagnostic.size();i++)
		delete energyDiagnostic[i];
	for (i=0;i<pointDiagnostic.size();i++)
		delete pointDiagnostic[i];
	for (i=0;i<boxDiagnostic.size();i++)
		delete boxDiagnostic[i];
	for (i=0;i<conductor.size();i++)
		delete conductor[i];
	for (i=0;i<clippingRegion.size();i++)
		delete clippingRegion[i];
	for (i=0;i<region.size();i++)
		delete region[i];

	for (i=0;i<computeTool.size();i++)
		delete computeTool[i];
	for (i=0;i<module.size();i++)
		delete module[i];

	if (uniformDeviate!=NULL)
		delete uniformDeviate;
	if (gaussianDeviate!=NULL)
		delete gaussianDeviate;

	if (dynamic_cast<std::ofstream*>(tw_out))
		((std::ofstream*)tw_out)->close();

	#ifdef USE_OPENCL
	if (waveBuffer!=NULL)
		clReleaseMemObject(waveBuffer);
	#endif
}

void Grid::Run()
{
	std::ofstream twstat;
	if (strip[0].Get_rank()==0)
	{
		twstat.open("twstat");
		twstat << "TurboWAVE is initializing.";
		twstat.close();
	}

	try
	{
		(*tw_out) << std::endl << "*** Prepare Simulation ***" << std::endl << std::endl;

		PrepareSimulation();

		(*tw_out) << std::endl << "*** Begin Simulation ***" << std::endl << std::endl;

		tw::Int startTime = GetSeconds();
		lastTime = startTime;

		if (GetSeconds()<0)
		{
			(*tw_out) << std::endl << "WARNING: System clock is not responding properly." << std::endl << std::endl;
		}

		(*tw_out) << "Current status can be viewed in 'twstat' file or by pressing enter key." << std::endl;
		(*tw_out) << "Enter 'help' for list of interactive commands." << std::endl << std::endl;

		while (stepNow <= stepsToTake && elapsedTime < elapsedTimeMax)
		{
			if ((GetSeconds() > lastTime + 5) && strip[0].Get_rank()==0)
			{
				twstat.open("twstat");
				InteractiveCommand("status",&twstat);
				twstat.close();
				lastTime = GetSeconds();
			}
			FundamentalCycle();
		}

		(*tw_out) << "Completed " << stepNow - 1 << " steps in " << GetSeconds() - startTime << " seconds." << std::endl;
		(*tw_out) << "Simulated elapsed time = " << elapsedTime << std::endl;
		if (strip[0].Get_rank()==0)
		{
			twstat.open("twstat");
			twstat << "Completed " << stepNow - 1 << " steps in " << GetSeconds() - startTime << " seconds." << std::endl;
			twstat << "Simulated elapsed time = " << elapsedTime << std::endl;
			twstat.close();
		}
	}
	catch (tw::FatalError& e)
	{
		(*tw_out) << "FATAL ERROR: " << e.what() << std::endl;
		(*tw_out) << "Simulation failed --- exiting now." << std::endl;
		#ifdef USE_TW_MPI
		if (tw_out != &std::cout)
		{
			std::cout << "FATAL ERROR: " << e.what() << std::endl;
			std::cout << "Simulation failed --- exiting now." << std::endl;
		}
		#endif
		if (strip[0].Get_rank()==0)
		{
			twstat.open("twstat");
			twstat << "The simulation failed. For more info see stdout." << std::endl;
			twstat.close();
		}
		completed = true;
		exit(1);
	}

	completed = true;
}

void Grid::SetupGeometry()
{
	// This routine assumes that MetricSpace::width, and MetricSpace::corner are valid
	switch (gridGeometry)
	{
		case cartesian:
			if (stepNow==1)
				(*tw_out) << "Using CARTESIAN Grid" << std::endl;
			SetCartesianGeometry();
			break;
		case cylindrical:
			if (stepNow==1)
				(*tw_out) << "Using CYLINDRICAL Grid" << std::endl;
			SetCylindricalGeometry();
			break;
		case spherical:
			if (stepNow==1)
				(*tw_out) << "Using SPHERICAL Grid" << std::endl;
			SetSphericalGeometry();
			break;
	}
}

void Grid::PrepareSimulation()
{
	std::ofstream twstat;
	tw::Int i;

	#ifdef USE_OPENCL
	PrintGPUInformation();
	#endif

	if (!restarted)
		GridFromInputFile();

	ReadInputFile();

	// Grid initialization done during reading of input file
	// because Module constructors are allowed to assume the Grid is fully specified

	// Initialize Regions

	if (!restarted)
		for (i=1;i<clippingRegion.size();i++)
			clippingRegion[i]->Initialize(*this,this);
	// region 0 is not saved in restart file
	clippingRegion[0]->Initialize(*this,this);

	// Initialize Injection Objects

	if (!restarted)
	{
		for (i=0;i<wave.size();i++)
			wave[i]->Initialize();
		for (i=0;i<pulse.size();i++)
			pulse[i]->Initialize();
		for (i=0;i<conductor.size();i++)
			conductor[i]->Initialize(*this);
	}
	#ifdef USE_OPENCL
	cl_int err;
	std::valarray<tw::Float> packed_waves(15*wave.size());
	for (i=0;i<wave.size();i++)
	{
		packed_waves[15*i] = wave[i]->a0;
		packed_waves[15*i+1] = wave[i]->w;
		for (tw::Int c=0;c<3;c++)
		{
			packed_waves[15*i+2+c] = wave[i]->laserFrame.u[c];
			packed_waves[15*i+5+c] = wave[i]->laserFrame.w[c];
			packed_waves[15*i+8+c] = wave[i]->focusPosition[c];
		}
		packed_waves[15*i+11] = wave[i]->pulseShape.t1;
		packed_waves[15*i+12] = wave[i]->pulseShape.t2;
		packed_waves[15*i+13] = wave[i]->pulseShape.t3;
		packed_waves[15*i+14] = wave[i]->pulseShape.t4;
	}
	waveBuffer = clCreateBuffer(context,CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,sizeof(tw::Float)*packed_waves.size(),&packed_waves[0],&err);
	#endif

	// Sort Modules

	ModuleComparator comparatorObject;
	std::sort(module.begin(),module.end(),comparatorObject);

	// Initialize Computational Tools
	// Must precede module initialization

	(*tw_out) << std::endl << "Initializing Compute Tools..." << std::endl << std::endl;

	for (i=0;i<computeTool.size();i++)
	{
		(*tw_out) << "Tool: " << computeTool[i]->name << std::endl;
		computeTool[i]->Initialize();
		computeTool[i]->WarningMessage(tw_out);
	}
	if (computeTool.size()==0)
		(*tw_out) << "(no tools)" << std::endl;

	// Initialize Modules

	(*tw_out) << std::endl << "Initializing Modules..." << std::endl << std::endl;

	for (i=0;i<module.size();i++)
		module[i]->ExchangeResources();

	for (i=0;i<module.size();i++)
	{
		(*tw_out) << "Module: " << module[i]->name << std::endl;
		module[i]->Initialize();
		module[i]->WarningMessage(tw_out);
	}
}

#ifdef USE_OPENCL

void Grid::PrintGPUInformation()
{
	cl_ulong ninfo;

	*tw_out << initMessage;

	*tw_out << "GPU INFORMATION" << std::endl;
	*tw_out << "--------------------------------------------" << std::endl;

	clGetDeviceInfo(gpu,CL_DEVICE_GLOBAL_MEM_SIZE,sizeof(ninfo),&ninfo,NULL);
	*tw_out << "Global memory: " << ninfo << std::endl;

	clGetDeviceInfo(gpu,CL_DEVICE_LOCAL_MEM_SIZE,sizeof(ninfo),&ninfo,NULL);
	*tw_out << "Local memory: " << ninfo << std::endl;

	clGetDeviceInfo(gpu,CL_DEVICE_MAX_WORK_GROUP_SIZE,sizeof(ninfo),&ninfo,NULL);
	*tw_out << "Maximum work group size: " << ninfo << std::endl;

	clGetDeviceInfo(gpu,CL_DEVICE_PREFERRED_VECTOR_WIDTH_FLOAT,sizeof(ninfo),&ninfo,NULL);
	*tw_out << "Float vector width: " << ninfo << std::endl;

	clGetDeviceInfo(gpu,CL_DEVICE_PREFERRED_VECTOR_WIDTH_DOUBLE,sizeof(ninfo),&ninfo,NULL);
	*tw_out << "Double vector width: " << ninfo << std::endl;

	*tw_out << "--------------------------------------------" << std::endl << std::endl;

}

#endif

void Grid::InteractiveCommand(const std::string& cmd,std::ostream *theStream)
{
	if (cmd=="help" || cmd=="?")
	{
		*theStream << "--- List of Interactive Commands ---" << std::endl;
		*theStream << "status or [enter] : print current step and other status indicators" << std::endl;
		*theStream << "metrics : print grid and time step metrics for this simulation" << std::endl;
		*theStream << "list : list modules and compute tools and their ID numbers" << std::endl;
		//*theStream << "peek [x] [y] [z] : print current data at cell x,y,z" << std::endl;
		*theStream << "Ctrl-C : abort the simulation" << std::endl;
		*theStream << std::endl;
	}
	if (cmd=="status" || cmd=="")
	{
		*theStream << "Current step: " << stepNow << std::endl;
		*theStream << "Current step size: " << dt << std::endl;
		*theStream << "Current elapsed time: " << elapsedTime << std::endl;
		for (tw::Int i=0;i<module.size();i++)
			module[i]->StatusMessage(theStream);
		for (tw::Int i=0;i<computeTool.size();i++)
			computeTool[i]->StatusMessage(theStream);
		*theStream << std::endl;
	}
	if (cmd=="list")
	{
		*theStream << "--- List of Modules ---" << std::endl;
		for (tw::Int i=0;i<module.size();i++)
			*theStream << i << " = " << module[i]->name << std::endl;
		*theStream << "--- List of Tools ---" << std::endl;
		for (tw::Int i=0;i<computeTool.size();i++)
			*theStream << i << " = " << computeTool[i]->name << std::endl;
		*theStream << std::endl;
	}
	if (cmd=="metrics")
	{
		*theStream << "Steps to take: " << stepsToTake << std::endl;
		*theStream << "Steps remaining: " << stepsToTake - stepNow << std::endl;
		*theStream << "Global grid size: " << globalCells[1] << "," << globalCells[2] << "," << globalCells[3] << std::endl;
		*theStream << "Local grid size: " << localCells[1] << "," << localCells[2] << "," << localCells[3] << std::endl;
		*theStream << "MPI Domains: " << domains[1] << "," << domains[2] << "," << domains[3] << std::endl;
		*theStream << std::endl;
	}
	if (cmd.find("peek")!=std::string::npos)
	{
		*theStream << "Not implemented yet.";
		*theStream << std::endl;
	}
}

void Grid::FundamentalCycle()
{
	tw::Int i;

	Diagnose();

	for (i=0;i<module.size();i++)
		module[i]->Reset();

	for (i=0;i<module.size();i++)
		module[i]->Update();

	elapsedTime += dt;
	signalPosition += signalSpeed*dt;
	antiSignalPosition -= signalSpeed*dt;
	stepNow++;

	if (adaptiveGrid)
		for (i=0;i<module.size();i++)
			module[i]->AdaptGrid();

	if (movingWindow && signalPosition>=(windowPosition + spacing.z) && dim[3]>1)
		MoveWindow();

	if (!movingWindow && antiSignalPosition<=(antiWindowPosition - spacing.z) && dim[3]>1)
		AntiMoveWindow();
}

bool Grid::MangleModuleName(std::string& name)
{
	bool trouble,did_mangle;
	tw::Int id = 2;
	std::string mangled(name);
	do
	{
		trouble = false;
		for (tw::Int i=0;i<module.size();i++)
			if (module[i]->name==mangled)
				trouble = true;
		if (trouble)
			mangled = name + std::to_string(id);
		id++;
	} while (trouble);
	did_mangle = (name!=mangled);
	name = mangled;
	return did_mangle;
}

Module* Grid::GetModule(const std::string& name)
{
	for (tw::Int i=0;i<module.size();i++)
		if (module[i]->name==name)
			return module[i];
	throw tw::FatalError("Could not find module: " + name);
	return NULL;
}

void Grid::MangleToolName(std::string& name)
{
	bool trouble,did_mangle;
	tw::Int id = 2;
	std::string mangled(name);
	do
	{
		trouble = false;
		for (tw::Int i=0;i<computeTool.size();i++)
			if (computeTool[i]->name==mangled)
				trouble = true;
		if (trouble)
			mangled = name + std::to_string(id);
		id++;
	} while (trouble);
	did_mangle = (name!=mangled);
	name = mangled;
	return did_mangle;
}

ComputeTool* Grid::CreateTool(const std::string& basename,tw::tool_type theType)
{
	std::string name(basename);
	MangleToolName(name);
	(*tw_out) << "Creating Tool " << name << "..." << std::endl;
	computeTool.push_back(ComputeTool::CreateObjectFromType(name,theType,this,this));
	computeTool.back()->refCount++;
	return computeTool.back();
}

ComputeTool* Grid::GetTool(const std::string& name)
{
	for (tw::Int i=0;i<computeTool.size();i++)
	{
		if (computeTool[i]->name==name)
		{
			computeTool[i]->refCount++;
			return computeTool[i];
		}
	}
	throw tw::FatalError("Could not find tool: " + name);
	return NULL;
}

ComputeTool* Grid::LoadRestartedTool(std::ifstream& inFile)
{
	// Read in the name and find the tool
	// No need to read data, it has already happened.
	std::string tmp;
	inFile >> tmp;
	inFile.ignore();
	return GetTool(tmp);
}

ComputeTool* Grid::ToolFromDirective(std::stringstream& inputString,const std::string& command)
{
	std::string word;
	tw::tool_type type;
	// Handle retrieval of named tools
	if (command=="get")
	{
		inputString >> word >> word >> word >> word >> word; // tool with name = [name]
		return GetTool(word);
	}
	// Handle creation of new tools on the fly
	type = ComputeTool::CreateTypeFromDirective(inputString,command);
	if (type!=tw::tool_type::nullTool)
		return CreateTool("tool",type);
	// If no tool return null pointer.  Modules have to be written to handle this.
	return NULL;
}

bool Grid::RemoveTool(ComputeTool *theTool)
{
	auto iter = std::find(computeTool.begin(),computeTool.end(),theTool);
	if (iter==computeTool.end())
		throw tw::FatalError("Attempt to remove a non-existant tool.");
	(*iter)->refCount--;
	if ((*iter)->refCount==0)
	{
		delete *iter;
		computeTool.erase(iter);
		return true;
	}
	return false;
}

void Grid::SetupTimeInfo(tw::Float dt)
{
	this->dt = dt;
	dth = 0.5*dt;

	tw::Int i;
	for (i=0;i<module.size();i++)
	{
		module[i]->dt = dt;
		module[i]->dth = dth;
		module[i]->dti = 1.0/dt;
	}
}

void Grid::MoveWindow()
{
	tw::Int i;
	windowPosition += spacing.z;
	corner.z += spacing.z;
	globalCorner.z += spacing.z;
	for (i=lb[3];i<=ub[3];i++)
		X(i,3) += spacing.z;

	for (i=0;i<clippingRegion.size();i++)
		if (clippingRegion[i]->moveWithWindow)
			clippingRegion[i]->Translate(tw::vec3(0,0,spacing.z));
		else
			clippingRegion[i]->Initialize(*this,this);

	for (i=0;i<module.size();i++)
		module[i]->MoveWindow();
}

void Grid::AntiMoveWindow()
{
	tw::Int i;
	antiWindowPosition -= spacing.z;

	for (i=0;i<module.size();i++)
		module[i]->AntiMoveWindow();
}

void Grid::ReadData(std::ifstream& inFile)
{
	tw::Int i;
	tw::Int num;

	Task::ReadData(inFile);
	MetricSpace::ReadData(inFile);
	inFile.read((char *)&gridGeometry,sizeof(tw_geometry));
	inFile.read((char *)&unitDensityCGS,sizeof(tw::Float));
	inFile.read((char *)&dt0,sizeof(tw::Float));
	inFile.read((char *)&dt,sizeof(tw::Float));
	inFile.read((char *)&dth,sizeof(tw::Float));
	inFile.read((char *)&dtMin,sizeof(tw::Float));
	inFile.read((char *)&dtMax,sizeof(tw::Float));
	inFile.read((char *)&elapsedTime,sizeof(tw::Float));
	inFile.read((char *)&elapsedTimeMax,sizeof(tw::Float));
	inFile.read((char *)&signalPosition,sizeof(tw::Float));
	inFile.read((char *)&signalSpeed,sizeof(tw::Float));
	inFile.read((char *)&windowPosition,sizeof(tw::Float));
	inFile.read((char *)&antiSignalPosition,sizeof(tw::Float));
	inFile.read((char *)&antiWindowPosition,sizeof(tw::Float));
	inFile.read((char *)&movingWindow,sizeof(bool));
	inFile.read((char *)&adaptiveTimestep,sizeof(bool));
	inFile.read((char *)&adaptiveGrid,sizeof(bool));
	inFile.read((char *)&appendMode,sizeof(bool));
	inFile.read((char *)&fullOutput,sizeof(bool));
	inFile.read((char *)&neutralize,sizeof(bool));
	inFile.read((char *)&smoothing,sizeof(tw::Int));
	inFile.read((char *)&compensation,sizeof(tw::Int));
	inFile.read((char *)&stepsToTake,sizeof(tw::Int));
	inFile.read((char *)&dumpPeriod,sizeof(tw::Int));
	inFile.read((char *)&binaryFormat,sizeof(tw::Int));
	inFile.read((char *)bc0,sizeof(tw_boundary_spec)*4);
	inFile.read((char *)bc1,sizeof(tw_boundary_spec)*4);
	inFile.read((char *)&radialProgressionFactor,sizeof(tw::Float));

	(*tw_out) << "Local Grid = " << dim[1] << "x" << dim[2] << "x" << dim[3] << std::endl;
	#ifdef USE_OPENCL
	InitializeMetricsBuffer(context,dt);
	#endif

	inFile.read((char *)&num,sizeof(tw::Int));
	for (i=0;i<num;i++)
	{
		(*tw_out) << "Add Nonuniform Region" << std::endl;
		region.push_back(new NonUniformRegion(1,2,1.0,1.0));
		inFile.read((char*)region.back(),sizeof(NonUniformRegion));
	}

	inFile.read((char *)&num,sizeof(tw::Int));
	for (i=1;i<num;i++) // don't read index 0, it is created by constructor
	{
		clippingRegion.push_back(Region::CreateObjectFromFile(clippingRegion,inFile));
		(*tw_out) << "Add Region " << clippingRegion.back()->name << std::endl;
	}

	if (uniformDeviate!=NULL) delete uniformDeviate;
	uniformDeviate = new UniformDeviate(1);
	uniformDeviate->ReadData(inFile);

	if (gaussianDeviate!=NULL) delete gaussianDeviate;
	gaussianDeviate = new GaussianDeviate(1);
	gaussianDeviate->ReadData(inFile);

	// Read ComputeTool objects

	inFile.read((char *)&num,sizeof(num));
	for (i=0;i<num;i++)
	{
		computeTool.push_back(ComputeTool::CreateObjectFromFile(inFile,this,this));
		(*tw_out) << "Installed Tool " << computeTool.back()->name << std::endl;
	}

	// Read Modules
	// Quasi-tools are read at the same time

	inFile.read((char *)&num,sizeof(tw::Int));
	for (i=0;i<num;i++)
	{
		module.push_back(Module::CreateObjectFromFile(inFile,this));
		// try to have Module::ReadData figure out to populate containment structures for both modules and tools.
		(*tw_out) << "Installed Module " << module.back()->name << std::endl;
	}

	// Read Grid managed objects
	// Probably most of them should be tools.

	inFile.read((char *)&num,sizeof(tw::Int));
	for (i=0;i<num;i++)
	{
		(*tw_out) << "Add Explicit Wave" << std::endl;
		wave.push_back(new Wave(gaussianDeviate));
		wave.back()->ReadData(inFile);
	}

	inFile.read((char *)&num,sizeof(tw::Int));
	for (i=0;i<num;i++)
	{
		(*tw_out) << "Add PGC Pulse" << std::endl;
		pulse.push_back(new Pulse(gaussianDeviate));
		pulse.back()->ReadData(inFile);
	}

	inFile.read((char *)&num,sizeof(tw::Int));
	for (i=0;i<num;i++)
	{
		(*tw_out) << "Add Conductor" << std::endl;
		conductor.push_back(new Conductor(clippingRegion));
		conductor.back()->ReadData(inFile);
	}

	// Read Grid managed diagnostics

	inFile.read((char *)&num,sizeof(tw::Int));
	for (i=0;i<num;i++)
	{
		(*tw_out) << "Add Energy Series" << std::endl;
		energyDiagnostic.push_back(new EnergySeriesDescriptor(clippingRegion));
		energyDiagnostic.back()->ReadData(inFile);
	}

	inFile.read((char *)&num,sizeof(tw::Int));
	for (i=0;i<num;i++)
	{
		(*tw_out) << "Add Point Series" << std::endl;
		pointDiagnostic.push_back(new PointSeriesDescriptor(clippingRegion));
		pointDiagnostic.back()->ReadData(inFile);
	}

	inFile.read((char *)&num,sizeof(tw::Int));
	for (i=0;i<num;i++)
	{
		(*tw_out) << "Add Box Diagnostic" << std::endl;
		boxDiagnostic.push_back(new GridDataDescriptor(clippingRegion));
		boxDiagnostic.back()->ReadData(inFile);
	}
}

void Grid::WriteData(std::ofstream& outFile)
{
	tw::Int i;

	Task::WriteData(outFile);
	MetricSpace::WriteData(outFile);
	outFile.write((char *)&gridGeometry,sizeof(tw_geometry));
	outFile.write((char *)&unitDensityCGS,sizeof(tw::Float));
	outFile.write((char *)&dt0,sizeof(tw::Float));
	outFile.write((char *)&dt,sizeof(tw::Float));
	outFile.write((char *)&dth,sizeof(tw::Float));
	outFile.write((char *)&dtMin,sizeof(tw::Float));
	outFile.write((char *)&dtMax,sizeof(tw::Float));
	outFile.write((char *)&elapsedTime,sizeof(tw::Float));
	outFile.write((char *)&elapsedTimeMax,sizeof(tw::Float));
	outFile.write((char *)&signalPosition,sizeof(tw::Float));
	outFile.write((char *)&signalSpeed,sizeof(tw::Float));
	outFile.write((char *)&windowPosition,sizeof(tw::Float));
	outFile.write((char *)&antiSignalPosition,sizeof(tw::Float));
	outFile.write((char *)&antiWindowPosition,sizeof(tw::Float));
	outFile.write((char *)&movingWindow,sizeof(bool));
	outFile.write((char *)&adaptiveTimestep,sizeof(bool));
	outFile.write((char *)&adaptiveGrid,sizeof(bool));
	outFile.write((char *)&appendMode,sizeof(bool));
	outFile.write((char *)&fullOutput,sizeof(bool));
	outFile.write((char *)&neutralize,sizeof(bool));
	outFile.write((char *)&smoothing,sizeof(tw::Int));
	outFile.write((char *)&compensation,sizeof(tw::Int));
	outFile.write((char *)&stepsToTake,sizeof(tw::Int));
	outFile.write((char *)&dumpPeriod,sizeof(tw::Int));
	outFile.write((char *)&binaryFormat,sizeof(tw::Int));
	outFile.write((char *)bc0,sizeof(tw_boundary_spec)*4);
	outFile.write((char *)bc1,sizeof(tw_boundary_spec)*4);
 	outFile.write((char *)&radialProgressionFactor,sizeof(tw::Float));

	i = region.size();
	outFile.write((char *)&i,sizeof(tw::Int));
	for (i=0;i<region.size();i++)
		outFile.write((char*)region[i],sizeof(NonUniformRegion));

	i = clippingRegion.size();
	outFile.write((char *)&i,sizeof(tw::Int));
	for (i=1;i<clippingRegion.size();i++) // don't write index 0, it is created by constructor
	clippingRegion[i]->WriteData(outFile);

	uniformDeviate->WriteData(outFile);
	gaussianDeviate->WriteData(outFile);

	i = module.size();
	outFile.write((char *)&i,sizeof(tw::Int));
	for (i=0;i<module.size();i++)
		module[i]->WriteData(outFile);

	i = wave.size();
	outFile.write((char *)&i,sizeof(tw::Int));
	for (i=0;i<wave.size();i++)
		wave[i]->WriteData(outFile);

	i = pulse.size();
	outFile.write((char *)&i,sizeof(tw::Int));
	for (i=0;i<pulse.size();i++)
		pulse[i]->WriteData(outFile);

	i = conductor.size();
	outFile.write((char *)&i,sizeof(tw::Int));
	for (i=0;i<conductor.size();i++)
		conductor[i]->WriteData(outFile);

	i = energyDiagnostic.size();
	outFile.write((char *)&i,sizeof(tw::Int));
	for (i=0;i<energyDiagnostic.size();i++)
		energyDiagnostic[i]->WriteData(outFile);

	i = pointDiagnostic.size();
	outFile.write((char *)&i,sizeof(tw::Int));
	for (i=0;i<pointDiagnostic.size();i++)
		pointDiagnostic[i]->WriteData(outFile);

	i = boxDiagnostic.size();
	outFile.write((char *)&i,sizeof(tw::Int));
	for (i=0;i<boxDiagnostic.size();i++)
		boxDiagnostic[i]->WriteData(outFile);
}


///////////////////////////
//  Read the Input File  //
///////////////////////////

void Grid::OpenInputFile(std::ifstream& inFile)
{
	std::string fileName;
	fileName = InputPathName() + "stdin";
	inFile.open(fileName.c_str());
	if (!(inFile.rdstate() & std::ios::failbit))
		return;
	fileName = InputPathName() + "stdin.txt";
	inFile.open(fileName.c_str());
	if (!(inFile.rdstate() & std::ios::failbit))
		return;
	throw tw::FatalError("couldn't open stdin");
}

std::string Grid::InputFileFirstPass()
{
	// The first pass is used to fully initialize the task

	try
	{
		Lock();

		bool foundGrid = false;
		bool foundRestart = false;
		std::stringstream messageOut,fileName;

		int numRanksProvided,worldRank;
		MPI_Comm_size(MPI_COMM_WORLD,&numRanksProvided);
		MPI_Comm_rank(MPI_COMM_WORLD,&worldRank);
		// world rank is suitable for reading task data from restart file
		// because this data is the same in every restart file

		std::ifstream inFile;
		std::stringstream inputString;

		OpenInputFile(inFile);

		PreprocessInputFile(inFile,inputString);
		inFile.close();

		inputString.seekg(0);

		std::string com1,com2,word;

		do
		{
			inputString >> com1;

			if (com1=="threads")
			{
				messageOut << "WARNING: threads directive no longer supported.  Use command line arguments instead." << std::endl;
			}

			if (com1=="affinity")
			{
				tw::input::ReadArray(affinityMask,inputString); // eg, affinity = { 0 2 4 6 }
			}

			if (com1=="hardware")
			{
				inputString >> word >> com1;
				if (com1=="device")
				{
					inputString >> com2;
					if (com2=="string")
						inputString >> word >> deviceSearchString; // eg, hardware acceleration device string = nvidia
					if (com2=="numbers")
						tw::input::ReadArray(deviceIDList,inputString); // eg, hardware acceleration device numbers = { 0 , 1 }
				}
				if (com1=="platform") // eg, hardware acceleration platform string = cuda
				{
					inputString >> word >> word >> platformSearchString;
				}
			}

			if (com1=="new")
			{
				inputString >> com1;
				if (com1=="grid")
				{
					foundGrid = true;

					do
					{
						inputString >> com2;

						if (com2=="dimensions") // eg, dimensions = 32 32 32
						{
							inputString >> word;
							inputString >> globalCells[1] >> globalCells[2] >> globalCells[3];
						}
						if (com2=="decomposition") // eg, decomposition = 8 4 1
						{
							inputString >> word;
							inputString >> domains[1] >> domains[2] >> domains[3];
						}
					} while (com2!="}");
				}
				else
					tw::input::ExitInputFileBlock(inputString);
			}

			if (com1=="generate")
				tw::input::ExitInputFileBlock(inputString);

			if (com1=="open")
			{
				inputString >> com2;
				if (com2=="restart")
				{
					foundRestart = true;
					restarted = true;
					inputString >> com2 >> com2;
					fileName.str("");
					fileName << InputPathName() << worldRank << "_" << com2;
					std::ifstream inFile(fileName.str().c_str());
					if (inFile.rdstate() & std::ios::failbit)
						throw tw::FatalError("could not open restart file");
					Task::ReadData(inFile);
					inFile.close();
				}
			}

			if (com1=="stdout") // eg, stdout = full
			{
				inputString >> com2 >> com2;
				if (com2=="full")
					fullOutput = true;
			}

			if (com1=="xboundary" || com1=="yboundary" || com1=="zboundary" ) // eg, xboundary = absorbing absorbing
			{
				tw::input::ReadBoundaryTerm(bc0,bc1,inputString,com1);
				periodic[1] = bc0[1]==cyclic ? 1 : 0;
				periodic[2] = bc0[2]==cyclic ? 1 : 0;
				periodic[3] = bc0[3]==cyclic ? 1 : 0;
			}

		} while (!inputString.eof());

		if (!foundGrid && !foundRestart)
			throw tw::FatalError("neither a grid directive nor a restart file was found");

		// Check integer viability
		int64_t totalCellsPerRank = int64_t(globalCells[1])*int64_t(globalCells[2])*int64_t(globalCells[3])/int64_t(numRanksProvided);
		if (totalCellsPerRank>=pow(2,31) && sizeof(tw::Int)==4)
			throw tw::FatalError("You must recompile turboWAVE with 64 bit integers to handle this many grid cells.");

		// Verify and (if necessary) correct decomposition
		if (NumTasks() != numRanksProvided)
		{
			messageOut << "WARNING: Bad decomposition ";
			tw::Int ax1=1,ax2=2,ax3=3; // to be sorted so ax1 is longest
			for (tw::Int i=1;i<=3;i++)
			{
				domains[i] = 1;
				if (globalCells[i]>=globalCells[1] && globalCells[i]>=globalCells[2] && globalCells[i]>=globalCells[3])
					ax1 = i;
			}
			for (tw::Int i=1;i<=3;i++)
				ax2 = i==ax1 ? ax2 : i;
			for (tw::Int i=1;i<=3;i++)
				ax3 = i==ax1 || i==ax2 ? ax3 : i;
			if (globalCells[ax2]<globalCells[ax3])
				std::swap(ax2,ax3);
			domains[ax1] = numRanksProvided;
			while (globalCells[ax1]%(domains[ax1]*2)!=0 && domains[ax1]>0)
			{
				domains[ax1] /= 2;
				domains[ax2] *= 2;
			}
			messageOut << "(defaulting to " << domains[1] << "x" << domains[2] << "x" << domains[3] << ")" << std::endl;
		}
		messageOut << NumTasks() << "-Way Decomposition" << std::endl;

		// Set up the domain decomposition
		// Grid/restart provide domains[] , globalCells[] , periodic[] as inputs
		// communicators, cornerCell[], localCells[], domainIndex[], n0[] , n1[] are computed
		for (tw::Int i=1;i<=3;i++)
		{
			if (globalCells[i]%domains[i]!=0)
				throw tw::FatalError("global number of cells is not divisible by number of domains along axis");
		}
		Task::Initialize(domains,globalCells,periodic);
		for (tw::Int i=1;i<=3;i++)
		{
			if (localCells[i]%2!=0 && globalCells[i]>1)
				throw tw::FatalError("local number of cells is not even along non-ignorable axis");
		}
		Resize(localCells[1],localCells[2],localCells[3],tw::vec3(0.0,0.0,0.0),tw::vec3(1.0,1.0,1.0),2);

		// Random numbers
		uniformDeviate = new UniformDeviate(1 + strip[0].Get_rank()*(MaxSeed()/numRanksProvided));
		gaussianDeviate = new GaussianDeviate(1 + strip[0].Get_rank()*(MaxSeed()/numRanksProvided) + MaxSeed()/(2*numRanksProvided));

		// Set up standard outputs
		std::stringstream stdoutString;
		stdoutString << strip[0].Get_rank() << "_stdout.txt";
		if (strip[0].Get_rank() == 0)
			tw_out = &std::cout;
		else
		{
			if (fullOutput)
				tw_out = new std::ofstream(stdoutString.str().c_str());
			else
				tw_out = new std::stringstream; // put output into a throw away string
		}

		Unlock();

		return messageOut.str();
	}

	catch (tw::FatalError& e)
	{
		std::cout << "FATAL ERROR: " << e.what() << std::endl;
		std::cout << "Could not start simulation --- exiting now." << std::endl;
		exit(1);
	}
}

void Grid::GridFromInputFile()
{
	Lock();

	tw::Int i;
	std::string com1,com2,word;

	// adaptive grid variables
	tw::Int i1,i2;
	tw::Float length;

	corner = globalCorner = size = globalSize = tw::vec3(0.0,0.0,0.0);

	(*tw_out) << "Extract Grid from Input File..." << std::endl << std::endl;

	std::ifstream inFile;
	std::stringstream inputString;
	OpenInputFile(inFile);
	PreprocessInputFile(inFile,inputString);
	inFile.close();

	inputString.seekg(0);

	do
	{
		inputString >> com1;

		if (com1=="new")
		{
			std::vector<std::string> preamble = tw::input::EnterInputFileBlock(inputString,"{");

			if (preamble[0]=="grid")
			{
				if (preamble.size()!=1)
					throw tw::FatalError("Ill-formed grid block.")
				do
				{
					inputString >> com2;
					if (com2=="corner") // eg, corner = 0.0 0.0 0.0
					{
						inputString >> word;
						inputString >> globalCorner.x >> globalCorner.y >> globalCorner.z;
					}
					if (com2=="cell") // eg, cell size = 0.5 0.5 0.5
					{
						inputString >> word >> word;
						inputString >> spacing.x >> spacing.y >> spacing.z;
					}
					if (com2=="geometry") // eg, geometry = cylindrical
					{
						gridGeometry = cartesian;
						inputString >> com2 >> com2;
						if (com2=="cylindrical")
							gridGeometry = cylindrical;
						if (com2=="spherical")
							gridGeometry = spherical;
					}
					if (com2=="radial") // eg, radial progression factor = 1.03
					{
						inputString >> com2 >> com2 >> com2 >> radialProgressionFactor;
					}
					if (com2=="region") // eg, region : start = 1 , end = 100 , length = 1e4
					{
						inputString >> com2 >> com2 >> com2 >> i1 >> com2 >> com2 >> i2 >> com2 >> com2 >> length;
						region.push_back(new NonUniformRegion(i1,i2,length,spacing.z));
					}

					if (com2=="adaptive") // eg, adaptive timestep = yes
					{
						inputString >> com2;
						if (com2=="timestep")
						{
							inputString >> com2 >> com2;
							adaptiveTimestep = (com2=="yes" || com2=="true" || com2=="on");
							(*tw_out) << "Adaptive timestep = " << adaptiveTimestep << std::endl;
						}
						if (com2=="grid")
						{
							inputString >> com2 >> com2;
							adaptiveGrid = (com2=="yes" || com2=="true" || com2=="on");
						}
					}
				} while (com2!="}");

				com2 = "grid";
			}
		}

		if (com1=="timestep") // eg, timestep = 1.0
		{
			inputString >> word;
			inputString >> dt0;
			SetupTimeInfo(dt0);
			(*tw_out) << "Timestep = " << dt << std::endl;
		}

		com1 = "???";

	} while (!inputString.eof());

	Unlock();

	(*tw_out) << "Allocate " << dim[1] << "x" << dim[2] << "x" << dim[3] << " Grid" << std::endl;
	size = spacing * tw::vec3(dim[1],dim[2],dim[3]);
	globalSize = spacing * tw::vec3(globalCells[1],globalCells[2],globalCells[3]);
	corner = globalCorner + tw::vec3(domainIndex[1],domainIndex[2],domainIndex[3]) * size;
	Resize(dim[1],dim[2],dim[3],corner,size,2);
	SetCellWidthsAndLocalSize();
	SetGlobalSizeAndLocalCorner();
	SetupGeometry();
	#ifdef USE_OPENCL
	InitializeMetricsBuffer(context,dt);
	#endif
}

void Grid::ReadSubmoduleBlock(std::stringstream& inputString,Module *sup)
{
	// To be called by supermodules that want to add submodules while
	// reading their own input file block.

	// Get the preamble = words that come between "new" and the opening brace
	std::vector<std::string> preamble = tw::input::EnterInputFileBlock(inputString,"{");
	// if an object has a name, it is expected to be the last string in the preamble
	std::string object_name(preamble.back());

	tw::module_type whichModule = Module::CreateTypeFromInput(preamble);
	if (whichModule!=tw::module_type::nullModule)
	{
		if (Module::SingularType(whichModule))
			if (module_type_exists(whichModule))
				throw tw::FatalError("Singular module type was created twice.  Check order of input file.");
		createdModuleTypes.push_back(whichModule);
		MangleModuleName(object_name);
		(*tw_out) << "Installing Module " << object_name << "..." << std::endl;
		module.push_back(Module::CreateObjectFromType(object_name,whichModule,this));
		module.back()->ReadInputFileBlock(inputString);
	}
	else
		throw tw::FatalError("Module type not recognized: " + preamble[0]);

	bool added = sup->AddSubmodule(module.back());
	if (!added)
		throw tw::FatalError("Unhandled " + preamble[0] + ". Check order of input file.");
}

void Grid::ReadInputFile()
{
	Lock();

	std::string com1,word;
	Profile* theProfile;
	std::ifstream inFile;
	std::stringstream inputString;

	(*tw_out) << std::endl << "Reading Input File..." << std::endl << std::endl;

	OpenInputFile(inFile);
	PreprocessInputFile(inFile,inputString);
	inFile.close();

	inputString.seekg(0);

	bool module_type_exists = [&] (tw::module_type whichType)
	{
		return std::find(createdModuleTypes.begin(),createdModuleTypes.end(),whichType) != createdModuleTypes.end();
	};

	do
	{
		inputString >> com1;

		if (com1=="new")
		{
			// Get the preamble = words that come between "new" and the opening brace
			std::vector<std::string> preamble = tw::input::EnterInputFileBlock(inputString,"{");
			// if an object has a name, it is expected to be the last string in the preamble
			std::string object_name(preamble.back());

			// Straight tool installation
			tw::tool_type whichTool = ComputeTool::CreateTypeFromInput(preamble);
			if (whichTool!=tw::tool_type::nullTool)
			{
				(*tw_out) << "Installing Tool " << object_name << "..." << std::endl;
				if (MangleToolName(object_name))
					throw tw::FatalError("Encountered duplicate tool name.");
				computeTool.push_back(ComputeTool::CreateObjectFromType(object_name,whichTool,this,this));
				computeTool.back()->ReadInputFileBlock(inputString);
			}

			// Handle modules whose creation is automatically triggered by another module
			// Must do this before reading in the submodule
			tw::module_type super_type = Module::CreateSupermoduleTypeFromSubmoduleKey(preamble[0]);
			if (super_type!=tw::module_type::nullModule)
				if (!module_type_exists(super_type) || !Module::SingularType(super_type))
				{
					(*tw_out) << "Installing super-Module automatically..." << std::endl;
					createdModuleTypes.push_back(super_type);
					object_name = "super";
					MangleModuleName(object_name);
					module.push_back(Module::CreateObjectFromType(object_name,super_type,this));
				}

			// Straight module installation
			tw::module_type whichModule = Module::CreateTypeFromInput(preamble);
			if (whichModule!=tw::module_type::nullModule)
			{
				if (Module::SingularType(whichModule))
					if (module_type_exists(whichModule))
						throw tw::FatalError("Singular module type was created twice.  Check order of input file.");
				createdModuleTypes.push_back(whichModule);
				MangleModuleName(object_name);
				(*tw_out) << "Installing Module " << object_name << "..." << std::endl;
				module.push_back(Module::CreateObjectFromType(object_name,whichModule,this));
				module.back()->ReadInputFileBlock(inputString);
			}

			// Handle containment of submodules
			// This is used for supermodules that were automatically created.
			// Submodules that are defined inside a supermodule block should never see this code.
			// The most recently created valid supermodule is selected.
			if (super_type!=tw::module_type::nullModule)
			{
				bool added = false;
				for (auto it=module.end();it>=module.begin();--it)
				{
					added = (*it)->AddSubmodule(module.back());
					if (added)
						break;
				}
				if (!added)
					throw tw::FatalError("Unhandled " + preamble[0] + ". Check order of input file.");
			}

			// Handle low level objects (not modules or tools) that should be owned by a module
			// The module must know how to read the block, or delegate it to the quasitool
			if (Module::QuasitoolNeedsModule(preamble))
			{
				bool processed = false;
				for (tw::Int i=0;i<module.size();i++)
					processed = processed || module[i]->ReadQuasitoolBlock(preamble,inputString)
				if (!processed)
					throw tw::FatalError("Unhandled " + preamble[0] + ". Check order of input file.");
			}

			// Need to add code to:
			// check if laser solvers installed a propagator
			// auto creation of groups (perhaps let chemistry object handle upon initialization)

			// The remaining objects are explicitly managed by Grid
			// Perhaps they should be repackaged as ComputeTool objects

			if (preamble[0]=="wave")
			{
				wave.push_back(new Wave(gaussianDeviate));
				wave.back()->ReadInputFile(inputString,preamble[0]);
				wave.back()->pulseShape.delay += dth;
			}

			if (preamble[0]=="pulse")
			{
				pulse.push_back(new Pulse(gaussianDeviate));
				pulse.back()->ReadInputFile(inputString,preamble[0]);
			}

			if (preamble[0]=="region")
			{
				clippingRegion.push_back(Region::CreateObjectFromString(clippingRegion,preamble[1]));
				clippingRegion.back()->name = preamble[2];
				clippingRegion.back()->ReadInputFileBlock(inputString);
			}

			if (preamble[0]=="conductor")
			{
				conductor.push_back(new Conductor(clippingRegion));
				conductor.back()->ReadInputFile(inputString,preamble[0]);
			}

			if (preamble[0]=="energy") // eg, new energy series { ... }
			{
				inputString >> word >> word;
				energyDiagnostic.push_back(new EnergySeriesDescriptor(clippingRegion));
				energyDiagnostic.back()->ReadInputFile(inputString);
			}

			if (preamble[0]=="point") // eg, new point series { ... }
			{
				inputString >> word	>> word;
				pointDiagnostic.push_back(new PointSeriesDescriptor(clippingRegion));
				pointDiagnostic.back()->ReadInputFile(inputString);
			}

			if (preamble[0]=="box") // eg, new box series { ... }
			{
				inputString >> word >> word;
				boxDiagnostic.push_back(new GridDataDescriptor(clippingRegion));
				boxDiagnostic.back()->ReadInputFile(inputString);
			}
		}

		if (com1=="generate")
		{
			std::string key,name;
			inputString >> key >> name >> word;
			theProfile = tw::input::GetProfile(this,name,key);

			if (theProfile!=NULL)
			{
				(*tw_out) << "Create " << key << " " << name << std::endl;
				theProfile->ReadInputFileBlock(inputString,neutralize);
				(*tw_out) << "   Clipping region = " << theProfile->theRgn->name << std::endl;
			}
			else
				(*tw_out) << "WARNING: Couldn't find " << name << std::endl;
		}

		// Outside declarations: must come after the above

		if (com1=="xboundary" || com1=="yboundary" || com1=="zboundary" ) // eg, xboundary = absorbing absorbing
		{
			// already done in first pass, but must take block off string again
			tw::input::ReadBoundaryTerm(bc0,bc1,inputString,com1);
		}

		if (com1=="open") // eg, open restart file dump1
		{
			std::stringstream fileName;
			std::ifstream restartFile;
			inputString >> word >> word >> word;
			fileName << InputPathName() << strip[0].Get_rank() << "_" << word;
			(*tw_out) << "Reading restart file " << fileName.str() << "..." << std::endl;
			restartFile.open(fileName.str().c_str());
			ReadData(restartFile);
			restartFile.close();
		}

		if (com1=="normalize" || com1=="unit") // eg, normalize density to 1e16, or unit density = 1e16
		{
			inputString >> word >> word;
			inputString >> unitDensityCGS;
			(*tw_out) << "Unit of density = " << unitDensityCGS << " cm^-3" << std::endl;
		}

		if (com1=="dtmax") // eg, dtmax = 100.0
		{
			inputString >> word;
			inputString >> dtMax;
			(*tw_out) << "Set Maximum Timestep = " << dtMax << std::endl;
		}

		if (com1=="dtmin") // eg, dtmin = 1.0
		{
			inputString >> word;
			inputString >> dtMin;
			(*tw_out) << "Set Minimum Timestep = " << dtMin << std::endl;
		}

		if (com1=="maxtime") // eg, maxtime = 1e4
		{
			inputString >> word;
			inputString >> elapsedTimeMax;
			(*tw_out) << "Set Maximum Elapsed Time = " << elapsedTimeMax << std::endl;
		}

		if (com1=="steps") // eg, steps = 10
		{
			inputString >> word;
			inputString >> stepsToTake;
			(*tw_out) << "Steps to Take = " << stepsToTake << std::endl;
		}

		if (com1=="dump") // eg, dump period = 1024
		{
			inputString >> word >> word >> dumpPeriod;
			(*tw_out) << "Dump Period = " << dumpPeriod << std::endl;
		}

		if (com1=="neutralize") // eg, neutralize = yes
		{
			inputString >> word >> word;
			neutralize = (word=="yes" || word=="true" || word=="on");
			(*tw_out) << "Full neutralization = " << neutralize << std::endl;
		}

		if (com1=="window") // eg, window speed = 1
		{
			inputString >> word >> word >> signalSpeed;
			(*tw_out) << "Window Speed = " << signalSpeed << std::endl;
		}

		if (com1=="moving") // eg, moving window = yes
		{
			inputString >> word >> word >> word;
			movingWindow = (word=="yes" || word=="on" || word=="true");
			(*tw_out) << "Moving Window = " << movingWindow << std::endl;
		}

		if (com1=="smoothing" || com1=="smoother") // eg, smoothing = on, or smoothing = 2 (smooth twice)
		{
			inputString >> word >> word;
			char *endptr;
			if (word=="yes" || word=="on" || word=="true")
				smoothing = 4;
			else
				smoothing = strtol(word.c_str(),&endptr,10);
			if (smoothing>0)
				compensation = 1;
		}

		if (com1=="compensation") // eg, compensation = 1
		{
			inputString >> word >> compensation;
		}

		if (com1=="binary")
		{
			inputString >> word >> word >> word;
			if (word=="2d")
				binaryFormat = 2;
			if (word=="3d")
				binaryFormat = 3;
		}

		if (com1=="append")
		{
			inputString >> word >> word >> word;
			appendMode = (word=="on" || word=="true" || word=="yes");
		}

		com1 = "???";

	} while (!inputString.eof());

	if (smoothing==0)
		(*tw_out) << "No Smoothing" << std::endl;
	else
		(*tw_out) << "Smoothing passes = " << smoothing << " , Compensation passes = " << compensation << std::endl;

	Unlock();
}
