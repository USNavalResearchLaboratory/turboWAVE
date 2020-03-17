#include "simulation.h"
#include "particles.h"
#include "fieldSolve.h"
#include "electrostatic.h"
#include "laserSolve.h"
#include "fluid.h"
#include "quantum.h"
#include "solidState.h"
//#include <unistd.h>

////////////////////////
//  SIMULATION CLASS  //
////////////////////////


Simulation::Simulation(const std::string& file_name)
{
	inputFileName = file_name;
	clippingRegion.push_back(new EntireRegion(clippingRegion));

	gridGeometry = tw::grid::cartesian;

	dt0 = 0.1;
	SetupTimeInfo(dt0);
	dtCritical = tw::small_pos;
	dtMin = tw::small_pos;
	dtMax = tw::big_pos;
	elapsedTime = 0.0;
	elapsedTimeMax = tw::big_pos;
	signalPosition = 0.0;
	windowPosition = 0.0;
	signalSpeed = 1.0;
	antiSignalPosition = 0.0;
	antiWindowPosition = 0.0;

	neutralize = true;
	movingWindow = false;
	restarted = false;
	completed = false;
	adaptiveTimestep = false;
	adaptiveGrid = false;

	stepNow = 1;
	stepsToTake = 32;
	lastTime = 0;
	outputLevel = 0;
	errorCheckingLevel = 0;
	dumpPeriod = 0;

	bc0[1] = tw::bc::par::periodic;
	bc1[1] = tw::bc::par::periodic;
	bc0[2] = tw::bc::par::periodic;
	bc1[2] = tw::bc::par::periodic;
	bc0[3] = tw::bc::par::absorbing;
	bc1[3] = tw::bc::par::absorbing;

	outerDirectives.Add("affinity",new tw::input::List<std::valarray<tw::Int>,tw::Int>(&affinityMask),false);
	outerDirectives.Add("hardware acceleration device string",new tw::input::String(&deviceSearchString),false);
	outerDirectives.Add("hardware acceleration device numbers",new tw::input::List<std::valarray<tw::Int>,tw::Int>(&deviceIDList),false);
	outerDirectives.Add("hardware acceleration platform string",new tw::input::String(&platformSearchString),false);
	outerDirectives.Add("timestep",new tw::input::Float(&dt0));
	outerDirectives.Add("output level",new tw::input::Int(&outputLevel),false);
	outerDirectives.Add("xboundary",new tw::input::Enums<tw::bc::par>(tw::bc::par_map(),&bc0[1],&bc1[1]));
	outerDirectives.Add("yboundary",new tw::input::Enums<tw::bc::par>(tw::bc::par_map(),&bc0[2],&bc1[2]));
	outerDirectives.Add("zboundary",new tw::input::Enums<tw::bc::par>(tw::bc::par_map(),&bc0[3],&bc1[3]));
	outerDirectives.Add("unit density",new tw::input::Float(&unitDensityCGS),false);
	outerDirectives.Add("dtmin",new tw::input::Float(&dtMin),false);
	outerDirectives.Add("dtmax",new tw::input::Float(&dtMax),false);
	outerDirectives.Add("dtcrit",new tw::input::Float(&dtCritical),false);
	outerDirectives.Add("maxtime",new tw::input::Float(&elapsedTimeMax),false);
	outerDirectives.Add("steps",new tw::input::Int(&stepsToTake));
	outerDirectives.Add("dump period",new tw::input::Int(&dumpPeriod),false);
	outerDirectives.Add("neutralize",new tw::input::Bool(&neutralize),false);
	outerDirectives.Add("window speed",new tw::input::Float(&signalSpeed),false);
	outerDirectives.Add("moving window",new tw::input::Bool(&movingWindow),false);
	outerDirectives.Add("error checking level",new tw::input::Int(&errorCheckingLevel),false);

	// allow user to choose any of 8 corners with strings like corner001, corner010, etc.
	tw::vec3 corner[2][2][2];
	auto corner_str = [&] (tw::Int i,tw::Int j,tw::Int k)
	{
		return "corner" + std::to_string(i) + std::to_string(j) + std::to_string(k);
	};
	for (tw::Int i=0;i<2;i++)
		for (tw::Int j=0;j<2;j++)
			for (tw::Int k=0;k<2;k++)
				gridDirectives.Add(corner_str(i,j,k),new tw::input::Numbers<tw::Float>(&corner[i][j][k][0],3),false);
	gridDirectives.Add("corner",new tw::input::Numbers<tw::Float>(&globalCorner[0],3),false);
	gridDirectives.Add("cell size",new tw::input::Numbers<tw::Float>(&spacing[0],3));
	gridDirectives.Add("adaptive timestep",new tw::input::Bool(&adaptiveTimestep),false);
	gridDirectives.Add("adaptive grid",new tw::input::Bool(&adaptiveGrid),false);
	gridDirectives.Add("dimensions",new tw::input::Numbers<tw::Int>(&globalCells[1],3));
	gridDirectives.Add("decomposition",new tw::input::Numbers<tw::Int>(&domains[1],3));
	std::map<std::string,tw::grid::geometry> geo = {{"cartesian",tw::grid::cartesian},{"cylindrical",tw::grid::cylindrical},{"spherical",tw::grid::spherical}};
	gridDirectives.Add("geometry",new tw::input::Enums<tw::grid::geometry>(geo,&gridGeometry),false);
}

Simulation::~Simulation()
{
	for (auto rgn : clippingRegion)
		delete rgn;

	for (auto m : module)
		delete m;
	// clean up after any modules that did not release tools
	for (auto tool : computeTool)
		delete tool;

	if (uniformDeviate!=NULL)
		delete uniformDeviate;
	if (gaussianDeviate!=NULL)
		delete gaussianDeviate;

	if (dynamic_cast<std::ofstream*>(tw_out))
		((std::ofstream*)tw_out)->close();
	if (dynamic_cast<std::stringstream*>(tw_out))
		delete tw_out;
}

void Simulation::SetCellWidthsAndLocalSize()
{
	// Overwrite uniform cell widths with requested variable widths.
	// Correct the local domain size to reflect the variable cell widths.

	for (auto tool : computeTool)
	{
		Warp *w = dynamic_cast<Warp*>(tool);
		if (w)
		{
			tw::Int ax = tw::grid::naxis(w->ax);
			for (tw::Int i=lfg[ax];i<=ufg[ax];i++)
				dX(i,ax) = spacing[ax-1] + w->AddedCellWidth(cornerCell[ax]-1+i);
		}
	}

	size = 0.0;
	for (tw::Int ax=1;ax<=3;ax++)
		for (tw::Int i=1;i<=dim[ax];i++)
			size[ax-1] += dX(i,ax);
}

void Simulation::SetGlobalSizeAndLocalCorner()
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

void Simulation::Run()
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

		#ifdef USE_HPC
		(*tw_out) << "Current status can be viewed in 'twstat' file." << std::endl;
		(*tw_out) << "This executable does not support interactive commands." << std::endl << std::endl;
		#endif
		#ifdef USE_DESKTOP
		(*tw_out) << "Current status can be viewed in 'twstat' file or by pressing enter key." << std::endl;
		(*tw_out) << "Enter 'help' for list of interactive commands." << std::endl << std::endl;
		#endif

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

void Simulation::SetupGeometry()
{
	// This routine assumes that MetricSpace::width, and MetricSpace::corner are valid
	switch (gridGeometry)
	{
		case tw::grid::cartesian:
			if (stepNow==1)
				(*tw_out) << "Using CARTESIAN Grid" << std::endl;
			SetCartesianGeometry();
			break;
		case tw::grid::cylindrical:
			if (stepNow==1)
				(*tw_out) << "Using CYLINDRICAL Grid" << std::endl;
			SetCylindricalGeometry();
			break;
		case tw::grid::spherical:
			if (stepNow==1)
				(*tw_out) << "Using SPHERICAL Grid" << std::endl;
			SetSphericalGeometry();
			break;
	}
}

void Simulation::PrepareSimulation()
{
	std::ofstream twstat;

	#ifdef USE_OPENCL
	PrintGPUInformation();
	#endif

	if (!restarted)
		SetupLocalGrid();

	ReadInputFile();

	// Attach clipping regions to tools
	for (auto tool : computeTool)
	{
		if (tool->region_name=="tw::entire")
			tool->theRgn = clippingRegion[0];
		else
			tool->theRgn = Region::FindRegion(clippingRegion,tool->region_name);
	}

	// The following is where Modules process the ComputeTool instances attached by the user.
	for (auto m : module)
		m->VerifyInput();

	// If a diagnostic tool is not attached to any module attach it to all modules
	for (auto tool : computeTool)
		if (dynamic_cast<Diagnostic*>(tool) && tool->refCount==0)
			for (auto m : module)
			{
				m->moduleTool.push_back(tool);
				tool->refCount++;
			}

	// The Task and MetricSpace inherited members are initialized during input file reading,
	// because Module constructors are allowed to assume the grid is fully specified.

	// Initialize Regions

	if (!restarted)
		for (tw::Int i=1;i<clippingRegion.size();i++)
			clippingRegion[i]->Initialize(*this,this);
	// region 0 is not saved in restart file
	clippingRegion[0]->Initialize(*this,this);

	// Sort Modules

	ModuleComparator comparatorObject;
	std::sort(module.begin(),module.end(),comparatorObject);

	// Initialize Computational Tools
	// Must precede module initialization

	(*tw_out) << std::endl << "Initializing Compute Tools..." << std::endl << std::endl;

	for (auto tool : computeTool)
	{
		(*tw_out) << "Tool: " << tool->name << std::endl;
		tool->Initialize();
		tool->WarningMessage(tw_out);
	}
	if (computeTool.size()==0)
		(*tw_out) << "(no tools)" << std::endl;

	// Initialize Modules

	(*tw_out) << std::endl << "Initializing Modules..." << std::endl << std::endl;

	for (auto m : module)
		m->ExchangeResources();

	for (auto m : module)
	{
		(*tw_out) << "Module: " << m->name << std::endl;
		m->Initialize();
		m->WarningMessage(tw_out);
	}
}

#ifdef USE_OPENCL

void Simulation::PrintGPUInformation()
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

void Simulation::InteractiveCommand(const std::string& cmd,std::ostream *theStream)
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
		for (auto m : module)
			m->StatusMessage(theStream);
		for (auto tool : computeTool)
			tool->StatusMessage(theStream);
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

void Simulation::FundamentalCycle()
{
	Diagnose();

	for (auto m : module)
		m->Reset();

		for (auto m : module)
		m->Update();

	elapsedTime += dt;
	signalPosition += signalSpeed*dt;
	antiSignalPosition -= signalSpeed*dt;
	stepNow++;

	if (adaptiveGrid)
		for (auto m : module)
			m->AdaptGrid();

	if (movingWindow && signalPosition>=(windowPosition + spacing.z) && dim[3]>1)
		MoveWindow();

	if (!movingWindow && antiSignalPosition<=(antiWindowPosition - spacing.z) && dim[3]>1)
		AntiMoveWindow();
}

bool Simulation::MangleModuleName(std::string& name)
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

bool Simulation::CheckModule(const std::string& name)
{
	for (tw::Int i=0;i<module.size();i++)
		if (module[i]->name==name)
			return true;
	return false;
}

Module* Simulation::GetModule(const std::string& name)
{
	for (tw::Int i=0;i<module.size();i++)
		if (module[i]->name==name)
			return module[i];
	throw tw::FatalError("Could not get module: <" + name+">");
	return NULL;
}

tw::Int Simulation::FindModule(const std::string& name)
{
	for (tw::Int i=0;i<module.size();i++)
		if (module[i]->name==name)
			return i;
	throw tw::FatalError("Could not find module: <" + name+">");
	return 0;
}

bool Simulation::MangleToolName(std::string& name)
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

ComputeTool* Simulation::CreateTool(const std::string& basename,tw::tool_type theType)
{
	std::string name(basename);
	MangleToolName(name);
	(*tw_out) << "Creating Tool <" << name << ">..." << std::endl;
	computeTool.push_back(ComputeTool::CreateObjectFromType(name,theType,this,this));
	computeTool.back()->refCount++;
	return computeTool.back();
}

ComputeTool* Simulation::GetTool(const std::string& name)
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

ComputeTool* Simulation::GetRestartedTool(std::ifstream& inFile)
{
	// Read in the name and find the tool
	// No need to read data, it has already happened.
	std::string tmp;
	inFile >> tmp;
	inFile.ignore();
	return GetTool(tmp);
}

void Simulation::ToolFromDirective(std::vector<ComputeTool*>& tool,std::stringstream& inputString,const std::string& command)
{
	// The first argument to this function is typically NOT the tool list owned by Simulation.
	// Instead it is the list owned by a module.

	std::string word;

	// Handle retrieval of named tools
	if (command=="get")
	{
		inputString >> word;
		if (word=="=")
			throw tw::FatalError("Expected a name after <get>, not the <=> separator.");
		tw::input::StripQuotes(word);
		if (CheckModule(word))
			throw tw::FatalError("Tried to <get> module "+word+", but <get> can only be used for tools.");
		tool.push_back(GetTool(word));
		return;
	}
}

bool Simulation::RemoveTool(ComputeTool *theTool)
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

void Simulation::MoveWindow()
{
	tw::Int i;
	windowPosition += spacing.z;
	corner.z += spacing.z;
	globalCorner.z += spacing.z;
	for (i=lfg[3];i<=ufg[3];i++)
		X(i,3) += spacing.z;

	for (i=0;i<clippingRegion.size();i++)
		if (clippingRegion[i]->moveWithWindow)
			clippingRegion[i]->Translate(tw::vec3(0,0,spacing.z));
		else
			clippingRegion[i]->Initialize(*this,this);

	for (i=0;i<module.size();i++)
		module[i]->MoveWindow();
}

void Simulation::AntiMoveWindow()
{
	tw::Int i;
	antiWindowPosition -= spacing.z;

	for (i=0;i<module.size();i++)
		module[i]->AntiMoveWindow();
}

void Simulation::ReadCheckpoint(std::ifstream& inFile)
{
	tw::Int i;
	tw::Int num;

	Task::ReadCheckpoint(inFile);
	MetricSpace::ReadCheckpoint(inFile);
	inFile.read((char *)&gridGeometry,sizeof(tw::grid::geometry));
	inFile.read((char *)&unitDensityCGS,sizeof(tw::Float));
	inFile.read((char *)&dt0,sizeof(tw::Float));
	inFile.read((char *)&dt,sizeof(tw::Float));
	inFile.read((char *)&dth,sizeof(tw::Float));
	inFile.read((char *)&dtCritical,sizeof(tw::Float));
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
	inFile.read((char *)&outputLevel,sizeof(outputLevel));
	inFile.read((char *)&errorCheckingLevel,sizeof(errorCheckingLevel));
	inFile.read((char *)&neutralize,sizeof(bool));
	inFile.read((char *)&stepsToTake,sizeof(tw::Int));
	inFile.read((char *)&dumpPeriod,sizeof(tw::Int));
	inFile.read((char *)bc0,sizeof(tw::bc::par)*4);
	inFile.read((char *)bc1,sizeof(tw::bc::par)*4);

	(*tw_out) << "Local Grid = " << dim[1] << "x" << dim[2] << "x" << dim[3] << std::endl;
	#ifdef USE_OPENCL
	InitializeMetricsBuffer(context,dt);
	#endif

	inFile.read((char *)&num,sizeof(tw::Int));
	for (i=1;i<num;i++) // don't read index 0, it is created by constructor
	{
		clippingRegion.push_back(Region::CreateObjectFromFile(clippingRegion,inFile));
		(*tw_out) << "Add Region " << clippingRegion.back()->name << std::endl;
	}

	if (uniformDeviate!=NULL) delete uniformDeviate;
	uniformDeviate = new UniformDeviate(1);
	uniformDeviate->ReadCheckpoint(inFile);

	if (gaussianDeviate!=NULL) delete gaussianDeviate;
	gaussianDeviate = new GaussianDeviate(1);
	gaussianDeviate->ReadCheckpoint(inFile);

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
		// Following calls Module::ReadCheckpoint, the base takes care of module containment, and populating the ComputeTool list.
		// Setup of ComputeTools is not completed until after input file processing, when Module::VerifyInput() is called.
		module.push_back(Module::CreateObjectFromFile(inFile,this));
		(*tw_out) << "Installed Module " << module.back()->name << std::endl;
	}
}

void Simulation::WriteCheckpoint(std::ofstream& outFile)
{
	tw::Int i;

	Task::WriteCheckpoint(outFile);
	MetricSpace::WriteCheckpoint(outFile);
	outFile.write((char *)&gridGeometry,sizeof(tw::grid::geometry));
	outFile.write((char *)&unitDensityCGS,sizeof(tw::Float));
	outFile.write((char *)&dt0,sizeof(tw::Float));
	outFile.write((char *)&dt,sizeof(tw::Float));
	outFile.write((char *)&dth,sizeof(tw::Float));
	outFile.write((char *)&dtCritical,sizeof(tw::Float));
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
	outFile.write((char *)&outputLevel,sizeof(outputLevel));
	outFile.write((char *)&errorCheckingLevel,sizeof(errorCheckingLevel));
	outFile.write((char *)&neutralize,sizeof(bool));
	outFile.write((char *)&stepsToTake,sizeof(tw::Int));
	outFile.write((char *)&dumpPeriod,sizeof(tw::Int));
	outFile.write((char *)bc0,sizeof(tw::bc::par)*4);
	outFile.write((char *)bc1,sizeof(tw::bc::par)*4);

	i = clippingRegion.size();
	outFile.write((char *)&i,sizeof(tw::Int));
	for (i=1;i<clippingRegion.size();i++) // don't write index 0, it is created by constructor
		clippingRegion[i]->WriteCheckpoint(outFile);

	uniformDeviate->WriteCheckpoint(outFile);
	gaussianDeviate->WriteCheckpoint(outFile);

	i = computeTool.size();
	outFile.write((char *)&i,sizeof(tw::Int));
	for (i=0;i<computeTool.size();i++)
		computeTool[i]->WriteCheckpoint(outFile);

	i = module.size();
	outFile.write((char *)&i,sizeof(tw::Int));
	for (i=0;i<module.size();i++)
		module[i]->WriteCheckpoint(outFile);
}


///////////////////////////
//  Read the Input File  //
///////////////////////////

void Simulation::OpenInputFile(std::ifstream& inFile)
{
	std::string inputPath;
	inputPath = InputPathName() + inputFileName;
	inFile.open(inputPath.c_str());
	if (inFile.rdstate() & std::ios::failbit)
		throw tw::FatalError("couldn't open input file " + inputPath);
}

std::string Simulation::InputFileFirstPass()
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

		tw::input::PreprocessInputFile(inFile,inputString);
		inFile.close();

		inputString.seekg(0);
		outerDirectives.Reset();
		gridDirectives.Reset();

		std::string com1,com2,word;

		do
		{
			com1 = outerDirectives.ReadNext(inputString);

			if (com1=="new")
			{
				tw::input::Preamble preamble = tw::input::EnterInputFileBlock(com1,inputString,"{=");
				if (preamble.end_token=="=")
				{
					// Pass through quasitools by searching for the next keyword.
					std::stringstream temp(inputString.str());
					temp.seekg(inputString.tellg());
					do
					{
						temp >> word;
						if (word!="new" && word!="generate" && word!="open")
							inputString >> word;
					} while (word!="new" && word!="generate" && word!="open");
				}
				else
				{
					if (preamble.words[0]=="grid")
					{
						if (preamble.words.size()!=1)
							throw tw::FatalError("Ill-formed grid block.");
						foundGrid = true;
						corner = globalCorner = size = globalSize = tw::vec3(0.0,0.0,0.0);
						// allow user to choose any of 8 corners with strings like corner001, corner010, etc.
						tw::vec3 corner[2][2][2];
						auto corner_str = [&] (tw::Int i,tw::Int j,tw::Int k)
						{
							return "corner" + std::to_string(i) + std::to_string(j) + std::to_string(k);
						};
						do
						{
							com2 = gridDirectives.ReadNext(inputString);
							if (com2=="tw::EOF")
								throw tw::FatalError("Encountered EOF while processing <grid>.");
							if (com2=="new" || com2=="generate" || com2=="get" || com2=="open")
								throw tw::FatalError("Keyword <"+com2+"> inside grid block is not allowed.");
						} while (com2!="}");
						gridDirectives.ThrowErrorIfMissingKeys("grid");
						tw::Int corners_given = gridDirectives.TestKey("corner") ? 1 : 0;
						for (tw::Int i=0;i<2;i++)
							for (tw::Int j=0;j<2;j++)
								for (tw::Int k=0;k<2;k++)
								{
									if (gridDirectives.TestKey(corner_str(i,j,k)))
									{
										corners_given++;
										if (corners_given>1)
											throw tw::FatalError("Grid geometry is overspecified.");
										globalCorner.x = corner[i][j][k].x - globalCells[1]*spacing.x*i;
										globalCorner.y = corner[i][j][k].y - globalCells[2]*spacing.y*j;
										globalCorner.z = corner[i][j][k].z - globalCells[3]*spacing.z*k;
									}
								}
					}
					else
						tw::input::ExitInputFileBlock(inputString,true);
				}
			}

			if (com1=="generate")
				tw::input::ExitInputFileBlock(inputString,false);

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
					Task::ReadCheckpoint(inFile);
					inFile.close();
				}
			}
		} while (!inputString.eof());

		if (outerDirectives.TestKey("timestep"))
			UpdateTimestep(dt0);
		else
			throw tw::FatalError("Could not find timestep directive.");
		periodic[1] = bc0[1]==tw::bc::par::periodic ? 1 : 0;
		periodic[2] = bc0[2]==tw::bc::par::periodic ? 1 : 0;
		periodic[3] = bc0[3]==tw::bc::par::periodic ? 1 : 0;

		if (!foundGrid && !foundRestart)
			throw tw::FatalError("neither a grid directive nor a restart file was found");

		if (outerDirectives.TestKey("unit density"))
			AttachUnits(unitDensityCGS);

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
		// Simulation/restart provide domains[] , globalCells[] , periodic[] as inputs
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
			if (outputLevel>0)
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

void Simulation::SetupLocalGrid()
{
	tw::Int Nx = localCells[1], Ny = localCells[2], Nz = localCells[3];
	(*tw_out) << "Allocate " << Nx << "x" << Ny << "x" << Nz << " Grid" << std::endl;
	size = spacing * tw::vec3(Nx,Ny,Nz);
	globalSize = spacing * tw::vec3(globalCells[1],globalCells[2],globalCells[3]);
	corner = globalCorner + tw::vec3(domainIndex[1],domainIndex[2],domainIndex[3]) * size;
	Resize(Nx,Ny,Nz,corner,size,2);
	SetCellWidthsAndLocalSize();
	SetGlobalSizeAndLocalCorner();
	SetupGeometry();
	#ifdef USE_OPENCL
	InitializeMetricsBuffer(context,dt);
	#endif
}

void Simulation::NestedDeclaration(const std::string& com,std::stringstream& inputString,Module *sup)
{
	// To be called by supermodules that want to add submodules or tools while
	// reading their own input file block.

	// Get the preamble = words that come between "new" and the opening brace
	tw::input::Preamble preamble = tw::input::EnterInputFileBlock(com,inputString,"{=");
	if (preamble.attaching)
		throw tw::FatalError(preamble.err_prefix+"keyword <for> is not allowed in a nested declaration.");

	auto module_type_exists = [&] (tw::module_type whichType)
	{
		return std::find(createdModuleTypes.begin(),createdModuleTypes.end(),whichType) != createdModuleTypes.end();
	};

	tw::module_type whichModule = Module::CreateTypeFromInput(preamble);
	tw::tool_type whichTool = ComputeTool::CreateTypeFromInput(preamble);
	if (whichModule==tw::module_type::none && whichTool==tw::tool_type::none)
		throw tw::FatalError(preamble.err_prefix+"key was not recognized.");
	if (whichModule!=tw::module_type::none && whichTool!=tw::tool_type::none)
		throw tw::FatalError(preamble.err_prefix+"key claimed by both Module and Tool, this is a bug in the code.");

	if (whichModule!=tw::module_type::none)
	{
		if (Module::SingularType(whichModule))
			if (module_type_exists(whichModule))
				throw tw::FatalError(preamble.err_prefix+"Singular module type was created twice.  Check order of input file.");
		createdModuleTypes.push_back(whichModule);
		MangleModuleName(preamble.obj_name);
		(*tw_out) << "   Attaching nested module <" << preamble.obj_name << ">..." << std::endl;
		module.push_back(Module::CreateObjectFromType(preamble.obj_name,whichModule,this));
		module.back()->ReadInputFileBlock(inputString);
		bool added = sup->AddSubmodule(module.back());
		if (!added)
			throw tw::FatalError(preamble.err_prefix+"parent module rejected the child.");
	}

	if (whichTool!=tw::tool_type::none)
	{
		ComputeTool *tool = CreateTool(preamble.obj_name,whichTool);
		(*tw_out) << "   Attaching nested tool <" << tool->name << ">..." << std::endl;
		tool->ReadInputFileBlock(inputString);
		sup->moduleTool.push_back(tool);
	}
}

void Simulation::ReadInputFile()
{
	Lock();

	std::string com1,word;
	Profile* theProfile;
	std::ifstream inFile;
	std::stringstream inputString;

	(*tw_out) << std::endl << "Reading Input File..." << std::endl << std::endl;

	OpenInputFile(inFile);
	tw::input::PreprocessInputFile(inFile,inputString);
	inFile.close();

	inputString.seekg(0);
	outerDirectives.Reset();
	gridDirectives.Reset();

	auto module_type_exists = [&] (tw::module_type whichType)
	{
		return std::find(createdModuleTypes.begin(),createdModuleTypes.end(),whichType) != createdModuleTypes.end();
	};

	auto find_super = [&] (Module *sub)
	{
		// used for submodules defined outside the supermodule's block
		bool added = false;
		for (auto it=module.end()-1;it>=module.begin();--it)
		{
			added = (*it)->AddSubmodule(sub);
			if (added)
				break;
		}
		if (added)
			(*tw_out) << "   (super=" << sub->super->name << ")" << std::endl;
		else
			(*tw_out) << "   (super=none)" << std::endl;
	};

	do
	{
		// Process outer directives again, after first pass.
		// This is redundant, but helps minimize keyword collisions with object names.
		com1 = outerDirectives.ReadNext(inputString);

		if (com1=="new" || com1=="generate")
		{
			bool processed = false;

			// Get the preamble = words that come between "new" and the opening token, and derived information
			tw::input::Preamble preamble = tw::input::EnterInputFileBlock(com1,inputString,"{=");

			// Install a pre or post declared tool
			tw::tool_type whichTool = ComputeTool::CreateTypeFromInput(preamble);
			if (whichTool!=tw::tool_type::none)
			{
				processed = true;
				if (preamble.attaching)
				{
					ComputeTool *tool = CreateTool(preamble.obj_name,whichTool);
					(*tw_out) << "Attaching <" << tool->name << "> to <" << preamble.owner_name << ">..." << std::endl;
					tool->ReadInputFileBlock(inputString);
					GetModule(preamble.owner_name)->moduleTool.push_back(tool);
				}
				else
				{
					bool duplicate = MangleToolName(preamble.obj_name);
					(*tw_out) << "Creating Tool <" << preamble.obj_name << ">..." << std::endl;
					if (duplicate && errorCheckingLevel>0)
						throw tw::FatalError(preamble.err_prefix+"duplicate tool name.");
					// Do not use CreateTool, do not want to increase refCount
					computeTool.push_back(ComputeTool::CreateObjectFromType(preamble.obj_name,whichTool,this,this));
					computeTool.back()->ReadInputFileBlock(inputString);
				}
			}

			// Handle modules whose creation is automatically triggered by another module
			// Must do this before reading in the submodule
			tw::module_type super_type = Module::CreateSupermoduleTypeFromSubmoduleKey(preamble.words[0]);
			if (super_type!=tw::module_type::none)
				if (!module_type_exists(super_type) || !Module::SingularType(super_type))
				{
					createdModuleTypes.push_back(super_type);
					std::string super_module_name = preamble.obj_name + "_sup";
					MangleModuleName(super_module_name);
					(*tw_out) << "Installing supermodule triggered by <" << preamble.obj_name << ">..." << std::endl;
					module.push_back(Module::CreateObjectFromType(super_module_name,super_type,this));
					find_super(module.back());
				}

			// Straight module installation
			tw::module_type whichModule = Module::CreateTypeFromInput(preamble);
			if (whichModule!=tw::module_type::none)
			{
				processed = true;
				if (Module::SingularType(whichModule))
					if (module_type_exists(whichModule))
						throw tw::FatalError(preamble.err_prefix + "singular module type was created twice.  Check order of input file.");
				createdModuleTypes.push_back(whichModule);
				MangleModuleName(preamble.obj_name);
				(*tw_out) << "Installing module <" << preamble.obj_name << ">..." << std::endl;
				module.push_back(Module::CreateObjectFromType(preamble.obj_name,whichModule,this));
				find_super(module.back()); // note next line might change module.back()
				module.back()->ReadInputFileBlock(inputString);
			}

			// Handle low level objects (not modules or tools) that should be owned by a module
			// The module must know how to read the block, or delegate it to the quasitool
			if (Module::QuasitoolNeedsModule(preamble))
			{
				(*tw_out) << "Processing quasitool <" << preamble.obj_name << ">..." << std::endl;
				for (tw::Int i=0;i<module.size();i++)
					processed = processed || module[i]->ReadQuasitoolBlock(preamble,inputString);
				if (!processed)
					throw tw::FatalError("Unhandled key <" + preamble.words[0] + ">. Check order of input file.");
			}

			// Regions are neither modules nor tools

			if (preamble.words[0]=="region")
			{
				processed = true;
				clippingRegion.push_back(Region::CreateObjectFromString(clippingRegion,preamble.words[1]));
				clippingRegion.back()->name = preamble.obj_name;
				clippingRegion.back()->ReadInputFileBlock(inputString);
			}

			if (preamble.words[0]=="grid") // drop out of grid block, it was already processed during first pass.
			{
				processed = true;
				tw::input::ExitInputFileBlock(inputString,true);
			}

			if (!processed)
				throw tw::FatalError(preamble.err_prefix+"keys were not understood.");
		}

		if (com1=="open") // eg, open restart file dump1
		{
			std::stringstream fileName;
			std::ifstream restartFile;
			inputString >> word >> word >> word;
			fileName << InputPathName() << strip[0].Get_rank() << "_" << word;
			(*tw_out) << "Reading restart file " << fileName.str() << "..." << std::endl;
			restartFile.open(fileName.str().c_str());
			ReadCheckpoint(restartFile);
			restartFile.close();
		}

	} while (!inputString.eof());

	outerDirectives.ThrowErrorIfMissingKeys("Simulation");

	Unlock();
}

void Simulation::Diagnose()
{
	tw::Int master = 0;
	tw::Int curr = strip[0].Get_rank();

	std::vector<Diagnostic*> diagnostic;
	for (auto tool : computeTool)
		if (dynamic_cast<Diagnostic*>(tool))
			diagnostic.push_back(dynamic_cast<Diagnostic*>(tool));

	// DIAGNOSTIC PREP
	// Main purpose of this step is to let modules transfer data from compute devices.
	// This has a high cost, so only do if necessary.

	bool doing_restart=false, doing_diagnostics=false;
	doing_restart = dumpPeriod>0 && stepNow%dumpPeriod==0;
	for (auto d : diagnostic)
		doing_diagnostics = doing_diagnostics || d->WriteThisStep(elapsedTime,dt,stepNow);
	if (doing_restart || doing_diagnostics)
	{
		for (auto m : module)
			m->StartDiagnostics();
	}
	else
	{
		return;
	}

	// RESTART MECHANISM

	Lock();
	if (doing_restart)
	{
		std::ofstream restartFile;
		std::string fileName = std::to_string(curr) + "_dump" + std::to_string(stepNow/dumpPeriod);
		restartFile.open(fileName.c_str(),std::ios::binary);
		WriteCheckpoint(restartFile);
		restartFile.close();
	}
	Unlock();

	// MAIN DIAGNOSTIC LOOP

	auto has_diagnostic = [&] (Module *m,Diagnostic *diagnostic)
	{
		for (auto d : m->moduleTool)
			if (d==diagnostic)
				return true;
		return false;
	};

	for (auto d : diagnostic)
	{
		if (d->WriteThisStep(elapsedTime,dt,stepNow))
		{
			d->Start();
			d->Float("time",elapsedTime,true);
			d->Float("dt",dt,true);
			d->Float("zwindow",windowPosition,true);
			for (auto m : module)
				if (has_diagnostic(m,d))
					m->Report(*d);
			d->Finish();
		}
	}
}
