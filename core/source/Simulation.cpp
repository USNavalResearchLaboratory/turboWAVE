#include "simulation.h"
#include "particles.h"
#include "fieldSolve.h"
#include "electrostatic.h"
#include "laserSolve.h"
#include "fluid.h"
#include "quantum.h"
#include "solidState.h"
//#include <unistd.h>

GridReader::GridReader(tw::UnitConverter& uc)
{
	globalCorner = 0.0;
	adaptiveGrid = adaptiveTimestep = false;
	geo = tw::grid::cartesian;

	// allow user to choose any of 8 corners with strings like corner001, corner010, etc.
	auto corner_str = [&] (tw::Int i,tw::Int j,tw::Int k)
	{
		return "corner" + std::to_string(i) + std::to_string(j) + std::to_string(k);
	};
	for (tw::Int i=0;i<2;i++)
		for (tw::Int j=0;j<2;j++)
			for (tw::Int k=0;k<2;k++)
				directives.Add(corner_str(i,j,k),new tw::input::Numbers<tw::Float>(&cornerSet[i][j][k][0],3),false);
	directives.Add("corner",new tw::input::Numbers<tw::Float>(&globalCorner[0],3),false);
	directives.Add("cell size",new tw::input::Numbers<tw::Float>(&spacing[0],3));
	directives.Add("adaptive timestep",new tw::input::Bool(&adaptiveTimestep),false);
	directives.Add("adaptive grid",new tw::input::Bool(&adaptiveGrid),false);
	directives.Add("dimensions",new tw::input::Numbers<tw::Int>(&req_dim[1],3));
	directives.Add("decomposition",new tw::input::Numbers<tw::Int>(&req_dom[1],3));
	std::map<std::string,tw::grid::geometry> geo_map = {{"cartesian",tw::grid::cartesian},{"cylindrical",tw::grid::cylindrical},{"spherical",tw::grid::spherical}};
	directives.Add("geometry",new tw::input::Enums<tw::grid::geometry>(geo_map,&geo),false);
	directives.AttachUnits(uc);
}

void GridReader::Read(std::stringstream& inputString,tw::input::Preamble& preamble)
{
	std::string com2;
	if (preamble.words.size()!=1)
		throw tw::FatalError("Ill-formed grid block.");
	// allow user to choose any of 8 corners with strings like corner001, corner010, etc.
	auto corner_str = [&] (tw::Int i,tw::Int j,tw::Int k)
	{
		return "corner" + std::to_string(i) + std::to_string(j) + std::to_string(k);
	};
	do
	{
		com2 = directives.ReadNext(inputString);
		if (com2=="tw::EOF")
			throw tw::FatalError("Encountered EOF while processing <grid>.");
		if (com2=="new" || com2=="generate" || com2=="get")
			throw tw::FatalError("Keyword <"+com2+"> inside grid block is not allowed.");
	} while (com2!="}");
	directives.ThrowErrorIfMissingKeys("grid");
	tw::Int corners_given = directives.TestKey("corner") ? 1 : 0;
	for (tw::Int i=0;i<2;i++)
		for (tw::Int j=0;j<2;j++)
			for (tw::Int k=0;k<2;k++)
			{
				if (directives.TestKey(corner_str(i,j,k)))
				{
					corners_given++;
					if (corners_given>1)
						throw tw::FatalError("Grid geometry is overspecified.");
					globalCorner.x = cornerSet[i][j][k].x - req_dim[1]*spacing.x*i;
					globalCorner.y = cornerSet[i][j][k].y - req_dim[2]*spacing.y*j;
					globalCorner.z = cornerSet[i][j][k].z - req_dim[3]*spacing.z*k;
				}
			}
}

void GridReader::UpdateTask(Task& tsk)
{
	for (tw::Int i=0;i<4;i++)
	{
		tsk.globalCells[i] = req_dim[i];
		tsk.domains[i] = req_dom[i];
	}
}


////////////////////////
//  SIMULATION CLASS  //
////////////////////////


Simulation::Simulation(const std::string& test_name,
	const std::string& file_name,
	const std::string& restart_name,
	const std::string& platform,
	const std::string& device,
	const tw::Int& outputLevel,
	const tw::Int& errorCheckingLevel)
{
	unitTest = test_name;
	inputFileName = file_name;
	restartFileName = restart_name;
	platformSearchString = platform;
	deviceSearchString = device;
	this->outputLevel = outputLevel;
	this->errorCheckingLevel = errorCheckingLevel;
	clippingRegion.push_back(new EntireRegion(clippingRegion));

	nativeUnits = tw::units::plasma;
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
	completed = false;

	stepNow = 0;
	stepsToTake = 32;
	lastTime = 0;
	dumpPeriod = 0;

	bc0[1] = tw::bc::par::periodic;
	bc1[1] = tw::bc::par::periodic;
	bc0[2] = tw::bc::par::periodic;
	bc1[2] = tw::bc::par::periodic;
	bc0[3] = tw::bc::par::absorbing;
	bc1[3] = tw::bc::par::absorbing;

	outerDirectives.Add("native units",new tw::input::Enums<tw::units>(tw::get_unit_map(),&nativeUnits),false);
	outerDirectives.Add("unit density",new tw::input::Float(&unitDensityCGS));
	outerDirectives.Add("timestep",new tw::input::Float(&dt0));
	outerDirectives.Add("xboundary",new tw::input::Enums<tw::bc::par>(tw::bc::par_map(),&bc0[1],&bc1[1]));
	outerDirectives.Add("yboundary",new tw::input::Enums<tw::bc::par>(tw::bc::par_map(),&bc0[2],&bc1[2]));
	outerDirectives.Add("zboundary",new tw::input::Enums<tw::bc::par>(tw::bc::par_map(),&bc0[3],&bc1[3]));
	outerDirectives.Add("dtmin",new tw::input::Float(&dtMin),false);
	outerDirectives.Add("dtmax",new tw::input::Float(&dtMax),false);
	outerDirectives.Add("dtcrit",new tw::input::Float(&dtCritical),false);
	outerDirectives.Add("maxtime",new tw::input::Float(&elapsedTimeMax),false);
	outerDirectives.Add("steps",new tw::input::Int(&stepsToTake));
	outerDirectives.Add("dump period",new tw::input::Int(&dumpPeriod),false);
	outerDirectives.Add("neutralize",new tw::input::Bool(&neutralize),false);
	outerDirectives.Add("window speed",new tw::input::Float(&signalSpeed),false);
	outerDirectives.Add("moving window",new tw::input::Bool(&movingWindow),false);
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

void Simulation::SetupIO()
{
	// Set up standard outputs
	if (strip[0].Get_rank()==0)
	{
		tw_out = &std::cout;
		tw_err = &std::cerr;
	}
	else
	{
		if (outputLevel>0)
		{
			tw_out = new std::ofstream(std::to_string(strip[0].Get_rank()) + "_stdout.txt");
			tw_err = new std::ofstream(std::to_string(strip[0].Get_rank()) + "_stderr.txt");
		}
		else
		{
			// put outputs into throw-away strings
			tw_out = new std::stringstream;
			tw_err = new std::stringstream;
		}
	}

	*tw_out << std::endl << "*** Starting turboWAVE Session ***" << std::endl << std::endl;
	*tw_out << "Floating point precision = " << sizeof(tw::Float)*8 << " bits" << std::endl;
	#ifdef USE_OPENMP
	*tw_out << "Maximum OpenMP threads = " << omp_get_max_threads() << std::endl;
	#endif
}

void Simulation::Run()
{
	#ifdef USING_TW_MPI
	// Wait until master thread is done
	// TODO: fold this into TW_MPI's MPI_init
	TW_MPI_Lock();
	TW_MPI_Unlock();
	#endif

	SetupIO();
	InputFileFirstPass();

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

		(*tw_out) << "Completed " << stepNow << " steps in " << GetSeconds() - startTime << " seconds." << std::endl;
		(*tw_out) << "Simulated elapsed time = " << elapsedTime << std::endl;
		if (strip[0].Get_rank()==0)
		{
			twstat.open("twstat");
			twstat << "Completed " << stepNow << " steps in " << GetSeconds() - startTime << " seconds." << std::endl;
			twstat << "Simulated elapsed time = " << elapsedTime << std::endl;
			twstat.close();
		}
	}
	catch (tw::FatalError& e)
	{
		(*tw_err) << "FATAL ERROR: " << e.what() << std::endl;
		(*tw_err) << "Simulation failed --- exiting now." << std::endl;
		#ifdef USE_TW_MPI
		if (tw_out != &std::cout)
		{
			std::cerr << "FATAL ERROR: " << e.what() << std::endl;
			std::cerr << "Simulation failed --- exiting now." << std::endl;
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

	*tw_out << std::endl << "*** TurboWAVE Session Completed ***" << std::endl;
	completed = true;
}

void Simulation::Test()
{
	#ifdef USING_TW_MPI
	// Wait until master thread is done
	// TODO: fold this into TW_MPI's MPI_init
	TW_MPI_Lock();
	TW_MPI_Unlock();
	#endif

	SetupIO();

	// Get basic MPI data
	int numRanksProvided,worldRank;
	MPI_Comm_size(MPI_COMM_WORLD,&numRanksProvided);
	MPI_Comm_rank(MPI_COMM_WORLD,&worldRank);

	// Random numbers
	uniformDeviate = new UniformDeviate(1 + strip[0].Get_rank()*(MaxSeed()/numRanksProvided));
	gaussianDeviate = new GaussianDeviate(1 + strip[0].Get_rank()*(MaxSeed()/numRanksProvided) + MaxSeed()/(2*numRanksProvided));

	// diverting output
	std::ostream *alt_out = new std::stringstream;
	std::ostream *save_out = tw_out;

	// Create a test grid
	if (numRanksProvided!=2)
		throw tw::FatalError("test must use -n 2");
	AttachUnits(tw::units::plasma,1e19);
	globalCells[1] = 4;
	globalCells[2] = 4;
	globalCells[3] = 4;
	periodic[1] = 1;
	periodic[2] = 1;
	periodic[3] = 0;
	domains[1] = 1;
	domains[2] = 1;
	domains[3] = 2;
	Task::Initialize(domains,globalCells,periodic);
	globalCorner = tw::vec3(0.0,0.0,0.0);
	spacing = tw::vec3(0.2,0.2,0.2);
	globalSize = spacing * tw::vec3(globalCells[1],globalCells[2],globalCells[3]);
	Resize(*this,globalCorner,globalSize,2);
	UpdateTimestep(0.1);

	// Verify grid creation
	tw::Float tolerance = 1e-10;
	for (tw::Int ax=1;ax<=3;ax++)
	{
		assert(size[ax-1]-spacing[ax-1]*dim[ax]<tolerance);
		assert(corner[ax-1]-(globalCorner[ax-1]+size[ax-1]*strip[ax].Get_rank())<tolerance);
	}
	tw::Int src,dst;
	strip[1].Shift(1,1,&src,&dst);
	assert(src==strip[1].Get_rank());
	assert(dst==strip[1].Get_rank());
	strip[2].Shift(1,1,&src,&dst);
	assert(src==strip[2].Get_rank());
	assert(dst==strip[2].Get_rank());
	strip[3].Shift(1,1,&src,&dst);
	assert(src==strip[3].Get_rank()-1 || src==MPI_PROC_NULL);
	assert(dst==strip[3].Get_rank()+1 || dst==MPI_PROC_NULL);

	bool some_success = false;
	ComputeTool *tool;
	Module *module;

	for (auto m : ComputeTool::Map())
	{
		if (unitTest==m.first || unitTest=="--all")
		{
			tw_out = save_out;
			*tw_out << std::endl << "testing " << term::bold << term::cyan << m.first << term::reset_all << std::endl;
			try
			{
				tw_out = alt_out;
				tool = CreateTool("test_tool",m.second);
			}
			catch (tw::FatalError& e)
			{
				tw_out = save_out;
				*tw_out << "    " << term::yellow << "tool rejected the environment" << term::reset_all << std::endl;
				tool = NULL;
			}
			tw_out = alt_out;
			if (tool!=NULL)
			{
				if (tool->Test())
				{
					tw_out = save_out;
					*tw_out << "    " << term::ok << " " << term::green << "success" << term::reset_all << std::endl;
					some_success = true;
				}
				else
				{
					tw_out = save_out;
					*tw_out << "    " << term::yellow << "no test available for this environment" << term::reset_all << std::endl;
				}
				RemoveTool(tool);
			}
		}
	}

	for (auto m : Module::Map())
	{
		if (unitTest==m.first || unitTest=="--all")
		{
			tw_out = save_out;
			*tw_out << std::endl << "testing " << term::bold << term::cyan << m.first << term::reset_all << std::endl;
			try
			{
				tw_out = alt_out;
				module = Module::CreateObjectFromType("test_module",m.second,this);
			}
			catch (tw::FatalError& e)
			{
				tw_out = save_out;
				*tw_out << "    " << term::yellow << "module rejected the environment" << term::reset_all << std::endl;
				module = NULL;
			}
			tw_out = alt_out;
			if (module!=NULL)
			{
				if (module->Test())
				{
					tw_out = save_out;
					*tw_out << "    " << term::ok << " " << term::green << "success" << term::reset_all << std::endl;
					some_success = true;
				}
				else
				{
					tw_out = save_out;
					*tw_out << "    " << term::yellow << "no test available for this environment" << term::reset_all << std::endl;
				}
				delete module;
			}
		}
	}

	tw_out = save_out;
	delete alt_out;
	assert(some_success);
	completed = true;
}

void Simulation::PrepareSimulation()
{
	std::ofstream twstat;

	#ifdef USE_OPENCL
	PrintGPUInformation();
	#endif

	ReadInputFile();

	// The following is where Modules process the ComputeTool instances attached by the user.
	for (auto m : module)
		m->VerifyInput();

	// Attach clipping regions to tools
	for (auto tool : computeTool)
	{
		if (tool->region_name=="tw::entire")
			tool->theRgn = clippingRegion[0];
		else
			tool->theRgn = Region::FindRegion(clippingRegion,tool->region_name);
	}

	// Start the metadata dictionary
	if (strip[0].Get_rank()==0)
	{
		std::ofstream metadata_file("tw_metadata.json");
		metadata_file << "{" << std::endl;
		metadata_file << "\"version\": \"" << TW_VERSION_STRING << "\"}";
		// n.b. closing brace must not be followed by anything due to append strategy
		metadata_file.close();
	}

	// If a diagnostic tool is not attached to any module attach it to all modules
	for (auto tool : computeTool)
		if (dynamic_cast<Diagnostic*>(tool) && tool->refCount==0)
			for (auto m : module)
			{
				m->moduleTool.push_back(tool);
				tool->refCount++;
			}

	// Check for duplicate filenames
	std::vector<std::string> fileName;
	for (auto tool : computeTool)
	{
		Diagnostic *d = dynamic_cast<Diagnostic*>(tool);
		if (d)
		{
			if (std::find(fileName.begin(),fileName.end(),d->filename)==fileName.end())
				fileName.push_back(d->filename);
			else
				throw tw::FatalError("Duplicate file name used in diagnostic <"+d->name+">.");
		}
	}

	// The Task and MetricSpace inherited members are initialized during input file reading,
	// because Module constructors are allowed to assume the grid is fully specified.

	// Initialize Regions

	for (auto rgn : clippingRegion)
		rgn->Initialize(*this,this);

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

	// Read checkpoint data

	if (restartFileName!="tw::none")
	{
		std::stringstream fileName;
		std::ifstream restartFile;
		fileName << strip[0].Get_rank() << "_dump.chk";
		(*tw_out) << std::endl << "Reading restart file " << fileName.str() << "..." << std::endl << std::endl;
		restartFile.open(fileName.str().c_str());
		ReadCheckpoint(restartFile);
		restartFile.close();
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

	// RESTART MECHANISM

	tw::Int curr = strip[0].Get_rank();
	bool doing_restart = dumpPeriod>0 && stepNow%dumpPeriod==0;
	// Lock();
	if (doing_restart)
	{
		std::ofstream restartFile;
		std::string fileName = std::to_string(curr) + "_dump.chk";
		restartFile.open(fileName.c_str(),std::ios::binary);
		WriteCheckpoint(restartFile);
		restartFile.close();
	}
	// Unlock();
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

ComputeTool* Simulation::GetTool(const std::string& name,bool attaching)
{
	for (tw::Int i=0;i<computeTool.size();i++)
	{
		if (computeTool[i]->name==name)
		{
			if (attaching)
				computeTool[i]->refCount++;
			return computeTool[i];
		}
	}
	throw tw::FatalError("Could not find tool: " + name);
	return NULL;
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
		tool.push_back(GetTool(word,true));
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
	std::string objectName;

	MetricSpace::ReadCheckpoint(inFile);
	inFile.read((char *)&stepNow,sizeof(tw::Int));
	stepNow++;
	stepsToTake += stepNow-1;
	inFile.read((char *)&elapsedTime,sizeof(tw::Float));
	inFile.read((char *)&signalPosition,sizeof(tw::Float));
	inFile.read((char *)&windowPosition,sizeof(tw::Float));
	inFile.read((char *)&antiSignalPosition,sizeof(tw::Float));
	inFile.read((char *)&antiWindowPosition,sizeof(tw::Float));

	if (uniformDeviate!=NULL)
		uniformDeviate->ReadCheckpoint(inFile);

	if (gaussianDeviate!=NULL)
		gaussianDeviate->ReadCheckpoint(inFile);

	// Read Region objects

	for (tw::Int i=0;i<clippingRegion.size();i++)
	{
		inFile >> objectName;
		inFile.ignore();
		(*tw_out) << "Read checkpoint data for region <" << objectName << ">..." << std::endl;
		Region::FindRegion(clippingRegion,objectName)->ReadCheckpoint(inFile);
	}

	// Read ComputeTool objects

	for (tw::Int i=0;i<computeTool.size();i++)
	{
		inFile >> objectName;
		inFile.ignore();
		(*tw_out) << "Read checkpoint data for tool <" << objectName << ">..." << std::endl;
		GetTool(objectName,false)->ReadCheckpoint(inFile);
	}

	// Read Module objects

	for (tw::Int i=0;i<module.size();i++)
	{
		inFile >> objectName;
		inFile.ignore();
		(*tw_out) << "Read checkpoint data for module <" << objectName << ">..." << std::endl;
		GetModule(objectName)->ReadCheckpoint(inFile);
	}
}

void Simulation::WriteCheckpoint(std::ofstream& outFile)
{
	MetricSpace::WriteCheckpoint(outFile);
	outFile.write((char *)&stepNow,sizeof(tw::Int));
	outFile.write((char *)&elapsedTime,sizeof(tw::Float));
	outFile.write((char *)&signalPosition,sizeof(tw::Float));
	outFile.write((char *)&windowPosition,sizeof(tw::Float));
	outFile.write((char *)&antiSignalPosition,sizeof(tw::Float));
	outFile.write((char *)&antiWindowPosition,sizeof(tw::Float));

	if (uniformDeviate!=NULL)
		uniformDeviate->WriteCheckpoint(outFile);

	if (gaussianDeviate!=NULL)
		gaussianDeviate->WriteCheckpoint(outFile);

	for (auto obj : clippingRegion)
	{
		(*tw_out) << "Checkpointing <" << obj->name << ">" << std::endl;
		obj->WriteCheckpoint(outFile);
	}

	for (auto obj : computeTool)
	{
		(*tw_out) << "Checkpointing <" << obj->name << ">" << std::endl;
		obj->WriteCheckpoint(outFile);
	}

	for (auto obj : module)
	{
		(*tw_out) << "Checkpointing <" << obj->name << ">" << std::endl;
		obj->WriteCheckpoint(outFile);
	}
}


///////////////////////////
//  Read the Input File  //
///////////////////////////


/// The first pass through the input file is used to fully initialize the `Task` and
/// `MetricSpace` parent classes.
void Simulation::InputFileFirstPass()
{
	int numRanksProvided,worldRank;
	MPI_Comm_size(MPI_COMM_WORLD,&numRanksProvided);
	MPI_Comm_rank(MPI_COMM_WORLD,&worldRank);
	// world rank is suitable for reading task data from restart file
	// because this data is the same in every restart file

	try
	{
		// Lock();

		bool foundGrid = false;
		std::stringstream fileName;
		std::stringstream inputString;

		tw::input::PreprocessInputFile(tw::input::FileEnv(inputFileName),inputString);
		AttachUnits(tw::input::GetNativeUnits(inputString.str()),tw::input::GetUnitDensityCGS(inputString.str()));

		inputString.seekg(0);
		outerDirectives.AttachUnits(units);
		outerDirectives.Reset();
		GridReader grid(units);

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
						if (word!="new" && word!="generate")
							inputString >> word;
					} while (word!="new" && word!="generate");
				}
				else
				{
					if (preamble.words[0]=="grid")
					{
						foundGrid = true;
						grid.Read(inputString,preamble);
						grid.UpdateTask(*this);
					}
					if (preamble.words[0]=="warp")
					{
						MangleToolName(preamble.obj_name);
						*tw_out << "Creating Tool <" << preamble.obj_name << ">..." << std::endl;
						// Do not use CreateTool, do not want to increase refCount
						computeTool.push_back(ComputeTool::CreateObjectFromType(preamble.obj_name,tw::tool_type::warp,this,this));
						computeTool.back()->ReadInputFileBlock(inputString);
						computeTool.back()->Initialize(); // OK and necessary to init here
						AttachWarp(dynamic_cast<Warp*>(computeTool.back()));
					}
					if (preamble.words[0]!="grid" && preamble.words[0]!="warp")
						tw::input::ExitInputFileBlock(inputString,true);
				}
			}

			if (com1=="generate")
				tw::input::ExitInputFileBlock(inputString,false);

		} while (!inputString.eof());
		outerDirectives.ThrowErrorIfMissingKeys("Simulation");
		if (!foundGrid)
			throw tw::FatalError("Grid directive was not found.");

		UpdateTimestep(dt0);
		periodic[1] = bc0[1]==tw::bc::par::periodic ? 1 : 0;
		periodic[2] = bc0[2]==tw::bc::par::periodic ? 1 : 0;
		periodic[3] = bc0[3]==tw::bc::par::periodic ? 1 : 0;

		// Check integer viability
		int64_t totalCellsPerRank = int64_t(globalCells[1])*int64_t(globalCells[2])*int64_t(globalCells[3])/int64_t(numRanksProvided);
		if (totalCellsPerRank>=pow(2,31) && sizeof(tw::Int)==4)
			throw tw::FatalError("You must recompile turboWAVE with 64 bit integers to handle this many grid cells.");

		// Verify and (if necessary) correct decomposition
		if (NumTasks() != numRanksProvided)
		{
			*tw_out << "WARNING: Bad decomposition ";
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
			*tw_out << "(defaulting to " << domains[1] << "x" << domains[2] << "x" << domains[3] << ")" << std::endl;
		}
		*tw_out << NumTasks() << "-Way Decomposition" << std::endl;

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
		Resize(*this,grid.GlobalCorner(),grid.GlobalSize(),2,grid.Geometry());

		// Random numbers
		uniformDeviate = new UniformDeviate(1 + strip[0].Get_rank()*(MaxSeed()/numRanksProvided));
		gaussianDeviate = new GaussianDeviate(1 + strip[0].Get_rank()*(MaxSeed()/numRanksProvided) + MaxSeed()/(2*numRanksProvided));

		// Unlock();
	}

	catch (tw::FatalError& e)
	{
		if (worldRank==0)
		{
			std::cerr << "FATAL ERROR: " << e.what() << std::endl;
			std::cerr << "Could not start simulation --- exiting now." << std::endl;
		}
		// Following acts as a barrier before exiting; we should implement MPI_Barrier in TW_MPI.
		char buf[255];
		MPI_Bcast(buf,1,MPI_BYTE,0,MPI_COMM_WORLD);
		exit(1);
	}
}

void Simulation::NestedDeclaration(const std::string& com,std::stringstream& inputString,Module *super)
{
	// To be called by supermodules that want to add submodules or tools while
	// reading their own input file block.
	// This function can be called recursively.

	// Get the preamble = words that come between "new" and the opening brace
	tw::input::Preamble preamble = tw::input::EnterInputFileBlock(com,inputString,"{=");
	if (preamble.attaching)
		throw tw::FatalError(preamble.err_prefix+"keyword <for> is not allowed in a nested declaration.");

	tw::module_type whichModule = Module::CreateTypeFromInput(preamble);
	tw::tool_type whichTool = ComputeTool::CreateTypeFromInput(preamble);
	if (whichModule==tw::module_type::none && whichTool==tw::tool_type::none)
		throw tw::FatalError(preamble.err_prefix+"key was not recognized.");
	if (whichModule!=tw::module_type::none && whichTool!=tw::tool_type::none)
		throw tw::FatalError(preamble.err_prefix+"key claimed by both Module and Tool, this is a bug in the code.");

	if (whichModule!=tw::module_type::none)
	{
		if (Module::SingularType(whichModule))
			if (module_map.find(whichModule)!=module_map.end())
				throw tw::FatalError(preamble.err_prefix+"Singular module type was created twice.  Check order of input file.");
		MangleModuleName(preamble.obj_name);
		(*tw_out) << "   Attaching nested module <" << preamble.obj_name << ">..." << std::endl;
		Module *sub = Module::CreateObjectFromType(preamble.obj_name,whichModule,this);
		module.push_back(sub);
		module_map[whichModule] = sub;
		// The following may lead to a recursive call of this function.
		// If so the module vector and map can be modified.
		sub->ReadInputFileBlock(inputString);
		super->AddSubmodule(sub);
	}

	if (whichTool!=tw::tool_type::none)
	{
		ComputeTool *tool = CreateTool(preamble.obj_name,whichTool);
		(*tw_out) << "   Attaching nested tool <" << tool->name << ">..." << std::endl;
		tool->ReadInputFileBlock(inputString);
		super->moduleTool.push_back(tool);
	}
}

Module* Simulation::RecursiveAutoSuper(tw::module_type reqType,const std::string& basename)
{
	// If it already exists we are done
	if (module_map.find(reqType)!=module_map.end())
			return module_map[reqType];

	// Automatic creation of supermodule
	if (!Module::AutoModuleType(reqType))
		throw tw::FatalError("Module <"+basename+"> requires a supermodule that cannot be created automatically.");
	std::string super_module_name = basename + "_sup";
	MangleModuleName(super_module_name);
	(*tw_out) << "Installing supermodule triggered by <" << basename << ">..." << std::endl;
	Module *super = Module::CreateObjectFromType(super_module_name,reqType,this);
	module.push_back(super);
	module_map[reqType] = super;

	// Handle recursion
	tw::module_type superSuperType = Module::RequiredSupermoduleType(reqType);
	if (superSuperType!=tw::module_type::none)
	{
		Module *superSuper = RecursiveAutoSuper(superSuperType,super_module_name);
		superSuper->AddSubmodule(super);
	}
	return super;
}

void Simulation::ReadInputFile()
{
	// Lock();

	std::string com1,word;
	Profile* theProfile;
	std::stringstream inputString;

	(*tw_out) << std::endl << "Reading Input File..." << std::endl << std::endl;

	tw::input::PreprocessInputFile(tw::input::FileEnv(inputFileName),inputString);

	inputString.seekg(0);
	outerDirectives.Reset();

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

			// Intercept items processed during the first pass
			if (preamble.words[0]=="grid")
			{
				preamble.words[0] = "tw::none";
				processed = true;
				tw::input::ExitInputFileBlock(inputString,true);
			}
			if (preamble.words[0]=="warp")
			{
				preamble.words[0] = "tw::none";
				processed = true;
				tw::input::ExitInputFileBlock(inputString,true);
			}

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

			// Module Installation
			tw::module_type whichModule = Module::CreateTypeFromInput(preamble);
			if (whichModule!=tw::module_type::none)
			{
				processed = true;
				if (Module::SingularType(whichModule))
					if (module_map.find(whichModule)!=module_map.end())
						throw tw::FatalError(preamble.err_prefix + "singular module type was created twice.  Check order of input file.");
				MangleModuleName(preamble.obj_name);
				(*tw_out) << "Installing module <" << preamble.obj_name << ">..." << std::endl;
				Module *sub = Module::CreateObjectFromType(preamble.obj_name,whichModule,this);
				module.push_back(sub);
				module_map[whichModule] = sub;
				sub->ReadInputFileBlock(inputString); // important to note this can change module vector and map if there are nested declarations
				if (preamble.attaching)
				{
					Module *super = GetModule(preamble.owner_name);
					(*tw_out) << "Attaching <" << preamble.obj_name << "> to <" << preamble.owner_name << ">..." << std::endl;
					super->AddSubmodule(sub);
				}
				else
				{
					// If not explicitly attaching, but supermodule is required, find or create one
					tw::module_type reqType = Module::RequiredSupermoduleType(whichModule);
					if (reqType!=tw::module_type::none)
					{
						Module *super = RecursiveAutoSuper(reqType,preamble.obj_name);
						super->AddSubmodule(sub);
					}
				}
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
				clippingRegion.back()->directives.AttachUnits(units);
				clippingRegion.back()->name = preamble.obj_name;
				clippingRegion.back()->ReadInputFileBlock(inputString);
			}

			if (!processed)
				throw tw::FatalError(preamble.err_prefix+"keys were not understood.");
		}

	} while (!inputString.eof());

	outerDirectives.ThrowErrorIfMissingKeys("Simulation");

	// Unlock();
}

void Simulation::Diagnose()
{
	std::vector<Diagnostic*> diagnostic;
	for (auto tool : computeTool)
		if (dynamic_cast<Diagnostic*>(tool))
			diagnostic.push_back(dynamic_cast<Diagnostic*>(tool));

	// DIAGNOSTIC PREP
	// Main purpose of this step is to let modules transfer data from compute devices.
	// This has a high cost, so only do if necessary.

	bool doing_diagnostics=false;
	for (auto d : diagnostic)
		doing_diagnostics = doing_diagnostics || d->WriteThisStep(elapsedTime,dt,stepNow);
	if (doing_diagnostics)
	{
		for (auto m : module)
			m->StartDiagnostics();
	}
	else
	{
		return;
	}

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
