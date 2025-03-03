module;

#include <algorithm>
#include <tree_sitter/api.h>
#include "tw_includes.h"

module twmodule;
import base;
import input;
import factory;

void Simulation::UpdateTimestep(tw::Float dt0)
{
	SetupTimeInfo(dt0);
	for (auto m : module)
		m->SetupTimeInfo(dt0);
}

GridReader::GridReader(tw::UnitConverter& uc)
{
	adaptiveGrid = adaptiveTimestep = found = false;
	geo = tw::grid::cartesian;
	absoluteRef = 0.0;

	directives.Add("origin",new tw::input::Numbers<tw::Float>(&relativeRef[0],4));
	directives.Add("shift",new tw::input::Numbers<tw::Float>(&absoluteRef[0],4),false);
	directives.Add("cell size",new tw::input::Numbers<tw::Float>(&spacing[0],4));
	directives.Add("adaptive timestep",new tw::input::Bool(&adaptiveTimestep),false);
	directives.Add("adaptive grid",new tw::input::Bool(&adaptiveGrid),false);
	directives.Add("dimensions",new tw::input::Numbers<tw::Int>(&req_dim[0],4));
	directives.Add("decomposition",new tw::input::Numbers<tw::Int>(&req_dom[0],4));
	std::map<std::string,tw::grid::geometry> geo_map = {{"cartesian",tw::grid::cartesian},{"cylindrical",tw::grid::cylindrical},{"spherical",tw::grid::spherical}};
	directives.Add("geometry",new tw::input::Enums<tw::grid::geometry>(geo_map,&geo),false);
	directives.AttachUnits(uc);
}

bool GridReader::Read(const TSTreeCursor *curs0,const std::string& src)
{
	TSTreeCursor curs = ts_tree_cursor_copy(curs0);
	if (tw::input::node_kind(&curs)=="new") {
		ts_tree_cursor_goto_first_child(&curs); // `new` token
		ts_tree_cursor_goto_next_sibling(&curs);
		if (tw::input::node_text(&curs,src)=="grid") {
			ts_tree_cursor_goto_next_sibling(&curs);
			if (tw::input::node_kind(&curs)=="block") {
				ts_tree_cursor_goto_first_child(&curs);
				tw::input::next_named_node(&curs,true);
				do {
					if (!directives.ReadNext(&curs,src)) {
						tw::input::ThrowParsingError(&curs,src,"only assignments allowed in grid block");
					}
				} while (tw::input::next_named_node(&curs,false));
				directives.ThrowErrorIfMissingKeys("grid");
				if (req_dom[0]!=1)
					tw::input::ThrowParsingError(&curs,src,"Time decomposition must be 1");
				found = true;
				return true;
			} else {
				tw::input::ThrowParsingError(&curs,src,"Ill-formed grid block");
			}
		}
	}
	return false;
}

tw::vec4 GridReader::GlobalSize()
{
	return tw::vec4(
		spacing[0] * req_dim[0],
		spacing[1] * req_dim[1],
		spacing[2] * req_dim[2],
		spacing[3] * req_dim[3]
	);
}

tw::vec4 GridReader::GlobalCorner()
{
	return absoluteRef - relativeRef*GlobalSize();
}

void GridReader::UpdateTask(Task& tsk)
{
	for (tw::Int i=0;i<4;i++)
	{
		tsk.globalCells[i] = req_dim[i];
		tsk.domains[i] = req_dom[i];
	}
}

void GridReader::UpdateSpace(MetricSpace& ms)
{
	// whatever not handled by Resize method
	ms.adaptiveGrid = adaptiveGrid;
	ms.adaptiveTimestep = adaptiveTimestep;
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

	stepNow = 0;
	lastTime = 0;
	dumpPeriod = 0;
	inputFilePass = 0;

	bc0[1] = tw::bc::par::periodic;
	bc1[1] = tw::bc::par::periodic;
	bc0[2] = tw::bc::par::periodic;
	bc1[2] = tw::bc::par::periodic;
	bc0[3] = tw::bc::par::absorbing;
	bc1[3] = tw::bc::par::absorbing;

	outerDirectives.Add("native units",new tw::input::Enums<tw::units>(tw::get_unit_map(),&nativeUnits),false);
	outerDirectives.Add("unit density",new tw::input::Float(&unitDensityCGS));
	outerDirectives.Add("xboundary",new tw::input::Enums<tw::bc::par>(tw::bc::par_map(),&bc0[1],&bc1[1]));
	outerDirectives.Add("yboundary",new tw::input::Enums<tw::bc::par>(tw::bc::par_map(),&bc0[2],&bc1[2]));
	outerDirectives.Add("zboundary",new tw::input::Enums<tw::bc::par>(tw::bc::par_map(),&bc0[3],&bc1[3]));
	outerDirectives.Add("dtmin",new tw::input::Float(&dtMin),false);
	outerDirectives.Add("dtmax",new tw::input::Float(&dtMax),false);
	outerDirectives.Add("dtcrit",new tw::input::Float(&dtCritical),false);
	outerDirectives.Add("maxtime",new tw::input::Float(&elapsedTimeMax),false);
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

	if (dynamic_cast<std::ofstream*>(tw_err))
		((std::ofstream*)tw_err)->close();
	if (dynamic_cast<std::stringstream*>(tw_err))
		delete tw_err;
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

	*tw_out << std::endl << term::green << term::bold << "Starting turboWAVE Session" << term::reset_all << std::endl << std::endl;
	*tw_out << "Floating point precision = " << sizeof(tw::Float)*8 << " bits" << std::endl;
	*tw_out << "Maximum OpenMP threads = " << omp_get_max_threads() << std::endl;
}

void Simulation::Run()
{
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
		(*tw_out) << std::endl << term::green << term::bold << "Prepare Simulation" << term::reset_all << std::endl << std::endl;

		PrepareSimulation();

		(*tw_out) << std::endl << term::green << term::bold << "Begin Simulation" << term::reset_all << std::endl << std::endl;

		tw::Int startTime = GetSeconds();
		lastTime = startTime;

		if (GetSeconds()<0)
		{
			(*tw_out) << std::endl << term::warning << ": System clock is not responding properly." << std::endl << std::endl;
		}

		(*tw_out) << "Current status can be viewed in 'twstat' file." << std::endl;
		(*tw_out) << "This executable does not support interactive commands." << std::endl << std::endl;

		while (stepNow <= dim[0] && elapsedTime < elapsedTimeMax)
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
		if (strip[0].Get_rank()==0)
		{
			twstat.open("twstat");
			twstat << "The simulation failed. For more info see stdout." << std::endl;
			twstat.close();
		}
		exit(1);
	}

	*tw_out << std::endl << term::green << term::bold << "TurboWAVE Session Completed" << term::reset_all << std::endl;
}

void Simulation::Test()
{
	outputLevel = 1;
	SetupIO();

	// Get basic MPI data
	int numRanksProvided,worldRank;
	MPI_Comm_size(MPI_COMM_WORLD,&numRanksProvided);
	MPI_Comm_rank(MPI_COMM_WORLD,&worldRank);

	// Random numbers
	uniformDeviate = new UniformDeviate(1 + strip[0].Get_rank()*(MaxSeed()/numRanksProvided));
	gaussianDeviate = new GaussianDeviate(1 + strip[0].Get_rank()*(MaxSeed()/numRanksProvided) + MaxSeed()/(2*numRanksProvided));

	// diverting output
	std::stringstream test_out;
	std::stringstream test_report;
	std::ostream *alt_out = new std::stringstream;
	std::ostream *save_out = tw_out;

	// each test must setup its own task and grid
	if (numRanksProvided!=2)
		throw tw::FatalError("test must use -n 2");
	AttachUnits(tw::units::plasma,1e19);

	tw::Float failed = 0;
	bool tested = false;
	success_count = 0;
	failure_count = 0;
	ComputeTool *tool;
	Module *module;

	// Testables that are not ComputeTool or Module
	// TODO: we need to unify handling of all testables.

	*tw_out << std::endl << "testing " << term::bold << term::cyan << "metric space" << term::reset_all << std::endl;
 	Initialize(tw::idx4(1,1,1,2).array,tw::idx4(1,4,1,4).array,tw::idx4(0,1,1,0).array);
	Resize(*this,tw::vec4(0,0,0,0),tw::vec4(0.1,0.8,0.2,0.8),2);
	tw::Int testId = 1;
	failed = 0;
	try {
		tested = MetricSpace::Test(testId);
	} catch (tw::FatalError& e) {
		failed = 1;
		test_report << failure_count+1 << ". metric space" << std::endl;
		test_report << e.what() << std::endl;
	}
	strip[0].AllSum(&failed,&failed,sizeof(tw::Float),0);
	if (tested && failed == 0) {
		success_count++;
		test_out << "    " << term::ok << " " << term::green << "metric space" << term::reset_all << std::endl;
	} else if (failed > 0) {
		failure_count++;
		test_out << "    " << term::err << " " << term::red << "metric space" << term::reset_all << std::endl;
	}
	*tw_out << test_out.str();
	test_out.str("");
	test_out.clear();

	// Test ComputeTools

	for (auto m : ComputeTool::Map())
	{
		if (unitTest==m.first || unitTest=="--all")
		{
			*tw_out << std::endl << "testing " << term::bold << term::cyan << m.first << term::reset_all << std::endl;
			tw::Int gridId = 1;
			while (ComputeTool::SetTestGrid(m.second,gridId,this,this))
			{
				std::string gridStr("    "+std::to_string(this->globalCells[1])+"x"+std::to_string(this->globalCells[2])+"x"+std::to_string(this->globalCells[3])+" grid");
				tw_out = alt_out;
				tw::Int testId = 1;
				do
				{
					try
					{
						tool = CreateTool("test_tool",m.second);
					}
					catch (tw::FatalError& e)
					{
						test_out << gridStr + "    " << term::yellow << "tool rejected the environment" << term::reset_all << std::endl;
						tool = NULL;
					}
					if (tool!=NULL)
					{
						failed = 0;
						try {
							tested = tool->Test(testId);
						} catch(tw::FatalError& e) {
							failed = 1;
							test_report << failure_count+1 << ". " << m.first << " grid " << gridId << std::endl;
							test_report << e.what() << std::endl;
						}
						strip[0].AllSum(&failed,&failed,sizeof(tw::Float),0);
						if (tested && failed == 0) {
							success_count++;
							test_out << gridStr << "    " << term::ok << " " << term::green << tool->testName << term::reset_all << std::endl;
						} else if (failed > 0) {
							failure_count++;
							test_out << gridStr << "    " << term::err << " " << term::red << tool->testName <<  term::reset_all << std::endl;
						}
						RemoveTool(tool);
					}
				} while (testId>1);
				MPI_Barrier(MPI_COMM_WORLD);
				tw_out = save_out;
				*tw_out << test_out.str();
				test_out.str("");
				test_out.clear();
				gridId++;
			}
		}
	}

	// Test Modules

	for (auto m : Module::Map())
	{
		if (unitTest==m.first || unitTest=="--all")
		{
			*tw_out << std::endl << "testing " << term::bold << term::cyan << m.first << term::reset_all << std::endl;
			tw::Int gridId = 1;
			while (Module::SetTestGrid(m.second,gridId,this))
			{
				std::string gridStr("    "+std::to_string(this->globalCells[1])+"x"+std::to_string(this->globalCells[2])+"x"+std::to_string(this->globalCells[3])+" grid");
				tw_out = alt_out;
				tw::Int testId = 1;
				do
				{
					try
					{
	 					module = factory::CreateModuleFromType("test_module",m.second,this);
					}
					catch (tw::FatalError& e)
					{
						test_out << "    " << term::yellow << "module rejected the environment" << term::reset_all << std::endl;
						module = NULL;
					}
					if (module!=NULL)
					{
						failed = 0;
						try {
							tested = module->Test(testId);
						} catch(tw::FatalError& e) {
							failed = 1;
							test_report << ++failure_count << ". " << m.first << " grid " << gridId << std::endl;
							test_report << e.what() << std::endl;
						}
						strip[0].AllSum(&failed,&failed,sizeof(tw::Float),0);
						if (tested && failed == 0) {
							success_count++;
							test_out << gridStr << "    " << term::ok << " " << term::green << module->testName << term::reset_all << std::endl;
						} else if (failed > 0) {
							failure_count++;
							test_out << gridStr << "    " << term::err << " " << term::red << module->testName <<  term::reset_all << std::endl;
						}
						delete module;
					}
				} while (testId>1);
				MPI_Barrier(MPI_COMM_WORLD);
				tw_out = save_out;
				*tw_out << test_out.str();
				test_out.str("");
				test_out.clear();
				gridId++;
			}
		}
	}

	*tw_out << std::endl;
	if (test_report.str().size()>0) {
		*tw_err << term::bold << term::red << "Unit Tests Failing" << term::reset_all << std::endl << std::endl << test_report.str();
	} else if (failure_count > 0) {
		*tw_err << term::bold << term::red << "Unit Tests Failing" << term::reset_all << std::endl;
		*tw_err << "Details may be in node-scoped output files" << std::endl;
	} else {
		*tw_out << term::bold << term::green << "Unit Tests Passing" << term::reset_all << " - " << success_count << " succeeded, " << failure_count << " failed" << std::endl;
	}
	delete alt_out;
	MPI_Barrier(MPI_COMM_WORLD);
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
		*theStream << "Current step size: " << spacing[0] << std::endl;
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
		*theStream << "Steps to take: " << dim[0] << std::endl;
		*theStream << "Steps remaining: " << dim[0] - stepNow << std::endl;
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

	// TODO: we don't advance the time corner because the particle cell is being
	// encoded with the time level as well.  It might make more sense to hold the
	// particle's index at 0 and advance the corner instead.  But then we might need
	// to update the time box for everything that inherits DiscreteSpace, such as Module.
	// There is further confusion in the fact that the dimension of the time box is
	// steps to take, not steps stored in memory.
	//corner[0] += spacing[0];

	elapsedTime += spacing[0];
	signalPosition += signalSpeed*spacing[0];
	antiSignalPosition -= signalSpeed*spacing[0];
	stepNow++;

	if (adaptiveGrid)
		for (auto m : module)
			m->AdaptGrid();

	if (movingWindow && signalPosition>=(windowPosition + spacing[3]) && dim[3]>1)
		MoveWindow();

	if (!movingWindow && antiSignalPosition<=(antiWindowPosition - spacing[3]) && dim[3]>1)
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
	computeTool.push_back(factory::CreateToolFromType(name,theType,this,this));
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

/// @brief handle retrieval of named tools
/// @param tool tool list owned by a module
/// @param curs should be on a `get` node
/// @param src source document
void Simulation::ToolFromDirective(std::vector<ComputeTool*>& tool,TSTreeCursor *curs,const std::string& src)
{
	auto word = tw::input::next_named_node_text(curs,src);
	tw::input::StripQuotes(word);
	if (CheckModule(word))
		tw::input::ThrowParsingError(curs,src,"Tried to <get> module, but <get> can only be used for tools.");
	tool.push_back(GetTool(word,true));
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
	windowPosition += spacing[3];
	corner[3] += spacing[3];
	globalCorner[3] += spacing[3];
	
	for (i=lfg[3];i<=ufg[3];i++)
		X(i,3) += spacing[3];

	for (i=0;i<clippingRegion.size();i++)
		if (clippingRegion[i]->moveWithWindow)
			clippingRegion[i]->Translate(tw::vec3(0,0,spacing[3]));
		else
			clippingRegion[i]->Initialize(*this,this);

	for (i=0;i<module.size();i++)
		module[i]->MoveWindow();
}

void Simulation::AntiMoveWindow()
{
	tw::Int i;
	antiWindowPosition -= spacing[3];

	for (i=0;i<module.size();i++)
		module[i]->AntiMoveWindow();
}

void Simulation::ReadCheckpoint(std::ifstream& inFile)
{
	std::string objectName;

	MetricSpace::ReadCheckpoint(inFile);
	inFile.read((char *)&stepNow,sizeof(tw::Int));
	stepNow++;
	dim[0] += stepNow-1;
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

/// @brief main visitor callback for AST walker
/// @param curs cursor, typically passed from walker
/// @param inputFilePass 0=syntax, 1=grid, 2=objects
/// @return action for the walker to take upon exit
tw::input::navigation Simulation::visit(TSTreeCursor *curs) {

	if (inputFilePass == 0) {
		TSNode node = ts_tree_cursor_current_node(curs);
		if (ts_node_has_error(node)) {
			return tw::input::navigation::gotoChild;
		} else if (ts_node_is_error(node)) {
			throw tw::FatalError("syntax error at " + tw::input::loc_str(curs));
		}
		return tw::input::navigation::gotoSibling;


	} else if (inputFilePass == 1) {
		if (tw::input::node_kind(curs) == "input_file") {
			return tw::input::navigation::gotoChild;
		} else if (tw::input::node_kind(curs) == "assignment") {
			// Process outer assignments
			outerDirectives.ReadNext(curs,src);
			return tw::input::navigation::gotoSibling;
		} else if (tw::input::node_kind(curs) == "new") {
			if (!outerDirectives.TestKey("native units")) {
				throw tw::FatalError("`new` encountered before `native units`");
			}
			if (!outerDirectives.TestKey("unit density")) {
				throw tw::FatalError("`new` encountered before `unit density`");
			}
			ts_tree_cursor_goto_first_child(curs);
			ts_tree_cursor_goto_next_sibling(curs);
			std::string key = tw::input::node_text(curs,src);
			ts_tree_cursor_goto_parent(curs);
			if (key == "grid") {
				gridReader->Read(curs,src);
				gridReader->UpdateTask(*this);
				gridReader->UpdateSpace(*this);
				return tw::input::navigation::gotoSibling;
			} else if (key == "warp") {
				TSTreeCursor saveCurs = ts_tree_cursor_copy(curs);
				tw::input::Preamble preamble = tw::input::GetPreamble(curs,src);
				MangleToolName(preamble.obj_name);
				*tw_out << "Creating Tool <" << preamble.obj_name << ">..." << std::endl;
				// Do not use CreateTool, do not want to increase refCount
				computeTool.push_back(factory::CreateToolFromType(preamble.obj_name,tw::tool_type::warp,this,this));
				computeTool.back()->ReadInputFileBlock(curs,src);
				computeTool.back()->Initialize(); // OK and necessary to init here
				AttachWarp(dynamic_cast<Warp*>(computeTool.back()));
				*curs = ts_tree_cursor_copy(&saveCurs);
				return tw::input::navigation::gotoSibling;
			} else {
				// not a grid or warp, so skip it
				return tw::input::navigation::gotoSibling;
			}
		} else {
			// not a basic `new` directive so skip it
			return tw::input::navigation::gotoSibling;
		}


	} else if (inputFilePass == 2) {
		TSNode node = ts_tree_cursor_current_node(curs);
		std::string typ = ts_node_type(node);
		if (tw::input::node_kind(curs) == "input_file") {
			return tw::input::navigation::gotoChild;
		} else if (typ == "comment") {
			return tw::input::navigation::gotoSibling;
		} else if (typ == "assignment") {
			// Process outer directives again, after first pass.
			// This is redundant, but helps minimize keyword collisions with object names.
			outerDirectives.ReadNext(curs,src);
			return tw::input::navigation::gotoSibling;
		} else if (typ == "new" || typ == "associative_new" || typ == "generate") {
			tw::input::Preamble preamble = tw::input::GetPreamble(curs,src);

			// Intercept items processed during the first pass
			if (preamble.obj_key=="grid" || preamble.obj_key=="warp")
				return tw::input::navigation::gotoParentSibling;

			// Install a pre or post declared tool
			tw::tool_type whichTool = ComputeTool::CreateTypeFromInput(preamble);
			if (whichTool!=tw::tool_type::none)
			{
				if (preamble.attaching)
				{
					ComputeTool *tool = CreateTool(preamble.obj_name,whichTool);
					(*tw_out) << "Attaching <" << tool->name << "> to <" << preamble.owner_name << ">..." << std::endl;
					tool->ReadInputFileBlock(curs,src);
					GetModule(preamble.owner_name)->moduleTool.push_back(tool);
				}
				else
				{
					bool duplicate = MangleToolName(preamble.obj_name);
					(*tw_out) << "Creating Tool <" << preamble.obj_name << ">..." << std::endl;
					if (duplicate && errorCheckingLevel>0)
						tw::input::ThrowParsingError(curs,src,"duplicate tool name.");
					// Do not use CreateTool, do not want to increase refCount
					computeTool.push_back(factory::CreateToolFromType(preamble.obj_name,whichTool,this,this));
					computeTool.back()->ReadInputFileBlock(curs,src);
				}
				return tw::input::navigation::gotoSibling;
			}

			// Module Installation
			tw::module_type whichModule = Module::CreateTypeFromInput(preamble);
			if (whichModule!=tw::module_type::none)
			{
				if (Module::SingularType(whichModule))
					if (module_map.find(whichModule)!=module_map.end())
						tw::input::ThrowParsingError(curs,src,"singular module type was created twice.  Check order of input file.");
				MangleModuleName(preamble.obj_name);
				(*tw_out) << "Installing module <" << preamble.obj_name << ">..." << std::endl;
				Module *sub = factory::CreateModuleFromType(preamble.obj_name,whichModule,this);
				module.push_back(sub);
				module_map[whichModule] = sub;
				sub->ReadInputFileBlock(curs,src); // important to note this can change module vector and map if there are nested declarations
				if (preamble.attaching && preamble.owner_name.size() > 0)
				{
					Module *super = GetModule(preamble.owner_name);
					(*tw_out) << "Attaching <" << preamble.obj_name << "> to <" << preamble.owner_name << ">..." << std::endl;
					super->AddSubmodule(sub);
				}
				else if (preamble.attaching && preamble.owner_name.size() == 0)
				{
					// If not explicitly attaching, but supermodule is required, find or create one
					tw::module_type reqType = Module::RequiredSupermoduleType(whichModule);
					if (reqType!=tw::module_type::none)
					{
						Module *super = RecursiveAutoSuper(reqType,preamble.obj_name);
						super->AddSubmodule(sub);
					}
				}
				return tw::input::navigation::gotoSibling;
			}

			// Regions are neither modules nor tools
			if (preamble.obj_key.substr(0,6)=="region")
			{
				std::string rgnType = tw::input::trim(preamble.obj_key.substr(6));
				clippingRegion.push_back(Region::CreateObjectFromString(clippingRegion,rgnType));
				clippingRegion.back()->directives.AttachUnits(units);
				clippingRegion.back()->name = preamble.obj_name;
				clippingRegion.back()->ReadInputFileBlock(curs,src);
				return tw::input::navigation::gotoSibling;
			}
			
			throw tw::FatalError("unknown key <" + preamble.obj_key + "> at " + tw::input::loc_str(curs));

		} else if (typ == "reaction" || typ == "collision" || typ == "excitation") {
			for (tw::Int i=0;i<module.size();i++)
				module[i]->ReadQuasitoolBlock(curs,src);
			return tw::input::navigation::gotoSibling;
		} else {
			tw::input::ThrowParsingError(curs,src,"unhandled directive");
			return tw::input::navigation::exit; // compiler doesn't know above line throws
		}


	} else {
		throw tw::FatalError("invalid number of input file passes");
	}
}

tw::input::navigation Simulation::descend(TSTreeCursor *curs) {
	return tw::input::navigation::exit;
}

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

		std::stringstream fileName;
		tw::input::FileEnv fenv(inputFileName);
		fenv.OpenDeck(src);
		TSTree *tree = tw::input::GetTree(src);

		// TODO: handle preprocessing
		//tw::input::PreprocessInputFile(tw::input::FileEnv(inputFileName),inputString);
		AttachUnits(tw::input::GetNativeUnits(tree,src),tw::input::GetUnitDensityCGS(tree,src));

		outerDirectives.AttachUnits(units);
		outerDirectives.Reset();
		gridReader = new GridReader(units);

		inputFilePass = 0;
		tw::input::WalkTree(tree,this);
		inputFilePass = 1;
		tw::input::WalkTree(tree,this);

		outerDirectives.ThrowErrorIfMissingKeys("Simulation");
		if (!gridReader->FoundGrid())
			throw tw::FatalError("Grid directive was not found.");

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
			*tw_out << term::warning << ": Bad decomposition ";
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
		Resize(*this,gridReader->GlobalCorner(),gridReader->GlobalSize(),2,gridReader->Geometry());
		delete gridReader;

		// Random numbers
		uniformDeviate = new UniformDeviate(1 + strip[0].Get_rank()*(MaxSeed()/numRanksProvided));
		gaussianDeviate = new GaussianDeviate(1 + strip[0].Get_rank()*(MaxSeed()/numRanksProvided) + MaxSeed()/(2*numRanksProvided));

		// Unlock();
	}

	catch (tw::FatalError& e)
	{
		if (worldRank==0)
		{
			std::cerr << term::red << "FATAL ERROR: " << term::reset_all << e.what() << std::endl;
			std::cerr << "Could not start simulation --- exiting now." << std::endl;
		}
		MPI_Barrier(MPI_COMM_WORLD);
		exit(1);
	}
}

/// @brief to be called by supermodules that want to add submodules or tools, can be recursive
/// @param curs should be on the outer directive node
/// @param src source document
/// @param super module that is adding the item
void Simulation::NestedDeclaration(TSTreeCursor *curs,const std::string& src,Module *super)
{
	tw::input::Preamble preamble = tw::input::GetPreamble(curs,src);
	if (preamble.attaching)
		tw::input::ThrowParsingError(curs,src,"keyword <for> is not allowed in a nested declaration.");

	tw::module_type whichModule = Module::CreateTypeFromInput(preamble);
	tw::tool_type whichTool = ComputeTool::CreateTypeFromInput(preamble);
	if (whichModule==tw::module_type::none && whichTool==tw::tool_type::none)
		tw::input::ThrowParsingError(curs,src,"key was not recognized.");
	if (whichModule!=tw::module_type::none && whichTool!=tw::tool_type::none)
		tw::input::ThrowParsingError(curs,src,"key claimed by both Module and Tool, this is a bug in the code.");

	if (whichModule!=tw::module_type::none)
	{
		if (Module::SingularType(whichModule))
			if (module_map.find(whichModule)!=module_map.end())
				tw::input::ThrowParsingError(curs,src,"Singular module type was created twice.  Check order of input file.");
		MangleModuleName(preamble.obj_name);
		(*tw_out) << "   Attaching nested module <" << preamble.obj_name << ">..." << std::endl;
		Module *sub = factory::CreateModuleFromType(preamble.obj_name,whichModule,this);
		module.push_back(sub);
		module_map[whichModule] = sub;
		// The following may lead to a recursive call of this function.
		// If so the module vector and map can be modified.
		sub->ReadInputFileBlock(curs,src);
		super->AddSubmodule(sub);
	}

	if (whichTool!=tw::tool_type::none)
	{
		ComputeTool *tool = CreateTool(preamble.obj_name,whichTool);
		(*tw_out) << "   Attaching nested tool <" << tool->name << ">..." << std::endl;
		tool->ReadInputFileBlock(curs,src);
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
	Module *super = factory::CreateModuleFromType(super_module_name,reqType,this);
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

	(*tw_out) << "Reading Input File..." << std::endl << std::endl;

	// TODO: handle preprocessing
	//tw::input::PreprocessInputFile(tw::input::FileEnv(inputFileName),inputString);

	outerDirectives.Reset();
	TSTree *tree = tw::input::GetTree(src);
	inputFilePass = 2;
	tw::input::WalkTree(tree,this);
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
		doing_diagnostics = doing_diagnostics || d->WriteThisStep(elapsedTime,spacing[0],stepNow);
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
		if (d->WriteThisStep(elapsedTime,spacing[0],stepNow))
		{
			d->Start();
			d->ReportNumber("time",elapsedTime,true);
			d->ReportNumber("dt",spacing[0],true);
			d->ReportNumber("zwindow",windowPosition,true);
			for (auto m : module)
				if (has_diagnostic(m,d))
					m->Report(*d);
			d->Finish();
		}
	}
}
