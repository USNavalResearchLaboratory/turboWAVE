module;

#include "tw_includes.h"
#include "tw_logger.h"
#ifndef USE_STD_MODULE
    #include <algorithm>
#endif

module twmodule;
import base;
import logger;
import input;
import factory;

Simulation::Simulation(const bool interactive,
	const std::string& test_name,
	const std::string& file_name,
	const std::string& restart_name,
	const std::string& platform,
	const std::string& device,
	const tw::Int& outputLevel,
	const tw::Int& errorCheckingLevel)
{
	this->interactive = interactive;
	unitTest = test_name;
	inputFileName = file_name;
	restartFileName = restart_name;
	platformSearchString = platform;
	deviceSearchString = device;
	this->outputLevel = outputLevel;
	this->errorCheckingLevel = errorCheckingLevel;
	clippingRegion.push_back(new EntireRegion(clippingRegion));

	nativeUnits = tw::units::plasma;

	neutralize = true;
	movingWindow = false;
	maxWindowPosition[0] = tw::big_pos;

	stepNow = 0;
	previous_timestamp = 0;
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
	outerDirectives.Add("dtmin",new tw::input::Float(&min_spacing[0]),false);
	outerDirectives.Add("dtmax",new tw::input::Float(&max_spacing[0]),false);
	outerDirectives.Add("dtcrit",new tw::input::Float(&critical_spacing[0]),false);
	outerDirectives.Add("maxtime",new tw::input::Float(&maxWindowPosition[0]),false);
	outerDirectives.Add("dump period",new tw::input::Int(&dumpPeriod),false);
	outerDirectives.Add("neutralize",new tw::input::Bool(&neutralize),false);
	outerDirectives.Add("window speed",new tw::input::Float(&solutionVelocity[3]),false);
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
}

void Simulation::UpdateTimestep(tw::Float dt0)
{
	DynSpace::SetupTimeInfo(dt0);
}

void Simulation::SetupIO()
{
	int worldRank;
	MPI_Comm_rank(MPI_COMM_WORLD,&worldRank);
	// Set up standard outputs
	if (worldRank!=0)
	{
		if (outputLevel>0)
		{
			auto tw_out = new std::ofstream(std::to_string(worldRank) + "_stdout.txt");
			std::cout.rdbuf(tw_out->rdbuf());
			auto tw_err = new std::ofstream(std::to_string(worldRank) + "_stderr.txt");
			std::cerr.rdbuf(tw_err->rdbuf());
		}
		else
		{
			// put outputs into throw-away strings
			auto tw_out = new std::stringstream;
			std::cout.rdbuf(tw_out->rdbuf());
			auto tw_err = new std::stringstream;
			std::cerr.rdbuf(tw_err->rdbuf());
		}
	}

	std::println(std::cout,"\n{}{}Starting turboWAVE Session{}\n",term::green,term::bold,term::reset_all);
	std::println(std::cout,"Floating point precision = {} bits",sizeof(tw::Float)*8);
	std::println(std::cout,"Maximum OpenMP threads = {}", omp_get_max_threads());
}

void Simulation::Run()
{
	std::ofstream twstat;
	SetupIO();

	try
	{
		InputFileFirstPass();

		if (strip[0].Get_rank()==0)
		{
			twstat.open("twstat");
			twstat << "TurboWAVE is initializing.";
			twstat.close();
		}

		std::println(std::cout,"\n{}{}Prepare Simulation{}\n",term::green,term::bold,term::reset_all);
		std::flush(std::cout);

		PrepareSimulation();

		std::println(std::cout,"\n{}{}Begin Simulation{}\n",term::green,term::bold,term::reset_all);
		std::flush(std::cout);

		tw::Int startTime = GetSeconds();
		previous_timestamp = startTime;

		if (GetSeconds()<0) {
			std::println(std::cout,"\n{}: System clock is not responding properly.\n",term::warning);
		}

		std::println(std::cout,"Current status can be viewed in 'twstat' file.");
		if (interactive) {
			std::println(std::cout,"Type `help` for interactive commands.");
		} else {
			std::println(std::cout,"interactive thread is disabled");
		}
		std::flush(std::cout);

		while (stepNow <= stepsToTake && windowPosition[0] < maxWindowPosition[0])
		{
			if ((GetSeconds() > previous_timestamp + 5) && strip[0].Get_rank()==0)
			{
				twstat.open("twstat");
				InteractiveCommand("status",&twstat);
				twstat.close();
				previous_timestamp = GetSeconds();
			}
			FundamentalCycle();
		}

		std::println(std::cout,"Completed {} steps in {} seconds.",stepNow,GetSeconds() - startTime);
		std::println(std::cout,"Simulated elapsed time = {}",windowPosition[0]);
		if (strip[0].Get_rank()==0)
		{
			twstat.open("twstat");
			twstat << "Completed " << stepNow << " steps in " << GetSeconds() - startTime << " seconds." << std::endl;
			twstat << "Simulated elapsed time = " << windowPosition[0] << std::endl;
			twstat.close();
		}
		failure_count = 0;
	}
	catch (tw::FatalError& e)
	{
		std::println(std::cerr,"{}: {}",term::error,e.what());
		if (strip[0].Get_rank()==0)
		{
			twstat.open("twstat");
			twstat << "The simulation failed. For more info see stdout." << std::endl;
			twstat.close();
		}
		failure_count = 1;
	}

	if (failure_count==0) {
		std::println(std::cout,"\n{}{}TurboWAVE Session Completed{}",term::green,term::bold,term::reset_all);
	} else {
		std::println(std::cerr,"\n{}{}TurboWAVE Session Failed{}",term::red,term::bold,term::reset_all);
	}
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

	// each test must setup its own task and grid
	if (numRanksProvided!=2)
		throw tw::FatalError("test must use -n 2");
	AttachUnits(tw::units::plasma,1e19);

	tw::Float failed = 0;
	success_count = 0;
	failure_count = 0;
	ComputeTool *tool;
	Module *module;
	
	logger::INFO("Test Compute Tools");

	for (auto [tool_key,theType] : ComputeTool::Map())
	{
		if (unitTest==tool_key || unitTest=="--all")
		{
			std::println(std::cout,"\ntesting {}{}{}{}",term::bold,term::cyan,tool_key,term::reset_all);
			tw::Int enviro = 1;
			while (ComputeTool::SetTestEnvironment(theType,enviro,this,this))
			{
				auto gridStr = std::format("    {}x{}x{} grid",this->GlobalDim(1),this->GlobalDim(2),this->GlobalDim(3));
				tool = CreateTool("test_tool",theType);
				logger::TRACE("register tests using tool prototype...");
				tool->RegisterTests();
				auto test_names = tool->test_names;
				auto test_cases = tool->tests;
				logger::TRACE("remove prototype tool...");
				RemoveTool(tool);
				for (auto id=0; id < test_cases.size(); id++) {
					logger::DEBUG(std::format("begin test case {}",id));
					try {
						tool = CreateTool("test_tool",theType);
						failed = 0;
						try {
							test_cases[id](*tool);
							success_count++;
							std::println(test_out,"{}    {} {}{}{}",gridStr,term::ok,term::green,test_names[id],term::reset_all);
						} catch(tw::FatalError& e) {
							failed = 1;
							std::println(test_out,"{}    {} {}{}{}",gridStr,term::err,term::red,test_names[id],term::reset_all);
							std::println(test_report,"{}. {}: environment {}",++failure_count,tool_key,enviro);
							std::println(test_report,"{}",e.what());
						}
						logger::TRACE("begin MPI sum of failures");
						strip[0].AllSum(&failed,&failed,sizeof(tw::Float),0);
						logger::TRACE("end MPI sum of failures");
						RemoveTool(tool);
					} catch (tw::FatalError& e) {
						std::println(test_out,"{}    {}tool rejected the environment{}",gridStr,term::yellow,term::reset_all);
					}
					logger::DEBUG(std::format("finish test case {}",id));
				}
				logger::TRACE("wait at MPI barrier");
				MPI_Barrier(MPI_COMM_WORLD);
				logger::TRACE("continue from MPI barrier");
				std::print(std::cout,"{}",test_out.str());
				test_out.str("");
				test_out.clear();
				enviro++;
			}
		} else {
			logger::INFO(std::format("skip test of tool {}",tool_key));
		}
	}

	logger::INFO("Test Modules");

	for (auto [mod_key,theType] : Module::Map())
	{
		if (unitTest==mod_key || unitTest=="--all")
		{
			std::println(std::cout,"\ntesting {}{}{}{}",term::bold,term::cyan,mod_key,term::reset_all);
			tw::Int enviro = 1;
			while (Module::SetTestEnvironment(theType,enviro,this))
			{
				auto gridStr = std::format("    {}x{}x{} grid",this->GlobalDim(1),this->GlobalDim(2),this->GlobalDim(3));
				module = factory::CreateModuleFromType("test_module",theType,this);
				logger::TRACE("register tests using module prototype...");
				module->RegisterTests();
				auto test_names = tool->test_names;
				auto test_cases = tool->tests;
				logger::TRACE("remove prototype module...");
				delete module;
				for (auto id = 0; id < test_cases.size(); id++) {
					logger::DEBUG(std::format("begin test case {}",id));
					try {
	 					module = factory::CreateModuleFromType("test_module",theType,this);
						failed = 0;
						try {
							test_cases[id](*module);
							success_count++;
							std::println(test_out,"{}    {} {}{}{}",gridStr,term::ok,term::green,test_names[id],term::reset_all);
						} catch(tw::FatalError& e) {
							failed = 1;
							std::println(test_out,"{}    {} {}{}{}",gridStr,term::err,term::red,test_names[id],term::reset_all);
							std::println(test_report,"{}. {}: environment {}",++failure_count,mod_key,enviro);
							std::println(test_report,"{}",e.what());
						}
						logger::TRACE("begin MPI sum of failures");
						strip[0].AllSum(&failed,&failed,sizeof(tw::Float),0);
						logger::TRACE("end MPI sum of failures");
						delete module;
					} catch (tw::FatalError& e) {
						std::println(test_out,"    {}module rejected the environment{}",term::yellow,term::reset_all);
					}
					logger::DEBUG(std::format("finish test case {}",id));
				}
				logger::TRACE("wait at MPI barrier");
				MPI_Barrier(MPI_COMM_WORLD);
				logger::TRACE("continue from MPI barrier");
				std::print(std::cout,"{}",test_out.str());
				test_out.str("");
				test_out.clear();
				enviro++;
			}
		} else {
			logger::INFO(std::format("skip test of module {}",mod_key));
		}
	}

	std::println(std::cout,"");
	if (test_report.str().size()>0) {
		std::println(std::cerr,"{}{}Unit Tests Failing{}\n",term::bold,term::red,term::reset_all);
		std::println(std::cout,"{}",test_report.str());
	} else if (failure_count > 0) {
		std::println(std::cerr,"{}{}Unit Tests Failing{}\n",term::bold,term::red,term::reset_all);
		std::println(std::cerr,"Details may be in node-scoped output files");
	} else {
		std::println(std::cout,"{}{}Unit Tests Passing{} - {} succeeded, {} failed",term::bold,term::green,term::reset_all,success_count,failure_count);
	}
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

	std::println(std::cout,"\nInitializing Compute Tools...\n");

	for (auto tool : computeTool)
	{
		std::println(std::cout,"Tool: {}",tool->name);
		tool->Initialize();
		tool->WarningMessage();
	}
	if (computeTool.size()==0)
		std::println(std::cout,"(no tools)");

	// Initialize Modules

	std::println(std::cout,"\nInitialize Modules...\n");

	for (auto m : module)
		m->ExchangeResources();

	for (auto m : module)
	{
		std::println(std::cout,"Module: {}",m->name);
		m->Initialize();
		m->WarningMessage();
	}

	// Read checkpoint data

	if (restartFileName!="tw::none")
	{
		std::stringstream fileName;
		std::ifstream restartFile;
		fileName << strip[0].Get_rank() << "_dump.chk";
		std::println(std::cout,"\nReading restart file {}...",fileName.str());
		restartFile.open(fileName.str().c_str());
		ReadCheckpoint(restartFile);
		restartFile.close();
	}
}

#ifdef USE_OPENCL

void Simulation::PrintGPUInformation()
{
	cl_ulong ninfo;

	std::println(initMessage);

	std::println(std::cout,"GPU INFORMATION");
	std::println(std::cout,"--------------------------------------------");

	clGetDeviceInfo(gpu,CL_DEVICE_GLOBAL_MEM_SIZE,sizeof(ninfo),&ninfo,NULL);
	std::println(std::cout,"Global memory: {}",ninfo);

	clGetDeviceInfo(gpu,CL_DEVICE_LOCAL_MEM_SIZE,sizeof(ninfo),&ninfo,NULL);
	std::println(std::cout,"Local memory: {}",ninfo);

	clGetDeviceInfo(gpu,CL_DEVICE_MAX_WORK_GROUP_SIZE,sizeof(ninfo),&ninfo,NULL);
	std::println(std::cout,"Maximum work group size: {}",ninfo);

	clGetDeviceInfo(gpu,CL_DEVICE_PREFERRED_VECTOR_WIDTH_FLOAT,sizeof(ninfo),&ninfo,NULL);
	std::println(std::cout,"Float vector width: {}",ninfo);

	clGetDeviceInfo(gpu,CL_DEVICE_PREFERRED_VECTOR_WIDTH_DOUBLE,sizeof(ninfo),&ninfo,NULL);
	std::println(std::cout,"Double vector width: {}",ninfo);

	std::println(std::cout,"--------------------------------------------\n");

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
		*theStream << "Current elapsed time: " << windowPosition[0] << std::endl;
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
		*theStream << "Global grid size: " << GlobalDim(1) << "," << GlobalDim(2) << "," << GlobalDim(3) << std::endl;
		*theStream << "Local grid size: " << Dim(1) << "," << Dim(2) << "," << Dim(3) << std::endl;
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
	logger::TRACE(std::format("step {}",stepNow));
	Diagnose();

	for (auto m : module)
		m->Reset();

	for (auto m : module)
		m->Update();

	solutionPosition += spacing[0]*solutionVelocity;
	altSolutionPosition += spacing[0]*(solutionVelocity - tw::vec4(0,0,0,1));
	stepNow++;

	if (adaptiveGrid)
		for (auto m : module)
			m->AdaptGrid();

	corner[0] += spacing[0];
	windowPosition[0] += spacing[0];
	
	if (movingWindow && solutionPosition[3]>=windowPosition[3] + spacing[3] && dim[3]>1) {
		logger::TRACE(std::format("move lab at {:.5} triggered by {:.5}",windowPosition[3],solutionPosition[3]));
		MoveWindow();
	} else {
		logger::TRACE(std::format("hold lab at {:.5} restrained by {:.5}",windowPosition[3],solutionPosition[3]));
	}

	if (!movingWindow && altSolutionPosition[3]<=altWindowPosition[3] - spacing[3] && dim[3]>1) {
		logger::TRACE(std::format("move alt at {:.5} triggered by {:.5}",altWindowPosition[3],altSolutionPosition[3]));
		AntiMoveWindow();
	} else {
		logger::TRACE(std::format("hold alt at {:.5} restrained by {:.5}",altWindowPosition[3],altSolutionPosition[3]));
	}

	// RESTART MECHANISM

	tw::Int curr = strip[0].Get_rank();
	bool doing_restart = dumpPeriod>0 && stepNow%dumpPeriod==0;
	if (doing_restart)
	{
		std::ofstream restartFile;
		std::string fileName = std::to_string(curr) + "_dump.chk";
		restartFile.open(fileName.c_str(),std::ios::binary);
		WriteCheckpoint(restartFile);
		restartFile.close();
	}
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
	logger::INFO(std::format("Creating Tool <{}>...",name));
	computeTool.push_back(factory::CreateToolFromType(name,theType,this,this));
	computeTool.back()->refCount++;
	return computeTool.back();
}

ComputeTool* Simulation::GetTool(const std::string& name,bool attaching)
{
	for (auto i=0;i<computeTool.size();i++)
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

bool Simulation::RemoveTool(ComputeTool *theTool)
{
	auto iter = std::find(computeTool.begin(),computeTool.end(),theTool);
	if (iter==computeTool.end())
		throw tw::FatalError("Attempt to remove a non-existant tool.");
	(*iter)->refCount--;
	logger::TRACE(std::format("remove tool to refcount {}",(*iter)->refCount));
	if ((*iter)->refCount==0)
	{
		logger::WARN("tool was not removed");
		//delete *iter;
		//computeTool.erase(iter);
		//return true;
	}
	return false;
}

void Simulation::MoveWindow()
{
	windowPosition[3] += spacing[3];
	corner[3] += spacing[3];
	globalCorner[3] += spacing[3];
	
	for (auto i=lfg[3];i<=ufg[3];i++)
		X(i,3) += spacing[3];

	for (auto i=0;i<clippingRegion.size();i++)
		if (clippingRegion[i]->moveWithWindow)
			clippingRegion[i]->Translate(tw::vec3(0,0,spacing[3]));
		else
			clippingRegion[i]->Initialize(*this,this);

	for (auto i=0;i<module.size();i++)
		module[i]->MoveWindow();
}

void Simulation::AntiMoveWindow()
{
	altWindowPosition[3] -= spacing[3];

	for (auto i=0;i<module.size();i++)
		module[i]->AntiMoveWindow();
}

void Simulation::ReadCheckpoint(std::ifstream& inFile)
{
	std::string objectName;

	MetricSpace::ReadCheckpoint(inFile);
	inFile.read((char *)&stepNow,sizeof(tw::Int));
	stepNow++;

	if (uniformDeviate!=NULL)
		uniformDeviate->ReadCheckpoint(inFile);

	if (gaussianDeviate!=NULL)
		gaussianDeviate->ReadCheckpoint(inFile);

	// Read Region objects

	for (tw::Int i=0;i<clippingRegion.size();i++)
	{
		inFile >> objectName;
		inFile.ignore();
		std::println(std::cout,"Read checkpoint data for region <{}>...",objectName);
		Region::FindRegion(clippingRegion,objectName)->ReadCheckpoint(inFile);
	}

	// Read ComputeTool objects

	for (tw::Int i=0;i<computeTool.size();i++)
	{
		inFile >> objectName;
		inFile.ignore();
		std::println(std::cout,"Read checkpoint data for tool <{}>...",objectName);
		GetTool(objectName,false)->ReadCheckpoint(inFile);
	}

	// Read Module objects

	for (tw::Int i=0;i<module.size();i++)
	{
		inFile >> objectName;
		inFile.ignore();
		std::println(std::cout,"Read checkpoint data for module <{}>...",objectName);
		GetModule(objectName)->ReadCheckpoint(inFile);
	}
}

void Simulation::WriteCheckpoint(std::ofstream& outFile)
{
	MetricSpace::WriteCheckpoint(outFile);
	outFile.write((char *)&stepNow,sizeof(tw::Int));

	if (uniformDeviate!=NULL)
		uniformDeviate->WriteCheckpoint(outFile);

	if (gaussianDeviate!=NULL)
		gaussianDeviate->WriteCheckpoint(outFile);

	for (auto obj : clippingRegion)
	{
		std::println(std::cout,"Checkpointing <{}>...",obj->name);
		obj->WriteCheckpoint(outFile);
	}

	for (auto obj : computeTool)
	{
		std::println(std::cout,"Checkpointing <{}>...",obj->name);
		obj->WriteCheckpoint(outFile);
	}

	for (auto obj : module)
	{
		std::println(std::cout,"Checkpointing <{}>...",obj->name);
		obj->WriteCheckpoint(outFile);
	}
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
		doing_diagnostics = doing_diagnostics || d->WriteThisStep(windowPosition[0],spacing[0],stepNow);
	if (doing_diagnostics)
	{
		logger::TRACE(std::format("prepare diagnostics"));
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
		logger::TRACE(std::format("diagnostic {}",d->name));
		if (d->WriteThisStep(windowPosition[0],spacing[0],stepNow))
		{
			d->Start();
			d->ReportNumber("time",windowPosition[0],true);
			d->ReportNumber("dt",spacing[0],true);
			d->ReportNumber("zwindow",windowPosition[3],true);
			for (auto m : module)
				if (has_diagnostic(m,d))
					m->Report(*d);
			d->Finish();
		}
	}
}
