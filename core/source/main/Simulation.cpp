module;

#include <tree_sitter/api.h>
#include "tw_includes.h"
#include "tw_logger.h"
#ifndef USE_STD_MODULE
    #include <algorithm>
#endif

export module simulation;
import base;
import metric_space;
import driver;
import logger;
import input;
import read_grid;
import factory;

/// This class is the root of the driver tree.
/// This basically orchestrates everything.
export struct Simulation: Driver, tw::input::Visitor
{
	bool interactive;
	tw::Int outputLevel,errorCheckingLevel,success_count,failure_count;
	tw::Int previous_timestamp;

	tw::Int dumpPeriod;
	tw::Float dummy_unit_density; // prevent unknown key error
	std::string dummy_native_units; // prevent unknown key error

	tw::Int inputFilePass;
	tw::input::DirectiveReader outerDirectives;
	std::unique_ptr<GridReader> gridReader;
	std::string raw_src; // raw input file
	std::string src; // expanded input file

	Simulation(const bool interactive,
		const std::string& unitTest,
		const std::string& inputFileName,
		const std::string& restartFileName,
		const std::string& platform,
		const std::string& device,
		const tw::Int& outputLevel,
		const tw::Int& errorCheckingLevel,
		MetricSpace *ms,
		Task *tsk);
	void SetupIO();
	virtual void Run();
	virtual void Test();
	void PrepareSimulation();
	void FundamentalCycle();
	void MoveWindow();
	void AntiMoveWindow();
	void Diagnose();

	tw::input::navigation visit(TSTreeCursor *curs);
	void InputFileFirstPass();
	void ParseNestedDeclaration(TSTreeCursor *curs,const std::string& src,Driver *sup);
	Driver* AutoCreateSupers(tw::tool_type reqType,const std::string& basename);
	void ReadInputFile();
	void ReadCheckpoint(std::ifstream& inFile);
	void WriteCheckpoint(std::ofstream& outFile);
	void InteractiveCommand(const std::string& cmd,std::ostream *theStream);

	#ifdef USE_OPENCL
	void PrintGPUInformation();
	void CellUpdateProtocol(cl_kernel k)
	{
		MetricSpace::CellUpdateProtocol(k,commandQueue);
	}
	void ElementUpdateProtocol(cl_kernel k)
	{
		MetricSpace::ElementUpdateProtocol(k,commandQueue);
	}
	void LocalUpdateProtocol(cl_kernel k)
	{
		MetricSpace::LocalUpdateProtocol(k,commandQueue);
	}
	void PointUpdateProtocol(cl_kernel k)
	{
		MetricSpace::PointUpdateProtocol(k,commandQueue);
	}
	void StripUpdateProtocol(cl_kernel k,tw::Int axis,tw::Int stripArgument)
	{
		MetricSpace::StripUpdateProtocol(k,commandQueue,axis,stripArgument);
	}
	#endif


};

Simulation::Simulation(const bool interactive,
	const std::string& test_name,
	const std::string& file_name,
	const std::string& restart_name,
	const std::string& platform,
	const std::string& device,
	const tw::Int& outputLevel,
	const tw::Int& errorCheckingLevel,
	MetricSpace *ms,
	Task *tsk) : Driver("simulation",ms,tsk)
{
	this->interactive = interactive;
	task->unitTest = test_name;
	task->inputFileName = file_name;
	task->restartFileName = restart_name;
	task->platformSearchString = platform;
	task->deviceSearchString = device;
	this->outputLevel = outputLevel;
	this->errorCheckingLevel = errorCheckingLevel;

	previous_timestamp = 0;
	dumpPeriod = 0;
	inputFilePass = 0;

	outerDirectives.Add("dump period",new tw::input::Int(&dumpPeriod),false);
	// this gets read directly before tree-walking, just throw it out during tree-walking
	outerDirectives.Add("unit density",new tw::input::Float(&dummy_unit_density),false);
	// this gets read directly before tree-walking, just throw it out during tree-walking
	outerDirectives.Add("native units",new tw::input::String(&dummy_native_units),false);
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

		if (task->strip[0].Get_rank()==0)
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

		while (space->IsNotFinished())
		{
			if ((GetSeconds() > previous_timestamp + 5) && task->strip[0].Get_rank()==0)
			{
				twstat.open("twstat");
				InteractiveCommand("status",&twstat);
				twstat.close();
				previous_timestamp = GetSeconds();
			}
			FundamentalCycle();
		}

		std::println(std::cout,"Completed {} steps in {} seconds.",space->StepNow(),GetSeconds() - startTime);
		std::println(std::cout,"Simulated elapsed time = {}",space->WindowPos(0));
		if (task->strip[0].Get_rank()==0)
		{
			twstat.open("twstat");
			twstat << "Completed " << space->StepNow() << " steps in " << GetSeconds() - startTime << " seconds." << std::endl;
			twstat << "Simulated elapsed time = " << space->WindowPos(0) << std::endl;
			twstat.close();
		}
		failure_count = 0;
	}
	catch (tw::FatalError& e)
	{
		std::println(std::cerr,"{}: {}",term::error,e.what());
		if (task->strip[0].Get_rank()==0)
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
	task->uniformDeviate = new UniformDeviate(1 + task->strip[0].Get_rank()*(MaxSeed()/numRanksProvided));
	task->gaussianDeviate = new GaussianDeviate(1 + task->strip[0].Get_rank()*(MaxSeed()/numRanksProvided) + MaxSeed()/(2*numRanksProvided));

	// diverting output
	std::stringstream test_out;
	std::stringstream test_report;

	// each test must setup its own task and grid
	if (numRanksProvided!=2)
		throw tw::FatalError("test must use -n 2");
	const tw::Float testUnitDensityCGS = 1e19;
	space->AttachUnits(tw::units::plasma,testUnitDensityCGS);

	tw::Float failed = 0;
	success_count = 0;
	failure_count = 0;
	
	std::println(std::cout,"\nTest Compute Tools");

	for (auto [tool_key,theType] : ComputeTool::Map())
	{
		if (tw::IsDriver(theType)) {
			continue;
		}
		if (task->unitTest==tool_key || task->unitTest=="--all")
		{
			std::println(std::cout,"\ntesting {}{}{}{}",term::bold,term::cyan,tool_key,term::reset_all);
			tw::Int enviro = 1;
			while (ComputeTool::SetTestEnvironment(theType,enviro,space,task,testUnitDensityCGS))
			{
				auto gridStr = std::format("    {}x{}x{} grid",space->GlobalDim(1),space->GlobalDim(2),space->GlobalDim(3));
				for (auto id=0; ; id++) {
					try {
						auto tool = CreateTool("test_tool",theType);
						tool->RegisterTests();
						if (id >= tool->tests.size() || id > 100) {
							break;
						}
						logger::DEBUG(std::format("begin test case {}",id));
						failed = 0;
						try {
							tool->tests[id](*tool);
							success_count++;
							std::println(test_out,"{}    {} {}{}{}",gridStr,term::ok,term::green,tool->test_names[id],term::reset_all);
						} catch(tw::FatalError& e) {
							failed = 1;
							std::println(test_out,"{}    {} {}{}{}",gridStr,term::err,term::red,tool->test_names[id],term::reset_all);
							std::println(test_report,"{}. {}: environment {}",++failure_count,tool_key,enviro);
							std::println(test_report,"{}",e.what());
						}
						logger::TRACE("begin MPI sum of failures");
						task->strip[0].AllSum(&failed,&failed,sizeof(tw::Float),0);
						logger::TRACE("end MPI sum of failures");
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

	std::println(std::cout,"\nTest Modules");

	for (auto [mod_key,theType] : ComputeTool::Map())
	{
		if (!tw::IsDriver(theType)) {
			continue;
		}
		if (task->unitTest==mod_key || task->unitTest=="--all")
		{
			std::println(std::cout,"\ntesting {}{}{}{}",term::bold,term::cyan,mod_key,term::reset_all);
			tw::Int enviro = 1;
			while (Driver::SetTestEnvironment(theType,enviro,space,task,testUnitDensityCGS))
			{
				auto gridStr = std::format("    {}x{}x{} grid",space->GlobalDim(1),space->GlobalDim(2),space->GlobalDim(3));
				for (auto id = 0; ; id++) {
					try {
	 					auto driver = factory::CreateDriverFromType("test_module",theType,space,task);
						driver->RegisterTests();
						if (id >= driver->tests.size() || id > 100) {
							delete driver;
							break;
						}
						logger::DEBUG(std::format("begin test case {}",id));
						failed = 0;
						try {
							driver->tests[id](*driver);
							logger::TRACE("test succeeded");
							success_count++;
							std::println(test_out,"{}    {} {}{}{}",gridStr,term::ok,term::green,driver->test_names[id],term::reset_all);
						} catch(tw::FatalError& e) {
							failed = 1;
							std::println(test_out,"{}    {} {}{}{}",gridStr,term::err,term::red,driver->test_names[id],term::reset_all);
							std::println(test_report,"{}. {}: environment {}",++failure_count,mod_key,enviro);
							std::println(test_report,"{}",e.what());
						}
						logger::TRACE("begin MPI sum of failures");
						task->strip[0].AllSum(&failed,&failed,sizeof(tw::Float),0);
						logger::TRACE("end MPI sum of failures");
						delete driver;
					} catch (tw::FatalError& e) {
						std::println(test_out,"    {}driver rejected the environment{}",term::yellow,term::reset_all);
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

	// Attach tools to engines (root level only)
	for (auto tool : tools) {
		auto engine = std::dynamic_pointer_cast<Engine>(tool);
		if (engine) {
			engine->AttachTools(tools);
		}
	}

	// Call VerifyInput on root-level drivers, sub-drivers have to be explicitly called within.
	// A common task is to verify that needed tools are attached.
	for (auto d : sub_drivers) {
		logger::DEBUG(std::format("verify input for {}",d->name));
		d->VerifyInput();
	}

	// Start the metadata dictionary
	if (task->strip[0].Get_rank()==0)
	{
		std::ofstream metadata_file("tw_metadata.json");
		metadata_file << "{" << std::endl;
		metadata_file << "\"version\": \"" << TW_VERSION_STRING << "\"}";
		// n.b. closing brace must not be followed by anything due to append strategy
		metadata_file.close();
	}

	// If a diagnostic tool is not attached to any module attach it to all modules
	for (auto tool : tools) {
		// n.b. loop variable increments use_count
		logger::DEBUG(std::format("tool {} used by {}",tool->name,tool.use_count()-1));
		if (std::dynamic_pointer_cast<Diagnostic>(tool)) {
			if (tool.use_count()==2) {
				for (auto d : sub_drivers) {
					d->tools.push_back(tool);
				}
			}
		}
	}

	// Check for duplicate filenames
	std::vector<std::string> fileName;
	for (auto tool : tools) {
		auto d = std::dynamic_pointer_cast<Diagnostic>(tool);
		if (d) {
			if (std::find(fileName.begin(),fileName.end(),d->filename)==fileName.end()) {
				fileName.push_back(d->filename);
			} else {
				throw tw::FatalError("Duplicate file name used in diagnostic <"+d->name+">.");
			}
		}
	}

	// Sort Drivers

	DriverComparator comparatorObject;
	std::sort(sub_drivers.begin(),sub_drivers.end(),comparatorObject);

	// Initialize tools that are on the root list.
	// Must precede module initialization.
	// If a tool is held privately it has to be initialized explicitly by the driver.
	// TODO: two passes, non-engines followed by engines, e.g. so regions are intialized first

	std::println(std::cout,"\nInitializing Shared Tools...\n");
	std::flush(std::cout);

	for (auto tool : tools)
	{
		std::println(std::cout,"Tool: {}",tool->name);
		std::flush(std::cout);
		tool->Initialize();
	}
	if (tools.size()==0)
		std::println(std::cout,"(no shared tools)");

	// Initialize root level drivers.
	// Super-drivers have to explicitly initialize sub-drivers.

	std::println(std::cout,"\nInitialize Drivers...\n");
	std::flush(std::cout);

	for (auto d : sub_drivers)
		d->ExchangeResources();

	for (auto d : sub_drivers)
	{
		std::println(std::cout,"Driver: {}",d->name);
		std::flush(std::cout);
		d->Initialize();
	}

	RecursiveTreeDisplay(&std::cout,0,16);

	// Read checkpoint data

	if (task->restartFileName!="tw::none")
	{
		std::stringstream fileName;
		std::ifstream restartFile;
		fileName << task->strip[0].Get_rank() << "_dump.chk";
		std::println(std::cout,"\nReading restart file {}...",fileName.str());
		std::flush(std::cout);
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
	if (cmd=="help" || cmd=="?") {
		*theStream << "--- List of Interactive Commands ---" << std::endl;
		*theStream << "status or /: print current step and other status indicators" << std::endl;
		*theStream << "metrics : print grid and time step metrics for this simulation" << std::endl;
		*theStream << "list : list modules and compute tools and their ID numbers" << std::endl;
		//*theStream << "peek [x] [y] [z] : print current data at cell x,y,z" << std::endl;
		*theStream << "Ctrl-C : abort the simulation" << std::endl;
		*theStream << std::endl;
	} else if (cmd=="status" || cmd=="/") {
		*theStream << "Current step: " << space->StepNow() << std::endl;
		*theStream << "Current step size: " << space->dx(0) << std::endl;
		*theStream << "Current elapsed time: " << space->WindowPos(0) << std::endl;
		for (auto d : sub_drivers)
			d->StatusMessage(theStream);
		for (auto tool : tools)
			tool->StatusMessage(theStream);
		*theStream << std::endl;
	} else if (cmd=="list") {
		*theStream << "--- Driver Tree ---" << std::endl;
		RecursiveTreeDisplay(theStream,0,16);
	} else if (cmd=="metrics") {
		*theStream << "Steps to take: " << space->StepsToTake() << std::endl;
		*theStream << "Steps remaining: " << space->StepsToTake() - space->StepNow() << std::endl;
		*theStream << "Global grid size: " << space->GlobalDim(1) << "," << space->GlobalDim(2) << "," << space->GlobalDim(3) << std::endl;
		*theStream << "Local grid size: " << space->Dim(1) << "," << space->Dim(2) << "," << space->Dim(3) << std::endl;
		*theStream << "MPI Domains: " << task->domains[1] << "," << task->domains[2] << "," << task->domains[3] << std::endl;
		*theStream << std::endl;
	} else if (cmd.find("peek")!=std::string::npos) {
		*theStream << "Not implemented yet.";
		*theStream << std::endl;
	} else {
		*theStream << "Unknown command.";
		*theStream << std::endl;
	}
}

void Simulation::FundamentalCycle()
{
	logger::DEBUG(std::format("step {}",space->StepNow()));
	Diagnose();

	for (auto d : sub_drivers) {
		logger::DEBUG(std::format("reset {}",d->name));
		d->Reset();
	}

	for (auto d : sub_drivers) {
		logger::DEBUG(std::format("update {}",d->name));
		d->Update();
	}

	if (space->IsAxisAdaptive(3))
		for (auto d : sub_drivers)
			d->AdaptGrid();

	auto shift = space->Advance(1);
	if (shift==-1) {
		AntiMoveWindow();
	} else if (shift==1) {
		MoveWindow();
	}

	// RESTART MECHANISM

	tw::Int curr = task->strip[0].Get_rank();
	bool doing_restart = dumpPeriod>0 && space->StepNow()%dumpPeriod==0;
	if (doing_restart)
	{
		logger::INFO(std::format("checkpoint {}",space->StepNow()/dumpPeriod));
		std::ofstream restartFile;
		std::string fileName = std::to_string(curr) + "_dump.chk";
		restartFile.open(fileName.c_str(),std::ios::binary);
		WriteCheckpoint(restartFile);
		restartFile.close();
	}
}

void Simulation::MoveWindow()
{
	// DynSpace did not advance this
	for (auto i=space->LFG(3); i<=space->UFG(3); i++) {
		space->X(i,3) += space->dx(3);
	}
	// TODO: clipping regions formerly advanced their coordinates here, now we will instead
	// apply a transformation to unwind moving window the same as rotations etc.
	// Also, there was a call to Initialize here, why?
	for (auto d : sub_drivers)
		d->MoveWindow();
}

void Simulation::AntiMoveWindow()
{
	for (auto d : sub_drivers)
		d->AntiMoveWindow();
}

void Simulation::ReadCheckpoint(std::ifstream& inFile)
{
	std::string objectName;

	space->ReadCheckpoint(inFile);

	if (task->uniformDeviate!=NULL)
		task->uniformDeviate->ReadCheckpoint(inFile);

	if (task->gaussianDeviate!=NULL)
		task->gaussianDeviate->ReadCheckpoint(inFile);

	// TODO: need to reed a more tree-like structure, whereas before we could rely on flat lists.
	// We probably need to forbid checkpointing anything in a shared tool.
}

void Simulation::WriteCheckpoint(std::ofstream& outFile)
{
	space->WriteCheckpoint(outFile);

	if (task->uniformDeviate!=NULL)
		task->uniformDeviate->WriteCheckpoint(outFile);

	if (task->gaussianDeviate!=NULL)
		task->gaussianDeviate->WriteCheckpoint(outFile);
}

void Simulation::Diagnose()
{
	std::vector<std::shared_ptr<Diagnostic>> diagnostics;
	for (auto tool : tools)
		if (std::dynamic_pointer_cast<Diagnostic>(tool))
			diagnostics.push_back(std::dynamic_pointer_cast<Diagnostic>(tool));

	// DIAGNOSTIC PREP
	// Main purpose of this step is to let modules transfer data from compute devices.
	// This has a high cost, so only do if necessary.

	bool doing_diagnostics=false;
	for (auto diag : diagnostics)
		doing_diagnostics = doing_diagnostics || diag->WriteThisStep();
	if (!doing_diagnostics) {
		return;
	}

	logger::TRACE(std::format("prepare diagnostics"));
	for (auto d : sub_drivers)
		d->StartDiagnostics();

	// MAIN DIAGNOSTIC LOOP

	auto has_diagnostic = [&] (Driver *d,std::shared_ptr<Diagnostic> diagnostic)
	{
		for (auto maybe_diag : d->tools)
			if (maybe_diag==diagnostic)
				return true;
		return false;
	};

	for (auto diag : diagnostics)
	{
		logger::TRACE(std::format("diagnostic {}",diag->name));
		if (diag->WriteThisStep())
		{
			diag->Start();
			diag->ReportNumber("time",space->WindowPos(0),true);
			diag->ReportNumber("dt",space->dx(0),true);
			diag->ReportNumber("zwindow",space->WindowPos(3),true);
			for (auto sd : sub_drivers) {
				if (has_diagnostic(sd,diag)) {
					logger::TRACE(std::format("{} is reporting",sd->name));
					sd->Report(*diag);
				}
			}
			diag->Finish();
		}
	}
}
