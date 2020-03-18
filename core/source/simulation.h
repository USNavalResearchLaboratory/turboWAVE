#include "meta_base.h"
#include "meta_tools.h"

struct Simulation;

#include "module.h"

struct Simulation:Task,MetricSpace
{
	std::string inputFileName,restartFileName;
	std::ostream *tw_out;
	tw::input::DirectiveReader outerDirectives,gridDirectives;

	tw::grid::geometry gridGeometry;
	tw::Float dt0,dtMin,dtMax,dtCritical,elapsedTime,elapsedTimeMax;
	tw::Float signalPosition,windowPosition,signalSpeed;
	tw::Float antiSignalPosition,antiWindowPosition;
	tw::Float unitDensityCGS;

	bool neutralize,movingWindow,adaptiveTimestep,adaptiveGrid;
	bool completed;
	tw::Int outputLevel,errorCheckingLevel;

	tw::Int stepsToTake,stepNow;
	tw::Int lastTime;

	tw::Int dumpPeriod;

	tw::bc::par bc0[4],bc1[4];

	std::vector<Region*> clippingRegion;
	std::vector<ComputeTool*> computeTool;
	std::vector<Module*> module;
	std::vector<tw::module_type> createdModuleTypes;

	Simulation(const std::string& inputFileName,const std::string& restartFileName);
	virtual ~Simulation();
	virtual void Run();
	void SetupGeometry();
	void PrepareSimulation();
	void FundamentalCycle();
	void MoveWindow();
	void AntiMoveWindow();
	void Diagnose();

	bool MangleModuleName(std::string& name);
	bool CheckModule(const std::string& name);
	Module* GetModule(const std::string& name);
	tw::Int FindModule(const std::string& name);

	bool MangleToolName(std::string& name);
	ComputeTool* CreateTool(const std::string& basename,tw::tool_type theType);
	ComputeTool* GetTool(const std::string& name,bool attaching);
	void ToolFromDirective(std::vector<ComputeTool*>& tool,std::stringstream& inputString,const std::string& command);
	bool RemoveTool(ComputeTool *theTool);

	void SetCellWidthsAndLocalSize();
	void SetGlobalSizeAndLocalCorner();

	tw::Float ToLab(tw::Float zeta,tw::Float relativeTime);
	tw::Float ToLight(tw::Float z,tw::Float relativeTime);
	template <class T,class U>
	U ValueOnLabGrid(T& A,tw::strip s,tw::Int k,tw::Float relativeTime);
	template <class T,class U>
	U ValueOnLightGrid(T& A,tw::strip s,tw::Int k,tw::Float relativeTime);

	void OpenInputFile(std::ifstream& inFile);
	std::string InputFileFirstPass();
	void SetupLocalGrid();
	void NestedDeclaration(const std::string& com,std::stringstream& inputString,Module *sup);
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

	void UpdateTimestep(tw::Float dt0)
	{
		SetupTimeInfo(dt0);
		for (auto m : module)
			m->SetupTimeInfo(dt0);
	}

	bool IsFirstStep()
	{
		return stepNow == 0;
	}

	bool IsLastStep()
	{
		// Actual steps taken = stepsToTake+1
		// This allows us to write both the initial condition and the last available data.
		// If there is a restart the initial step is stepsToTake+1 and 1 less step is taken.
		return stepNow == stepsToTake;
	}

	bool Completed()
	{
		return completed;
	}
};

inline tw::Float Simulation::ToLab(tw::Float zeta,tw::Float relativeTime)
{
	return zeta + antiWindowPosition - windowPosition + signalSpeed*(elapsedTime + relativeTime - dth);
}

inline tw::Float Simulation::ToLight(tw::Float z,tw::Float relativeTime)
{
	return z + windowPosition - antiWindowPosition - signalSpeed*(elapsedTime + relativeTime - dth);
}

template <class T,class U>
U Simulation::ValueOnLabGrid(T& A,tw::strip s,tw::Int k,tw::Float relativeTime)
{
	// Take a quantity known on the light grid and get its value in a cell of the lab grid
	tw::Int klight;
	tw::Float z,zeta,w;
	z = Pos(s,k).z - corner.z;
	zeta = ToLight(z,relativeTime);
	klight = MyFloor(zeta*freq.z + 0.5001);
	w = 0.5 - zeta*freq.z + tw::Float(klight);
	return w*A(s,klight) + (one - w)*A(s,klight+1);
}

template <class T,class U>
U Simulation::ValueOnLightGrid(T& A,tw::strip s,tw::Int k,tw::Float relativeTime)
{
	// Take a quantity known on the lab grid and get its value in a cell of the light grid
	tw::Int klab;
	tw::Float z,zeta,w;
	zeta = Pos(s,k).z - corner.z;
	z = ToLab(zeta,relativeTime);
	klab = MyFloor(z*freq.z + 0.4999);
	w = 0.5 - z*freq.z + tw::Float(klab);
	return w*A(s,klab) + (one - w)*A(s,klab+1);
}
