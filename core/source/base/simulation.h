#include "meta_base.h"
#include "meta_tools.h"

struct Simulation;

#include "module.h"

/// This class reads the grid block from the input file.  The grid block
/// contains data pertinent to both `Task` and `MetricSpace`, hence the
/// need for a broker class.
class GridReader
{
	tw::input::DirectiveReader directives;
	tw::grid::geometry geo;
	tw::Int req_dim[4],req_dom[4];
	tw::vec4 relativeRef,absoluteRef,spacing;
	bool adaptiveTimestep,adaptiveGrid,found;
	public:
	GridReader(tw::UnitConverter& uc);
	bool Read(const TSTreeCursor *curs,const std::string& src);
	void UpdateTask(Task& tsk);
	void UpdateSpace(MetricSpace& ms);
	tw::vec4 GlobalCorner();
	tw::vec4 GlobalSize();
	tw::grid::geometry Geometry() { return geo; }
	bool FoundGrid() { return found; }
};

struct Simulation:Task,MetricSpace,tw::input::Visitor
{
	std::ostream *tw_out,*tw_err;
	tw::input::DirectiveReader outerDirectives;
	GridReader *gridReader;

	tw::Float dtMin,dtMax,dtCritical,elapsedTime,elapsedTimeMax;
	tw::Float signalPosition,windowPosition,signalSpeed;
	tw::Float antiSignalPosition,antiWindowPosition;
	tw::Float unitDensityCGS;
	tw::units nativeUnits;

	bool neutralize,movingWindow;
	bool completed;
	tw::Int outputLevel,errorCheckingLevel,inputFilePass;
	std::string src; // current input document

	tw::Int stepNow;
	tw::Int lastTime;

	tw::Int dumpPeriod;

	tw::bc::par bc0[4],bc1[4];

	std::vector<Region*> clippingRegion;
	std::vector<ComputeTool*> computeTool;
	std::vector<Module*> module;
	// Map of the most recently created module of a given type
	std::map<tw::module_type,Module*> module_map;

	Simulation(const std::string& unitTest,
		const std::string& inputFileName,
		const std::string& restartFileName,
		const std::string& platform,
		const std::string& device,
		const tw::Int& outputLevel,
		const tw::Int& errorCheckingLevel);
	virtual ~Simulation();
	void SetupIO();
	virtual void Run();
	virtual void Test();
	void PrepareSimulation();
	void FundamentalCycle();
	void MoveWindow();
	void AntiMoveWindow();
	void Diagnose();

	bool MangleModuleName(std::string& name);
	bool CheckModule(const std::string& name);
	Module* GetModule(const std::string& name);

	bool MangleToolName(std::string& name);
	ComputeTool* CreateTool(const std::string& basename,tw::tool_type theType);
	ComputeTool* GetTool(const std::string& name,bool attaching);
	void ToolFromDirective(std::vector<ComputeTool*>& tool,TSTreeCursor *curs,const std::string& src);
	bool RemoveTool(ComputeTool *theTool);

	tw::Float ToLab(tw::Float zeta,tw::Float relativeTime);
	tw::Float ToLight(tw::Float z,tw::Float relativeTime);
	template <class T,class U>
	U ValueOnLabGrid(T& A,tw::strip s,tw::Int k,tw::Float relativeTime);
	template <class T,class U>
	U ValueOnLightGrid(T& A,tw::strip s,tw::Int k,tw::Float relativeTime);

	tw::input::navigation visit(TSTreeCursor *curs);
	tw::input::navigation descend(TSTreeCursor *curs);
	void InputFileFirstPass();
	void NestedDeclaration(TSTreeCursor *curs,const std::string& src,Module *sup);
	Module* RecursiveAutoSuper(tw::module_type reqType,const std::string& basename);
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
		return stepNow == dim[0];
	}

	bool Completed()
	{
		return completed;
	}
};

inline tw::Float Simulation::ToLab(tw::Float zeta,tw::Float relativeTime)
{
	return zeta + antiWindowPosition - windowPosition + signalSpeed*(elapsedTime + relativeTime - 0.5*spacing[0]);
}

inline tw::Float Simulation::ToLight(tw::Float z,tw::Float relativeTime)
{
	return z + windowPosition - antiWindowPosition - signalSpeed*(elapsedTime + relativeTime - 0.5*spacing[0]);
}

template <class T,class U>
U Simulation::ValueOnLabGrid(T& A,tw::strip s,tw::Int k,tw::Float relativeTime)
{
	// Take a quantity known on the light grid and get its value in a cell of the lab grid
	tw::Int klight;
	tw::Float z,zeta,w;
	z = Pos(s,k).z - corner[3];
	zeta = ToLight(z,relativeTime);
	klight = MyFloor(zeta*freq[3] + 0.5001);
	w = 0.5 - zeta*freq[3] + tw::Float(klight);
	return w*A(s,klight) + (one - w)*A(s,klight+1);
}

template <class T,class U>
U Simulation::ValueOnLightGrid(T& A,tw::strip s,tw::Int k,tw::Float relativeTime)
{
	// Take a quantity known on the lab grid and get its value in a cell of the light grid
	tw::Int klab;
	tw::Float z,zeta,w;
	zeta = Pos(s,k).z - corner[3];
	z = ToLab(zeta,relativeTime);
	klab = MyFloor(z*freq[3] + 0.4999);
	w = 0.5 - z*freq[3] + tw::Float(klab);
	return w*A(s,klab) + (one - w)*A(s,klab+1);
}
