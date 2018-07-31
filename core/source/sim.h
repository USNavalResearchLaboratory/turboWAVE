#include "definitions.h"
#include "tasks.h"
#include "ctools.h"
#include "3dmath.h"
#include "metricSpace.h"
#include "3dfields.h"
#include "region.h"
#include "numerics.h"
#include "computeTool.h"
#include "parabolic.h"
#include "elliptic.h"
#include "hyperbolic.h"
#include "fft.h"
#include "functions.h"
#include "injection.h"
#include "physics.h"

struct Grid;

#include "input.h"
#include "diagnostics.h"
#include "module.h"

void ReduceInputFile(std::ifstream& inputFile,std::stringstream& out);

void PreprocessInputFile(std::ifstream& inputFile,std::stringstream& out);

void ExitInputFileBlock(std::stringstream& inputString);

struct NonUniformRegion
{
	tw::Int i1,i2,ih,N;
	tw::Float L,gridSum,dz,aux;

	NonUniformRegion(tw::Int first,tw::Int last,tw::Float length,tw::Float dz0);
	tw::Float AddedCellWidth(tw::Int globalCell);
	tw::Float ACoefficient(tw::Float length);
};

struct Grid:Task,MetricSpace
{
	std::ostream *tw_out;

	std::vector<NonUniformRegion*> region;
	tw::Float radialProgressionFactor;

	tw_geometry gridGeometry;
	tw::Float dt0,dt,dth,dtMin,dtMax,elapsedTime,elapsedTimeMax;
	tw::Float signalPosition,windowPosition,signalSpeed;
	tw::Float antiSignalPosition,antiWindowPosition;
	tw::Float unitDensityCGS;

	tw::Int smoothing,compensation;
	bool neutralize,movingWindow,adaptiveTimestep,adaptiveGrid;
	bool restarted,appendMode,fullOutput,completed;

	tw::Int stepsToTake,stepNow;
	tw::Int lastTime;

	tw::Int dumpPeriod,sortPeriod,sortX,sortY,sortZ;
	tw::Int binaryFormat;

	tw_boundary_spec bc0[4],bc1[4];

	std::vector<ComputeTool*> computeTool;
	std::vector<Module*> module;

	std::vector<Region*> clippingRegion;
	std::vector<Wave*> wave;
	std::vector<Pulse*> pulse;
	std::vector<Conductor*> conductor;

	std::vector<EnergySeriesDescriptor*> energyDiagnostic;
	std::vector<PointSeriesDescriptor*> pointDiagnostic;
	std::vector<GridDataDescriptor*> boxDiagnostic;

	UniformDeviate *uniformDeviate;
	GaussianDeviate *gaussianDeviate;

	#ifdef USE_OPENCL
	cl_mem waveBuffer;
	#endif

	Grid();
	virtual ~Grid();
	virtual void Run();
	void SetupGeometry();
	void PrepareSimulation();
	void FundamentalCycle();
	void MoveWindow();
	void AntiMoveWindow();
	void Diagnose();
	void EmergencyDump();
	void GetGlobalBoxDataIndexing(GridDataDescriptor* theBox,tw::Int pts[4],tw::Int glob[6],tw::Int skip[4]);
	void GetLocalBoxDataIndexing(GridDataDescriptor* theBox,const tw::Int pts[4],const tw::Int glob[6],const tw::Int skip[4],tw::Int loc[6],const tw::Int coords[4]);
	void WriteBoxDataHeader(const std::string& quantity,GridDataDescriptor* theBox);
	void WriteMomentumDataHeader(const std::string& quantity,GridDataDescriptor* theBox);
	void WriteBoxData(const std::string& quantity,GridDataDescriptor* theBox,tw::Float* theData,const tw::Int *stride);
	void WriteMomentumData(const std::string& quantity,GridDataDescriptor* theBox,tw::Float* theData,const tw::Int *stride);
	void WriteCellDataHeader(GridDataDescriptor* theBox);
	void WriteCellData(GridDataDescriptor* theBox);

	ComputeTool* AddPrivateTool(tw_tool whichTool);
	ComputeTool* AddSharedTool(tw_tool whichTool);
	bool RemoveTool(ComputeTool *theTool);

	void SetCellWidthsAndLocalSize();
	void SetGlobalSizeAndLocalCorner();

	tw::Float ToLab(tw::Float zeta,tw::Float relativeTime);
	tw::Float ToLight(tw::Float z,tw::Float relativeTime);
	template <class T,class U>
	U ValueOnLabGrid(T& A,tw::Int i,tw::Int j,tw::Int k,tw::Float relativeTime);
	template <class T,class U>
	U ValueOnLightGrid(T& A,tw::Int i,tw::Int j,tw::Int k,tw::Float relativeTime);

	void OpenInputFile(std::ifstream& inFile);
	std::string InputFileFirstPass();
	void GridFromInputFile();
	void ReadInputFile();
	void ReadData(std::ifstream& inFile);
	void WriteData(std::ofstream& outFile);
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

	void SetupTimeInfo(tw::Float dt);

	bool IsFirstStep()
	{
		return stepNow == 1;
	}

	bool IsLastStep()
	{
		return stepNow == stepsToTake;
	}

	bool Completed()
	{
		return completed;
	}
};

inline tw::Float Grid::ToLab(tw::Float zeta,tw::Float relativeTime)
{
	return zeta + antiWindowPosition - windowPosition + signalSpeed*(elapsedTime + relativeTime - dth);
}

inline tw::Float Grid::ToLight(tw::Float z,tw::Float relativeTime)
{
	return z + windowPosition - antiWindowPosition - signalSpeed*(elapsedTime + relativeTime - dth);
}

template <class T,class U>
U Grid::ValueOnLabGrid(T& A,tw::Int i,tw::Int j,tw::Int k,tw::Float relativeTime)
{
	// Take a quantity known on the light grid and get its value in a cell of the lab grid
	tw::Int cell;
	tw::Float z,zeta,w;
	z = Pos(i,j,k).z - corner.z;
	zeta = ToLight(z,relativeTime);
	cell = MyFloor(zeta*freq.z + 0.5001);
	w = 0.5 - zeta*freq.z + tw::Float(cell);
	return w*A(i,j,cell) + (one - w)*A(i,j,cell+1);
}

template <class T,class U>
U Grid::ValueOnLightGrid(T& A,tw::Int i,tw::Int j,tw::Int k,tw::Float relativeTime)
{
	// Take a quantity known on the lab grid and get its value in a cell of the light grid
	tw::Int cell;
	tw::Float z,zeta,w;
	zeta = Pos(i,j,k).z - corner.z;
	z = ToLab(zeta,relativeTime);
	cell = MyFloor(z*freq.z + 0.4999);
	w = 0.5 - z*freq.z + tw::Float(cell);
	return w*A(i,j,cell) + (one - w)*A(i,j,cell+1);
}
