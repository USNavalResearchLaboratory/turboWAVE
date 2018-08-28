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

#include "physics.h" // ASHER_MOD

//////////////////////////
//                      //
//  COMPUTE TOOL CLASS  //
//                      //
//////////////////////////


ComputeTool::ComputeTool(MetricSpace *ms,Task *tsk,bool shared)
{
	// Do not create ComputeTools directly.
	// Instances are created through Grid::AddSharedTool or Grid::AddPrivateTool
	space = ms;
	task = tsk;
	typeCode = nullTool;
	sharedTool = false;
	refCount = 0;
	name = "tw_tool";
	programFilename = "";
	buildLog = "";
}

ComputeTool::~ComputeTool()
{
	#ifdef USE_OPENCL

	if (programFilename!="")
		clReleaseProgram(program);

	#endif
}

void ComputeTool::Initialize()
{
	// Nothing to do.  Override if there is some initialization that depends
	// on the creation of all modules.
}

void ComputeTool::InitializeCLProgram(const std::string& filename)
{
	#ifdef USE_OPENCL
	task->InitializeCLProgram(program,filename,buildLog);
	programFilename = filename;
	#endif
}

ComputeTool* ComputeTool::CreateObjectFromType(tw_tool theType,MetricSpace *ms,Task *tsk,bool shared)
{
	// Called by Grid::AddSharedTool or Grid::AddPrivateTool
	ComputeTool *ans;
	switch (theType)
	{
		case nullTool:
			break;
		case eigenmodePropagator:
			ans = new EigenmodePropagator(ms,tsk,shared);
			break;
		case adiPropagator:
			ans = new ADIPropagator(ms,tsk,shared);
			break;
		//case conservativePropagator:
	//		ans = new ConservativePropagator(ms,tsk,shared);
	//		break;
		case isotropicPropagator:
			ans = new IsotropicPropagator(ms,tsk,shared);
			break;
		case generalParabolicPropagator:
			ans = new ParabolicSolver(ms,tsk,shared);
			break;
		case schroedingerPropagator:
			ans = new SchroedingerPropagator(ms,tsk,shared);
			break;
		case iterativePoissonSolver:
			ans = new IterativePoissonSolver(ms,tsk,shared);
			break;
		case ellipticSolver1D:
			ans = new EllipticSolver1D(ms,tsk,shared);
			break;
		case facrPoissonSolver:
			ans = new PoissonSolver(ms,tsk,shared);
			break;
		case eigenmodePoissonSolver:
			ans = new EigenmodePoissonSolver(ms,tsk,shared);
			break;
		case yeePropagatorPML:
			ans = new YeePropagatorPML(ms,tsk,shared);
			break;
		case lorentzPropagator:
			ans = new LorentzPropagator(ms,tsk,shared);
			break;
		case eosDataTool: // ASHER_MOD
			ans = new EOSDataTool(ms,tsk,shared);
			break;
		case eosIdealGas:
			ans = new EOSIdealGas(ms,tsk,shared);
			break;
		case eosHotElectrons:
			ans = new EOSHotElectrons(ms,tsk,shared);
			break;
		case eosMixture:
			ans = new EOSMixture(ms,tsk,shared);
			break;
		case eosIdealGasMix:
			ans = new EOSIdealGasMix(ms,tsk,shared);
			break;
		case eosMieGruneisen:
			ans = new EOSMieGruneisen(ms,tsk,shared);
			break;
		case eosMieGruneisen2:
			ans = new EOSMieGruneisen2(ms,tsk,shared);
			break;
	}
	return ans;
}

void ComputeTool::WarningMessage(std::ostream *theStream)
{
	// Save any messages from OpenCL compiler for later output
	if (buildLog.size()>4)
	{
		*theStream << "WARNING : Build log for " << programFilename << " is not empty:" << std::endl;
		*theStream << buildLog << std::endl;
	}
}
