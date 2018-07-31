#include "sim.h"
#include "fieldSolve.h"
#include "electrostatic.h"
#include "laserSolve.h"
#include "fluid.h"
#include "quantum.h"
#include "particles.h"
#include "solidState.h"


////////////////////////
//                    //
//    MODULE CLASS    //
//                    //
////////////////////////


Module::Module(Grid* theGrid)
{
	name = "Generic";
	owner = theGrid;
	programFilename = "";
	buildLog = "";
	dt = theGrid->dt;
	dth = theGrid->dth;
	dti = 1.0/dt;
	updateOrderingIndex = 1;
	typeCode = nullModule;
	suppressNextUpdate = false;
	DiscreteSpace::operator=(*theGrid);
}

Module::~Module()
{
	#ifdef USE_OPENCL

	if (programFilename!="")
		clReleaseProgram(program);

	#endif
}

void Module::InitializeCLProgram(const std::string& filename)
{
	#ifdef USE_OPENCL
	owner->InitializeCLProgram(program,filename,buildLog);
	programFilename = filename;
	#endif
}

void Module::PublishResource(void* resource,const std::string& description)
{
	tw::Int i;
	for (i=0;i<owner->module.size();i++)
		owner->module[i]->InspectResource(resource,description);
}

bool Module::InspectResource(void* resource,const std::string& description)
{
	return false;
}

void Module::Initialize()
{
	tw::Int i;
	if (!owner->restarted)
		for (i=0;i<profile.size();i++)
			profile[i]->Initialize(owner);
}

void Module::ReadData(std::ifstream& inFile)
{
	inFile >> name;
	inFile.ignore();
	DiscreteSpace::ReadData(inFile);
	inFile.read((char *)&dt,sizeof(tw::Float));
	inFile.read((char *)&dth,sizeof(tw::Float));
	dti = 1.0/dt;

	tw::Int i,num;
	inFile.read((char *)&num,sizeof(tw::Int));
	for (i=0;i<num;i++)
		profile.push_back(Profile::CreateObjectFromFile(owner->clippingRegion,inFile));
}

void Module::WriteData(std::ofstream& outFile)
{
	outFile.write((char *)&typeCode,sizeof(tw_module));
	outFile << name << " ";
	DiscreteSpace::WriteData(outFile);
	outFile.write((char *)&dt,sizeof(tw::Float));
	outFile.write((char *)&dth,sizeof(tw::Float));

	tw::Int i = profile.size();
	outFile.write((char *)&i,sizeof(tw::Int));
	for (i=0;i<profile.size();i++)
		profile[i]->WriteData(outFile);
}

void Module::ReadInputFileBlock(std::stringstream& inputString)
{
	std::string com;
	do
	{
		inputString >> com;
		ReadInputFileTerm(inputString,com);
	} while (com!="}");
}

void Module::ReadInputFileTerm(std::stringstream& inputString,std::string& command)
{
	std::string word;

	if (command=="name") // eg, name = He
	{
		inputString >> word;
		inputString >> name;
	}
}

void Module::StartDiagnostics()
{
}

Module* Module::CreateObjectFromFile(std::ifstream& inFile,Grid *theGrid)
{
	tw_module typeCode;
	Module *ans;
	inFile.read((char*)&typeCode,sizeof(tw_module));
	switch (typeCode)
	{
		case nullModule:
			break;
		case curvilinearDirectSolver:
			ans = new CurvilinearDirectSolver(theGrid);
			break;
		case electrostatic:
			ans = new Electrostatic(theGrid);
			break;
		case coulombSolver:
			ans = new CoulombSolver(theGrid);
			break;
		case directSolver:
			ans = new DirectSolver(theGrid);
			break;
		case qsLaser:
			ans = new QSSolver(theGrid);
			break;
		case pgcLaser:
			ans = new PGCSolver(theGrid);
			break;
		case boundElectrons:
			ans = new BoundElectrons(theGrid);
			break;
		case fluidFields:
			ans = new Fluid(theGrid);
			break;
		case equilibriumGroup:
			ans = new EquilibriumGroup(theGrid);
			break;
		case chemical:
			ans = new Chemical(theGrid);
			break;
		case chemistry:
			ans = new Chemistry(theGrid);
			break;
		case species:
			ans = new Species(theGrid);
			break;
		case kinetics:
			ans = new Kinetics(theGrid);
			break;
		case atomicPhysics:
			ans = new AtomicPhysics(theGrid);
			break;
		case schroedingerModule:
			ans = new Schroedinger(theGrid);
			break;
		case pauliModule:
			ans = new Pauli(theGrid);
			break;
		case kleinGordon:
			ans = new KleinGordon(theGrid);
			break;
		case dirac:
			ans = new Dirac(theGrid);
			break;
	}
	ans->ReadData(inFile);
	return ans;
}

void Module::WarningMessage(std::ostream *theStream)
{
	if (buildLog.size()>4)
	{
		*theStream << "WARNING : Build log for " << programFilename << " is not empty:" << std::endl;
		*theStream << buildLog << std::endl;
	}
	if (owner->movingWindow)
		for (tw::Int i=0;i<profile.size();i++)
			if (profile[i]->theRgn->moveWithWindow==true && profile[i]->theRgn->name!="entire")
				*theStream << "WARNING : region " << profile[i]->theRgn->name << " in motion in module " << name << std::endl;
}
