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


Module::Module(const std::string& name,Grid* theGrid)
{
	this->name = name;
	owner = theGrid;
	super = NULL;
	programFilename = "";
	buildLog = "";
	dt = theGrid->dt;
	dth = theGrid->dth;
	dti = 1.0/dt;
	updateSequencePriority = 1;
	typeCode = tw::module_type::nullModule;
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

bool Module::ValidSubmodule(Module* sub)
{
	return false;
}

bool Module::AddSubmodule(Module* sub)
{
	if (ValidSubmodule(sub))
	{
		submodule.push_back(sub);
		sub->super = this;
		return true;
	}
	else
		return false;
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
	DiscreteSpace::ReadData(inFile);
	inFile.read((char *)&dt,sizeof(tw::Float));
	inFile.read((char *)&dth,sizeof(tw::Float));
	dti = 1.0/dt;

	tw::Int i,num;
	inFile.read((char *)&num,sizeof(tw::Int));
	for (i=0;i<num;i++)
		profile.push_back(Profile::CreateObjectFromFile(owner->clippingRegion,inFile));

	// Find supermodule and setup containment
	std::string super_name;
	inFile >> super_name;
	inFile.ignore();
	owner->GetModule(super_name)->AddSubmodule(this);
}

void Module::WriteData(std::ofstream& outFile)
{
	outFile.write((char *)&typeCode,sizeof(typeCode));
	outFile << name << " ";
	DiscreteSpace::WriteData(outFile);
	outFile.write((char *)&dt,sizeof(tw::Float));
	outFile.write((char *)&dth,sizeof(tw::Float));

	tw::Int i = profile.size();
	outFile.write((char *)&i,sizeof(tw::Int));
	for (i=0;i<profile.size();i++)
		profile[i]->WriteData(outFile);

	// Write the name of the supermodule
	// This can be used to reconstruct the whole hierarchy, assuming modules are sorted correctly
	outFile << super->name << " ";
}

void Module::ReadInputFileBlock(std::stringstream& inputString)
{
	std::string com;
	do
	{
		inputString >> com;
		ReadInputFileDirective(inputString,com);
	} while (com!="}");
}

bool Module::ReadQuasitoolBlock(const std::vector<std::string>& preamble,std::stringstream& inputString)
{
	return false;
}

void Module::ReadInputFileDirective(std::stringstream& inputString,const std::string& command)
{
}

void Module::StartDiagnostics()
{
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

////// STATIC MEMBERS FOLLOW

bool Module::SingularType(tw::module_type theType)
{
	return theType==tw::module_type::kinetics || theType==tw::module_type::chemistry;
}

tw::module_type Module::CreateSupermoduleTypeFromSubmoduleKey(const std::string& key)
{
	if (key=="species")
		return tw::module_type::kinetics;
	if (key=="chemical")
		return tw::module_type::equilibriumGroup;
	if (key=="group")
		return tw::module_type::chemistry;
	return tw::module_type::nullModule;
}

bool Module::QuasitoolNeedsModule(const std::vector<std::string>& preamble)
{
	bool ans = false;
	ans = ans || (preamble[0]=="phase");
	ans = ans || (preamble[0]=="orbit");
	ans = ans || (preamble[0]=="detector");
	ans = ans || (preamble[0]=="reaction");
	ans = ans || (preamble[0]=="collision");
	ans = ans || (preamble[0]=="excitation");
	return ans;
}

tw::module_type Module::CreateTypeFromInput(const std::vector<std::string>& preamble)
{
	if (preamble[0]=="curvilinear" && preamble[1]=="direct")
		return tw::module_type::curvilinearDirectSolver;
	if (preamble[0]=="coulomb")
		return tw::module_type::coulombSolver;
	if (preamble[0]=="electrostatic")
		return tw::module_type::electrostatic;
	if (preamble[0]=="direct")
		return tw::module_type::directSolver;
	if (preamble[0]=="quasistatic")
		return tw::module_type::qsLaser;
	if (preamble[0]=="pgc" || preamble[0]=="PGC")
		return tw::module_type::pgcLaser;
	if (preamble[0]=="bound")
		return tw::module_type::boundElectrons;
	if (preamble[0]=="schroedinger" || preamble[0]=="atomic")
		return tw::module_type::schroedinger;
	if (preamble[0]=="pauli")
		return tw::module_type::pauli;
	if (preamble[0]=="klein")
		return tw::module_type::kleinGordon;
	if (preamble[0]=="dirac")
		return tw::module_type::dirac;
	if (preamble[0]=="fluid")
		return tw::module_type::fluidFields;
	if (preamble[0]=="chemistry")
		return tw::module_type::chemistry;
	if (preamble[0]=="group")
		return tw::module_type::equilibriumGroup;
	if (preamble[0]=="chemical")
		return tw::module_type::chemical;
	if (preamble[0]=="species")
		return tw::module_type::species;
	if (preamble[0]=="kinetics")
		return tw::module_type::kinetics;
	return tw::module_type::nullModule;
}

Module* Module::CreateObjectFromType(const std::string& name,tw::module_type theType,Grid* theGrid)
{
	Module *ans;
	switch (theType)
	{
		case tw::module_type::nullModule:
			ans = NULL;
			break;
		case tw::module_type::curvilinearDirectSolver:
			ans = new CurvilinearDirectSolver(name,theGrid);
			break;
		case tw::module_type::electrostatic:
			ans = new Electrostatic(name,theGrid);
			break;
		case tw::module_type::coulombSolver:
			ans = new CoulombSolver(name,theGrid);
			break;
		case tw::module_type::directSolver:
			ans = new DirectSolver(name,theGrid);
			break;
		case tw::module_type::qsLaser:
			ans = new QSSolver(name,theGrid);
			break;
		case tw::module_type::pgcLaser:
			ans = new PGCSolver(name,theGrid);
			break;
		case tw::module_type::boundElectrons:
			ans = new BoundElectrons(name,theGrid);
			break;
		case tw::module_type::fluidFields:
			ans = new Fluid(name,theGrid);
			break;
		case tw::module_type::equilibriumGroup:
			ans = new EquilibriumGroup(name,theGrid);
			break;
		case tw::module_type::chemical:
			ans = new Chemical(name,theGrid);
			break;
		case tw::module_type::chemistry:
			ans = new Chemistry(name,theGrid);
			break;
		case tw::module_type::species:
			ans = new Species(name,theGrid);
			break;
		case tw::module_type::kinetics:
			ans = new Kinetics(name,theGrid);
			break;
		case tw::module_type::schroedinger:
			ans = new Schroedinger(name,theGrid);
			break;
		case tw::module_type::pauli:
			ans = new Pauli(name,theGrid);
			break;
		case tw::module_type::kleinGordon:
			ans = new KleinGordon(name,theGrid);
			break;
		case tw::module_type::dirac:
			ans = new Dirac(name,theGrid);
			break;
	}
	return ans;
}

Module* Module::CreateObjectFromFile(std::ifstream& inFile,Grid *theGrid)
{
	// This might work, but the containment heierarchy better be sorted.
	tw::module_type theType;
	std::string name;
	Module *ans;
	inFile.read((char*)&theType,sizeof(tw::module_type));
	inFile >> name;
	inFile.ignore();
	ans = CreateObjectFromType(name,theType,theGrid);
	ans->ReadData(inFile);
	return ans;
}
