module;

#include "tw_includes.h"

export module laser_solve;
import input;
import twmodule;
import compute_tool;
import parabolic;
import fields;
import diagnostics;
import numerics;
import logger;
#include "tw_logger.h"

using namespace tw::bc;

export struct LaserSolver:Module
{
	tw::Float laserFreq;
	tw::Int resolution;
	tw_polarization_type polarizationType;
	ComplexField a0,a1; // vector potential
	ComplexField chi; // defined by j = chi*a

	// fields that are tracked at higher resolution than the baseline grid
	MetricSpace HRSpace;
	ComplexField HRa0,HRa1;
	ComplexField HRchi;
	LaserPropagator *propagator;
	BoxDiagnostic *HRBoxDiagnostic;
	GlobalSpline<tw::Complex> *spliner;

	bool debug; // usually used to suppress envelope evolution

	LaserSolver(const std::string& name,Simulation* sim);
	virtual ~LaserSolver();
	virtual void ExchangeResources();
	virtual void Initialize();
	virtual void Reset();

	virtual void VerifyInput();
	virtual void ReadCheckpoint(std::ifstream& inFile);
	virtual void WriteCheckpoint(std::ofstream& outFile);

	void Downsample(const ComplexField& hiRes,ComplexField& loRes);
	void Upsample(ComplexField& hiRes,ComplexField& loRes);
	virtual void Update();
	tw::vec3 GetIonizationKick(const tw::Float& a2,const tw::Float& q0,const tw::Float& m0);
};

export struct QSSolver:LaserSolver
{
	QSSolver(const std::string& name,Simulation* sim);
};

export struct PGCSolver:LaserSolver
{
	Field F;

	PGCSolver(const std::string& name,Simulation* sim);
	virtual void ExchangeResources();
	virtual void Initialize();

	virtual void MoveWindow();
	virtual void AntiMoveWindow();

	virtual void Update();
	virtual void ComputeFinalFields();

	virtual void Report(Diagnostic&);
};

//////////////////////////////
//                          //
//     LASER SOLVER BASE    //
//                          //
//////////////////////////////

LaserSolver::LaserSolver(const std::string& name,Simulation* sim):Module(name,sim)
{
	if (native.native!=tw::units::plasma)
		throw tw::FatalError("LaserSolver module requires <native units = plasma>");

	updateSequencePriority = tw::priority::field;
	laserFreq = 10.0;
	polarizationType = linearPolarization;
	propagator = NULL;
	HRBoxDiagnostic = NULL;
	debug = false;
	resolution = 1;

	a0.Initialize(*this,owner);
	a1.Initialize(*this,owner);
	chi.Initialize(*this,owner);

	spliner = new GlobalSpline<tw::Complex>(&sim->strip[3],Num(1)*Num(2),Dim(3));

	directives.Add("carrier frequency",new tw::input::Float(&laserFreq));
	std::map<std::string,tw_polarization_type> pol = {{"linear",linearPolarization},{"circular",circularPolarization},{"radial",radialPolarization}};
	directives.Add("polarization",new tw::input::Enums<tw_polarization_type>(pol,&polarizationType),false);
	directives.Add("debug",new tw::input::Bool(&debug),false);
	directives.Add("resolution",new tw::input::Int(&resolution),false);
}

LaserSolver::~LaserSolver()
{
	if (propagator!=NULL)
		owner->RemoveTool(propagator);
	if (HRBoxDiagnostic!=NULL)
		owner->RemoveTool(HRBoxDiagnostic);
	delete spliner;
}

void LaserSolver::ExchangeResources()
{
	PublishResource(&a0,"laser:a0");
	PublishResource(&a1,"laser:a1");
	PublishResource(&chi,"laser:chi");
	PublishResource(&laserFreq,"laser:carrierFrequency");
	PublishResource(&polarizationType,"laser:polarizationType");
}

void LaserSolver::Initialize()
{
	logger::TRACE("initialize laser base");

	tw::vec3 pos;
	tw::Float polarizationFactor;
	const tw::Float dt = dx(0);
	const tw::Float dth = 0.5*dt;

	Module::Initialize();

	for (auto tool : moduleTool) {
		auto diag = dynamic_cast<BoxDiagnostic*>(tool);
		if (diag!=NULL && diag!=HRBoxDiagnostic) {
			diag->no_reports.push_back("a_real");
			diag->no_reports.push_back("a_imag");
			diag->no_reports.push_back("j1_real");
			diag->no_reports.push_back("j1_imag");
			HRBoxDiagnostic->CopyParams(*diag);
		}
	}

	HRa0.Initialize(HRSpace,owner);
	HRa1.Initialize(HRSpace,owner);
	HRchi.Initialize(HRSpace,owner);

	propagator->SetData(laserFreq,dt,polarizationType,owner->movingWindow,&HRSpace);
	propagator->SetBoundaryConditions(HRa0,HRa1,HRchi);
	propagator->SetBoundaryConditions(a0,a1,chi);

	if (polarizationType==circularPolarization) {
		polarizationFactor = 1.414;
	} else {
		polarizationFactor = 1.0;
	}

	for (auto cell : EntireCellRange(HRSpace)) {
		for (auto pulse : wave) {
			pos = HRSpace.Pos(cell);
			pos.z = HRSpace.ToLab(pos.z,-dth);
			HRa0(cell) += polarizationFactor*pulse->VectorPotentialEnvelope(-dth,pos,laserFreq);
			pos = HRSpace.Pos(cell);
			pos.z = HRSpace.ToLab(pos.z,dth);
			HRa1(cell) += polarizationFactor*pulse->VectorPotentialEnvelope(dth,pos,laserFreq);
		}
	}

	Downsample(HRa0,a0);
	Downsample(HRa1,a1);
}

tw::vec3 LaserSolver::GetIonizationKick(const tw::Float& a2,const tw::Float& q0,const tw::Float& m0)
{
	tw::Float phase;
	tw::vec3 ans;
	if (polarizationType==circularPolarization)
	{
		phase = owner->uniformDeviate->Next()*2.0*pi;
		// remember "a" has been multiplied by std::sqrt(2) at the beginning
		ans.x = q0*std::sqrt(0.5*a2)*std::cos(phase);
		ans.y = q0*std::sqrt(0.5*a2)*std::sin(phase);
		ans.z = 0.25*q0*q0*a2/m0;
	}
	else
	{
		// For linear polarization, assume phase is at zero of vector potential (peak of field).
		// Expression for pz is valid even for a>1, only caveat is this assumes the electron is
		// seeing plane wave fields during its first time level of freedom.
		ans.x = 0.0;
		ans.y = 0.0;
		ans.z = 0.25*q0*q0*a2/m0;
	}
	return ans;
}

void LaserSolver::Downsample(const ComplexField& hiRes,ComplexField& loRes) {
	if (resolution==1) {
		#pragma omp parallel
		{
			for (auto cell : EntireCellRange(loRes)) {
				loRes(cell,0) = hiRes(cell,0);
				loRes(cell,1) = hiRes(cell,1);
			}
		}
		return;
	}
	#pragma omp parallel
	{
		StripRange loRange(loRes,3,strongbool::yes);
		StripRange hiRange(hiRes,3,strongbool::yes);
		auto loStrip = loRange.begin();
		auto hiStrip = hiRange.begin();
		do {
			for (auto s=1;s<=loRes.Dim(3);s++) {
				// |       x       |
				// |   x   |   x   |
				// | x | x | x | x |
				tw::Int l = 1 + (s-1)*resolution;
				loRes(*loStrip,s,0) = 0.0;
				loRes(*loStrip,s,1) = 0.0;
				for (auto i=0;i<resolution;i++) {
					loRes(*loStrip,s,0) += hiRes(*hiStrip,l+i,0);
					loRes(*loStrip,s,1) += hiRes(*hiStrip,l+i,1);
				}
				loRes(*loStrip,s,0) /= resolution;
				loRes(*loStrip,s,1) /= resolution;
			}
			++loStrip;
			++hiStrip;
		} while (loStrip!=loRange.end() && hiStrip!=hiRange.end());
	}
	loRes.DownwardCopy(tw::grid::z,Element(0,1),1);
	loRes.UpwardCopy(tw::grid::z,Element(0,1),1);
	loRes.ApplyBoundaryCondition(Element(0,1));
}

void LaserSolver::Upsample(ComplexField& hiRes,ComplexField& loRes) {
	#pragma omp parallel
	{
		StripRange rng(loRes,3,strongbool::yes);
		for (auto it=rng.begin(); it!=rng.end(); ++it) {
			spliner->SetStrip(it.global_count(),&loRes(*it,0),loRes.Stride(3)/2);
		}
	}
	spliner->Solve();
	#pragma omp parallel
	{
		StripRange loRange(loRes,3,strongbool::yes);
		StripRange hiRange(hiRes,3,strongbool::yes);
		auto loStrip = loRange.begin();
		auto hiStrip = hiRange.begin();
		do {
			for (auto i=1;i<=hiRes.Dim(3);i++) {
				tw::Float x = 0.5 + tw::Float(i-0.5)/resolution;
				hiRes(*hiStrip,i) = spliner->Interpolate(x,loStrip.global_count(),&loRes(*loStrip,0),loRes.Stride(3)/2);
			}
			++loStrip;
			++hiStrip;
		} while (loStrip!=loRange.end() && hiStrip!=hiRange.end());
	}
	hiRes.DownwardCopy(tw::grid::z,Element(0,1),1);
	hiRes.UpwardCopy(tw::grid::z,Element(0,1),1);
	hiRes.ApplyBoundaryCondition(Element(0,1));
}

void LaserSolver::Update()
{
	logger::TRACE("start laser update");
	if (!debug) {
		Upsample(HRchi,chi);
		propagator->Advance(HRa0,HRa1,HRchi);
		Downsample(HRa0,a0);
		Downsample(HRa1,a1);
	}
}

void LaserSolver::Reset()
{
	chi = tw::Complex(0,0);
}

void LaserSolver::VerifyInput()
{
	Module::VerifyInput();
	for (auto tool : moduleTool) {
		propagator = dynamic_cast<LaserPropagator*>(tool);
		if (propagator!=NULL)
			break;
	}
	if (propagator==NULL) {
		propagator = (LaserPropagator*)owner->CreateTool("default_adi",tw::tool_type::adiPropagator);
		moduleTool.push_back(propagator);
	}
	tw::Int HRGlobalCells[4];
	HRGlobalCells[0] = owner->GlobalDim(0); 
	HRGlobalCells[1] = owner->GlobalDim(1); 
	HRGlobalCells[2] = owner->GlobalDim(2); 
	HRGlobalCells[3] = owner->GlobalDim(3)*resolution;
	logger::DEBUG(std::format("creating high resolution space x{}",resolution));
	HRSpace.Resize(owner,HRGlobalCells,owner->GlobalCorner(),owner->GlobalPhysicalSize(),owner->Layers(3),owner->gridGeometry);
	HRBoxDiagnostic = (BoxDiagnostic*)owner->CreateTool("hr_box",tw::tool_type::boxDiagnostic);
	moduleTool.push_back(HRBoxDiagnostic);
	HRBoxDiagnostic->filename = "refined";
	HRBoxDiagnostic->space = &HRSpace;
	HRBoxDiagnostic->reports.push_back("a_real");
	HRBoxDiagnostic->reports.push_back("a_imag");
	HRBoxDiagnostic->reports.push_back("j1_real");
	HRBoxDiagnostic->reports.push_back("j1_imag");
}

void LaserSolver::ReadCheckpoint(std::ifstream& inFile)
{
	Module::ReadCheckpoint(inFile);
	a0.ReadCheckpoint(inFile);
	a1.ReadCheckpoint(inFile);
	HRa0.ReadCheckpoint(inFile);
	HRa1.ReadCheckpoint(inFile);
}

void LaserSolver::WriteCheckpoint(std::ofstream& outFile)
{
	Module::WriteCheckpoint(outFile);
	a0.WriteCheckpoint(outFile);
	a1.WriteCheckpoint(outFile);
	HRa0.WriteCheckpoint(outFile);
	HRa1.WriteCheckpoint(outFile);
}


//////////////////////////////
//                          //
// Quasistatic Laser Solver //
//                          //
//////////////////////////////


QSSolver::QSSolver(const std::string& name,Simulation* sim):LaserSolver(name,sim)
{
}


//////////////////////////////
//                          //
//     PGC LASER SOLVER     //
//                          //
//////////////////////////////


PGCSolver::PGCSolver(const std::string& name,Simulation* sim):LaserSolver(name,sim)
{
	F.Initialize(8,*this,owner);
}

void PGCSolver::ExchangeResources()
{
	LaserSolver::ExchangeResources();
	PublishResource(&F,"laser:F");
}

void PGCSolver::Initialize()
{
	LaserSolver::Initialize();
	logger::TRACE("initialize pgc");

	// Here we deal with boundary conditions particular to PGC
	// B.C.'s for laser amplitude and sources are dealt with in propagator objects

	// longitudinal boundary conditions
	if (owner->movingWindow)
	{
		F.SetBoundaryConditions(Element(0,5),tw::grid::z,fld::neumannWall,fld::dirichletWall);
	}
	else
	{
		F.SetBoundaryConditions(Element(0,5),tw::grid::z,fld::neumannWall,fld::neumannWall);
		F.SetBoundaryConditions(Element(2),tw::grid::z,fld::dirichletWall,fld::dirichletWall);
		F.SetBoundaryConditions(Element(5),tw::grid::z,fld::dirichletWall,fld::dirichletWall);
	}
	F.SetBoundaryConditions(Element(6,7),tw::grid::z,fld::neumannWall,fld::neumannWall);

	// transverse boundary conditions
	if (owner->bc0[1]==par::axisymmetric || owner->bc0[1]==par::reflecting)
	{
		F.SetBoundaryConditions(Element(0,5),tw::grid::x, fld::neumannWall, fld::neumannWall);
		F.SetBoundaryConditions(Element(0),tw::grid::x,fld::dirichletWall,fld::neumannWall);
		F.SetBoundaryConditions(Element(3),tw::grid::x,fld::dirichletWall,fld::neumannWall);
	}
	else
	{
		F.SetBoundaryConditions(Element(0,5),tw::grid::x, fld::neumannWall, fld::neumannWall);
	}
	if (owner->bc0[2]==par::axisymmetric || owner->bc0[2]==par::reflecting)
	{
		F.SetBoundaryConditions(Element(0,5),tw::grid::y, fld::neumannWall, fld::neumannWall);
		F.SetBoundaryConditions(Element(1),tw::grid::y,fld::dirichletWall,fld::neumannWall);
		F.SetBoundaryConditions(Element(4),tw::grid::y,fld::dirichletWall,fld::neumannWall);
	}
	else
	{
		F.SetBoundaryConditions(Element(0,5),tw::grid::y, fld::neumannWall, fld::neumannWall);
	}
	F.SetBoundaryConditions(Element(6,7),tw::grid::x,fld::neumannWall,fld::neumannWall);
	F.SetBoundaryConditions(Element(6,7),tw::grid::y,fld::neumannWall,fld::neumannWall);

	ComputeFinalFields();
}

void PGCSolver::MoveWindow()
{
	LaserSolver::MoveWindow();
	logger::TRACE("field shift");
	for (auto s : StripRange(*this,3,strongbool::yes))
		F.Shift(s,-1,0.0);
	F.DownwardCopy(tw::grid::z,1);
}

void PGCSolver::AntiMoveWindow()
{
	const tw::Float dth = 0.5*dx(0);
	for (auto s : StripRange(*this,3,strongbool::yes))
	{
		tw::Float polarizationFactor = polarizationType==circularPolarization ? 1.414 : 1.0;
		tw::Complex incoming0(0,0);
		tw::Complex incoming1(0,0);
		for (auto pulse : wave)
		{
			tw::vec3 pos = owner->Pos(s,0);
			pos.z = owner->ToLab(pos.z,-dth);
			incoming0 += polarizationFactor*pulse->VectorPotentialEnvelope(owner->WindowPos(0)-dth,pos,laserFreq);
			pos = owner->Pos(s,0);
			pos.z = owner->ToLab(pos.z,dth);
			incoming1 += polarizationFactor*pulse->VectorPotentialEnvelope(owner->WindowPos(0)+dth,pos,laserFreq);
		}
		a0.Shift(s,1,incoming0);
		a1.Shift(s,1,incoming1);
	}
	a0.UpwardCopy(tw::grid::z,1);
	a1.UpwardCopy(tw::grid::z,1);
}

void PGCSolver::Update()
{
	logger::TRACE("start PGC update");
	chi.DepositFromNeighbors();
	chi.ApplyFoldingCondition();
	chi.DivideCellVolume(*owner);
	chi.ApplyBoundaryCondition();
	chi.Smooth(*owner,smoothing,compensation);
	#pragma omp parallel
	{
		const tw::Float dth = 0.5*dx(0);
		for (auto s : StripRange(*this,3,strongbool::yes))
		{
			for (tw::Int k=1;k<=dim[3];k++)
				chi(s,k) = owner->ValueOnLightGrid<ComplexField,tw::Complex>(chi,s,k,dth);
		}
	}
	chi.DownwardCopy(tw::grid::z,1);
	chi.UpwardCopy(tw::grid::z,1);
	chi.ApplyBoundaryCondition();

	LaserSolver::Update();
	ComputeFinalFields();
}

void PGCSolver::ComputeFinalFields()
{
	logger::TRACE("compute PGC forces");
	#pragma omp parallel
	{
		const tw::Float dth = 0.5*dx(0);
		for (auto s : StripRange(*this,3,strongbool::yes))
		{
			for (tw::Int k=1;k<=dim[3];k++)
			{
				F(s,k,7) = norm(owner->ValueOnLabGrid<ComplexField,tw::Complex>(a1,s,k,dth));
				F(s,k,6) = norm(owner->ValueOnLabGrid<ComplexField,tw::Complex>(a0,s,k,-dth));
				F(s,k,6) = 0.5*(F(s,k,6) + F(s,k,7));
			}
		}
	}

	F.DownwardCopy(tw::grid::z,Element(6,7),1);
	F.UpwardCopy(tw::grid::z,Element(6,7),1);

	#pragma omp parallel
	{
		for (auto cell : InteriorCellRange(*this))
		{
			F(cell,0) = F(cell,6,1);
			F(cell,1) = F(cell,6,2);
			F(cell,2) = F(cell,6,3);

			F(cell,3) = F(cell,7,1);
			F(cell,4) = F(cell,7,2);
			F(cell,5) = F(cell,7,3);
		}
	}

	F.CopyFromNeighbors(Element(0,5));
	F.ApplyBoundaryCondition();
}

void PGCSolver::Report(Diagnostic& diagnostic)
{
	LaserSolver::Report(diagnostic);

	if (diagnostic.name == "hr_box") {
		diagnostic.ReportField("a_real",HRa1,0,tw::dims::vector_potential,"$\\Re A$");
		diagnostic.ReportField("a_imag",HRa1,1,tw::dims::vector_potential,"$\\Im A$");
		diagnostic.ReportField("j1_real",HRchi,0,tw::dims::current_density,"$\\Re j$");
		diagnostic.ReportField("j1_imag",HRchi,1,tw::dims::current_density,"$\\Im j$");
		return;
	}

	ScalarField temp;
	temp.Initialize(*this,owner);

	const tw::Float dti = dk(0);
	const tw::Float dth = 0.5*dx(0);
	for (auto cell : InteriorCellRange(*this))
	{
		const tw::Complex aNow = half*(a0(cell)+a1(cell));
		const tw::Complex dtau = dti*(a1(cell)-a0(cell));
		const tw::Complex dzeta = half*(a0(cell,0,3) + a1(cell,0,3)) + ii*half*(a0(cell,1,3) + a1(cell,1,3));
		const tw::Complex eNow = ii*laserFreq*aNow - (dtau-dzeta);
		const tw::Complex bNow = ii*laserFreq*aNow + dzeta;
		temp(cell) = 0.25*(norm(eNow) + norm(bNow));
	}
	diagnostic.VolumeIntegral("LaserEnergy",temp,0);

	for (auto cell : InteriorCellRange(*this))
	{
		const tw::Complex aNow = half*(a0(cell)+a1(cell));
		//const tw::Complex dtau = dti*(a1(cell)-a0(cell));
		const tw::Complex dzeta = half*(a0(cell,0,3) + a1(cell,0,3)) + ii*half*(a0(cell,1,3) + a1(cell,1,3));
		//const tw::Complex eNow = ii*laserFreq*aNow - (dtau-dzeta);
		const tw::Complex bNow = ii*laserFreq*aNow + dzeta;
		temp(cell) = imag( conj(aNow)*bNow - aNow*conj(bNow) );
	}
	diagnostic.VolumeIntegral("WaveAction",temp,0);

	for (auto s : StripRange(*this,3,strongbool::no))
	{
		for (tw::Int k=1;k<=dim[3];k++)
		{
			tw::Complex dadt = dti*(owner->ValueOnLabGrid<ComplexField,tw::Complex>(a1,s,k,dth) - owner->ValueOnLabGrid<ComplexField,tw::Complex>(a0,s,k,-dth));
			tw::Complex anow = tw::Float(0.5)*(owner->ValueOnLabGrid<ComplexField,tw::Complex>(a0,s,k,-dth) + owner->ValueOnLabGrid<ComplexField,tw::Complex>(a1,s,k,dth));
			temp(s,k) = -real(dadt - ii*laserFreq*anow);
		}
	}
	diagnostic.ReportField("e_real",temp,0,tw::dims::electric_field,"$\\Re E$");

	for (auto s : StripRange(*this,3,strongbool::no))
	{
		for (tw::Int k=1;k<=dim[3];k++)
		{
			tw::Complex dadt = dti*(owner->ValueOnLabGrid<ComplexField,tw::Complex>(a1,s,k,dth) - owner->ValueOnLabGrid<ComplexField,tw::Complex>(a0,s,k,-dth));
			tw::Complex anow = tw::Float(0.5)*(owner->ValueOnLabGrid<ComplexField,tw::Complex>(a0,s,k,-dth) + owner->ValueOnLabGrid<ComplexField,tw::Complex>(a1,s,k,dth));
			temp(s,k) = -imag(dadt - ii*laserFreq*anow);
		}
	}
	diagnostic.ReportField("e_imag",temp,0,tw::dims::electric_field,"$\\Im E$");

	diagnostic.ReportField("a2",F,7,tw::dims::none,"$a^2$");
	diagnostic.ReportField("chi_real",chi,0,tw::dims::none,"$\\Re \\chi$");
	diagnostic.ReportField("chi_imag",chi,1,tw::dims::none,"$\\Im \\chi$");
}
