module;

#include "tw_includes.h"
#include "tw_logger.h"

export module laser_solve;
import input;
import twmodule;
import compute_tool;
import parabolic;
import fields;
import diagnostics;
import numerics;
import injection;
import logger;

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

	tw::Waves waves;
	std::shared_ptr<LaserPropagator> propagator;
	std::shared_ptr<BoxDiagnostic> HRBoxDiagnostic;
	std::unique_ptr<GlobalSpline<tw::Complex>> spliner;

	bool debug; // usually used to suppress envelope evolution

	LaserSolver(const std::string& name,Simulation* sim);
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
	debug = false;
	resolution = 1;

	a0.Initialize(*this,owner);
	a1.Initialize(*this,owner);
	chi.Initialize(*this,owner);

	spliner = std::make_unique<GlobalSpline<tw::Complex>>(&sim->strip[3],Num(1)*Num(2),Dim(3));

	directives.Add("carrier frequency",new tw::input::Float(&laserFreq));
	std::map<std::string,tw_polarization_type> pol = {{"linear",linearPolarization},{"circular",circularPolarization},{"radial",radialPolarization}};
	directives.Add("polarization",new tw::input::Enums<tw_polarization_type>(pol,&polarizationType),false);
	directives.Add("debug",new tw::input::Bool(&debug),false);
	directives.Add("resolution",new tw::input::Int(&resolution),false);
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

	for (auto tool : tools) {
		auto diag = std::dynamic_pointer_cast<BoxDiagnostic>(tool);
		if (diag && diag!=HRBoxDiagnostic) {
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

	for (auto cell : EntireCellRange(HRSpace,1)) {
		for (auto pulse : waves) {
			pos = HRSpace.Pos(cell);
			pos.z = HRSpace.ToLab(pos.z,-dth);
			HRa0.Pack(cell, HRa0(cell) + polarizationFactor*pulse->VectorPotentialEnvelope(-dth,pos,laserFreq));
			pos = HRSpace.Pos(cell);
			pos.z = HRSpace.ToLab(pos.z,dth);
			HRa1.Pack(cell, HRa1(cell) + polarizationFactor*pulse->VectorPotentialEnvelope(dth,pos,laserFreq));
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
			for (auto cell : EntireCellRange(loRes,1)) {
				loRes(cell,0) = hiRes(cell,0);
				loRes(cell,1) = hiRes(cell,1);
			}
		}
		return;
	}
	#pragma omp parallel
	{
		StripRange loRange(loRes,3,0,1,strongbool::yes);
		StripRange hiRange(hiRes,3,0,1,strongbool::yes);
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
	loRes.DownwardCopy(tw::grid::z,1);
	loRes.UpwardCopy(tw::grid::z,1);
	loRes.ApplyBoundaryCondition();
}

void LaserSolver::Upsample(ComplexField& hiRes,ComplexField& loRes) {
	#pragma omp parallel
	{
		StripRange rng(loRes,3,0,1,strongbool::yes);
		for (auto it=rng.begin(); it!=rng.end(); ++it) {
			spliner->SetStrip(it.global_count(),&loRes(*it,0,0),loRes.Stride(3),loRes.Stride(4));
		}
	}
	spliner->Solve();
	#pragma omp parallel
	{
		StripRange loRange(loRes,3,0,1,strongbool::yes);
		StripRange hiRange(hiRes,3,0,1,strongbool::yes);
		auto loStrip = loRange.begin();
		auto hiStrip = hiRange.begin();
		do {
			for (auto i=1;i<=hiRes.Dim(3);i++) {
				tw::Float x = 0.5 + tw::Float(i-0.5)/resolution;
				hiRes.Pack(*hiStrip,i, spliner->Interpolate(x,loStrip.global_count()));
			}
			++loStrip;
			++hiStrip;
		} while (loStrip!=loRange.end() && hiStrip!=hiRange.end());
	}
	hiRes.DownwardCopy(tw::grid::z,1);
	hiRes.UpwardCopy(tw::grid::z,1);
	hiRes.ApplyBoundaryCondition();
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
	for (auto tool : tools) {
		if (std::dynamic_pointer_cast<LaserPropagator>(tool)) {
			propagator = std::dynamic_pointer_cast<LaserPropagator>(tool);
		} else if (std::dynamic_pointer_cast<Wave>(tool)) {
			waves.push_back(std::dynamic_pointer_cast<Wave>(tool));
		}
	}
	if (!propagator) {
		auto name = owner->CreateTool("default_adi",tw::tool_type::adiPropagator);
		propagator = std::dynamic_pointer_cast<LaserPropagator>(owner->UseTool(name));
	}
	tw::node5 HRGlobalCells {
		owner->GlobalDim(0),
		owner->GlobalDim(1),
		owner->GlobalDim(2),
		owner->GlobalDim(3)*resolution,
		1,
	};
	tw::node4 layers {
		owner->Layers(0),
		owner->Layers(1),
		owner->Layers(2),
		owner->Layers(3)
	};
	logger::DEBUG(std::format("creating high resolution space x{}",resolution));
	HRSpace.Resize(owner,HRGlobalCells,owner->GlobalCorner(),owner->GlobalPhysicalSize(),std_packing,layers,owner->gridGeometry);
	auto name = owner->CreateTool("hr_box",tw::tool_type::boxDiagnostic);
	HRBoxDiagnostic = std::dynamic_pointer_cast<BoxDiagnostic>(owner->UseTool(name));
	HRBoxDiagnostic->filename = "refined";
	HRBoxDiagnostic->space = &HRSpace;
	HRBoxDiagnostic->reports.push_back("a_real");
	HRBoxDiagnostic->reports.push_back("a_imag");
	HRBoxDiagnostic->reports.push_back("j1_real");
	HRBoxDiagnostic->reports.push_back("j1_imag");
	tools.push_back(HRBoxDiagnostic);
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
	F.Initialize(8,*this,sim);
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
		F.SetBoundaryConditions(Rng(0,6),tw::grid::z,fld::neumannWall,fld::dirichletWall);
	}
	else
	{
		F.SetBoundaryConditions(Rng(0,6),tw::grid::z,fld::neumannWall,fld::neumannWall);
		F.SetBoundaryConditions(Rng(2),tw::grid::z,fld::dirichletWall,fld::dirichletWall);
		F.SetBoundaryConditions(Rng(5),tw::grid::z,fld::dirichletWall,fld::dirichletWall);
	}
	F.SetBoundaryConditions(Rng(6,8),tw::grid::z,fld::neumannWall,fld::neumannWall);

	// transverse boundary conditions
	if (owner->bc0[1]==par::axisymmetric || owner->bc0[1]==par::reflecting)
	{
		F.SetBoundaryConditions(Rng(0,6),tw::grid::x, fld::neumannWall, fld::neumannWall);
		F.SetBoundaryConditions(Rng(0),tw::grid::x,fld::dirichletWall,fld::neumannWall);
		F.SetBoundaryConditions(Rng(3),tw::grid::x,fld::dirichletWall,fld::neumannWall);
	}
	else
	{
		F.SetBoundaryConditions(Rng(0,6),tw::grid::x, fld::neumannWall, fld::neumannWall);
	}
	if (owner->bc0[2]==par::axisymmetric || owner->bc0[2]==par::reflecting)
	{
		F.SetBoundaryConditions(Rng(0,6),tw::grid::y, fld::neumannWall, fld::neumannWall);
		F.SetBoundaryConditions(Rng(1),tw::grid::y,fld::dirichletWall,fld::neumannWall);
		F.SetBoundaryConditions(Rng(4),tw::grid::y,fld::dirichletWall,fld::neumannWall);
	}
	else
	{
		F.SetBoundaryConditions(Rng(0,6),tw::grid::y, fld::neumannWall, fld::neumannWall);
	}
	F.SetBoundaryConditions(Rng(6,8),tw::grid::x,fld::neumannWall,fld::neumannWall);
	F.SetBoundaryConditions(Rng(6,8),tw::grid::y,fld::neumannWall,fld::neumannWall);

	ComputeFinalFields();
}

void PGCSolver::MoveWindow()
{
	LaserSolver::MoveWindow();
	logger::TRACE("field shift");
	for (auto s : StripRange(*this,3,0,1,strongbool::yes))
		F.Shift(Rng(0,8),s,-1,0.0);
	F.DownwardCopy(Rng(0,8),tw::grid::z,1);
}

void PGCSolver::AntiMoveWindow()
{
	// assumes resolution = 1
	const tw::Float dth = 0.5*dx(0);
	for (auto s : StripRange(*this,3,0,1,strongbool::yes))
	{
		tw::Float polarizationFactor = polarizationType==circularPolarization ? 1.414 : 1.0;
		tw::Complex incoming0(0,0);
		tw::Complex incoming1(0,0);
		for (auto pulse : waves)
		{
			tw::vec3 pos = owner->Pos(s,0);
			pos.z = owner->ToLab(pos.z,-dth);
			incoming0 += polarizationFactor*pulse->VectorPotentialEnvelope(owner->WindowPos(0)-dth,pos,laserFreq);
			pos = owner->Pos(s,0);
			pos.z = owner->ToLab(pos.z,dth);
			incoming1 += polarizationFactor*pulse->VectorPotentialEnvelope(owner->WindowPos(0)+dth,pos,laserFreq);
		}
		HRa0.Shift(Rng(0,2),s,1,(tw::Float*)&incoming0);
		HRa1.Shift(Rng(0,2),s,1,(tw::Float*)&incoming1);
	}
	HRa0.UpwardCopy(tw::grid::z,1);
	HRa1.UpwardCopy(tw::grid::z,1);
}

void PGCSolver::Update()
{
	logger::TRACE("start PGC update");
	chi.DepositFromNeighbors();
	chi.ApplyFoldingCondition();
	chi.DivideCellVolume(Rng(0,2),*owner);
	chi.ApplyBoundaryCondition();
	chi.Smooth(Rng(0,2),*owner,smoothing,compensation);
	#pragma omp parallel
	{
		const tw::Float dth = 0.5*dx(0);
		for (auto s : StripRange(*this,3,0,1,strongbool::yes))
		{
			for (tw::Int k=1;k<=dim[3];k++)
				chi.Pack(s,k, owner->ValueOnLightGrid<ComplexField,tw::Complex>(chi,s,k,dth));
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
	logger::TRACE("compute PGC potentials");
	#pragma omp parallel
	{
		const tw::Float dth = 0.5*dx(0);
		for (auto s : StripRange(*this,3,0,1,strongbool::yes))
		{
			for (tw::Int k=1;k<=dim[3];k++)
			{
				F(s,k,7) = norm(owner->ValueOnLabGrid<ComplexField,tw::Complex>(a1,s,k,dth));
				F(s,k,6) = norm(owner->ValueOnLabGrid<ComplexField,tw::Complex>(a0,s,k,-dth));
				F(s,k,6) = 0.5*(F(s,k,6) + F(s,k,7));
			}
		}
	}

	F.DownwardCopy(Rng(6,8),tw::grid::z,1);
	F.UpwardCopy(Rng(6,8),tw::grid::z,1);

	logger::TRACE("compute PGC forces");
	#pragma omp parallel
	{
		for (auto cell : InteriorCellRange(*this,1))
		{
			F(cell,0) = F.d1(cell,6,1);
			F(cell,1) = F.d1(cell,6,2);
			F(cell,2) = F.d1(cell,6,3);

			F(cell,3) = F.d1(cell,7,1);
			F(cell,4) = F.d1(cell,7,2);
			F(cell,5) = F.d1(cell,7,3);
		}
	}

	F.CopyFromNeighbors(Rng(0,6));
	F.ApplyBoundaryCondition(Rng(0,6));
}

void PGCSolver::Report(Diagnostic& diagnostic)
{
	LaserSolver::Report(diagnostic);

	if (diagnostic.name == "hr_box") {
		diagnostic.ReportField("a_real",HRa1,1,0,tw::dims::vector_potential,"$\\Re A$");
		diagnostic.ReportField("a_imag",HRa1,1,1,tw::dims::vector_potential,"$\\Im A$");
		diagnostic.ReportField("j1_real",HRchi,1,0,tw::dims::current_density,"$\\Re j$");
		diagnostic.ReportField("j1_imag",HRchi,1,1,tw::dims::current_density,"$\\Im j$");
		return;
	}

	ScalarField temp;
	temp.Initialize(*this,owner);

	const tw::Float dti = dk(0);
	const tw::Float dth = 0.5*dx(0);
	for (auto cell : InteriorCellRange(*this,1))
	{
		const tw::Complex aNow = half*(a0(cell)+a1(cell));
		const tw::Complex dtau = dti*(a1(cell)-a0(cell));
		const tw::Complex dzeta = half*(a0.d1(cell,0,3) + a1.d1(cell,0,3)) + ii*half*(a0.d1(cell,1,3) + a1.d1(cell,1,3));
		const tw::Complex eNow = ii*laserFreq*aNow - (dtau-dzeta);
		const tw::Complex bNow = ii*laserFreq*aNow + dzeta;
		temp(cell) = 0.25*(norm(eNow) + norm(bNow));
	}
	diagnostic.VolumeIntegral("LaserEnergy",temp,1,0);

	for (auto cell : InteriorCellRange(*this,1))
	{
		const tw::Complex aNow = half*(a0(cell)+a1(cell));
		//const tw::Complex dtau = dti*(a1(cell)-a0(cell));
		const tw::Complex dzeta = half*(a0.d1(cell,0,3) + a1.d1(cell,0,3)) + ii*half*(a0.d1(cell,1,3) + a1.d1(cell,1,3));
		//const tw::Complex eNow = ii*laserFreq*aNow - (dtau-dzeta);
		const tw::Complex bNow = ii*laserFreq*aNow + dzeta;
		temp(cell) = imag( conj(aNow)*bNow - aNow*conj(bNow) );
	}
	diagnostic.VolumeIntegral("WaveAction",temp,1,0);

	for (auto s : StripRange(*this,3,0,1,strongbool::no))
	{
		for (tw::Int k=1;k<=dim[3];k++)
		{
			tw::Complex dadt = dti*(owner->ValueOnLabGrid<ComplexField,tw::Complex>(a1,s,k,dth) - owner->ValueOnLabGrid<ComplexField,tw::Complex>(a0,s,k,-dth));
			tw::Complex anow = tw::Float(0.5)*(owner->ValueOnLabGrid<ComplexField,tw::Complex>(a0,s,k,-dth) + owner->ValueOnLabGrid<ComplexField,tw::Complex>(a1,s,k,dth));
			temp(s,k) = -real(dadt - ii*laserFreq*anow);
		}
	}
	diagnostic.ReportField("e_real",temp,1,0,tw::dims::electric_field,"$\\Re E$");

	for (auto s : StripRange(*this,3,0,1,strongbool::no))
	{
		for (tw::Int k=1;k<=dim[3];k++)
		{
			tw::Complex dadt = dti*(owner->ValueOnLabGrid<ComplexField,tw::Complex>(a1,s,k,dth) - owner->ValueOnLabGrid<ComplexField,tw::Complex>(a0,s,k,-dth));
			tw::Complex anow = tw::Float(0.5)*(owner->ValueOnLabGrid<ComplexField,tw::Complex>(a0,s,k,-dth) + owner->ValueOnLabGrid<ComplexField,tw::Complex>(a1,s,k,dth));
			temp(s,k) = -imag(dadt - ii*laserFreq*anow);
		}
	}
	diagnostic.ReportField("e_imag",temp,1,0,tw::dims::electric_field,"$\\Im E$");

	diagnostic.ReportField("a2",F,1,7,tw::dims::none,"$a^2$");
	diagnostic.ReportField("chi_real",chi,1,0,tw::dims::none,"$\\Re \\chi$");
	diagnostic.ReportField("chi_imag",chi,1,1,tw::dims::none,"$\\Im \\chi$");
}
