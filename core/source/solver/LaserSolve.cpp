module;

#include "tw_includes.h"
#include "tw_logger.h"

export module laser_solve;
import input;
import driver;
import parabolic;
import fields;
import diagnostics;
import numerics;
import injection;
import logger;

using namespace tw::bc;

export struct LaserSolver:Driver
{
	tw::Float laserFreq;
	tw::Int resolution;
	tw::Int N0; // time levels stored
	tw_polarization_type polarizationType;
	ComplexField a; // vector potential
	ComplexField chi; // defined by j = chi*a

	// fields that are tracked at higher resolution than the baseline grid
	MetricSpace HRSpace;
	ComplexField HRa;
	ComplexField HRchi;

	tw::Waves waves;
	std::shared_ptr<ForwardPropagator> propagator;
	std::unique_ptr<GlobalSpline<tw::Complex>> spliner;

	bool debug; // usually used to suppress envelope evolution

	LaserSolver(const std::string& name,MetricSpace *ms,Task *tsk);
	virtual void ExchangeResources();
	virtual void Initialize();
	virtual void Reset();

	virtual void VerifyInput();
	virtual void ReadCheckpoint(std::ifstream& inFile);
	virtual void WriteCheckpoint(std::ofstream& outFile);

	void Downsample(const ComplexField& hiRes,ComplexField& loRes,tw::Int n);
	void Upsample(ComplexField& hiRes,ComplexField& loRes,tw::Int n);
	virtual void Update();
	tw::vec3 GetIonizationKick(const tw::Float& a2,const tw::Float& q0,const tw::Float& m0);
};

export struct QSSolver:LaserSolver
{
	QSSolver(const std::string& name,MetricSpace *ms,Task *tsk);
};

export struct PGCSolver:LaserSolver
{
	Field F;

	PGCSolver(const std::string& name,MetricSpace *ms,Task *tsk);
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

LaserSolver::LaserSolver(const std::string& name,MetricSpace *ms,Task *tsk):Driver(name,ms,tsk)
{
	if (native.unit_system!=tw::units::plasma)
		throw tw::FatalError("LaserSolver module requires <native units = plasma>");

	updateSequencePriority = tw::priority::field;
	laserFreq = 10.0;
	polarizationType = linearPolarization;
	debug = false;
	resolution = 1;
	N0 = 3;

	a.Initialize(space->ax0(N0),task);
	chi.Initialize(space->ax0(N0),task);

	spliner = std::make_unique<GlobalSpline<tw::Complex>>(&task->strip[3],Num(1)*Num(2),Dim(3));

	directives.Add("carrier frequency",new tw::input::Float(&laserFreq));
	std::map<std::string,tw_polarization_type> pol = {{"linear",linearPolarization},{"circular",circularPolarization},{"radial",radialPolarization}};
	directives.Add("polarization",new tw::input::Enums<tw_polarization_type>(pol,&polarizationType),false);
	directives.Add("debug",new tw::input::Bool(&debug),false);
	directives.Add("resolution",new tw::input::Int(&resolution),false);
}

void LaserSolver::ExchangeResources()
{
	PublishResource(&chi,"laser:chi");
	PublishResource(&laserFreq,"laser:carrierFrequency");
	PublishResource(&polarizationType,"laser:polarizationType");
}

void LaserSolver::VerifyInput()
{
	Driver::VerifyInput();
	for (auto tool : tools) {
		if (std::dynamic_pointer_cast<ForwardPropagator>(tool)) {
			propagator = std::dynamic_pointer_cast<ForwardPropagator>(tool);
		} else if (std::dynamic_pointer_cast<Wave>(tool)) {
			waves.push_back(std::dynamic_pointer_cast<Wave>(tool));
		}
	}
	if (!propagator) {
		auto new_tool = CreateTool("default_fwd",tw::tool_type::forwardPropagator);
		AddTool(new_tool);
		propagator = std::dynamic_pointer_cast<ForwardPropagator>(new_tool);
	}
	tw::node5 HRGlobalCells {
		space->GlobalDim(0),
		space->GlobalDim(1),
		space->GlobalDim(2),
		space->GlobalDim(3)*resolution,
		1,
	};
	tw::node4 layers {
		space->Layers(0),
		space->Layers(1),
		space->Layers(2),
		space->Layers(3)
	};
	logger::DEBUG(std::format("creating high resolution space x{}",resolution));
	HRSpace.Resize(task,HRGlobalCells,space->GlobalCorner(),space->GlobalPhysicalSize(),std_packing,layers,space->gridGeometry);
}

void LaserSolver::Initialize()
{
	logger::DEBUG("initialize laser base");
	for (auto tool : tools) {
		if (std::dynamic_pointer_cast<BoxDiagnostic>(tool)) {
			logger::DEBUG(std::format("Add grid variant to {}",tool->name));
			std::dynamic_pointer_cast<BoxDiagnostic>(tool)->AddVariant(&HRSpace);
		}
	}

	tw::vec3 pos;
	tw::Float polarizationFactor;
	const tw::Float dt = dx(0);
	const tw::Float dth = 0.5*dt;

	Driver::Initialize();

	HRa.Initialize(HRSpace.ax0(N0),task);
	HRchi.Initialize(HRSpace.ax0(N0),task);

	propagator->SetData(laserFreq,dt,polarizationType,space->IsStdMovingWindow(),&HRSpace);
	propagator->SetBoundaryConditions(HRa,HRchi);
	propagator->SetBoundaryConditions(a,chi);

	if (polarizationType==circularPolarization) {
		polarizationFactor = 1.414;
	} else {
		polarizationFactor = 1.0;
	}

	for (auto n=1; n<=2; n++) {
		for (auto cell : EntireCellRange(HRa,n)) {
			for (auto pulse : waves) {
				pos = HRSpace.Pos(cell);
				pos.z = HRSpace.ToLab(pos.z,(1.5-n)*dt);
				HRa.Pack(cell, HRa(cell) + polarizationFactor*pulse->VectorPotentialEnvelope((1.5-n)*dt,pos,laserFreq));
			}
		}
		Downsample(HRa,a,n);
	}
}

tw::vec3 LaserSolver::GetIonizationKick(const tw::Float& a2,const tw::Float& q0,const tw::Float& m0)
{
	tw::Float phase;
	tw::vec3 ans;
	if (polarizationType==circularPolarization)
	{
		phase = task->uniformDeviate->Next()*2.0*pi;
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

void LaserSolver::Downsample(const ComplexField& hiRes,ComplexField& loRes,tw::Int n) {
	if (resolution==1) {
		#pragma omp parallel
		{
			for (auto cell : EntireCellRange(loRes,n)) {
				loRes(cell,0) = hiRes(cell,0);
				loRes(cell,1) = hiRes(cell,1);
			}
		}
		return;
	}
	#pragma omp parallel
	{
		StripRange loRange(loRes,3,0,n,strongbool::yes);
		StripRange hiRange(hiRes,3,0,n,strongbool::yes);
		for (auto loStrip=loRange.begin(),hiStrip=hiRange.begin(); loStrip!=loRange.end() && hiStrip!=hiRange.end(); ++loStrip,++hiStrip) {
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
		}
	}
	auto r = Rng04(n,n+1,0,2);
	loRes.Field::DownwardCopy(r,tw::grid::z,1);
	loRes.Field::UpwardCopy(r,tw::grid::z,1);
	loRes.Field::ApplyBoundaryCondition(r);
}

void LaserSolver::Upsample(ComplexField& hiRes,ComplexField& loRes,tw::Int n) {
	#pragma omp parallel
	{
		StripRange rng(loRes,3,0,n,strongbool::yes);
		for (auto it=rng.begin(); it!=rng.end(); ++it) {
			spliner->SetStrip(it.global_count(),&loRes(*it,0,0),loRes.Stride(3),loRes.Stride(4));
		}
	}
	spliner->Solve();
	#pragma omp parallel
	{
		StripRange loRange(loRes,3,0,n,strongbool::yes);
		StripRange hiRange(hiRes,3,0,n,strongbool::yes);
		for (auto loStrip=loRange.begin(),hiStrip=hiRange.begin(); loStrip!=loRange.end() && hiStrip!=hiRange.end(); ++loStrip,++hiStrip) {
			for (auto i=1;i<=hiRes.Dim(3);i++) {
				tw::Float x = 0.5 + tw::Float(i-0.5)/resolution;
				hiRes.Pack(*hiStrip,i, spliner->Interpolate(x,loStrip.global_count()));
			}
		}
	}
	auto r = Rng04(n,n+1,0,2);
	hiRes.Field::DownwardCopy(r,tw::grid::z,1);
	hiRes.Field::UpwardCopy(r,tw::grid::z,1);
	hiRes.Field::ApplyBoundaryCondition(r);
}

void LaserSolver::Update()
{
	logger::TRACE("start laser update");
	if (!debug) {
		Upsample(HRchi,chi,1);
		propagator->Advance(HRa,HRchi);
		Downsample(HRa,a,1);
		Downsample(HRa,a,2);
	}
}

void LaserSolver::Reset()
{
	chi = tw::Complex(0,0);
}

void LaserSolver::ReadCheckpoint(std::ifstream& inFile)
{
	Driver::ReadCheckpoint(inFile);
	a.ReadCheckpoint(inFile);
	HRa.ReadCheckpoint(inFile);
}

void LaserSolver::WriteCheckpoint(std::ofstream& outFile)
{
	Driver::WriteCheckpoint(outFile);
	a.WriteCheckpoint(outFile);
	HRa.WriteCheckpoint(outFile);
}


//////////////////////////////
//                          //
// Quasistatic Laser Solver //
//                          //
//////////////////////////////


QSSolver::QSSolver(const std::string& name,MetricSpace *ms,Task *tsk):LaserSolver(name,ms,tsk)
{
}


//////////////////////////////
//                          //
//     PGC LASER SOLVER     //
//                          //
//////////////////////////////


PGCSolver::PGCSolver(const std::string& name,MetricSpace *ms,Task *tsk):LaserSolver(name,ms,tsk)
{
	F.Initialize(8,*space,task);
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
	if (space->IsStdMovingWindow())
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
	if (space->bc0[1]==par::axisymmetric || space->bc0[1]==par::reflecting)
	{
		F.SetBoundaryConditions(Rng(0,6),tw::grid::x, fld::neumannWall, fld::neumannWall);
		F.SetBoundaryConditions(Rng(0),tw::grid::x,fld::dirichletWall,fld::neumannWall);
		F.SetBoundaryConditions(Rng(3),tw::grid::x,fld::dirichletWall,fld::neumannWall);
	}
	else
	{
		F.SetBoundaryConditions(Rng(0,6),tw::grid::x, fld::neumannWall, fld::neumannWall);
	}
	if (space->bc0[2]==par::axisymmetric || space->bc0[2]==par::reflecting)
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
	auto rng1 = StripRange(HRa,3,0,1,strongbool::yes);
	auto rng2 = StripRange(HRa,3,0,2,strongbool::yes);
	for (auto nxt=rng1.begin(),prv=rng2.begin(); nxt!=rng1.end() && prv!=rng1.end(); ++nxt,++prv)
	// for (auto [nxt,prv] : std::views::zip(
	// 	StripRange(HRa,3,0,1,strongbool::yes),
	// 	StripRange(HRa,3,0,2,strongbool::yes)))
	{
		tw::Float polarizationFactor = polarizationType==circularPolarization ? 1.414 : 1.0;
		tw::Complex older(0,0);
		tw::Complex newer(0,0);
		for (auto pulse : waves)
		{
			tw::vec3 pos = space->Pos(*nxt,0);
			pos.z = space->ToLab(pos.z,-dth);
			older += polarizationFactor*pulse->VectorPotentialEnvelope(space->WindowPos(0)-dth,pos,laserFreq);
			pos = space->Pos(*nxt,0);
			pos.z = space->ToLab(pos.z,dth);
			newer += polarizationFactor*pulse->VectorPotentialEnvelope(space->WindowPos(0)+dth,pos,laserFreq);
		}
		HRa.Shift(Rng(0,2),*nxt,1,(tw::Float*)&newer);
		HRa.Shift(Rng(0,2),*prv,1,(tw::Float*)&older);
	}
	HRa.Field::UpwardCopy(Rng04(1,3,0,2),tw::grid::z,1);
}

void PGCSolver::Update()
{
	logger::TRACE("start PGC update");
	chi.DepositFromNeighbors();
	chi.ApplyFoldingCondition();
	chi.DivideCellVolume(Rng(0,2),*space);
	chi.ApplyBoundaryCondition();
	chi.Smooth(Rng(0,2),*space,smoothing,compensation);
	#pragma omp parallel
	{
		const tw::Float dth = 0.5*dx(0);
		for (auto s : StripRange(*this,3,0,1,strongbool::yes))
		{
			for (tw::Int k=1;k<=dim[3];k++)
				chi.Pack(s,k, space->ValueOnLightGrid<ComplexField,tw::Complex>(chi,s,k,dth));
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
		for (auto nxt : StripRange(*this,3,0,1,strongbool::yes))
		{
			auto prv = tw::strip(nxt,2);
			for (tw::Int k=1;k<=dim[3];k++)
			{
				F(nxt,k,7) = norm(space->ValueOnLabGrid<ComplexField,tw::Complex>(a,nxt,k,dth));
				F(nxt,k,6) = norm(space->ValueOnLabGrid<ComplexField,tw::Complex>(a,prv,k,-dth));
				F(nxt,k,6) = 0.5*(F(nxt,k,6) + F(nxt,k,7));
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

	diagnostic.SwitchVariant(1);
	diagnostic.ReportField("a_real",HRa,1,0,tw::dims::vector_potential,"$\\Re A$");
	diagnostic.ReportField("a_imag",HRa,1,1,tw::dims::vector_potential,"$\\Im A$");
	diagnostic.ReportField("j1_real",HRchi,1,0,tw::dims::current_density,"$\\Re j$");
	diagnostic.ReportField("j1_imag",HRchi,1,1,tw::dims::current_density,"$\\Im j$");
	diagnostic.SwitchVariant(0);

	ComplexField temp;
	temp.Initialize(*space,task);

	const tw::Float dti = dk(0);
	const tw::Float dth = 0.5*dx(0);
	for (auto [d,nxt,prv] : std::views::zip(
		InteriorCellRange(temp,1),
		InteriorCellRange(a,1),
		InteriorCellRange(a,2)))
	{
		const tw::Complex aNow = half*(a(prv)+a(nxt));
		const tw::Complex dtau = dti*(a(nxt)-a(prv));
		const tw::Complex dzeta = half*(a.d1(prv,0,3) + a.d1(nxt,0,3)) + ii*half*(a.d1(prv,1,3) + a.d1(nxt,1,3));
		const tw::Complex eNow = ii*laserFreq*aNow - (dtau-dzeta);
		const tw::Complex bNow = ii*laserFreq*aNow + dzeta;
		temp(d,0) = 0.25*(norm(eNow) + norm(bNow));
		temp(d,1) = imag( conj(aNow)*bNow - aNow*conj(bNow) );
	}
	diagnostic.VolumeIntegral("LaserEnergy",temp,1,0);
	diagnostic.VolumeIntegral("WaveAction",temp,1,1);

	for (auto [d,nxt,prv] : std::views::zip(
		StripRange(temp,3,0,1,strongbool::no),
		StripRange(a,3,0,1,strongbool::no),
		StripRange(a,3,0,2,strongbool::no)))
	{
		for (tw::Int k=1;k<=dim[3];k++)
		{
			tw::Complex dadt = dti*(space->ValueOnLabGrid<ComplexField,tw::Complex>(a,nxt,k,dth) - space->ValueOnLabGrid<ComplexField,tw::Complex>(a,prv,k,-dth));
			tw::Complex anow = tw::Float(0.5)*(space->ValueOnLabGrid<ComplexField,tw::Complex>(a,prv,k,-dth) + space->ValueOnLabGrid<ComplexField,tw::Complex>(a,nxt,k,dth));
			temp.Pack(d, k, ii*laserFreq*anow - dadt);
		}
	}
	diagnostic.ReportField("e_real",temp,1,0,tw::dims::electric_field,"$\\Re E$");
	diagnostic.ReportField("e_imag",temp,1,1,tw::dims::electric_field,"$\\Im E$");

	diagnostic.ReportField("a2",F,1,7,tw::dims::none,"$a^2$");
	diagnostic.ReportField("chi_real",chi,1,0,tw::dims::none,"$\\Re \\chi$");
	diagnostic.ReportField("chi_imag",chi,1,1,tw::dims::none,"$\\Im \\chi$");
}
