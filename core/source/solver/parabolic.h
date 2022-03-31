enum tw_polarization_type {linearPolarization,circularPolarization,radialPolarization};

struct LaserPropagator:ComputeTool
{
	ComplexField *p0;
	tw::Float w0,dt;
	tw_polarization_type polarization;
	bool movingWindow;

	LaserPropagator(const std::string& name,MetricSpace *m,Task *tsk);
	virtual ~LaserPropagator();
	virtual void SetData(tw::Float w0,tw::Float dt,tw_polarization_type pol,bool mov);
	void SetBoundaryConditions(ComplexField& a0,ComplexField& a1,ComplexField& chi);
	virtual void Advance(ComplexField& a0,ComplexField& a1,ComplexField& chi) = 0;
};

struct EigenmodePropagator:LaserPropagator
{
	// following selects number of radial modes to keep and absorbing layers
	tw::Int modes,layers;
	// causality = 0.0 gives centered differencing
	// causality = -1.0 uses forward (downwind) differencing
	// causality = 1.0 uses backward (upwind) differencing
	// non-integral values are legal, and give a linear combination
	tw::Float causality,dampingTime;
	std::valarray<tw::Float> eigenvalue,hankel,inverseHankel;
	GlobalIntegrator<tw::Complex>* globalIntegrator;

	EigenmodePropagator(const std::string& name,MetricSpace *m,Task *tsk);
	virtual ~EigenmodePropagator();
	virtual void Initialize();
	virtual void SetData(tw::Float w0,tw::Float dt,tw_polarization_type pol,bool mov);
	virtual void Advance(ComplexField& a0,ComplexField& a1,ComplexField& chi);
};

struct ADIPropagator:LaserPropagator
{
	GlobalIntegrator<tw::Complex>* xGlobalIntegrator;
	GlobalIntegrator<tw::Complex>* yGlobalIntegrator;
	std::valarray<tw::Complex> X1,X2,X3;
	std::valarray<tw::Complex> Y1,Y2,Y3;
	bool evenTime;

	ADIPropagator(const std::string& name,MetricSpace *m,Task *tsk);
	virtual ~ADIPropagator();
	virtual void Advance(ComplexField& a0,ComplexField& a1,ComplexField& chi);
};

struct ConservativePropagator:LaserPropagator
{
};

struct SchroedingerPropagator:ComputeTool
{
	std::vector<GlobalIntegrator<tw::Complex>*> globalIntegrator;

	SchroedingerPropagator(const std::string& name,MetricSpace *m,Task *tsk);
	virtual ~SchroedingerPropagator();
	virtual void DepositCurrent(const tw::grid::axis& axis,ComplexField& psi0,ComplexField& psi1,Field& A4,Field& J4,tw::Complex dt);
	virtual void ApplyNumerator(const tw::grid::axis& axis,ComplexField& psi,Field& A4,bool keepA2Term,tw::Complex dt);
	virtual void ApplyDenominator(const tw::grid::axis& axis,ComplexField& psi,Field& A4,bool keepA2Term,tw::Complex dt);
	virtual void UpdateSpin(ComplexField& psi,ComplexField& chi,Field& A4,tw::Float adt);
};

struct ParabolicSolver:BoundedTool
{
	// inhomogeneous boundary conditions are handled as described in elliptic.h
	std::vector<GlobalIntegrator<tw::Float>*> globalIntegrator;

	ParabolicSolver(const std::string& name,MetricSpace *m,Task *tsk);
	virtual ~ParabolicSolver();

	void FixTemperature(Field& T,const Element& e,Region* theRegion,const tw::Float& T0);
	void FormOperatorStencil(tw::Float *D1,tw::Float *D2,const ScalarField& fluxMask,Field *coeff,tw::Int c,const tw::strip& s,tw::Int i);
	virtual void Advance(const tw::grid::axis& axis,ScalarField& psi,ScalarField& fluxMask,tw::Float coeff,tw::Float dt);
	virtual void Advance(ScalarField& psi,ScalarField& fluxMask,tw::Float coeff,tw::Float dt);
	virtual void Advance(	const tw::grid::axis& axis,
							Field& psi,
							tw::Int psi_idx,
							ScalarField& fluxMask,
							Field *coeff1,
							tw::Int c1_idx,
							Field *coeff2,
							tw::Int c2_idx,
							tw::Float dt);
	virtual void Advance(	Field& psi,
							tw::Int psi_idx,
							ScalarField& fluxMask,
							Field *coeff1,
							tw::Int c1_idx,
							Field *coeff2,
							tw::Int c2_idx,
							tw::Float dt);
};

struct IsotropicPropagator:ComputeTool
{
	GlobalIntegrator<tw::Complex>* zGlobalIntegrator;
	std::valarray<tw::Complex> Z1,Z2,Z3;

	IsotropicPropagator(const std::string& name,MetricSpace *m,Task *tsk);
	virtual ~IsotropicPropagator();

	void SetupIncomingWaveLeft(const tw::strip& s,ComplexField& amplitude,tw::Complex a0,tw::Complex a1,tw::Complex w0) const
	{
		amplitude(s,0) = a1 - a0 + ii*w0*space->dl(s,1,3)*half*(a0+a1);
	}

	void SetupIncomingWaveRight(const tw::strip& s,ComplexField& amplitude,tw::Complex aN,tw::Complex aN1,tw::Complex w0) const
	{
		amplitude(s,space->Dim(s.Axis())+1) = aN1 - aN - ii*w0*space->dl(s,space->Dim(s.Axis())+1,3)*half*(aN+aN1);
	}

	virtual void Advance(ComplexField& amplitude, // this is the array to update
		ComplexField& refractiveIndex, // index of refraction array
		ScalarField& nu_e, // electron collision frequency array
		tw::Float laserFrequency,tw::Float dt);
};
