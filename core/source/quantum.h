namespace qo
{
	enum waveEquationType { schroedinger,pauli,klein_gordon,dirac };
	enum waveFunctionType { random,free,helicity,bound,lookup };
	enum potentialType { softCore,bachelet };

	struct State
	{
		qo::waveEquationType waveEquation;
		qo::waveFunctionType waveType;
		bool cylindricalAtom;
		tw::Complex amplitude;
		tw::Float nr; // radial quantum number
		tw::Float Lam,Jam,jzam; // angular momentum: parity(orbital), total, total-z
		// N.b. quantum numbers have different meanings depending on type of state.
		// Scalars: Jam = [0,1,...], Lam = Jam, jzam = [-Lam,...,Lam]
		// Spinors: Jam = [0.5,1.5,...], Lam = Jam +- 1/2, jzam = [-Jam,...,Jam]
		// Cylindrical Scalars: Lam = jzam , is signed, and Jam is ignored.
		// Cylindrical Spinors: Lam = jzam+-1/2, is signed, and Jam is ignored.
		tw::vec3 k,size; // free or random state parameters
		bool negativeEnergy;
		tw::Float helicity;
		tw::vec3 cell_size; // cell size used in lookup table
		tw::Float qnuc,qorb,morb,energy,normalizationConstant;
		MetricSpace *space;
		UniformDeviate *ud;
		std::valarray<tw::Complex> radialFunction;
		tw::Int components;

		State(MetricSpace *m,UniformDeviate *u);
		void Initialize(qo::waveEquationType waveEq,tw::Float nuclearCharge,tw::Float orbitingCharge,tw::Float orbitingMass);
		bool GoodQuantumNumbers();

		tw::Complex Random(const tw::vec3& r,const tw::Float& t,const tw::Int& comp) const;
		tw::Complex Free(bool helicityState,const tw::vec3& r,const tw::Float& t,const tw::Int& comp) const;
		tw::Complex Bound(const tw::vec3& r,const tw::Float& t,const tw::Int& comp) const;
		tw::Complex Lookup(const tw::vec3& r,const std::valarray<tw::Complex>& fr,const tw::Float& dr,const tw::Float& t,const tw::Int& comp) const;

		tw::Complex Amplitude(const tw::vec3& r,const tw::Float& t,const tw::Int& comp) const;

		void ReadInputFileBlock(std::stringstream& inputString);
		void ReadCheckpoint(std::ifstream& inFile);
		void WriteCheckpoint(std::ofstream& outFile);
	};
}

struct AtomicPhysics:Module
{
	tw::Float alpha; // fine structure constant

	// Specification of potential
	qo::potentialType potentialTypeSpec;
	tw::Float residualCharge;
	tw::Float nuclearRadius; // soft core
	tw::Float c1,c2,a1,a2,Asr[6],asr[3]; // Bachelet
	tw::vec3 B0; // constant magnetic field

	// Electronic state
	std::vector<qo::State> waveFunction;
	std::vector<qo::State> refState;
	tw::Float q0,m0;
	tw::Float groundStateEnergy;

	// Solution method
	LorentzPropagator *photonPropagator;
	bool keepA2Term,dipoleApproximation;
	tw::Float timeRelaxingToGround;

	// Data
	Field psi_r,psi_i; // real and imaginary parts of wavefunction
	Field A4,Ao4; // EM 4-potential and old 4-potential
	Field J4; // EM 4-current

	AtomicPhysics(const std::string& name,Simulation* sim);
	virtual ~AtomicPhysics();
	void ExchangeResources();
	tw::Float GetSphericalPotential(tw::Float r) const;
	virtual void Initialize();
	void FormPotentials(tw::Float t);
	void FormGhostCellPotentials(tw::Float t);
	tw::vec4 GetA4AtOrigin();

	tw::Complex Psi(const tw::Int& i,const tw::Int& j,const tw::Int& k,const tw::Int& c) const
	{
		return tw::Complex(psi_r(i,j,k,c),psi_i(i,j,k,c));
	}
	tw::Complex Psi(const tw::cell& cell,const tw::Int& c) const
	{
		return tw::Complex(psi_r(cell,c),psi_i(cell,c));
	}
	tw::Complex Psi(const tw::xstrip<1>& v,const tw::Int& i,const tw::Int& c) const
	{
		return tw::Complex(psi_r(v,i,c),psi_i(v,i,c));
	}
	virtual void VerifyInput();
	virtual void ReadInputFileDirective(std::stringstream& inputString,const std::string& command);
	virtual void ReadCheckpoint(std::ifstream& inFile);
	virtual void WriteCheckpoint(std::ofstream& outFile);
};

struct Schroedinger:AtomicPhysics
{
	ComplexField psi0,psi1; // old and new wavefunction
	ComplexField scratch,v,w; // for distributed parallelism
	#ifdef USE_OPENCL
	cl_kernel applyNumerator,applyDenominator,chargeKernel,currentKernel,stitchKernel;
	#endif
	SchroedingerPropagator *propagator;

	Schroedinger(const std::string& name,Simulation* sim);
	virtual ~Schroedinger();
	virtual void Initialize();
	virtual void Update();
	virtual void VerifyInput();

	virtual void UpdateJ4();
	virtual void Normalize();

	virtual void StartDiagnostics();
	virtual void Report(Diagnostic&);
};

struct Pauli:AtomicPhysics
{
	ComplexField psi0,psi1,chi0,chi1; // old and new wavefunction
	ComplexField scratch,v,w; // for distributed parallelism
	SchroedingerPropagator *propagator;

	Pauli(const std::string& name,Simulation* sim);
	virtual ~Pauli();
	virtual void Initialize();
	virtual void Update();
	virtual void VerifyInput();

	//virtual void UpdateJ4();
	virtual void Normalize();

	virtual void StartDiagnostics();
	virtual void Report(Diagnostic&);
};

struct KleinGordon:AtomicPhysics
{
	#ifdef USE_OPENCL
	cl_kernel updatePsi;
	cl_kernel updateChi;
	#endif

	KleinGordon(const std::string& name,Simulation* sim);
	virtual ~KleinGordon();
	virtual void Initialize();
	virtual void Update();

	tw::Float ComputeRho(const tw::cell& cell);
	// The following give the Feshbach-Villars decomposition
	// If sgn = 1.0, returns electron part, if sgn=-1.0 returns positron part
	tw::Complex FV(const tw::Int& i,const tw::Int& j,const tw::Int& k,const tw::Float& sgn)
	{
		return (Psi(i,j,k,0) + sgn*Psi(i,j,k,1))/root2;
	}
	tw::Complex FV(const tw::cell& cell,const tw::Float sgn)
	{
		return (Psi(cell,0) + sgn*Psi(cell,1))/root2;
	}
	tw::Complex FV(const tw::xstrip<1>& v,const tw::Int i,const tw::Float& sgn)
	{
		return (Psi(v,i,0) + sgn*Psi(v,i,1))/root2;
	}
	virtual void UpdateJ4();
	virtual void Normalize();

	virtual void StartDiagnostics();
	virtual void Report(Diagnostic&);
};

struct Dirac:AtomicPhysics
{
	#ifdef USE_OPENCL
	cl_kernel leapFrog;
	#endif

	Dirac(const std::string& name,Simulation* sim);
	virtual ~Dirac();
	virtual void Initialize();
	template <tw::Int OUT1,tw::Int OUT2,tw::Int IN1,tw::Int IN2>
	void LeapFrog(tw::Float sgn);
	virtual void Update();

	tw::Float ComputeRho(const tw::cell& cell);
	virtual void UpdateJ4();
	virtual void Normalize();

	virtual void StartDiagnostics();
	virtual void Report(Diagnostic&);
};

template <tw::Int OUT1,tw::Int OUT2,tw::Int IN1,tw::Int IN2>
void Dirac::LeapFrog(tw::Float sgn)
{
	// Leapfrog electron/positron component over positron/electron component
	// The spinor (OUT1,OUT2) is leapfrogged over the spinor (IN1,IN2)
	// If (OUT1,OUT2) is the electron then sgn=1.0 (positive energy)
	// If (OUT1,OUT2) is the positron then sgn=-1.0 (negative energy)

	static const tw::Int AB = tw::vec_align_bytes;
	#pragma omp parallel firstprivate(sgn)
	{
		alignas(AB) tw::Float Ur[dim[1]],Ui[dim[1]];
		alignas(AB) tw::Float D1r[dim[1]],D1i[dim[1]];
		alignas(AB) tw::Float D2r[dim[1]],D2i[dim[1]];

		for (auto v : VectorStripRange<1>(*this,false))
		{
			// unitary operator of time translation for diagonal part of Hamiltonian
			#pragma omp simd aligned(Ur,Ui:AB)
			for (tw::Int i=1;i<=dim[1];i++)
			{
				const tw::Float dq = dt*(sgn*m0+q0*A4(v,i,0));
				Ur[i-1] = cos(dq);
				Ui[i-1] = -sin(dq);
			}

			#pragma omp simd aligned(D1r,D2r,D1i,D2i:AB)
			for (tw::Int i=1;i<=dim[1];i++)
			{
				D1r[i-1] = -psi_r(v,i,IN1,3) - q0*A4(v,i,3)*psi_i(v,i,IN1);
				D1i[i-1] = -psi_i(v,i,IN1,3) + q0*A4(v,i,3)*psi_r(v,i,IN1);
				D2r[i-1] = -psi_r(v,i,IN1,1) + psi_i(v,i,IN1,2) - q0*A4(v,i,1)*psi_i(v,i,IN1) - q0*A4(v,i,2)*psi_r(v,i,IN1);
				D2i[i-1] = -psi_i(v,i,IN1,1) - psi_r(v,i,IN1,2) + q0*A4(v,i,1)*psi_r(v,i,IN1) - q0*A4(v,i,2)*psi_i(v,i,IN1);
			}
			#pragma omp simd aligned(D1r,D2r,D1i,D2i:AB)
			for (tw::Int i=1;i<=dim[1];i++)
			{
				D1r[i-1] += -psi_r(v,i,IN2,1) - psi_i(v,i,IN2,2) - q0*A4(v,i,1)*psi_i(v,i,IN2) + q0*A4(v,i,2)*psi_r(v,i,IN2);
				D1i[i-1] += -psi_i(v,i,IN2,1) + psi_r(v,i,IN2,2) + q0*A4(v,i,1)*psi_r(v,i,IN2) + q0*A4(v,i,2)*psi_i(v,i,IN2);
				D2r[i-1] += psi_r(v,i,IN2,3) + q0*A4(v,i,3)*psi_i(v,i,IN2);
				D2i[i-1] += psi_i(v,i,IN2,3) - q0*A4(v,i,3)*psi_r(v,i,IN2);
			}

			#pragma omp simd aligned(Ur,Ui,D1r,D1i:AB)
			for (tw::Int i=1;i<=dim[1];i++)
			{
				psi_r(v,i,OUT1) += dth*D1r[i-1];
				psi_i(v,i,OUT1) += dth*D1i[i-1];
				complex_multiply_assign(psi_r(v,i,OUT1),psi_i(v,i,OUT1),Ur[i-1],Ui[i-1]);
				psi_r(v,i,OUT1) += dth*D1r[i-1];
				psi_i(v,i,OUT1) += dth*D1i[i-1];
			}

			#pragma omp simd aligned(Ur,Ui,D2r,D2i:AB)
			for (tw::Int i=1;i<=dim[1];i++)
			{
				psi_r(v,i,OUT2) += dth*D2r[i-1];
				psi_i(v,i,OUT2) += dth*D2i[i-1];
				complex_multiply_assign(psi_r(v,i,OUT2),psi_i(v,i,OUT2),Ur[i-1],Ui[i-1]);
				psi_r(v,i,OUT2) += dth*D2r[i-1];
				psi_i(v,i,OUT2) += dth*D2i[i-1];
			}
		}
	}
}
