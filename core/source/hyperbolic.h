struct YeePropagatorPML:ComputeTool
{
	YeePropagatorPML(const std::string& name,MetricSpace *m,Task *tsk);

	#ifdef USE_OPENCL
	virtual ~YeePropagatorPML();
	cl_kernel k_advanceE,k_prepCenteredFields,k_advanceB,k_centeredFields;
	void SetupComputeKernels(Field& F,Field& A,Field& PMLx,Field& PMLy,Field& PMLz,Field& j4);
	#endif

	void AdvanceE(Field& A,Field& PMLx,Field& PMLy,Field& PMLz,Field& j4);
	void AdvanceB(Field& A,Field& PMLx,Field& PMLy,Field& PMLz);
	void PrepCenteredFields(Field& F,Field& A);
	void CenteredFields(Field& F,Field& A);
	void UpdateInteriorBoundaryE(Field& A,const ScalarField& conductor);
	void UpdateInteriorBoundaryB(Field& A,const ScalarField& conductor);
	void UpdateExteriorBoundary(Field& A,Field& PMLx,Field& PMLy,Field& PMLz);
};

struct LorentzPropagator:ComputeTool
{
	LorentzPropagator(const std::string& name,MetricSpace *m,Task *tsk);

	#ifdef USE_OPENCL
	virtual ~LorentzPropagator();
	cl_kernel k_advance,k_swap,k_midstep,k_undoMidstep;
	void SetupComputeKernels(Field& A4,Field& Ao4,Field& j4);
	#endif

	void Advance(Field& A4,Field& Ao4,Field& j4,const tw::Float mult,const tw::Float dt);
	void MidstepEstimate(Field& A4,Field& Ao4);
	void UndoMidstepEstimate(Field& A4,Field& Ao4);
};
