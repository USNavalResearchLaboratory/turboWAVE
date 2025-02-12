struct BoundElectrons:Module
{
	tw::Float q0,m0; // charge, rest mass
	ScalarField dens;
	Field R0,R1;
	tw::vec3 resFreq,dampFreq,oscStrength;
	// following are rows of contracted anharmonic tensor, indexed from 1
	std::valarray<tw::Float> a1,a2,a3; // second order anharmonic coefficients
	std::valarray<tw::Float > packet; // packet to send to compute device
	tw::Float b,d; // third/fifth order anharmonic coefficient (isotropic)

	#ifdef USE_OPENCL
	cl_kernel k_update;
	cl_mem packet_buffer;
	#endif

	tw::Float theta,phi;
	tw::basis crystalBasis;

	Field* EM;
	Field* sources;
	Field* laser;
	ComplexField* chi;
	tw::Float* carrierFrequency;
	bool* circular;
	ScalarField* fixed;
	Vec3Field* ESField;

	BoundElectrons(const std::string& name,Simulation* sim);
	virtual ~BoundElectrons();
	virtual void Initialize();
	virtual bool InspectResource(void* resource,const std::string& description);
	virtual void MoveWindow();
	virtual void Update();

	virtual bool ReadInputFileDirective(const TSTreeCursor *curs,const std::string& src);
	virtual void ReadCheckpoint(std::ifstream& inFile);
	virtual void WriteCheckpoint(std::ofstream& outFile);

	virtual void StartDiagnostics();
	virtual void Report(Diagnostic&);
};
