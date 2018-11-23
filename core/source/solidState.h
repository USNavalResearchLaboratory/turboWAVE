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

	BoundElectrons(const std::string& name,Grid* theGrid);
	~BoundElectrons();
	virtual void Initialize();
	virtual bool InspectResource(void* resource,const std::string& description);
	virtual void MoveWindow();
	virtual void Update();

	virtual void ReadInputFileDirective(std::stringstream& inputString,const std::string& command);
	virtual void ReadData(std::ifstream& inFile);
	virtual void WriteData(std::ofstream& outFile);

	virtual void StartDiagnostics();
	virtual void EnergyHeadings(std::ofstream& outFile);
	virtual void EnergyColumns(std::vector<tw::Float>& cols,std::vector<bool>& avg,const Region& theRgn);
	virtual void BoxDiagnosticHeader(GridDataDescriptor*);
	virtual void BoxDiagnose(GridDataDescriptor*);
	virtual void PointDiagnosticHeader(std::ofstream& outFile);
	virtual void PointDiagnose(std::ofstream& outFile,const weights_3D& w);
};
