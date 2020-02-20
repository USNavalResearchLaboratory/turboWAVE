void WriteDVHeader(std::ofstream& outFile,tw::Int version,tw::Int xDim,tw::Int yDim,tw::Int zDim,float x0,float x1,float y0,float y1,float z0,float z1);

struct DiagnosticDescriptor
{
	Region *theRgn;
	std::vector<Region*>& rgnList;
	std::string filename;
	tw::Int period;
	tw::Float tRef,t0,t1,timePeriod;
	bool moveWithWindow;

	DiagnosticDescriptor(std::vector<Region*>& rgnList);
	bool WriteThisStep(tw::Float elapsedTime,tw::Float dt,tw::Int stepNow);
	void ReadData(std::ifstream& inFile);
	void WriteData(std::ofstream& outFile);
};

struct EnergySeriesDescriptor : DiagnosticDescriptor
{
	tw::Int numSigFigs;

	EnergySeriesDescriptor(std::vector<Region*>& rgnList);
	void ReadInputFile(std::stringstream& inputString);
	void ReadData(std::ifstream& inFile);
	void WriteData(std::ofstream& outFile);
};

struct PointSeriesDescriptor : DiagnosticDescriptor
{
	tw::vec3 thePoint;

	PointSeriesDescriptor(std::vector<Region*>& rgnList);
	void ReadInputFile(std::stringstream& inputString);
	void ReadData(std::ifstream& inFile);
	void WriteData(std::ofstream& outFile);
};

struct GridDataDescriptor : DiagnosticDescriptor
{
	bool average;
	tw::Int xSkip,ySkip,zSkip;

	GridDataDescriptor(std::vector<Region*>& rgnList);
	void ReadInputFile(std::stringstream& inputString);
	void ReadData(std::ifstream& inFile);
	void WriteData(std::ofstream& outFile);
};

struct PhaseSpaceDescriptor : DiagnosticDescriptor
{
	tw::vec3 min,max;
	std::string hAxis,vAxis;
	tw::Int hDim,vDim;
	std::valarray<tw::Float> horVolFactor,verVolFactor;
	std::ofstream outFile;

	PhaseSpaceDescriptor(std::vector<Region*>& rgnList);
	void SetupGeometry(tw::dom::geometry theGeometry);
	void ReadInputFile(std::stringstream& inputString);
	void ReadData(std::ifstream& inFile);
	void WriteData(std::ofstream& outFile);
};

struct ParticleDetectorDescriptor : DiagnosticDescriptor
{
	tw::basis orientation;
	tw::vec3 position;
	tw::Float minGamma;
	std::ofstream outFile;

	ParticleDetectorDescriptor(std::vector<Region*>& rgnList);
	void ReadInputFile(std::stringstream& inputString);
	void ReadData(std::ifstream& inFile);
	void WriteData(std::ofstream& outFile);
	void WriteRecords(std::valarray<float>& records);
};

struct ParticleOrbitDescriptor : DiagnosticDescriptor
{
	tw::Float minGamma;
	std::ofstream outFile;

	ParticleOrbitDescriptor(std::vector<Region*>& rgnList);
	void ReadInputFile(std::stringstream& inputString);
	void ReadData(std::ifstream& inFile);
	void WriteData(std::ofstream& outFile);
};

struct FarFieldDetectorDescriptor : DiagnosticDescriptor
{
	tw::Float radius,theta0,theta1,phi0,phi1;
	tw::Int thetaPts,phiPts,timePts;
	std::ofstream AthetaFile,AphiFile;
	MetricSpace *space;
	Task *task;
	Vec3Field A; // index as (t,theta,phi)

	FarFieldDetectorDescriptor(MetricSpace *space,Task *task,std::vector<Region*>& rgnList);
	void ReadInputFile(std::stringstream& inputString);
	void ReadData(std::ifstream& inFile);
	void WriteData(std::ofstream& outFile);
	void AccumulateField(const tw::Float& elapsedTime,Field& J4);
};
