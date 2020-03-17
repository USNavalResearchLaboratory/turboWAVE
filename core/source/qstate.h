namespace qo
{
	enum waveEquation { schroedinger,pauli,klein_gordon,dirac };
	enum potentialType { softCore,bachelet };
}

struct HamiltonianParameters
{
	qo::waveEquation form;
	tw::Float qorb,morb,qnuc,rnuc;
	tw::Float c1,c2,a1,a2; // Bachelet
	tw::vec3 B0; // Magnetic field
};

class QState : public ComputeTool
{
protected:
	bool cylindricalAtom;
	tw::Complex amplitude;

public:
	QState(const std::string& name,MetricSpace *ms,Task *tsk);
	virtual tw::Complex Amplitude(const HamiltonianParameters& H,const tw::vec3& r,const tw::Float& t,const tw::Int& comp) const;
	virtual bool GoodQuantumNumbers(const HamiltonianParameters& H) const;
	virtual tw::Float Energy(const HamiltonianParameters& H) const;
	virtual tw::Float NormalizationConstant(const HamiltonianParameters& H) const;
	virtual void ReadCheckpoint(std::ifstream& inFile);
	virtual void WriteCheckpoint(std::ofstream& outFile);
};

class RandomState : public QState
{
	tw::vec3 size;
public:
	RandomState(const std::string& name,MetricSpace *ms,Task *tsk);
	virtual tw::Complex Amplitude(const HamiltonianParameters& H,const tw::vec3& r,const tw::Float& t,const tw::Int& comp) const;
	virtual void ReadCheckpoint(std::ifstream& inFile);
	virtual void WriteCheckpoint(std::ofstream& outFile);
};

class FreeState : public QState
{
	tw::vec3 size;
	tw::vec4 k4;
	tw::vec3 spin;
public:
	FreeState(const std::string& name,MetricSpace *ms,Task *tsk);
	virtual tw::Complex Amplitude(const HamiltonianParameters& H,const tw::vec3& r,const tw::Float& t,const tw::Int& comp) const;
	virtual void ReadCheckpoint(std::ifstream& inFile);
	virtual void WriteCheckpoint(std::ofstream& outFile);
};

class BoundState : public QState
{
	tw::vec4 qnumbers;
	tw::Float& nr = qnumbers[0]; // radial quantum number
	tw::Float& Jam = qnumbers[1]; // total angular momentum
	tw::Float& Lam = qnumbers[2]; // parity(orbital)
	tw::Float& jzam = qnumbers[3]; // total angular momentum component
	// N.b. quantum numbers have different meanings depending on type of state.
	// Scalars: Jam = [0,1,...], Lam = Jam, jzam = [-Lam,...,Lam]
	// Spinors: Jam = [0.5,1.5,...], Lam = Jam +- 1/2, jzam = [-Jam,...,Jam]
	// Cylindrical Scalars: Lam = jzam , is signed, and Jam is ignored.
	// Cylindrical Spinors: Lam = jzam+-1/2, is signed, and Jam is ignored.
public:
	BoundState(const std::string& name,MetricSpace *ms,Task *tsk);
	virtual tw::Complex Amplitude(const HamiltonianParameters& H,const tw::vec3& r,const tw::Float& t,const tw::Int& comp) const;
	virtual bool GoodQuantumNumbers(const HamiltonianParameters& H) const;
	virtual tw::Float Energy(const HamiltonianParameters& H) const;
	virtual tw::Float NormalizationConstant(const HamiltonianParameters& H) const;
	virtual void ReadCheckpoint(std::ifstream& inFile);
	virtual void WriteCheckpoint(std::ofstream& outFile);
};

class TabulatedState : public QState
{
	tw::Float energy,nr,Lam,Jam,jzam;
	tw::Int components;
	std::string filename;
	std::valarray<tw::Complex> radialFunction;
	tw::vec3 cell_size; // cell size used in lookup table
public:
	TabulatedState(const std::string& name,MetricSpace *ms,Task *tsk);
	virtual void Initialize();
	virtual tw::Complex Amplitude(const HamiltonianParameters& H,const tw::vec3& r,const tw::Float& t,const tw::Int& comp) const;
	virtual void ReadCheckpoint(std::ifstream& inFile);
	virtual void WriteCheckpoint(std::ofstream& outFile);
};
