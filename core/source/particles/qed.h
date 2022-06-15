struct QED:ComputeTool
{
	QED(const std::string& name,MetricSpace *m,Task *tsk);
	virtual void Initialize();
	virtual tw::Float CalculateRate(tw::Float,tw::Float);
	virtual tw::Float NewQParameter(tw::Float,tw::Float);

	std::string photon_name, electron_name, positron_name;

	std::string tablesPath;
	std::vector<tw::Float> eta_vector, chi_vector, frac_vector;
	std::vector<std::vector<tw::Float>> LT1, LT2, LT3, LT4, LT5;

	tw::Float photon_cutoff_energy = 0.1f; // default, ~51.1 keV
	const tw::Float alpha_f = 0.0072973525693; // fine-structure constant
	const tw::Float lambda_C = 2.0*pi*cgs::hbar/cgs::me/cgs::c; // Compton wavelength (cm)
};

struct PhotonGenerator:QED
{
	PhotonGenerator(const std::string& name,MetricSpace *m,Task *tsk);
	virtual void Initialize();
	virtual tw::Float CalculateRate(tw::Float eta,tw::Float gamma);
	virtual tw::Float NewQParameter(tw::Float eta,tw::Float Py);
};

struct PairCreator:QED
{
	PairCreator(const std::string& name,MetricSpace *m,Task *tsk);
	virtual void Initialize();
	virtual tw::Float CalculateRate(tw::Float chi,tw::Float energy);
	virtual tw::Float NewQParameter(tw::Float chi,tw::Float Pf);
};
