module;

#include "tw_includes.h"

export module qed;
import input;
import compute_tool;
import fields;

export struct QED:ComputeTool
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

export struct PhotonGenerator:QED
{
	PhotonGenerator(const std::string& name,MetricSpace *m,Task *tsk);
	virtual void Initialize();
	virtual tw::Float CalculateRate(tw::Float eta,tw::Float gamma);
	virtual tw::Float NewQParameter(tw::Float eta,tw::Float Py);
};

export struct PairCreator:QED
{
	PairCreator(const std::string& name,MetricSpace *m,Task *tsk);
	virtual void Initialize();
	virtual tw::Float CalculateRate(tw::Float chi,tw::Float energy);
	virtual tw::Float NewQParameter(tw::Float chi,tw::Float Pf);
};

QED::QED(const std::string& name,MetricSpace *m,Task *tsk) : ComputeTool(name,m,tsk)
{
	ComputeTool::Initialize();
	directives.Add("tables path",new tw::input::String(&tablesPath),false);
}

PhotonGenerator::PhotonGenerator(const std::string& name,MetricSpace *m,Task *tsk) : QED(name,m,tsk)
{
	directives.Add("photon species",new tw::input::String(&photon_name),true);
	directives.Add("cutoff energy",new tw::input::Float(&photon_cutoff_energy),false);
}

PairCreator::PairCreator(const std::string& name,MetricSpace *m,Task *tsk) : QED(name,m,tsk)
{
	directives.Add("electron species",new tw::input::String(&electron_name),true);
	directives.Add("positron species",new tw::input::String(&positron_name),true);
}

void QED::Initialize() {}

void PhotonGenerator::Initialize()
{
	std::ifstream table;
	std::string tableName, line;

	std::vector<std::vector<tw::Float>> temp_matrix;
	std::vector<tw::Float> temp_row;
	tw::Float num;

	if (!tablesPath.empty() && tablesPath.back()!='/')
		tablesPath.append("/");

	/* Key
		Quantum parameters
		  eta: fermions, chi: photons

		Functions
		  Py(eta,chi): Cumulative probability of emitting a chi-photon by an eta-fermion.
		  h(eta): Synchrotron function used to calculate the photo-emission rate.

		Look-up tables
			LT1: element (i,j): chi[j] corresponding to eta[i]
			LT2: row 1: eta[i], row 2: h(eta[i])
			LT3: element (i,j): Py(eta[i],chi[j])

		Notes
		- Some of the raw table values are log10; we invert them upon loading here.
		- The first row of the raw LT3 file contains log10(eta[i]) values, which we store separately.
	*/

	// Load LT1

	tableName = tablesPath + "mmc1.table";
	table.open(tableName);

	if (!table.is_open())
		throw tw::FatalError("couldn't open QED table <" + tableName + ">");

	while (getline(table,line))
	{
		std::istringstream ss(line);
		temp_row.clear();
		while (ss >> num)
			temp_row.push_back(std::pow(10,num));
		LT1.push_back(temp_row);
	} table.close();

	// Load LT2

	tableName = tablesPath + "mmc2.table";
	table.open(tableName);

	if (!table.is_open())
		throw tw::FatalError("couldn't open QED table <" + tableName + ">");

	temp_matrix.clear();

	while (getline(table,line))
	{
		std::istringstream ss(line);
		temp_row.clear();
		while (ss >> num)
			temp_row.push_back(std::pow(10,num));
		temp_matrix.push_back(temp_row);
	} table.close();

	for (int j=0;j<temp_matrix[0].size();j++)
	{
		temp_row.clear();
		for (int i=0;i<temp_matrix.size();i++)
			temp_row.push_back(temp_matrix[i][j]);
		LT2.push_back(temp_row);
	}

	// Load LT3

	tableName = tablesPath + "mmc3.table";
	table.open(tableName);

	if (!table.is_open())
		throw tw::FatalError("couldn't open QED table <" + tableName + ">");

	temp_matrix.clear();

	while (getline(table,line))
	{
		std::istringstream ss(line);
		temp_row.clear();
		while (ss >> num)
			temp_row.push_back(num);
		temp_matrix.push_back(temp_row);
	} table.close();

	for (int j=0;j<temp_matrix[0].size();j++)
		eta_vector.push_back(std::pow(10,temp_matrix[0][j]));

	LT3 = temp_matrix;
	LT3.erase(LT3.begin());
}

void PairCreator::Initialize()
{
	std::ifstream table;
	std::string tableName, line;

	std::vector<std::vector<tw::Float>> temp_matrix;
	std::vector<tw::Float> temp_row;
	tw::Float num;

	if (!tablesPath.empty() && tablesPath.back()!='/')
		tablesPath.append("/");

	/* Key
		Quantum parameters
		  eta: fermions, chi: photons

		Functions
		  Pf(chi,f): Cumulative probability of a chi-photon giving fractions f & 1-f of
		  	its energy to the electron & positron, respectively, of a newly-generated pair.
		  T(chi): Emissivity function used to calculate the pair-creation rate.

		Look-up tables
			LT4: element (i,j): Pf(chi[i],f[j])
			LT5: row 1: chi[i], row2: T(chi[i])

		Notes
		- Some of the raw table values are log10; we invert them upon loading here.
		- The first row of the raw LT4 file contains log10(chi[i]) values, the
		  	second f[j], and the rest is the Pf(chi[i],f[j]) table.
		- The first two rows of the LT4 file are stored separately.
	*/

	// Load LT4

	tableName = tablesPath + "mmc4.table";
	table.open(tableName);

	if (!table.is_open())
		throw tw::FatalError("couldn't open QED table <" + tableName + ">");

	temp_matrix.clear();

	while (getline(table,line))
	{
		std::istringstream ss(line);
		temp_row.clear();
		while (ss >> num)
			temp_row.push_back(num);
		temp_matrix.push_back(temp_row);
	} table.close();

	for (int j=0;j<temp_matrix[0].size();j++)
		chi_vector.push_back(std::pow(10,temp_matrix[0][j]));
	for (int j=0;j<temp_matrix[1].size();j++)
		frac_vector.push_back(temp_matrix[1][j]);

	LT4 = temp_matrix;
	LT4.erase(LT4.begin());
	LT4.erase(LT4.begin());

	// Load LT5

	tableName = tablesPath + "mmc5.table";
	table.open(tableName);

	if (!table.is_open())
		throw tw::FatalError("couldn't open QED table <" + tableName + ">");

	temp_matrix.clear();

	while (getline(table,line))
	{
		std::istringstream ss(line);
		temp_row.clear();
		while (ss >> num)
			temp_row.push_back(std::pow(10,num));
		temp_matrix.push_back(temp_row);
	} table.close();

	for (int j=0;j<temp_matrix[0].size();j++)
	{
		temp_row.clear();
		for (int i=0;i<temp_matrix.size();i++)
			temp_row.push_back(temp_matrix[i][j]);
		LT5.push_back(temp_row);
	}
}

tw::Float QED::CalculateRate(tw::Float a,tw::Float b) {return 0.0f;}

tw::Float PhotonGenerator::CalculateRate(tw::Float eta,tw::Float gamma)
{
	tw::Float rate, synch;

	// Interpolate synchrotron function value
	if (eta <= LT2[0].front())
		synch = LT2[1][0];
	else if (eta >= LT2[0].back())
		synch = LT2[1].back();
	else
	{
		tw::Float x[2], y[2];
		int Idx = NearestIdx(LT2[0],eta);

		x[0] = LT2[0][Idx-1]; x[1] = LT2[0][Idx];
		y[0] = LT2[1][Idx-1]; y[1] = LT2[1][Idx];

		synch = y[0] + (eta-x[0])*(y[1]-y[0])/(x[1]-x[0]);
	}

	rate = std::sqrt(3.0)*(alpha_f/lambda_C)*cgs::c*eta*synch/gamma;
	rate = rate * tw::dims::frequency >> cgs >> native;

	return rate;
}

tw::Float PairCreator::CalculateRate(tw::Float chi,tw::Float energy)
{
	tw::Float rate, emiss;

	// Interpolate pair emissivity value
	if (chi <= LT5[0].front())
		return 0.0f;
	else if (chi >= LT5[0].back())
		emiss = LT5[1].back();
	else
	{
		tw::Float x[2], y[2];
		int Idx = NearestIdx(LT5[0],chi);

		x[0] = LT5[0][Idx-1]; x[1] = LT5[0][Idx];
		y[0] = LT5[1][Idx-1]; y[1] = LT5[1][Idx];

		emiss = y[0] + (chi-x[0])*(y[1]-y[0])/(x[1]-x[0]);
	}

	rate = 2.0*pi*(alpha_f/lambda_C)*cgs::me*cub(cgs::c)*chi*emiss/energy;
	rate = rate * tw::dims::frequency >> cgs >> native;

	return rate;
}

tw::Float QED::NewQParameter(tw::Float a,tw::Float b) {return 0.0f;}

// Returns a new photon quantum parameter given that of the generating fermion
tw::Float PhotonGenerator::NewQParameter(tw::Float eta,tw::Float Py)
{
	tw::Float chi;

	// Fermion quantum parameter out of table range
	if (eta <= eta_vector.front())
		return 0.0f;

	// Look up cumulative emission probability & quantum parameter values from table
	int i = (eta < eta_vector.back()) ? NearestIdx(eta_vector,eta) : eta_vector.size()-1;
	int j1 = NearestIdx(LT3[i-1],Py);
	int j2 = NearestIdx(LT3[i],Py);

	// RNG gave bad Py value; tell caller to re-generate it
	if (j1==0 || j2==0)
		return -1.0f;

	tw::Float eta1 = eta_vector[i-1];
	tw::Float eta2 = eta_vector[i];

	tw::Float chi11 = LT1[i-1][j1-1];
	tw::Float chi12 = LT1[i-1][j1];
	tw::Float chi21 = LT1[i][j2-1];
	tw::Float chi22 = LT1[i][j2];

	tw::Float P11 = LT3[i-1][j1-1];
	tw::Float P12 = LT3[i-1][j1];
	tw::Float P21 = LT3[i][j2-1];
	tw::Float P22 = LT3[i][j2];

	// Interpolate new photon quantum parameter
	chi = (eta2-eta1)*Py + (eta-eta1)*(chi21*P22-chi22*P21)/(chi22-chi21)
		+ (eta2-eta)*(chi11*P12-chi12*P11)/(chi12-chi11);
	chi /= (eta-eta1)*(P22-P21)/(chi22-chi21) + (eta2-eta)*(P12-P11)/(chi12-chi11);

	return chi;
}

// Returns the energy fraction of the generating photon that a fermion will receive
tw::Float PairCreator::NewQParameter(tw::Float chi,tw::Float Pf)
{
	tw::Float frac;

	// Photon quantum parameter out of table range
	if (chi <= chi_vector.front())
		return 0.5f;

	// Look up cumulative creation probability & energy fraction values from table
	int i = (chi < chi_vector.back()) ? NearestIdx(chi_vector,chi) : chi_vector.size()-1;
	int j1 = NearestIdx(LT4[i-1],Pf);
	int j2 = NearestIdx(LT4[i],Pf);

	// RNG gave bad Pf value; tell caller to re-generate it
	if (j1==0 || j2==0)
		return -1.0f;

	tw::Float chi1 = chi_vector[i-1];
	tw::Float chi2 = chi_vector[i];

	tw::Float f11 = frac_vector[j1-1];
	tw::Float f12 = frac_vector[j1];
	tw::Float f21 = frac_vector[j2-1];
	tw::Float f22 = frac_vector[j2];

	tw::Float P11 = LT4[i-1][j1-1];
	tw::Float P12 = LT4[i-1][j1];
	tw::Float P21 = LT4[i][j2-1];
	tw::Float P22 = LT4[i][j2];

	// Interpolate energy fraction
	frac = (chi2-chi1)*Pf + (chi-chi1)*(f21*P22-f22*P21)/(f22-f21)
		+ (chi2-chi)*(f11*P12-f12*P11)/(f12-f11);
	frac /= (chi-chi1)*(P22-P21)/(f22-f21) + (chi2-chi)*(P12-P11)/(f12-f11);

	return frac;
}
