module;

#include "tw_includes.h"

export module qstate;
import input;
import driver;
import metric_space;
import functions;

export namespace qo
{
	enum waveEquation { schroedinger,pauli,klein_gordon,dirac };
	enum potentialType { softCore,bachelet };
}

export struct HamiltonianParameters
{
	qo::waveEquation form;
	tw::Float qorb,morb,qnuc,rnuc;
	tw::Float c1,c2,a1,a2; // Bachelet
	tw::vec3 B0; // Magnetic field
};

export class QState : public ComputeTool
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
};

export class RandomState : public QState
{
	tw::vec3 size;
public:
	RandomState(const std::string& name,MetricSpace *ms,Task *tsk);
	virtual tw::Complex Amplitude(const HamiltonianParameters& H,const tw::vec3& r,const tw::Float& t,const tw::Int& comp) const;
};

export class FreeState : public QState
{
	tw::vec3 size;
	tw::vec4 k4;
	tw::vec3 spin;
public:
	FreeState(const std::string& name,MetricSpace *ms,Task *tsk);
	virtual tw::Complex Amplitude(const HamiltonianParameters& H,const tw::vec3& r,const tw::Float& t,const tw::Int& comp) const;
};

export class BoundState : public QState
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
};

export class TabulatedState : public QState
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
};

QState::QState(const std::string& name,MetricSpace *ms,Task *tsk) : ComputeTool(name,ms,tsk)
{
	amplitude = tw::Complex(1.0,0.0);
	if (ms->Dim(3)==1)
		cylindricalAtom = true;
	else
		cylindricalAtom = false;
	directives.Add("amplitude",new tw::input::Complex(&amplitude),false);
	directives.Add("cylindrical",new tw::input::Bool(&cylindricalAtom),false);
}

bool QState::GoodQuantumNumbers(const HamiltonianParameters& H) const
{
	return true;
}

tw::Float QState::Energy(const HamiltonianParameters& H) const
{
	return 0.0;
}

tw::Float QState::NormalizationConstant(const HamiltonianParameters& H) const
{
	return 1.0;
}

tw::Complex QState::Amplitude(const HamiltonianParameters& H,const tw::vec3& r,const tw::Float& t,const tw::Int& comp) const
{
	return tw::Complex(0.0,0.0);
}

RandomState::RandomState(const std::string& name,MetricSpace *ms,Task *tsk) : QState(name,ms,tsk)
{
	directives.Add("size",new tw::input::Vec3(&size));
}

tw::Complex RandomState::Amplitude(const HamiltonianParameters& H,const tw::vec3& r,const tw::Float& t,const tw::Int& comp) const
{
	tw::Complex ans;
	ans = task->uniformDeviate->Next() * std::exp( -sqr(r.x/size.x) - sqr(r.y/size.y) - sqr(r.z/size.z) );
	ans *= amplitude * std::exp( ii*task->uniformDeviate->Next()*two*pi );
	return ans;
}


FreeState::FreeState(const std::string& name,MetricSpace *ms,Task *tsk) : QState(name,ms,tsk)
{
	directives.Add("size",new tw::input::Vec3(&size));
	directives.Add("k4",new tw::input::Vec4(&k4));
	directives.Add("spin",new tw::input::Vec3(&spin));
}

tw::Complex FreeState::Amplitude(const HamiltonianParameters& H,const tw::vec3& r,const tw::Float& t,const tw::Int& comp) const
{
	const tw::spinor w(spin/Magnitude(spin),0.5);
	const tw::vec3 k = k4.spatial();
	const bool relativistic = H.form==qo::klein_gordon || H.form==qo::dirac;
	const tw::Float energy = (relativistic ? std::sqrt((k^k)+H.morb*H.morb) : 0.5*(k^k)/H.morb);
	tw::Complex ans(amplitude);
	ans *= std::exp( -sqr(r.x/size.x) - sqr(r.y/size.y) - sqr(r.z/size.z) );
	ans *= std::exp( ii*((k^r) - energy*t) );
	if (H.form==qo::pauli)
		ans *= w[comp];
	if (H.form==qo::dirac)
	{
		const tw::Float signedEnergy = k4[0]<0.0 ? -energy : energy;
		// Create an electron bispinor coefficient in the rest frame, L&L QED 23
		tw::bispinor u(w);
		// Boost to desired momentum
		u.Boost(tw::vec4(signedEnergy,k));
		// If negative energy apply charge conjugation operator
		// This reverses momentum and spin
		if (signedEnergy<0.0)
		{
			ans = std::conj(ans);
			u.ChargeConjugate();
		}
		ans *= u[comp];
	}
	return ans;
}


BoundState::BoundState(const std::string& name,MetricSpace *ms,Task *tsk) : QState(name,ms,tsk)
{
	directives.Add("nr_j_l_m",new tw::input::Vec4(&qnumbers));
}

tw::Float BoundState::Energy(const HamiltonianParameters& H) const
{
	tw::Float energy;
	const tw::Float qnuc = H.qnuc;
	const tw::Float qorb = H.qorb;
	const tw::Float morb = H.morb;
	const tw::Float n = nr + Lam + 1.0;
	switch (H.form)
	{
		case qo::schroedinger:
			energy = -0.5*morb*sqr(qnuc*qorb/n);
			break;
		case qo::pauli:
			energy = -0.5*morb*sqr(qnuc*qorb/n);
			break;
		case qo::klein_gordon:
			if (cylindricalAtom)
				energy = morb/std::sqrt(1.0 + sqr(qnuc*qorb)/sqr(nr+0.5+std::sqrt(sqr(jzam)-sqr(qnuc*qorb))));
			else
				energy = morb/std::sqrt(1.0 + sqr(qnuc*qorb)/sqr(nr+0.5+std::sqrt(sqr(Lam+0.5)-sqr(qnuc*qorb))));
			break;
		case qo::dirac:
			if (cylindricalAtom)
				energy = morb/std::sqrt(1.0 + sqr(qnuc*qorb)/sqr(nr+std::sqrt(sqr(jzam)-sqr(qnuc*qorb))));
			else
				energy = morb/std::sqrt(1.0 + sqr(qnuc*qorb)/sqr(nr+std::sqrt(sqr(Jam+0.5)-sqr(qnuc*qorb))));
			break;
	}
	return energy;
}

tw::Float BoundState::NormalizationConstant(const HamiltonianParameters& H) const
{
	tw::Float normalizationConstant;
	const tw::Float n = nr + Lam + 1.0;
	const tw::Float kr = -H.morb*H.qnuc*H.qorb/n;
	switch (H.form)
	{
		case qo::schroedinger:
			normalizationConstant = std::pow(kr,1.5)*std::sqrt(Factorial(n+Lam)/Factorial(nr))/(Factorial(2.0*Lam+1.0)*std::sqrt(n*pi));
			break;
		case qo::pauli:
			normalizationConstant = std::pow(kr,1.5)*std::sqrt(Factorial(n+Lam)/Factorial(nr))/(Factorial(2.0*Lam+1.0)*std::sqrt(n*pi));
			break;
		default:
			normalizationConstant = 1.0;
	}
	return normalizationConstant;
}

bool BoundState::GoodQuantumNumbers(const HamiltonianParameters& H) const
{
	tw::Float int_part;
	auto good_lam = [&] (tw::Float test) { return std::modf(test,&int_part)==0.0; };
	auto good_jam = [&] (tw::Float test) { return std::abs(std::modf(test,&int_part))==0.5; };
	if (H.form==qo::schroedinger || H.form==qo::klein_gordon)
	{
		if (cylindricalAtom)
			return good_lam(jzam) && jzam==Lam;
		return Jam>=0.0 && good_lam(Jam) && good_lam(jzam) && Jam==Lam && jzam>=-Jam && jzam<=Jam;
	}
	if (H.form==qo::pauli || H.form==qo::dirac)
	{
		if (cylindricalAtom)
			return good_jam(jzam) && std::abs(Lam-jzam)==0.5 && ((nr==0.0 && (Lam-jzam)*jzam<0.0) || nr>0.0);
		return Jam>0.0 && good_jam(Jam) && good_jam(jzam) && std::abs(Jam-Lam)==0.5 && jzam>=-Jam && jzam<=Jam && ((nr==0.0 && Lam<Jam) || nr>0.0);
	}
	return true;
}

tw::Complex BoundState::Amplitude(const HamiltonianParameters& H,const tw::vec3& r,const tw::Float& t,const tw::Int& comp) const
{
	tw::Complex ans;
	tw::vec3 R = r;
	const tw::Float energy = Energy(H);
	const tw::Float normalizationConstant = NormalizationConstant(H);
	const tw::Float qnuc = H.qnuc;
	const tw::Float qorb = H.qorb;
	const tw::Float morb = H.morb;

	if (H.form==qo::schroedinger)
	{
		if (cylindricalAtom)
		{
			space->CurvilinearToCylindrical(&R);
			ans = tw::Complex(0,0);
		}
		else
		{
			space->CurvilinearToSpherical(&R);
			tw::Float kr = std::sqrt(-two*morb*energy);
			ans = ConfluentHypergeometric(-tw::Int(nr),two*(Lam+one),two*kr*R.x)*std::pow(two*kr*R.x,Lam)*std::exp(-kr*R.x);
			ans *= SphericalHarmonic(Lam,jzam,R.y,R.z);
		}
	}

	if (H.form==qo::pauli)
	{
		if (cylindricalAtom)
		{
			space->CurvilinearToCylindrical(&R);
			ans = tw::Complex(0,0);
		}
		else
		{
			space->CurvilinearToSpherical(&R);
			tw::Float kr = std::sqrt(-two*morb*energy);
			ans = ConfluentHypergeometric(-tw::Int(nr),two*(Lam+one),two*kr*R.x)*std::pow(two*kr*R.x,Lam)*std::exp(-kr*R.x);
			ans *= SphericalHarmonic(Lam,jzam,R.y,R.z);
			if (comp==0)
				ans *= Jam-Lam+0.5;
			if (comp==1)
				ans *= 0.5-Jam+Lam;
		}
	}

	if (H.form==qo::klein_gordon)
	{
		tw::Float kr = std::sqrt(morb*morb - sqr(energy));
		if (cylindricalAtom)
		{
			space->CurvilinearToCylindrical(&R);
			const tw::Float gam = std::sqrt(sqr(jzam) - sqr(qnuc*qorb));
			ans = ConfluentHypergeometric(-tw::Int(nr),two*gam+one,two*kr*R.x);
			ans *= std::pow(R.x,gam)*std::exp(-kr*R.x);
			ans *= std::exp(ii*jzam*R.y);
			if (comp==1)
				ans *= (energy - qnuc*qorb/R.x)/morb;
		}
		else
		{
			space->CurvilinearToSpherical(&R);
			const tw::Float gam = std::sqrt(sqr(Lam+0.5)-sqr(qnuc*qorb));
			ans = ConfluentHypergeometric(-tw::Int(nr),two*gam+one,two*kr*R.x);
			ans *= std::pow(R.x,gam-0.5)*std::exp(-kr*R.x);
			ans *= SphericalHarmonic(Lam,jzam,R.y,R.z);
			if (comp==1)
				ans *= (energy - qnuc*qorb/R.x)/morb;
		}
	}

	if (H.form==qo::dirac)
	{
		tw::Complex Q1,Q2;
		tw::Float kr = std::sqrt(morb*morb - sqr(energy));
		if (cylindricalAtom)
		{
			space->CurvilinearToCylindrical(&R);
			const tw::Float szam = jzam-Lam;
			const tw::Float kappa = -2*szam*jzam;
			const tw::Float gam = std::sqrt(sqr(kappa) - sqr(qnuc*qorb));
			Q1 = ConfluentHypergeometric(-tw::Int(nr),two*gam+one,two*kr*R.x);
			Q2 = ConfluentHypergeometric(1-tw::Int(nr),two*gam+one,two*kr*R.x);
			Q2 *= -(kappa - qnuc*qorb*morb/kr) / (gam - qnuc*qorb*energy/kr);
			tw::Complex f = tw::Float(std::sqrt(morb+energy)*std::exp(-kr*R.x)*std::pow(R.x,gam-0.5))*(Q1+Q2);
			tw::Complex g = ii*tw::Float(std::sqrt(morb-energy)*std::exp(-kr*R.x)*std::pow(R.x,gam-0.5))*(Q1-Q2);
			if (comp==0)
				ans = f*(szam+0.5)*std::exp(ii*(jzam-half)*R.y);
			if (comp==1)
				ans = f*(0.5-szam)*std::exp(ii*(jzam+half)*R.y);
			if (comp==2)
				ans = g*(0.5-szam)*std::exp(ii*(jzam-half)*R.y);
			if (comp==3)
				ans = g*(szam+0.5)*std::exp(ii*(jzam+half)*R.y);
		}
		else
		{
			space->CurvilinearToSpherical(&R);
			const tw::Float Sam = Jam-Lam;
			const tw::Float Lamp = 2*Jam - Lam;
			const tw::Float kappa = -2*Sam*(Jam+0.5);
			const tw::Float gam = std::sqrt(sqr(kappa)-sqr(qorb*qnuc));
			Q1 = ConfluentHypergeometric(-tw::Int(nr),two*gam+one,two*kr*R.x);
			Q2 = ConfluentHypergeometric(1-tw::Int(nr),two*gam+one,two*kr*R.x);
			Q2 *= -(gam - qnuc*qorb*energy/kr)/(kappa - qnuc*qorb*morb/kr);
			tw::Complex f = tw::Float(std::sqrt(morb+energy)*std::exp(-kr*R.x)*std::pow(2.0*kr*R.x,gam-1.0))*(Q1+Q2);
			tw::Complex g = tw::Float(std::sqrt(morb-energy)*std::exp(-kr*R.x)*std::pow(2.0*kr*R.x,gam-1.0))*(Q1-Q2);
			ans = 0.0;
			if (comp==0)
				ans = f*SphericalHarmonicSpinor(Jam,Lam,jzam,R.y,R.z)[0];
			if (comp==1)
				ans = f*SphericalHarmonicSpinor(Jam,Lam,jzam,R.y,R.z)[1];
			if (comp==2)
				ans = tw::Float(std::pow(-1.0,0.5*(1.0+Lam-Lamp)))*g*SphericalHarmonicSpinor(Jam,Lamp,jzam,R.y,R.z)[0];
			if (comp==3)
				ans = tw::Float(std::pow(-1.0,0.5*(1.0+Lam-Lamp)))*g*SphericalHarmonicSpinor(Jam,Lamp,jzam,R.y,R.z)[1];
		}
	}

	return normalizationConstant*amplitude*std::exp(-ii*energy*t)*ans;
}


TabulatedState::TabulatedState(const std::string& name,MetricSpace *ms,Task *tsk) : QState(name,ms,tsk)
{
	directives.Add("filename",new tw::input::String(&filename));
}

void TabulatedState::Initialize()
{
	QState::Initialize();
	tw::Int num;
	std::string word;
	std::stringstream contents;
	tw::Float realAmplitude,imagAmplitude;
	tw::input::FileEnv file_env(task->inputFileName);
	if (!file_env.FindAndOpen(filename,contents))
		throw tw::FatalError("couldn't open state file <" + filename + ">");
	do
	{
		contents >> word;
		//if (word=="soft_core_radius")
		//if (word=="nuclear_charge")
		//if (word=="Bz")
		if (word=="energy")
			contents >> word >> energy;
		if (word=="pts")
			contents >> word >> num;
		if (word=="components")
			contents >> word >> components;
		if (word=="cell_width")
			contents >> word >> cell_size.x;
		if (word=="nr_j_l_m")
			contents >> word >> nr >> Jam >> Lam >> jzam;
		if (word=="cylindrical")
		{
			contents >> word >> word;
			cylindricalAtom = (word=="true" || word=="yes" || word=="on") ? true : false;
		}
	} while (word!="start_data");

	radialFunction.resize(components*num);
	for (tw::Int i=0;i<components*num;i++)
	{
		contents >> realAmplitude >> imagAmplitude;
		radialFunction[i] = tw::Complex(realAmplitude,imagAmplitude);
	}
}

tw::Complex TabulatedState::Amplitude(const HamiltonianParameters& H,const tw::vec3& r,const tw::Float& t,const tw::Int& comp) const
{
	tw::Complex ans;
	tw::vec3 R = r;
	const tw::Float dr = cell_size.x;
	if (cylindricalAtom)
		space->CurvilinearToCylindrical(&R);
	else
		space->CurvilinearToSpherical(&R);
	tw::Int dim = radialFunction.size()/components;
	tw::Int s = MyFloor((R.x-0.5*dr)/dr);
	tw::Float w = (R.x-0.5*dr)/dr - tw::Float(s);
	if (s>=0 && s<dim-1)
		ans = radialFunction[comp*dim+s]*(one-w) + radialFunction[comp*dim+s+1]*w;
	if (s<0)
		ans = radialFunction[comp*dim];
	if (s==dim-1)
		ans = radialFunction[comp*dim+dim-1];
	if (cylindricalAtom)
	{
		if (H.form==qo::schroedinger || H.form==qo::klein_gordon)
		{
			ans *= std::exp(ii*jzam*R.y);
		}
		if (H.form==qo::pauli || H.form==qo::dirac)
		{
			if (comp==0) ans *= std::exp(ii*(jzam-half)*R.y);
			if (comp==1) ans *= std::exp(ii*(jzam+half)*R.y);
			if (comp==2) ans *= std::exp(ii*(jzam-half)*R.y);
			if (comp==3) ans *= std::exp(ii*(jzam+half)*R.y);
		}
	}
	else
	{
		if (H.form==qo::schroedinger || H.form==qo::klein_gordon)
		{
			ans *= SphericalHarmonic(Lam,jzam,R.y,R.z);
		}
		if (H.form==qo::pauli || H.form==qo::dirac)
		{
			const tw::Float Lamp = 2*Jam - Lam;
			if (comp==0)
				ans *= SphericalHarmonicSpinor(Jam,Lam,jzam,R.y,R.z)[0];
			if (comp==1)
				ans *= SphericalHarmonicSpinor(Jam,Lam,jzam,R.y,R.z)[1];
			if (comp==2)
				ans *= SphericalHarmonicSpinor(Jam,Lamp,jzam,R.y,R.z)[0];
			if (comp==3)
				ans *= SphericalHarmonicSpinor(Jam,Lamp,jzam,R.y,R.z)[1];
		}
	}
	return amplitude*std::exp(-ii*energy*t)*ans;
}
