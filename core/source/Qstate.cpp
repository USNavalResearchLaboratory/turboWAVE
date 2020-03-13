#include "meta_base.h"
#include "computeTool.h"
#include "qstate.h"

QState::QState(const std::string& name,MetricSpace *ms,Task *tsk) : ComputeTool(name,ms,tsk)
{
	typeCode = tw::tool_type::none;
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

tw::Complex QState::Amplitude(const HamiltonianParameters& H,const tw::vec3& r,const tw::Float& t,const tw::Int& comp) const
{
	return tw::Complex(0.0,0.0);
}

void QState::ReadCheckpoint(std::ifstream& inFile)
{
	ComputeTool::ReadCheckpoint(inFile);
	inFile.read((char*)&amplitude,sizeof(tw::Complex));
	inFile.read((char*)&cylindricalAtom,sizeof(bool));
}

void QState::WriteCheckpoint(std::ofstream& outFile)
{
	ComputeTool::WriteCheckpoint(outFile);
	outFile.write((char*)&amplitude,sizeof(tw::Complex));
	outFile.write((char*)&cylindricalAtom,sizeof(bool));
}


RandomState::RandomState(const std::string& name,MetricSpace *ms,Task *tsk) : QState(name,ms,tsk)
{
	typeCode = tw::tool_type::randomState;
	directives.Add("size",new tw::input::Vec3(&size));
}

tw::Complex RandomState::Amplitude(const HamiltonianParameters& H,const tw::vec3& r,const tw::Float& t,const tw::Int& comp) const
{
	tw::Complex ans;
	ans = task->uniformDeviate->Next() * exp( -sqr(r.x/size.x) - sqr(r.y/size.y) - sqr(r.z/size.z) );
	ans *= amplitude * std::exp( ii*task->uniformDeviate->Next()*two*pi );
	return ans;
}

void RandomState::ReadCheckpoint(std::ifstream& inFile)
{
	QState::ReadCheckpoint(inFile);
	inFile.read((char*)&size,sizeof(tw::vec3));
}

void RandomState::WriteCheckpoint(std::ofstream& outFile)
{
	QState::WriteCheckpoint(outFile);
	outFile.write((char*)&size,sizeof(tw::vec3));
}


FreeState::FreeState(const std::string& name,MetricSpace *ms,Task *tsk) : QState(name,ms,tsk)
{
	typeCode = tw::tool_type::freeState;
	directives.Add("size",new tw::input::Vec3(&size));
	directives.Add("k4",new tw::input::Vec4(&k4));
	directives.Add("spin",new tw::input::Vec3(&spin));
}

tw::Complex FreeState::Amplitude(const HamiltonianParameters& H,const tw::vec3& r,const tw::Float& t,const tw::Int& comp) const
{
	const tw::spinor w(spin/Magnitude(spin),0.5);
	const tw::vec3 k = k4.spatial();
	const bool relativistic = H.form==qo::klein_gordon || H.form==qo::dirac;
	const tw::Float energy = (relativistic ? sqrt((k^k)+H.morb*H.morb) : 0.5*(k^k)/H.morb);
	tw::Complex ans(amplitude);
	ans *= exp( -sqr(r.x/size.x) - sqr(r.y/size.y) - sqr(r.z/size.z) );
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

void FreeState::ReadCheckpoint(std::ifstream& inFile)
{
	QState::ReadCheckpoint(inFile);
	inFile.read((char*)&k4,sizeof(tw::vec4));
	inFile.read((char*)&size,sizeof(tw::vec3));
	inFile.read((char*)&spin,sizeof(spin));
}

void FreeState::WriteCheckpoint(std::ofstream& outFile)
{
	QState::WriteCheckpoint(outFile);
	outFile.write((char*)&k4,sizeof(tw::vec4));
	outFile.write((char*)&size,sizeof(tw::vec3));
	outFile.write((char*)&spin,sizeof(spin));
}


BoundState::BoundState(const std::string& name,MetricSpace *ms,Task *tsk) : QState(name,ms,tsk)
{
	typeCode = tw::tool_type::boundState;
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
				energy = morb/sqrt(1.0 + sqr(qnuc*qorb)/sqr(nr+0.5+sqrt(sqr(jzam)-sqr(qnuc*qorb))));
			else
				energy = morb/sqrt(1.0 + sqr(qnuc*qorb)/sqr(nr+0.5+sqrt(sqr(Lam+0.5)-sqr(qnuc*qorb))));
			break;
		case qo::dirac:
			if (cylindricalAtom)
				energy = morb/sqrt(1.0 + sqr(qnuc*qorb)/sqr(nr+sqrt(sqr(jzam)-sqr(qnuc*qorb))));
			else
				energy = morb/sqrt(1.0 + sqr(qnuc*qorb)/sqr(nr+sqrt(sqr(Jam+0.5)-sqr(qnuc*qorb))));
			break;
	}
	return energy;
}

tw::Float BoundState::NormalizationConstant(const HamiltonianParameters& H) const
{
	tw::Float normalizationConstant;
	const tw::Float n = nr + Lam + 1.0;
	const tw::Float kr = H.morb*H.qnuc*H.qorb/n;
	switch (H.form)
	{
		case qo::schroedinger:
			normalizationConstant = pow(kr,1.5)*sqrt(Factorial(n+Lam)/Factorial(nr))/(Factorial(2.0*Lam+1.0)*sqrt(n*pi));
			break;
		case qo::pauli:
			normalizationConstant = pow(kr,1.5)*sqrt(Factorial(n+Lam)/Factorial(nr))/(Factorial(2.0*Lam+1.0)*sqrt(n*pi));
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
			tw::Float kr = sqrt(-two*morb*energy);
			ans = ConfluentHypergeometric(-tw::Int(nr),two*(Lam+one),two*kr*R.x)*pow(two*kr*R.x,Lam)*exp(-kr*R.x);
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
			tw::Float kr = sqrt(-two*morb*energy);
			ans = ConfluentHypergeometric(-tw::Int(nr),two*(Lam+one),two*kr*R.x)*pow(two*kr*R.x,Lam)*exp(-kr*R.x);
			ans *= SphericalHarmonic(Lam,jzam,R.y,R.z);
			if (comp==0)
				ans *= Jam-Lam+0.5;
			if (comp==1)
				ans *= 0.5-Jam+Lam;
		}
	}

	if (H.form==qo::klein_gordon)
	{
		tw::Float kr = sqrt(morb*morb - sqr(energy));
		if (cylindricalAtom)
		{
			space->CurvilinearToCylindrical(&R);
			const tw::Float gam = sqrt(sqr(jzam) - sqr(qnuc*qorb));
			ans = ConfluentHypergeometric(-tw::Int(nr),two*gam+one,two*kr*R.x);
			ans *= pow(R.x,gam)*exp(-kr*R.x);
			ans *= std::exp(ii*jzam*R.y);
			if (comp==1)
				ans *= (energy - qnuc*qorb/R.x)/morb;
		}
		else
		{
			space->CurvilinearToSpherical(&R);
			const tw::Float gam = sqrt(sqr(Lam+0.5)-sqr(qnuc*qorb));
			ans = ConfluentHypergeometric(-tw::Int(nr),two*gam+one,two*kr*R.x);
			ans *= pow(R.x,gam-0.5)*exp(-kr*R.x);
			ans *= SphericalHarmonic(Lam,jzam,R.y,R.z);
			if (comp==1)
				ans *= (energy - qnuc*qorb/R.x)/morb;
		}
	}

	if (H.form==qo::dirac)
	{
		tw::Complex Q1,Q2;
		tw::Float kr = sqrt(morb*morb - sqr(energy));
		if (cylindricalAtom)
		{
			space->CurvilinearToCylindrical(&R);
			const tw::Float szam = jzam-Lam;
			const tw::Float kappa = -2*szam*jzam;
			const tw::Float gam = sqrt(sqr(kappa) - sqr(qnuc*qorb));
			Q1 = ConfluentHypergeometric(-tw::Int(nr),two*gam+one,two*kr*R.x);
			Q2 = ConfluentHypergeometric(1-tw::Int(nr),two*gam+one,two*kr*R.x);
			Q2 *= -(kappa - qnuc*qorb*morb/kr) / (gam - qnuc*qorb*energy/kr);
			tw::Complex f = tw::Float(sqrt(morb+energy)*exp(-kr*R.x)*pow(R.x,gam-0.5))*(Q1+Q2);
			tw::Complex g = ii*tw::Float(sqrt(morb-energy)*exp(-kr*R.x)*pow(R.x,gam-0.5))*(Q1-Q2);
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
			const tw::Float gam = sqrt(sqr(kappa)-sqr(qorb*qnuc));
			Q1 = ConfluentHypergeometric(-tw::Int(nr),two*gam+one,two*kr*R.x);
			Q2 = ConfluentHypergeometric(1-tw::Int(nr),two*gam+one,two*kr*R.x);
			Q2 *= -(gam - qnuc*qorb*energy/kr)/(kappa - qnuc*qorb*morb/kr);
			tw::Complex f = tw::Float(sqrt(morb+energy)*exp(-kr*R.x)*pow(2.0*kr*R.x,gam-1.0))*(Q1+Q2);
			tw::Complex g = tw::Float(sqrt(morb-energy)*exp(-kr*R.x)*pow(2.0*kr*R.x,gam-1.0))*(Q1-Q2);
			ans = 0.0;
			if (comp==0)
				ans = f*SphericalHarmonicSpinor(Jam,Lam,jzam,R.y,R.z)[0];
			if (comp==1)
				ans = f*SphericalHarmonicSpinor(Jam,Lam,jzam,R.y,R.z)[1];
			if (comp==2)
				ans = tw::Float(pow(-1.0,0.5*(1.0+Lam-Lamp)))*g*SphericalHarmonicSpinor(Jam,Lamp,jzam,R.y,R.z)[0];
			if (comp==3)
				ans = tw::Float(pow(-1.0,0.5*(1.0+Lam-Lamp)))*g*SphericalHarmonicSpinor(Jam,Lamp,jzam,R.y,R.z)[1];
		}
	}

	return normalizationConstant*amplitude*std::exp(-ii*energy*t)*ans;
}

void BoundState::ReadCheckpoint(std::ifstream& inFile)
{
	QState::ReadCheckpoint(inFile);
	inFile.read((char*)&Lam,sizeof(tw::Float));
	inFile.read((char*)&Jam,sizeof(tw::Float));
	inFile.read((char*)&jzam,sizeof(tw::Float));
	inFile.read((char*)&nr,sizeof(tw::Float));
}

void BoundState::WriteCheckpoint(std::ofstream& outFile)
{
	QState::WriteCheckpoint(outFile);
	outFile.write((char*)&Lam,sizeof(tw::Float));
	outFile.write((char*)&Jam,sizeof(tw::Float));
	outFile.write((char*)&jzam,sizeof(tw::Float));
	outFile.write((char*)&nr,sizeof(tw::Float));
}


TabulatedState::TabulatedState(const std::string& name,MetricSpace *ms,Task *tsk) : QState(name,ms,tsk)
{
	typeCode = tw::tool_type::tabulatedState;
	directives.Add("filename",new tw::input::String(&filename));
}

void TabulatedState::Initialize()
{
	QState::Initialize();
	tw::Int num;
	std::string word;
	tw::Float realAmplitude,imagAmplitude;
	std::ifstream inFile(filename.c_str());
	if (inFile.rdstate() & std::ios::failbit)
		throw tw::FatalError("couldn't open state file <" + filename + ">");
	do
	{
		inFile >> word;
		//if (word=="soft_core_radius")
		//if (word=="nuclear_charge")
		//if (word=="Bz")
		if (word=="energy")
			inFile >> word >> energy;
		if (word=="pts")
			inFile >> word >> num;
		if (word=="components")
			inFile >> word >> components;
		if (word=="cell_width")
			inFile >> word >> cell_size.x;
		if (word=="nr_j_l_m")
			inFile >> word >> nr >> Jam >> Lam >> jzam;
		if (word=="cylindrical")
		{
			inFile >> word >> word;
			cylindricalAtom = (word=="true" || word=="yes" || word=="on") ? true : false;
		}
	} while (word!="start_data");

	radialFunction.resize(components*num);
	for (tw::Int i=0;i<components*num;i++)
	{
		inFile >> realAmplitude >> imagAmplitude;
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

void TabulatedState::ReadCheckpoint(std::ifstream& inFile)
{
	QState::ReadCheckpoint(inFile);
	inFile >> filename;
	inFile.ignore();
}

void TabulatedState::WriteCheckpoint(std::ofstream& outFile)
{
	QState::WriteCheckpoint(outFile);
	outFile << filename << " ";
}
