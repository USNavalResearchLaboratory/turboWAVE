#include "simulation.h"
#include "fieldSolve.h"
#include "quantum.h"
using namespace tw::bc;

qo::State::State(MetricSpace *m,UniformDeviate *u)
{
	space = m;
	ud = u;
	waveType = qo::bound;
	waveEquation = qo::schroedinger;
	Lam = 0.0;
	Jam = jzam = 0.0;
	nr = 0.0;
	amplitude = tw::Complex(1.0,0.0);
	cell_size = tw::vec3(1.0,1.0,1.0);
	size = tw::vec3(1.0,1.0,1.0);
	k = tw::vec3(0.0,0.0,1.0);
	negativeEnergy = false;
	helicity = 0.5;
	qnuc = 1.0;
	qorb = -1.0;
	morb = 1.0;
	if (m->Dim(3)==1)
		cylindricalAtom = true;
	else
		cylindricalAtom = false;
	energy = -0.5;
	normalizationConstant = 1.0;
	components = 1;
}

bool qo::State::GoodQuantumNumbers()
{
	tw::Float int_part;
	auto good_lam = [&] (tw::Float test) { return std::modf(test,&int_part)==0.0; };
	auto good_jam = [&] (tw::Float test) { return std::abs(std::modf(test,&int_part))==0.5; };
	if (waveEquation==qo::schroedinger || waveEquation==qo::klein_gordon)
	{
		if (cylindricalAtom)
			return good_lam(jzam) && jzam==Lam;
		return Jam>=0.0 && good_lam(Jam) && good_lam(jzam) && Jam==Lam && jzam>=-Jam && jzam<=Jam;
	}
	if (waveEquation==qo::pauli || waveEquation==qo::dirac)
	{
		if (cylindricalAtom)
			return good_jam(jzam) && std::abs(Lam-jzam)==0.5 && ((nr==0.0 && (Lam-jzam)*jzam<0.0) || nr>0.0);
		return Jam>0.0 && good_jam(Jam) && good_jam(jzam) && std::abs(Jam-Lam)==0.5 && jzam>=-Jam && jzam<=Jam && ((nr==0.0 && Lam<Jam) || nr>0.0);
	}
	return true;
}

void qo::State::Initialize(qo::waveEquationType waveEq,tw::Float nuclearCharge,tw::Float orbitingCharge,tw::Float orbitingMass)
{
	waveEquation = waveEq;
	qnuc = nuclearCharge;
	qorb = orbitingCharge;
	morb = orbitingMass;
	const tw::Float n = nr + Lam + 1.0;
	const tw::Float kr = morb*qnuc*qorb/n;

	switch (waveType)
	{
		case qo::random:
			energy = 0.0;
			break;
		case qo::lookup:
			if (!GoodQuantumNumbers())
				throw tw::FatalError("Quantum numbers are invalid.");
			break;
		case qo::bound:
			if (!GoodQuantumNumbers())
				throw tw::FatalError("Quantum numbers are invalid.");
			switch (waveEquation)
			{
				case schroedinger:
					energy = -0.5*morb*sqr(qnuc*qorb/n);
					normalizationConstant = pow(kr,1.5)*sqrt(Factorial(n+Lam)/Factorial(nr))/(Factorial(2.0*Lam+1.0)*sqrt(n*pi));
					break;
				case pauli:
					energy = -0.5*morb*sqr(qnuc*qorb/n);
					normalizationConstant = pow(kr,1.5)*sqrt(Factorial(n+Lam)/Factorial(nr))/(Factorial(2.0*Lam+1.0)*sqrt(n*pi));
					break;
				case klein_gordon:
					if (cylindricalAtom)
						energy = morb/sqrt(1.0 + sqr(qnuc*qorb)/sqr(nr+0.5+sqrt(sqr(jzam)-sqr(qnuc*qorb))));
					else
						energy = morb/sqrt(1.0 + sqr(qnuc*qorb)/sqr(nr+0.5+sqrt(sqr(Lam+0.5)-sqr(qnuc*qorb))));
					break;
				case dirac:
					if (cylindricalAtom)
						energy = morb/sqrt(1.0 + sqr(qnuc*qorb)/sqr(nr+sqrt(sqr(jzam)-sqr(qnuc*qorb))));
					else
						energy = morb/sqrt(1.0 + sqr(qnuc*qorb)/sqr(nr+sqrt(sqr(Jam+0.5)-sqr(qnuc*qorb))));
					break;
			}
			break;
		case qo::free:
		case qo::helicity:
			switch (waveEquation)
			{
				case schroedinger:
					energy = 0.5*(k^k)/morb;
					break;
				case pauli:
					energy = 0.5*(k^k)/morb;
					break;
				case klein_gordon:
					energy = sqrt((k^k) + morb*morb);
					break;
				case dirac:
					energy = sqrt((k^k) + morb*morb);
					break;
			}
			break;
	}
}

tw::Complex qo::State::Random(const tw::vec3& r,const tw::Float& t,const tw::Int& comp) const
{
	tw::Complex ans;
	ans = ud->Next() * exp( -sqr(r.x/size.x) - sqr(r.y/size.y) - sqr(r.z/size.z) );
	ans *= amplitude * std::exp( ii*ud->Next()*two*pi );
	return normalizationConstant*ans;
}

tw::Complex qo::State::Lookup(const tw::vec3& r,const std::valarray<tw::Complex>& fr,const tw::Float& dr,const tw::Float& t,const tw::Int& comp) const
{
	tw::Complex ans;
	tw::vec3 R = r;
	if (cylindricalAtom)
		space->CurvilinearToCylindrical(&R);
	else
		space->CurvilinearToSpherical(&R);
	tw::Int dim = fr.size()/components;
	tw::Int s = MyFloor((R.x-0.5*dr)/dr);
	tw::Float w = (R.x-0.5*dr)/dr - tw::Float(s);
	if (s>=0 && s<dim-1)
		ans = fr[comp*dim+s]*(one-w) + fr[comp*dim+s+1]*w;
	if (s<0)
		ans = fr[comp*dim];
	if (s==dim-1)
		ans = fr[comp*dim+dim-1];
	if (cylindricalAtom)
	{
		if (waveEquation==schroedinger || waveEquation==klein_gordon)
		{
			ans *= std::exp(ii*jzam*R.y);
		}
		if (waveEquation==pauli || waveEquation==dirac)
		{
			if (comp==0) ans *= std::exp(ii*(jzam-half)*R.y);
			if (comp==1) ans *= std::exp(ii*(jzam+half)*R.y);
			if (comp==2) ans *= std::exp(ii*(jzam-half)*R.y);
			if (comp==3) ans *= std::exp(ii*(jzam+half)*R.y);
		}
	}
	else
	{
		if (waveEquation==schroedinger || waveEquation==klein_gordon)
		{
			ans *= SphericalHarmonic(Lam,jzam,R.y,R.z);
		}
		if (waveEquation==pauli || waveEquation==dirac)
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
	return normalizationConstant*amplitude*std::exp(-ii*energy*t)*ans;
}

tw::Complex qo::State::Free(bool helicityState,const tw::vec3& r,const tw::Float& t,const tw::Int& comp) const
{
	tw::spinor w;
	if (helicityState)
		w = tw::spinor(k/Magnitude(k),helicity);
	else
		w = tw::spinor(0.5+helicity,0.5-helicity);
	tw::Complex ans(amplitude);
	ans *= exp( -sqr(r.x/size.x) - sqr(r.y/size.y) - sqr(r.z/size.z) );
	ans *= std::exp( ii*((k^r) - energy*t) );
	ans *= normalizationConstant;

	switch (waveEquation)
	{
		case schroedinger:
			// nothing to do
			break;
		case pauli:
			ans *= w[comp];
			break;
		case klein_gordon:
			// nothing to do for positive energy
			// need to put something here for negative energy
			break;
		case dirac:
			// Create an electron bispinor coefficient in the rest frame, L&L QED 23
			tw::bispinor u(w);
			// Boost to desired momentum
			u.Boost(tw::vec4(energy,k));
			// If negative energy apply charge conjugation operator
			// This reverses momentum and spin
			if (negativeEnergy)
			{
				ans = std::conj(ans);
				u.ChargeConjugate();
			}
			ans *= u[comp];
			break;
	}
	return ans;
}

tw::Complex qo::State::Bound(const tw::vec3& r,const tw::Float& t,const tw::Int& comp) const
{
	tw::Complex ans;
	tw::vec3 R = r;

	if (waveEquation==schroedinger)
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

	if (waveEquation==pauli)
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

	if (waveEquation==klein_gordon)
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

	if (waveEquation==dirac)
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

tw::Complex qo::State::Amplitude(const tw::vec3& r,const tw::Float& t,const tw::Int& comp) const
{
	switch (waveType)
	{
		case qo::lookup:
			return Lookup(r,radialFunction,cell_size.x,t,comp);
			break;
		case qo::bound:
			return Bound(r,t,comp);
			break;
		case qo::free:
			return Free(false,r,t,comp);
			break;
		case qo::helicity:
			return Free(true,r,t,comp);
			break;
		case qo::random:
			return Random(r,t,comp);
			break;
	}
	return tw::Complex(0.0,0.0);
}

void qo::State::ReadInputFileBlock(std::stringstream& inputString)
{
	std::string command,word;
	tw::Float realAmplitude,imagAmplitude;
	do
	{
		inputString >> command;
		if (command=="type")
		{
			inputString >> word >> word;
			if (word=="lookup")
				waveType = qo::lookup;
			if (word=="free")
				waveType = qo::free;
			if (word=="helicity")
				waveType = qo::helicity;
			if (word=="bound")
				waveType = qo::bound;
			if (word=="random")
				waveType = qo::random;
		}
		if (command=="file")
		{
			inputString >> word >> word;
			std::ifstream inFile(word.c_str());
			tw::Int i,num;
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
			for (i=0;i<components*num;i++)
			{
				inFile >> realAmplitude >> imagAmplitude;
				radialFunction[i] = tw::Complex(realAmplitude,imagAmplitude);
			}
		}
		if (command=="dr")
			throw tw::FatalError("state specification 'dr' no longer supported, specify in file");
		if (command=="amplitude")
		{
			inputString >> word >> realAmplitude >> imagAmplitude;
			amplitude = tw::Complex(realAmplitude,imagAmplitude);
		}
		if (command=="nr_j_l_m")
		{
			inputString >> word >> nr >> Jam >> Lam >> jzam;
		}
		if (command=="nrlz" || command=="nrllz" || command=="lam" || command=="jam" || command=="zam" || command=="nr")
		{
			throw tw::FatalError("Quantum state specification is deprecated.");
		}
		if (command=="size")
		{
			inputString >> word >> size.x >> size.y >> size.z;
		}
		if (command=="k")
			throw tw::FatalError("Free state specification is deprecated.");
		if (command=="k_e_s")
		{
			inputString >> word >> k.x >> k.y >> k.z >> energy >> helicity;
			negativeEnergy = (energy<0.0);
		}
		if (command=="cylindrical")
		{
			inputString >> word >> word;
			cylindricalAtom = (word=="true" || word=="yes" || word=="on") ? true : false;
		}
	} while (command!="}");
}

void qo::State::ReadData(std::ifstream& inFile)
{
	tw::Int i,num;
	inFile.read((char*)&waveType,sizeof(waveFunctionType));
	inFile.read((char*)&waveEquation,sizeof(waveEquationType));
	inFile.read((char*)&amplitude,sizeof(tw::Complex));
	inFile.read((char*)&Lam,sizeof(tw::Float));
	inFile.read((char*)&Jam,sizeof(tw::Float));
	inFile.read((char*)&jzam,sizeof(tw::Float));
	inFile.read((char*)&nr,sizeof(tw::Float));
	inFile.read((char*)&k,sizeof(tw::vec3));
	inFile.read((char*)&negativeEnergy,sizeof(bool));
	inFile.read((char*)&helicity,sizeof(tw::Float));
	inFile.read((char*)&size,sizeof(tw::vec3));
	inFile.read((char*)&cell_size,sizeof(tw::vec3));
	inFile.read((char*)&qnuc,sizeof(tw::Float));
	inFile.read((char*)&qorb,sizeof(tw::Float));
	inFile.read((char*)&morb,sizeof(tw::Float));
	inFile.read((char*)&energy,sizeof(tw::Float));
	inFile.read((char*)&normalizationConstant,sizeof(tw::Float));
	inFile.read((char*)&cylindricalAtom,sizeof(bool));

	inFile.read((char*)&num,sizeof(tw::Int));
	inFile.read((char*)&components,sizeof(tw::Int));
	radialFunction.resize(num);
	for (i=0;i<num;i++)
		inFile.read((char*)&radialFunction[i],sizeof(tw::Float));
}

void qo::State::WriteData(std::ofstream& outFile)
{
	tw::Int i,num;
	outFile.write((char*)&waveType,sizeof(waveFunctionType));
	outFile.write((char*)&waveEquation,sizeof(waveEquationType));
	outFile.write((char*)&amplitude,sizeof(tw::Complex));
	outFile.write((char*)&Lam,sizeof(tw::Float));
	outFile.write((char*)&Jam,sizeof(tw::Float));
	outFile.write((char*)&jzam,sizeof(tw::Float));
	outFile.write((char*)&nr,sizeof(tw::Float));
	outFile.write((char*)&k,sizeof(tw::vec3));
	outFile.write((char*)&negativeEnergy,sizeof(bool));
	outFile.write((char*)&helicity,sizeof(tw::Float));
	outFile.write((char*)&size,sizeof(tw::vec3));
	outFile.write((char*)&cell_size,sizeof(tw::vec3));
	outFile.write((char*)&qnuc,sizeof(tw::Float));
	outFile.write((char*)&qorb,sizeof(tw::Float));
	outFile.write((char*)&morb,sizeof(tw::Float));
	outFile.write((char*)&energy,sizeof(tw::Float));
	outFile.write((char*)&normalizationConstant,sizeof(tw::Float));
	outFile.write((char*)&cylindricalAtom,sizeof(bool));

	num = radialFunction.size();
	outFile.write((char*)&num,sizeof(tw::Int));
	outFile.write((char*)&components,sizeof(tw::Int));
	for (i=0;i<num;i++)
		outFile.write((char*)&radialFunction[i],sizeof(tw::Float));
}


////////////////////////////////
//                            //
//   ATOMIC PHYSICS MODULE    //
//                            //
////////////////////////////////


AtomicPhysics::AtomicPhysics(const std::string& name,Simulation* sim):Module(name,sim)
{
	typeCode = tw::module_type::nullModule;
	keepA2Term = true;
	dipoleApproximation = true;
	alpha = 0.0072973525664;
	timeRelaxingToGround = 0.0;

	potentialTypeSpec = qo::softCore;
	nuclearRadius = 0.01;
	residualCharge = 1.0;
	B0 = 0.0;
	m0 = 1.0;
	q0 = -1.0;

	A4.Initialize(4,*this,owner,tw::dom::xAxis);
	Ao4.Initialize(4,*this,owner,tw::dom::xAxis);
	J4.Initialize(4,*this,owner,tw::dom::xAxis);

	photonPropagator = NULL;

	#ifdef USE_OPENCL
	InitializeCLProgram("quantum.cl");
	A4.InitializeComputeBuffer();
	Ao4.InitializeComputeBuffer();
	J4.InitializeComputeBuffer();
	#endif

	directives.Add("orbiting charge",new tw::input::Float(&q0));
	directives.Add("orbiting mass",new tw::input::Float(&m0));
	directives.Add("B0",new tw::input::Vec3(&B0));
	directives.Add("keep a2 term",new tw::input::Bool(&keepA2Term));
	directives.Add("dipole approximation",new tw::input::Bool(&dipoleApproximation));
	directives.Add("relaxation time",new tw::input::Float(&timeRelaxingToGround));
	directives.Add("soft core potential charge",new tw::input::Custom);
	directives.Add("bachelet potential",new tw::input::Custom);
}

AtomicPhysics::~AtomicPhysics()
{
	if (photonPropagator!=NULL)
		owner->RemoveTool(photonPropagator);
	// release any base class OpenCL kernels here
}

tw::Float AtomicPhysics::GetSphericalPotential(tw::Float r) const
{
	tw::Int j;
	tw::Float ans;
	if (potentialTypeSpec==qo::softCore)
		ans = residualCharge/sqrt(sqr(nuclearRadius) + sqr(r));
	if (potentialTypeSpec==qo::bachelet)
	{
		ans = (residualCharge / r) * (c1*erf(sqrt(a1)*r) + c2*erf(sqrt(a2)*r));
		for (j=0;j<3;j++)
			ans -= (Asr[j] + r*r*Asr[j+3])*exp(-asr[j]*r*r);
	}
	return ans;
}

void AtomicPhysics::Initialize()
{
	Module::Initialize();
	#ifdef USE_OPENCL
	photonPropagator->SetupComputeKernels(A4,Ao4,J4);
	#endif

	// Boundary conditions should preserve hermiticity
	// One way is to have A = 0 and grad(psi)=0 for components normal to boundary

	tw::bc::fld psiDefaultBC,A4DefaultBC;
	psiDefaultBC = fld::neumannWall;
	A4DefaultBC = fld::dirichletWall;
	psi_r.SetBoundaryConditions(tw::dom::xAxis,psiDefaultBC,psiDefaultBC);
	psi_r.SetBoundaryConditions(tw::dom::yAxis,psiDefaultBC,psiDefaultBC);
	psi_r.SetBoundaryConditions(tw::dom::zAxis,psiDefaultBC,psiDefaultBC);
	psi_i.SetBoundaryConditions(tw::dom::xAxis,psiDefaultBC,psiDefaultBC);
	psi_i.SetBoundaryConditions(tw::dom::yAxis,psiDefaultBC,psiDefaultBC);
	psi_i.SetBoundaryConditions(tw::dom::zAxis,psiDefaultBC,psiDefaultBC);
	A4.SetBoundaryConditions(tw::dom::xAxis,A4DefaultBC,A4DefaultBC);
	A4.SetBoundaryConditions(tw::dom::yAxis,A4DefaultBC,A4DefaultBC);
	A4.SetBoundaryConditions(tw::dom::zAxis,A4DefaultBC,A4DefaultBC);
	J4.SetBoundaryConditions(tw::dom::xAxis,A4DefaultBC,A4DefaultBC);
	J4.SetBoundaryConditions(tw::dom::yAxis,A4DefaultBC,A4DefaultBC);
	J4.SetBoundaryConditions(tw::dom::zAxis,A4DefaultBC,A4DefaultBC);

	switch (owner->gridGeometry)
	{
		case tw::dom::cartesian:
			break;
		case tw::dom::cylindrical:
			A4.SetBoundaryConditions(Element(0),tw::dom::xAxis,fld::neumannWall,fld::dirichletWall);
			A4.SetBoundaryConditions(Element(3),tw::dom::xAxis,fld::neumannWall,fld::dirichletWall);
			J4.SetBoundaryConditions(Element(0),tw::dom::xAxis,fld::neumannWall,fld::dirichletWall);
			J4.SetBoundaryConditions(Element(3),tw::dom::xAxis,fld::neumannWall,fld::dirichletWall);
			break;
		case tw::dom::spherical:
			A4.SetBoundaryConditions(Element(0),tw::dom::yAxis,fld::neumannWall,fld::neumannWall);
			A4.SetBoundaryConditions(Element(1),tw::dom::yAxis,fld::neumannWall,fld::neumannWall);
			J4.SetBoundaryConditions(Element(0),tw::dom::yAxis,fld::neumannWall,fld::neumannWall);
			J4.SetBoundaryConditions(Element(1),tw::dom::yAxis,fld::neumannWall,fld::neumannWall);
			break;
	}

	// Wavefunction and Potentials are initialized in sub classes
}

void AtomicPhysics::ExchangeResources()
{
	PublishResource(&J4,"qo:j4");
}

void AtomicPhysics::FormPotentials(tw::Float t)
{
	// residual charge has to be in desired units (it is Z for atomic units, Z*sqrt(alpha) for natural)
	// vector potential in wave block must also be in desired units, unlike earlier versions which used mc^2/e in every case

	#pragma omp parallel firstprivate(t)
	{
		tw::Float phiNow,r;
		tw::vec3 A0,A1,r_curv,r_cart;
		for (auto cell : EntireCellRange(*this))
		{
			r = owner->SphericalRadius(owner->Pos(cell));
			phiNow = GetSphericalPotential(r);
			r_cart = dipoleApproximation ? tw::vec3(0,0,0) : owner->Pos(cell);
			owner->CurvilinearToCartesian(&r_cart);
			A0 = A1 = tw::vec3(-0.5*r_cart.y*B0.z,0.5*r_cart.x*B0.z,0.0);
			for (tw::Int s=0;s<owner->wave.size();s++)
			{
				A0 += owner->wave[s]->VectorPotential(t-dt,r_cart);
				A1 += owner->wave[s]->VectorPotential(t,r_cart);
			}
			r_curv = owner->Pos(cell);
			owner->TangentVectorToCurvilinear(&A0,r_curv);
			owner->TangentVectorToCurvilinear(&A1,r_curv);

			Ao4(cell,0) = phiNow;
			Ao4(cell,1) = A0.x;
			Ao4(cell,2) = A0.y;
			Ao4(cell,3) = A0.z;

			A4(cell,0) = phiNow;
			A4(cell,1) = A1.x;
			A4(cell,2) = A1.y;
			A4(cell,3) = A1.z;
		}
	}

	Ao4.CopyFromNeighbors();
	//Ao4.ApplyBoundaryCondition();
	A4.CopyFromNeighbors();
	//A4.ApplyBoundaryCondition();
}

void AtomicPhysics::FormGhostCellPotentials(tw::Float t)
{
	for (tw::Int ax=1;ax<=3;ax++)
		if (A4.Dim(ax)>1)
			#pragma omp parallel firstprivate(t,ax)
			{
				for (auto s : StripRange(*this,ax,strongbool::no))
					for (tw::Int ghostCell=0;ghostCell<=Dim(s.Axis())+1;ghostCell+=Dim(s.Axis())+1)
					{
						tw::vec3 pos(owner->Pos(s,ghostCell));
						tw::vec3 A3(-0.5*pos.y*B0.z,0.5*pos.x*B0.z,0.0);
						for (tw::Int wv=0;wv<owner->wave.size();wv++)
							A3 += owner->wave[wv]->VectorPotential(t,pos);
						if ((ghostCell==0 && owner->n0[ax]==MPI_PROC_NULL) || (ghostCell!=0 && owner->n1[ax]==MPI_PROC_NULL))
						{
							A4(s,ghostCell,1) = A3.x;
							A4(s,ghostCell,2) = A3.y;
							A4(s,ghostCell,3) = A3.z;
						}
					}
			}
}

tw::vec4 AtomicPhysics::GetA4AtOrigin()
{
	tw::Int s;
	tw::vec3 aNow,r_cart;
	tw::vec4 A;

	aNow = 0.0;
	r_cart = tw::vec3(0,0,0);
	A[0] = GetSphericalPotential(0.0);
	for (s=0;s<owner->wave.size();s++)
		aNow += owner->wave[s]->VectorPotential(owner->elapsedTime,r_cart);
	A[1] = aNow.x;
	A[2] = aNow.y;
	A[3] = aNow.z;
	return A;
}

void AtomicPhysics::VerifyInput()
{
	Module::VerifyInput();
	if (owner->gridGeometry!=tw::dom::cartesian)
		if (B0.x!=0.0 || B0.y!=0.0 || B0.z!=0.0)
			throw tw::FatalError("Static B field assumes Cartesian geometry.");
	for (auto tool : moduleTool)
	{
		photonPropagator = dynamic_cast<LorentzPropagator*>(tool);
		if (photonPropagator!=NULL)
			break;
	}
	if (photonPropagator==NULL)
		photonPropagator = (LorentzPropagator*)owner->CreateTool("default_photons",tw::tool_type::lorentzPropagator);
}

void AtomicPhysics::ReadInputFileDirective(std::stringstream& inputString,const std::string& command)
{
	std::string word;
	// Must intercept new before Module::ReadInputFileDirective,
	// otherwise submodule interpretation will be made.
	if (command=="new")
	{
		inputString >> word;
		if (word=="wavefunction")
		{
			waveFunction.push_back(qo::State(owner,owner->uniformDeviate));
			waveFunction.back().ReadInputFileBlock(inputString);
		}
		if (word=="reference")
		{
			refState.push_back(qo::State(owner,owner->uniformDeviate));
			refState.back().ReadInputFileBlock(inputString);
		}
	}
	else
		Module::ReadInputFileDirective(inputString,command);
	// note: examples of charge are geared toward atomic units
	// if using natural units, unit of charge is sqrt(alpha) ~ 0.085
	// at present unit conversions are not performed in quantum modules---you have to enter it as it will be used
	if (command=="soft core potential charge") // eg, soft core potential , charge = 1.0 , radius = 0.01
	{
		potentialTypeSpec = qo::softCore;
		inputString >> word >> residualCharge >> word >> word >> nuclearRadius;
	}
	if (command=="bachelet potential") // eg, bachelet potential = ...
	{
		potentialTypeSpec = qo::bachelet;
		inputString >> word;
		inputString >> residualCharge;
		inputString >> c1 >> c2 >> a1 >> a2;
		inputString >> Asr[0] >> Asr[1] >> Asr[2] >> Asr[3] >> Asr[4] >> Asr[5];
		inputString >> asr[0] >> asr[1] >> asr[2];
	}
}

void AtomicPhysics::ReadData(std::ifstream& inFile)
{
	tw::Int num;

	Module::ReadData(inFile);
	inFile.read((char *)&potentialTypeSpec,sizeof(qo::potentialType));
	inFile.read((char *)&residualCharge,sizeof(tw::Float));
	inFile.read((char *)&nuclearRadius,sizeof(tw::Float));
	inFile.read((char *)&c1,sizeof(tw::Float));
	inFile.read((char *)&c2,sizeof(tw::Float));
	inFile.read((char *)&a1,sizeof(tw::Float));
	inFile.read((char *)&a2,sizeof(tw::Float));
	inFile.read((char *)Asr,6*sizeof(tw::Float));
	inFile.read((char *)asr,3*sizeof(tw::Float));
	inFile.read((char *)&B0,sizeof(tw::vec3));
	inFile.read((char *)&q0,sizeof(tw::Float));
	inFile.read((char *)&m0,sizeof(tw::Float));
	inFile.read((char *)&keepA2Term,sizeof(bool));
	inFile.read((char *)&dipoleApproximation,sizeof(bool));
	inFile.read((char *)&timeRelaxingToGround,sizeof(tw::Float));
	inFile.read((char *)&groundStateEnergy,sizeof(tw::Float));

	inFile.read((char *)&num,sizeof(tw::Int));
	for (tw::Int i=0;i<num;i++)
	{
		waveFunction.push_back(qo::State(owner,owner->uniformDeviate));
		inFile.read((char *)&waveFunction.back(),sizeof(qo::State));
		waveFunction.back().space = owner;
		waveFunction.back().ud = owner->uniformDeviate;
	}

	inFile.read((char *)&num,sizeof(tw::Int));
	for (tw::Int i=0;i<num;i++)
	{
		refState.push_back(qo::State(owner,owner->uniformDeviate));
		inFile.read((char *)&refState.back(),sizeof(qo::State));
		refState.back().space = owner;
		refState.back().ud = owner->uniformDeviate;
	}

	psi_r.ReadData(inFile);
	psi_i.ReadData(inFile);
	J4.ReadData(inFile);
	Ao4.ReadData(inFile);
	A4.ReadData(inFile);
}

void AtomicPhysics::WriteData(std::ofstream& outFile)
{
	tw::Int num;

	Module::WriteData(outFile);
	outFile.write((char *)&potentialTypeSpec,sizeof(qo::potentialType));
	outFile.write((char *)&residualCharge,sizeof(tw::Float));
	outFile.write((char *)&nuclearRadius,sizeof(tw::Float));
	outFile.write((char *)&c1,sizeof(tw::Float));
	outFile.write((char *)&c2,sizeof(tw::Float));
	outFile.write((char *)&a1,sizeof(tw::Float));
	outFile.write((char *)&a2,sizeof(tw::Float));
	outFile.write((char *)Asr,6*sizeof(tw::Float));
	outFile.write((char *)asr,3*sizeof(tw::Float));
	outFile.write((char *)&B0,sizeof(tw::vec3));
	outFile.write((char *)&q0,sizeof(tw::Float));
	outFile.write((char *)&m0,sizeof(tw::Float));
	outFile.write((char *)&keepA2Term,sizeof(bool));
	outFile.write((char *)&dipoleApproximation,sizeof(bool));
	outFile.write((char *)&timeRelaxingToGround,sizeof(tw::Float));
	outFile.write((char *)&groundStateEnergy,sizeof(tw::Float));

	num = waveFunction.size();
	outFile.write((char *)&num,sizeof(tw::Int));
	for (tw::Int i=0;i<num;i++)
		outFile.write((char *)&waveFunction[i],sizeof(qo::State));

	num = refState.size();
	outFile.write((char *)&num,sizeof(tw::Int));
	for (tw::Int i=0;i<num;i++)
		outFile.write((char *)&refState[i],sizeof(qo::State));

	psi_r.WriteData(outFile);
	psi_i.WriteData(outFile);
	J4.WriteData(outFile);
	Ao4.WriteData(outFile);
	A4.WriteData(outFile);
}


///////////////////////////////////////
//                                   //
// NON-RELATIVISTIC PARTICLE (TDSE)  //
//                                   //
///////////////////////////////////////


Schroedinger::Schroedinger(const std::string& name,Simulation* sim):AtomicPhysics(name,sim)
{
	// Should move OpenCL stuff into the propagator tool
	typeCode = tw::module_type::schroedinger;
	#ifndef USE_OPENCL
	propagator = NULL;
	#endif
	psi0.Initialize(*this,owner);
	psi1.Initialize(*this,owner);
	v.Initialize(*this,owner);
	w.Initialize(*this,owner);
	scratch.Initialize(*this,owner);

	#ifdef USE_OPENCL

	cl_int err;
	applyNumerator = clCreateKernel(program,"ApplyNumerator",&err);
	applyDenominator = clCreateKernel(program,"ApplyDenominator",&err);
	chargeKernel = clCreateKernel(program,"DepositCharge",&err);
	currentKernel = clCreateKernel(program,"DepositCurrent",&err);
	stitchKernel = clCreateKernel(program,"Stitch",&err);

	scratch.InitializeComputeBuffer();
	v.InitializeComputeBuffer();
	w.InitializeComputeBuffer();
	psi0.InitializeComputeBuffer();
	psi1.InitializeComputeBuffer();
	J4.InitializeComputeBuffer();

	clSetKernelArg(applyNumerator,0,sizeof(cl_mem),&psi1.computeBuffer);
	clSetKernelArg(applyNumerator,1,sizeof(cl_mem),&A4.computeBuffer);
	clSetKernelArg(applyNumerator,2,sizeof(cl_mem),&owner->metricsBuffer);

	clSetKernelArg(applyDenominator,0,sizeof(cl_mem),&psi1.computeBuffer);
	clSetKernelArg(applyDenominator,1,sizeof(cl_mem),&A4.computeBuffer);
	clSetKernelArg(applyDenominator,2,sizeof(cl_mem),&owner->metricsBuffer);
	clSetKernelArg(applyDenominator,6,sizeof(cl_mem),&scratch.computeBuffer);
	clSetKernelArg(applyDenominator,7,sizeof(cl_mem),&v.computeBuffer);
	clSetKernelArg(applyDenominator,8,sizeof(cl_mem),&w.computeBuffer);

	clSetKernelArg(chargeKernel,0,sizeof(cl_mem),&psi1.computeBuffer);
	clSetKernelArg(chargeKernel,1,sizeof(cl_mem),&J4.computeBuffer);

	clSetKernelArg(currentKernel,0,sizeof(cl_mem),&psi0.computeBuffer);
	clSetKernelArg(currentKernel,1,sizeof(cl_mem),&psi1.computeBuffer);
	clSetKernelArg(currentKernel,2,sizeof(cl_mem),&A4.computeBuffer);
	clSetKernelArg(currentKernel,3,sizeof(cl_mem),&J4.computeBuffer);
	clSetKernelArg(currentKernel,4,sizeof(cl_mem),&owner->metricsBuffer);

	clSetKernelArg(stitchKernel,0,sizeof(cl_mem),&psi1.computeBuffer);
	clSetKernelArg(stitchKernel,2,sizeof(cl_mem),&v.computeBuffer);
	clSetKernelArg(stitchKernel,3,sizeof(cl_mem),&w.computeBuffer);

	#endif
}

Schroedinger::~Schroedinger()
{
	#ifndef USE_OPENCL
	if (propagator!=NULL)
		owner->RemoveTool(propagator);
	#endif
	#ifdef USE_OPENCL
	clReleaseKernel(applyNumerator);
	clReleaseKernel(applyDenominator);
	clReleaseKernel(chargeKernel);
	clReleaseKernel(currentKernel);
	clReleaseKernel(stitchKernel);
	#endif
}

void Schroedinger::VerifyInput()
{
	AtomicPhysics::VerifyInput();
	#ifndef USE_OPENCL
	if (propagator==NULL)
		propagator = (SchroedingerPropagator*)owner->CreateTool("TDSE",tw::tool_type::schroedingerPropagator);
	#endif
}

void Schroedinger::Initialize()
{
	AtomicPhysics::Initialize();

	const tw::bc::fld psiDefaultBC = fld::neumannWall;
	psi0.SetBoundaryConditions(tw::dom::xAxis,psiDefaultBC,psiDefaultBC);
	psi0.SetBoundaryConditions(tw::dom::yAxis,psiDefaultBC,psiDefaultBC);
	psi0.SetBoundaryConditions(tw::dom::zAxis,psiDefaultBC,psiDefaultBC);
	psi1.SetBoundaryConditions(tw::dom::xAxis,psiDefaultBC,psiDefaultBC);
	psi1.SetBoundaryConditions(tw::dom::yAxis,psiDefaultBC,psiDefaultBC);
	psi1.SetBoundaryConditions(tw::dom::zAxis,psiDefaultBC,psiDefaultBC);

	if (!owner->restarted)
	{
		// Solve for the lowest energy s-state on a spherical grid.
		// This is used only to print the numerical ground state energy level.
		const tw::Float maxR = owner->SphericalRadius(GlobalCorner(*owner)+GlobalPhysicalSize(*owner));
		const tw::Float dr = dx(*owner) * owner->ScaleFactor(1,tw::vec3(tw::small_pos,0.0,0.0));
		const tw::Float r = maxR>30.0 ? 30.0 : maxR;
		const tw::Int dim = MyCeil(r/dr);
		std::valarray<tw::Float> eigenvector(dim),phi_r(dim);
		for (tw::Int i=0;i<dim;i++)
			phi_r[i] = GetSphericalPotential((tw::Float(i)+0.5)*dr);
		groundStateEnergy = GetSphericalGroundState(eigenvector,phi_r,dr);
		(*owner->tw_out) << "Numerical ground state energy = " << groundStateEnergy << std::endl;

		for (tw::Int s=0;s<waveFunction.size();s++)
		{
			waveFunction[s].Initialize(qo::schroedinger,residualCharge,q0,m0);
			(*owner->tw_out) << "Wavefunction " << s << " : energy = " << waveFunction[s].energy << std::endl;
		}
		for (tw::Int s=0;s<refState.size();s++)
		{
			refState[s].Initialize(qo::schroedinger,residualCharge,q0,m0);
			(*owner->tw_out) << "Reference " << s << " : energy = " << refState[s].energy << std::endl;
		}

		#pragma omp parallel
		{
			for (auto cell : EntireCellRange(*this))
				for (tw::Int s=0;s<waveFunction.size();s++)
					psi1(cell) += waveFunction[s].Amplitude(owner->Pos(cell),0.0,0);
		}
		Normalize();
		psi0 = psi1;
		#pragma omp parallel
		{
			for (auto cell : EntireCellRange(*this))
				J4(cell,0) = norm(psi1(cell));
		}
	}

	FormPotentials(owner->elapsedTime);

	#ifdef USE_OPENCL
	psi1.SendToComputeBuffer();
	A4.SendToComputeBuffer();
	#endif
}

#ifdef USE_OPENCL
void Schroedinger::Update()
{
	Field mpi_packet;
	tw::Int i,j,k;
	tw::Float partitionFactor,relax,totalProbability;
	partitionFactor = 1.0/tw::Float(owner->Dimensionality());
	relax = owner->elapsedTime < timeRelaxingToGround ? 1.0 : 0.0;

	// Clear source vector
	J4.MADDComputeBuffer(0.0,0.0);

	// Update vector potential
	tw::vec4 A40 = GetA4AtOrigin();
	A4.FillComputeBufferVec4(Element(1,3),A40);

	// Set up compute arguments
	// argument 3 is the strip argument
	clSetKernelArg(applyNumerator,4,sizeof(tw::Float),&partitionFactor);
	clSetKernelArg(applyNumerator,5,sizeof(tw::Float),&relax);
	clSetKernelArg(applyDenominator,4,sizeof(tw::Float),&partitionFactor); // 3 is the strip
	clSetKernelArg(applyDenominator,5,sizeof(tw::Float),&relax);

	// Start Charge Deposition
	owner->CellUpdateProtocol(chargeKernel);

	// Loop over Dimensions Updating Wavefunction and Currents
	for (i=1;i<=3;i++)
	{
		if (owner->localCells[i]>1)
		{
			CopyComputeBuffer(psi0,psi1);
			j = i > 2 ? i-2 : i+1; // 2 , 3 , 1
			k = i > 1 ? i-1 : i+2; // 3 , 1 , 2
			DiscreteSpace mpi_layout;
			mpi_layout.Resize(owner->localCells[j],owner->localCells[k],1,corner,size);
			mpi_packet.Initialize(16,mpi_layout,owner);
			mpi_packet.InitializeComputeBuffer();
			clSetKernelArg(applyDenominator,9,sizeof(cl_mem),&mpi_packet.computeBuffer);
			owner->StripUpdateProtocol(applyNumerator,i,3);
			owner->StripUpdateProtocol(applyDenominator,i,3);

			// Global Integration
			mpi_packet.ReceiveFromComputeBuffer();
			ComputeAlphasAndBetas<tw::Complex>(&owner->strip[i],owner->localCells2[j]*owner->localCells2[k],(tw::Complex*)&mpi_packet(0,0,0,0));
			mpi_packet.SendToComputeBuffer();
			clSetKernelArg(stitchKernel,1,sizeof(cl_mem),&owner->stripBuffer[i]);
			clSetKernelArg(stitchKernel,4,sizeof(cl_mem),&mpi_packet.computeBuffer);
			owner->LocalUpdateProtocol(stitchKernel);
			psi1.UpdateGhostCellsInComputeBuffer();

			if (owner->elapsedTime > timeRelaxingToGround)
				owner->StripUpdateProtocol(currentKernel,i,5);
		}
	}

	// Finish Charge Deposition
	owner->CellUpdateProtocol(chargeKernel);

	// Normalize After Relaxation Step
	if (owner->elapsedTime < timeRelaxingToGround)
	{
		CopyComputeBuffer(scratch,psi1);
		scratch.DestructiveComplexMod2ComputeBuffer();
		scratch.ZeroGhostCellsInComputeBuffer();
		scratch.WeightComputeBufferByVolume(*owner,0.0);
		totalProbability = scratch.DestructiveSumComputeBuffer();
		owner->strip[0].AllSum(&totalProbability,&totalProbability,sizeof(tw::Float),0);
		psi1.MADDComputeBuffer(1.0/sqrt(totalProbability),0.0);
	}

	// Get J4 for Bohmian Trajectories
	J4.UpdateGhostCellsInComputeBuffer();
	J4.ReceiveFromComputeBuffer();
}
#else
void Schroedinger::Update()
{
	tw::Complex dtc = owner->elapsedTime < timeRelaxingToGround ? -ii*dt : dt;
	FormPotentials(owner->elapsedTime);
	J4 = 0.0;

	propagator->DepositCurrent(tw::dom::tAxis,psi0,psi1,A4,J4,dtc);

	psi0 = psi1;
	propagator->ApplyNumerator(tw::dom::xAxis,psi1,A4,keepA2Term,dtc);
	propagator->ApplyDenominator(tw::dom::xAxis,psi1,A4,keepA2Term,dtc);
	propagator->DepositCurrent(tw::dom::xAxis,psi0,psi1,A4,J4,dtc);

	psi0 = psi1;
	propagator->ApplyNumerator(tw::dom::yAxis,psi1,A4,keepA2Term,dtc);
	propagator->ApplyDenominator(tw::dom::yAxis,psi1,A4,keepA2Term,dtc);
	propagator->DepositCurrent(tw::dom::yAxis,psi0,psi1,A4,J4,dtc);

	psi0 = psi1;
	propagator->ApplyNumerator(tw::dom::zAxis,psi1,A4,keepA2Term,dtc);
	propagator->ApplyDenominator(tw::dom::zAxis,psi1,A4,keepA2Term,dtc);
	propagator->DepositCurrent(tw::dom::zAxis,psi0,psi1,A4,J4,dtc);

	propagator->DepositCurrent(tw::dom::tAxis,psi0,psi1,A4,J4,dtc);

	J4.CopyFromNeighbors();
	J4.ApplyBoundaryCondition();
	if (owner->elapsedTime < timeRelaxingToGround)
		Normalize();
}
// void Schroedinger::Update()
// {
// 	// This advances the wavefunction conservatively, but uses the naive
// 	// centered differenced current deposition
// 	tw::Complex dtc = owner->elapsedTime < timeRelaxingToGround ? -ii*dt : dt;
// 	FormPotentials(owner->elapsedTime);
//
// 	psi0 = psi1;
// 	propagator->ApplyNumerator(tw::dom::xAxis,psi1,A4,keepA2Term,dtc);
// 	propagator->ApplyDenominator(tw::dom::xAxis,psi1,A4,keepA2Term,dtc);
// 	propagator->ApplyNumerator(tw::dom::yAxis,psi1,A4,keepA2Term,dtc);
// 	propagator->ApplyDenominator(tw::dom::yAxis,psi1,A4,keepA2Term,dtc);
// 	propagator->ApplyNumerator(tw::dom::zAxis,psi1,A4,keepA2Term,dtc);
// 	propagator->ApplyDenominator(tw::dom::zAxis,psi1,A4,keepA2Term,dtc);
// 	if (owner->elapsedTime < timeRelaxingToGround)
// 		Normalize();
// 	else
// 		UpdateJ4();
// }
#endif

void Schroedinger::UpdateJ4()
{
	// Here is the naive centered differencing to determine the non-relativistic probability current.
	// The results are usually pathological.  The conservative scheme should be used instead.
	tw::Int i,j,k;
	tw::Complex psiNow,psi_x,psi_y,psi_z;
	if (owner->elapsedTime < timeRelaxingToGround)
		return;

	for (k=1;k<=dim[3];k++)
		for (j=1;j<=dim[2];j++)
			for (i=1;i<=dim[1];i++)
			{
				psiNow = psi1(i,j,k);
				psi_x = (psi1(i+1,j,k) - psi1(i-1,j,k))/owner->dL(i,j,k,1);
				psi_y = (psi1(i,j+1,k) - psi1(i,j-1,k))/owner->dL(i,j,k,2);
				psi_z = (psi1(i,j,k+1) - psi1(i,j,k-1))/owner->dL(i,j,k,3);

				J4(i,j,k,0) = norm(psiNow);
				J4(i,j,k,1) = -real((half*ii/m0) * (conj(psiNow)*psi_x - conj(psi_x)*psiNow)) - q0*A4(i,j,k,1)*norm(psiNow)/m0;
				J4(i,j,k,2) = -real((half*ii/m0) * (conj(psiNow)*psi_y - conj(psi_y)*psiNow)) - q0*A4(i,j,k,2)*norm(psiNow)/m0;
				J4(i,j,k,3) = -real((half*ii/m0) * (conj(psiNow)*psi_z - conj(psi_z)*psiNow)) - q0*A4(i,j,k,3)*norm(psiNow)/m0;
			}
	J4.CopyFromNeighbors();
	J4.ApplyBoundaryCondition();
}

void Schroedinger::Normalize()
{
	tw::Float totalProbability = 0.0;
	for (auto cell : InteriorCellRange(*this))
		totalProbability += norm(psi1(cell)) * owner->dS(cell,0);
	owner->strip[0].AllSum(&totalProbability,&totalProbability,sizeof(tw::Float),0);
	psi1 *= 1.0/sqrt(totalProbability);
	psi1.CopyFromNeighbors();
	psi1.ApplyBoundaryCondition();
}

void Schroedinger::EnergyHeadings(std::ofstream& outFile)
{
	outFile << "TotalProb Ex Ey Ez Ax Ay Az Px Py Pz ";

	tw::Int i;
	for (i=0;i<refState.size();i++)
		outFile << "real<ref" << i << "|psi> imag<ref" << i << "|psi> ";
}

void Schroedinger::EnergyColumns(std::vector<tw::Float>& cols,std::vector<bool>& avg,const Region& theRgn)
{
	tw::Int i,j,k,s;
	tw::Int x0,x1,y0,y1,z0,z1;
	tw::Float dV;
	tw::vec3 pos,ENow,ANow;
	std::valarray<tw::Complex> overlap(refState.size());
	tw::Complex psiNow;
	tw::Float totalProb = 0.0;
	tw::vec3 dipoleMoment = 0.0;
	ENow = 0.0;
	ANow = 0.0;
	overlap = tw::Complex(0,0);

	for (s=0;s<owner->wave.size();s++)
	{
		ANow += owner->wave[s]->VectorPotential(owner->elapsedTime,tw::vec3(0,0,0));
		ENow -= dti*owner->wave[s]->VectorPotential(owner->elapsedTime+0.5*dt,tw::vec3(0,0,0));
		ENow += dti*owner->wave[s]->VectorPotential(owner->elapsedTime-0.5*dt,tw::vec3(0,0,0));
	}

	theRgn.GetLocalCellBounds(&x0,&x1,&y0,&y1,&z0,&z1);
	for (k=z0;k<=z1;k++)
		for (j=y0;j<=y1;j++)
			for (i=x0;i<=x1;i++)
			{
				pos = owner->Pos(i,j,k);
				if (theRgn.Inside(pos,*owner))
				{
					dV = owner->dS(i,j,k,0);
					psiNow = psi1(i,j,k);
					totalProb += norm(psiNow)*dV;
					dipoleMoment += norm(psiNow)*pos*dV;
					// Gauge transformation (remove uniform A)
					psiNow *= std::exp(-ii*q0*(ANow^pos));
					for (s=0;s<refState.size();s++)
						overlap[s] += conj(refState[s].Amplitude(pos,owner->elapsedTime,0))*psiNow*dV;
				}
			}

	cols.push_back(totalProb); avg.push_back(false);
	cols.push_back(ENow.x); avg.push_back(true);
	cols.push_back(ENow.y); avg.push_back(true);
	cols.push_back(ENow.z); avg.push_back(true);
	cols.push_back(ANow.x); avg.push_back(true);
	cols.push_back(ANow.y); avg.push_back(true);
	cols.push_back(ANow.z); avg.push_back(true);
	cols.push_back(dipoleMoment.x); avg.push_back(false);
	cols.push_back(dipoleMoment.y); avg.push_back(false);
	cols.push_back(dipoleMoment.z); avg.push_back(false);
	for (s=0;s<refState.size();s++)
	{
		cols.push_back(real(overlap[s])); avg.push_back(false);
		cols.push_back(imag(overlap[s])); avg.push_back(false);
	}
}

void Schroedinger::BoxDiagnosticHeader(GridDataDescriptor* box)
{
	owner->WriteBoxDataHeader("psi_r",box);
	owner->WriteBoxDataHeader("psi_i",box);
	owner->WriteBoxDataHeader("rho",box);
	owner->WriteBoxDataHeader("Jx",box);
	owner->WriteBoxDataHeader("Jy",box);
	owner->WriteBoxDataHeader("Jz",box);
	owner->WriteBoxDataHeader("phi",box);
	owner->WriteBoxDataHeader("Ax",box);
	owner->WriteBoxDataHeader("Ay",box);
	owner->WriteBoxDataHeader("Az",box);
	owner->WriteBoxDataHeader("divJ",box);
}

void Schroedinger::BoxDiagnose(GridDataDescriptor* box)
{
	owner->WriteBoxData("psi_r",box,&psi1(0,0,0,0),psi1.Stride());
	owner->WriteBoxData("psi_i",box,&psi1(0,0,0,1),psi1.Stride());

	owner->WriteBoxData("rho",box,&J4(0,0,0,0),J4.Stride());
	owner->WriteBoxData("Jx",box,&J4(0,0,0,1),J4.Stride());
	owner->WriteBoxData("Jy",box,&J4(0,0,0,2),J4.Stride());
	owner->WriteBoxData("Jz",box,&J4(0,0,0,3),J4.Stride());

	owner->WriteBoxData("phi",box,&A4(0,0,0,0),A4.Stride());
	owner->WriteBoxData("Ax",box,&A4(0,0,0,1),A4.Stride());
	owner->WriteBoxData("Ay",box,&A4(0,0,0,2),A4.Stride());
	owner->WriteBoxData("Az",box,&A4(0,0,0,3),A4.Stride());

	tw::Int i,j,k;
	ScalarField divJ;
	divJ.Initialize(*this,owner);
	for (k=1;k<=dim[3];k++)
		for (j=1;j<=dim[2];j++)
			for (i=1;i<=dim[1];i++)
				divJ(i,j,k) = div<1,2,3>(J4,i,j,k,*owner);
	owner->WriteBoxData("divJ",box,&divJ(0,0,0),divJ.Stride());
}


void Schroedinger::PointDiagnosticHeader(std::ofstream& outFile)
{
	outFile << "psi_r psi_i Ax Ay Az ";
}

void Schroedinger::PointDiagnose(std::ofstream& outFile,const weights_3D& w)
{
	tw::Int s;
	tw::vec3 ANow = 0.0;
	tw::Complex psiNow;
	psi1.Interpolate(&psiNow,w);
	for (s=0;s<owner->wave.size();s++)
		ANow += owner->wave[s]->VectorPotential(owner->elapsedTime,tw::vec3(0,0,0));
	outFile << real(psiNow) << " " << imag(psiNow) << " " << ANow.x << " " << ANow.y << " " << ANow.z << " ";
}

void Schroedinger::StartDiagnostics()
{
	#ifdef USE_OPENCL
	psi1.ReceiveFromComputeBuffer();
	A4.ReceiveFromComputeBuffer();
	#endif
}


/////////////////////////////////////////
//                                     //
// Dual TDSE Coupled by Pauli Matrices //
//                                     //
/////////////////////////////////////////

Pauli::Pauli(const std::string& name,Simulation* sim):AtomicPhysics(name,sim)
{
	throw tw::FatalError("Pauli module is not supported in this version of TW.");
	typeCode = tw::module_type::pauli;
	#ifndef USE_OPENCL
	propagator = NULL;
	#endif
	psi0.Initialize(*this,owner);
	psi1.Initialize(*this,owner);
	chi0.Initialize(*this,owner);
	chi1.Initialize(*this,owner);
	v.Initialize(*this,owner);
	w.Initialize(*this,owner);
	scratch.Initialize(*this,owner);
	#ifdef USE_OPENCL
	scratch.InitializeComputeBuffer();
	v.InitializeComputeBuffer();
	w.InitializeComputeBuffer();
	psi0.InitializeComputeBuffer();
	psi1.InitializeComputeBuffer();
	chi0.InitializeComputeBuffer();
	chi1.InitializeComputeBuffer();
	J4.InitializeComputeBuffer();
	#endif
}

Pauli::~Pauli()
{
	#ifndef USE_OPENCL
	if (propagator!=NULL)
		owner->RemoveTool(propagator);
	#endif
}

void Pauli::VerifyInput()
{
	AtomicPhysics::VerifyInput();
	#ifndef USE_OPENCL
	if (propagator==NULL)
		propagator = (SchroedingerPropagator*)owner->CreateTool("TDSE",tw::tool_type::schroedingerPropagator);
	#endif
}

void Pauli::Initialize()
{
	AtomicPhysics::Initialize();
}

#ifdef USE_OPENCL
void Pauli::Update()
{
	cl_int err;
	cl_kernel spinKernel = clCreateKernel(program,"UpdateSpin",&err);

	// Spin
	clSetKernelArg(spinKernel,0,sizeof(cl_mem),&psi1.computeBuffer);
	clSetKernelArg(spinKernel,1,sizeof(cl_mem),&chi1.computeBuffer);
	clSetKernelArg(spinKernel,2,sizeof(cl_mem),&A4.computeBuffer);
	clSetKernelArg(spinKernel,3,sizeof(cl_mem),&owner->metricsBuffer);
	owner->LocalUpdateProtocol(spinKernel);

	// ADD: something to handle normalization after a relaxation step

	clReleaseKernel(spinKernel);
}
#else
void Pauli::Update()
{
	// KNOWN PROBLEM: we have to add spin term to current density

	tw::Complex dtc = owner->elapsedTime < timeRelaxingToGround ? -ii*dt : dt;
	FormPotentials(owner->elapsedTime);
	J4 = 0.0;

	propagator->DepositCurrent(tw::dom::tAxis,psi0,psi1,A4,J4,dtc);
	propagator->DepositCurrent(tw::dom::tAxis,chi0,chi1,A4,J4,dtc);

	psi0 = psi1;
	chi0 = chi1;
	propagator->ApplyNumerator(tw::dom::xAxis,psi1,A4,keepA2Term,dtc);
	propagator->ApplyDenominator(tw::dom::xAxis,psi1,A4,keepA2Term,dtc);
	propagator->DepositCurrent(tw::dom::xAxis,psi0,psi1,A4,J4,dtc);
	propagator->ApplyNumerator(tw::dom::xAxis,chi1,A4,keepA2Term,dtc);
	propagator->ApplyDenominator(tw::dom::xAxis,chi1,A4,keepA2Term,dtc);
	propagator->DepositCurrent(tw::dom::xAxis,chi0,chi1,A4,J4,dtc);

	psi0 = psi1;
	chi0 = chi1;
	propagator->ApplyNumerator(tw::dom::yAxis,psi1,A4,keepA2Term,dtc);
	propagator->ApplyDenominator(tw::dom::yAxis,psi1,A4,keepA2Term,dtc);
	propagator->DepositCurrent(tw::dom::yAxis,psi0,psi1,A4,J4,dtc);
	propagator->ApplyNumerator(tw::dom::yAxis,chi1,A4,keepA2Term,dtc);
	propagator->ApplyDenominator(tw::dom::yAxis,chi1,A4,keepA2Term,dtc);
	propagator->DepositCurrent(tw::dom::yAxis,chi0,chi1,A4,J4,dtc);

	psi0 = psi1;
	chi0 = chi1;
	propagator->ApplyNumerator(tw::dom::zAxis,psi1,A4,keepA2Term,dtc);
	propagator->ApplyDenominator(tw::dom::zAxis,psi1,A4,keepA2Term,dtc);
	propagator->DepositCurrent(tw::dom::zAxis,psi0,psi1,A4,J4,dtc);
	propagator->ApplyNumerator(tw::dom::zAxis,chi1,A4,keepA2Term,dtc);
	propagator->ApplyDenominator(tw::dom::zAxis,chi1,A4,keepA2Term,dtc);
	propagator->DepositCurrent(tw::dom::zAxis,chi0,chi1,A4,J4,dtc);

	propagator->DepositCurrent(tw::dom::tAxis,psi0,psi1,A4,J4,dtc);
	propagator->DepositCurrent(tw::dom::tAxis,chi0,chi1,A4,J4,dtc);

   	propagator->UpdateSpin(psi1,chi1,A4,alpha*dt);

	if (owner->elapsedTime < timeRelaxingToGround)
		Normalize();
}
#endif

void Pauli::Normalize()
{
	tw::Float totalProbability = 0.0;
	for (auto cell : InteriorCellRange(*this))
		totalProbability += (norm(psi1(cell))+norm(chi1(cell))) * owner->dS(cell,0);
	owner->strip[0].AllSum(&totalProbability,&totalProbability,sizeof(tw::Float),0);
	psi1 *= 1.0/sqrt(totalProbability);
	chi1 *= 1.0/sqrt(totalProbability);
	psi1.CopyFromNeighbors();
	psi1.ApplyBoundaryCondition();
	chi1.CopyFromNeighbors();
	chi1.ApplyBoundaryCondition();
}

void Pauli::EnergyHeadings(std::ofstream& outFile)
{
	outFile << "TotalProb Ex Ey Ez Ax Ay Az Px Py Pz Sz ";

	tw::Int i;
	for (i=0;i<refState.size();i++)
		outFile << "real<ref" << i << "|psi> imag<ref" << i << "|psi> ";
}

void Pauli::EnergyColumns(std::vector<tw::Float>& cols,std::vector<bool>& avg,const Region& theRgn)
{
	tw::Int i,j,k,s;
	tw::Int x0,x1,y0,y1,z0,z1;
	tw::Float dV;
	tw::vec3 pos,ENow,ANow;
	std::valarray<tw::Complex> overlap(refState.size());

	tw::Float totalProb = 0.0, Sz = 0.0;
	tw::vec3 dipoleMoment = 0.0;
	ENow = 0.0;
	ANow = 0.0;
	overlap = tw::Complex(0,0);

	theRgn.GetLocalCellBounds(&x0,&x1,&y0,&y1,&z0,&z1);
	for (k=z0;k<=z1;k++)
		for (j=y0;j<=y1;j++)
			for (i=x0;i<=x1;i++)
			{
				pos = owner->Pos(i,j,k);
				if (theRgn.Inside(pos,*owner))
				{
					dV = owner->dS(i,j,k,0);
					totalProb += (norm(psi1(i,j,k)) + norm(chi1(i,j,k)))*dV;
					dipoleMoment += (norm(psi1(i,j,k)) + norm(chi1(i,j,k)))*pos*dV;
					Sz += (norm(psi1(i,j,k)) - norm(chi1(i,j,k)))*dV;
					for (s=0;s<refState.size();s++)
						overlap[s] += 0.0; // need to revisit
				}
			}

	for (s=0;s<owner->wave.size();s++)
	{
		ANow += owner->wave[s]->VectorPotential(owner->elapsedTime,tw::vec3(0,0,0));
		ENow -= dti*owner->wave[s]->VectorPotential(owner->elapsedTime+0.5*dt,tw::vec3(0,0,0));
		ENow += dti*owner->wave[s]->VectorPotential(owner->elapsedTime-0.5*dt,tw::vec3(0,0,0));
	}

	cols.push_back(totalProb); avg.push_back(false);
	cols.push_back(ENow.x); avg.push_back(true);
	cols.push_back(ENow.y); avg.push_back(true);
	cols.push_back(ENow.z); avg.push_back(true);
	cols.push_back(ANow.x); avg.push_back(true);
	cols.push_back(ANow.y); avg.push_back(true);
	cols.push_back(ANow.z); avg.push_back(true);
	cols.push_back(dipoleMoment.x); avg.push_back(false);
	cols.push_back(dipoleMoment.y); avg.push_back(false);
	cols.push_back(dipoleMoment.z); avg.push_back(false);
	cols.push_back(Sz); avg.push_back(false);
	for (s=0;s<refState.size();s++)
	{
		cols.push_back(real(overlap[s])); avg.push_back(false);
		cols.push_back(imag(overlap[s])); avg.push_back(false);
	}
}

void Pauli::BoxDiagnosticHeader(GridDataDescriptor* box)
{
	owner->WriteBoxDataHeader("psi_r",box);
	owner->WriteBoxDataHeader("psi_i",box);
	owner->WriteBoxDataHeader("chi_r",box);
	owner->WriteBoxDataHeader("chi_i",box);
	owner->WriteBoxDataHeader("rho",box);
	owner->WriteBoxDataHeader("Jx",box);
	owner->WriteBoxDataHeader("Jy",box);
	owner->WriteBoxDataHeader("Jz",box);
	owner->WriteBoxDataHeader("phi",box);
	owner->WriteBoxDataHeader("Ax",box);
	owner->WriteBoxDataHeader("Ay",box);
	owner->WriteBoxDataHeader("Az",box);
	owner->WriteBoxDataHeader("Sz",box);
}

void Pauli::BoxDiagnose(GridDataDescriptor* box)
{
	owner->WriteBoxData("psi_r",box,&psi1(0,0,0,0),psi1.Stride());
	owner->WriteBoxData("psi_i",box,&psi1(0,0,0,1),psi1.Stride());
	owner->WriteBoxData("chi_r",box,&chi1(0,0,0,0),chi1.Stride());
	owner->WriteBoxData("chi_i",box,&chi1(0,0,0,1),chi1.Stride());

	owner->WriteBoxData("rho",box,&J4(0,0,0,0),J4.Stride());
	owner->WriteBoxData("Jx",box,&J4(0,0,0,1),J4.Stride());
	owner->WriteBoxData("Jy",box,&J4(0,0,0,2),J4.Stride());
	owner->WriteBoxData("Jz",box,&J4(0,0,0,3),J4.Stride());

	owner->WriteBoxData("phi",box,&A4(0,0,0,0),A4.Stride());
	owner->WriteBoxData("Ax",box,&A4(0,0,0,1),A4.Stride());
	owner->WriteBoxData("Ay",box,&A4(0,0,0,2),A4.Stride());
	owner->WriteBoxData("Az",box,&A4(0,0,0,3),A4.Stride());

	ScalarField Sz;
	Sz.Initialize(*this,owner);
	for (auto cell : InteriorCellRange(*this))
		Sz(cell) = norm(psi1(cell))-norm(chi1(cell));
	owner->WriteBoxData("Sz",box,&Sz(0,0,0),Sz.Stride());
}

void Pauli::StartDiagnostics()
{
	#ifdef USE_OPENCL
	psi0.ReceiveFromComputeBuffer();
	psi1.ReceiveFromComputeBuffer();
	chi0.ReceiveFromComputeBuffer();
	chi1.ReceiveFromComputeBuffer();
	A4.ReceiveFromComputeBuffer();
	J4.ReceiveFromComputeBuffer();
	#endif
}


///////////////////////////////////////
//                                   //
//  RELATIVISTIC PARTICLE (KG EQN)   //
//                                   //
///////////////////////////////////////


KleinGordon::KleinGordon(const std::string& name,Simulation* sim) : AtomicPhysics(name,sim)
{
	// Wavefunction is in Hamiltonian 2-component representation.
	// This is the preliminary Hamiltonian form from Feshbach-Villars, Eq. 2.12.
	// Unlike the symmetric form, this allows a leap-frog scheme to be employed.
	// In this representation, component 0 is the usual scalar wavefunction.
	// Component 1 is psi1 = (id/dt - q*phi)*psi0/m
	// The symmetric form is obtained from (psi0+psi1)/sqrt(2) , (psi0-psi1)/sqrt(2)
	typeCode = tw::module_type::kleinGordon;
	m0 = 1.0;
	q0 = -sqrt(alpha);
	residualCharge = sqrt(alpha);
	dipoleApproximation = false;
	psi_r.Initialize(2,*this,owner,tw::dom::xAxis);
	psi_i.Initialize(2,*this,owner,tw::dom::xAxis);

	#ifdef USE_OPENCL
	cl_int err;
	updatePsi = clCreateKernel(program,"KGAdvance1",&err);
	updateChi = clCreateKernel(program,"KGAdvance2",&err);

	psi_r.InitializeComputeBuffer();
	psi_i.InitializeComputeBuffer();

	clSetKernelArg(updatePsi,0,sizeof(cl_mem),&psi_r.computeBuffer);
	clSetKernelArg(updatePsi,1,sizeof(cl_mem),&psi_i.computeBuffer);
	clSetKernelArg(updatePsi,2,sizeof(cl_mem),&A4.computeBuffer);
	clSetKernelArg(updatePsi,3,sizeof(cl_mem),&owner->metricsBuffer);

	clSetKernelArg(updateChi,0,sizeof(cl_mem),&psi_r.computeBuffer);
	clSetKernelArg(updateChi,1,sizeof(cl_mem),&psi_i.computeBuffer);
	clSetKernelArg(updateChi,2,sizeof(cl_mem),&A4.computeBuffer);
	clSetKernelArg(updateChi,3,sizeof(cl_mem),&owner->metricsBuffer);
	#endif
}

KleinGordon::~KleinGordon()
{
	#ifdef USE_OPENCL
	clReleaseKernel(updatePsi);
	clReleaseKernel(updateChi);
	#endif
}

void KleinGordon::Initialize()
{
	AtomicPhysics::Initialize();

	if (!owner->restarted)
	{
		for (tw::Int s=0;s<waveFunction.size();s++)
		{
			waveFunction[s].Initialize(qo::klein_gordon,residualCharge,q0,m0);
			(*owner->tw_out) << "Wavefunction " << s << " : energy = " << waveFunction[s].energy << std::endl;
		}
		for (tw::Int s=0;s<refState.size();s++)
		{
			refState[s].Initialize(qo::klein_gordon,residualCharge,q0,m0);
			(*owner->tw_out) << "Reference " << s << " : energy = " << refState[s].energy << std::endl;
		}

		for (auto cell : InteriorCellRange(*this))
		{
			const tw::vec3 pos = owner->Pos(cell);
			for (tw::Int s=0;s<waveFunction.size();s++)
			{
				psi_r(cell,0) += real(waveFunction[s].Amplitude(pos,0.0,0));
				psi_i(cell,0) += imag(waveFunction[s].Amplitude(pos,0.0,0));
				psi_r(cell,1) += real(waveFunction[s].Amplitude(pos,dth,1));
				psi_i(cell,1) += imag(waveFunction[s].Amplitude(pos,dth,1));
			}
		}

		psi_r.CopyFromNeighbors();
		psi_i.CopyFromNeighbors();
		FormPotentials(owner->elapsedTime);
		Normalize();
		UpdateJ4();
	}

	#ifdef USE_OPENCL

	psi_r.SendToComputeBuffer();
	psi_i.SendToComputeBuffer();
	Ao4.SendToComputeBuffer();
	A4.SendToComputeBuffer();
	clSetKernelArg(updatePsi,4,sizeof(m0),&m0);
	clSetKernelArg(updatePsi,5,sizeof(q0),&q0);
	clSetKernelArg(updateChi,4,sizeof(m0),&m0);
	clSetKernelArg(updateChi,5,sizeof(q0),&q0);

	#endif
}

tw::Float KleinGordon::ComputeRho(const tw::cell& cell)
{
	return q0*(norm(FV(cell,1.0)) - norm(FV(cell,-1.0)));
}

void KleinGordon::UpdateJ4()
{
	#pragma omp parallel
	{
		for (auto v : VectorStripRange<1>(*this,false))
		{
			for (tw::Int i=1;i<=dim[1];i++)
			{
				J4(v,i,0) = q0*(norm(FV(v,i,1.0)) - norm(FV(v,i,-1.0)));
				J4(v,i,1) = 0.0;
				J4(v,i,2) = 0.0;
				J4(v,i,3) = 0.0;
			}
		}
	}
	J4.CopyFromNeighbors();
	J4.ApplyBoundaryCondition();
}

void KleinGordon::Normalize()
{
	tw::Float totalCharge = 0.0;
	for (auto cell : InteriorCellRange(*this))
		totalCharge += ComputeRho(cell) * owner->dS(cell,0);
	owner->strip[0].AllSum(&totalCharge,&totalCharge,sizeof(tw::Float),0);
	psi_r *= sqrt(fabs(q0/totalCharge));
	psi_i *= sqrt(fabs(q0/totalCharge));
	psi_r.CopyFromNeighbors();
	psi_r.ApplyBoundaryCondition();
	psi_i.CopyFromNeighbors();
	psi_i.ApplyBoundaryCondition();
}

#ifdef USE_OPENCL
void KleinGordon::Update()
{
	// Update wavefunction
	owner->LocalUpdateProtocol(updatePsi);
	// Boundary conditions and communications
	psi_r.UpdateGhostCellsInComputeBuffer(Element(0));
	psi_i.UpdateGhostCellsInComputeBuffer(Element(0));

	photonPropagator->Advance(A4,Ao4,J4,0.0,dt);
	FormGhostCellPotentials(owner->elapsedTime+dt);
	A4.SendGhostCellsToComputeBuffer();
	photonPropagator->MidstepEstimate(A4,Ao4);

	// Update auxiliary wavefunction
	owner->LocalUpdateProtocol(updateChi);
	// Boundary conditions and communications
	psi_r.UpdateGhostCellsInComputeBuffer(Element(1));
	psi_i.UpdateGhostCellsInComputeBuffer(Element(1));

	photonPropagator->UndoMidstepEstimate(A4,Ao4);
}
#else
void KleinGordon::Update()
{
	// Solve in the Hamiltonian leap-frog representation
	// dpsi/dt = -i*q*phi*psi - i*m*chi
	// dchi/dt = -i*q*phi*chi - i*m*psi + i*(Dk)^2(psi)/m
	// Leapfrog Scheme:
	// At start we know psi(n-1/2),chi(n)
	// Advance psi to n+1/2 using chi(n),A(n)
	// Advance chi to n+1 using psi(n+1/2),A(n+1/2)

	static const tw::Int AB = tw::vec_align_bytes;
	#pragma omp parallel
	{
		alignas(AB) tw::Float Ur[dim[1]],Ui[dim[1]],Dr[dim[1]],Di[dim[1]];
		// Update psi
		for (auto v : VectorStripRange<1>(*this,false))
		{
			#pragma omp simd aligned(Ur,Ui:AB)
			for (tw::Int i=1;i<=dim[1];i++)
			{
				// unitary operator of time translation for diagonal part of Hamiltonian
				const tw::Float dq = dt*q0*A4(v,i,0);
				Ur[i-1] = cos(dq);
				Ui[i-1] = -sin(dq);
			}
			#pragma omp simd aligned(Dr,Di:AB)
			for (tw::Int i=1;i<=dim[1];i++)
			{
				Dr[i-1] = m0*psi_i(v,i,1);
				Di[i-1] = -m0*psi_r(v,i,1);
			}
			#pragma omp simd aligned(Ur,Ui,Dr,Di:AB)
			for (tw::Int i=1;i<=dim[1];i++)
			{
				psi_r(v,i,0) += dth*Dr[i-1];
				psi_i(v,i,0) += dth*Di[i-1];
				complex_multiply_assign(psi_r(v,i,0),psi_i(v,i,0),Ur[i-1],Ui[i-1]);
				psi_r(v,i,0) += dth*Dr[i-1];
				psi_i(v,i,0) += dth*Di[i-1];
			}
		}
	}

	psi_r.CopyFromNeighbors(Element(0));
	psi_r.ApplyBoundaryCondition(Element(0));
	psi_i.CopyFromNeighbors(Element(0));
	psi_i.ApplyBoundaryCondition(Element(0));

	photonPropagator->Advance(A4,Ao4,J4,0.0,dt);
	FormGhostCellPotentials(owner->elapsedTime+dt);
	photonPropagator->MidstepEstimate(A4,Ao4);

	#pragma omp parallel
	{
		alignas(AB) tw::Float Ur[dim[1]],Ui[dim[1]],Dr[dim[1]],Di[dim[1]];
		// Update chi
		for (auto v : VectorStripRange<1>(*this,false))
		{
			#pragma omp simd aligned(Ur,Ui:AB)
			for (tw::Int i=1;i<=dim[1];i++)
			{
				// unitary operator of time translation for diagonal part of Hamiltonian
				const tw::Float dq = dt*q0*A4(v,i,0);
				Ur[i-1] = cos(dq);
				Ui[i-1] = -sin(dq);
			}
			#pragma omp simd aligned(Dr,Di:AB)
			for (tw::Int i=1;i<=dim[1];i++)
			{
				// Evaluate i(Dk)^2 assuming div(A)=0
				// This adds i*del^2(psi)
				Dr[i-1] = -(psi_i.d2(v,i,0,1) + psi_i.d2(v,i,0,2) + psi_i.d2(v,i,0,3));
				Di[i-1] = (psi_r.d2(v,i,0,1) + psi_r.d2(v,i,0,2) + psi_r.d2(v,i,0,3));
				// This adds 2q*div(A*psi)
				Dr[i-1] += 2*q0*freq.x*(psi_r.sfwd(v,i,0,1)*A4.sfwd(v,i,1,1) - psi_r.sbak(v,i,0,1)*A4.sbak(v,i,1,1));
				Dr[i-1] += 2*q0*freq.y*(psi_r.sfwd(v,i,0,2)*A4.sfwd(v,i,2,2) - psi_r.sbak(v,i,0,2)*A4.sbak(v,i,2,2));
				Dr[i-1] += 2*q0*freq.z*(psi_r.sfwd(v,i,0,3)*A4.sfwd(v,i,3,3) - psi_r.sbak(v,i,0,3)*A4.sbak(v,i,3,3));
				Di[i-1] += 2*q0*freq.x*(psi_i.sfwd(v,i,0,1)*A4.sfwd(v,i,1,1) - psi_i.sbak(v,i,0,1)*A4.sbak(v,i,1,1));
				Di[i-1] += 2*q0*freq.y*(psi_i.sfwd(v,i,0,2)*A4.sfwd(v,i,2,2) - psi_i.sbak(v,i,0,2)*A4.sbak(v,i,2,2));
				Di[i-1] += 2*q0*freq.z*(psi_i.sfwd(v,i,0,3)*A4.sfwd(v,i,3,3) - psi_i.sbak(v,i,0,3)*A4.sbak(v,i,3,3));
				// This adds -i*q^2*A^2*psi
				Dr[i-1] += q0*q0*(sqr(A4(v,i,1))+sqr(A4(v,i,2))+sqr(A4(v,i,3)))*psi_i(v,i,0);
				Di[i-1] -= q0*q0*(sqr(A4(v,i,1))+sqr(A4(v,i,2))+sqr(A4(v,i,3)))*psi_r(v,i,0);
				// Finish by incorporating mass
				Dr[i-1] = Dr[i-1]/m0 + m0*psi_i(v,i,0);
				Di[i-1] = Di[i-1]/m0 - m0*psi_r(v,i,0);
			}
			#pragma omp simd aligned(Ur,Ui,Dr,Di:AB)
			for (tw::Int i=1;i<=dim[1];i++)
			{
				// Perform the integration
				psi_r(v,i,1) += dth*Dr[i-1];
				psi_i(v,i,1) += dth*Di[i-1];
				complex_multiply_assign(psi_r(v,i,1),psi_i(v,i,1),Ur[i-1],Ui[i-1]);
				psi_r(v,i,1) += dth*Dr[i-1];
				psi_i(v,i,1) += dth*Di[i-1];
			}
		}
	}

	psi_r.CopyFromNeighbors(Element(1));
	psi_r.ApplyBoundaryCondition(Element(1));
	psi_i.CopyFromNeighbors(Element(1));
	psi_i.ApplyBoundaryCondition(Element(1));

	photonPropagator->UndoMidstepEstimate(A4,Ao4);
}
#endif

void KleinGordon::EnergyHeadings(std::ofstream& outFile)
{
	outFile << "TotalCharge Ex Ey Ez Ax Ay Az Px Py Pz ";

	tw::Int i;
	for (i=0;i<refState.size();i++)
		outFile << "real<ref" << i << "|psi> imag<ref" << i << "|psi> ";
}

void KleinGordon::EnergyColumns(std::vector<tw::Float>& cols,std::vector<bool>& avg,const Region& theRgn)
{
	tw::Int i,j,k,s;
	tw::Int x0,x1,y0,y1,z0,z1;
	tw::Float dV,chargeDensity;
	tw::vec3 pos,ENow,ANow,dipoleMoment;
	tw::Complex estate,pstate,psi_ref,chi_ref;
	std::valarray<tw::Complex> overlap(refState.size());

	tw::Float totalCharge = 0.0;
	ENow = 0.0;
	ANow = 0.0;
	dipoleMoment = 0.0;
	overlap = tw::Complex(0,0);

	for (s=0;s<owner->wave.size();s++)
	{
		ANow += owner->wave[s]->VectorPotential(owner->elapsedTime,tw::vec3(0,0,0));
		ENow -= dti*owner->wave[s]->VectorPotential(owner->elapsedTime+0.5*dt,tw::vec3(0,0,0));
		ENow += dti*owner->wave[s]->VectorPotential(owner->elapsedTime-0.5*dt,tw::vec3(0,0,0));
	}

	theRgn.GetLocalCellBounds(&x0,&x1,&y0,&y1,&z0,&z1);
	for (k=z0;k<=z1;k++)
		for (j=y0;j<=y1;j++)
			for (i=x0;i<=x1;i++)
			{
				pos = owner->Pos(i,j,k);
				if (theRgn.Inside(pos,*owner))
				{
					dV = owner->dS(i,j,k,0);
					chargeDensity = J4(i,j,k,0);
					totalCharge += chargeDensity*dV;
					dipoleMoment += chargeDensity*pos*dV;
					// Feshbach-Villars decomposition for overlaps
					estate = FV(i,j,k,1.0);
					pstate = FV(i,j,k,-1.0);
					// Gauge transformation (remove uniform A)
					estate *= std::exp(-ii*q0*(ANow^pos));
					pstate *= std::exp(-ii*q0*(ANow^pos));
					for (s=0;s<refState.size();s++)
					{
						psi_ref = refState[s].Amplitude(pos,0.0,0)/root2;
						chi_ref = refState[s].Amplitude(pos,0.0,1)/root2;
						overlap[s] += conj(psi_ref+chi_ref)*estate*dV;
						overlap[s] -= conj(psi_ref-chi_ref)*pstate*dV;
					}
				}
			}

	cols.push_back(totalCharge); avg.push_back(false);
	cols.push_back(ENow.x); avg.push_back(true);
	cols.push_back(ENow.y); avg.push_back(true);
	cols.push_back(ENow.z); avg.push_back(true);
	cols.push_back(ANow.x); avg.push_back(true);
	cols.push_back(ANow.y); avg.push_back(true);
	cols.push_back(ANow.z); avg.push_back(true);
	cols.push_back(dipoleMoment.x); avg.push_back(false);
	cols.push_back(dipoleMoment.y); avg.push_back(false);
	cols.push_back(dipoleMoment.z); avg.push_back(false);
	for (s=0;s<refState.size();s++)
	{
		cols.push_back(real(overlap[s])); avg.push_back(false);
		cols.push_back(imag(overlap[s])); avg.push_back(false);
	}
}

void KleinGordon::BoxDiagnosticHeader(GridDataDescriptor* box)
{
	owner->WriteBoxDataHeader("rho",box);
	owner->WriteBoxDataHeader("Jx",box);
	owner->WriteBoxDataHeader("Jy",box);
	owner->WriteBoxDataHeader("Jz",box);
	owner->WriteBoxDataHeader("phi",box);
	owner->WriteBoxDataHeader("Ax",box);
	owner->WriteBoxDataHeader("Ay",box);
	owner->WriteBoxDataHeader("Az",box);
	owner->WriteBoxDataHeader("psi0_r",box);
	owner->WriteBoxDataHeader("psi1_r",box);
}

void KleinGordon::BoxDiagnose(GridDataDescriptor* box)
{
	owner->WriteBoxData("rho",box,&J4(0,0,0,0),J4.Stride());
	owner->WriteBoxData("Jx",box,&J4(0,0,0,1),J4.Stride());
	owner->WriteBoxData("Jy",box,&J4(0,0,0,2),J4.Stride());
	owner->WriteBoxData("Jz",box,&J4(0,0,0,3),J4.Stride());

	owner->WriteBoxData("phi",box,&A4(0,0,0,0),A4.Stride());
	owner->WriteBoxData("Ax",box,&A4(0,0,0,1),A4.Stride());
	owner->WriteBoxData("Ay",box,&A4(0,0,0,2),A4.Stride());
	owner->WriteBoxData("Az",box,&A4(0,0,0,3),A4.Stride());

	owner->WriteBoxData("psi0_r",box,&psi_r(0,0,0,0),psi_r.Stride());
	owner->WriteBoxData("psi1_r",box,&psi_r(0,0,0,1),psi_r.Stride());
}

void KleinGordon::StartDiagnostics()
{
	#ifdef USE_OPENCL
	psi_r.ReceiveFromComputeBuffer();
	psi_i.ReceiveFromComputeBuffer();
	A4.ReceiveFromComputeBuffer();
	#endif
	UpdateJ4();
}


///////////////////////////////////////
//                                   //
//   RELATIVISTIC PARTICLE (DIRAC)   //
//                                   //
///////////////////////////////////////


Dirac::Dirac(const std::string& name,Simulation* sim) : AtomicPhysics(name,sim)
{
	// Wavefunction is in standard representation
	// i.e., gamma^0 = diag(1,1,-1,-1)
	typeCode = tw::module_type::dirac;
	m0 = 1.0;
	q0 = -sqrt(alpha);
	residualCharge = sqrt(alpha);
	dipoleApproximation = false;
	psi_r.Initialize(4,*this,owner,tw::dom::xAxis);
	psi_i.Initialize(4,*this,owner,tw::dom::xAxis);

	#ifdef USE_OPENCL
	cl_int err;
	leapFrog = clCreateKernel(program,"DiracAdvance",&err);

	psi_r.InitializeComputeBuffer();
	psi_i.InitializeComputeBuffer();

	#endif
}

Dirac::~Dirac()
{
	#ifdef USE_OPENCL
	clReleaseKernel(leapFrog);
	#endif
}

void Dirac::Initialize()
{
	AtomicPhysics::Initialize();
	psi_r.SetBoundaryConditions(tw::dom::xAxis,fld::neumannWall,fld::none);
	psi_i.SetBoundaryConditions(tw::dom::xAxis,fld::neumannWall,fld::none);

	if (!owner->restarted)
	{
		for (tw::Int s=0;s<waveFunction.size();s++)
		{
			waveFunction[s].Initialize(qo::dirac,residualCharge,q0,m0);
			(*owner->tw_out) << "Wavefunction " << s << " : energy = " << waveFunction[s].energy << std::endl;
		}
		for (tw::Int s=0;s<refState.size();s++)
		{
			refState[s].Initialize(qo::dirac,residualCharge,q0,m0);
			(*owner->tw_out) << "Reference " << s << " : energy = " << refState[s].energy << std::endl;
		}

		#pragma omp parallel
		{
			for (auto cell : InteriorCellRange(*this))
			{
				tw::vec3 pos = owner->Pos(cell);
				for (tw::Int s=0;s<waveFunction.size();s++)
				{
					psi_r(cell,0) += real(waveFunction[s].Amplitude(pos,0.0,0));
					psi_r(cell,1) += real(waveFunction[s].Amplitude(pos,0.0,1));
					psi_r(cell,2) += real(waveFunction[s].Amplitude(pos,dth,2));
					psi_r(cell,3) += real(waveFunction[s].Amplitude(pos,dth,3));
					psi_i(cell,0) += imag(waveFunction[s].Amplitude(pos,0.0,0));
					psi_i(cell,1) += imag(waveFunction[s].Amplitude(pos,0.0,1));
					psi_i(cell,2) += imag(waveFunction[s].Amplitude(pos,dth,2));
					psi_i(cell,3) += imag(waveFunction[s].Amplitude(pos,dth,3));
				}
			}
		}

		psi_r.CopyFromNeighbors();
		psi_i.CopyFromNeighbors();
		FormPotentials(owner->elapsedTime);
		Normalize();
		UpdateJ4();
	}

	#ifdef USE_OPENCL

	psi_r.SendToComputeBuffer();
	psi_i.SendToComputeBuffer();
	Ao4.SendToComputeBuffer();
	A4.SendToComputeBuffer();
	clSetKernelArg(leapFrog,0,sizeof(cl_mem),&psi_r.computeBuffer);
	clSetKernelArg(leapFrog,1,sizeof(cl_mem),&psi_i.computeBuffer);
	clSetKernelArg(leapFrog,2,sizeof(cl_mem),&A4.computeBuffer);
	clSetKernelArg(leapFrog,3,sizeof(cl_mem),&owner->metricsBuffer);
	clSetKernelArg(leapFrog,4,sizeof(m0),&m0);
	clSetKernelArg(leapFrog,5,sizeof(q0),&q0);

	#endif
}

tw::Float Dirac::ComputeRho(const tw::cell& cell)
{
	tw::Float ans = 0.0;
	for (tw::Int c=0;c<4;c++)
		ans += sqr(psi_r(cell,c)) + sqr(psi_i(cell,c));
	return q0*ans;
}

void Dirac::UpdateJ4()
{
	#pragma omp parallel
	{
		tw::Complex z0,z1,z2,z3;
		for (auto cell : InteriorCellRange(*this))
		{
			z0 = tw::Complex(psi_r(cell,0),psi_i(cell,0));
			z1 = tw::Complex(psi_r(cell,1),psi_i(cell,1));
			z2 = tw::Complex(psi_r(cell,2),psi_i(cell,2));
			z3 = tw::Complex(psi_r(cell,3),psi_i(cell,3));
			J4(cell,0) = q0*(norm(z0)+norm(z1)+norm(z2)+norm(z3));
			J4(cell,1) = two*q0*real(conj(z0)*z3 + z1*conj(z2));
			J4(cell,2) = two*q0*imag(conj(z0)*z3 + z1*conj(z2));
			J4(cell,3) = two*q0*real(z0*conj(z2) - z1*conj(z3));
		}
	}
	J4.CopyFromNeighbors();
	J4.ApplyBoundaryCondition();
}

void Dirac::Normalize()
{
	tw::Float totalCharge = 0.0;
	for (auto cell : InteriorCellRange(*this))
		totalCharge += ComputeRho(cell) * owner->dS(cell,0);
	owner->strip[0].AllSum(&totalCharge,&totalCharge,sizeof(tw::Float),0);
	psi_r *= sqrt(fabs(q0/totalCharge));
	psi_i *= sqrt(fabs(q0/totalCharge));
	psi_r.CopyFromNeighbors();
	psi_r.ApplyBoundaryCondition();
	psi_i.CopyFromNeighbors();
	psi_i.ApplyBoundaryCondition();
}

#ifdef USE_OPENCL
void Dirac::Update()
{
	// Solve in the standard representation
	// [gamma^mu(p_mu - qA_mu) - m]psi = 0
	// gamma^0 = | 1  0 |     gamma^i = |  0  t |
	//           | 0 -1 |               | -t  0 |
	// where t is the i^th Pauli matrix
	// Leapfrog Scheme:
	// At start we know psi0(n-1/2),psi1(n-1/2),psi2(n),psi3(n),A(n-1),A(n)
	// Advance psi0,psi1 to n+1/2 using psi2(n),psi3(n),A(n)
	// Advance psi2,psi3 to n+1 using psi0(n+1/2),psi1(n+1/2),A(n+1/2)

	cl_int components[4] = { 0,1,2,3 };
	tw::Float mass[2] = { m0 , -m0 };

	// Update upper pair in the bispinor
	clSetKernelArg(leapFrog,4,sizeof(tw::Float),&mass[0]);
	clSetKernelArg(leapFrog,6,sizeof(cl_int),&components[0]);
	clSetKernelArg(leapFrog,7,sizeof(cl_int),&components[1]);
	clSetKernelArg(leapFrog,8,sizeof(cl_int),&components[2]);
	clSetKernelArg(leapFrog,9,sizeof(cl_int),&components[3]);
	owner->LocalUpdateProtocol(leapFrog);
	// Boundary conditions and communications
	psi_r.UpdateGhostCellsInComputeBuffer(Element(0,1));
	psi_i.UpdateGhostCellsInComputeBuffer(Element(0,1));

	photonPropagator->Advance(A4,Ao4,J4,0.0,dt);
	FormGhostCellPotentials(owner->elapsedTime+dt);
	A4.SendGhostCellsToComputeBuffer();
	photonPropagator->MidstepEstimate(A4,Ao4);

	// Update lower pair in the bispinor
	clSetKernelArg(leapFrog,4,sizeof(tw::Float),&mass[1]);
	clSetKernelArg(leapFrog,6,sizeof(cl_int),&components[2]);
	clSetKernelArg(leapFrog,7,sizeof(cl_int),&components[3]);
	clSetKernelArg(leapFrog,8,sizeof(cl_int),&components[0]);
	clSetKernelArg(leapFrog,9,sizeof(cl_int),&components[1]);
	owner->LocalUpdateProtocol(leapFrog);
	// Boundary conditions and communications
	psi_r.UpdateGhostCellsInComputeBuffer(Element(2,3));
	psi_i.UpdateGhostCellsInComputeBuffer(Element(2,3));

	photonPropagator->UndoMidstepEstimate(A4,Ao4);
}
#else
void Dirac::Update()
{
	// Solve in the standard representation
	// [gamma^mu(p_mu - qA_mu) - m]psi = 0
	// gamma^0 = | 1  0 |     gamma^i = |  0  t |
	//           | 0 -1 |               | -t  0 |
	// where t is the i^th Pauli matrix
	// Leapfrog Scheme:
	// At start we know psi0(n-1/2),psi1(n-1/2),psi2(n),psi3(n),A(n-1),A(n)
	// Advance psi0,psi1 to n+1/2 using psi2(n),psi3(n),A(n)
	// Advance psi2,psi3 to n+1 using psi0(n+1/2),psi1(n+1/2),A(n+1/2)

	LeapFrog<0,1,2,3>(1.0);
	psi_r.CopyFromNeighbors(Element(0,1));
	psi_r.ApplyBoundaryCondition(Element(0,1));
	psi_i.CopyFromNeighbors(Element(0,1));
	psi_i.ApplyBoundaryCondition(Element(0,1));

	photonPropagator->Advance(A4,Ao4,J4,0.0,dt);
	FormGhostCellPotentials(owner->elapsedTime+dt);
	photonPropagator->MidstepEstimate(A4,Ao4);

	LeapFrog<2,3,0,1>(-1.0);
	psi_r.CopyFromNeighbors(Element(2,3));
	psi_r.ApplyBoundaryCondition(Element(2,3));
	psi_i.CopyFromNeighbors(Element(2,3));
	psi_i.ApplyBoundaryCondition(Element(2,3));

	photonPropagator->UndoMidstepEstimate(A4,Ao4);
}
#endif

void Dirac::EnergyHeadings(std::ofstream& outFile)
{
	outFile << "TotalCharge Ex Ey Ez Ax Ay Az Px Py Pz ";

	tw::Int i;
	for (i=0;i<refState.size();i++)
		outFile << "real<ref" << i << "|psi> imag<ref" << i << "|psi> ";
}

void Dirac::EnergyColumns(std::vector<tw::Float>& cols,std::vector<bool>& avg,const Region& theRgn)
{
	tw::Int i,j,k,s;
	tw::Int x0,x1,y0,y1,z0,z1;
	tw::Float dV,chargeDensity;
	tw::vec3 pos,ENow,ANow,dipoleMoment;
	std::valarray<tw::Complex> overlap(refState.size());

	tw::Float totalCharge = 0.0;
	ENow = 0.0;
	ANow = 0.0;
	dipoleMoment = 0.0;
	overlap = tw::Complex(0,0);

	for (s=0;s<owner->wave.size();s++)
	{
		ANow += owner->wave[s]->VectorPotential(owner->elapsedTime,tw::vec3(0,0,0));
		ENow -= dti*owner->wave[s]->VectorPotential(owner->elapsedTime+0.5*dt,tw::vec3(0,0,0));
		ENow += dti*owner->wave[s]->VectorPotential(owner->elapsedTime-0.5*dt,tw::vec3(0,0,0));
	}

	theRgn.GetLocalCellBounds(&x0,&x1,&y0,&y1,&z0,&z1);
	for (k=z0;k<=z1;k++)
		for (j=y0;j<=y1;j++)
			for (i=x0;i<=x1;i++)
			{
				pos = owner->Pos(i,j,k);
				if (theRgn.Inside(pos,*owner))
				{
					dV = owner->dS(i,j,k,0);
					chargeDensity = J4(i,j,k,0);
					totalCharge += chargeDensity*dV;
					dipoleMoment += chargeDensity*pos*dV;
					// Need to revisit statistical interpretation
					for (s=0;s<refState.size();s++)
						overlap[s] += 0.0;
				}
			}

	cols.push_back(totalCharge); avg.push_back(false);
	cols.push_back(ENow.x); avg.push_back(true);
	cols.push_back(ENow.y); avg.push_back(true);
	cols.push_back(ENow.z); avg.push_back(true);
	cols.push_back(ANow.x); avg.push_back(true);
	cols.push_back(ANow.y); avg.push_back(true);
	cols.push_back(ANow.z); avg.push_back(true);
	cols.push_back(dipoleMoment.x); avg.push_back(false);
	cols.push_back(dipoleMoment.y); avg.push_back(false);
	cols.push_back(dipoleMoment.z); avg.push_back(false);
	for (s=0;s<refState.size();s++)
	{
		cols.push_back(real(overlap[s])); avg.push_back(false);
		cols.push_back(imag(overlap[s])); avg.push_back(false);
	}
}

void Dirac::BoxDiagnosticHeader(GridDataDescriptor* box)
{
	owner->WriteBoxDataHeader("rho",box);
	owner->WriteBoxDataHeader("Jx",box);
	owner->WriteBoxDataHeader("Jy",box);
	owner->WriteBoxDataHeader("Jz",box);
	owner->WriteBoxDataHeader("phi",box);
	owner->WriteBoxDataHeader("Ax",box);
	owner->WriteBoxDataHeader("Ay",box);
	owner->WriteBoxDataHeader("Az",box);
	owner->WriteBoxDataHeader("psi0_r",box);
	owner->WriteBoxDataHeader("psi1_r",box);
	owner->WriteBoxDataHeader("psi2_r",box);
	owner->WriteBoxDataHeader("psi3_r",box);
}

void Dirac::BoxDiagnose(GridDataDescriptor* box)
{
	owner->WriteBoxData("rho",box,&J4(0,0,0,0),J4.Stride());
	owner->WriteBoxData("Jx",box,&J4(0,0,0,1),J4.Stride());
	owner->WriteBoxData("Jy",box,&J4(0,0,0,2),J4.Stride());
	owner->WriteBoxData("Jz",box,&J4(0,0,0,3),J4.Stride());

	owner->WriteBoxData("phi",box,&A4(0,0,0,0),A4.Stride());
	owner->WriteBoxData("Ax",box,&A4(0,0,0,1),A4.Stride());
	owner->WriteBoxData("Ay",box,&A4(0,0,0,2),A4.Stride());
	owner->WriteBoxData("Az",box,&A4(0,0,0,3),A4.Stride());

	owner->WriteBoxData("psi0_r",box,&psi_r(0,0,0,0),psi_r.Stride());
	owner->WriteBoxData("psi1_r",box,&psi_r(0,0,0,1),psi_r.Stride());
	owner->WriteBoxData("psi2_r",box,&psi_r(0,0,0,2),psi_r.Stride());
	owner->WriteBoxData("psi3_r",box,&psi_r(0,0,0,3),psi_r.Stride());
}

void Dirac::PointDiagnosticHeader(std::ofstream& outFile)
{
	outFile << "phi Ax Ay Az psi0_r psi0_i psi1_r psi1_i psi2_r psi2_i psi3_r psi3_i ";
}

void Dirac::PointDiagnose(std::ofstream& outFile,const weights_3D& w)
{
	std::valarray<tw::Float> ANow(4),psiNow_r(4),psiNow_i(4);
	A4.Interpolate(ANow,w);
	psi_r.Interpolate(psiNow_r,w);
	psi_i.Interpolate(psiNow_i,w);
	for (tw::Int c=0;c<4;c++)
		outFile << ANow[c] << " ";
	for (tw::Int c=0;c<4;c++)
		outFile << psiNow_r[c] << " " << psiNow_i[c] << " ";
}

void Dirac::StartDiagnostics()
{
	#ifdef USE_OPENCL
	psi_r.ReceiveFromComputeBuffer();
	psi_i.ReceiveFromComputeBuffer();
	A4.ReceiveFromComputeBuffer();
	#endif
	UpdateJ4();
}
