#include "sim.h"
#include "fieldSolve.h"
#include "fluid.h"

///////////////////////
//                   //
//   FLUID MODULE    //
//                   //
///////////////////////


Fluid::Fluid(Grid* theGrid):Module(theGrid)
{
	name = "Fluid";
	typeCode = fluidFields;
	charge = -1.0;
	mass = 1.0;
	thermalMomentum = 0.0;
	enCrossSection = 0.0;
	initialIonizationFraction = 1.0;
	coulombCollisions = false;

	state0.Initialize(4,*this,owner);
	state1.Initialize(4,*this,owner);
	fixed.Initialize(*this,owner);
	gas.Initialize(*this,owner);

	vel.Initialize(4,*this,owner);

	EM = J4 = NULL;
	laser = NULL;
	chi = NULL;
	carrierFrequency = NULL;
}

bool Fluid::InspectResource(void* resource,const std::string& description)
{
	if (description=="electromagnetic:F")
	{
		EM = (Field*)resource;
		return true;
	}

	if (description=="electromagnetic:sources")
	{
		J4 = (Field*)resource;
		return true;
	}

	if (description=="laser:F")
	{
		laser = (Field*)resource;
		return true;
	}

	if (description=="laser:chi")
	{
		chi = (ComplexField*)resource;
		return true;
	}

	if (description=="laser:carrierFrequency")
	{
		carrierFrequency = (tw::Float*)resource;
		return true;
	}

	return false;
}

void Fluid::Initialize()
{
	Module::Initialize();
	ionization.Initialize(owner->unitDensityCGS,carrierFrequency);

	if (owner->restarted)
		return;

	boundarySpec xNormal,xParallel,other;
	xNormal = owner->bc0[1]==axisymmetric || owner->bc0[1]==reflecting ? dirichletWall : neumannWall;
	xParallel = neumannWall;
	other = neumannWall;

	if (owner->movingWindow)
	{
		state0.SetBoundaryConditions(xAxis,xParallel,other);
		state0.SetBoundaryConditions(yAxis,other,other);
		state0.SetBoundaryConditions(zAxis,other,none);
		state0.SetBoundaryConditions(Element(1),xAxis,xNormal,other);

		state1.SetBoundaryConditions(xAxis,xParallel,other);
		state1.SetBoundaryConditions(yAxis,other,other);
		state1.SetBoundaryConditions(zAxis,other,none);
		state1.SetBoundaryConditions(Element(1),xAxis,xNormal,other);

		gas.SetBoundaryConditions(xAxis,xParallel,other);
		gas.SetBoundaryConditions(yAxis,other,other);
		gas.SetBoundaryConditions(zAxis,other,none);

		vel.SetBoundaryConditions(xAxis,xParallel,other);
		vel.SetBoundaryConditions(yAxis,other,other);
		vel.SetBoundaryConditions(zAxis,other,other);
	}
	else
	{
		state0.SetBoundaryConditions(xAxis,neumannWall,neumannWall);
		state0.SetBoundaryConditions(yAxis,neumannWall,neumannWall);
		state0.SetBoundaryConditions(zAxis,neumannWall,neumannWall);
		state0.SetBoundaryConditions(Element(1),xAxis,xNormal,neumannWall);
		state0.SetBoundaryConditions(Element(3),zAxis,dirichletWall,dirichletWall);

		state1.SetBoundaryConditions(xAxis,neumannWall,neumannWall);
		state1.SetBoundaryConditions(yAxis,neumannWall,neumannWall);
		state1.SetBoundaryConditions(zAxis,neumannWall,neumannWall);
		state1.SetBoundaryConditions(Element(1),xAxis,xNormal,neumannWall);
		state1.SetBoundaryConditions(Element(3),zAxis,dirichletWall,dirichletWall);

		gas.SetBoundaryConditions(xAxis,neumannWall,neumannWall);
		gas.SetBoundaryConditions(yAxis,neumannWall,neumannWall);
		gas.SetBoundaryConditions(zAxis,neumannWall,neumannWall);

		vel.SetBoundaryConditions(xAxis,neumannWall,neumannWall);
		vel.SetBoundaryConditions(yAxis,neumannWall,neumannWall);
		vel.SetBoundaryConditions(zAxis,neumannWall,neumannWall);
	}

	#pragma omp parallel
	{
		tw::vec3 pos,A0,A1;
		tw::Float density;
		UnitConverter uc(owner->unitDensityCGS);

		for (CellIterator cell(*this,true);cell<cell.end();++cell)
		{
			pos = owner->Pos(cell);

			for (tw::Int s=0;s<profile.size();s++)
			{
				for (tw::Int c=1;c<=3;c++)
				{
					state0(cell,c) = profile[s]->driftMomentum[c-1];
					state1(cell,c) = profile[s]->driftMomentum[c-1];
				}
				density = profile[s]->GetValue(pos,*owner);
				gas(cell) += (1.0 - initialIonizationFraction)*density;
				state0(cell,0) += initialIonizationFraction*density; // ionization fraction should not be zero
				state1(cell,0) += initialIonizationFraction*density;
				if (owner->neutralize)
					fixed(cell) += initialIonizationFraction*density;
			}

			A0 = A1 = tw::vec3(0,0,0);
			for (tw::Int s=0;s<owner->wave.size();s++)
			{
				A0 += owner->wave[s]->VectorPotential(-dth,pos);
				A1 += owner->wave[s]->VectorPotential(0.0,pos);
			}
			state0(cell,1) += A0.x;
			state0(cell,2) += A0.y;
			state0(cell,3) += 0.5*(A0.x*A0.x + A0.y*A0.y);
			state1(cell,1) += A1.x;
			state1(cell,2) += A1.y;
			state1(cell,3) += 0.5*(A1.x*A1.x + A1.y*A1.y);
			for (tw::Int s=0;s<owner->pulse.size();s++)
			{
				state0(cell,3) += 0.25*norm(owner->pulse[s]->VectorPotentialEnvelope(-dth,pos));
				state1(cell,3) += 0.25*norm(owner->pulse[s]->VectorPotentialEnvelope(0.0,pos));
			}
		}
	}

	// temperature is set to last profile's temperature
	thermalMomentum = profile.back()->thermalMomentum.x;
	if (profile.back()->temperature!=0.0)
		thermalMomentum = sqrt(profile.back()->temperature*mass); // appropriate for exp(-v^2/(2*vth^2)) convention
	if (thermalMomentum==0.0)
		throw tw::FatalError("Fluid module requires temperature specification.");
	if (initialIonizationFraction<=0.0 || initialIonizationFraction>1.0)
		throw tw::FatalError("Fluid module reports initial ionization fraction out of range.");

	state0.CopyFromNeighbors();
	state0.ApplyBoundaryCondition();
	state1.CopyFromNeighbors();
	state1.ApplyBoundaryCondition();
	gas.CopyFromNeighbors();
	gas.ApplyBoundaryCondition();
}

void Fluid::MoveWindow()
{
	#pragma omp parallel
	{
		for (StripIterator s(*this,3,strongbool::yes);s<s.end();++s)
		{
			tw::Int k = s.Dim()+1;
			tw::vec3 pos,A0,A1;
			tw::Float incomingGas,incomingPlasma[4];
			pos = owner->Pos(s,k);
			incomingGas = incomingPlasma[0] = incomingPlasma[1] = incomingPlasma[2] = incomingPlasma[3] = 0.0;
			for (tw::Int p=0;p<profile.size();p++)
			{
				if (ionization.ionizationModel==noIonization)
					incomingPlasma[0] += profile[p]->GetValue(pos,*owner);
				else
				{
					incomingGas += profile[p]->GetValue(pos,*owner);
					incomingPlasma[0] += 1e-6*profile[p]->GetValue(pos,*owner); // add a little plasma
				}
			}
			state0.Shift(s,-1,incomingPlasma);
			state1.Shift(s,-1,incomingPlasma);
			fixed.Shift(s,-1,0.0);
			gas.Shift(s,-1,incomingGas);
			if (owner->neutralize)
				fixed(s,k) += incomingPlasma[0];

			A0 = A1 = tw::vec3(0,0,0);
			for (tw::Int w=0;w<owner->wave.size();w++)
			{
				A0 += owner->wave[w]->VectorPotential(owner->elapsedTime-dth,pos);
				A1 += owner->wave[w]->VectorPotential(owner->elapsedTime,pos);
			}
			state0(s,k,1) += A0.x;
			state0(s,k,2) += A0.y;
			state0(s,k,3) += 0.5*(A0.x*A0.x + A0.y*A0.y);
			state1(s,k,1) += A1.x;
			state1(s,k,2) += A1.y;
			state1(s,k,3) += 0.5*(A1.x*A1.x + A1.y*A1.y);
		}
	}

	state0.DownwardCopy(zAxis,1);
	state1.DownwardCopy(zAxis,1);
	fixed.DownwardCopy(zAxis,1);
	gas.DownwardCopy(zAxis,1);
}

void Fluid::AddDensity(tw::Float densToAdd,tw::Int i,tw::Int j,tw::Int k)
{
	if (densToAdd>0.0)
	{
		tw::Float scaleFactor = state1(i,j,k,0) / (state1(i,j,k,0) + densToAdd);
		state1(i,j,k,0) += densToAdd;
		state1(i,j,k,1) *= scaleFactor;
		state1(i,j,k,2) *= scaleFactor;
		state1(i,j,k,3) *= scaleFactor;
		if (owner->neutralize)
			fixed(i,j,k) += densToAdd;
	}
}

void Fluid::Update()
{
	boundarySpec bc;
	tw::Float field,ionizedDensity;
	UnitConverter uc(owner->unitDensityCGS);

	const tw::Float q0 = charge;
	const tw::Float m0 = mass;
	const tw::Float m0i = 1.0/m0;
	const tw::Float kT = sqr(thermalMomentum)/m0; // appropriate for exp(-v^2/(2*vth^2)) convention
	const tw::Float coulomb = coulombCollisions ? 1.0 : 0.0;

	// Time centering information upon entry:
	// qty		time		qty		time
	// ---		----		---		----
	// n0		-1/2		n1		0
	// 						p1		-1/2
	// a0		-1/2		a1		1/2
	// aa0		0			aa1		1/2
	// E		0			B		1/2

	// Time centering information on output:
	// qty		time		qty					time
	// ---		----		---					----
	// rho		1			chi					1/2
	// J		1/2

	// Advance the momentum to t=1/2

	#pragma omp parallel
	{
		// Apply half of E-field impulse and put relativistic mass into vel(cell,0)
		for (VectorizingIterator<3> v(*this,false);v<v.end();++v)
		{
			#pragma omp simd
			for (tw::Int k=1;k<=dim[3];k++)
			{
				state1(v,k,1) += dth*q0*(*EM).sfwd(v,k,0,1);
				state1(v,k,2) += dth*q0*(*EM).sfwd(v,k,1,2);
				state1(v,k,3) += dth*q0*(*EM).sfwd(v,k,2,3);
				vel(v,k,0) = m0*m0 + sqr(state1(v,k,1)) + sqr(state1(v,k,2)) + sqr(state1(v,k,3));
			}
			if (laser)
			{
				#pragma omp simd
				for (tw::Int k=1;k<=dim[3];k++)
					vel(v,k,0) += 0.5*q0*q0*(*laser)(v,k,6);
			}
			#pragma omp simd
			for (tw::Int k=1;k<=dim[3];k++)
				vel(v,k,0) = sqrt(vel(v,k,0));
		}
	}

	vel.CopyFromNeighbors(Element(0));
	vel.ApplyBoundaryCondition(Element(0));
	#pragma omp parallel
	{
		tw::Float kT_eff,temp;
		std::valarray<tw::Float> nuColl(dim[3]+1);
		for (VectorizingIterator<3> v(*this,false);v<v.end();++v)
		{
			#pragma omp simd
			for (tw::Int k=1;k<=dim[3];k++)
			{
				kT_eff = kT + (vel(v,k,0) - m0);
				nuColl[k] = gas(v,k) * enCrossSection * sqrt(kT_eff*m0i);
				temp = coulomb * 1e-5 * uc.SimToCGS(density_dim,fixed(v,k));
				temp *= pow(uc.sim_to_eV(kT_eff),-1.5);
				temp *= uc.CGSValue(time_dim);
				nuColl[k] += temp;
			}
			#pragma omp simd
			for (tw::Int k=1;k<=dim[3];k++)
			{
				state1(v,k,1) -= dt*vel(v,k,0,1);
				state1(v,k,2) -= dt*vel(v,k,0,2);
				state1(v,k,3) -= dt*vel(v,k,0,3);
				state1(v,k,1) += dth*q0*(*EM).sfwd(v,k,0,1);
				state1(v,k,2) += dth*q0*(*EM).sfwd(v,k,1,2);
				state1(v,k,3) += dth*q0*(*EM).sfwd(v,k,2,3);
				state1(v,k,1) /= 1.0 + nuColl[k]*dt;
				state1(v,k,2) /= 1.0 + nuColl[k]*dt;
				state1(v,k,3) /= 1.0 + nuColl[k]*dt;
			}
		}
	}
	state1.CopyFromNeighbors();
	state1.ApplyBoundaryCondition();

	// Compute 1/mass and velocity at t = 1/2

	#pragma omp parallel
	{
		for (VectorizingIterator<3> v(*this,true);v<v.end();++v)
		{
			#pragma omp simd
			for (tw::Int k=lb[3];k<=ub[3];k++)
				vel(v,k,0) = m0*m0 + sqr(state1(v,k,1)) + sqr(state1(v,k,2)) + sqr(state1(v,k,3));
			if (laser)
			{
				#pragma omp simd
				for (tw::Int k=lb[3];k<=ub[3];k++)
					vel(v,k,0) += 0.5*q0*q0*(*laser)(v,k,7);
			}
			#pragma omp simd
			for (tw::Int k=lb[3];k<=ub[3];k++)
			{
				vel(v,k,0) = 1.0/sqrt(vel(v,k,0));
				vel(v,k,1) = state1(v,k,1)*vel(v,k,0);
				vel(v,k,2) = state1(v,k,2)*vel(v,k,0);
				vel(v,k,3) = state1(v,k,3)*vel(v,k,0);
			}
		}
	}

	// FCT Update of density - trial step

	state0 = state1; // dens0(t=0) , dens1(t=0)
	FCT_Driver convector(&state0,&state1,&vel,NULL,owner);
	convector.SetDensityElements(Element(0));
	if (dim[1]>1)
	{
		convector.SetVelocityElement(1);
		convector.Convect(xAxis,dirichletCell,dirichletCell,dth);
	}
	if (dim[2]>1)
	{
		convector.SetVelocityElement(2);
		convector.Convect(yAxis,dirichletCell,dirichletCell,dth);
	}
	if (dim[3]>1)
	{
		bc = owner->movingWindow ? neumannWall : dirichletCell;
		convector.SetVelocityElement(3);
		convector.Convect(zAxis,bc,bc,dth);
	}

	// FCT Update of density - full step

	Swap(state0,state1); // dens0(t=0) , dens1(t=1/2)
	if (dim[1]>1)
	{
		convector.SetVelocityElement(1);
		convector.Convect(xAxis,dirichletCell,dirichletCell,dt);
		convector.GetTrueFlux(vel,Element(1),Element(0));
	}
	if (dim[2]>1)
	{
		convector.SetVelocityElement(2);
		convector.Convect(yAxis,dirichletCell,dirichletCell,dt);
		convector.GetTrueFlux(vel,Element(2),Element(0));
	}
	if (dim[3]>1)
	{
		bc = owner->movingWindow ? neumannWall : dirichletCell;
		convector.SetVelocityElement(3);
		convector.Convect(zAxis,bc,bc,dt);
		convector.GetTrueFlux(vel,Element(3),Element(0));
	}
	Swap(state0,state1); // dens0(t=1/2) , dens1(t=1)

	// Ionization

	if (ionization.ionizationModel!=noIonization)
	{
		for (tw::Int i=lb[1];i<=ub[1];i++)
			for (tw::Int j=lb[2];j<=ub[2];j++)
				for (tw::Int k=lb[3];k<=ub[3];k++)
				{
					field = sqrt(sqr((*EM)(i,j,k,0))+sqr((*EM)(i,j,k,1))+sqr((*EM)(i,j,k,2)));
					ionizedDensity = gas(i,j,k)*dt*ionization.Rate(field,0.0);
					if (ionizedDensity > gas(i,j,k))
						ionizedDensity = gas(i,j,k);
					gas(i,j,k) = fabs(gas(i,j,k) - ionizedDensity);
					AddDensity(ionizedDensity,i,j,k);
				}
	}

	// Deposit current at t = 1/2 and charge at t = 1
	// Get conservative current based on actual FCT fluxes which are stored
	// in vel(cell,i) ; these are currents integrated over a cell wall.

	if (J4)
	{
		// This is a deposition routine.
		// Writing to ghost cells would double count the sources.
		#pragma omp parallel
		{
			// Need special treatment for current in ignorable directions.
			// Here cell walls and centers need not be distinguished (dim=1).
			if (dim[1]==1)
			{
				for (VectorizingIterator<3> v(*this,false);v<v.end();++v)
					#pragma omp simd
					for (tw::Int k=1;k<=dim[3];k++)
						vel(v,k,1) *= owner->dS(v,k,1) * dt * state0(v,k,0);
			}
			if (dim[2]==1)
			{
				for (VectorizingIterator<3> v(*this,false);v<v.end();++v)
					#pragma omp simd
					for (tw::Int k=1;k<=dim[3];k++)
						vel(v,k,2) *= owner->dS(v,k,2) * dt * state0(v,k,0);
			}
			if (dim[3]==1)
			{
				for (VectorizingIterator<3> v(*this,false);v<v.end();++v)
					#pragma omp simd
					for (tw::Int k=1;k<=dim[3];k++)
						vel(v,k,3) *= owner->dS(v,k,3) * dt * state0(v,k,0);
			}
		}
		#pragma omp parallel
		{
			for (VectorizingIterator<3> v(*this,false);v<v.end();++v)
			{
				#pragma omp simd
				for (tw::Int k=1;k<=dim[3];k++)
				{
					(*J4)(v,k,0) += owner->dS(v,k,0) * q0 * (state1(v,k,0) - fixed(v,k));
					(*J4)(v,k,1) += q0 * vel(v,k,1) * dti;
					(*J4)(v,k,2) += q0 * vel(v,k,2) * dti;
					(*J4)(v,k,3) += q0 * vel(v,k,3) * dti;
				}
			}
		}
	}

	if (chi)
	{
		#pragma omp parallel
		{
			for (CellIterator cell(*this,false);cell<cell.end();++cell)
				(*chi)(cell) -= owner->dS(cell,0)*q0*q0*state0(cell,0)/(m0*vel(cell,0));
		}
	}
}

void Fluid::ReadInputFileTerm(std::stringstream& inputString,std::string& command)
{
	std::string word;
	UnitConverter uc(owner->unitDensityCGS);

	Module::ReadInputFileTerm(inputString,command);
	ionization.ReadInputFileTerm(inputString,command);
	if (command=="charge") // eg, charge = 1.0
	{
		inputString >> word;
		inputString >> charge;
	}
	if (command=="mass") // eg, mass = 1800.0
	{
		inputString >> word;
		inputString >> mass;
	}
	if (command=="collision")
	{
		throw tw::FatalError("collision frequency variable no longer supported.");
	}
	if (command=="initial") // eg, initial ionization fraction = 1e-3
	{
		inputString >> word >> word >> word;
		inputString >> initialIonizationFraction;
	}
	if (command=="neutral") // eg, neutral cross section = 1e-15 // cm^2
	{
		inputString >> word >> word >> word;
		inputString >> enCrossSection;
		enCrossSection = uc.CGSToSim(cross_section_dim,enCrossSection);
	}
	if (command=="coulomb") // eg, coulomb collisions = true
	{
		inputString >> word >> word >> word;
		coulombCollisions = (word=="true" || word=="on" || word=="yes");
	}
}

void Fluid::ReadData(std::ifstream& inFile)
{
	Module::ReadData(inFile);
	inFile.read((char *)&charge,sizeof(tw::Float));
	inFile.read((char *)&mass,sizeof(tw::Float));
	inFile.read((char *)&thermalMomentum,sizeof(tw::Float));
	inFile.read((char *)&enCrossSection,sizeof(tw::Float));
	inFile.read((char *)&initialIonizationFraction,sizeof(tw::Float));
	inFile.read((char *)&ionization,sizeof(IonizationData));
	inFile.read((char *)&coulombCollisions,sizeof(bool));

	state0.ReadData(inFile);
	state1.ReadData(inFile);
	gas.ReadData(inFile);
	fixed.ReadData(inFile);
}

void Fluid::WriteData(std::ofstream& outFile)
{
	Module::WriteData(outFile);
	outFile.write((char *)&charge,sizeof(tw::Float));
	outFile.write((char *)&mass,sizeof(tw::Float));
	outFile.write((char *)&thermalMomentum,sizeof(tw::Float));
	outFile.write((char *)&enCrossSection,sizeof(tw::Float));
	outFile.write((char *)&initialIonizationFraction,sizeof(tw::Float));
	outFile.write((char *)&ionization,sizeof(IonizationData));
	outFile.write((char *)&coulombCollisions,sizeof(bool));

	state0.WriteData(outFile);
	state1.WriteData(outFile);
	gas.WriteData(outFile);
	fixed.WriteData(outFile);
}

void Fluid::BoxDiagnosticHeader(GridDataDescriptor* box)
{
	owner->WriteBoxDataHeader(name+"_e",box);
	owner->WriteBoxDataHeader(name+"_n",box);
}

void Fluid::BoxDiagnose(GridDataDescriptor* box)
{
	owner->WriteBoxData(name+"_e",box,&state1(0,0,0,0),state1.Stride());
	owner->WriteBoxData(name+"_n",box,&gas(0,0,0),gas.Stride());
}

void Fluid::PointDiagnosticHeader(std::ofstream& outFile)
{
	outFile << name << " ";
}

void Fluid::PointDiagnose(std::ofstream& outFile,const weights_3D& w)
{
	std::valarray<tw::Float> dens(4);
	state1.Interpolate(dens,w);
	outFile << dens[0]  << " ";
}


/////////////////////
//                 //
//  INTERACTIONS   //
//                 //
/////////////////////


void SubReaction::ReadData(std::ifstream& inFile)
{
	tw::Int i,num,data;

	inFile.read((char *)&num,sizeof(tw::Int));
	for (i=0;i<num;i++)
	{
		inFile.read((char *)&data,sizeof(tw::Int));
		reactant.push_back(data);
	}

	inFile.read((char *)&num,sizeof(tw::Int));
	for (i=0;i<num;i++)
	{
		inFile.read((char *)&data,sizeof(tw::Int));
		product.push_back(data);
	}

	inFile.read((char *)&heat,sizeof(tw::Float));
	inFile.read((char *)&vheat,sizeof(tw::Float));
}

void SubReaction::WriteData(std::ofstream& outFile)
{
	tw::Int i,num,data;

	num = reactant.size();
	outFile.write((char *)&num,sizeof(tw::Int));
	for (i=0;i<num;i++)
	{
		data = reactant[i];
		outFile.write((char *)&data,sizeof(tw::Int));
	}

	num = product.size();
	outFile.write((char *)&num,sizeof(tw::Int));
	for (i=0;i<num;i++)
	{
		data = product[i];
		outFile.write((char *)&data,sizeof(tw::Int));
	}

	outFile.write((char *)&heat,sizeof(tw::Float));
	outFile.write((char *)&vheat,sizeof(tw::Float));
}

void Reaction::ReadData(std::ifstream& inFile)
{
	tw::Int i,num;

	inFile.read((char *)&catalyst,sizeof(tw::Int));
	inFile.read((char *)&numBodies,sizeof(tw::Int));
	inFile.read((char *)&c1,sizeof(tw::Float));
	inFile.read((char *)&c2,sizeof(tw::Float));
	inFile.read((char *)&c3,sizeof(tw::Float));
	inFile.read((char *)&T0,sizeof(tw::Float));
	inFile.read((char *)&T1,sizeof(tw::Float));
	inFile.read((char *)b,sizeof(tw::Float)*9);

	inFile.read((char *)&num,sizeof(tw::Int));
	for (i=0;i<num;i++)
	{
		sub.push_back(new SubReaction);
		sub.back()->ReadData(inFile);
	}
}

void Reaction::WriteData(std::ofstream& outFile)
{
	tw::Int i,num;

	outFile.write((char *)&catalyst,sizeof(tw::Int));
	outFile.write((char *)&numBodies,sizeof(tw::Int));
	outFile.write((char *)&c1,sizeof(tw::Float));
	outFile.write((char *)&c2,sizeof(tw::Float));
	outFile.write((char *)&c3,sizeof(tw::Float));
	outFile.write((char *)&T0,sizeof(tw::Float));
	outFile.write((char *)&T1,sizeof(tw::Float));
	outFile.write((char *)b,sizeof(tw::Float)*9);

	num = sub.size();
	outFile.write((char *)&num,sizeof(tw::Int));
	for (i=0;i<num;i++)
		sub[i]->WriteData(outFile);
}

void Excitation::ReadData(std::ifstream& inFile)
{
	inFile.read((char *)&exciter,sizeof(tw::Int));
	inFile.read((char *)&excitee,sizeof(tw::Int));
	inFile.read((char *)&c1,sizeof(tw::Float));
	inFile.read((char *)&c2,sizeof(tw::Float));
	inFile.read((char *)&c3,sizeof(tw::Float));
	inFile.read((char *)&T0,sizeof(tw::Float));
	inFile.read((char *)&T1,sizeof(tw::Float));
	inFile.read((char *)b,sizeof(tw::Float)*9);
	inFile.read((char *)&level,sizeof(tw::Float));
}

void Excitation::WriteData(std::ofstream& outFile)
{
	outFile.write((char *)&exciter,sizeof(tw::Int));
	outFile.write((char *)&excitee,sizeof(tw::Int));
	outFile.write((char *)&c1,sizeof(tw::Float));
	outFile.write((char *)&c2,sizeof(tw::Float));
	outFile.write((char *)&c3,sizeof(tw::Float));
	outFile.write((char *)&T0,sizeof(tw::Float));
	outFile.write((char *)&T1,sizeof(tw::Float));
	outFile.write((char *)b,sizeof(tw::Float)*9);
	outFile.write((char *)&level,sizeof(tw::Float));
}

void Collision::ReadData(std::ifstream& inFile)
{
	inFile.read((char *)&type,sizeof(sparc::collisionType));
	inFile.read((char *)&chem1,sizeof(tw::Int));
	inFile.read((char *)&chem2,sizeof(tw::Int));
	inFile.read((char *)&crossSection,sizeof(tw::Float));
	inFile.read((char *)&ks,sizeof(tw::Float));
	inFile.read((char *)&T_ref,sizeof(tw::Float));
	inFile.read((char *)&n_ref,sizeof(tw::Float));
}

void Collision::WriteData(std::ofstream& outFile)
{
	outFile.write((char *)&type,sizeof(sparc::collisionType));
	outFile.write((char *)&chem1,sizeof(tw::Int));
	outFile.write((char *)&chem2,sizeof(tw::Int));
	outFile.write((char *)&crossSection,sizeof(tw::Float));
	outFile.write((char *)&ks,sizeof(tw::Float));
	outFile.write((char *)&T_ref,sizeof(tw::Float));
	outFile.write((char *)&n_ref,sizeof(tw::Float));
}


/////////////////////
//                 //
//   GROUP MODULE  //
//                 //
/////////////////////


EquilibriumGroup::EquilibriumGroup(Grid* theGrid):Module(theGrid)
{
	name = "Group";
	typeCode = equilibriumGroup;
	forceFilter = 1.0;

	// ASHER_MOD
	//eosMixData = (EOSMixture *)owner->AddPrivateTool(eosIdealGasMix); // For bug testing with original EOS calc
	eosMixData = (EOSMixture *)owner->AddPrivateTool(eosMixture);
	//
}

// ASHER_MOD -- EquilibriumGroup now needs a destructor because of eosMixData
EquilibriumGroup::~EquilibriumGroup()
{
	owner->RemoveTool(eosMixData); 
}

void EquilibriumGroup::ReadOneChemical(std::stringstream& inputString)
{
	Chemical *temp = new Chemical(owner);
	temp->chemBoss = chemBoss;
	temp->group = this;
	chemical.push_back(temp);
	chemBoss->chemical.push_back(temp);
	owner->module.push_back(temp);
	temp->name = name;
	name += "(group)";
	temp->ReadInputFileBlock(inputString);
	// ASHER_MOD
	// Select from a number of EOS models based on input
	if (temp->eosName=="IdealGas") {
		if (temp->mass==1.0) {
		temp->eosData = (EOSDataTool*)owner->AddPrivateTool(eosIdealGasElectrons);
		} else {
		temp->eosData = (EOSDataTool*)owner->AddPrivateTool(eosIdealGas);
		}
	}
		
	if (temp->eosName=="MieGruneisen") {
	temp->eosData = (EOSDataTool*)owner->AddPrivateTool(eosMieGruneisen);
	// let the EOS Object know what the input Gruneisen Parameter is
	((EOSMieGruneisen*)(temp->eosData))->GRUN = temp->GRUN;
	}	
	if (temp->eosName=="MieGruneisen2") {
	temp->eosData = (EOSDataTool*)owner->AddPrivateTool(eosMieGruneisen2);
	// let the EOS Object know what the input Gruneisen Parameter is
	((EOSMieGruneisen2*)(temp->eosData))->GRUN = temp->GRUN;
	((EOSMieGruneisen2*)(temp->eosData))->n0 = temp->n0;
	((EOSMieGruneisen2*)(temp->eosData))->c0 = temp->c0;
	((EOSMieGruneisen2*)(temp->eosData))->S1 = temp->S1;
	}	
	eosMixData->eosComponents.push_back(temp->eosData);
}

void EquilibriumGroup::ReadInputFileBlock(std::stringstream& inputString)
{
	Chemical *temp;
	std::string word;
	do {
		inputString >> word;
		if (word=="new")
		{
			inputString >> word;
			temp = new Chemical(owner);
			temp->chemBoss = chemBoss;
			temp->group = this;
			chemical.push_back(temp);
			chemBoss->chemical.push_back(temp);
			owner->module.push_back(temp);
			inputString >> temp->name;
			temp->ReadInputFileBlock(inputString);
			// ASHER_MOD
			// Select from a number of EOS models based on input
			if (temp->eosName=="IdealGas") {
				if (temp->mass==1.0) {
				temp->eosData = (EOSDataTool*)owner->AddPrivateTool(eosIdealGasElectrons);
				} else {
				temp->eosData = (EOSDataTool*)owner->AddPrivateTool(eosIdealGas);
				}
			}
			if (temp->eosName=="MieGruneisen") {
				temp->eosData = (EOSDataTool*)owner->AddPrivateTool(eosMieGruneisen);
				// let the EOS Object know what the input Gruneisen Parameter is
				((EOSMieGruneisen*)(temp->eosData))->GRUN = temp->GRUN;
			}	
			if (temp->eosName=="MieGruneisen2") {
			temp->eosData = (EOSDataTool*)owner->AddPrivateTool(eosMieGruneisen2);
				// let the EOS Object know what the input Gruneisen Parameter is
				((EOSMieGruneisen2*)(temp->eosData))->GRUN = temp->GRUN;
				((EOSMieGruneisen2*)(temp->eosData))->n0 = temp->n0;
				((EOSMieGruneisen2*)(temp->eosData))->c0 = temp->c0;
				((EOSMieGruneisen2*)(temp->eosData))->S1 = temp->S1;
			}	
			eosMixData->eosComponents.push_back(temp->eosData);

		}
		if (word=="mobile")
		{
			inputString >> word >> word;
			forceFilter = (word=="yes" || word=="on" || word=="true") ? 1.0 : 0.0;
		}
	} while (word!="}");
}

void EquilibriumGroup::Initialize()
{
	tw::Int c;

	Module::Initialize();
	mass.resize(chemical.size());
	charge.resize(chemical.size());
	cvm.resize(chemical.size());
	excitationEnergy.resize(chemical.size());
	thermo_cond_cvm.resize(chemical.size());
	k_visc_m.resize(chemical.size());

	for (c=0;c<chemical.size();c++)
	{
		mass[c] = chemical[c]->mass;
		charge[c] = chemical[c]->charge;
		cvm[c] = chemical[c]->cvm;
		excitationEnergy[c] = chemical[c]->excitationEnergy;
		thermo_cond_cvm[c] = chemical[c]->thermometricConductivity * cvm[c];
		k_visc_m[c] = chemical[c]->kinematicViscosity * mass[c];
	}
}

void EquilibriumGroup::ReadData(std::ifstream& inFile)
{
	tw::Int num;
	Module::ReadData(inFile);
	inFile.read((char *)&num,sizeof(tw::Int));
	mass.resize(num);
	charge.resize(num);
	cvm.resize(num);
	excitationEnergy.resize(num);
	thermo_cond_cvm.resize(num);
	k_visc_m.resize(num);

	inFile.read((char *)&mass[0],sizeof(tw::Float)*num);
	inFile.read((char *)&charge[0],sizeof(tw::Float)*num);
	inFile.read((char *)&cvm[0],sizeof(tw::Float)*num);
	inFile.read((char *)&excitationEnergy[0],sizeof(tw::Float)*num);
	inFile.read((char *)&thermo_cond_cvm[0],sizeof(tw::Float)*num);
	inFile.read((char *)&k_visc_m[0],sizeof(tw::Float)*num);
	inFile.read((char *)&forceFilter,sizeof(tw::Float));

	inFile.read((char *)&e,sizeof(Element));
	inFile.read((char *)&npx,sizeof(tw::Int));
	inFile.read((char *)&npy,sizeof(tw::Int));
	inFile.read((char *)&npz,sizeof(tw::Int));
	inFile.read((char *)&U,sizeof(tw::Int));
	inFile.read((char *)&Xi,sizeof(tw::Int));

	inFile.read((char *)&T,sizeof(tw::Int));
	inFile.read((char *)&Tv,sizeof(tw::Int));
	inFile.read((char *)&P,sizeof(tw::Int));
	inFile.read((char *)&K,sizeof(tw::Int));
	inFile.read((char *)&visc,sizeof(tw::Int));
}

void EquilibriumGroup::WriteData(std::ofstream& outFile)
{
	tw::Int num = chemical.size();
	Module::WriteData(outFile);
	outFile.write((char *)&num,sizeof(tw::Int));

	outFile.write((char *)&mass[0],sizeof(tw::Float)*num);
	outFile.write((char *)&charge[0],sizeof(tw::Float)*num);
	outFile.write((char *)&cvm[0],sizeof(tw::Float)*num);
	outFile.write((char *)&excitationEnergy[0],sizeof(tw::Float)*num);
	outFile.write((char *)&thermo_cond_cvm[0],sizeof(tw::Float)*num);
	outFile.write((char *)&k_visc_m[0],sizeof(tw::Float)*num);
	outFile.write((char *)&forceFilter,sizeof(tw::Float));

	outFile.write((char *)&e,sizeof(Element));
	outFile.write((char *)&npx,sizeof(tw::Int));
	outFile.write((char *)&npy,sizeof(tw::Int));
	outFile.write((char *)&npz,sizeof(tw::Int));
	outFile.write((char *)&U,sizeof(tw::Int));
	outFile.write((char *)&Xi,sizeof(tw::Int));

	outFile.write((char *)&T,sizeof(tw::Int));
	outFile.write((char *)&Tv,sizeof(tw::Int));
	outFile.write((char *)&P,sizeof(tw::Int));
	outFile.write((char *)&K,sizeof(tw::Int));
	outFile.write((char *)&visc,sizeof(tw::Int));
}

void EquilibriumGroup::BoxDiagnosticHeader(GridDataDescriptor* box)
{
	std::string filename;
	filename = "T_" + name;
	owner->WriteBoxDataHeader(filename,box);
	filename = "P_" + name;                  // ASHER_MOD -- add pressure to diagnostic
	owner->WriteBoxDataHeader(filename,box);
	filename = "Cv_" + name;                 // ASHER_MOD -- add heat capacitance to the diagnostic
	owner->WriteBoxDataHeader(filename,box);
	filename = "K_" + name;
	owner->WriteBoxDataHeader(filename,box);
	filename = "X_" + name;
	owner->WriteBoxDataHeader(filename,box);
	filename = "vx_" + name;
	owner->WriteBoxDataHeader(filename,box);
	filename = "vy_" + name;
	owner->WriteBoxDataHeader(filename,box);
	filename = "vz_" + name;
	owner->WriteBoxDataHeader(filename,box);
}

void EquilibriumGroup::BoxDiagnose(GridDataDescriptor* box)
{
	tw::Int ax;
	std::string filename;
	char ax_lab[3] = { 'x' , 'y' , 'z' };

	filename = "T_" + name;
	CopyFieldData(chemBoss->scratch,Element(0),chemBoss->eos1,Element(T));
	chemBoss->scratch *= chemBoss->fluxMask;
	owner->WriteBoxData(filename,box,&chemBoss->scratch(0,0,0),chemBoss->scratch.Stride());

	filename = "P_" + name;
	CopyFieldData(chemBoss->scratch,Element(0),chemBoss->eos1,Element(P));
	chemBoss->scratch *= chemBoss->fluxMask;
	owner->WriteBoxData(filename,box,&chemBoss->scratch(0,0,0),chemBoss->scratch.Stride());

	filename = "Cv_" + name;
	CopyFieldData(chemBoss->scratch,Element(0),chemBoss->eos1,Element(Cv));
	chemBoss->scratch *= chemBoss->fluxMask;
	owner->WriteBoxData(filename,box,&chemBoss->scratch(0,0,0),chemBoss->scratch.Stride());

	filename = "K_" + name;
	CopyFieldData(chemBoss->scratch,Element(0),chemBoss->eos1,Element(K));
	chemBoss->scratch *= chemBoss->fluxMask;
	owner->WriteBoxData(filename,box,&chemBoss->scratch(0,0,0),chemBoss->scratch.Stride());

	filename = "X_" + name;
	CopyFieldData(chemBoss->scratch,Element(0),chemBoss->eos1,Element(Tv));
	chemBoss->scratch *= chemBoss->fluxMask;
	owner->WriteBoxData(filename,box,&chemBoss->scratch(0,0,0),chemBoss->scratch.Stride());

	for (ax=1;ax<=3;ax++)
	{
		filename = ax_lab[ax-1];
		filename = "v" + filename + "_" + name;
		LoadVelocity(chemBoss->scratch,chemBoss->state1,ax);
		chemBoss->scratch *= chemBoss->fluxMask;
		owner->WriteBoxData(filename,box,&chemBoss->scratch(0,0,0),chemBoss->scratch.Stride());
		//owner->WriteBoxData(filename,box,&chemBoss->state1(0,0,0,npx+ax-1),chemBoss->state1.Stride());
	}
}

void EquilibriumGroup::PointDiagnosticHeader(std::ofstream& outFile)
{
	outFile << ("T_"+name) << " " << ("X_"+name) << " ";
}

void EquilibriumGroup::PointDiagnose(std::ofstream& outFile,const weights_3D& w)
{
	std::valarray<tw::Float> data(chemBoss->state1.Components());
	chemBoss->state1.Interpolate(data,Element(U,Xi),w);
	outFile << data[U] << " " << data[Xi] << " ";
}


/////////////////////
//                 //
// CHEMICAL MODULE //
//                 //
/////////////////////


Chemical::Chemical(Grid* theGrid):Module(theGrid)
{
	name = "Chemical";
	typeCode = chemical;
	charge = -1.0;
	mass = 1.0;
	cvm = 1.5;
	excitationEnergy = 0.0;
	thermometricConductivity = kinematicViscosity = 0.0;
	effectiveMass = 1.0; // for electrons moving though this chemical
	transitionDensity = 1.0; // effective mass comes in for density >> transitionDensity
	permittivity = tw::Complex(1.0,0.0);
	indexInGroup = 0;

	// ASHER_MOD
	eosName = "IdealGas"; // default unless specified otherwise
	GRUN    = 2.0;        // currently it defaults to the value of Cu at room temperature (0.1 for H20),
	n0 = 3.3e3;           // not that there is any reason why it should.
	c0 = 1.3248e-5;       // Note these parameters are meaningless if you're using the IdealGas model.
	S1 = 1.5;
}

// ASHER_MOD -- Chemical now needs a destructor because of the eosData
Chemical::~Chemical() 
{
	// ASHER_MOD
	owner->RemoveTool(eosData); 
}

void Chemical::Initialize()
{
	Module::Initialize();
	if (ionization.ionizationModel!=noIonization)
		ionization.Initialize(owner->unitDensityCGS,&chemBoss->laserFrequency);
	GenerateFluid(chemBoss->state1,true);
	// can't generate in Chemistry::Initialize since profiles would not be initialized
}

bool Chemical::GenerateFluid(Field& f,bool init)
{
	tw::Int i,j,k,s;
	tw::vec3 pos,np;
	tw::Float add,kT,dens,kinetic,vibrational;

	bool timeGate,didGenerate=false;

	const tw::Int ns = indexInState;
	const tw::Int npx = group->npx;
	const tw::Int npy = group->npy;
	const tw::Int npz = group->npz;
	const tw::Int U = group->U;
	const tw::Int Xi = group->Xi;

	for (s=0;s<profile.size();s++)
	{
		switch (profile[s]->timingMethod)
		{
			case triggeredProfile:
				timeGate = owner->elapsedTime>=profile[s]->t0 && !profile[s]->wasTriggered;
				add = 1.0;
				break;
			case maintainedProfile:
				timeGate = owner->elapsedTime>=profile[s]->t0 && owner->elapsedTime<=profile[s]->t1;
				add = 0.0;
				break;
			default:
				timeGate = false;
		}
		if ( timeGate )
		{
			profile[s]->wasTriggered = true;
			didGenerate = true;
			for (k=lb[3];k<=ub[3];k++)
				for (j=lb[2];j<=ub[2];j++)
					for (i=lb[1];i<=ub[1];i++)
					{
						pos = owner->Pos(i,j,k);
						dens = profile[s]->GetValue(pos,*owner);
						if (profile[s]->whichQuantity==densityProfile && dens>0.0)
						{
							kT = sqr(profile[s]->thermalMomentum.x)/mass; // appropriate for exp(-v^2/(2*vth^2)) convention
							if (profile[s]->temperature!=0.0)
								kT = profile[s]->temperature;
							np.x = dens*profile[s]->driftMomentum.x;
							np.y = dens*profile[s]->driftMomentum.y;
							np.z = dens*profile[s]->driftMomentum.z;
							kinetic = 0.5*Norm(np)/(tw::small_pos + mass*dens);
							vibrational = dens*excitationEnergy/(fabs(exp(excitationEnergy/kT) - 1.0) + tw::small_pos);

							f(i,j,k,ns) = add*f(i,j,k,ns) + dens;
							f(i,j,k,npx) = add*f(i,j,k,npx) + np.x;
							f(i,j,k,npy) = add*f(i,j,k,npy) + np.y;
							f(i,j,k,npz) = add*f(i,j,k,npz) + np.z;
							f(i,j,k,U) = add*f(i,j,k,U) + cvm*dens*kT + kinetic + vibrational;
							f(i,j,k,Xi) = add*f(i,j,k,Xi) + vibrational;
						}
						if (profile[s]->whichQuantity==energyProfile)
						{
							f(i,j,k,U) = add*f(i,j,k,U) + dens;
						}
						if (profile[s]->whichQuantity==pxProfile)
						{
							f(i,j,k,npx) = add*f(i,j,k,npx) + dens;
						}
						if (profile[s]->whichQuantity==pyProfile)
						{
							f(i,j,k,npy) = add*f(i,j,k,npy) + dens;
						}
						if (profile[s]->whichQuantity==pzProfile)
						{
							f(i,j,k,npz) = add*f(i,j,k,npz) + dens;
						}
					}
		}
	}
	f.ApplyBoundaryCondition(Element(ns,Xi));
	return didGenerate;
}

void Chemical::ReadInputFileTerm(std::stringstream& inputString,std::string& command)
{
	tw::Int i;
	std::string word;
	UnitConverter uc(owner->unitDensityCGS);
	Module::ReadInputFileTerm(inputString,command);
	ionization.ReadInputFileTerm(inputString,command);
	if (command=="ion") // eg, ion species = N3
	{
		inputString >> word >> word >> word;
		for (i=0;i<chemBoss->chemical.size();i++)
			if (chemBoss->chemical[i]->name==word)
				ionization.ionSpecies = i;
	}
	if (command=="electron") // eg, electron species = electrons
	{
		inputString >> word >> word >> word;
		for (i=0;i<chemBoss->chemical.size();i++)
			if (chemBoss->chemical[i]->name==word)
				ionization.electronSpecies = i;
	}
	if (command=="charge") // eg, charge = 1.0
	{
		inputString >> word;
		inputString >> charge;
	}
	if (command=="mass") // eg, mass = 1800.0
	{
		inputString >> word;
		inputString >> mass;
	}
	if (command=="effective") // eg, effective mass = 1.8 for density >> 1.0
	{
		inputString >> word >> word >> effectiveMass >> word >> word >> word >> transitionDensity;
	}
	if (command=="permittivity") // eg, permittivity = 5.0 , 0.0
	{
		tw::Float tempr,tempi;
		inputString >> word;
		inputString >> tempr >> tempi;
		permittivity = tw::Complex(tempr,tempi);
	}
	if (command=="thermometric") // eg, thermometric conductivity = 1.0 // cm^2/s
	{
		inputString >> word >> word >> thermometricConductivity;
		thermometricConductivity *= uc.CGSValue(time_dim)/sqr(uc.CGSValue(length_dim));
	}
	if (command=="kinematic") // eg, kinematic viscosity = 1.0 // cm^2/s
	{
		inputString >> word >> word >> kinematicViscosity;
		kinematicViscosity *= uc.CGSValue(time_dim)/sqr(uc.CGSValue(length_dim));
	}
	if (command=="vibrational" || command=="vibration" || command=="excitation") // eg, vibrational energy = 0.3
	{
		inputString >> word >> word >> excitationEnergy;
		excitationEnergy = uc.eV_to_sim(excitationEnergy);
	}
	if (command=="rotational") // eg, rotational degrees of freedom = 2
	{
		inputString >> word >> word >> word >> word >> cvm;
		cvm = 0.5*(3.0 + cvm);
	}
	if (command=="cv") // eg, cv = 2.5 (this is really m*cv, or m*cv/kB)
	{
		inputString >> word >> cvm;
	}
	// ASHER_MOD
	// EOS Model to use
	if (command=="EOS") {
		inputString >> word >> word;
		eosName = word;
	}
	// dimensionless Gruneisen parameter (if Gruneisen Model is used)
	if (command=="GRUN") {
		inputString >> word;
		inputString >> GRUN;
	}

	// Reference Density
	if (command=="n0") {
		inputString >> word;
		inputString >> n0;
	}

	// y - intercept of Hugoniot fit (usually appriximately speed of sound)
	if (command=="c0") {
		inputString >> word;
		inputString >> c0;
	}

	// coefficient of linear fit of Hugoniot data
	if (command=="S1") {
		inputString >> word;
		inputString >> S1;
	}
}

void Chemical::ReadData(std::ifstream& inFile)
{
	Module::ReadData(inFile);
	inFile.read((char *)&charge,sizeof(tw::Float));
	inFile.read((char *)&mass,sizeof(tw::Float));
	inFile.read((char *)&permittivity,sizeof(tw::Complex));
	inFile.read((char *)&effectiveMass,sizeof(tw::Float));
	inFile.read((char *)&transitionDensity,sizeof(tw::Float));
	inFile.read((char *)&excitationEnergy,sizeof(tw::Float));
	inFile.read((char *)&cvm,sizeof(tw::Float));
	inFile.read((char *)&thermometricConductivity,sizeof(tw::Float));
	inFile.read((char *)&kinematicViscosity,sizeof(tw::Float));
	inFile.read((char *)&ionization,sizeof(IonizationData));
}

void Chemical::WriteData(std::ofstream& outFile)
{
	Module::WriteData(outFile);
	outFile.write((char *)&charge,sizeof(tw::Float));
	outFile.write((char *)&mass,sizeof(tw::Float));
	outFile.write((char *)&permittivity,sizeof(tw::Complex));
	outFile.write((char *)&effectiveMass,sizeof(tw::Float));
	outFile.write((char *)&transitionDensity,sizeof(tw::Float));
	outFile.write((char *)&excitationEnergy,sizeof(tw::Float));
	outFile.write((char *)&cvm,sizeof(tw::Float));
	outFile.write((char *)&thermometricConductivity,sizeof(tw::Float));
	outFile.write((char *)&kinematicViscosity,sizeof(tw::Float));
	outFile.write((char *)&ionization,sizeof(IonizationData));
}

void Chemical::BoxDiagnosticHeader(GridDataDescriptor* box)
{
	owner->WriteBoxDataHeader(name,box);
}

void Chemical::BoxDiagnose(GridDataDescriptor* box)
{
	CopyFieldData(chemBoss->scratch,Element(0),chemBoss->state1,Element(indexInState));
	chemBoss->scratch *= chemBoss->fluxMask;
	owner->WriteBoxData(name,box,&chemBoss->scratch(0,0,0),chemBoss->scratch.Stride());
}

void Chemical::PointDiagnosticHeader(std::ofstream& outFile)
{
	outFile << name << " ";
}

void Chemical::PointDiagnose(std::ofstream& outFile,const weights_3D& w)
{
	std::valarray<tw::Float> densNow(chemBoss->state1.Components());
	chemBoss->state1.Interpolate(densNow,Element(indexInState),w);
	outFile << densNow[indexInState] << " ";
}

///////////////////
//               //
//   CHEMISTRY   //
//               //
///////////////////


Chemistry::Chemistry(Grid* theGrid):Module(theGrid)
{
	name = "Chemistry";
	typeCode = chemistry;
	epsilonFactor = 1e-4;
	tolerance = 1e-6;
	overrelaxation = 1.9;
	maxIterations = 10000;
	laserFrequency = 1.0;
	phi_lbc = phi_rbc = 0.0;

	// some arrays must be initialized later since we don't know how many elements they have yet
	scratch.Initialize(*this,owner);
	scratch2.Initialize(*this,owner);
	fluxMask.Initialize(*this,owner);
	rho0.Initialize(*this,owner);
	rho.Initialize(*this,owner);
	phi.Initialize(*this,owner);
	nu_e.Initialize(*this,owner);
	me_eff.Initialize(*this,owner);
	laserAmplitude.Initialize(*this,owner);
	radiativeLosses.Initialize(*this,owner);
	radiationIntensity.Initialize(*this,owner);
	refractiveIndex.Initialize(*this,owner);
	gres.resize(owner->globalCells[3]+2);

	parabolicSolver = (ParabolicSolver*)owner->AddPrivateTool(generalParabolicPropagator);
	if (owner->Dimensionality()==1)
		ellipticSolver = (EllipticSolver*)owner->AddPrivateTool(ellipticSolver1D);
	else
		ellipticSolver = (EllipticSolver*)owner->AddPrivateTool(iterativePoissonSolver);
	laserPropagator = (IsotropicPropagator*)owner->AddPrivateTool(isotropicPropagator);
	radModel = sparc::noRadiation;
	lasModel = sparc::vacuum;
	plasModel = sparc::neutral;
	electrons = NULL;
}

Chemistry::~Chemistry()
{
	tw::Int i;
	owner->RemoveTool(parabolicSolver);
	owner->RemoveTool(ellipticSolver);
	owner->RemoveTool(laserPropagator);
	for (i=0;i<reaction.size();i++)
		delete reaction[i];
	for (i=0;i<excitation.size();i++)
		delete excitation[i];
	for (i=0;i<collision.size();i++)
		delete collision[i];
}

void Chemistry::Initialize()
{
	tw::Int i,j,k,s,totalElements;

	Module::Initialize();

	ellipticSolver->tolerance = tolerance;
	ellipticSolver->overrelaxation = overrelaxation;
	ellipticSolver->maxIterations = maxIterations;
	ellipticSolver->lbc = phi_lbc;
	ellipticSolver->rbc = phi_rbc;
	for (i=0;i<owner->conductor.size();i++)
		ellipticSolver->FixPotential(phi,owner->conductor[i]->theRgn,0.0);

	// Initialize laser frequency and refractive index
	if (owner->pulse.size())
	{
		laserFrequency = 0.0;
		for (s=0;s<owner->pulse.size();s++)
			laserFrequency += owner->pulse[s]->w;
		laserFrequency /= tw::Float(owner->pulse.size());
	}
	refractiveIndex = tw::Complex(1.0,0.0);

	// Electrostatic boundary conditions
	ellipticSolver->SetBoundaryConditions(phi,owner->bc0[1]==axisymmetric ? neumannWall : dirichletCell,dirichletCell,dirichletCell,dirichletCell,dirichletCell,neumannWall);
	rho.SetBoundaryConditions(xAxis,neumannWall,neumannWall);
	rho.SetBoundaryConditions(yAxis,neumannWall,neumannWall);
	rho.SetBoundaryConditions(zAxis,neumannWall,neumannWall);
	scratch.SetBoundaryConditions(xAxis,neumannWall,neumannWall);
	scratch.SetBoundaryConditions(yAxis,neumannWall,neumannWall);
	scratch.SetBoundaryConditions(zAxis,neumannWall,neumannWall);

	// Find number of state variables characterizing the whole system
	// For each group we have:
	// density1,density2,...,npx,npy,npz,U,Xi
	// The 5 components per group of the EOS array (T,Tv,P,K,visc) are not counted
	totalElements = 0;
	for (i=0;i<group.size();i++)
		totalElements += group[i]->chemical.size() + 6; // ASHER_MOD -- now that I'm including Cv in the EOS

	if (!owner->restarted)
	{
		state0.Initialize(totalElements,*this,owner);
		state1.Initialize(totalElements,*this,owner);
		creationRate.Initialize(totalElements,*this,owner);
		destructionRate.Initialize(totalElements,*this,owner);
		eos0.Initialize(group.size()*6,*this,owner); // ASHER_MOD -- now that I'm including Cv in the EOS
		eos1.Initialize(group.size()*6,*this,owner); // there are 6 elements inside these (5->6)
	}

	// variables defining normal component boundary conditions
	boundarySpec bc0[4],bc1[4];
	for (tw::Int ax=1;ax<=3;ax++)
	{
		bc0[ax] = (owner->bc0[ax] == reflecting || owner->bc0[ax] == axisymmetric) ? dirichletWall : neumannWall;
		bc1[ax] = (owner->bc1[ax] == reflecting || owner->bc1[ax] == axisymmetric) ? dirichletWall : neumannWall;
	}

	// setup the flux mask used to make conductors impermeable
	// also used for reflecting boundary conditions at simulation walls
	fluxMask = 1.0;
	for (CellIterator cell(*this,true);cell<cell.end();++cell)
		for (s=0;s<owner->conductor.size();s++)
			if (owner->conductor[s]->theRgn->Inside(owner->Pos(cell),*owner))
				fluxMask(cell) = 0.0;
	for (tw::Int ax=1;ax<=3;ax++)
		for (StripIterator strip(*this,ax,strongbool::yes);strip<strip.end();++strip)
		{
			if (bc0[ax]==dirichletWall && owner->n0[ax]==MPI_PROC_NULL && dim[ax]>1)
				fluxMask(strip,lb[ax]) = fluxMask(strip,0) = 0.0;
			if (bc1[ax]==dirichletWall && owner->n1[ax]==MPI_PROC_NULL && dim[ax]>1)
				fluxMask(strip,ub[ax]) = fluxMask(strip,dim[ax]+1) = 0.0;
		}

	// vector elements of state are modified further in species loop below
	state0.SetBoundaryConditions(xAxis,neumannWall,neumannWall);
	state0.SetBoundaryConditions(yAxis,neumannWall,neumannWall);
	state0.SetBoundaryConditions(zAxis,neumannWall,neumannWall);
	state1.SetBoundaryConditions(xAxis,neumannWall,neumannWall);
	state1.SetBoundaryConditions(yAxis,neumannWall,neumannWall);
	state1.SetBoundaryConditions(zAxis,neumannWall,neumannWall);
	eos0.SetBoundaryConditions(xAxis,neumannWall,neumannWall);
	eos0.SetBoundaryConditions(yAxis,neumannWall,neumannWall);
	eos0.SetBoundaryConditions(zAxis,neumannWall,neumannWall);
	eos1.SetBoundaryConditions(xAxis,neumannWall,neumannWall);
	eos1.SetBoundaryConditions(yAxis,neumannWall,neumannWall);
	eos1.SetBoundaryConditions(zAxis,neumannWall,neumannWall);
	creationRate.SetBoundaryConditions(xAxis,neumannWall,neumannWall);
	creationRate.SetBoundaryConditions(yAxis,neumannWall,neumannWall);
	creationRate.SetBoundaryConditions(zAxis,neumannWall,neumannWall);
	destructionRate.SetBoundaryConditions(xAxis,neumannWall,neumannWall);
	destructionRate.SetBoundaryConditions(yAxis,neumannWall,neumannWall);
	destructionRate.SetBoundaryConditions(zAxis,neumannWall,neumannWall);
	nu_e.SetBoundaryConditions(xAxis,neumannWall,neumannWall);
	nu_e.SetBoundaryConditions(yAxis,neumannWall,neumannWall);
	nu_e.SetBoundaryConditions(zAxis,neumannWall,neumannWall);

	// Set up the indices that point to the conglomerate field elements
	// Also fix the normal component boundary conditions
	j = 0; // running index within EOS array
	k = 0; // running index within state array
	for (i=0;i<group.size();i++)
	{
		group[i]->T = j++;
		group[i]->Tv = j++;
		group[i]->P = j++;
		group[i]->K = j++;
		group[i]->visc = j++;

		group[i]->Cv  = j++; // ASHER_MOD -- let Cv have its own index

		group[i]->e.low = k;
		group[i]->e.high = k + group[i]->chemical.size() - 1;
		k += group[i]->chemical.size();
		for (s=0;s<group[i]->chemical.size();s++)
		{
			group[i]->chemical[s]->indexInGroup = s;
			group[i]->chemical[s]->indexInState = group[i]->e.low + s;
		}

		// X-Component
		state0.SetBoundaryConditions(Element(k),xAxis,bc0[1],bc1[1]);
		state1.SetBoundaryConditions(Element(k),xAxis,bc0[1],bc1[1]);
		creationRate.SetBoundaryConditions(Element(k),xAxis,bc0[1],bc1[1]);
		destructionRate.SetBoundaryConditions(Element(k),xAxis,bc0[1],bc1[1]);
		group[i]->npx = k++;
		// Y-Component
		state0.SetBoundaryConditions(Element(k),yAxis,bc0[2],bc1[2]);
		state1.SetBoundaryConditions(Element(k),yAxis,bc0[2],bc1[2]);
		creationRate.SetBoundaryConditions(Element(k),yAxis,bc0[2],bc1[2]);
		destructionRate.SetBoundaryConditions(Element(k),yAxis,bc0[2],bc1[2]);
		group[i]->npy = k++;
		// Z-Component
		state0.SetBoundaryConditions(Element(k),zAxis,bc0[3],bc1[3]);
		state1.SetBoundaryConditions(Element(k),zAxis,bc0[3],bc1[3]);
		creationRate.SetBoundaryConditions(Element(k),zAxis,bc0[3],bc1[3]);
		destructionRate.SetBoundaryConditions(Element(k),zAxis,bc0[3],bc1[3]);
		group[i]->npz = k++;

		group[i]->U = k++;
		group[i]->Xi = k++;
	}

	// find electrons
	for (i=0;i<chemical.size();i++)
		if (chemical[i]->mass==1.0)
		{
			electrons = chemical[i];
			ie = chemical[i]->indexInState;
			electrons->group->forceFilter = 0.0;
		}
	me_eff = 1.0;
	nu_e = 1.0;
}

void Chemistry::Reset()
{
	tw::Int i,c;
	bool didGenerate = false;
	for (i=0;i<group.size();i++)
		for (c=0;c<group[i]->chemical.size();c++)
			didGenerate += group[i]->chemical[c]->GenerateFluid(state1,false);
	if (didGenerate || owner->IsFirstStep())
		ApplyEOS(state1,eos1);
	state0 = state1;
	eos0 = eos1;
	rho0 = rho;
}

tw::Float Chemistry::CoulombCrossSection(tw::Float q1,tw::Float q2,tw::Float m12,tw::Float v12,tw::Float N1,tw::Float N2,tw::Float T1,tw::Float T2) const
{
	UnitConverter uc(owner->unitDensityCGS);
	tw::Float rmin,rmax,rmin_alt,coulombLog,hbar;
	hbar = uc.hbar / (uc.MKSValue(time_dim) * uc.MKSValue(energy_dim));
	rmin = fabs(q1*q2)/(tw::small_pos + 4.0*pi*m12*v12*v12*uc.MKSValue(number_dim));
	rmin_alt = hbar/(tw::small_pos + 2.0*m12*v12);
	if (rmin_alt<rmin) rmin = rmin_alt;
	rmax = 1.0/sqrt(tw::small_pos + N1*q1*q1/T1 + N2*q2*q2/T2);
	coulombLog = log(rmax/rmin);
	if (coulombLog<1.0) coulombLog = 1.0;
	return (32.0/pi)*pow(v12,tw::Float(-4.0))*sqr(q1*q2/(4*pi*m12))*coulombLog/uc.MKSValue(number_dim);
}

tw::Float Chemistry::ElectronPhononRateCoeff(tw::Float Ti,tw::Float EFermi,tw::Float ks,tw::Float nref) const
{
	// here, the rate coefficient is collision frequency divided by a reference density
	// this allows us to multiply away the collisions in regions where the metallic density is at "background level"
	UnitConverter uc(owner->unitDensityCGS);
	// get quantities into cgs
	tw::Float vF = uc.c*100.0*sqrt(2.0*EFermi);
	tw::Float qe = uc.SimToCGS(charge_dim,1.0);
	tw::Float kB_Ti = uc.SimToCGS(energy_dim,Ti);
	tw::Float hbar = 1e7 * uc.hbar;
	// collision frequency in real units
	tw::Float nu = 2.0*ks*qe*qe*kB_Ti/(sqr(hbar)*vF);
	// return normalized rate coefficient
	return (uc.CGSValue(time_dim)/nref) * nu;
}

void Chemistry::ComputeElectronCollisionFrequency()
{
	#pragma omp parallel
	{
		tw::Float T1,T2,m1,m2,m12,v12,q1,q2,N1,N2,Ni,Ti,sigma;
		tw::Float aggregateCollFreq,nlattice,nratio,weight,weightSum;
		tw::Float coulombCollFreq,phononCollFreq;
		Chemical *chem1,*chem2;
		for (CellIterator cell(*this,false);cell<cell.end();++cell)
		{
			aggregateCollFreq = 0.0;

			// Compute effective mass in this cell
			me_eff(cell) = 0.0;
			weightSum = 0.0;
			for (tw::Int s=0;s<chemical.size();s++)
			{
				nlattice = state1(cell,chemical[s]->indexInState);
				nratio = nlattice / chemical[s]->transitionDensity;
				weight = chemical[s]->mass*nlattice;
				weightSum += weight;
				me_eff(cell) += weight*(1.0 + nratio*chemical[s]->effectiveMass)/(1.0 + nratio);
			}
			me_eff(cell) /= weightSum;

			// Compute collision frequency
			for (tw::Int s=0;s<collision.size();s++)
			{
				chem1 = chemical[collision[s]->chem1];
				chem2 = chemical[collision[s]->chem2];
				if (chem1==electrons || chem2==electrons)
				{
					N1 = state1(cell,chem1->indexInState);
					N2 = state1(cell,chem2->indexInState);
					T1 = eos1(cell,chem1->group->T);
					T2 = eos1(cell,chem2->group->T);
					m1 = chem1->mass;
					m2 = chem2->mass;
					q1 = chem1->charge;
					q2 = chem2->charge;
					m12 = m1*m2/(m1+m2);
					v12 = sqrt(8.0*(T1/m1 + T2/m2)/pi);
					Ni = (chem1==electrons ? N2 : N1);
					Ti = (chem1==electrons ? T2 : T1);

					if (collision[s]->type==sparc::hard_sphere)
					{
						sigma = collision[s]->crossSection;
						aggregateCollFreq += (4.0/3.0) * sigma * Ni * v12;
					}

					if (collision[s]->type==sparc::coulomb)
					{
						sigma = CoulombCrossSection(q1,q2,m12,v12,N1,N2,T1,T2);
						aggregateCollFreq += (4.0/3.0) * sigma * Ni * v12;
					}

					if (collision[s]->type==sparc::metallic)
					{
						phononCollFreq = Ni * ElectronPhononRateCoeff(Ti,collision[s]->T_ref,collision[s]->ks,collision[s]->n_ref);
						coulombCollFreq = (4.0/3.0) * CoulombCrossSection(q1,q2,m12,v12,N1,N2,T1,T2) * Ni * v12;
						aggregateCollFreq += coulombCollFreq*phononCollFreq / (coulombCollFreq + phononCollFreq);
					}
				}
			}
			nu_e(cell) = aggregateCollFreq;
		}
	}
}

void Chemistry::ComputeCollisionalSources()
{
	#pragma omp parallel
	{
		tw::Int r,s,c,i1,i2;
		tw::vec3 forceDensity,forceNow;
		tw::Float powerDensity,powerNow,vibrationalDensity,vibrationsNow;

		tw::Float logT,T1,T2,m1,m2,m12,v12,q1,q2,N1,N2,sigma;
		tw::Float rate1,rate2,rateNow,weight;
		tw::vec3 v1_v2;
		Chemical *chem,*chem1,*chem2;

		tw::Float Te,Tv,Ti,energy,level,n0;
		tw::Float Xv;

		UnitConverter uc(owner->unitDensityCGS);

		// REACTIONS

		for (r=0;r<reaction.size();r++)
			for (CellIterator cell(*this,false);cell<cell.end();++cell)
			{
				// Work out the primitive reaction rate (number density per unit time)
				Te = eos1(cell,chemical[reaction[r]->catalyst]->group->T);
				if (Te<reaction[r]->T0 || Te>reaction[r]->T1)
					rateNow = 0.0;
				else
				{
					if (reaction[r]->c1==0.0)
					{
						logT = log(tw::small_pos + uc.sim_to_eV(fabs(Te)));
						rateNow = 0.0;
						for (s=0;s<9;s++)
							rateNow += reaction[r]->b[s]*pow(logT,tw::Float(s));
						rateNow = exp(rateNow) * uc.CGSValue(time_dim) * std::pow(uc.CGSValue(density_dim),tw::Float(reaction[r]->numBodies));
					}
					else
						rateNow = reaction[r]->c1 * pow(Te,reaction[r]->c2) * exp(-reaction[r]->c3/Te);
				}
				for (s=0;s<reaction[r]->sub.size();s++)
					for (c=0;c<reaction[r]->sub[s]->reactant.size();c++)
						rateNow *= state1(cell,chemical[reaction[r]->sub[s]->reactant[c]]->indexInState);

				// Update creation and destruction arrays
				if (rateNow>0.0)
				{
					for (s=0;s<reaction[r]->sub.size();s++)
					{
						powerDensity = rateNow*(reaction[r]->sub[s]->heat + reaction[r]->sub[s]->vheat);
						vibrationalDensity = rateNow*(reaction[r]->sub[s]->vheat);
						forceDensity = 0.0;
						for (c=0;c<reaction[r]->sub[s]->reactant.size();c++)
						{
							chem = chemical[reaction[r]->sub[s]->reactant[c]];
							powerNow = rateNow*state1(cell,chem->group->U)/(tw::small_pos + chem->group->DensitySum(state1,cell));
							forceNow = rateNow*chem->mass*chem->group->Velocity(state1,cell);
							vibrationsNow = rateNow*state1(cell,chem->group->Xi)/(tw::small_pos + chem->group->ConditionalDensitySum(state1,chem->group->excitationEnergy,cell));

							DestroyMass(cell,chem->indexInGroup,rateNow,chem->group);
							DestroyTotalEnergy(cell,powerNow,chem->group);
							DestroyMomentum(cell,1,forceNow.x,chem->group);
							DestroyMomentum(cell,2,forceNow.y,chem->group);
							DestroyMomentum(cell,3,forceNow.z,chem->group);
							DestroyVibrations(cell,vibrationsNow,chem->group);

							powerDensity += powerNow;
							forceDensity += forceNow;
							vibrationalDensity += vibrationsNow;
						}
						// parcel out conserved quantities weighted by particle number
						// loses vibrational energy if any products atomic (total is still conserved)
						weight = 1.0/tw::Float(reaction[r]->sub[s]->product.size());
						for (c=0;c<reaction[r]->sub[s]->product.size();c++)
						{
							chem = chemical[reaction[r]->sub[s]->product[c]];
							CreateMass(cell,chem->indexInGroup,rateNow,chem->group);
							CreateTotalEnergy(cell,weight*powerDensity,chem->group);
							CreateMomentum(cell,1,weight*forceDensity.x,chem->group);
							CreateMomentum(cell,2,weight*forceDensity.y,chem->group);
							CreateMomentum(cell,3,weight*forceDensity.z,chem->group);
							if (chem->excitationEnergy>0.0)
								CreateVibrations(cell,weight*vibrationalDensity,chem->group);
						}
					}
				}
			}

		// COLLISIONS (MOMENTUM TRANSFER, THERMAL ENERGY TRANSFER)

		for (s=0;s<collision.size();s++)
			for (CellIterator cell(*this,false);cell<cell.end();++cell)
			{
				chem1 = chemical[collision[s]->chem1];
				chem2 = chemical[collision[s]->chem2];
				N1 = state1(cell,chem1->indexInState);
				N2 = state1(cell,chem2->indexInState);
				T1 = eos1(cell,chem1->group->T);
				T2 = eos1(cell,chem2->group->T);
				m1 = chem1->mass;
				m2 = chem2->mass;
				v1_v2 = chem1->group->Velocity(state1,cell) - chem2->group->Velocity(state1,cell);
				q1 = chem1->charge;
				q2 = chem2->charge;
				m12 = m1*m2/(m1+m2);
				v12 = sqrt(8*(T1/m1 + T2/m2)/pi);
				Ti = chem1==electrons ? T2 : T1;

				if (collision[s]->type==sparc::hard_sphere)
				{
					sigma = collision[s]->crossSection;
					rateNow = 4.0 * sigma * N1 * N2 * v12 * m12;
				}

				if (collision[s]->type==sparc::coulomb)
				{
					sigma = CoulombCrossSection(q1,q2,m12,v12,N1,N2,T1,T2);
					rateNow = 4.0 * sigma * N1 * N2 * v12 * m12;
				}

				if (collision[s]->type==sparc::metallic)
				{
					rate1 = 4.0*CoulombCrossSection(q1,q2,m12,v12,N1,N2,T1,T2)*v12;
					rate2 = ElectronPhononRateCoeff(Ti,collision[s]->T_ref,collision[s]->ks,collision[s]->n_ref);
					rateNow = N1*N2*m12/(1.0/rate1 + 1.0/rate2);
		 		}

				for (tw::Int ax=1;ax<=3;ax++)
				{
					CreateMomentum( cell , ax , rateNow * v1_v2[ax] / 3.0 , chem2->group);
					DestroyMomentum( cell , ax , rateNow * v1_v2[ax] / 3.0 , chem1->group);
				}
				if (T1>T2)
				{
					DestroyTotalEnergy( cell , rateNow * (T1 - T2) / (m1 + m2) , chem1->group);
					CreateTotalEnergy( cell ,rateNow * (T1 - T2) / (m1 + m2) , chem2->group);
				}
				else
				{
					CreateTotalEnergy( cell , rateNow * (T2 - T1) / (m1 + m2) , chem1->group);
					DestroyTotalEnergy( cell ,rateNow * (T2 - T1) / (m1 + m2) , chem2->group);
				}
			}

		// VIBRATIONS

		for (s=0;s<excitation.size();s++)
			for (CellIterator cell(*this,false);cell<cell.end();++cell)
			{
				chem1 = chemical[excitation[s]->exciter];
				chem2 = chemical[excitation[s]->excitee];
				i1 = chem1->indexInState;
				i2 = chem2->indexInState;
				Te = eos1(cell,chem1->group->T);
				Tv = eos1(cell,chem2->group->Tv);
				energy = chem2->excitationEnergy;
				level = excitation[s]->level;

				// compute reaction rate for this excitation
				if (excitation[s]->c1==0.0)
				{
					logT = log(tw::small_pos + uc.sim_to_eV(fabs(Te)));
					Xv = 0.0;
					for (c=0;c<9;c++)
						Xv += excitation[s]->b[c]*pow(logT,tw::Float(c));
					Xv = uc.CGSToSim(rate_coefficient_2_dim,exp(Xv));
				}
				else
					Xv = excitation[s]->c1 * pow(Te,excitation[s]->c2) * exp(-excitation[s]->c3/Te);

				// compute rate of change of vibrational energy
				if (level>0)
				{
					n0 = state1(cell,i2) * (1.0 - exp(-energy/Tv));
					rateNow = energy * level * Xv * state1(cell,i1) * n0 * (1.0 - exp(energy*level/Te - energy*level/Tv));
				}
				else
				{
					rateNow = energy * Xv * state1(cell,i1) * state1(cell,i2) * (1.0 - exp(energy/Te - energy/Tv));
				}

				CreateTotalAndVibrational(cell,rateNow,chem2->group);
				DestroyTotalEnergy(cell,rateNow,chem1->group);
			}
		}
}

void Chemistry::ComputeRadiativeSources()
{
	#pragma omp parallel
	{
		tw::Int s;
		Chemical *s1,*s2;
		tw::Float photoRate,Emag;
		UnitConverter uc(owner->unitDensityCGS);

		// Ohmic heating due to laser fields
		if (electrons)
			for (CellIterator cell(*this,false);cell<cell.end();++cell)
				CreateTotalEnergy(cell,0.5*nu_e(cell)*state1(cell,ie)*norm(laserAmplitude(cell))/(sqr(laserFrequency) + sqr(nu_e(cell))),electrons->group);

		// Photoionization
		for (CellIterator cell(*this,false);cell<cell.end();++cell)
			if (radiationIntensity(cell)>0.0)
				for (s=0;s<chemical.size();s++)
				{
					IonizationData& ionization = chemical[s]->ionization;
					if (ionization.ionizationModel!=noIonization)
					{
						Emag = sqrt(norm(laserAmplitude(cell)));
						s1 = chemical[ionization.ionSpecies];
						s2 = chemical[ionization.electronSpecies];
						photoRate = state1(cell,chemical[s]->indexInState)*ionization.Rate(0.0,Emag);
						CreateMass(cell,s1->indexInGroup,photoRate,s1->group);
						CreateMass(cell,s2->indexInGroup,photoRate,s2->group);
						DestroyMass(cell,chemical[s]->indexInGroup,photoRate,chemical[s]->group);
					}
				}

		// Compute radiative losses
		// Assume optically thin, estimate mean free path from Zel'dovich table 5.2
		// Strictly the formula only works in LTE, but when it is significant we probably are in LTE anyway
		// Hence we form equilibrium temperature from (total pressure / total density)
		if (radModel==sparc::thin)
		{
			for (CellIterator cell(*this,false);cell<cell.end();++cell)
			{
				tw::Float stef_boltz = 5.67e-8; // W/m^2/K^4
				tw::Float Ptot=0.0,ntot=0.0,TK,lossNow,meanFreePath;
				for (s=0;s<group.size();s++)
				{
					Ptot += eos1(cell,group[s]->P);
					ntot += group[s]->DensitySum(state1,cell);
				}
				TK = uc.SimToMKS(temperature_dim,Ptot/(tw::small_pos + ntot));
				meanFreePath = 8.0e-14 * sqr(TK); // m
				radiativeLosses(cell) = uc.MKSToSim(power_density_dim,4.0*stef_boltz*pow(TK,tw::Float(4.0))/(tw::small_pos + meanFreePath));
				for (s=0;s<group.size();s++)
				{
					lossNow = eos1(cell,group[s]->P)*radiativeLosses(cell)/(tw::small_pos + Ptot);
					DestroyTotalEnergy(cell,lossNow,group[s]);
				}
			}
		}
	}
}

tw::vec3 Chemistry::ComputeForceOnBody(tw::Int i,tw::Int j,tw::Int k)
{
	tw::vec3 ans;
	tw::Int s;
	EquilibriumGroup* g;

	for (s=0;s<group.size();s++)
	{
		g = group[s];
		const tw::Float Ax0 = owner->dS(i,j,k,1);
		const tw::Float Ax1 = owner->dS(i+1,j,k,1);
		const tw::Float Ay0 = owner->dS(i,j,k,2);
		const tw::Float Ay1 = owner->dS(i,j+1,k,2);
		const tw::Float Az0 = owner->dS(i,j,k,3);
		const tw::Float Az1 = owner->dS(i,j,k+1,3);

		const tw::Float Pc = eos1(i,j,k,g->P);

		if (fluxMask(i-1,j,k)==0.0 && fluxMask(i,j,k)==1.0)
			ans.x -= Pc*Ax0;

		if (fluxMask(i,j,k)==1.0 && fluxMask(i+1,j,k)==0.0)
			ans.x += Pc*Ax1;

		if (fluxMask(i,j-1,k)==0.0 && fluxMask(i,j,k)==1.0)
			ans.y -= Pc*Ay0;

		if (fluxMask(i,j,k)==1.0 && fluxMask(i,j+1,k)==0.0)
			ans.y += Pc*Ay1;

		if (fluxMask(i,j,k-1)==0.0 && fluxMask(i,j,k)==1.0)
			ans.z -= Pc*Az0;

		if (fluxMask(i,j,k)==1.0 && fluxMask(i,j,k+1)==0.0)
			ans.z += Pc*Az1;
	}
	return ans;
}

void Chemistry::ComputeHydroSources()
{
	#pragma omp parallel
	{
		EquilibriumGroup* g;

		for (tw::Int gidx=0;gidx<group.size();gidx++)
		{
			g = group[gidx];
			if (g->chemical[0]==electrons)
			{
				// electrons are not advanced using hydro
				// charge and current are computed implicitly in field advance
				// momentum density is fixed to give the heavy particle velocity (see ApplyEOS)
				// the velocity is then used to convect the energy density
				// mass and momentum density are convected, but then reset to restore quasineutrality and heavy particle velocity

				// take away any collisional sources of momentum
				for (CellIterator cell(*this,false);cell<cell.end();++cell)
				{
					creationRate(cell,g->npx) = 0.0;
					creationRate(cell,g->npy) = 0.0;
					creationRate(cell,g->npz) = 0.0;
					destructionRate(cell,g->npx) = 0.0;
					destructionRate(cell,g->npy) = 0.0;
					destructionRate(cell,g->npz) = 0.0;
				}
			}
			else
			{
				for (tw::Int ax=1;ax<=3;ax++)
				{
					#pragma omp barrier
					g->LoadVelocity(scratch,state1,ax);
					#pragma omp barrier
					for (CellIterator cell(*this,false);cell<cell.end();++cell)
					{
						tw::Float dV,dS0,dS1,dl0,dl1,P0,P1,v0,v1;
						tw::Float forceDensity = 0.0;
						tw::Float powerDensity = 0.0;
						owner->GetCellMetrics(cell,ax,&dV,&dS0,&dS1,&dl0,&dl1);

						const tw::Float E1 = (phi.bak(cell,0,ax) - phi.fwd(cell,0,ax)) / (dl0 + dl1);
						const tw::Float nm = g->DensityWeightedSum(state1,g->mass,cell);
						const tw::Float nq = g->DensityWeightedSum(state1,g->charge,cell);
						const tw::Float Pc = eos1(cell,g->P);
						const tw::Float vc = scratch(cell);
						const tw::Float f0 = fluxMask.bak(cell,0,ax);
						const tw::Float fc = fluxMask(cell);
						const tw::Float f1 = fluxMask.fwd(cell,0,ax);

						P0 = eos1.bak(cell,g->P,ax);
						P1 = eos1.fwd(cell,g->P,ax);
						v0 = scratch.bak(cell,0,ax);
						v1 = scratch.fwd(cell,0,ax);

						P0 = (P0*f0 + Pc*fc) / (f0 + fc + tw::small_pos);
						P1 = (P1*f1 + Pc*fc) / (f1 + fc + tw::small_pos);
						v0 = 0.5*f0*fc*(v0 + vc);
						v1 = 0.5*f1*fc*(v1 + vc);

						// external forces

						forceDensity += nq * E1;
						powerDensity += nq * (vc * E1);

						// pressure force

						forceDensity += (dS0*P0 - dS1*P1)/dV;

						// work done by pressure

						powerDensity += (v0*P0*dS0 - v1*P1*dS1)/dV;

						// N3 force from walls

						if (f0==0.0 && fc==1.0 && vc<0.0)
							forceDensity += nm*vc*vc/dl0;

						if (fc==1.0 && f1==0.0 && vc>0.0)
							forceDensity -= nm*vc*vc/dl1;

						// put everything into the sources

						CreateMomentum(cell,ax,fc*forceDensity,g);
						if (powerDensity>0.0)
							CreateTotalEnergy(cell,fc*powerDensity,g);
						else
							DestroyTotalEnergy(cell,-fc*powerDensity,g);
					}
				} // end loop over axes

				// Undifferentiated tensor divergence terms
				//#pragma omp barrier
				if (owner->gridGeometry==cylindrical)
					for (CellIterator cell(*this,false);cell<cell.end();++cell)
					{
						const tw::Float nm = g->DensityWeightedSum(state1,g->mass,cell);
						const tw::Float Pc = eos1(cell,g->P);
						const tw::vec3 vc = g->Velocity(state1,cell);
						const tw::vec3 pos = owner->Pos(cell);
						CreateMomentum(cell,1,fluxMask(cell)*(nm*sqr(vc.y) + Pc)/pos.x,g);
						DestroyMomentum(cell,2,fluxMask(cell)*nm*vc.x*vc.y/pos.x,g);
					}
				if (owner->gridGeometry==spherical)
					for (CellIterator cell(*this,false);cell<cell.end();++cell)
					{
						const tw::Float nm = g->DensityWeightedSum(state1,g->mass,cell);
						const tw::Float Pc = eos1(cell,g->P);
						const tw::vec3 vc = g->Velocity(state1,cell);
						const tw::vec3 pos = owner->Pos(cell);
						const tw::Float tanz = tan(pos.z);
						CreateMomentum(cell,1,fluxMask(cell)*(nm*(sqr(vc.y) + sqr(vc.z)) + 2.0*Pc)/pos.x,g);
						DestroyMomentum(cell,2,fluxMask(cell)*(nm*(vc.x*vc.z - vc.y*vc.y/tanz) - Pc/tanz)/pos.x,g);
						DestroyMomentum(cell,3,fluxMask(cell)*nm*(vc.x*vc.y + vc.y*vc.z/tanz)/pos.x,g);
					}

			} // end else
		} // end loop over groups
	} // end parallel region
}

void Chemistry::ComputeSources()
{
	creationRate = 0.0;
	destructionRate = 0.0;

	ComputeCollisionalSources();
	ComputeRadiativeSources();
	ComputeHydroSources();

	creationRate.CopyFromNeighbors();
	creationRate.ApplyBoundaryCondition();

	destructionRate.CopyFromNeighbors();
	destructionRate.ApplyBoundaryCondition();
}

void Chemistry::LaserAdvance(tw::Float dt)
{
	// Currently SPARC ignores differences in frequency between injected pulses.
	// The average frequency of all the pulses is used throughout (see also Chemistry::Initialize).

	if (owner->pulse.size() && lasModel==sparc::vacuum) // use a prescribed field
	{
		#pragma omp parallel
		{
			for (CellIterator cell(*this,false);cell<cell.end();++cell)
			{
				laserAmplitude(cell) = 0.0;
				for (tw::Int s=0;s<owner->pulse.size();s++)
					laserAmplitude(cell) += laserFrequency*owner->pulse[s]->VectorPotentialEnvelope(owner->elapsedTime,owner->Pos(cell));
				radiationIntensity(cell) = 0.5*norm(laserAmplitude(cell));
			}
		}
	}
	if (electrons && owner->pulse.size() && lasModel==sparc::isotropic && dim[3]>1)
	{
		#pragma omp parallel
		{
			tw::Complex a0,a1; // amplitudes of incoming waves on left
			tw::Complex aN,aN1; // amplitudes of incoming waves on right
			for (StripIterator strip(*this,3,strongbool::no);strip<strip.end();++strip)
			{
				a0 = a1 = aN = aN1 = 0.0;
				for (tw::Int s=0;s<owner->pulse.size();s++)
				{
					if (owner->pulse[s]->direction.z > 0.0)
					{
						a0 += std::exp(ii*laserFrequency*owner->Pos(strip,0).z)*laserFrequency*owner->pulse[s]->VectorPotentialEnvelope(owner->elapsedTime,owner->Pos(strip,0));
						a1 += std::exp(ii*laserFrequency*owner->Pos(strip,1).z)*laserFrequency*owner->pulse[s]->VectorPotentialEnvelope(owner->elapsedTime,owner->Pos(strip,1));
					}
					else
					{
						aN += std::exp(ii*laserFrequency*owner->Pos(strip,0).z)*laserFrequency*owner->pulse[s]->VectorPotentialEnvelope(owner->elapsedTime,owner->Pos(strip,dim[3]));
						aN1 += std::exp(ii*laserFrequency*owner->Pos(strip,1).z)*laserFrequency*owner->pulse[s]->VectorPotentialEnvelope(owner->elapsedTime,owner->Pos(strip,dim[3]+1));
					}
				}
				laserPropagator->SetupIncomingWaveLeft(strip,laserAmplitude,a0,a1,laserFrequency);
				laserPropagator->SetupIncomingWaveRight(strip,laserAmplitude,aN,aN1,laserFrequency);

				for (tw::Int k=1;k<=dim[3];k++)
				{
					refractiveIndex(strip,k) = 0.0;
					// Add dispersionless part of susceptibility
					for (tw::Int s=0;s<chemical.size();s++)
						refractiveIndex(strip,k) += state1(strip,k,chemical[s]->indexInState) * (chemical[s]->permittivity-one);
					// Add plasma contribution to susceptibility
					refractiveIndex(strip,k) -= state1(strip,k,ie)/sqr(laserFrequency)*(one - ii*nu_e(strip,k)/laserFrequency)/(one + sqr(nu_e(strip,k)/laserFrequency));
					// Convert susceptibility to refractive index
					refractiveIndex(strip,k) = std::sqrt(one + refractiveIndex(strip,k));
					// Following may be used to populate laser amplitude with vacuum field due to first pulse
					// In such case, laserPropagator->Advance should be commented out
					//laserAmplitude(strip,k) = laserFrequency*owner->pulse[0]->VectorPotentialEnvelope(owner->elapsedTime,owner->Pos(strip,k));
				}
			}
		}

		laserPropagator->Advance(laserAmplitude,refractiveIndex,nu_e,laserFrequency,dt);

		#pragma omp parallel
		{
			for (CellIterator cell(*this,false);cell<cell.end();++cell)
				radiationIntensity(cell) = real(refractiveIndex(cell))*0.5*norm(laserAmplitude(cell));
		}
	}
}

tw::Float Chemistry::EstimateTimeStep()
{
	std::vector<tw::Float> dtMax(tw::GetOMPMaxThreads());
	tw::Float dtMaxAllThreads,dtMaxAllNodes;

	#pragma omp parallel
	{
		const tw::Int tid = tw::GetOMPThreadNum();
		const tw::Float sqrt_eps = sqrt(epsilonFactor);
		const tw::Float creationDominance = 100.0;
		tw::Float dts;

		if (owner->IsFirstStep())
			dtMax[tid] = owner->dt0;
		else
			dtMax[tid] = owner->dtMax;

		// Courant condition

		for (CellIterator cell(*this,false);cell<cell.end();++cell)
			for (tw::Int s=0;s<group.size();s++)
			{
				const tw::vec3 vel = group[s]->Velocity(state1,cell);
				if (Norm(vel)!=0.0)
				{
					const tw::Float dxi2 = sqr(0.5*owner->dL(cell,1));
					const tw::Float deta2 = sqr(0.5*owner->dL(cell,2));
					const tw::Float dzeta2 = sqr(0.5*owner->dL(cell,3));
					dts = sqr(vel.x)/dxi2 + sqr(vel.y)/deta2 + sqr(vel.z)/dzeta2;
					dts = 0.9/sqrt(dts);
					if (dts < dtMax[tid])
					{
						dtMax[tid] = dts;
						//statusMessage.str("");
						//statusMessage << "Limited by Courant condition : " << group[s]->name << " v = " << Magnitude(vel) << std::endl;
					}
				}
			}

		// Reaction rates, collision rates, other source terms

		for (CellIterator cell(*this,false);cell<cell.end();++cell)
			for (tw::Int s=0;s<chemical.size();s++)
			{
				const tw::Int c = chemical[s]->indexInState;
				if (creationRate(cell,c) > creationDominance*destructionRate(cell,c))
					dts = sqrt_eps*fabs( (tw::small_pos+state1(cell,c)) / destructionRate(cell,c) );
				else
					dts = sqrt_eps*fabs( (tw::small_pos+state1(cell,c)) / (creationRate(cell,c) - destructionRate(cell,c)) );
				if (dts < dtMax[tid])
				{
					dtMax[tid] = dts;
					//statusMessage.str("");
					//statusMessage << "Limited by dN/dt : " << chemical[s]->name << " : N=" << state1(cell,c) << " , dN/dt=" << creationRate(cell,c) << "-" << destructionRate(cell,c) << std::endl;
				}
			}

		for (CellIterator cell(*this,false);cell<cell.end();++cell)
			for (tw::Int s=0;s<group.size();s++)
			{
				const tw::Int c = group[s]->U;
				if (creationRate(cell,c) > creationDominance*destructionRate(cell,c))
					dts = sqrt_eps*fabs( (tw::small_pos+state1(cell,c)) / destructionRate(cell,c) );
				else
					dts = sqrt_eps*fabs( (tw::small_pos+state1(cell,c)) / (creationRate(cell,c) - destructionRate(cell,c)) );
				if (dts < dtMax[tid])
				{
					dtMax[tid] = dts;
					//statusMessage.str("");
					//statusMessage << "Limited by dU/dt : " << group[s]->name << " : U=" << state1(cell,c) << " , dU/dt=" << creationRate(cell,c) << "-" << destructionRate(cell,c) << std::endl;
				}
			}

		for (CellIterator cell(*this,false);cell<cell.end();++cell)
			for (tw::Int s=0;s<group.size();s++)
			{
				const tw::Int c = group[s]->Xi;
				if (creationRate(cell,c) > creationDominance*destructionRate(cell,c))
					dts = sqrt_eps*fabs( (tw::small_pos+state1(cell,c)) / destructionRate(cell,c) );
				else
					dts = sqrt_eps*fabs( (tw::small_pos+state1(cell,c)) / (creationRate(cell,c) - destructionRate(cell,c)) );
				if (dts < dtMax[tid])
				{
					dtMax[tid] = dts;
					//statusMessage.str("");
					//statusMessage << "Limited by dX/dt : " << group[s]->name << " : X=" << state1(cell,c) << " , dX/dt=" << creationRate(cell,c) << "-" << destructionRate(cell,c) << std::endl;
				}
			}
	}

	// Choose the smallest maximum step from all the threads
	dtMaxAllThreads = *std::min_element(std::begin(dtMax),std::end(dtMax));
	// Choose the smallest maximum step from all the nodes
	dtMaxAllNodes = owner->strip[0].GetMin(dtMaxAllThreads);
	// Don't let it fall below the minimum time step
	return dtMaxAllNodes < owner->dtMin ? owner->dtMin : dtMaxAllNodes;
}

void Chemistry::DiffusionAdvance(tw::Float dt)
{
	tw::Int s,ax;
	EquilibriumGroup *g;

	for (s=0;s<group.size();s++)
	{
		g = group[s];

		// HEAT CONDUCTION

		g->LoadMassDensityCv(scratch,state1);
		g->LoadMassDensity(scratch2,state1);
		parabolicSolver->Advance(eos1,g->T,fluxMask,&scratch,0,&eos1,g->K,dt);

		#pragma omp parallel
		{
			for (CellIterator cell(*this,false);cell<cell.end();++cell)
			{
				state1(cell,g->U) = 0.5*scratch2(cell)*Norm(g->Velocity(state1,cell));
				state1(cell,g->U) += scratch(cell)*eos1(cell,g->T);
			}
		}

		// VISCOSITY

		for (ax=1;ax<=3;ax++)
		{
			g->LoadVelocity(scratch,state1,ax);
			CopyBoundaryConditions(scratch,0,state1,g->npx+ax-1);
			parabolicSolver->Advance(scratch,0,fluxMask,&scratch2,0,&eos1,g->visc,dt);
			#pragma omp parallel
			{
				for (CellIterator cell(*this,false);cell<cell.end();++cell)
					state1(cell,g->npx+ax-1) = scratch2(cell)*scratch(cell);
			}
		}
	}

	state1.CopyFromNeighbors();
	state1.ApplyBoundaryCondition();
	eos1.CopyFromNeighbors();
	eos1.ApplyBoundaryCondition();
}

void Chemistry::FieldAdvance(tw::Float dt)
{
	if (!electrons || plasModel==sparc::neutral)
		return;

	const tw::Float q0 = electrons->charge;
	const tw::Float m0 = electrons->mass;
	const tw::Int Pe = electrons->group->P;

	// set up source and coefficient
	#pragma omp parallel
	{
		tw::Float D1,D2,P0,P1,mu0,mu1,dV,dS0,dS1,dl0,dl1;
		for (CellIterator cell(*this,false);cell<cell.end();++cell)
		{
			rho(cell) = rho0(cell);
			for (tw::Int ax=1;ax<=3;ax++)
			{
				owner->GetCellMetrics(cell,ax,&dV,&dS0,&dS1,&dl0,&dl1);
				P0 = eos1.bak(cell,Pe,ax);
				P1 = eos1.fwd(cell,Pe,ax);
				mu0 = 2.0*(q0/m0)/(nu_e(cell)+nu_e.bak(cell,0,ax));
				mu1 = 2.0*(q0/m0)/(nu_e(cell)+nu_e.fwd(cell,0,ax));
				// Laplacian coefficients weighted by mobility
				// ( for calculation of div(mu*grad(P)) )
				D1 = mu0*dS0/(dl0*dV);
				D2 = mu1*dS1/(dl1*dV);
				// re-use rho to hold rho_eff  = rho(t=0) + dt*div(mu*grad(P))
				rho(cell) += dt*(D1*P0 - (D1+D2)*eos1(cell,Pe) + D2*P1);
			}
			// coefficients multiplying grad(phi) go into scratch array
			scratch(cell) = 1.0 + dt*q0*state1(cell,ie)*(q0/m0)/nu_e(cell);
		}
	}
	scratch.CopyFromNeighbors();
	scratch.ApplyBoundaryCondition();

	// Solve the elliptical equation --- div(scratch*grad(phi)) = -rho_eff
	ellipticSolver->SetCoefficients(&scratch);
	ellipticSolver->Solve(phi,rho,-1.0);

	// Compute the new charge density
	#pragma omp parallel for collapse(3) schedule(static)
	for (tw::Int i=1;i<=dim[1];i++)
		for (tw::Int j=1;j<=dim[2];j++)
			for (tw::Int k=1;k<=dim[3];k++)
			{
				rho(i,j,k) = (phi(i,j,k) - phi(i-1,j,k))*owner->dS(i,j,k,1) / owner->dl(i,j,k,1);
				rho(i,j,k) += (phi(i,j,k) - phi(i+1,j,k))*owner->dS(i+1,j,k,1) / owner->dl(i+1,j,k,1);
				rho(i,j,k) += (phi(i,j,k) - phi(i,j-1,k))*owner->dS(i,j,k,2) / owner->dl(i,j,k,2);
				rho(i,j,k) += (phi(i,j,k) - phi(i,j+1,k))*owner->dS(i,j+1,k,2) / owner->dl(i,j+1,k,2);
				rho(i,j,k) += (phi(i,j,k) - phi(i,j,k-1))*owner->dS(i,j,k,3) / owner->dl(i,j,k,3);
				rho(i,j,k) += (phi(i,j,k) - phi(i,j,k+1))*owner->dS(i,j,k+1,3) / owner->dl(i,j,k+1,3);
				rho(i,j,k) /= owner->dS(i,j,k,0);
			}
	rho.CopyFromNeighbors();
	rho.ApplyBoundaryCondition();
}

void Chemistry::HydroAdvance(const axisSpec& axis,tw::Float dt)
{
	// Convect the fluid

	tw::Int ax = naxis(axis);

	if (dim[ax] > 1)
	{
		FCT_Driver convector(&state0,&state1,&scratch,&fluxMask,owner);
		boundarySpec bc0 = (owner->bc0[ax] == reflecting || owner->bc0[ax] == axisymmetric) ? dirichletCell : neumannWall;
		boundarySpec bc1 = owner->bc1[ax] == reflecting ? dirichletCell : neumannWall;

		for (tw::Int c=0;c<group.size();c++)
		{
			convector.SetDensityElements(Element(group[c]->e.low,group[c]->e.high+6));
			convector.SetVelocityElement(0);
			#pragma omp parallel
			{
				group[c]->LoadVelocity(scratch,state1,ax);
				for (CellIterator cell(*this,true);cell<cell.end();++cell)
					scratch(cell) *= fluxMask(cell);
			}
			convector.Convect(axis,bc0,bc1,dt);
		}
		state0.ApplyBoundaryCondition();
	}
}

void Chemistry::ChemAdvance(tw::Float dt)
{
	creationRate *= dt;
	state0 += creationRate;

	destructionRate *= dt;
	state0 -= destructionRate;

	state0.CopyFromNeighbors();
	state0.ApplyBoundaryCondition();
}

void Chemistry::ApplyEOS(Field& hydro,Field& eos)
{
	// Load (P,T) into eos using (n,np,u) from hydro
	// Reverse calculation is applied following diffusion steps
	// Here we also compute transport coefficients and impose assumed conditions such as quasineutrality
	// Finally, we throw an error if numerical failure is detected

	// ASHER_MOD -- loop over groups and use eos object functions
	for (tw::Int c=0;c<group.size();c++)
	{
		if (group[c]->chemical[0]!=electrons)
		{
//			((EOSIdealGasMix *)(group[c]->eosMixData))->ApplyEOS( *group[c],group[c]->mass,group[c]->charge,group[c]->cvm,
			group[c]->eosMixData->ApplyEOS( group[c]->mass,group[c]->charge,group[c]->cvm,
											group[c]->excitationEnergy,group[c]->thermo_cond_cvm,group[c]->k_visc_m,
											group[c]->e,group[c]->T,group[c]->Tv,group[c]->P,group[c]->K,
											group[c]->visc,group[c]->Cv,group[c]->npx,group[c]->npy,group[c]->npz,
											group[c]->U,group[c]->Xi,ie,nu_e,hydro,eos);
		}
	}


	#pragma omp parallel
	{
		tw::Float ntot,ne,ionChargeDensity,nm,nmcv,epsvn,nv,KE,IE;
		tw::vec3 ionVelocity;

		for (CellIterator cell(*this,false);cell<cell.end();++cell)
		{
			ionChargeDensity = 0.0;
			ionVelocity = 0.0;
			for (tw::Int c=0;c<group.size();c++)
			{
				if (group[c]->chemical[0]!=electrons)
				{
					// ASHER_MOD : Comments for original EOS calculations here for reference

					// ntot = group[c]->DensitySum(hydro,cell);
					// nm = group[c]->DensityWeightedSum(hydro,group[c]->mass,cell);
					// nmcv = group[c]->DensityWeightedSum(hydro,group[c]->cvm,cell);
					// epsvn = group[c]->DensityWeightedSum(hydro,group[c]->excitationEnergy,cell);
					// nv = group[c]->ConditionalDensitySum(hydro,group[c]->excitationEnergy,cell);
					// KE = 0.5*nm*Norm(group[c]->Velocity(hydro,cell));
					// IE = hydro(cell,group[c]->U) - KE;
					// if (IE<=0.0)
					// {
					// 	hydro(cell,group[c]->U) = KE*1.00001 + tw::small_pos;
					// 	IE = hydro(cell,group[c]->U) - KE;
					// }
					// eos(cell,group[c]->T) = IE/(tw::small_pos + nmcv);
					// eos(cell,group[c]->Tv) = (epsvn/(nv+tw::small_pos))/log(1.0001 + epsvn/(hydro(cell,group[c]->Xi)+tw::small_pos));
					// eos(cell,group[c]->P) = ntot * eos(cell,group[c]->T);
					// eos(cell,group[c]->K) = group[c]->DensityWeightedSum(hydro,group[c]->thermo_cond_cvm,cell); // thermometricConductivity * nmcv;
					// eos(cell,group[c]->visc) = group[c]->DensityWeightedSum(hydro,group[c]->k_visc_m,cell); // kinematicViscosity * nm;

					ionChargeDensity += group[c]->DensityWeightedSum(hydro,group[c]->charge,cell);
					ionVelocity += group[c]->Velocity(hydro,cell);
				}
			}
			if (electrons)
			{
				// Impose quasineutrality and assume electron velocity = average heavy particle velocity
				ne = -ionChargeDensity/electrons->charge;
				hydro(cell,ie) = ne;
				hydro(cell,electrons->group->npx) = ne*electrons->mass*ionVelocity.x/tw::Float(group.size()-1);
				hydro(cell,electrons->group->npy) = ne*electrons->mass*ionVelocity.y/tw::Float(group.size()-1);
				hydro(cell,electrons->group->npz) = ne*electrons->mass*ionVelocity.z/tw::Float(group.size()-1);

				// // EOS for electrons
				// nmcv = 1.5*ne;
				// eos(cell,electrons->group->T) = hydro(cell,electrons->group->U)/(tw::small_pos + nmcv);
				// eos(cell,electrons->group->Tv) = 0.0;
				// eos(cell,electrons->group->P) = ne*eos(cell,electrons->group->T);
				// eos(cell,electrons->group->K) = 3.2*ne*eos(cell,electrons->group->T)/(electrons->mass*nu_e(cell));
				// eos(cell,electrons->group->visc) = 0.0;
				// // Braginskii has for e-viscosity 0.73*ne*eos(cell,electrons->group->T)/nu_e(cell)
				// // However, we are forcing electrons to move with ions and so should not diffuse velocity field
			}
		}
	}


	// ASHER_MOD -- electron EOSs have the be calculated after the quasineutrality-ensuring hydro advance
	for (tw::Int c=0;c<group.size();c++)
	{
		if (group[c]->chemical[0]==electrons)
		{
			// so far only does an ideal gas case
			electrons->eosData->ApplyEOS(   electrons->mass,electrons->charge,electrons->cvm,				
											electrons->excitationEnergy,electrons->thermometricConductivity,electrons->kinematicViscosity,
											electrons->group->e,electrons->group->T,electrons->group->Tv,electrons->group->P,electrons->group->K,
											electrons->group->visc,group[c]->Cv,electrons->group->npx,electrons->group->npy,electrons->group->npz,
											electrons->group->U,electrons->group->Xi,ie,nu_e,hydro,eos);
		}
	}


	// Check for numerical failure, defined by NaN in the hydro state vector
	tw::Int badCells = 0;
	for (CellIterator cell(*this,true);cell<cell.end();++cell)
		for (tw::Int c=0;c<hydro.Components();c++)
			badCells += std::isnan(hydro(cell,c));
	if (badCells)
		throw tw::FatalError("Encountered NaN in hydrodynamic state");

	hydro.CopyFromNeighbors();
	hydro.ApplyBoundaryCondition();
	eos.CopyFromNeighbors();
	eos.ApplyBoundaryCondition();
}

void Chemistry::FirstOrderAdvance(tw::Float dt,bool computeSources)
{
	if (computeSources)
	{
		ComputeElectronCollisionFrequency();
		ComputeSources();
	}
	HydroAdvance(xAxis,dt);
	HydroAdvance(yAxis,dt);
	HydroAdvance(zAxis,dt);
	ChemAdvance(dt);
	ApplyEOS(state0,eos0);

	Swap(state0,state1);
	Swap(eos0,eos1);

	DiffusionAdvance(dt);
	FieldAdvance(dt);
}

void Chemistry::Update()
{
	tw::Float dts;

	if (!owner->adaptiveTimestep)
	{
		FirstOrderAdvance(0.5*dt,true);
		FirstOrderAdvance(dt,true);
	}
	else
	{
		ComputeElectronCollisionFrequency();
		ComputeSources();
		dts = EstimateTimeStep();
		owner->SetupTimeInfo(dts);
		FirstOrderAdvance(0.5*dt,false);
		FirstOrderAdvance(dt,true);
	}

	LaserAdvance(dt);
}

void Chemistry::ParseReaction(std::stringstream& inputString)
{
	tw::Int i,j;
	std::string word;
	tw::Float sign,energy;
	tw::Int lastChem = -1;
	bool rhs = false;
	bool foundSpecies = false;

	reaction.push_back(new Reaction);
	reaction.back()->numBodies = 0;
	inputString >> word >> word >> word; // take off "=" and "{"

	UnitConverter uc(owner->unitDensityCGS);

	// Read in reactants and products
	do
	{
		reaction.back()->sub.push_back(new SubReaction);
		reaction.back()->sub.back()->heat = 0.0;
		reaction.back()->sub.back()->vheat = 0.0;
		rhs = false;
		do
		{
			if (word=="->")
			{
				rhs = true;
				inputString >> word;
			}
			if (word=="+")
			{
				sign = 1.0;
				inputString >> word;
			}
			if (word=="-")
			{
				sign = -1.0;
				inputString >> word;
			}
			if (word==":")
				inputString >> word;

			// see if this is a heat of reaction
			// if no number is found, encode this result by setting energy=0
			// WARNING: this fails if a species name begins with inf or nan
			try { energy = uc.eV_to_sim(std::stod(word,NULL)); }
			catch (std::invalid_argument) { energy = 0.0; }
			if (energy!=0.0) // it is
			{
				if (lastChem==-1)
					throw tw::FatalError("heat of reaction read before chemical name");
				if (word.back()=='v')
					reaction.back()->sub.back()->vheat = energy*sign;
				else
					reaction.back()->sub.back()->heat = energy*sign;
			}
			else
			{
				foundSpecies = false;
				for (i=0;i<chemical.size();i++)
				{
					if (chemical[i]->name==word)
					{
						foundSpecies = true;
						lastChem = i;
						if (rhs)
						{
							reaction.back()->sub.back()->product.push_back(lastChem);
						}
						else
						{
							reaction.back()->sub.back()->reactant.push_back(lastChem);
							reaction.back()->numBodies++;
						}
					}
				}
				if (!foundSpecies)
				{
					(*owner->tw_out) << "ERROR: couldn't find chemical " << word << std::endl;
					abort();
				}
			}
			inputString >> word;
		} while (word!=":" && word!="}");
	} while (word!="}");

	// Get reaction rate data and normalize
	inputString >> word;
	if (word=="rate")
	{
		inputString >> word >> reaction.back()->c1 >> reaction.back()->c2 >> reaction.back()->c3;
		reaction.back()->c1 *= pow(uc.eV_to_sim(1.0),-reaction.back()->c2);
		reaction.back()->c1 *= uc.CGSValue(time_dim) * pow(uc.CGSValue(density_dim),tw::Float(reaction.back()->numBodies-1));
		reaction.back()->c3 = uc.eV_to_sim(reaction.back()->c3);
	}
	if (word=="janev_rate")
	{
		reaction.back()->c1 = 0.0;
		inputString >> word;
		for (i=0;i<9;i++)
			inputString >> reaction.back()->b[i];
		// leave janev rates in cm^3/s and eV
	}

	// Get temperature range and catalyst
	inputString >> word;
	foundSpecies = false;
	for (i=0;i<chemical.size();i++)
	{
		if (chemical[i]->name==word)
		{
			foundSpecies = true;
			reaction.back()->catalyst = i;
		}
	}
	if (!foundSpecies)
	{
		(*owner->tw_out) << "ERROR: couldn't find chemical " << word << std::endl;
		abort();
	}
	inputString >> word;
	tw::input::PythonRange(word,&reaction.back()->T0,&reaction.back()->T1);
	reaction.back()->T0 = uc.eV_to_sim(reaction.back()->T0);
	reaction.back()->T1 = uc.eV_to_sim(reaction.back()->T1);

	// Print message
	(*owner->tw_out) << "Reaction: ";
	for (i=0;i<reaction.back()->sub.size();i++)
		for (j=0;j<reaction.back()->sub[i]->reactant.size();j++)
			(*owner->tw_out) << chemical[reaction.back()->sub[i]->reactant[j]]->name << " ";
	(*owner->tw_out) << "-> ";
	for (i=0;i<reaction.back()->sub.size();i++)
		for (j=0;j<reaction.back()->sub[i]->product.size();j++)
			(*owner->tw_out) << chemical[reaction.back()->sub[i]->product[j]]->name << " ";
	(*owner->tw_out) << "   Rate = " << reaction.back()->c1;
	if (reaction.back()->c2!=0.0)
		(*owner->tw_out) << " * T(" << chemical[reaction.back()->catalyst]->name << ")^" << reaction.back()->c2;
	if (reaction.back()->c3!=0.0)
		(*owner->tw_out) << " * exp[-" << reaction.back()->c3 << "/T(" << chemical[reaction.back()->catalyst]->name << ")]";
	for (i=0;i<reaction.back()->sub.size();i++)
		(*owner->tw_out) << " : " << reaction.back()->sub[i]->heat << " " << reaction.back()->sub[i]->vheat << "v";
	(*owner->tw_out) << " : " << reaction.back()->T0 << " " << reaction.back()->T1;
	(*owner->tw_out) << std::endl;
}

void Chemistry::ParseExcitation(std::stringstream& inputString)
{
	tw::Int i;
	std::string word,species;
	UnitConverter uc(owner->unitDensityCGS);

	excitation.push_back(new Excitation);

	inputString >> word >> word; // take off "="
	for (i=0;i<chemical.size();i++)
		if (chemical[i]->name==word)
			excitation.back()->exciter = i;

	inputString >> word >> word;
	for (i=0;i<chemical.size();i++)
		if (chemical[i]->name==word)
			excitation.back()->excitee = i;

	inputString >> word >> word >> excitation.back()->level;

	inputString >> word;
	if (word=="rate")
	{
		inputString >> word >> excitation.back()->c1 >> excitation.back()->c2 >> excitation.back()->c3;
		excitation.back()->c1 = pow(uc.eV_to_sim(1.0),-excitation.back()->c2) * uc.CGSToSim(rate_coefficient_2_dim,excitation.back()->c1);
		excitation.back()->c3 = uc.eV_to_sim(excitation.back()->c3);
	}
	if (word=="janev_rate")
	{
		excitation.back()->c1 = 0.0;
		inputString >> word;
		for (i=0;i<9;i++)
			inputString >> excitation.back()->b[i];
		// leave janev rates in cm^3/s and eV
	}

	(*owner->tw_out) << "Excitation: " << chemical[excitation.back()->exciter]->name << " " << chemical[excitation.back()->excitee]->name << std::endl;
}

void Chemistry::ParseCollision(std::stringstream& inputString)
{
	// hard sphere example
	// new collision = e <-> N , cross section = 5.0

	// coulomb collision example
	// new collision = e <-> N[+] , coulomb

	// metallic collision example
	// new collision = e <-> Cu[+] , metallic , ks = 2.4 , fermi_energy_ev = 7.0 , ref_density = 3000

	tw::Int i;
	std::string word,species;
	UnitConverter uc(owner->unitDensityCGS);

	collision.push_back(new Collision);

	inputString >> word >> word; // take off "="
	for (i=0;i<chemical.size();i++)
		if (chemical[i]->name==word)
			collision.back()->chem1 = i;

	inputString >> word >> word;
	for (i=0;i<chemical.size();i++)
		if (chemical[i]->name==word)
			collision.back()->chem2 = i;

	inputString >> word;

	if (word=="cross")
	{
		collision.back()->type = sparc::hard_sphere;
		inputString >> word >> word >> collision.back()->crossSection;
		(*owner->tw_out) << "Hard sphere ";
	}

	if (word=="coulomb")
	{
		collision.back()->type = sparc::coulomb;
		(*owner->tw_out) << "Coulomb ";
	}

	if (word=="metallic")
	{
		collision.back()->type = sparc::metallic;
		inputString >> word >> word >> collision.back()->ks;
		inputString >> word >> word >> collision.back()->T_ref;
		inputString >> word >> word >> collision.back()->n_ref;
		collision.back()->T_ref *= uc.eV_to_sim(1.0);
		(*owner->tw_out) << "Metallic ";
	}

	(*owner->tw_out) << "collision: " << chemical[collision.back()->chem1]->name << " " << chemical[collision.back()->chem2]->name << std::endl;
}


void Chemistry::ReadInputFileTerm(std::stringstream& inputString,std::string& command)
{
	std::string word;

	Module::ReadInputFileTerm(inputString,command);
	if (command=="epsilon") // eg epsilon factor = 0.1
		inputString >> word >> word >> epsilonFactor;
	if (command=="radiation") // eg radiation model = thin
	{
		inputString >> word >> word >> word;
		if (word=="thin")
			radModel = sparc::thin;
		if (word=="thick")
			radModel = sparc::thick;
	}
	if (command=="laser") // eg laser propagator = vacuum
	{
		inputString >> word >> word >> word;
		if (word=="vacuum")
			lasModel = sparc::vacuum;
		if (word=="isotropic")
			lasModel = sparc::isotropic;
	}
	if (command=="plasma") // eg plasma model = neutral
	{
		inputString >> word >> word >> word;
		if (word=="neutral")
			plasModel = sparc::neutral;
		if (word=="quasineutral")
			plasModel = sparc::quasineutral;
	}
	if (command=="tolerance") // eg tolerance = 1e-6
		inputString >> word >> tolerance;
	if (command=="overrelaxation") // eg overrelaxation = 1.9
		inputString >> word >> overrelaxation;
	if (command=="iterations") // eg iterations = 1000
		inputString >> word >> maxIterations;
	if (command=="reaction")
		ParseReaction(inputString);
	if (command=="excitation")
		ParseExcitation(inputString);
	if (command=="collision")
		ParseCollision(inputString);
	if (command=="dipole") // eg, dipole center = 0 0 0
		inputString >> word >> word >> dipoleCenter.x >> dipoleCenter.y >> dipoleCenter.z;
	if (command=="external") // eg, external potential = ( 0.0 , 1.0 )
	{
		inputString >> word >> word >> phi_lbc >> phi_rbc;
	}
}

void Chemistry::ReadData(std::ifstream& inFile)
{
	tw::Int i,num,eos_size;
	Module::ReadData(inFile);
	inFile.read((char *)&radModel,sizeof(sparc::radiationModel));
	inFile.read((char *)&lasModel,sizeof(sparc::laserModel));
	inFile.read((char *)&plasModel,sizeof(sparc::plasmaModel));
	inFile.read((char *)&epsilonFactor,sizeof(tw::Float));
	inFile.read((char *)&tolerance,sizeof(tw::Float));
	inFile.read((char *)&overrelaxation,sizeof(tw::Float));
	inFile.read((char *)&maxIterations,sizeof(tw::Int));
	inFile.read((char *)&dipoleCenter,sizeof(tw::vec3));
	inFile.read((char *)&phi_lbc,sizeof(tw::Float));
	inFile.read((char *)&phi_rbc,sizeof(tw::Float));

	inFile.read((char *)&num,sizeof(tw::Int));
	inFile.read((char *)&eos_size,sizeof(tw::Int));

	state0.Initialize(num,*this,owner);
	state1.Initialize(num,*this,owner);
	creationRate.Initialize(num,*this,owner);
	destructionRate.Initialize(num,*this,owner);
	eos0.Initialize(eos_size,*this,owner);
	eos1.Initialize(eos_size,*this,owner);

	eos1.ReadData(inFile);
	state1.ReadData(inFile);
	rho.ReadData(inFile);
	phi.ReadData(inFile);

	inFile.read((char *)&num,sizeof(tw::Int));
	for (i=0;i<num;i++)
	{
		(*owner->tw_out) << "Add Reaction" << std::endl;
		reaction.push_back(new Reaction);
		reaction.back()->ReadData(inFile);
	}

	inFile.read((char *)&num,sizeof(tw::Int));
	for (i=0;i<num;i++)
	{
		(*owner->tw_out) << "Add Excitation" << std::endl;
		excitation.push_back(new Excitation);
		excitation.back()->ReadData(inFile);
	}

	inFile.read((char *)&num,sizeof(tw::Int));
	for (i=0;i<num;i++)
	{
		(*owner->tw_out) << "Add Collision" << std::endl;
		collision.push_back(new Collision);
		collision.back()->ReadData(inFile);
	}
}

void Chemistry::WriteData(std::ofstream& outFile)
{
	tw::Int i;
	Module::WriteData(outFile);
	outFile.write((char *)&radModel,sizeof(sparc::radiationModel));
	outFile.write((char *)&lasModel,sizeof(sparc::laserModel));
	outFile.write((char *)&plasModel,sizeof(sparc::plasmaModel));
	outFile.write((char *)&epsilonFactor,sizeof(tw::Float));
	outFile.write((char *)&tolerance,sizeof(tw::Float));
	outFile.write((char *)&overrelaxation,sizeof(tw::Float));
	outFile.write((char *)&maxIterations,sizeof(tw::Int));
	outFile.write((char *)&dipoleCenter,sizeof(tw::vec3));
	outFile.write((char *)&phi_lbc,sizeof(tw::Float));
	outFile.write((char *)&phi_rbc,sizeof(tw::Float));


	i = state1.Components();
	outFile.write((char *)&i,sizeof(tw::Int));
	i = group.size()*6; // ASHER_MOD -- there are 6 EOSs now...
	outFile.write((char *)&i,sizeof(tw::Int));

	eos1.WriteData(outFile);
	state1.WriteData(outFile);
	rho.WriteData(outFile);
	phi.WriteData(outFile);

	i = reaction.size();
	outFile.write((char *)&i,sizeof(tw::Int));
	for (i=0;i<reaction.size();i++)
		reaction[i]->WriteData(outFile);

	i = excitation.size();
	outFile.write((char *)&i,sizeof(tw::Int));
	for (i=0;i<excitation.size();i++)
		excitation[i]->WriteData(outFile);

	i = collision.size();
	outFile.write((char *)&i,sizeof(tw::Int));
	for (i=0;i<collision.size();i++)
		collision[i]->WriteData(outFile);
}

void Chemistry::EnergyHeadings(std::ofstream& outFile)
{
	outFile << "Mass Charge Energy Px Py Pz Fx Fy Fz Dx Dy Dz ";
}

void Chemistry::EnergyColumns(std::vector<tw::Float>& cols,std::vector<bool>& avg,const Region& theRgn)
{
	tw::Int i,j,k,s;
	CellIterator cell(*this,false);
	tw::Int x0,x1,y0,y1,z0,z1;
	tw::Float dV;
	tw::vec3 r1,r2,pos,bodyForce,fluidMomentum,dipoleMoment;
	tw::Float totalEnergy = 0.0;
	tw::Float totalMass = 0.0;
	tw::Float totalCharge = 0.0;
	theRgn.GetLocalCellBounds(&x0,&x1,&y0,&y1,&z0,&z1);
	for (k=z0;k<=z1;k++)
		for (j=y0;j<=y1;j++)
			for (i=x0;i<=x1;i++)
			{
				cell.SetCell(i,j,k);
				pos = owner->Pos(i,j,k);
				if (theRgn.Inside(pos,*owner))
				{
					dV = owner->dS(i,j,k,0);
					owner->CurvilinearToCartesian(&(r1=dipoleCenter));
					owner->CurvilinearToCartesian(&(r2=pos));
					for (s=0;s<group.size();s++)
					{
						totalMass += group[s]->DensityWeightedSum(state1,group[s]->mass,cell) * dV;
						totalEnergy += state1(i,j,k,group[s]->U) * dV;
						totalCharge += rho(i,j,k) * dV;
						fluidMomentum.x += state1(i,j,k,group[s]->npx) * dV;
						fluidMomentum.y += state1(i,j,k,group[s]->npy) * dV;
						fluidMomentum.z += state1(i,j,k,group[s]->npz) * dV;
						dipoleMoment += (r2 - r1) * rho(i,j,k) * dV;
					}
					bodyForce += ComputeForceOnBody(i,j,k);
				}
			}
	cols.push_back(totalMass); avg.push_back(false);
	cols.push_back(totalCharge); avg.push_back(false);
	cols.push_back(totalEnergy); avg.push_back(false);
	cols.push_back(fluidMomentum.x); avg.push_back(false);
	cols.push_back(fluidMomentum.y); avg.push_back(false);
	cols.push_back(fluidMomentum.z); avg.push_back(false);
	cols.push_back(bodyForce.x); avg.push_back(false);
	cols.push_back(bodyForce.y); avg.push_back(false);
	cols.push_back(bodyForce.z); avg.push_back(false);
	cols.push_back(dipoleMoment.x); avg.push_back(false);
	cols.push_back(dipoleMoment.y); avg.push_back(false);
	cols.push_back(dipoleMoment.z); avg.push_back(false);
}

void Chemistry::BoxDiagnosticHeader(GridDataDescriptor* box)
{
	owner->WriteBoxDataHeader("impermeable",box);
	owner->WriteBoxDataHeader("collisionFreq",box);
	owner->WriteBoxDataHeader("massdensity",box);
	owner->WriteBoxDataHeader("rad-intensity",box);
	owner->WriteBoxDataHeader("rad-ereal",box);
	owner->WriteBoxDataHeader("rad-eimag",box);
	owner->WriteBoxDataHeader("rad-losses",box);
	owner->WriteBoxDataHeader("rad-nreal",box);
	owner->WriteBoxDataHeader("rad-nimag",box);
	owner->WriteBoxDataHeader("me_eff",box);
	if (electrons)
	{
		/*owner->WriteBoxDataHeader("chem-Ex",box);
		owner->WriteBoxDataHeader("chem-Ey",box);
		owner->WriteBoxDataHeader("chem-Ez",box);*/
		owner->WriteBoxDataHeader("chem-phi",box);
		owner->WriteBoxDataHeader("chem-rho",box);
		owner->WriteBoxDataHeader("chem-jz",box);
	}
}

void Chemistry::BoxDiagnose(GridDataDescriptor* box)
{
	tw::Int i,j,k,s;

	// If first step we need to apply EOS so that EquilibriumGroups can write out T

	if (owner->IsFirstStep())
	{
		ApplyEOS(state1,eos1);
		ComputeElectronCollisionFrequency();
	}

	// Impermeable Region

	owner->WriteBoxData("impermeable",box,&fluxMask(0,0,0),fluxMask.Stride());

	// Collision Diagnostic

	owner->WriteBoxData("collisionFreq",box,&nu_e(0,0,0),nu_e.Stride());

	// Mass Density Diagnostic

	scratch = 0.0;
	for (s=0;s<group.size();s++)
	{
		group[s]->LoadMassDensity(scratch2,state1);
		scratch2 *= fluxMask;
		scratch += scratch2;
	}
	owner->WriteBoxData("massdensity",box,&scratch(0,0,0),scratch.Stride());

	// Radiation Diagnostic

	owner->WriteBoxData("rad-intensity",box,&radiationIntensity(0,0,0),radiationIntensity.Stride());
	owner->WriteBoxData("rad-ereal",box,&laserAmplitude(0,0,0,0),laserAmplitude.Stride());
	owner->WriteBoxData("rad-eimag",box,&laserAmplitude(0,0,0,1),laserAmplitude.Stride());
	owner->WriteBoxData("rad-losses",box,&radiativeLosses(0,0,0),radiativeLosses.Stride());
	owner->WriteBoxData("rad-nreal",box,&refractiveIndex(0,0,0,0),refractiveIndex.Stride());
	owner->WriteBoxData("rad-nimag",box,&refractiveIndex(0,0,0,1),refractiveIndex.Stride());
	owner->WriteBoxData("me_eff",box,&me_eff(0,0,0),me_eff.Stride());

	// Electrostatics

	if (electrons)
	{
		owner->WriteBoxData("chem-phi",box,&phi(0,0,0),phi.Stride());
		owner->WriteBoxData("chem-rho",box,&rho(0,0,0),rho.Stride());

		tw::Float q0 = electrons->charge;
		tw::Float m0 = electrons->mass;
		tw::Int Pe = electrons->group->P;
		for (k=1;k<=dim[3];k++)
			for (j=1;j<=dim[2];j++)
				for (i=1;i<=dim[1];i++)
				{
					scratch(i,j,k) = q0*q0*state1(i,j,k,ie)*(phi(i,j,k-1)-phi(i,j,k+1))/owner->dL(i,j,k,3);
					scratch(i,j,k) += q0*(eos1(i,j,k-1,Pe)-eos1(i,j,k+1,Pe))/owner->dL(i,j,k,3);
					scratch(i,j,k) /= m0*nu_e(i,j,k);
				}
		owner->WriteBoxData("chem-jz",box,&scratch(0,0,0),scratch.Stride());
	}
}

void Chemistry::CustomDiagnose()
{
	/*const tw::Int customPeriod = 10;

	tw::Int i;
	float data;
	LocalDomain *master;
	ofstream outFile;

	master = (*owner->global)(0,0,0);

	if (owner->IsFirstStep() && !(owner->appendMode&&owner->restarted) && owner->local==master)
	{
		outFile.open("resolution.dvdat",ios::binary);
		WriteDVHeader(outFile,2,1,1,owner->global->dim[3],0.0,1.0,0.0,1.0,0.0,owner->global->dim[3]);
		outFile.close();
	}

	if (owner->stepNow%customPeriod==0 && owner->local==master)
	{
		outFile.open("resolution.dvdat",ios::binary | ios::app);
		for (i=1;i<=owner->global->dim[3];i++)
		{
			data = gres[i];
			WriteBigEndian((char *)&data,sizeof(float),0,outFile);
		}
		outFile.close();
	}*/
}

void Chemistry::StatusMessage(std::ostream *theStream)
{
	*theStream << statusMessage.str();
}
