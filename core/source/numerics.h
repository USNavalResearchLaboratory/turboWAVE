void Transform(tw::Float * array,tw::Int pts,tw::Int modes,tw::Int interval,std::valarray<tw::Float>& transform);
void ReverseTransform(tw::Float * array,tw::Int pts,tw::Int modes,tw::Int interval,std::valarray<tw::Float>& rev_transform);
tw::Float GetSphericalGroundState(std::valarray<tw::Float>& vec,std::valarray<tw::Float>& phi,tw::Float dr);
tw::Float GetCylindricalGroundState(std::valarray<tw::Float>& vec,std::valarray<tw::Float>& phi,tw::Float dr);
void ComputeTransformMatrices(tw::bc::fld radial_bc,std::valarray<tw::Float>& eigenvalue,std::valarray<tw::Float>& fwd,std::valarray<tw::Float>& rev,MetricSpace *space,Task *tsk);

template <class T,class U>
void TriDiagonal(std::valarray<T>& phi,std::valarray<T>& rho,U a,U b,U c);
template <class T,class U>
void TriDiagonal(std::valarray<T>& phi,std::valarray<T>& rho,U a,std::valarray<U>& b,U c);
template <class T,class U>
void TriDiagonal(std::valarray<T>& phi,std::valarray<T>& rho,std::valarray<U>& a,std::valarray<U>& b,std::valarray<U>& c);

// operations on a general matrix and eigenvectors
void LeftHalfRotation(ScalarField& A,tw::Int p,tw::Int q,tw::Float c,tw::Float s,tw::Int def1,tw::Int def2);
void RightHalfRotation(ScalarField& A,tw::Int p,tw::Int q,tw::Float c,tw::Float s,tw::Int def1,tw::Int def2);
void QRSweep(ScalarField& A,tw::Int deflation);
void QLSweep(ScalarField& A,tw::Int deflation);
void GetEigenvector(tw::Float eigenvalue,std::valarray<tw::Float>& vec,std::valarray<tw::Float>& a,std::valarray<tw::Float>& b,std::valarray<tw::Float>& c);
void NormalizeLeftRight(std::valarray<tw::Float>& left,std::valarray<tw::Float>& right);

// symmetric matrices
void SymmetricTridiagonalEigensystem(std::valarray<tw::Float>& eigenvalues,tw::Float *revTransform,std::valarray<tw::Float>& T1,std::valarray<tw::Float>& T2);
void SortEigensystem(std::valarray<tw::Float>& eigenvalues,tw::Float *revTransform);

template <class T>
void ComputeAlphasAndBetas(tw::comm *strip,tw::Int systems,T* k);

template <class T>
class GlobalIntegrator
{
	// the domains and systems are indexed from 0
	// the data is expected to include ghost cells, so the number of pts per "row" is N+2

	tw::comm *strip;
	tw::Int N,systems;
	std::valarray<T> v,w;
	std::valarray<T*> theData;
	std::valarray<tw::Int> stride; // units of T (beware of passing strides from Field objects)

	public:
	GlobalIntegrator(tw::comm *strip,tw::Int systems,tw::Int cells);
	void SetData(tw::Int system,T* theData,tw::Int stride);
	void SetMatrix(tw::Int system,T a,T b,T c,T theta,T eta);
	void SetMatrix(tw::Int system,T a,std::valarray<T>& b,T c,T theta,T eta);
	void SetMatrix(tw::Int system,std::valarray<T>& a,std::valarray<T>& b,std::valarray<T>& c);
	void Parallelize();
};

struct FCT_Engine
{
	// cells is the number of grid cells, not counting the 2 ghost cells
	// The V and A arrays have to be set up to contain cell volumes and low-side cell wall areas
	tw::Int cells;
	std::valarray<tw::Float> V,A,scratch;

	FCT_Engine(tw::Int ax,const MetricSpace& m);
	void Reset(const tw::strip& s,const MetricSpace& m,ScalarField *fluxMask);
	void Transport(std::valarray<tw::Float>& vel,std::valarray<tw::Float>& rho,std::valarray<tw::Float>& rho1,std::valarray<tw::Float>& diff,std::valarray<tw::Float>& flux,tw::Float dt);
	void Diffuse(std::valarray<tw::Float>& vel,std::valarray<tw::Float>& rho,std::valarray<tw::Float>& diff,tw::Float dt);
	void Limiter(tw::Float& adiff,const tw::Float& maxLow,const tw::Float& maxHigh);
	void Clip(std::valarray<tw::Float>& rho,std::valarray<tw::Float>& adiff,tw::Float rho00);
	void AntiDiffuse(std::valarray<tw::Float>& rho,std::valarray<tw::Float>& adiff,std::valarray<tw::Float>& flux);
};

struct FCT_Driver
{
	Element en; // indices of density components
	tw::Int vi; // index to velocity component
	Field *diff,*net_flux;

	// externally owned data
	Field* rho;
	Field* rho1;
	Field* vel;
	ScalarField* fluxMask;
	MetricSpace* ms;

	FCT_Driver(Field *rho,Field *rho1,Field *vel,ScalarField *fluxMask,MetricSpace *ms);
	~FCT_Driver();
	void SetDensityElements(const Element& e) { en = e; }
	void SetVelocityElement(tw::Int v) { vi = v; }
	void Convect(const tw::grid::axis& axis,tw::bc::fld low,tw::bc::fld high,tw::Float dt);
	void GetTrueFlux(Field& flux,const Element& dst,const Element& src)
	{
		CopyFieldData(flux,dst,*net_flux,src);
	}
};


/////////////////
//             //
//  TEMPLATES  //
//             //
/////////////////


///////////////////
//               //
//  TRIDIAGONAL  //
//               //
///////////////////


template <class T,class U>
void TriDiagonal(std::valarray<T>& phi,std::valarray<T>& rho,U a,U b,U c)
{
	// invert A * phi = rho
	// The rows of A are (b,c,0,...)(a,b,c,0,...)(0,a,b,c,0,...)...(0,...,a,b,c)(0,...,a,b)

	tw::Int n = phi.size();
	std::valarray<U> gam(n);

	if (n==1)
	{
		phi[0] = rho[0]/b;
		return;
	}

	tw::Int i;
	U bet;
	bet = b;
	phi[0] = rho[0]/bet;
	for (i=1;i<=n-1;i++)
	{
		gam[i] = c/bet;
		bet = b - a*gam[i];
		phi[i] = (rho[i] - a*phi[i-1])/bet;
	}

	for (i=n-2;i>=0;i--)
		phi[i] -= gam[i+1]*phi[i+1];
}

template <class T,class U>
void TriDiagonal(std::valarray<T>& phi,std::valarray<T>& rho,U a,std::valarray<U>& b,U c)
{
	// invert A * phi = rho
	// The rows of A are (b[0],c,0,...)(a,b[1],c,0,...)(0,a,b[2],c,0,...)...(0,...,a,b[N-2],c)(0,...,a,b[N-1])

	tw::Int n = phi.size();
	std::valarray<U> gam(n);

	if (n==1)
	{
		phi[0] = rho[0]/b[0];
		return;
	}

	tw::Int i;
	U bet;
	bet = b[0];
	phi[0] = rho[0]/bet;
	for (i=1;i<=n-1;i++)
	{
		gam[i] = c/bet;
		bet = b[i] - a*gam[i];
		phi[i] = (rho[i] - a*phi[i-1])/bet;
	}

	for (i=n-2;i>=0;i--)
		phi[i] -= gam[i+1]*phi[i+1];
}

template <class T,class U>
void TriDiagonal(std::valarray<T>& phi,std::valarray<T>& rho,std::valarray<U>& a,std::valarray<U>& b,std::valarray<U>& c)
{
	// invert A * phi = rho
	// The rows of A are (b[0],c[0],0,...)(a[1],b[1],c[1],0,...)(0,a,b,c,0,...)...(0,...,a,b,c)(0,...,a,b)

	tw::Int n = phi.size();
	std::valarray<U> gam(n);

	if (n==1)
	{
		phi[0] = rho[0]/b[0];
		return;
	}

	tw::Int i;
	U bet;
	bet = b[0];
	phi[0] = rho[0]/bet;
	for (i=1;i<=n-1;i++)
	{
		gam[i] = c[i-1]/bet;
		bet = b[i] - a[i]*gam[i];
		phi[i] = (rho[i] - a[i]*phi[i-1])/bet;
	}

	for (i=n-2;i>=0;i--)
		phi[i] -= gam[i+1]*phi[i+1];
}


/////////////////////////
//                     //
//  GLOBAL INTEGRATOR  //
//                     //
/////////////////////////

template <class T>
GlobalIntegrator<T>::GlobalIntegrator(tw::comm *strip,tw::Int systems,tw::Int cells)
{
	this->strip = strip;
	this->systems = systems;
	N = cells;

	v.resize(systems*(N+2));
	w.resize(systems*(N+2));

	theData.resize(systems);
	stride.resize(systems);
}

template <class T>
void GlobalIntegrator<T>::SetData(tw::Int system,T* row,tw::Int stride)
{
	// N.b. the stride is in units of T
	// When calling from TW beware of Field strides in units of tw::Float
	// and DiscreteSpace strides in units of cells
	this->theData[system] = row;
	this->stride[system] = stride;
}

template <class T>
void GlobalIntegrator<T>::SetMatrix(tw::Int system,T a,T b,T c,T theta,T eta)
{
	std::valarray<T> u(N),basisVector(N);

	// Set up the inversion vectors for one system

	basisVector = T(0.0);
	basisVector[0] = theta;
	TriDiagonal<T,T>(u,basisVector,a,b,c);
	v[std::slice(system*(N+2)+1,N,1)] = u;

	basisVector[0] = 0.0;
	basisVector[N-1] = eta;
	TriDiagonal<T,T>(u,basisVector,a,b,c);
	w[std::slice(system*(N+2)+1,N,1)] = u;
}

template <class T>
void GlobalIntegrator<T>::SetMatrix(tw::Int system,T a,std::valarray<T>& b,T c,T theta,T eta)
{
	std::valarray<T> u(N),basisVector(N);

	// Set up the inversion vectors for one system

	basisVector = T(0.0);
	basisVector[0] = theta;
	TriDiagonal<T,T>(u,basisVector,a,b,c);
	v[std::slice(system*(N+2)+1,N,1)] = u;

	basisVector[0] = 0.0;
	basisVector[N-1] = eta;
	TriDiagonal<T,T>(u,basisVector,a,b,c);
	w[std::slice(system*(N+2)+1,N,1)] = u;
}

template <class T>
void GlobalIntegrator<T>::SetMatrix(tw::Int system,std::valarray<T>& a,std::valarray<T>& b,std::valarray<T>& c)
{
	std::valarray<T> u(N),basisVector(N);

	// Set up the inversion vectors for one system

	basisVector = T(0.0);
	basisVector[0] = a[0];
	TriDiagonal<T,T>(u,basisVector,a,b,c);
	v[std::slice(system*(N+2)+1,N,1)] = u;

	basisVector[0] = 0.0;
	basisVector[N-1] = c[N-1];
	TriDiagonal<T,T>(u,basisVector,a,b,c);
	w[std::slice(system*(N+2)+1,N,1)] = u;
}

template <class T>
void ComputeAlphasAndBetas(tw::comm *strip,tw::Int systems,T* mpi_packet)
{
	tw::Int i,j,index,dsn,ds=8*systems,ss=8;
	tw::Int domains = strip->Get_size();
	tw::Int L = strip->Get_size() - 1;
	tw::Int currDomain = strip->Get_rank();
	//bool periodic[2];
	//strip->Get_periods(1,periodic);
	tw::Int unknowns = L;//periodic[1] ? L+1 : L;
	std::valarray<T> source(unknowns),ans(unknowns),a(unknowns),b(unknowns),c(unknowns);
	std::valarray<T> k(ds*domains);

	strip->Gather(mpi_packet,&k[0],sizeof(T)*ds,0);

	for (j=0;j<systems;j++)
	{
		// APERIODIC CASE

		if (currDomain==0 && domains>1)
		//if (currDomain==0 && domains>1 && !periodic[1])
		{
			// For alphas domain 0 is not used
			// indices to tridiagonal arguments are therefore offset by -1
			for (i=1;i<L;i++)
			{
				index = ss*j + ds*i;
				source[i-1] = k[index-ds+5] - k[index-ds+2]*k[index] + k[index-ds+2]*k[index+4]*k[index+5]/(tw::small_pos + k[index+2]);
				a[i-1] = k[index-ds+1];
				b[i-1] = one - k[index-ds+2]*k[index+3] + k[index-ds+2]*k[index+1]*k[index+4]/(tw::small_pos + k[index+2]);
				c[i-1] = k[index-ds+2]*k[index+4]/(tw::small_pos + k[index+2]);
			}
			source[L-1] = k[ss*j+ds*(L-1)+5] - k[ss*j+ds*(L-1)+2]*k[ss*j+ds*L];
			a[L-1] = k[ss*j+ds*(L-1)+1];
			b[L-1] = one - k[ss*j+ds*(L-1)+2]*k[ss*j+ds*L+3];
			TriDiagonal<T,T>(ans,source,a,b,c);
			for (i=1;i<=L;i++)
				k[ss*j + ds*i + 6] = ans[i-1]; // alphas

			// For betas domain L is not used
			source[0] = k[ss*j+ds] - k[ss*j+ds+3]*k[ss*j+5];
			b[0] = one - k[ss*j+2]*k[ss*j+ds+3];
			c[0] = k[ss*j+ds+4];
			for (i=1;i<L;i++)
			{
				index = ss*j + ds*i;
				source[i] = k[index+ds] - k[index+ds+3]*k[index+5] + k[index+ds+3]*k[index+1]*k[index]/(tw::small_pos + k[index+3]);
				a[i] = k[index+ds+3]*k[index+1]/(tw::small_pos + k[index+3]);
				b[i] = one - k[index+2]*k[index+ds+3] + k[index+ds+3]*k[index+1]*k[index+4]/(tw::small_pos + k[index+3]);
				c[i] = k[index+ds+4];
			}
			TriDiagonal<T,T>(ans,source,a,b,c);
			for (i=0;i<L;i++)
				k[ss*j + ds*i + 7] = ans[i]; // betas
		}

		// PERIODIC CASE (doesn't work)

// 		if (currDomain==0 && periodic[1])
// 		{
// 			for (i=0;i<=L;i++)
// 			{
// 				index = ss*j + ds*i;
// 				dsn = i==0 ? -ds*L : ds;
// 				source[i] = k[index-dsn+5] - k[index-dsn+2]*k[index] + k[index-dsn+2]*k[index+4]*k[index+5]/(tw::small_pos + k[index+2]);
// 				a[i] = k[index-dsn+1];
// 				b[i] = one - k[index-dsn+2]*k[index+3] + k[index-dsn+2]*k[index+1]*k[index+4]/(tw::small_pos + k[index+2]);
// 				c[i] = k[index-dsn+2]*k[index+4]/(tw::small_pos + k[index+2]);
// 			}
// 			TriDiagonal<T,T>(ans,source,a,b,c);
// 			for (i=0;i<=L;i++)
// 				k[ss*j + ds*i + 6] = ans[i]; // alphas
//
// 			for (i=0;i<=L;i++)
// 			{
// 				index = ss*j + ds*i;
// 				dsn = i==L ? -ds*L : ds;
// 				source[i] = k[index+dsn] - k[index+dsn+3]*k[index+5] + k[index+dsn+3]*k[index+1]*k[index]/(tw::small_pos + k[index+3]);
// 				a[i] = k[index+dsn+3]*k[index+1]/(tw::small_pos + k[index+3]);
// 				b[i] = one - k[index+2]*k[index+dsn+3] + k[index+dsn+3]*k[index+1]*k[index+4]/(tw::small_pos + k[index+3]);
// 				c[i] = k[index+dsn+4];
// 			}
// 			TriDiagonal<T,T>(ans,source,a,b,c);
// 			for (i=0;i<=L;i++)
// 				k[ss*j + ds*i + 7] = ans[i]; // betas
// 		}
	}

	strip->Scatter(&k[0],mpi_packet,sizeof(T)*ds,0);
}

template <class T>
void GlobalIntegrator<T>::Parallelize()
{
	int i,j;

	// make temporary space for calculation of alpha and beta
	// packet layout : ss=system stride ; domain stride = ss*systems
	// layout must match assumptions in ComputeAlphasAndBetas
	tw::Int ss=8;
	std::valarray<T> mpi_packet(ss*systems);

	for (j=0;j<systems;j++)
	{
		mpi_packet[0 + ss*j] = theData[j][stride[j]];
		mpi_packet[1 + ss*j] = v[j*(N+2) + N];
		mpi_packet[2 + ss*j] = w[j*(N+2) + N];
		mpi_packet[3 + ss*j] = v[j*(N+2) + 1];
		mpi_packet[4 + ss*j] = w[j*(N+2) + 1];
		mpi_packet[5 + ss*j] = theData[j][N*stride[j]];
		mpi_packet[6 + ss*j] = 0.0; // alpha, to be computed
		mpi_packet[7 + ss*j] = 0.0; // beta, to be computed
	}

	ComputeAlphasAndBetas<T>(strip,systems,&mpi_packet[0]);

	// APPLY THE CORRECTION FACTORS

	for (j=0;j<systems;j++)
	{
		theData[j][0] = mpi_packet[6 + ss*j];
		for (i=1;i<=N;i++)
			theData[j][i*stride[j]] -= mpi_packet[6 + ss*j]*v[j*(N+2)+i] + mpi_packet[7 + ss*j]*w[j*(N+2)+i];
		theData[j][(N+1)*stride[j]] = mpi_packet[7 + ss*j];

		// N.b. exterior ghost cell calculation only valid for periodic or dirichlet boundary conditions
		// (interior domains OK for any boundary conditions)
		// Currently dirichlet b.c. is hard coded as 0
		// Other cases must be handled outside this routine
	}
}
