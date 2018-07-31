#include "definitions.h"
#include "tasks.h"
#include "ctools.h"
#include "3dmath.h"
#include "metricSpace.h"
#include "3dfields.h"
#include "numerics.h"

void Transform(tw::Float *array,tw::Int num,tw::Int interval,std::valarray<tw::Float>& transform)
{
	tw::Int j,n;
	std::valarray<tw::Float> temp(num);

	for (j=0;j<num;j++)
		temp[j] = 0.0;
	for (n=0;n<num;n++)
		for (j=0;j<num;j++)
			temp[j] += transform[j*num+n] * array[n*interval];
	for (j=0;j<num;j++)
		array[j*interval] = temp[j];
}

void LeftHalfRotation(ScalarField& A,tw::Int p,tw::Int q,tw::Float c,tw::Float s,tw::Int def1,tw::Int def2)
{
	tw::Int i;
	tw::Float temp;

	// First index is column (opposite usual subscript notation)

	for (i=1+def1;i<=A.Dim(1)-def2;i++)
	{
		temp = c*A(i,p,1) - s*A(i,q,1);
		A(i,q,1) = s*A(i,p,1) + c*A(i,q,1);
		A(i,p,1) = temp;
	}
}

void RightHalfRotation(ScalarField& A,tw::Int p,tw::Int q,tw::Float c,tw::Float s,tw::Int def1,tw::Int def2)
{
	tw::Int j;
	tw::Float temp;

	// First index is column (opposite usual subscript notation)

	for (j=1+def1;j<=A.Dim(2)-def2;j++)
	{
		temp = c*A(p,j,1) - s*A(q,j,1);
		A(q,j,1) = s*A(p,j,1) + c*A(q,j,1);
		A(p,j,1) = temp;
	}
}

void QRSweep(ScalarField& A,tw::Int deflation)
{
	tw::Int j,p,q,n;
	tw::Float temp;
	n = A.Dim(2);
	std::valarray<tw::Float> s(n+1),c(n+1);

	// Form the R matrix
	for (j=2+deflation;j<=n;j++)
	{
		p = j-1;
		q = j;
		temp = A(p,p,1)/A(p,q,1);
		s[j] = pow(one + temp*temp,tw::Float(-half));
		c[j] = sqrt(one - s[j]*s[j]);
		if (temp>0.0)
			s[j] *= -1.0;
		LeftHalfRotation(A,p,q,c[j],s[j],deflation,0);
	}

	// Form RQ
	for (j=2+deflation;j<=n;j++)
	{
		p = j-1;
		q = j;
		RightHalfRotation(A,p,q,c[j],s[j],deflation,0);
	}
}

void QLSweep(ScalarField& A,tw::Int deflation)
{
	tw::Int i,p,q,n;
	tw::Float temp;
	n = A.Dim(1);
	std::valarray<tw::Float> s(n+1),c(n+1);

	// Form the L matrix
	for (i=2;i<=n-deflation;i++)
	{
		p = i-1;
		q = i;
		temp = A(p,p,1)/A(q,p,1);
		s[i] = pow(one + temp*temp,tw::Float(-half));
		c[i] = sqrt(one - s[i]*s[i]);
		if (temp>0.0)
			s[i] *= -1.0;
		RightHalfRotation(A,p,q,c[i],s[i],0,deflation);
	}

	// Form LQ
	for (i=2;i<=n-deflation;i++)
	{
		p = i-1;
		q = i;
		LeftHalfRotation(A,p,q,c[i],s[i],0,deflation);
	}
}

void GetEigenvector(tw::Float eigenvalue,std::valarray<tw::Float>& vec,std::valarray<tw::Float>& a,std::valarray<tw::Float>& b,std::valarray<tw::Float>& c)
{
	tw::Int i,iter;
	tw::Float norm;
	bool anotherIteration;
	std::valarray<tw::Float> ans(vec.size());
	UniformDeviate ud(1); // must have same seed on all nodes

	const tw::Int maxIterations = 100;
	const tw::Int randomizePeriod = 5;

	b -= eigenvalue;

	iter = 0;
	do
	{
		if (iter%randomizePeriod==0) // generate new random guess
		{
			norm = 0.0;
			for (i=0;i<vec.size();i++)
			{
				vec[i] = ud.Next();
				norm += vec[i]*vec[i];
			}
			vec *= 1.0/sqrt(norm);
		}

		TriDiagonal<tw::Float,tw::Float>(ans,vec,a,b,c);
		norm = 0.0;
		for (i=0;i<vec.size();i++)
			norm += ans[i]*ans[i];
		ans *= 1.0/sqrt(norm);

		anotherIteration = false;
		for (i=0;i<vec.size();i++)
			if (fabs((ans[i]-vec[i])/vec[i])>1e-10)
				anotherIteration = true;

		vec = ans;
		iter++;
	} while (anotherIteration && iter<maxIterations);

	b += eigenvalue;
}

void NormalizeLeftRight(std::valarray<tw::Float>& left,std::valarray<tw::Float>& right)
{
	tw::Int i;
	tw::Float norm = 0.0;
	for (i=0;i<left.size();i++)
		norm += left[i]*right[i];
	left *= 1.0/norm;
}

void SymmetricTridiagonalEigensystem(std::valarray<tw::Float>& eigenvalues,tw::Float *revTransform,std::valarray<tw::Float>& T1,std::valarray<tw::Float>& T2)
{
	// this is basically "tqli" from "numerical recipes" with zero offset inputs
	// revTransform must contain identity matrix, or can be NULL if eigenvectors not needed
	tw::Int m,l,i,n,k;
	tw::Float s,r,p,g,f,c,b,dd;
	std::valarray<tw::Float> d,e;
	n = T1.size();

	d.resize(n+1);
	e.resize(n+1);
	for (i=1;i<=n;i++)
	{
		d[i] = T2[i-1];
		e[i] = T1[i-1];
	}
	for (i=2;i<=n;i++)
		e[i-1] = e[i];
	e[n] = 0.0;

	for (l=1;l<=n;l++)
	{
		do
		{
			for (m=l;m<=n-1;m++)
			{
				dd = fabs(d[m]) + fabs(d[m+1]);
				if (fabs(e[m]) + dd == dd)
					break;
			}
			if (m!=l)
			{
				g = (d[l+1]-d[l])/(2.0*e[l]);
				r = pythag(g,1.0);
				g = d[m]-d[l]+e[l]/(g+(g>=0.0 ? fabs(r) : -fabs(r)));
				s = c = 1.0;
				p = 0.0;
				for (i=m-1;i>=l;i--)
				{
					f = s*e[i];
					b = c*e[i];
					e[i+1] = (r=pythag(f,g));
					if (r==0.0)
					{
						d[i+1] -= p;
						e[m] = 0.0;
						break;
					}
					s = f/r;
					c = g/r;
					g = d[i+1] - p;
					r = (d[i]-g)*s + 2.0*c*b;
					d[i+1] = g + (p=s*r);
					g = c*r - b;
					if (revTransform!=NULL)
					{
						for (k=0;k<n;k++)
						{
							f = revTransform[k*n+i];
							revTransform[k*n+i] = s*revTransform[k*n+i-1] + c*f;
							revTransform[k*n+i-1] = c*revTransform[k*n+i-1] - s*f;
						}
					}
				}
				if (r==0.0 && i>=l) continue;
				d[l] -= p;
				e[l] = g;
				e[m] = 0.0;
			}
		} while (m!=l);
	}
	for (i=0;i<n;i++)
		eigenvalues[i] = d[i+1];
}

tw::Float GetSphericalGroundState(std::valarray<tw::Float>& vec,std::valarray<tw::Float>& phi,tw::Float dr)
{
	// Diagonalize Hamiltonion : H0 = -0.5*del^2 - 1/sqrt(coreRadius^2 + r^2)
	tw::Int i,dim;
	dim = phi.size();
	tw::Float dr1,dr3,A1,A3,V2;
	tw::Float r0,r1,r2,r3,r4;
	tw::Float eigenvalue,normalization;

	std::valarray<tw::Float> T1,T2,T3,temp,Lambda;

	vec.resize(dim);
	T1.resize(dim);
	T2.resize(dim);
	T3.resize(dim);
	temp.resize(dim);
	Lambda.resize(dim);

	// Form the tridiagonal matrix to diagonalize
	// 0,1,2,3,4 = center, wall, center, wall, center
	r0 = -0.5*dr;
	r1 = 0.0;
	r2 = 0.5*dr;
	for (i=0;i<dim;i++)
	{
		r3 = r2 + 0.5*dr;
		r4 = r3 + 0.5*dr;
		dr1 = r2 - r0;
		dr3 = r4 - r2;
		A1 = 4.0*pi*sqr(r1);
		A3 = 4.0*pi*sqr(r3);
		V2 = (4.0/3.0)*pi*(cub(r3) - cub(r1));
		T1[i] = -0.5*A1/(dr1*V2);
		T2[i] = 0.5*(A1/dr1 + A3/dr3)/V2 - phi[i];
		T3[i] = -0.5*A3/(dr3*V2);
		r0 = r2;
		r1 = r3;
		r2 = r4;
		Lambda[i] = sqrt(V2);
	}
	// set neumann condition at r = 0 (this is actually redundant since T1[0]=0)
	T2[0] += T1[0];
	// set dirichlet condition at r = maxRadius
	T2[dim-1] -= T3[dim-1];

	// Symmetrize the matrix ; i.e., form Diag[Lambda] * Tri[T1,T2,T3] * Diag[Lambda]^-1
	// Note that if matrix were complex, this would make it Hermitian
	for (i=0;i<dim;i++)
	{
		T1[i] *= Lambda[i];
		T2[i] *= Lambda[i];
		T3[i] *= Lambda[i];
	}
	T2[0] /= Lambda[0];
	T1[1] /= Lambda[0];
	for (i=1;i<dim-1;i++)
	{
		T1[i+1] /= Lambda[i];
		T2[i] /= Lambda[i];
		T3[i-1] /= Lambda[i];
	}
	T2[dim-1] /= Lambda[dim-1];
	T3[dim-2] /= Lambda[dim-1];

	// Get least eigenvalue of symmetrized matrix
	SymmetricTridiagonalEigensystem(temp,NULL,T1,T2);
	eigenvalue = temp[0];
	for (i=1;i<dim;i++)
		if (temp[i]<eigenvalue)
			eigenvalue = temp[i];

	// Get the right eigenvector, re-scale, and normalize
	GetEigenvector(eigenvalue,vec,T1,T2,T3);
	for (i=0;i<dim;i++)
		vec[i] /= Lambda[i];
	normalization = 0.0;
	for (i=0;i<dim;i++)
		normalization += vec[i]*vec[i]*(4.0/3.0)*pi*(cub(dr*tw::Float(i+1)) - cub(dr*tw::Float(i)));
	vec /= sqrt(normalization);
	return eigenvalue;
}

tw::Float GetCylindricalGroundState(std::valarray<tw::Float>& vec,std::valarray<tw::Float>& phi,tw::Float dr)
{
	// Diagonalize Hamiltonion : H0 = -0.5*del^2 - 1/sqrt(coreRadius^2 + r^2)
	tw::Int i,dim;
	dim = phi.size();
	tw::Float dr1,dr3,A1,A3,V2;
	tw::Float r0,r1,r2,r3,r4;
	tw::Float eigenvalue,normalization;

	std::valarray<tw::Float> T1,T2,T3,temp,Lambda;

	vec.resize(dim);
	T1.resize(dim);
	T2.resize(dim);
	T3.resize(dim);
	temp.resize(dim);
	Lambda.resize(dim);

	// Form the tridiagonal matrix to diagonalize
	// 0,1,2,3,4 = center, wall, center, wall, center
	r0 = -0.5*dr;
	r1 = 0.0;
	r2 = 0.5*dr;
	for (i=0;i<dim;i++)
	{
		r3 = r2 + 0.5*dr;
		r4 = r3 + 0.5*dr;
		dr1 = r2 - r0;
		dr3 = r4 - r2;
		A1 = 2.0*pi*r1;
		A3 = 2.0*pi*r3;
		V2 = pi*(sqr(r3) - sqr(r1));
		T1[i] = -0.5*A1/(dr1*V2);
		T2[i] = 0.5*(A1/dr1 + A3/dr3)/V2 - phi[i];
		T3[i] = -0.5*A3/(dr3*V2);
		r0 = r2;
		r1 = r3;
		r2 = r4;
		Lambda[i] = sqrt(V2);
	}
	// set neumann condition at r = 0 (this is actually redundant since T1[0]=0)
	T2[0] += T1[0];
	// set dirichlet condition at r = maxRadius
	T2[dim-1] -= T3[dim-1];

	// Symmetrize the matrix ; i.e., form Diag[Lambda] * Tri[T1,T2,T3] * Diag[Lambda]^-1
	// Note that if matrix were complex, this would make it Hermitian
	for (i=0;i<dim;i++)
	{
		T1[i] *= Lambda[i];
		T2[i] *= Lambda[i];
		T3[i] *= Lambda[i];
	}
	T2[0] /= Lambda[0];
	T1[1] /= Lambda[0];
	for (i=1;i<dim-1;i++)
	{
		T1[i+1] /= Lambda[i];
		T2[i] /= Lambda[i];
		T3[i-1] /= Lambda[i];
	}
	T2[dim-1] /= Lambda[dim-1];
	T3[dim-2] /= Lambda[dim-1];

	// Get least eigenvalue of symmetrized matrix
	SymmetricTridiagonalEigensystem(temp,NULL,T1,T2);
	eigenvalue = temp[0];
	for (i=1;i<dim;i++)
		if (temp[i]<eigenvalue)
			eigenvalue = temp[i];

	// Get the right eigenvector, re-scale, and normalize
	GetEigenvector(eigenvalue,vec,T1,T2,T3);
	for (i=0;i<dim;i++)
		vec[i] /= Lambda[i];
	normalization = 0.0;
	for (i=0;i<dim;i++)
		normalization += vec[i]*vec[i]*pi*(sqr(dr*tw::Float(i+1)) - sqr(dr*tw::Float(i)));
	vec /= sqrt(normalization);
	return eigenvalue;
}

void ComputeTransformMatrices(std::valarray<tw::Float>& eigenvalue,std::valarray<tw::Float>& fwd,std::valarray<tw::Float>& rev,MetricSpace *space,Task *task)
{
	tw::Int i,j,dim;
	tw::Float dr1,dr2;
	std::valarray<tw::Float> T1,T2,T3,vec,temp,Lambda;
	std::valarray<tw::Float> dr,A,V;

	dim = task->globalCells[1];
	eigenvalue.resize(dim);
	fwd.resize(dim*dim);
	rev.resize(dim*dim);
	T1.resize(dim);
	T2.resize(dim);
	T3.resize(dim);
	vec.resize(dim);
	temp.resize(dim);
	Lambda.resize(dim);
	dr.resize(dim+2);
	A.resize(dim+2);
	V.resize(dim+2);

	fwd = 0.0;
	rev = 0.0;
	for (i=0;i<dim;i++)
	{
		fwd[i*dim + i] = 1.0;
		rev[i*dim + i] = 1.0;
	}

	// Find parameters describing radial grid spacings
	for (i=0;i<=task->localCells[1]+1;i++)
	{
		dr[i] = space->dX(i,1);
		A[i] = space->dS(i,1,1,1);
		V[i] = space->dS(i,1,1,0);
	}
	task->strip[1].Gather(&dr[1],&dr[1],task->localCells[1]*sizeof(tw::Float),0);
	task->strip[1].Bcast(&dr[0],(dim+2)*sizeof(tw::Float),0);
	task->strip[1].Gather(&A[1],&A[1],task->localCells[1]*sizeof(tw::Float),0);
	task->strip[1].Bcast(&A[0],(dim+2)*sizeof(tw::Float),0);
	task->strip[1].Gather(&V[1],&V[1],task->localCells[1]*sizeof(tw::Float),0);
	task->strip[1].Bcast(&V[0],(dim+2)*sizeof(tw::Float),0);
	dr[dim+1] = dr[dim];
	A[dim+1] = A[dim];
	V[dim+1] = V[dim];

	// Form the tridiagonal matrix to diagonalize
	for (i=1;i<=dim;i++)
	{
		dr1 = 0.5*dr[i] + 0.5*dr[i-1];
		dr2 = 0.5*dr[i] + 0.5*dr[i+1];
		T1[i-1] = A[i]/(dr1*V[i]);
		T2[i-1] = - (A[i]/dr1 + A[i+1]/dr2)/V[i];
		T3[i-1] = A[i+1]/(dr2*V[i]);
		Lambda[i-1] = sqrt(V[i]);
	}
	// set neumann condition at r = 0 (this is actually redundant since T1[0]=0)
	T2[0] += T1[0];
	// set dirichlet condition at r = R
	T2[dim-1] -= T3[dim-1];

	// Symmetrize the matrix ; i.e., form Diag[Lambda] * Tri[T1,T2,T3] * Diag[Lambda]^-1
	// Note that if matrix were complex, this would make it Hermitian
	for (i=0;i<dim;i++)
	{
		T1[i] *= Lambda[i];
		T2[i] *= Lambda[i];
		T3[i] *= Lambda[i];
	}
	T2[0] /= Lambda[0];
	T1[1] /= Lambda[0];
	for (i=1;i<dim-1;i++)
	{
		T1[i+1] /= Lambda[i];
		T2[i] /= Lambda[i];
		T3[i-1] /= Lambda[i];
	}
	T2[dim-1] /= Lambda[dim-1];
	T3[dim-2] /= Lambda[dim-1];

	// Get eigenvalues and reverse transform (columns of right eigenvectors)
	SymmetricTridiagonalEigensystem(eigenvalue,&rev[0],T1,T2);

	// Form forward transform by transposing reverse transform
	for (j=0;j<dim;j++)
		for (i=0;i<dim;i++)
			fwd[i*dim+j] = rev[j*dim+i];

	// Transform back to density representation ; i.e., form Diag[Lambda]^-1 * xfrm * Diag[Lambda]
	for (j=0;j<dim;j++)
		for (i=0;i<dim;i++)
		{
			fwd[j*dim+i] /= Lambda[j];
			rev[j*dim+i] /= Lambda[j];
		}
	for (j=0;j<dim;j++)
		for (i=0;i<dim;i++)
		{
			fwd[j*dim+i] *= Lambda[i];
			rev[j*dim+i] *= Lambda[i];
		}
}


////////////////////
//                //
//   FCT ENGINE   //
//                //
////////////////////


FCT_Engine::FCT_Engine(const StripIterator& s,const MetricSpace& m,ScalarField *fluxMask)
{
	tw::Int i,ax=s.Axis();
	cells = s.Dim();
	V.resize(cells+2);
	A.resize(cells+2);
	scratch.resize(cells+2);
	// The engine uses V[0] in the clipping stage
	// The engine never uses A[0]
	for (i=0;i<=cells+1;i++)
		V[i] = m.dS(s,i,0);
	for (i=1;i<=cells+1;i++)
		A[i] = m.dS(s,i,ax);
	if (fluxMask!=NULL)
		for (i=1;i<=cells+1;i++)
			A[i] *= 1.0 - tw::Float( (*fluxMask)(s,i-1) + (*fluxMask)(s,i) == 1.0 );
}

void FCT_Engine::Reset(const StripIterator& s,const MetricSpace& m,ScalarField *fluxMask)
{
	tw::Int i,ax=s.Axis();
	cells = s.Dim();
	// The engine uses V[0] in the clipping stage
	// The engine never uses A[0]
	for (i=0;i<=cells+1;i++)
		V[i] = m.dS(s,i,0);
	for (i=1;i<=cells+1;i++)
		A[i] = m.dS(s,i,ax);
	if (fluxMask!=NULL)
		for (i=1;i<=cells+1;i++)
			A[i] *= 1.0 - tw::Float( (*fluxMask)(s,i-1) + (*fluxMask)(s,i) == 1.0 );
}

void FCT_Engine::Transport(std::valarray<tw::Float>& vel,
	std::valarray<tw::Float>& rho,
	std::valarray<tw::Float>& rho1,
	std::valarray<tw::Float>& diff,
	std::valarray<tw::Float>& flux,
	tw::Float dt)
{
	// First step in FCT algorithm
	// vel, rho, rho1, flux, and diff are indexed from 0..cells+1
	// vel contains velocity in the direction of transport
	// vel is known at t = 0 for a 1rst order push, at t = 1/2 for a second order push
	// diff and flux are given at the low-side cell walls, all other quantities at the cell centers
	// rho is density at t = 0
	// rho1 is density at t = 0 for a 1rst order push (rho1 = rho), or at t = 1/2 for a 2nd order push
	// diff is an output giving the diffused "mass" into the low side of each cell
	// flux is an output giving the transported and diffused "mass" into the low side of each cell
	// after calling this, the values of the low and high ghost cells for rho must be supplied

	tw::Int i;

	// compute diffused mass
	#pragma omp simd
	for (i=1;i<=cells+1;i++)
		diff[i] = -0.25*A[i]*dt*fabs(vel[i]+vel[i-1])*(rho[i] - rho[i-1]);

	// compute transported mass
	#pragma omp simd
	for (i=1;i<=cells+1;i++)
		flux[i] = A[i]*dt*0.25*(vel[i]+vel[i-1])*(rho1[i]+rho1[i-1]);

	// compute transported density
	#pragma omp simd
	for (i=1;i<=cells;i++)
		rho[i] += (flux[i] - flux[i+1])/V[i];

	// compute transported and diffused mass
	#pragma omp simd
	for (i=1;i<=cells+1;i++)
		flux[i] += diff[i];
}

void FCT_Engine::Diffuse(std::valarray<tw::Float>& vel,
	std::valarray<tw::Float>& rho,
	std::valarray<tw::Float>& diff,
	tw::Float dt)
{
	// vel is the same as for FCT_Engine::Transport
	// rho is the output from FCT_Engine::Transport, with appropriate ghost cell values
	// diff is the array output from FCT_Engine::Transport, with appropriate ghost cell values
	// on output, diff is replaced by the anti-diffused mass into the low-side of each cell
	// after calling this, the values of the low/high ghost cells for rho must be supplied
	// HOWEVER, the GLOBAL ghost cells for rho should not be updated (leave them with the undiffused values)

	tw::Int i;
	scratch = rho;

	// compute transported and diffused density
	#pragma omp simd
	for (i=1;i<=cells;i++)
		rho[i] += (diff[i] - diff[i+1])/V[i];

	// compute anti-diffused mass (re-using diffused mass array)
	#pragma omp simd
	for (i=1;i<=cells+1;i++)
		diff[i] = 0.25*A[i]*dt*fabs(vel[i]+vel[i-1])*(scratch[i] - scratch[i-1]);
}

void FCT_Engine::Limiter(tw::Float& adiff,const tw::Float& maxLow,const tw::Float& maxHigh)
{
	tw::Float pos_channel=adiff,neg_channel=adiff,test;

	tw::Float maxHigh_pos = tw::Float(maxHigh>0.0);
	tw::Float maxLow_pos = tw::Float(maxLow>0.0);

	test = tw::Float(pos_channel>maxHigh);
	pos_channel = maxHigh_pos*(test*maxHigh + (1.0-test)*pos_channel);
	test = tw::Float(pos_channel>maxLow);
	pos_channel = maxLow_pos*(test*maxLow + (1.0-test)*pos_channel);

	test = tw::Float(neg_channel<maxHigh);
	neg_channel = (1.0-maxHigh_pos)*(test*maxHigh + (1.0-test)*neg_channel);
	test = tw::Float(neg_channel<maxLow);
	neg_channel = (1.0-maxLow_pos)*(test*maxLow + (1.0-test)*neg_channel);

	adiff = pos_channel*tw::Float(adiff>0.0) + neg_channel*tw::Float(adiff<=0.0);
}

void FCT_Engine::Clip(std::valarray<tw::Float>& rho,std::valarray<tw::Float>& adiff,tw::Float rho00)
{
	// Limit the antidiffusive fluxes given in adiff to ensure stability, positivity, etc.
	// rho00 gives the value in the cell below the low-side ghost cell
	// rho and adiff are the outputs from FCT_Engine::Diffuse
	// after calling, high ghost cell for adiff must be supplied

	tw::Int i;
	tw::Float maxAntiDiffusionLow,maxAntiDiffusionHigh;

	maxAntiDiffusionHigh = (rho[2]-rho[1])*V[1];
	maxAntiDiffusionLow = (rho[0] - rho00)*V[0];
	Limiter(adiff[1],maxAntiDiffusionLow,maxAntiDiffusionHigh);
	#pragma omp simd
	for (i=2;i<=cells;i++)
	{
		maxAntiDiffusionHigh = (rho[i+1]-rho[i])*V[i];
		maxAntiDiffusionLow = (rho[i-1]-rho[i-2])*V[i-1];
		Limiter(adiff[i],maxAntiDiffusionLow,maxAntiDiffusionHigh);
	}
}

void FCT_Engine::AntiDiffuse(std::valarray<tw::Float>& rho,std::valarray<tw::Float>& adiff,std::valarray<tw::Float>& flux)
{
	// adiff are the clipped fluxes from FCT_Engine::Clip
	// after calling, low and high ghost cells for rho must be supplied
	// flux is an ouptut giving the true flux through the cell walls
	// after calling, the low ghost cell for flux must be supplied

	tw::Int i;
	const tw::Float safety_factor = 0.999999;
	#pragma omp simd
	for (i=1;i<=cells;i++)
		rho[i] += safety_factor*(adiff[i] - adiff[i+1])/V[i];
	#pragma omp simd
	for (i=1;i<=cells+1;i++)
		flux[i] += safety_factor*adiff[i];
}



////////////////////
//                //
//   FCT DRIVER   //
//                //
////////////////////


FCT_Driver::FCT_Driver(Field *rho,Field *rho1,Field *vel,ScalarField *fluxMask,MetricSpace *ms)
{
	en = All(*rho);
	vi = 0;
	this->rho = rho;
	this->rho1 = rho1;
	this->vel = vel;
	this->fluxMask = fluxMask;
	this->ms = ms;
	diff = NULL;
	net_flux = NULL;
}

FCT_Driver::~FCT_Driver()
{
	if (diff!=NULL)
		delete diff;
	if (net_flux!=NULL)
		delete net_flux;
}

void FCT_Driver::Convect(const axisSpec& axis,boundarySpec low,boundarySpec high,tw::Float dt)
{
	const tw::Int ax = naxis(axis);
	const tw::Int N = ms->Dim(ax);

	if (diff!=NULL)
		delete diff;
	if (net_flux!=NULL)
		delete net_flux;
	diff = new Field;
	net_flux = new Field;
	diff->Initialize(en.Components(),*rho,rho->task);
	net_flux->Initialize(en.Components(),*rho,rho->task);
	diff->SetBoundaryConditions(axis,low,high);
	net_flux->SetBoundaryConditions(axis,low,high);

	#pragma omp parallel
	{
		tw::Int c;
		std::valarray<tw::Float> va_rho(N+2),va_rho1(N+2),va_vel(N+2),va_diff(N+2),va_flux(N+2);
		tw::Float va_rho00;
		StripIterator s(*ms,ax,strongbool::yes);
		FCT_Engine engine(s,*ms,fluxMask);

		// TRANSPORT

		for (c=en.low;c<=en.high;c++)
		{
			for (s=0;s<s.end();++s)
			{
				engine.Reset(s,*ms,fluxMask);
				rho->GetStrip(va_rho,s,c);
				rho1->GetStrip(va_rho1,s,c);
				vel->GetStrip(va_vel,s,vi);
				engine.Transport(va_vel,va_rho,va_rho1,va_diff,va_flux,dt);
				rho->SetStrip(va_rho,s,c);
				diff->SetStrip(va_diff,s,c-en.low);
				net_flux->SetStrip(va_flux,s,c-en.low);
			}
		}
		#pragma omp barrier
		#pragma omp single
		{
			rho->DownwardCopy(axis,en,1);
			rho->UpwardCopy(axis,en,1);
			rho->ApplyBoundaryCondition(en);
		}

		// DIFFUSE

		for (c=en.low;c<=en.high;c++)
		{
			for (s=0;s<s.end();++s)
			{
				engine.Reset(s,*ms,fluxMask);
				rho->GetStrip(va_rho,s,c);
				diff->GetStrip(va_diff,s,c-en.low);
				vel->GetStrip(va_vel,s,vi);
				engine.Diffuse(va_vel,va_rho,va_diff,dt);
				rho->SetStrip(va_rho,s,c);
				diff->SetStrip(va_diff,s,c-en.low);
			}
		}
		#pragma omp barrier
		#pragma omp single
		{
			rho->DownwardCopy(axis,en,1);
			rho->UpwardCopy(axis,en,2); // need 2 low-side ghost cells for clipping operation
			// leave ghost cell with transported un-diffused value (no call to rho->ApplyBoundaryCondition)
		}

		// CLIP

		for (c=en.low;c<=en.high;c++)
		{
			for (s=0;s<s.end();++s)
			{
				engine.Reset(s,*ms,fluxMask);
				rho->GetStrip(va_rho,s,c);
				diff->GetStrip(va_diff,s,c-en.low);
				va_rho00 = (*rho)(s,-1,c);
				engine.Clip(va_rho,va_diff,va_rho00);
				diff->SetStrip(va_diff,s,c-en.low);
			}
		}
		#pragma omp barrier
		#pragma omp single
		{
			diff->DownwardCopy(axis,1);
			diff->ApplyBoundaryCondition();
		}

		// ANTI-DIFFUSE

		for (c=en.low;c<=en.high;c++)
		{
			for (s=0;s<s.end();++s)
			{
				engine.Reset(s,*ms,fluxMask);
				rho->GetStrip(va_rho,s,c);
				diff->GetStrip(va_diff,s,c-en.low);
				net_flux->GetStrip(va_flux,s,c-en.low);
				engine.AntiDiffuse(va_rho,va_diff,va_flux);
				rho->SetStrip(va_rho,s,c);
				net_flux->SetStrip(va_flux,s,c-en.low);
			}
		}
		#pragma omp barrier
		#pragma omp single
		{
			rho->DownwardCopy(axis,en,1);
			rho->UpwardCopy(axis,en,1);
			rho->ApplyBoundaryCondition(en);
			net_flux->UpwardCopy(axis,1);
			net_flux->ApplyBoundaryCondition();
		}
	}
}
