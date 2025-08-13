module;

#include "tw_includes.h"

export module numerics:serial;
import dyn_space;
import metric_space;

/// @brief invert A * phi = rho, rows of A are (b,c,0,...)(a,b,c,0,...)(0,a,b,c,0,...)...(0,...,a,b,c)(0,...,a,b)
/// @tparam T type of the elements of the source and solution vectors
/// @tparam U type of the coefficients
/// @param phi solution vector
/// @param rho source vector
/// @param a lower subdiagonal
/// @param b diagonal
/// @param c upper subdiagonal
export template <class T, class U>
void TriDiagonal(std::valarray<T>& phi, std::valarray<T>& rho, U a, U b, U c)
{
	tw::Int n = phi.size();
	std::valarray<U> gam(n);

	if (n == 1)
	{
		phi[0] = rho[0] / b;
		return;
	}

	tw::Int i;
	U bet;
	bet = b;
	phi[0] = rho[0] / bet;
	for (i = 1; i <= n - 1; i++)
	{
		gam[i] = c / bet;
		bet = b - a * gam[i];
		phi[i] = (rho[i] - a * phi[i - 1]) / bet;
	}

	for (i = n - 2; i >= 0; i--)
		phi[i] -= gam[i + 1] * phi[i + 1];
}

/// @brief invert A * phi = rho, rows of A are (b,c,0,...)(a,b,c,0,...)(0,a,b,c,0,...)...(0,...,a,b,c)(0,...,a,b)
/// @tparam T type of the elements of the source and solution vectors
/// @tparam U type of the coefficients
/// @param phi solution vector
/// @param rho source vector
/// @param a lower subdiagonal
/// @param b diagonal
/// @param c upper subdiagonal
export template <class T, class U>
void TriDiagonal(std::valarray<T>& phi, std::valarray<T>& rho, U a, std::valarray<U>& b, U c)
{
	tw::Int n = phi.size();
	std::valarray<U> gam(n);

	if (n == 1)
	{
		phi[0] = rho[0] / b[0];
		return;
	}

	tw::Int i;
	U bet;
	bet = b[0];
	phi[0] = rho[0] / bet;
	for (i = 1; i <= n - 1; i++)
	{
		gam[i] = c / bet;
		bet = b[i] - a * gam[i];
		phi[i] = (rho[i] - a * phi[i - 1]) / bet;
	}

	for (i = n - 2; i >= 0; i--)
		phi[i] -= gam[i + 1] * phi[i + 1];
}

/// @brief invert A * phi = rho, rows of A are (b,c,0,...)(a,b,c,0,...)(0,a,b,c,0,...)...(0,...,a,b,c)(0,...,a,b)
/// @tparam T type of the elements of the source and solution vectors
/// @tparam U type of the coefficients
/// @param phi solution vector
/// @param rho source vector
/// @param a lower subdiagonal
/// @param b diagonal
/// @param c upper subdiagonal
export template <class T, class U>
void TriDiagonal(std::valarray<T>& phi, std::valarray<T>& rho, std::valarray<U>& a, std::valarray<U>& b, std::valarray<U>& c)
{
	tw::Int n = phi.size();
	std::valarray<U> gam(n);

	if (n == 1)
	{
		phi[0] = rho[0] / b[0];
		return;
	}

	tw::Int i;
	U bet;
	bet = b[0];
	phi[0] = rho[0] / bet;
	for (i = 1; i <= n - 1; i++)
	{
		gam[i] = c[i - 1] / bet;
		bet = b[i] - a[i] * gam[i];
		phi[i] = (rho[i] - a[i] * phi[i - 1]) / bet;
	}

	for (i = n - 2; i >= 0; i--)
		phi[i] -= gam[i + 1] * phi[i + 1];
}

export void Transform(tw::Float *array,tw::Int pts,tw::Int modes,tw::Int interval,std::valarray<tw::Float>& transform)
{
	// Truncation of modes assumes eigenvalues sorted with increasing *absolute* value
	tw::Int p,m;
	std::valarray<tw::Float> temp(modes);

	for (m=0;m<modes;m++)
		temp[m] = 0.0;
	for (p=0;p<pts;p++)
		for (m=0;m<modes;m++)
			temp[m] += transform[m*pts+p] * array[p*interval];
	for (m=0;m<modes;m++)
		array[m*interval] = temp[m];
	for (m=modes;m<pts;m++)
		array[m*interval] = 0.0;
}

export void ReverseTransform(tw::Float *array,tw::Int pts,tw::Int modes,tw::Int interval,std::valarray<tw::Float>& rev_transform)
{
	// Truncation of modes assumes eigenvalues sorted with increasing *absolute* value
	// This is the same operation as Transform only if pts=modes
	tw::Int p,m;
	std::valarray<tw::Float> temp(pts);

	for (p=0;p<pts;p++)
		temp[p] = 0.0;
	for (m=0;m<modes;m++)
		for (p=0;p<pts;p++)
			temp[p] += rev_transform[p*pts+m] * array[m*interval];
	for (p=0;p<pts;p++)
		array[p*interval] = temp[p];
}

// void LeftHalfRotation(ScalarField& A,tw::Int p,tw::Int q,tw::Float c,tw::Float s,tw::Int def1,tw::Int def2)
// {
// 	tw::Int i;
// 	tw::Float temp;

// 	// First index is column (opposite usual subscript notation)

// 	for (i=1+def1;i<=A.Dim(1)-def2;i++)
// 	{
// 		temp = c*A(i,p,1) - s*A(i,q,1);
// 		A(i,q,1) = s*A(i,p,1) + c*A(i,q,1);
// 		A(i,p,1) = temp;
// 	}
// }

// void RightHalfRotation(ScalarField& A,tw::Int p,tw::Int q,tw::Float c,tw::Float s,tw::Int def1,tw::Int def2)
// {
// 	tw::Int j;
// 	tw::Float temp;

// 	// First index is column (opposite usual subscript notation)

// 	for (j=1+def1;j<=A.Dim(2)-def2;j++)
// 	{
// 		temp = c*A(p,j,1) - s*A(q,j,1);
// 		A(q,j,1) = s*A(p,j,1) + c*A(q,j,1);
// 		A(p,j,1) = temp;
// 	}
// }

// void QRSweep(ScalarField& A,tw::Int deflation)
// {
// 	tw::Int j,p,q,n;
// 	tw::Float temp;
// 	n = A.Dim(2);
// 	std::valarray<tw::Float> s(n+1),c(n+1);

// 	// Form the R matrix
// 	for (j=2+deflation;j<=n;j++)
// 	{
// 		p = j-1;
// 		q = j;
// 		temp = A(p,p,1)/A(p,q,1);
// 		s[j] = pow(one + temp*temp,tw::Float(-half));
// 		c[j] = std::sqrt(one - s[j]*s[j]);
// 		if (temp>0.0)
// 			s[j] *= -1.0;
// 		LeftHalfRotation(A,p,q,c[j],s[j],deflation,0);
// 	}

// 	// Form RQ
// 	for (j=2+deflation;j<=n;j++)
// 	{
// 		p = j-1;
// 		q = j;
// 		RightHalfRotation(A,p,q,c[j],s[j],deflation,0);
// 	}
// }

// void QLSweep(ScalarField& A,tw::Int deflation)
// {
// 	tw::Int i,p,q,n;
// 	tw::Float temp;
// 	n = A.Dim(1);
// 	std::valarray<tw::Float> s(n+1),c(n+1);

// 	// Form the L matrix
// 	for (i=2;i<=n-deflation;i++)
// 	{
// 		p = i-1;
// 		q = i;
// 		temp = A(p,p,1)/A(q,p,1);
// 		s[i] = pow(one + temp*temp,tw::Float(-half));
// 		c[i] = std::sqrt(one - s[i]*s[i]);
// 		if (temp>0.0)
// 			s[i] *= -1.0;
// 		RightHalfRotation(A,p,q,c[i],s[i],0,deflation);
// 	}

// 	// Form LQ
// 	for (i=2;i<=n-deflation;i++)
// 	{
// 		p = i-1;
// 		q = i;
// 		LeftHalfRotation(A,p,q,c[i],s[i],0,deflation);
// 	}
// }

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
			vec *= 1.0/std::sqrt(norm);
		}

		TriDiagonal<tw::Float,tw::Float>(ans,vec,a,b,c);
		norm = 0.0;
		for (i=0;i<vec.size();i++)
			norm += ans[i]*ans[i];
		ans *= 1.0/std::sqrt(norm);

		anotherIteration = false;
		for (i=0;i<vec.size();i++)
			if (std::fabs((ans[i]-vec[i])/vec[i])>1e-10)
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

void SortEigensystem(std::valarray<tw::Float>& eigenvalues,tw::Float *revTransform)
{
	// Sort according to absolute value using simple insertion sort
	// Also sort revTransform matrix if it is not NULL
	tw::Int n = eigenvalues.size();
	for (tw::Int i=1;i<n;i++)
	{
		tw::Int j = i;
		while (std::fabs(eigenvalues[j-1]) > std::fabs(eigenvalues[j]))
		{
			std::swap(eigenvalues[j],eigenvalues[j-1]);
			if (revTransform!=NULL)
				for (tw::Int k=0;k<n;k++)
					std::swap(revTransform[k*n+j],revTransform[k*n+j-1]);
			j--;
			if (j==0)
				break;
		}
	}
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
				dd = std::fabs(d[m]) + std::fabs(d[m+1]);
				if (std::fabs(e[m]) + dd == dd)
					break;
			}
			if (m!=l)
			{
				g = (d[l+1]-d[l])/(2.0*e[l]);
				r = pythag(g,1.0);
				g = d[m]-d[l]+e[l]/(g+(g>=0.0 ? std::fabs(r) : -std::fabs(r)));
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
	SortEigensystem(eigenvalues,revTransform);
}

export tw::Float GetSphericalGroundState(std::valarray<tw::Float>& vec,std::valarray<tw::Float>& phi,tw::Float dr)
{
	// Diagonalize Hamiltonion : H0 = -0.5*del^2 - 1/std::sqrt(coreRadius^2 + r^2)
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
		Lambda[i] = std::sqrt(V2);
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
	vec /= std::sqrt(normalization);
	return eigenvalue;
}

export tw::Float GetCylindricalGroundState(std::valarray<tw::Float>& vec,std::valarray<tw::Float>& phi,tw::Float dr)
{
	// Diagonalize Hamiltonion : H0 = -0.5*del^2 - 1/std::sqrt(coreRadius^2 + r^2)
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
		Lambda[i] = std::sqrt(V2);
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
	vec /= std::sqrt(normalization);
	return eigenvalue;
}

export void ComputeTransformMatrices(tw::bc::fld radial_bc,std::valarray<tw::Float>& eigenvalue,std::valarray<tw::Float>& fwd,std::valarray<tw::Float>& rev,MetricSpace *space,Task *task)
{
	tw::Int i,j,dim;
	tw::Float dr1,dr2;
	std::valarray<tw::Float> T1,T2,T3,vec,temp,Lambda;
	std::valarray<tw::Float> dr,A,V;

	dim = space->GlobalDim(1);
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
	for (i=0;i<=space->Dim(1)+1;i++)
	{
		dr[i] = space->dX(i,1);
		A[i] = space->dS(i,1,1,1);
		V[i] = space->dS(i,1,1,0);
	}
	task->strip[1].Gather(&dr[1],&dr[1],space->Dim(1)*sizeof(tw::Float),0);
	task->strip[1].Bcast(&dr[0],(dim+2)*sizeof(tw::Float),0);
	task->strip[1].Gather(&A[1],&A[1],space->Dim(1)*sizeof(tw::Float),0);
	task->strip[1].Bcast(&A[0],(dim+2)*sizeof(tw::Float),0);
	task->strip[1].Gather(&V[1],&V[1],space->Dim(1)*sizeof(tw::Float),0);
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
		Lambda[i-1] = std::sqrt(V[i]);
	}
	// set neumann condition at r = 0 (this is actually redundant since T1[0]=0)
	T2[0] += T1[0];
	// set boundary condition at r = R
	switch (radial_bc)
	{
		case tw::bc::fld::dirichletCell:
			break;
		case tw::bc::fld::dirichletWall:
			T2[dim-1] -= T3[dim-1];
			break;
		case tw::bc::fld::neumannWall:
			T2[dim-1] += T3[dim-1];
			break;
		default:
			throw tw::FatalError("Unhandled boundary condition during diagonalization.");
			break;
	}

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
