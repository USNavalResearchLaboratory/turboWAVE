module;

#include "tw_includes.h"
#include <memory>
#include "tw_logger.h"
export module numerics:parallel;
import :serial;
import dyn_space;
import metric_space;
import logger;

/// @brief This solves a tridiagonal system accounting for domain decomposition both along and across the solution vector
/// @tparam T the type of the elements of the solution vector
export template <class T>
class GlobalIntegrator
{
	tw::comm* strip;
	tw::Int N, systems, components;
	std::valarray<T> v, w;
	std::valarray<tw::Float*> theData;
	std::valarray<tw::Int> posStride,compStride;

public:
	/// @brief create an object that inverts tridiagonals no matter the decomposition
	/// @param strip the communicator for the strip, usually task->strip[axis], where axis = 1, 2, or 3
	/// @param systems how many strips are being solved
	/// @param cells how many cells per strip, not including ghost cells
	GlobalIntegrator(tw::comm* strip, tw::Int systems, tw::Int cells) {
		this->components = sizeof(T)/sizeof(tw::Float);
		this->strip = strip;
		this->systems = systems;
		N = cells;

		v.resize(systems * (N + 2));
		w.resize(systems * (N + 2));

		theData.resize(systems);
		posStride.resize(systems);
		compStride.resize(systems);
	}
	/// @brief set the data for the given system
	/// @param system index of the system being solved starting at zero, a system is typically one of the strips being solved
	/// @param theData pointer to the first element of data (the near ghost cell) in this system
	/// @param posStride stride in units of tw::Float from start of `T` to next start of `T`
	/// @param compStride stride in units of tw::Float from first `T` component to next `T` component
	void SetData(tw::Int system, tw::Float* row, tw::Int posStride, tw::Int compStride) {
		this->theData[system] = row;
		this->posStride[system] = posStride;
		this->compStride[system] = compStride;
	}
	void Pack(const T& val,const tw::Int& system, const tw::Int& idx) {
		const tw::Float *ptr = (tw::Float*)&val;
		for (auto i=0; i< components; i++) {
			theData[system][i*compStride[system]+idx*posStride[system]] = ptr[i];
		}
	}
	T Unpack(const tw::Int& system, const tw::Int& idx) {
		T ans;
		tw::Float *ptr = (tw::Float*)&ans;
		for (auto i=0; i<components; i++) {
			ptr[i] = theData[system][i*compStride[system]+idx*posStride[system]];
		}
		return ans;
	}
	void SetMatrix(tw::Int system, T a, T b, T c, T theta, T eta) {
		std::valarray<T> u(N), basisVector(N);

		// Set up the inversion vectors for one system

		basisVector = T(0.0);
		basisVector[0] = theta;
		TriDiagonal<T, T>(u, basisVector, a, b, c);
		v[std::slice(system * (N + 2) + 1, N, 1)] = u;

		basisVector[0] = 0.0;
		basisVector[N - 1] = eta;
		TriDiagonal<T, T>(u, basisVector, a, b, c);
		w[std::slice(system * (N + 2) + 1, N, 1)] = u;
	}
	void SetMatrix(tw::Int system, T a, std::valarray<T>& b, T c, T theta, T eta) {
		std::valarray<T> u(N), basisVector(N);

		// Set up the inversion vectors for one system

		basisVector = T(0.0);
		basisVector[0] = theta;
		TriDiagonal<T, T>(u, basisVector, a, b, c);
		v[std::slice(system * (N + 2) + 1, N, 1)] = u;

		basisVector[0] = 0.0;
		basisVector[N - 1] = eta;
		TriDiagonal<T, T>(u, basisVector, a, b, c);
		w[std::slice(system * (N + 2) + 1, N, 1)] = u;
	}
	void SetMatrix(tw::Int system, std::valarray<T>& a, std::valarray<T>& b, std::valarray<T>& c) {
		std::valarray<T> u(N), basisVector(N);

		// Set up the inversion vectors for one system

		basisVector = T(0.0);
		basisVector[0] = a[0];
		TriDiagonal<T, T>(u, basisVector, a, b, c);
		v[std::slice(system * (N + 2) + 1, N, 1)] = u;

		basisVector[0] = 0.0;
		basisVector[N - 1] = c[N - 1];
		TriDiagonal<T, T>(u, basisVector, a, b, c);
		w[std::slice(system * (N + 2) + 1, N, 1)] = u;
	}
	void ComputeAlphasAndBetas(tw::comm* strip, tw::Int systems, T* mpi_packet) {
		tw::Int i, j, index, ds = 8 * systems, ss = 8;
		tw::Int domains = strip->Get_size();
		tw::Int L = strip->Get_size() - 1;
		tw::Int currDomain = strip->Get_rank();
		tw::Int unknowns = L;
		std::valarray<T> source(unknowns), ans(unknowns), a(unknowns), b(unknowns), c(unknowns);
		std::valarray<T> k(ds * domains);

		strip->Gather(mpi_packet, &k[0], sizeof(T) * ds, 0);

		for (j = 0; j < systems; j++)
		{
			// APERIODIC CASE

			if (currDomain == 0 && domains > 1)
			{
				// For alphas domain 0 is not used
				// indices to tridiagonal arguments are therefore offset by -1
				for (i = 1; i < L; i++)
				{
					index = ss * j + ds * i;
					source[i - 1] = k[index - ds + 5] - k[index - ds + 2] * k[index] + k[index - ds + 2] * k[index + 4] * k[index + 5] / (tw::small_pos + k[index + 2]);
					a[i - 1] = k[index - ds + 1];
					b[i - 1] = one - k[index - ds + 2] * k[index + 3] + k[index - ds + 2] * k[index + 1] * k[index + 4] / (tw::small_pos + k[index + 2]);
					c[i - 1] = k[index - ds + 2] * k[index + 4] / (tw::small_pos + k[index + 2]);
				}
				source[L - 1] = k[ss * j + ds * (L - 1) + 5] - k[ss * j + ds * (L - 1) + 2] * k[ss * j + ds * L];
				a[L - 1] = k[ss * j + ds * (L - 1) + 1];
				b[L - 1] = one - k[ss * j + ds * (L - 1) + 2] * k[ss * j + ds * L + 3];
				TriDiagonal<T, T>(ans, source, a, b, c);
				for (i = 1; i <= L; i++)
					k[ss * j + ds * i + 6] = ans[i - 1]; // alphas

				// For betas domain L is not used
				source[0] = k[ss * j + ds] - k[ss * j + ds + 3] * k[ss * j + 5];
				b[0] = one - k[ss * j + 2] * k[ss * j + ds + 3];
				c[0] = k[ss * j + ds + 4];
				for (i = 1; i < L; i++)
				{
					index = ss * j + ds * i;
					source[i] = k[index + ds] - k[index + ds + 3] * k[index + 5] + k[index + ds + 3] * k[index + 1] * k[index] / (tw::small_pos + k[index + 3]);
					a[i] = k[index + ds + 3] * k[index + 1] / (tw::small_pos + k[index + 3]);
					b[i] = one - k[index + 2] * k[index + ds + 3] + k[index + ds + 3] * k[index + 1] * k[index + 4] / (tw::small_pos + k[index + 3]);
					c[i] = k[index + ds + 4];
				}
				TriDiagonal<T, T>(ans, source, a, b, c);
				for (i = 0; i < L; i++)
					k[ss * j + ds * i + 7] = ans[i]; // betas
			}
		}

		strip->Scatter(&k[0], mpi_packet, sizeof(T) * ds, 0);
	}
	void Parallelize() {
		// make temporary space for calculation of alpha and beta
		// packet layout : ss=system stride ; domain stride = ss*systems
		// layout must match assumptions in ComputeAlphasAndBetas
		tw::Int ss = 8;
		std::valarray<T> mpi_packet(ss * systems);

		for (auto j = 0; j < systems; j++)
		{
			T val;

			mpi_packet[0 + ss * j] = Unpack(j,1);
			mpi_packet[1 + ss * j] = v[j * (N + 2) + N];
			mpi_packet[2 + ss * j] = w[j * (N + 2) + N];
			mpi_packet[3 + ss * j] = v[j * (N + 2) + 1];
			mpi_packet[4 + ss * j] = w[j * (N + 2) + 1];
			mpi_packet[5 + ss * j] = Unpack(j,N);
			mpi_packet[6 + ss * j] = 0.0; // alpha, to be computed
			mpi_packet[7 + ss * j] = 0.0; // beta, to be computed
		}

		ComputeAlphasAndBetas(strip, systems, &mpi_packet[0]);

		// APPLY THE CORRECTION FACTORS

		for (auto j = 0; j < systems; j++)
		{
			Pack(mpi_packet[6 + ss * j], j, 0);
			for (auto i = 1; i <= N; i++) {
				T val = Unpack(j,i) - mpi_packet[6 + ss * j] * v[j * (N + 2) + i] - mpi_packet[7 + ss * j] * w[j * (N + 2) + i];
				Pack(val,j,i);
			}
			Pack( mpi_packet[7 + ss * j], j, N+1);

			// N.b. exterior ghost cell calculation only valid for periodic or dirichlet boundary conditions
			// (interior domains OK for any boundary conditions)
			// Currently dirichlet b.c. is hard coded as 0
			// Other cases must be handled outside this routine
		}
	}
};

/// @brief Parallel cubic spline, serial part is from Stoer and Burlisch
/// @tparam T type we are splining
export template <class T>
class GlobalSpline {
	tw::Int dim,components;
	tw::comm *strip_comm;
	std::valarray<tw::Int> posStride,compStride;
	std::valarray<tw::Float*> theData;
	std::valarray<T> moments;
	std::unique_ptr<GlobalIntegrator<T>> gintegrator;
public:
	GlobalSpline(tw::comm* strip_comm, tw::Int strips, tw::Int cells) {
		this->components = sizeof(T)/sizeof(tw::Float);
		this->strip_comm = strip_comm;
		dim = cells;
		gintegrator = std::make_unique<GlobalIntegrator<T>>(strip_comm,strips,cells);
		moments.resize(strips*(cells+2));
		posStride.resize(strips);
		compStride.resize(strips);
		theData.resize(strips);
	}
	void Pack(const T& val,const tw::Int& strip, const tw::Int& idx) {
		const tw::Float *ptr = (tw::Float*)&val;
		for (auto i=0; i< components; i++) {
			theData[strip][i*compStride[strip]+idx*posStride[strip]] = ptr[i];
		}
	}
	T Unpack(const tw::Int& strip, const tw::Int& idx) {
		T ans;
		tw::Float *ptr = (tw::Float*)&ans;
		for (auto i=0; i<components; i++) {
			ptr[i] = theData[strip][i*compStride[strip]+idx*posStride[strip]];
		}
		return ans;
	}
	/// @brief thread safe function to setup one strip
	/// @param ax axis of the strip
	/// @param strip index of the strip
	/// @param theData pointer to near ghost cell of the strip's data
	/// @param posStride stride of the strip's positions
	/// @param compStride stride of the strip's components
	void SetStrip(tw::Int strip, tw::Float* theData, tw::Int posStride, tw::Int compStride) {
		std::valarray<T> phi(dim);
		std::valarray<T> rho(dim);
		std::valarray<T> a(dim);
		std::valarray<T> b(dim);
		std::valarray<T> c(dim);

		this->theData[strip] = theData;
		this->posStride[strip] = posStride;
		this->compStride[strip] = compStride;
		
		for (auto i = 1; i <= dim; i++) {
			rho[i-1] = tw::Float(3)*(Unpack(strip,i-1) - two*Unpack(strip,i) + Unpack(strip,i+1));
			a[i-1] = 0.5;
			b[i-1] = 2.0;
			c[i-1] = 0.5;
		}

		tw::Int n0,n1;
		strip_comm->Shift(1,1,&n0,&n1);
		if (n0==MPI_PROC_NULL) {
			c[0] = 0.0;
			rho[0] = 0.0;
		}
		if (n1==MPI_PROC_NULL) {
			a[dim-1] = 0.0;
			rho[dim-1] = 0.0;
		}

		TriDiagonal<T,T>(phi,rho,a,b,c);
		for (auto i = 1; i <= dim; i++) {
			moments[(dim+2)*strip + i] = phi[i-1];
		}
		gintegrator->SetData(strip,(tw::Float*)&moments[strip*(dim+2)],components,1);
		gintegrator->SetMatrix(strip,a,b,c);
	}
	void Solve() {
		gintegrator->Parallelize();
	}
	/// @brief get value of the spline at a point
	/// @param strip index of the strip
	/// @param x fractional point, interpoint spacing is always normalized
	/// @param theData pointer to near ghost cell of the strip's data
	/// @param posStride stride of the strip's positions
	/// @param compStride stride of the strip's components
	/// @return interpolated value
	T Interpolate(tw::Float x,tw::Int strip) {
		tw::Int i = x;
		tw::Float w = x - i;
		T m1 = moments[(dim+2)*strip + i];
		T m2 = moments[(dim+2)*strip + i + 1];
		T y1 = Unpack(strip,i);
		T y2 = Unpack(strip,i+1);
		return y1 + w*(y2-y1-(two*m1+m2)/tw::Float(6)) + half*w*w*m1 + w*w*w*(m2-m1)/tw::Float(6);
	}
};
