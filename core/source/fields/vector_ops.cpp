module;

#include "tw_includes.h"

export module fields:vector_ops;
import :base;

import base;

export template <tw::Int X, tw::Int Y, tw::Int Z>
void add_const_vec(const tw::Int& n, Field& vf, const tw::vec3& v0)
{
#pragma omp parallel
	{
		for (auto v : VectorStripRange<3>(vf,0,n,true))
		{
#pragma omp simd
			for (tw::Int k = 0; k <= vf.UNG(3); k++)
				vf(v, k, X) += v0.x;
#pragma omp simd
			for (tw::Int k = 0; k <= vf.UNG(3); k++)
				vf(v, k, Y) += v0.y;
#pragma omp simd
			for (tw::Int k = 0; k <= vf.UNG(3); k++)
				vf(v, k, Z) += v0.z;
		}
	}
}

export template<tw::Int T, tw::Int X, tw::Int Y, tw::Int Z>
void conserved_current_to_dens(const tw::Int& n, Field& current, const MetricSpace& m)
{
#pragma omp parallel
	{
		for (auto cell : EntireCellRange(m,n))
		{
			current(cell, T) /= m.dS(cell, 0);
			current(cell, X) /= m.dS(cell, 1) + tw::small_pos;
			current(cell, Y) /= m.dS(cell, 2) + tw::small_pos;
			current(cell, Z) /= m.dS(cell, 3) + tw::small_pos;
		}
	}
}

export template <tw::Int C, tw::Int X, tw::Int Y, tw::Int Z>
void assign_grad(const tw::Int& n, const Field& sf, Field& vf, const MetricSpace& m, const tw::Float& scaleFactor)
{
	// assign the gradient of a scalar field (sf) to a vector field (vf)
	const tw::Int xN1 = m.UNG(1), yN1 = m.UNG(2), zN1 = m.UNG(3);
#pragma omp parallel for collapse(3) schedule(static)
	for (tw::Int i = 1; i <= xN1; i++)
		for (tw::Int j = 1; j <= yN1; j++)
			for (tw::Int k = 1; k <= zN1; k++)
			{
				vf(n, i, j, k, X) = scaleFactor * (sf(n, i, j, k, C) - sf(n, i - 1, j, k, C)) / m.dl(i, j, k, 1);
				vf(n, i, j, k, Y) = scaleFactor * (sf(n, i, j, k, C) - sf(n, i, j - 1, k, C)) / m.dl(i, j, k, 2);
				vf(n, i, j, k, Z) = scaleFactor * (sf(n, i, j, k, C) - sf(n, i, j, k - 1, C)) / m.dl(i, j, k, 3);
			}
}

export template <tw::Int C, tw::Int X, tw::Int Y, tw::Int Z>
void add_grad(const tw::Int& n, const Field& sf, Field& vf, const MetricSpace& m, const tw::Float& scaleFactor)
{
	// add the gradient of a scalar field (sf) to a vector field (vf)
	const tw::Int xN1 = m.UNG(1), yN1 = m.UNG(2), zN1 = m.UNG(3);
#pragma omp parallel for collapse(3) schedule(static)
	for (tw::Int i = 1; i <= xN1; i++)
		for (tw::Int j = 1; j <= yN1; j++)
			for (tw::Int k = 1; k <= zN1; k++)
			{
				vf(n, i, j, k, X) += scaleFactor * (sf(n, i, j, k, C) - sf(n, i - 1, j, k, C)) / m.dl(i, j, k, 1);
				vf(n, i, j, k, Y) += scaleFactor * (sf(n, i, j, k, C) - sf(n, i, j - 1, k, C)) / m.dl(i, j, k, 2);
				vf(n, i, j, k, Z) += scaleFactor * (sf(n, i, j, k, C) - sf(n, i, j, k - 1, C)) / m.dl(i, j, k, 3);
			}
}

export template <tw::Int U, tw::Int V, tw::Int W, tw::Int X, tw::Int Y, tw::Int Z>
void add_curlB(const tw::Int& n, const Field& src, Field& dst, const MetricSpace& m, const tw::Float& scaleFactor)
{
	// curl of B-field on generalized Yee mesh
	const tw::Int xDim = m.Dim(1), yDim = m.Dim(2), zDim = m.Dim(3);
	const tw::Int xN1 = m.UNG(1), yN1 = m.UNG(2), zN1 = m.UNG(3);
#pragma omp parallel for collapse(2) schedule(static)
	for (tw::Int i = 1; i <= xN1; i++)
		for (tw::Int j = 1; j <= yDim; j++)
		{
			tw::xstrip<3> v(m, 0, {n,i,j,0});
			tw::xstrip<3> vj(m, 0, {n,i,j+1,0});
#pragma omp simd
			for (tw::Int k = 1; k <= zDim; k++)
			{
				dst(v, k, X) += scaleFactor * (src(v, k, V) * m.dlh(v, k, 2) - src(v, k + 1, V) * m.dlh(v, k + 1, 2)) / m.dS(v, k, 1);
				dst(v, k, X) += scaleFactor * (src(vj, k, W) * m.dlh(vj, k, 3) - src(v, k, W) * m.dlh(v, k, 3)) / m.dS(v, k, 1);
			}
		}
#pragma omp parallel for collapse(2) schedule(static)
	for (tw::Int i = 1; i <= xDim; i++)
		for (tw::Int j = 1; j <= yN1; j++)
		{
			tw::xstrip<3> v(m, 0, {n,i,j,0});
			tw::xstrip<3> vi(m, 0, {n,i+1,j,0});
#pragma omp simd
			for (tw::Int k = 1; k <= zDim; k++)
			{
				dst(v, k, Y) -= scaleFactor * (src(v, k, U) * m.dlh(v, k, 1) - src(v, k + 1, U) * m.dlh(v, k + 1, 1)) / m.dS(v, k, 2);
				dst(v, k, Y) -= scaleFactor * (src(vi, k, W) * m.dlh(vi, k, 3) - src(v, k, W) * m.dlh(v, k, 3)) / m.dS(v, k, 2);
			}
		}
#pragma omp parallel for collapse(2) schedule(static)
	for (tw::Int i = 1; i <= xDim; i++)
		for (tw::Int j = 1; j <= yDim; j++)
		{
			tw::xstrip<3> v(m, 0, {n,i,j,0});
			tw::xstrip<3> vi(m, 0, {n,i+1,j,0});
			tw::xstrip<3> vj(m, 0, {n,i,j+1,0});
#pragma omp simd
			for (tw::Int k = 1; k <= zN1; k++)
			{
				dst(v, k, Z) += scaleFactor * (src(v, k, U) * m.dlh(v, k, 1) - src(vj, k, U) * m.dlh(vj, k, 1)) / m.dS(v, k, 3);
				dst(v, k, Z) += scaleFactor * (src(vi, k, V) * m.dlh(vi, k, 2) - src(v, k, V) * m.dlh(v, k, 2)) / m.dS(v, k, 3);
			}
		}
}

export template <tw::Int U, tw::Int V, tw::Int W, tw::Int X, tw::Int Y, tw::Int Z>
void add_curlE(const tw::Int& n, const Field& src, Field& dst, const MetricSpace& m, const tw::Float& scaleFactor)
{
	// curl of E-field on generalized Yee mesh
	const tw::Int xN1 = m.UNG(1), yN1 = m.UNG(2), zN1 = m.UNG(3);
#pragma omp parallel for collapse(2) schedule(static)
	for (tw::Int i = 1; i <= xN1; i++)
		for (tw::Int j = 1; j <= yN1; j++)
		{
			tw::xstrip<3> v(m, 0, {n,i,j,0});
			tw::xstrip<3> vi(m, 0, {n,i-1,j,0});
			tw::xstrip<3> vj(m, 0, {n,i,j-1,0});
#pragma omp simd
			for (tw::Int k = 1; k <= zN1; k++)
			{
				dst(v, k, X) += scaleFactor * (src(v, k - 1, V) * m.dl(v, k - 1, 2) - src(v, k, V) * m.dl(v, k, 2)) / m.dSh(v, k, 1);
				dst(v, k, X) += scaleFactor * (src(v, k, W) * m.dl(v, k, 3) - src(vj, k, W) * m.dl(vj, k, 3)) / m.dSh(v, k, 1);
				dst(v, k, Y) -= scaleFactor * (src(v, k - 1, U) * m.dl(v, k - 1, 1) - src(v, k, U) * m.dl(v, k, 1)) / m.dSh(v, k, 2);
				dst(v, k, Y) -= scaleFactor * (src(v, k, W) * m.dl(v, k, 3) - src(vi, k, W) * m.dl(vi, k, 3)) / m.dSh(v, k, 2);
				dst(v, k, Z) += scaleFactor * (src(vj, k, U) * m.dl(vj, k, 1) - src(v, k, U) * m.dl(v, k, 1)) / m.dSh(v, k, 3);
				dst(v, k, Z) += scaleFactor * (src(v, k, V) * m.dl(v, k, 2) - src(vi, k, V) * m.dl(vi, k, 2)) / m.dSh(v, k, 3);
			}
		}
}

export template <tw::Int X, tw::Int Y, tw::Int Z>
tw::Float divE(const Field& vf, tw::Int n, tw::Int i, tw::Int j, tw::Int k, const MetricSpace& m)
{
	// divergence of E-field on generalized Yee mesh
	tw::Float ans, vol;
	vol = m.dS(i, j, k, 0);
	ans = vf(n, i + 1, j, k, X) * m.dS(i + 1, j, k, 1);
	ans -= vf(n, i, j, k, X) * m.dS(i, j, k, 1);
	ans += vf(n, i, j + 1, k, Y) * m.dS(i, j + 1, k, 2);
	ans -= vf(n, i, j, k, Y) * m.dS(i, j, k, 2);
	ans += vf(n, i, j, k + 1, Z) * m.dS(i, j, k + 1, 3);
	ans -= vf(n, i, j, k, Z) * m.dS(i, j, k, 3);
	return ans / vol;
}

template <tw::Int X, tw::Int Y, tw::Int Z>
inline tw::Float div(const Field& vf, tw::Int n, tw::Int i, tw::Int j, tw::Int k, const MetricSpace& m)
{
	// divergence of a vector field
	tw::Float ans, vol;
	vol = m.dS(i, j, k, 0);
	ans = (vf(n, i + 1, j, k, X) + vf(n, i, j, k, X)) * m.dS(i + 1, j, k, 1);
	ans -= (vf(n, i - 1, j, k, X) + vf(n, i, j, k, X)) * m.dS(i, j, k, 1);
	ans += (vf(n, i, j + 1, k, Y) + vf(n, i, j, k, Y)) * m.dS(i, j + 1, k, 2);
	ans -= (vf(n, i, j - 1, k, Y) + vf(n, i, j, k, Y)) * m.dS(i, j, k, 2);
	ans += (vf(n, i, j, k + 1, Z) + vf(n, i, j, k, Z)) * m.dS(i, j, k + 1, 3);
	ans -= (vf(n, i, j, k - 1, Z) + vf(n, i, j, k, Z)) * m.dS(i, j, k, 3);
	return 0.5 * ans / vol;
}

template <tw::Int C, tw::Int X, tw::Int Y, tw::Int Z>
inline tw::Float div(const Field& coeff, const Field& vf, tw::Int n, tw::Int i, tw::Int j, tw::Int k, const MetricSpace& m)
{
	// divergence of a scalar field times a vector field
	tw::Float ans, vol;
	vol = m.dS(i, j, k, 0);
	ans = (coeff(n, i + 1, j, k, C) + coeff(n, i, j, k, C)) * (vf(n, i + 1, j, k, X) + vf(n, i, j, k, X)) * m.dS(i + 1, j, k, 1);
	ans -= (coeff(n, i - 1, j, k, C) + coeff(n, i, j, k, C)) * (vf(n, i - 1, j, k, X) + vf(n, i, j, k, X)) * m.dS(i, j, k, 1);
	ans += (coeff(n, i, j + 1, k, C) + coeff(n, i, j, k, C)) * (vf(n, i, j + 1, k, Y) + vf(n, i, j, k, Y)) * m.dS(i, j + 1, k, 2);
	ans -= (coeff(n, i, j - 1, k, C) + coeff(n, i, j, k, C)) * (vf(n, i, j - 1, k, Y) + vf(n, i, j, k, Y)) * m.dS(i, j, k, 2);
	ans += (coeff(n, i, j, k + 1, C) + coeff(n, i, j, k, C)) * (vf(n, i, j, k + 1, Z) + vf(n, i, j, k, Z)) * m.dS(i, j, k + 1, 3);
	ans -= (coeff(n, i, j, k - 1, C) + coeff(n, i, j, k, C)) * (vf(n, i, j, k - 1, Z) + vf(n, i, j, k, Z)) * m.dS(i, j, k, 3);
	return 0.25 * ans / vol;
}

template <tw::Int C, tw::Int S>
inline tw::Float div(const Field& coeff, const Field& sf, tw::Int n, tw::Int i, tw::Int j, tw::Int k, const MetricSpace& m)
{
	// divergence of a scalar field times the gradient of another scalar field
	tw::Float ans, vol;
	vol = m.dS(i, j, k, 0);
	ans = (coeff(n, i + 1, j, k, C) + coeff(n, i, j, k, C)) * (sf(n, i + 1, j, k, S) - sf(n, i, j, k, S)) * m.dS(i + 1, j, k, 1) / m.dl(i + 1, j, k, 1);
	ans -= (coeff(n, i - 1, j, k, C) + coeff(n, i, j, k, C)) * (sf(n, i, j, k, S) - sf(n, i - 1, j, k, S)) * m.dS(i, j, k, 1) / m.dl(i, j, k, 1);
	ans += (coeff(n, i, j + 1, k, C) + coeff(n, i, j, k, C)) * (sf(n, i, j + 1, k, S) - sf(n, i, j, k, S)) * m.dS(i, j + 1, k, 2) / m.dl(i, j + 1, k, 2);
	ans -= (coeff(n, i, j - 1, k, C) + coeff(n, i, j, k, C)) * (sf(n, i, j, k, S) - sf(n, i, j - 1, k, S)) * m.dS(i, j, k, 2) / m.dl(i, j, k, 2);
	ans += (coeff(n, i, j, k + 1, C) + coeff(n, i, j, k, C)) * (sf(n, i, j, k + 1, S) - sf(n, i, j, k, S)) * m.dS(i, j, k + 1, 3) / m.dl(i, j, k + 1, 3);
	ans -= (coeff(n, i, j, k - 1, C) + coeff(n, i, j, k, C)) * (sf(n, i, j, k, S) - sf(n, i, j, k - 1, S)) * m.dS(i, j, k, 3) / m.dl(i, j, k, 3);
	return 0.5 * ans / vol;
}
