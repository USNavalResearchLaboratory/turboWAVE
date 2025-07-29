module;

#include "tw_includes.h"

/// This module provides easy access to aggregated interal dimensions.
/// As a simple example, a Vec3Field
/// is a Field with internal dimension 3, that works directly with `tw::vec3` types.
/// Accessors are also provided that assume dim[0] = 1.
export module fields:aggregates;
import :base;

import static_space;
import dyn_space;
import metric_space;
import tw_iterator;
import numerics;
import fft;
import logger;
#include "tw_logger.h"


export struct ScalarField: Field
{
	private:
	// hide the internal dimension initializer
	void Initialize(const tw::Int& components,const StaticSpace& ss,Task *task) {}

	public:
	// force the internal dimension
	void Initialize(const StaticSpace& ss,Task *task) {
		Field::Initialize(1,ss,task);
	}

	// gather and scatter

	void Interpolate(tw::Float *val, const weights_3D& weights) const
	{
		std::valarray<tw::Float> temp(*val,1);
		Field::Interpolate(Rng(0),temp,weights);
		*val = temp[0];
	}

	void InterpolateOnto(const tw::Float& val, const weights_3D& weights)
	{
		std::valarray<tw::Float> temp(val,1);
		Field::InterpolateOnto(Rng(0),temp,weights);
	}

	// integral transform

	void AxialSineTransform(const DynSpace& ds) {
		Field::RealAxialSineTransform(Rng(0),ds);
	}
	void InverseAxialSineTransform(const DynSpace& ds) {
		Field::RealInverseAxialSineTransform(Rng(0),ds);
	}
	void TransverseSineTransform(const DynSpace& ds) {
		Field::RealTransverseSineTransform(Rng(0),ds);
	}
	void InverseTransverseSineTransform(const DynSpace& ds) {
		Field::RealInverseTransverseSineTransform(Rng(0),ds);
	}
	void TransverseCosineTransform(const DynSpace& ds) {
		Field::RealTransverseCosineTransform(Rng(0),ds);
	}
	void InverseTransverseCosineTransform(const DynSpace& ds) {
		Field::RealInverseTransverseCosineTransform(Rng(0),ds);
	}
	void TransverseFFT(const DynSpace& ds) {
		Field::RealTransverseFFT(Rng(0),ds);
	}
	void InverseTransverseFFT(const DynSpace& ds) {
		Field::RealInverseTransverseFFT(Rng(0),ds);
	}
	void Hankel(tw::Int modes, std::valarray<tw::Float>& matrix) {
		Field::Hankel(Rng(0), modes, matrix);
	}
	void InverseHankel(tw::Int modes, std::valarray<tw::Float>& matrix) {
		Field::InverseHankel(Rng(0), modes, matrix);
	}

	// Boundaries and messages

	void SetBoundaryConditions(const tw::grid::axis& axis, tw::bc::fld low, tw::bc::fld high) {
		Field::SetBoundaryConditions(Rng(0), axis, low, high);
	}
	void ApplyFoldingCondition() {
		Field::ApplyFoldingCondition(Rng(0));
	}
	void ApplyBoundaryCondition(bool homogeneous = true) {
		Field::ApplyBoundaryCondition(Rng(0),homogeneous);
	}
	template <class T>
	void AdjustTridiagonalForBoundaries(const tw::grid::axis& axis, const tw::grid::side& side, std::valarray<T>& T1, std::valarray<T>& T2, std::valarray<T>& T3, std::valarray<T>& source, T val) {
		Field::AdjustTridiagonalForBoundaries(Rng(0), axis, side, T1, T2, T3, source, val);
	}
	void ZeroGhostCells() {
		Field::ZeroGhostCells(Rng(0));
	}
	void DownwardCopy(const tw::grid::axis& axis, tw::Int cells) {
		Field::DownwardCopy(Rng(0), axis, cells);
	}
	void UpwardCopy(const tw::grid::axis& axis, tw::Int cells) {
		Field::UpwardCopy(Rng(0), axis, cells);
	}
	void DownwardDeposit(const tw::grid::axis& axis, tw::Int cells) {
		Field::DownwardDeposit(Rng(0), axis, cells);
	}
	void UpwardDeposit(const tw::grid::axis& axis, tw::Int cells) {
		Field::UpwardDeposit(Rng(0), axis, cells);
	}
	void CopyFromNeighbors() {
		Field::CopyFromNeighbors(Rng(0));
	}
	void DepositFromNeighbors() {
		Field::DepositFromNeighbors(Rng(0));
	}

	// Access by value

	tw::Float operator () (const tw::Int& i, const tw::Int& j, const tw::Int& k) const
	{
		return Field::operator () (1,i,j,k,0);
	}
	tw::Float operator () (const tw::Int& n,const tw::Int& i, const tw::Int& j, const tw::Int& k) const
	{
		return Field::operator () (n,i,j,k,0);
	}
	tw::Float operator () (const tw::cell& cell) const
	{
		return Field::operator () (cell,0);
	}
	tw::Float operator () (const tw::strip& strip, const tw::Int& s) const
	{
		return Field::operator () (strip,s,0);
	}
	tw::Float operator () (const tw::xstrip<1>& v, const tw::Int& i) const
	{
		return Field::operator () (v,i,0);
	}
	tw::Float operator () (const tw::xstrip<3>& v, const tw::Int& k) const
	{
		return Field::operator () (v,k,0);
	}

	// Access by reference (only possible for ScalarField)

	tw::Float& operator () (const tw::Int& i, const tw::Int& j, const tw::Int& k)
	{
		return Field::operator () (0,i,j,k,0);
	}
	tw::Float& operator () (const tw::Int& n,const tw::Int& i, const tw::Int& j, const tw::Int& k)
	{
		return Field::operator () (n,i,j,k,0);
	}
	tw::Float& operator () (const tw::cell& cell)
	{
		return Field::operator () (cell,0);
	}
	tw::Float& operator () (const tw::strip& strip, const tw::Int& s)
	{
		return Field::operator () (strip,s,0);
	}
	tw::Float& operator () (const tw::xstrip<1>& v, const tw::Int& i)
	{
		return Field::operator () (v,i,0);
	}
	tw::Float& operator () (const tw::xstrip<3>& v, const tw::Int& k)
	{
		return Field::operator () (v,k,0);
	}

	// Copy assignment

	ScalarField& operator = (ScalarField& A)
	{
		for (auto cell : EntireCellRange(A,0))
			(*this)(cell) = A(cell);
		return *this;
	}
	ScalarField& operator += (ScalarField& A)
	{
		for (auto cell : EntireCellRange(A,0))
			(*this)(cell) += A(cell);
		return *this;
	}
	ScalarField& operator -= (ScalarField& A)
	{
		for (auto cell : EntireCellRange(A,0))
			(*this)(cell) -= A(cell);
		return *this;
	}
	ScalarField& operator *= (ScalarField& A)
	{
		for (auto cell : EntireCellRange(A,0))
			(*this)(cell) *= A(cell);
		return *this;
	}

	// Fill assignment

	ScalarField& operator = (tw::Float a)
	{
		for (auto cell : EntireCellRange(*this,0))
			(*this)(cell) = a;
		return *this;
	}
	ScalarField& operator += (tw::Float a)
	{
		for (auto cell : EntireCellRange(*this,0))
			(*this)(cell) += a;
		return *this;
	}
	ScalarField& operator -= (tw::Float a)
	{
		for (auto cell : EntireCellRange(*this,0))
			(*this)(cell) -= a;
		return *this;
	}
	ScalarField& operator *= (tw::Float a)
	{
		for (auto cell : EntireCellRange(*this,0))
			(*this)(cell) *= a;
		return *this;
	}
};

export struct ComplexField: Field
{
	private:
	// hide the internal dimension initializer
	void Initialize(const tw::Int& components,const StaticSpace& ss,Task *task) {}

	public:
	// force the internal dimension
	void Initialize(const StaticSpace& ss,Task *task) {
		Field::Initialize(2,ss,task);
	}

	// integral transform

	void FFT(const DynSpace& ds) {
		Field::ComplexFFT(Rng(0,2),ds);
	}
	void InverseFFT(const DynSpace& ds) {
		Field::ComplexInverseFFT(Rng(0,2),ds);
	}
	void TransverseFFT(const DynSpace& ds) {
		Field::ComplexTransverseFFT(Rng(0,2),ds);
	}
	void InverseTransverseFFT(const DynSpace& ds) {
		Field::ComplexInverseTransverseFFT(Rng(0,2),ds);
	}
	void Hankel(tw::Int modes, std::valarray<tw::Float>& matrix) {
		Field::Hankel(Rng(0,2), modes, matrix);
	}
	void InverseHankel(tw::Int modes, std::valarray<tw::Float>& matrix) {
		Field::InverseHankel(Rng(0,2), modes, matrix);
	}

	// Boundaries and messages

	void SetBoundaryConditions(const tw::grid::axis& axis, tw::bc::fld low, tw::bc::fld high) {
		Field::SetBoundaryConditions(Rng(0,2), axis, low, high);
	}
	void ApplyFoldingCondition() {
		Field::ApplyFoldingCondition(Rng(0,2));
	}
	void ApplyBoundaryCondition(bool homogeneous = true) {
		Field::ApplyBoundaryCondition(Rng(0,2),homogeneous);
	}
	template <class T>
	void AdjustTridiagonalForBoundaries(const tw::grid::axis& axis, const tw::grid::side& side, std::valarray<T>& T1, std::valarray<T>& T2, std::valarray<T>& T3, std::valarray<T>& source, T val) {
		Field::AdjustTridiagonalForBoundaries(Rng(0,2), axis, side, T1, T2, T3, source, val);
	}
	void ZeroGhostCells() {
		Field::ZeroGhostCells(Rng(0,2));
	}
	void DownwardCopy(const tw::grid::axis& axis, tw::Int cells) {
		Field::DownwardCopy(Rng(0,2), axis, cells);
	}
	void UpwardCopy(const tw::grid::axis& axis, tw::Int cells) {
		Field::UpwardCopy(Rng(0,2), axis, cells);
	}
	void DownwardDeposit(const tw::grid::axis& axis, tw::Int cells) {
		Field::DownwardDeposit(Rng(0,2), axis, cells);
	}
	void UpwardDeposit(const tw::grid::axis& axis, tw::Int cells) {
		Field::UpwardDeposit(Rng(0,2), axis, cells);
	}
	void CopyFromNeighbors() {
		Field::CopyFromNeighbors(Rng(0,2));
	}
	void DepositFromNeighbors() {
		Field::DepositFromNeighbors(Rng(0,2));
	}

	// mirror superclass reference access

	tw::Float& operator () (const tw::Int& n,const tw::Int& i, const tw::Int& j, const tw::Int& k, const tw::Int& c) {
		return Field::operator() (n,i,j,k,c);
	}
	tw::Float& operator () (const tw::cell& cell, const tw::Int& c) {
		return Field::operator () (cell,c);
	}
	tw::Float& operator()(const tw::strip &strip, const tw::Int &s, const tw::Int& c) {
		return Field::operator()(strip, s, c);
	}
	tw::Float& operator()(const tw::xstrip<1> &v, const tw::Int &i, const tw::Int& c) {
		return Field::operator()(v, i, c);
	}
	tw::Float& operator()(const tw::xstrip<3> &v, const tw::Int &k, const tw::Int& c) {
		return Field::operator()(v, k, c);
	}

	// mirror superclass value access

	tw::Float operator () (const tw::Int& n,const tw::Int& i, const tw::Int& j, const tw::Int& k, const tw::Int& c) const {
		return Field::operator() (n,i,j,k,c);
	}
	tw::Float operator () (const tw::cell& cell, const tw::Int& c) const {
		return Field::operator () (cell,c);
	}
	tw::Float operator()(const tw::strip &strip, const tw::Int &s, const tw::Int& c) const {
		return Field::operator()(strip, s, c);
	}
	tw::Float operator()(const tw::xstrip<1> &v, const tw::Int &i, const tw::Int& c) const {
		return Field::operator()(v, i, c);
	}
	tw::Float operator()(const tw::xstrip<3> &v, const tw::Int &k, const tw::Int& c) const {
		return Field::operator()(v, k, c);
	}

	// Access by value

	tw::Complex operator () (const tw::Int& i, const tw::Int& j, const tw::Int& k) const {
		return tw::Complex(Field::operator () (1,i,j,k,0), Field::operator () (0,i,j,k,1));
	}
	tw::Complex operator () (const tw::Int& n,const tw::Int& i, const tw::Int& j, const tw::Int& k) const {
		return tw::Complex(Field::operator () (n,i,j,k,0), Field::operator () (n,i,j,k,1));
	}
	tw::Complex operator () (const tw::cell& cell) const {
		return tw::Complex(Field::operator () (cell,0), Field::operator () (cell,1));
	}
	tw::Complex operator()(const tw::strip &strip, const tw::Int &s) const {
		return tw::Complex(Field::operator()(strip, s, 0), Field::operator()(strip, s, 1));
	}
	tw::Complex operator()(const tw::xstrip<1> &v, const tw::Int &i) const {
		return tw::Complex(Field::operator()(v, i, 0), Field::operator()(v, i, 1));
	}
	tw::Complex operator()(const tw::xstrip<3> &v, const tw::Int &k) const {
		return tw::Complex(Field::operator()(v, k, 0), Field::operator()(v, k, 1));
	}

	// copy assignment

	ComplexField& operator = (ComplexField& A) {
		this->array = A.array;
		return *this;
	}

	// fill assignment

	ComplexField& operator = (tw::Complex a) {
		for (auto cell : EntireCellRange(*this,1)) {
			(*this)(cell,0) = a.real();
			(*this)(cell,1) = a.imag();
		}
		return *this;
	}
};

export struct Vec3Field: Field
{
	private:
	// hide the internal dimension initializer
	void Initialize(const tw::Int& components,const StaticSpace& ss,Task *task) {}

	public:
	// force the internal dimension
	void Initialize(const StaticSpace& ss,Task *task) {
		Field::Initialize(3,ss,task);
	}

	// Boundaries and messages

	void ApplyBoundaryCondition(bool homogeneous = true) {
		Field::ApplyBoundaryCondition(Rng(0,3),homogeneous);
	}
	void DownwardCopy(const tw::grid::axis& axis, tw::Int cells) {
		Field::DownwardCopy(Rng(0,3), axis, cells);
	}
	void UpwardCopy(const tw::grid::axis& axis, tw::Int cells) {
		Field::UpwardCopy(Rng(0,3), axis, cells);
	}
	void CopyFromNeighbors() {
		Field::CopyFromNeighbors(Rng(0,3));
	}
	void DepositFromNeighbors() {
		Field::DepositFromNeighbors(Rng(0,3));
	}

	// mirror superclass reference access

	tw::Float& operator () (const tw::Int& n,const tw::Int& i, const tw::Int& j, const tw::Int& k, const tw::Int& c) {
		return Field::operator() (n,i,j,k,c);
	}
	tw::Float& operator () (const tw::cell& cell, const tw::Int& c) {
		return Field::operator () (cell,c);
	}
	tw::Float& operator()(const tw::strip &strip, const tw::Int &s, const tw::Int& c) {
		return Field::operator()(strip, s, c);
	}
	tw::Float& operator()(const tw::xstrip<1> &v, const tw::Int &i, const tw::Int& c) {
		return Field::operator()(v, i, c);
	}
	tw::Float& operator()(const tw::xstrip<3> &v, const tw::Int &k, const tw::Int& c) {
		return Field::operator()(v, k, c);
	}

	// mirror superclass value access

	tw::Float operator () (const tw::Int& n,const tw::Int& i, const tw::Int& j, const tw::Int& k, const tw::Int& c) const {
		return Field::operator() (n,i,j,k,c);
	}
	tw::Float operator () (const tw::cell& cell, const tw::Int& c) const {
		return Field::operator () (cell,c);
	}
	tw::Float operator()(const tw::strip &strip, const tw::Int &s, const tw::Int& c) const {
		return Field::operator()(strip, s, c);
	}
	tw::Float operator()(const tw::xstrip<1> &v, const tw::Int &i, const tw::Int& c) const {
		return Field::operator()(v, i, c);
	}
	tw::Float operator()(const tw::xstrip<3> &v, const tw::Int &k, const tw::Int& c) const {
		return Field::operator()(v, k, c);
	}

	// Access by value

	tw::vec3 operator () (const tw::Int& i, const tw::Int& j, const tw::Int& k) const {
		return tw::vec3((*this)(1,i,j,k,0), (*this)(0,i,j,k,1), (*this)(0,i,j,k,2));
	}
	tw::vec3 operator () (const tw::Int& n,const tw::Int& i, const tw::Int& j, const tw::Int& k) const {
		return tw::vec3((*this)(n,i,j,k,0), (*this)(n,i,j,k,1), (*this)(n,i,j,k,2));
	}
	tw::vec3 operator () (const tw::cell& cell) const {
		return tw::vec3((*this)(cell,0), (*this)(cell,1), (*this)(cell,2));
	}
	tw::vec3 operator()(const tw::strip &strip, const tw::Int &s) const {
		return tw::vec3((*this)(strip,s,0), (*this)(strip,s,1), (*this)(strip,s,2));
	}
	tw::vec3 operator()(const tw::xstrip<1> &v, const tw::Int &i) const {
		return tw::vec3((*this)(v,i,0), (*this)(v,i,1), (*this)(v,i,2));
	}
	tw::vec3 operator()(const tw::xstrip<3> &v, const tw::Int &k) const {
		return tw::vec3((*this)(v,k,0), (*this)(v,k,1), (*this)(v,k,2));
	}
};
