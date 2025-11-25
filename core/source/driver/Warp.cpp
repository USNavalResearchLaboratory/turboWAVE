module;

#include "tw_includes.h"

export module driver:warp;
import :tool;

import functions;

export struct Warp : ComputeTool,warp_base
{
	bool increasing;
	tw::Int rng[2];
	tw::Float L,gridSum;

	Warp(const std::string& name,MetricSpace *m,Task *tsk);
	virtual void Initialize();
	virtual void ReadCheckpoint(std::ifstream& inFile);
	virtual void WriteCheckpoint(std::ofstream& outFile);
	tw::Float AddedCellWidth(tw::Int globalCell);
};

Warp::Warp(const std::string& name,MetricSpace *m,Task *tsk) : ComputeTool(name,m,tsk)
{
	ax = tw::grid::z;
	increasing = true;
	directives.Add("axis",new tw::input::Enums<tw::grid::axis>(tw::grid::axis_map(),&ax));
	directives.Add("increasing",new tw::input::Bool(&increasing));
	directives.Add("index range",new tw::input::Numbers<tw::Int>(&rng[0],2));
	directives.Add("length",new tw::input::Float(&L));
}

void Warp::Initialize()
{
	const tw::Int N = rng[1] - rng[0] + 1;
	gridSum = 0.0;
	for (tw::Int i=1;i<=N;i++)
		gridSum += QuinticRise(tw::Float(i-1)/tw::Float(N-1));
}

tw::Float Warp::AddedCellWidth(tw::Int globalCell)
{
	const tw::Int N = rng[1] - rng[0] + 1;
	const tw::Float h = space->dx(tw::grid::naxis(ax));
	const tw::Float A = (1.0/gridSum)*(L/h - N);
	if (globalCell>=rng[0] && globalCell<=rng[1])
	{
		if (increasing)
			return h*A*QuinticRise(tw::Float(globalCell-rng[0])/tw::Float(N-1));
		else
			return h*A*QuinticFall(tw::Float(globalCell-rng[0])/tw::Float(N-1));
	}
	else
		return 0.0;
}

void Warp::ReadCheckpoint(std::ifstream& inFile)
{
	ComputeTool::ReadCheckpoint(inFile);
	inFile.read((char*)&ax,sizeof(ax));
	inFile.read((char*)&increasing,sizeof(increasing));
	inFile.read((char*)&rng[0],sizeof(rng));
	inFile.read((char*)&L,sizeof(L));
	inFile.read((char*)&gridSum,sizeof(gridSum));
}

void Warp::WriteCheckpoint(std::ofstream& outFile)
{
	ComputeTool::WriteCheckpoint(outFile);
	outFile.write((char*)&ax,sizeof(ax));
	outFile.write((char*)&increasing,sizeof(increasing));
	outFile.write((char*)&rng[0],sizeof(rng));
	outFile.write((char*)&L,sizeof(L));
	outFile.write((char*)&gridSum,sizeof(gridSum));
}
