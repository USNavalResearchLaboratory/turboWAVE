module;

#include "tw_includes.h"

export module fields:tools;
import :base;

import base;
import input;
import region;
import compute_tool;

/// Derivative of ComputeTool providing basic access to Field boundary conditions
export struct BoundedTool : ComputeTool
{
	tw::bc::fld x0,x1,y0,y1,z0,z1;
	tw::bc::fld x0s,x1s,y0s,y1s,z0s,z1s; // saved BC's

	BoundedTool(const std::string& name,MetricSpace *ms,Task *tsk);
	void SetBoundaryConditions(tw::bc::fld x0,tw::bc::fld x1,tw::bc::fld y0,tw::bc::fld y1,tw::bc::fld z0,tw::bc::fld z1);
	void SaveBoundaryConditions();
	void RestoreBoundaryConditions();
	void SetFieldsBoundaryConditions(Field& F,const Rng& r);
};

BoundedTool::BoundedTool(const std::string& name,MetricSpace *m,Task *tsk) : ComputeTool(name,m,tsk)
{
	x0 = x1 = y0 = y1 = z0 = z1 = tw::bc::fld::neumannWall;
	directives.Add("xboundary",new tw::input::Enums<tw::bc::fld>(tw::bc::fld_map(),&x0,&x1),false);
	directives.Add("yboundary",new tw::input::Enums<tw::bc::fld>(tw::bc::fld_map(),&y0,&y1),false);
	directives.Add("zboundary",new tw::input::Enums<tw::bc::fld>(tw::bc::fld_map(),&z0,&z1),false);
}

void BoundedTool::SetBoundaryConditions(tw::bc::fld x0,tw::bc::fld x1,tw::bc::fld y0,tw::bc::fld y1,tw::bc::fld z0,tw::bc::fld z1)
{
	this->x0 = x0;
	this->y0 = y0;
	this->z0 = z0;
	this->x1 = x1;
	this->y1 = y1;
	this->z1 = z1;
}

void BoundedTool::SaveBoundaryConditions()
{
	x0s = x0;
	x1s = x1;
	y0s = y0;
	y1s = y1;
	z0s = z0;
	z1s = z1;
}

void BoundedTool::RestoreBoundaryConditions()
{
	x0 = x0s;
	x1 = x1s;
	y0 = y0s;
	y1 = y1s;
	z0 = z0s;
	z1 = z1s;
}

void BoundedTool::SetFieldsBoundaryConditions(Field& F,const Rng& r)
{
	F.SetBoundaryConditions(r,tw::grid::x,x0,x1);
	F.SetBoundaryConditions(r,tw::grid::y,y0,y1);
	F.SetBoundaryConditions(r,tw::grid::z,z0,z1);
}
