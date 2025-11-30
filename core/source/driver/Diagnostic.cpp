module;

#include <tree_sitter/api.h>
#include "tw_includes.h"

export module driver:diagnostic;
export import :region;
export import :engine;


import input;
import fields;

export struct Diagnostic : Engine
{
	std::string filename;
	tw::Int skip[4];
	tw::Float t,tRef,t0,t1,timePeriod,gammaBoost;
	tw::vec3 vGalileo;
	bool headerWritten;
	const MetricSpace *ms; // may point at usual `space` or to an alternate

	Diagnostic(const std::string& name,MetricSpace *ms,Task *tsk) : Engine(name,ms,tsk) {
		skip[0] = 0;
		skip[1] = skip[2] = skip[3] = 1;
		t = 0.0;
		t0 = 0.0;
		t1 = tw::big_pos;
		tRef = tw::big_neg;
		timePeriod = 0.0;
		filename = "diagnostic";
		vGalileo = 0.0;
		gammaBoost = 1.0;
		headerWritten = false;
		directives.Add("period",new tw::input::Int(&skip[0]),false);
		directives.Add("time period",new tw::input::Float(&timePeriod),false);
		directives.Add("skip",new tw::input::Numbers<tw::Int>(&skip[1],3),false);
		directives.Add("filename",new tw::input::String(&filename),false);
		directives.Add("t0",new tw::input::Float(&t0),false);
		directives.Add("t1",new tw::input::Float(&t1),false);
		directives.Add("galilean velocity",new tw::input::Vec3(&vGalileo),false);
		directives.Add("boosted frame gamma",new tw::input::Float(&gammaBoost),false);
	}
	/// @brief copy base class parameters except for MetricSpace and filename
	void CopyParams(const Diagnostic& src) {
		for (auto i=0;i<4;i++)
			this->skip[i] = src.skip[i];
		this->t = src.t;
		this->tRef = src.tRef;
		this->t0 = src.t0;
		this->t1 = src.t1;
		this->timePeriod = src.timePeriod;
		this->gammaBoost = src.gammaBoost;
		this->vGalileo = src.vGalileo;
	}
	bool WriteThisStep() {
		tw::Float elapsedTime = space->WindowPos(0);
		tw::Float dt = space->dx(0);
		tw::Int stepNow = space->StepNow();
		t = elapsedTime; // save for use in reports

		if (timePeriod!=0.0)
		{
			if (tRef==elapsedTime)
				return true;
			if (elapsedTime >= t0 && elapsedTime <= t1)
			{
				tRef = elapsedTime;
				t0 += timePeriod;
				return true;
			}
			return false;
		}

		if (skip[0]!=0)
		{
			tw::Int startStep = tw::Int(t0/dt);
			if ( (elapsedTime>=t0) && (elapsedTime<=t1) && ((stepNow - startStep) % skip[0] == 0) )
				return true;
			else
				return false;
		}

		return false;
	}
	void StartGridFile(std::ofstream& grid) {
		std::string xname;
		if (filename=="tw::none")
			xname = "grid_warp.txt";
		else
			xname = filename + "_grid_warp.txt";
		if (headerWritten)
			grid.open(xname.c_str(),std::ios::app);
		else
			grid.open(xname.c_str());
		grid << "t = " << t << std::endl;
	}
	virtual void Start(const MetricSpace *alt = NULL) {
		this->ms = alt == NULL ? space : alt;
	}
	virtual void Finish() {;}
	virtual void ReportNumber(const std::string& label,tw::Float val,bool avg) {;}
	virtual void ReportField(const std::string& fieldName,const Field& F,const tw::Int n,const tw::Int c,
		const tw::dims unit = tw::dims::none,const std::string& pretty = "tw::none") {;}
	virtual void ReportParticle(const Particle& par,tw::Float m0) {;}
	virtual tw::Float VolumeIntegral(const std::string& fieldName,const Field& F,const tw::Int n,const tw::Int c) {
		tw::Float ans = 0.0;
		tw::Int loc[6];
		theRgn->GetLocalCellBounds(loc);
		for (tw::Int k=loc[4];k<=loc[5];k++)
			for (tw::Int j=loc[2];j<=loc[3];j++)
				for (tw::Int i=loc[0];i<=loc[1];i++)
					if (theRgn->Inside(space->Pos(i,j,k),0))
						ans += F(n,i,j,k,c) * space->dS(i,j,k,0);
		return ans;
	}
	virtual tw::Float FirstMoment(const std::string& fieldName,const Field& F,
		const tw::Int n,const tw::Int c,const tw::vec3& r0,const tw::grid::axis axis) {
		tw::Float ans = 0.0;
		const tw::Int ax = tw::grid::naxis(axis);
		tw::Int loc[6];
		theRgn->GetLocalCellBounds(loc);
		for (tw::Int k=loc[4];k<=loc[5];k++)
			for (tw::Int j=loc[2];j<=loc[3];j++)
				for (tw::Int i=loc[0];i<=loc[1];i++)
				{
					const tw::vec3 pos = space->Pos(i,j,k);
					if (theRgn->Inside(pos,0))
					{
						tw::vec3 r1 = r0;
						tw::vec3 r2 = pos;
						space->CurvilinearToCartesian(&r1);
						space->CurvilinearToCartesian(&r2);
						ans += F(n,i,j,k,c) * (r2[ax-1]-r1[ax-1]) * space->dS(i,j,k,0);
					}
				}
		return ans;
	}
	virtual void ReadCheckpoint(std::ifstream& inFile) {
		Engine::ReadCheckpoint(inFile);
		inFile.read((char *)&timePeriod,sizeof(tw::Float));
		inFile.read((char *)&t,sizeof(tw::Float));
		inFile.read((char *)&tRef,sizeof(tw::Float));
		inFile.read((char *)&t0,sizeof(tw::Float));
		inFile.read((char *)&t1,sizeof(tw::Float));
		inFile.read((char *)&headerWritten,sizeof(headerWritten));
	}
	virtual void WriteCheckpoint(std::ofstream& outFile) {
		Engine::WriteCheckpoint(outFile);
		outFile.write((char *)&timePeriod,sizeof(tw::Float));
		outFile.write((char *)&t,sizeof(tw::Float));
		outFile.write((char *)&tRef,sizeof(tw::Float));
		outFile.write((char *)&t0,sizeof(tw::Float));
		outFile.write((char *)&t1,sizeof(tw::Float));
		outFile.write((char *)&headerWritten,sizeof(headerWritten));
	}
};
