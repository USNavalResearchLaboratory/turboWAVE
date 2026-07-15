module;

#include <tree_sitter/api.h>
#include "tw_includes.h"

export module driver:diagnostic;
export import :region;
export import :engine;


import input;
import fields;

export std::string grid_file_name(const std::string& diag_name,tw::Int variant) {
	std::string xname("");
	if (diag_name!="tw::none") {
		xname = diag_name + "_";
	}
	xname += "grid_";
	if (variant > 0) {
		xname += std::to_string(variant) + "_";
	}
	xname += "warp.txt";
	return xname;
}

export struct Diagnostic : Engine
{
	std::shared_ptr<Region> theRgn;
	std::string filename;
	tw::Int skip[4];
	tw::Float t,tRef,t0,t1,timePeriod;
	tw::vec3 vGalileo,boost3;
	tw::vec4 boost;
	bool headerWritten;
	std::vector<MetricSpace*> variants;
	tw::Int variant;

	Diagnostic(const std::string& name,MetricSpace *ms,Task *tsk) : Engine(name,ms,tsk) {
		variants.push_back(ms);
		variant = 0;
		skip[0] = 0;
		skip[1] = skip[2] = skip[3] = 1;
		t = 0.0;
		t0 = 0.0;
		t1 = tw::max_pos;
		tRef = tw::max_neg;
		timePeriod = 0.0;
		filename = "diagnostic";
		vGalileo = 0.0;
		boost3 = tw::vec3(0,0,0);
		boost = tw::vec4(1,0,0,0);
		headerWritten = false;
		directives.Add("period",new tw::input::Int(&skip[0]),false);
		directives.Add("time period",new tw::input::Float(&timePeriod),false);
		directives.Add("skip",new tw::input::Numbers<tw::Int>(&skip[1],3),false);
		directives.Add("filename",new tw::input::String(&filename),false);
		directives.Add("t0",new tw::input::Float(&t0),false);
		directives.Add("t1",new tw::input::Float(&t1),false);
		directives.Add("galilean velocity",new tw::input::Vec3(&vGalileo),false);
		directives.Add("boost",new tw::input::Vec3(&boost3),false);
	}
	virtual void Initialize() {
		boost = tw::vec4(std::sqrt(1+Norm(boost3)),boost3);
		if (region.use_count()==0) {
			theRgn = std::make_shared<SimpleRegion>("default_entire",space,task,std::make_unique<EntireRegion>("entire",space,task));
		} else {
			theRgn = std::dynamic_pointer_cast<Region>(region);
			if (!theRgn) {
				throw tw::FatalError("dynamic cast to Region failed");
			}
		}
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
	void StartGridFile(std::ofstream& grid,tw::Int variant) {
		auto xname = grid_file_name(filename,variant);
		if (headerWritten)
			grid.open(xname.c_str(),std::ios::app);
		else
			grid.open(xname.c_str());
		grid << "t = " << t << std::endl;
	}
	/// Add a variant on the simulation's metric space to be used with this diagnostic.
	/// Returns the id number of this variant.
	/// Variant id = 0 will always be the simulation's primary metric space.
	tw::Int AddVariant(MetricSpace *ms) {
		variants.push_back(ms);
		return variants.size() - 1;
	}
	void SwitchVariant(tw::Int id) {
		this->variant = id;
		this->space = variants[id];
		this->theRgn->space = variants[id];
	}
	virtual void Start() {;}
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
					if (theRgn->Inside(space->Pos4(1,i,j,k),0))
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
					const auto pos = space->Pos4(1,i,j,k);
					if (theRgn->Inside(pos,0))
					{
						tw::vec3 r1 = r0;
						tw::vec3 r2 = pos.spatial();
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
