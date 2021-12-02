#include "simulation.h"
#include <iomanip>

////////////////////////////////////////////////////////////////
// Dictionary for storing metadata associated with .npy files //
////////////////////////////////////////////////////////////////

meta_writer::meta_writer(const tw::UnitConverter& units)
{
	native = units;
	mks = tw::UnitConverter(tw::units::mks,native);
	cgs = tw::UnitConverter(tw::units::cgs,native);
	plasma = tw::UnitConverter(tw::units::plasma,native);
	atomic = tw::UnitConverter(tw::units::atomic,native);
	natural = tw::UnitConverter(tw::units::natural,native);
}

std::string meta_writer::s(const std::string& raw)
{
	std::stringstream ans;
	ans << "\"";
	for (auto c : raw)
		if (c=='\\')
			ans << "\\\\";
		else
			ans << c;
	ans << "\"";
	return ans.str();
}

void meta_writer::start_entry(const std::string& name,const std::string& diagnostic_name)
{
	std::fstream outFile("tw_metadata.json",std::ios::out | std::ios::in); // we have to use input and output mode to avoid default truncation
	std::map<tw::units,std::string> m = tw::get_unit_map_r();
	// first back up and replace the closing brace with a comma
	outFile.seekp(-1,std::ios::end);
	outFile << "," << std::endl;
	outFile << s(name) << ": {" << std::endl;
	outFile << "\t\"native units\": " << s(m[native.native]) << "," << std::endl;
	if (diagnostic_name=="tw::none")
		outFile << "\t\"grid\": " << "\"grid_warp.txt\"," << std::endl;
	else
		outFile << "\t\"grid\": " << s(diagnostic_name + "_grid_warp.txt") << "," << std::endl;
	outFile << "\t\"axes\": " << "{" << std::endl;
	outFile.close();
}

void meta_writer::define_axis(const std::string& name,tw::Int ax,const std::string& label,tw::dims d,bool last)
{
	std::ofstream outFile("tw_metadata.json",std::ios::app);
	outFile << "\t\t\"" << ax << "\": {" << std::endl;
	outFile << "\t\t\t" << "\"label\" : " << s(label) << "," << std::endl;
	outFile << "\t\t\t" << "\"mks label\" : " << s(tw::mks_label(d)) << "," << std::endl;
	outFile << "\t\t\t" << "\"mks multiplier\" : " << (1.0*d>>native>>mks) << "," << std::endl;
	outFile << "\t\t\t" << "\"cgs label\" : " << s(tw::cgs_label(d)) << "," << std::endl;
	outFile << "\t\t\t" << "\"cgs multiplier\" : " << (1.0*d>>native>>cgs) << "," << std::endl;
	outFile << "\t\t\t" << "\"plasma label\" : " << s(tw::plasma_label(d)) << "," << std::endl;
	outFile << "\t\t\t" << "\"plasma multiplier\" : " << (1.0*d>>native>>plasma) << "," << std::endl;
	outFile << "\t\t\t" << "\"atomic label\" : " << s(tw::atomic_label(d)) << "," << std::endl;
	outFile << "\t\t\t" << "\"atomic multiplier\" : " << (1.0*d>>native>>atomic) << "," << std::endl;
	outFile << "\t\t\t" << "\"natural label\" : " << s(tw::natural_label(d)) << "," << std::endl;
	outFile << "\t\t\t" << "\"natural multiplier\" : " << (1.0*d>>native>>natural) << std::endl;
	outFile << "\t\t}" << (last==false ? "," : "") << std::endl;
	outFile.close();
}

std::string meta_writer::refine_label(const std::string& label,const tw::vec3& vGalileo,const tw::grid::geometry& geo)
{
	if (label=="$x$")
	{
		if (geo==tw::grid::cylindrical)
			return "$\\varrho$";
		if (geo==tw::grid::spherical)
			return "$r$";
	}
	if (label=="$y$" && geo!=tw::grid::cartesian)
		return "$\\varphi$";
	if (label=="$z$")
	{
		if (geo==tw::grid::spherical)
			return "$\\theta$";
		if (vGalileo.z==1.0)
			return "$(z-ct)$";
		if (vGalileo.z!=0.0 && vGalileo.z!=1.0)
			return "$\\zeta$";
	}
	return label;
}

void meta_writer::finish_entry()
{
	std::ofstream outFile("tw_metadata.json",std::ios::app);
	outFile << "\t}" << std::endl; // close axes dictionary
	outFile << "}}"; // close file and main dictionary, n.b. whitespace would spoil append strategy
}

//////////////////////////////////////////////////
// Writer class for four dimensional .npy files //
//////////////////////////////////////////////////

std::string npy_writer::form_header(tw::Int shape[4])
{
	// Form a string containing the header for a standard .npy file.
	// Shape elements are padded so we get the same length header every time.
	uint32_t HEADER_LEN;
	const size_t type_size = sizeof(HEADER_LEN);
	std::string vers("\x93NUMPY\x02");
	vers += '\0'; // null character cannot be put in the C-string above
	std::stringstream dict;
	dict << "{'descr': '<f4', 'fortran_order': False, 'shape': (";
	for (tw::Int i=0;i<4;i++)
		dict << std::setw(10) << std::setfill(' ') << shape[i] << ",";
	dict << "), }";
	HEADER_LEN = dict.str().size();
	while ((HEADER_LEN+8+type_size+1)%64)
	{
		HEADER_LEN++;
		dict << " ";
	}
	HEADER_LEN++;
	dict << "\n";
	// Form bytes for the little endian integer giving the header length
	char hbuff[type_size];
	memcpy(hbuff,(char*)&HEADER_LEN,type_size);
	BufferLittleEndian(hbuff,type_size);
	// Now put it all together in a string
	std::stringstream ans;
	ans << vers;
	for (tw::Int i=0;i<type_size;i++)
		ans << hbuff[i];
	ans << dict.str();
	return ans.str();
}

void npy_writer::write_header(const std::string& name,tw::Int shape[4])
{
	// Write the header for a standard .npy file.
	std::ofstream outFile;
	outFile.open(name.c_str(),std::ios::binary | std::ios::trunc);
	std::string header(form_header(shape));
	outFile.write(header.data(),header.size());
	outFile.close();
}

void npy_writer::update_shape(const std::string& name,tw::Int shape[4])
{
	// Updates the time dimension by evaluating the size of the current file.
	std::ifstream inFile;
	std::string header;
	inFile.open(name.c_str(),std::ios::binary);
	inFile.seekg(0,std::ios::end);
	size_t length = inFile.tellg();
	inFile.close();
	header = form_header(shape);
	shape[0] = (length - header.size())/(sizeof(float)*shape[1]*shape[2]*shape[3]);
	header = form_header(shape);
	// Now write the changes by simply overwriting the whole header
	std::fstream outFile; // we have to use input and output mode to avoid default truncation
	outFile.open(name.c_str(),std::ios::binary | std::ios::out | std::ios::in);
	outFile.seekp(0,std::ios::beg);
	outFile.write(header.data(),header.size());
	outFile.close();
}

void npy_writer::get_frame(const std::string& name,char *gData,tw::Int shape[4],tw::Int frame)
{
	std::ifstream inFile;
	std::string header(form_header(shape));
	size_t frameSize = sizeof(float)*shape[1]*shape[2]*shape[3];
	inFile.open(name.c_str(),std::ios::binary);
	inFile.seekg(header.size()+frame*frameSize,std::ios::beg);
	ReadLittleEndian(gData,frameSize,sizeof(float),inFile);
	inFile.close();
}

void npy_writer::set_frame(const std::string& name,const char *gData,tw::Int shape[4],tw::Int frame)
{
	std::fstream outFile;
	std::string header(form_header(shape));
	size_t frameSize = sizeof(float)*shape[1]*shape[2]*shape[3];
	outFile.open(name.c_str(),std::ios::binary | std::ios::out | std::ios::in);
	outFile.seekp(header.size()+frame*frameSize,std::ios::beg);
	WriteLittleEndian(gData,frameSize,sizeof(float),outFile);
	outFile.close();
}

void npy_writer::add_frame(const std::string& name,const char *gData,tw::Int shape[4])
{
	std::fstream outFile;
	outFile.open(name.c_str(),std::ios::binary | std::ios::out | std::ios::app);
	WriteLittleEndian(gData,sizeof(float)*shape[1]*shape[2]*shape[3],sizeof(float),outFile);
	outFile.close();
}

//////////////////////
// DIAGNOSTIC TOOLS //
//////////////////////


Diagnostic::Diagnostic(const std::string& name,MetricSpace *ms,Task *tsk) : ComputeTool(name,ms,tsk)
{
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

bool Diagnostic::WriteThisStep(tw::Float elapsedTime,tw::Float dt,tw::Int stepNow)
{
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

void Diagnostic::StartGridFile(std::ofstream& grid)
{
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

void Diagnostic::Start()
{
}

void Diagnostic::Finish()
{
}

void Diagnostic::Float(const std::string& label,tw::Float val,bool average)
{
}

tw::Float Diagnostic::VolumeIntegral(const std::string& fieldName,const struct Field& F,const tw::Int c)
{
	tw::Float ans = 0.0;
	tw::Int x0,x1,y0,y1,z0,z1;
	theRgn->GetLocalCellBounds(&x0,&x1,&y0,&y1,&z0,&z1);
	for (tw::Int k=z0;k<=z1;k++)
		for (tw::Int j=y0;j<=y1;j++)
			for (tw::Int i=x0;i<=x1;i++)
				if (theRgn->Inside(space->Pos(i,j,k),*space))
					ans += F(i,j,k,c) * space->dS(i,j,k,0);
	return ans;
}

tw::Float Diagnostic::FirstMoment(const std::string& fieldName,const struct Field& F,const tw::Int c,const tw::vec3& r0,const tw::grid::axis axis)
{
	tw::Float ans = 0.0;
	const tw::Int ax = tw::grid::naxis(axis);
	tw::Int x0,x1,y0,y1,z0,z1;
	theRgn->GetLocalCellBounds(&x0,&x1,&y0,&y1,&z0,&z1);
	for (tw::Int k=z0;k<=z1;k++)
		for (tw::Int j=y0;j<=y1;j++)
			for (tw::Int i=x0;i<=x1;i++)
			{
				const tw::vec3 pos = space->Pos(i,j,k);
				if (theRgn->Inside(pos,*space))
				{
					tw::vec3 r1 = r0;
					tw::vec3 r2 = pos;
					space->CurvilinearToCartesian(&r1);
					space->CurvilinearToCartesian(&r2);
					ans += F(i,j,k,c) * (r2[ax-1]-r1[ax-1]) * space->dS(i,j,k,0);
				}
			}
	return ans;
}

void Diagnostic::Field(const std::string& fieldName,const struct Field& F,const tw::Int c,tw::dims unit,const std::string& pretty)
{
}

void Diagnostic::Particle(const struct Particle& par,tw::Float m0,tw::Float tp)
{
}

void Diagnostic::ReadCheckpoint(std::ifstream& inFile)
{
	ComputeTool::ReadCheckpoint(inFile);
	inFile.read((char *)&timePeriod,sizeof(tw::Float));
	inFile.read((char *)&t,sizeof(tw::Float));
	inFile.read((char *)&tRef,sizeof(tw::Float));
	inFile.read((char *)&t0,sizeof(tw::Float));
	inFile.read((char *)&t1,sizeof(tw::Float));
	inFile.read((char *)&headerWritten,sizeof(headerWritten));
}

void Diagnostic::WriteCheckpoint(std::ofstream& outFile)
{
	ComputeTool::WriteCheckpoint(outFile);
	outFile.write((char *)&timePeriod,sizeof(tw::Float));
	outFile.write((char *)&t,sizeof(tw::Float));
	outFile.write((char *)&tRef,sizeof(tw::Float));
	outFile.write((char *)&t0,sizeof(tw::Float));
	outFile.write((char *)&t1,sizeof(tw::Float));
	outFile.write((char *)&headerWritten,sizeof(headerWritten));
}

TextTableBase::TextTableBase(const std::string& name,MetricSpace *ms,Task *tsk) : Diagnostic(name,ms,tsk)
{
	numSigFigs = 6;
	filename = "table";
	directives.Add("precision",new tw::input::Int(&numSigFigs),false);
}

void TextTableBase::Start()
{
	Diagnostic::Start();
	labels.clear();
	values.clear();
	avg.clear();
}

void TextTableBase::Finish()
{
	Diagnostic::Finish();
	if (values.size())
	{
		const tw::Int master = 0;
		const tw::Int curr = task->strip[0].Get_rank();
		const tw::Int numRanks = task->strip[0].Get_size();
		std::valarray<tw::Float> buff;
		buff.resize(values.size());
		for (tw::Int i=0;i<values.size();i++)
			buff[i] = values[i]/(avg[i]?numRanks:1.0);
		task->strip[0].Sum(&buff[0],&buff[0],values.size()*sizeof(tw::Float),master);
		for (tw::Int i=0;i<values.size();i++)
			values[i] = buff[i];
		if (curr==master)
		{
			std::string xname(filename+".txt");
			std::ofstream file;
			if (!headerWritten)
			{
				file.open(xname.c_str());
				for (auto s : labels)
					file << s << " ";
				file << std::endl;
				headerWritten = true;
			}
			else
				file.open(xname.c_str(),std::ios::app);
			file.precision(numSigFigs);
			for (auto v : values)
				file << v << " ";
			file << std::endl;
			file.close();
		}
	}
}

void TextTableBase::Float(const std::string& label,tw::Float val,bool average)
{
	Diagnostic::Float(label,val,average);
	labels.push_back(label);
	values.push_back(val);
	avg.push_back(average);
}

VolumeDiagnostic::VolumeDiagnostic(const std::string& name,MetricSpace *ms,Task *tsk) : TextTableBase(name,ms,tsk)
{
	filename = "energy";
}

tw::Float VolumeDiagnostic::VolumeIntegral(const std::string& fieldName,const struct Field& F,const tw::Int c)
{
	const tw::Float ans = Diagnostic::VolumeIntegral(fieldName,F,c);
	labels.push_back(fieldName);
	values.push_back(ans);
	avg.push_back(false);
	return ans;
}

tw::Float VolumeDiagnostic::FirstMoment(const std::string& fieldName,const struct Field& F,const tw::Int c,const tw::vec3& r0,const tw::grid::axis axis)
{
	const tw::Float ans = Diagnostic::FirstMoment(fieldName,F,c,r0,axis);
	labels.push_back(fieldName);
	values.push_back(ans);
	avg.push_back(false);
	return ans;
}

PointDiagnostic::PointDiagnostic(const std::string& name,MetricSpace *ms,Task *tsk) : TextTableBase(name,ms,tsk)
{
	thePoint = tw::vec3(1.0);
	filename = "pt1";
	directives.Add("point",new tw::input::Vec3(&thePoint));
}

void PointDiagnostic::Field(const std::string& fieldName,const struct Field& F,const tw::Int c,tw::dims unit,const std::string& pretty)
{
	tw::vec3 r = thePoint + vGalileo*t;
	if (space->IsPointValid(r)) // assumes uniform grid
	{
		std::valarray<tw::Float> ans(1);
		weights_3D w;
		space->GetWeights(&w,r);
		F.Interpolate(ans,Element(c),w);
		labels.push_back(fieldName);
		values.push_back(ans[0]);
		avg.push_back(false);
	}
	else
	{
		labels.push_back(fieldName);
		values.push_back(0.0);
		avg.push_back(false);
	}
}

void PointDiagnostic::ReadCheckpoint(std::ifstream& inFile)
{
	TextTableBase::ReadCheckpoint(inFile);
	inFile.read((char *)&thePoint,sizeof(tw::vec3));
}

void PointDiagnostic::WriteCheckpoint(std::ofstream& outFile)
{
	TextTableBase::WriteCheckpoint(outFile);
	outFile.write((char *)&thePoint,sizeof(tw::vec3));
}

BoxDiagnostic::BoxDiagnostic(const std::string& name,MetricSpace *ms,Task *tsk) : Diagnostic(name,ms,tsk)
{
	average = false;
	filename = "tw::none";
	directives.Add("average",new tw::input::Bool(&average),false);
	directives.Add("reports",new tw::input::List<std::vector<std::string>>(&reports),false);
}

void BoxDiagnostic::GetGlobalIndexing(tw::Int pts[4],tw::Int glb[6])
{
	// On output, pts has the number of cells, glb has the global index bounds
	// skip is corrected for ignorable dimensions
	theRgn->GetGlobalCellBounds(glb);

	for (tw::Int ax=1;ax<=3;ax++)
	{
		const tw::Int lb = 2*ax-2;
		const tw::Int ub = 2*ax-1;
		skip[ax] = space->Dim(ax)==1 ? 1 : skip[ax];
		// Force global bounds into skipping sequence
		// The upper bound is adjusted to respect the lower bound
		pts[ax] = (glb[ub] - glb[lb] + 1)/skip[ax];
		glb[ub] = glb[lb] + skip[ax]*(pts[ax]-1);
	}
}

void BoxDiagnostic::GetLocalIndexing(const tw::Int pts[4],const tw::Int glb[6],tw::Int loc[6],const tw::Int coords[4])
{
	// On output, pts has the number of cells, glb has the global index bounds, loc has the local index bounds
	// coords is an input with the cartesian MPI domain indices (not necessarily the domain of execution)
	tw::Int i0;
	for (tw::Int ax=1;ax<=3;ax++)
	{
		const tw::Int cell0 = coords[ax]*task->localCells[ax];
		const tw::Int lb = 2*ax-2;
		const tw::Int ub = 2*ax-1;
		// Get the lower local cell bounds and force into global skipping sequence
		i0 = 1 + cell0;
		i0 = i0 < glb[lb] ? glb[lb] : i0;
		i0 -= cell0 + ((i0 - glb[lb]) % skip[ax]);
		if (i0<1) i0 += skip[ax];
		loc[lb] = i0;
		// Get the upper local cell bounds and force into global skipping sequence
		i0 = task->localCells[ax] + cell0;
		i0 = i0 > glb[ub] ? glb[ub] : i0;
		i0 -= cell0 + ((i0 - glb[lb]) % skip[ax]);
		loc[ub] = i0;
	}
}

void BoxDiagnostic::Field(const std::string& fieldName,const struct Field& F,const tw::Int c,tw::dims unit,const std::string& pretty)
{
	if (reports.size()>0)
		if (std::find(reports.begin(),reports.end(),fieldName)==reports.end())
			return;

	tw::Int buffSize,ready,i0,i1;
	std::valarray<float> buffer,gData;
	std::string xname;
	meta_writer meta(native);
	npy_writer writer;
	tw::Int coords[4],pts[4],loc[6],glb[6],dim[4],s[4];
	tw::Int thisNode,curr,master;

	thisNode = task->strip[0].Get_rank();
	master = 0;

	if (filename=="tw::none")
		xname = fieldName + ".npy";
	else
		xname = filename + "_" + fieldName + ".npy";

	GetGlobalIndexing(pts,glb);
	for (tw::Int i=1;i<=3;i++) dim[i] = space->Dim(i);
	for (tw::Int i=1;i<=3;i++) s[i] = skip[i];

	if (thisNode==master)
		gData.resize(pts[1]*pts[2]*pts[3]);

	for (tw::Int kdom=(glb[4]-1)/dim[3];kdom<=(glb[5]-1)/dim[3];kdom++)
		for (tw::Int jdom=(glb[2]-1)/dim[2];jdom<=(glb[3]-1)/dim[2];jdom++)
			for (tw::Int idom=(glb[0]-1)/dim[1];idom<=(glb[1]-1)/dim[1];idom++)
			{
				curr = task->strip[0].Cart_rank(idom,jdom,kdom); // domain being written out
				coords[1] = idom; coords[2] = jdom; coords[3] = kdom;
				GetLocalIndexing(pts,glb,loc,coords);

				if (loc[0]<=dim[1] && loc[1]>=1 && loc[2]<=dim[2] && loc[3]>=1 && loc[4]<=dim[3] && loc[5]>=1)
				{
					if (curr==thisNode)
					{
						if (curr==master)
						{
							for (tw::Int k=loc[4];k<=loc[5];k+=s[3])
								for (tw::Int j=loc[2];j<=loc[3];j+=s[2])
									for (tw::Int i=loc[0];i<=loc[1];i+=s[1])
									{
										i0 = pts[2]*pts[3]*(task->GlobalCellIndex(i,1) - glb[0])/s[1];
										i0 += pts[3]*(task->GlobalCellIndex(j,2) - glb[2])/s[2];
										i0 += (task->GlobalCellIndex(k,3) - glb[4])/s[3];
										gData[i0] = F(i,j,k,c);
									}
						}
						else
						{
							task->strip[0].Recv(&ready,sizeof(tw::Int),master);
							buffSize = (loc[1]-loc[0]+s[1])*(loc[3]-loc[2]+s[2])*(loc[5]-loc[4]+s[3])/(s[1]*s[2]*s[3]);
							buffer.resize(buffSize);
							for (tw::Int k=loc[4];k<=loc[5];k+=s[3])
								for (tw::Int j=loc[2];j<=loc[3];j+=s[2])
									for (tw::Int i=loc[0];i<=loc[1];i+=s[1])
									{
										i0 = (loc[3]-loc[2]+s[2])*(loc[5]-loc[4]+s[3])*(i - loc[0])/(s[1]*s[2]*s[3]);
										i0 += (loc[5]-loc[4]+s[3])*(j - loc[2])/(s[2]*s[3]);
										i0 += (k - loc[4])/s[3];
										buffer[i0] = F(i,j,k,c);
									}
							task->strip[0].Send(&buffer[0],sizeof(float)*buffSize,master);
						}
					}
					if (thisNode==master && curr!=master)
					{
						task->strip[0].Send(&ready,sizeof(tw::Int),curr);
						buffSize = (loc[1]-loc[0]+s[1])*(loc[3]-loc[2]+s[2])*(loc[5]-loc[4]+s[3])/(s[1]*s[2]*s[3]);
						buffer.resize(buffSize);
						task->strip[0].Recv(&buffer[0],sizeof(float)*buffSize,curr);
						for (tw::Int k=loc[4];k<=loc[5];k+=s[3])
							for (tw::Int j=loc[2];j<=loc[3];j+=s[2])
								for (tw::Int i=loc[0];i<=loc[1];i+=s[1])
								{
									i0 = (loc[3]-loc[2]+s[2])*(loc[5]-loc[4]+s[3])*(i - loc[0])/(s[1]*s[2]*s[3]);
									i0 += (loc[5]-loc[4]+s[3])*(j - loc[2])/(s[2]*s[3]);
									i0 += (k - loc[4])/s[3];
									i1 = pts[2]*pts[3]*(coords[1]*task->localCells[1]+i - glb[0])/s[1];
									i1 += pts[3]*(coords[2]*task->localCells[2]+j - glb[2])/s[2];
									i1 += (coords[3]*task->localCells[3]+k - glb[4])/s[3];
									gData[i1] = buffer[i0];
								}
					}
				}
			}

	if (thisNode==master)
	{
		if (!headerWritten)
		{
			pts[0] = 0;
			writer.write_header(xname,pts);
			meta.start_entry(xname,filename);
			meta.define_axis(xname,0,"$t$",tw::dims::time);
			meta.define_axis(xname,1,meta.refine_label("$x$",vGalileo,space->geo),tw::dims::length);
			meta.define_axis(xname,2,meta.refine_label("$y$",vGalileo,space->geo),tw::dims::length);
			meta.define_axis(xname,3,meta.refine_label("$z$",vGalileo,space->geo),tw::dims::length);
			if (pretty=="tw::none")
				meta.define_axis(xname,4,fieldName,unit,true);
			else
				meta.define_axis(xname,4,pretty,unit,true);
			meta.finish_entry();
		}
		writer.add_frame(xname,(char*)&gData[0],pts);
		writer.update_shape(xname,pts);
	}
}

void BoxDiagnostic::Finish()
{
	// Write out the grid data

	std::ofstream gridFile;
	const tw::Int master = 0;
	const tw::Int curr_global = task->strip[0].Get_rank();
	tw::Int pts[4],glb[6];
	GetGlobalIndexing(pts,glb);
	if (curr_global==master)
		StartGridFile(gridFile);
	if (!headerWritten) // assuming static grid no need to keep writing spatial points
	{
		for (tw::Int ax=1;ax<=3;ax++)
		{
			// Message passing is needed to get the spatial points
			std::valarray<tw::Float> X(task->globalCells[ax]);
			const tw::Int offset = space->Dim(ax)*task->strip[ax].Get_rank();
			for (tw::Int i=1;i<=space->Dim(ax);i++)
				X[i-1+offset] = space->X(i,ax);
			task->strip[ax].Gather(&X[offset],&X[offset],task->localCells[ax]*sizeof(tw::Float),master);
			if (curr_global==master)
			{
				const tw::Int lb = glb[2*ax-2];
				const tw::Int ub = glb[2*ax-1];
				gridFile << "axis" << ax << " = ";
				for (tw::Int i=lb;i<ub;i+=skip[ax])
					gridFile << X[i-1] << " ";
				gridFile << X[ub-1] << std::endl;
			}
		}
	}
	if (curr_global==master)
		gridFile.close();
	headerWritten = true;
}

ParticleOrbits::ParticleOrbits(const std::string& name,MetricSpace *ms,Task *tsk) : Diagnostic(name,ms,tsk)
{
	filename = "par";
	minGamma = 1.0;
	directives.Add("minimum gamma",new tw::input::Float(&minGamma));
}

void ParticleOrbits::Start()
{
	Diagnostic::Start();
	parData.clear();
}

void ParticleOrbits::Finish()
{
	const tw::Int master = 0;
	std::valarray<float> parBuffer;
	parBuffer.resize(parData.size());
	for (tw::Int i=0;i<parData.size();i++)
		parBuffer[i] = parData[i];
	tw::Int pts = parData.size();
	if (task->strip[0].Get_rank()!=master)
	{
		task->strip[0].Send(&pts,sizeof(tw::Int),master);
		task->strip[0].Send(&parBuffer[0],sizeof(float)*pts,master);
	}
	if (task->strip[0].Get_rank()==master)
	{
		npy_writer writer;
		tw::Int shape[4] = { 0 , 1 , 1 , 8 };
		tw::Int false_shape[4] = { 0 , 1 , 1 , 8 };
		std::string xname = filename + ".npy";
		if (!headerWritten)
		{
			writer.write_header(xname,shape);
			headerWritten = true;
		}
		false_shape[2] = pts;
		writer.add_frame(xname,(char*)&parBuffer[0],false_shape);
		for (tw::Int i=0;i<task->strip[0].Get_size();i++)
		{
			if (i!=master)
			{
				task->strip[0].Recv(&pts,sizeof(tw::Int),i);
				parBuffer.resize(pts);
				task->strip[0].Recv(&parBuffer[0],sizeof(float)*pts,i);
				false_shape[2] = pts;
				writer.add_frame(xname,(char*)&parBuffer[0],false_shape);
			}
		}
		float separator[8] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
		writer.add_frame(xname,(char*)&separator[0],shape);
		writer.update_shape(xname,shape);
	}
}

void ParticleOrbits::Particle(const struct Particle& par,tw::Float m0,tw::Float tp)
{
	tw::vec4 x(tp,space->PositionFromPrimitive(par.q));
	tw::vec4 p(par.p);
	// Boosts will work in Cartesian or cylindrical, but not spherical.
	x.zBoost(gammaBoost,1.0);
	p.zBoost(gammaBoost,1.0);
	tw::vec3 x3 = x.spatial();
	if (theRgn->Inside(x3,*space))
		if (p[0] >= m0*minGamma)
		{
			space->CurvilinearToCartesian(&x3);
			x3 -= vGalileo*x[0];
			x = tw::vec4(x[0],x3);
			parData.push_back(x[1]);
			parData.push_back(p[1]);
			parData.push_back(x[2]);
			parData.push_back(p[2]);
			parData.push_back(x[3]);
			parData.push_back(p[3]);
			// TODO: support tags
			//parData.push_back(par.aux1);
			//parData.push_back(par.aux2);
		}
}

void ParticleOrbits::ReadCheckpoint(std::ifstream& inFile)
{
	Diagnostic::ReadCheckpoint(inFile);
	inFile.read((char *)&minGamma,sizeof(tw::Float));
}

void ParticleOrbits::WriteCheckpoint(std::ofstream& outFile)
{
	Diagnostic::WriteCheckpoint(outFile);
	outFile.write((char *)&minGamma,sizeof(tw::Float));
}

PhaseSpaceDiagnostic::PhaseSpaceDiagnostic(const std::string& name,MetricSpace *ms,Task *tsk) : Diagnostic(name,ms,tsk)
{
	filename = "xpx";
	ax[0] = tw::grid::t;
	ax[1] = tw::grid::x;
	ax[2] = tw::grid::px;
	ax[3] = tw::grid::py;
	bounds[0] = bounds[2] = bounds[4] = 0.0;
	bounds[1] = bounds[3] = bounds[5] = 1.0;
	dims[0] = 0;
	dims[1] = 100;
	dims[2] = 100;
	dims[3] = 1;
	accumulate = false;
	directives.Add("axes",new tw::input::Enums<tw::grid::axis>(tw::grid::axis_map(),&ax[1],&ax[2],&ax[3]));
	directives.Add("dimensions",new tw::input::Numbers<tw::Int>(&dims[1],3));
	directives.Add("bounds",new tw::input::Numbers<tw::Float>(&bounds[0],6));
	directives.Add("accumulate",new tw::input::Bool(&accumulate),false);
}

void PhaseSpaceDiagnostic::Start()
{
	Diagnostic::Start();

	if (!headerWritten && task->strip[0].Get_rank()==0)
	{
		std::map<tw::grid::axis,tw::dims> m = {{tw::grid::t,tw::dims::time},{tw::grid::x,tw::dims::length},{tw::grid::y,tw::dims::length},{tw::grid::z,tw::dims::length},
			{tw::grid::mass,tw::dims::energy},{tw::grid::px,tw::dims::momentum},{tw::grid::py,tw::dims::momentum},{tw::grid::pz,tw::dims::momentum},
			{tw::grid::g,tw::dims::none},{tw::grid::gbx,tw::dims::none},{tw::grid::gby,tw::dims::none},{tw::grid::gbz,tw::dims::none}};
		npy_writer writer;
		std::string xname = filename + ".npy";
		writer.write_header(xname,dims);
		meta_writer meta(native);
		meta.start_entry(xname,filename);
		meta.define_axis(xname,0,meta.refine_label(tw::grid::pretty_axis_label(ax[0]),vGalileo,space->geo),m[ax[0]]);
		meta.define_axis(xname,1,meta.refine_label(tw::grid::pretty_axis_label(ax[1]),vGalileo,space->geo),m[ax[1]]);
		meta.define_axis(xname,2,meta.refine_label(tw::grid::pretty_axis_label(ax[2]),vGalileo,space->geo),m[ax[2]]);
		meta.define_axis(xname,3,meta.refine_label(tw::grid::pretty_axis_label(ax[3]),vGalileo,space->geo),m[ax[3]]);
		meta.define_axis(xname,4,"$f({\\bf r},{\\bf p})$ (arb.)",tw::dims::none,true);
		meta.finish_entry();
		// don't set headerWritten until end of Finish() due to grid file
	}

	tw::vec3 phaseSpaceMin(bounds[0],bounds[2],bounds[4]);
	tw::vec3 phaseSpaceSize(bounds[1]-bounds[0],bounds[3]-bounds[2],bounds[5]-bounds[4]);
	plot_layout.Resize(dims[1],dims[2],dims[3],phaseSpaceMin,phaseSpaceSize,1);
	fxp.Initialize(plot_layout,task);
	fxp.SetBoundaryConditions(tw::grid::x,tw::bc::fld::dirichletCell,tw::bc::fld::dirichletCell);
	fxp.SetBoundaryConditions(tw::grid::y,tw::bc::fld::dirichletCell,tw::bc::fld::dirichletCell);
	fxp.SetBoundaryConditions(tw::grid::z,tw::bc::fld::dirichletCell,tw::bc::fld::dirichletCell);
}

void PhaseSpaceDiagnostic::Finish()
{
	tw::Int master = 0;
	task->strip[0].Sum(&fxp(0,0,0),&fxp(0,0,0),fxp.TotalBytes(),master);
	if (task->strip[0].Get_rank()==master)
	{
		npy_writer writer;
		fxp.ApplyBoundaryCondition();
		std::string xname = filename + ".npy";
		std::valarray<float> gData(dims[1]*dims[2]*dims[3]);
		if (accumulate && dims[0]!=0)
			writer.get_frame(xname,(char*)&gData[0],dims,0);
		else
			gData = 0.0f;
		for (tw::Int i=1;i<=dims[1];i++)
			for (tw::Int j=1;j<=dims[2];j++)
				for (tw::Int k=1;k<=dims[3];k++)
					gData[(i-1)*dims[2]*dims[3] + (j-1)*dims[3] + (k-1)] += fxp(i,j,k);
		if (accumulate && dims[0]!=0)
			writer.set_frame(xname,(char*)&gData[0],dims,0);
		else
			writer.add_frame(xname,(char*)&gData[0],dims);
		writer.update_shape(xname,dims);
	}

	// Write out the phase space grid data

	std::ofstream gridFile;
	if (task->strip[0].Get_rank()==0)
	{
		if (!accumulate || !headerWritten)
			StartGridFile(gridFile);
		if (!headerWritten) // assuming static grid no need to keep writing spatial points
		{
			for (tw::Int ax=1;ax<=3;ax++)
			{
				// Unlike box diagnostic the grid is always uniform, no message passing needed
				gridFile << "axis" << ax << " = ";
				const tw::Float ds = (bounds[2*ax-1] - bounds[2*ax-2])/tw::Float(dims[ax]);
				const tw::Float s1 = bounds[2*ax-2] + 0.5*ds;
				const tw::Float sN = bounds[2*ax-1] - 0.5*ds;
				for (tw::Int i=0;i<dims[ax]-1;i++)
					gridFile << s1 + ds*tw::Float(i) << " ";
				gridFile << sN << std::endl;
			}
		}
		gridFile.close();
	}
	headerWritten = true;
}

void PhaseSpaceDiagnostic::Particle(const struct Particle& par,tw::Float m0,tw::Float tp)
{
	weights_3D weights;
	tw::vec4 x(tp,space->PositionFromPrimitive(par.q));
	tw::vec4 v(par.p/m0);
	x.zBoost(gammaBoost,1.0);
	v.zBoost(gammaBoost,1.0);
	const tw::Float dV = (bounds[1]-bounds[0])*(bounds[3]-bounds[2])*(bounds[5]-bounds[4])/(dims[1]*dims[2]*dims[3]);
	tw::vec3 x3 = x.spatial();
	if (theRgn->Inside(x3,*space))
	{
		space->CurvilinearToCartesian(&x3);
		x3 -= vGalileo*x[0];
		x = tw::vec4(x[0],x3);
		std::map<tw::grid::axis,tw::Float> m =
			{{tw::grid::t,x[0]},{tw::grid::x,x[1]},{tw::grid::y,x[2]},{tw::grid::z,x[3]},
			{tw::grid::mass,m0*v[0]},{tw::grid::px,m0*v[1]},{tw::grid::py,m0*v[2]},{tw::grid::pz,m0*v[3]},
			{tw::grid::g,v[0]},{tw::grid::gbx,v[1]},{tw::grid::gby,v[2]},{tw::grid::gbz,v[3]}};
		tw::vec3 q(m[ax[1]],m[ax[2]],m[ax[3]]);
		if (dims[1]==1) q.x = 0.5*(bounds[0]+bounds[1]);
		if (dims[2]==1) q.y = 0.5*(bounds[2]+bounds[3]);
		if (dims[3]==1) q.z = 0.5*(bounds[4]+bounds[5]);
		if (q.x>=bounds[0] && q.x<=bounds[1] && q.y>=bounds[2] && q.y<=bounds[3] && q.z>=bounds[4] && q.z<=bounds[5])
		{
			fxp.GetWeights(&weights,q);
			fxp.InterpolateOnto( par.number/dV, weights );
		}
	}
}
