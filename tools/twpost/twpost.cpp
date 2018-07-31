#include "definitions.h"
#include "fft.h"
#include "twpost.h"

using namespace std;

void Usage()
{
	cout << "quit - exit program (eqv: stop, end, bye)" << endl;
	cout << "help - this (eqv: h , ?)" << endl;
	cout << "------In Memory Operations------" << endl;
	cout << "r1 - read data into slot 1" << endl;
	cout << "r2 - read data into slot 2" << endl;
	cout << "r3 - read data into slot 3" << endl;
	cout << "w1 - write data from slot 1" << endl;
	cout << "w2 - write data from slot 2" << endl;
	cout << "w3 - write data from slot 3" << endl;
	cout << "wv - write data from slot 1, 2, and 3 into as vector components in vtk" << endl;
	cout << "ex - select subset of data" << endl;
	cout << "mv - move data from one slot to another" << endl;
	cout << "add - put sum of slots 1 and 2 into slot 3" << endl;
	cout << "mul - put product of slots 1 and 2 into slot 3" << endl;
	cout << "int - put sum of squares of slots 1 and 2 into slot 3" << endl;
	cout << "rad - put radial value of slots 1 (X) and 2 (Y) into slot3" << endl;
	cout << "azm - put azimuthal value of slots 1 (X) and 2 (Y) into slot3" << endl;
	cout << "div - put (slot 1 / slot 2) into slot 3" << endl;
	cout << "gp - report grid parameters" << endl;
	cout << "lx - linear transform" << endl;
	cout << "as - put angular spectrum from slot 1 into slot 3" << endl;
	cout << "em - calculate emittance" << endl;
	cout << "se - put slice emittance into slot 3" << endl;
	cout << "ds - downsample" << endl;
	cout << "sm - smooth" << endl;
	cout << "ov - overlap integral" << endl;
	cout << "------In Place Frame-by-Frame Operations------" << endl;
	cout << "clip - clip data to reduce file size" << endl;
	cout << "------Notes------" << endl;
	cout << "For IN-MEMORY-OPERATIONS all the data must fit in main memory" << endl;
	cout << "For FRAME-by-FRAME OPERATIONS each frame must fit in main memory" << endl;
}

void ExtractExtension(string* extension,const string& pathname)
{
	int i = pathname.length();
	while (i>0)
	{
		i--;
		if (pathname[i]=='.')
		{
			*extension = pathname.substr(i+1);
			return;
		}
		if (pathname[i]=='/')
		{
			*extension = "";
			return;
		}
	}
	*extension = "";
}

void ExtractNameWithoutExtension(string* name,const string& pathname)
{
	int i = pathname.length();
	while (i>0)
	{
		i--;
		if (pathname[i]=='.')
		{
			*name = pathname.substr(0,i);
			return;
		}
		if (pathname[i]=='/')
		{
			*name = pathname.substr(i+1);
			return;
		}
	}
	*name = pathname;
}

dvdat::dvdat()
{
}

dvdat::~dvdat()
{
}

void dvdat::ReadFile(ifstream* inFile)
{
	int64 i,j,k,s;
	int temp32;
	char header[17];
	int64 headerSize;
	float data;
	
	inFile->read(header,16);
	headerSize = 16*sizeof(char) + 3*sizeof(int) + 6*sizeof(float);
	
	ReadBigEndian(inFile,(char*)&temp32,sizeof(int),sizeof(int)); xDim = temp32;
	ReadBigEndian(inFile,(char*)&temp32,sizeof(int),sizeof(int)); yDim = temp32;
	ReadBigEndian(inFile,(char*)&temp32,sizeof(int),sizeof(int)); zDim = temp32;
	ReadBigEndian(inFile,(char*)&x0,sizeof(float),sizeof(float));
	ReadBigEndian(inFile,(char*)&x1,sizeof(float),sizeof(float));
	ReadBigEndian(inFile,(char*)&y0,sizeof(float),sizeof(float));
	ReadBigEndian(inFile,(char*)&y1,sizeof(float),sizeof(float));
	ReadBigEndian(inFile,(char*)&z0,sizeof(float),sizeof(float));
	ReadBigEndian(inFile,(char*)&z1,sizeof(float),sizeof(float));
	t0 = 0.0;
	t1 = 1.0;
	
	cout << "dims: " << xDim << " x " << yDim << " x " << zDim << endl;
	
	inFile->seekg(0,ios::end);
	length = inFile->tellg();
	fDim = (length - headerSize)/(sizeof(float)*xDim*yDim*zDim);
	
	cout << "frames: " << fDim << endl;
	
	theData.resize(xDim*yDim*zDim*fDim);
	inFile->seekg(headerSize,ios::beg);
	for (s=0;s<fDim;s++)
	{
		cout << "frame " << s << endl;
		for (k=0;k<zDim;k++)
			for (j=0;j<yDim;j++)
				for (i=0;i<xDim;i++)
				{
					ReadBigEndian(inFile,(char*)&data,sizeof(float),sizeof(float));
					(*this)(i,j,k,s) = data;
				}
	}
}

void dvdat::ReadFrame(ifstream* inFile,int64 theFrame)
{
	int64 i,j,k,s;
	int temp32;
	char header[17];
	int64 headerSize;
	float data;
	
	inFile->seekg(0,ios::beg);
	inFile->read(header,16);
	headerSize = 16*sizeof(char) + 3*sizeof(int) + 6*sizeof(float);
	
	ReadBigEndian(inFile,(char*)&temp32,sizeof(int),sizeof(int)); xDim = temp32;
	ReadBigEndian(inFile,(char*)&temp32,sizeof(int),sizeof(int)); yDim = temp32;
	ReadBigEndian(inFile,(char*)&temp32,sizeof(int),sizeof(int)); zDim = temp32;
	ReadBigEndian(inFile,(char*)&x0,sizeof(float),sizeof(float));
	ReadBigEndian(inFile,(char*)&x1,sizeof(float),sizeof(float));
	ReadBigEndian(inFile,(char*)&y0,sizeof(float),sizeof(float));
	ReadBigEndian(inFile,(char*)&y1,sizeof(float),sizeof(float));
	ReadBigEndian(inFile,(char*)&z0,sizeof(float),sizeof(float));
	ReadBigEndian(inFile,(char*)&z1,sizeof(float),sizeof(float));
	t0 = 0.0;
	t1 = 1.0;
	
	cout << "dims: " << xDim << " x " << yDim << " x " << zDim << endl;
	
	inFile->seekg(0,ios::end);
	length = inFile->tellg();
	fDim = (length - headerSize)/(sizeof(float)*xDim*yDim*zDim);
	
	cout << "frames: " << fDim << endl;
	
	theData.resize(xDim*yDim*zDim);
	inFile->seekg(headerSize + sizeof(float)*xDim*yDim*zDim*theFrame,ios::beg);
	cout << "read frame " << theFrame << endl;
	for (k=0;k<zDim;k++)
		for (j=0;j<yDim;j++)
			for (i=0;i<xDim;i++)
			{
				ReadBigEndian(inFile,(char*)&data,sizeof(float),sizeof(float));
				(*this)(i,j,k,0) = data;
			}
}

void dvdat::WriteFile(ofstream* outFile)
{
	int64 i,j,k,s;
	int temp32;
	
	string version;
	float data;
	
	version = "DataViewer 2.0.0";
	outFile->write(version.data(),16);
	temp32 = xDim; WriteBigEndian(outFile,(char*)&temp32,sizeof(int),sizeof(int));
	temp32 = yDim; WriteBigEndian(outFile,(char*)&temp32,sizeof(int),sizeof(int));
	temp32 = zDim; WriteBigEndian(outFile,(char*)&temp32,sizeof(int),sizeof(int));
	WriteBigEndian(outFile,(char*)&x0,sizeof(float),sizeof(float));
	WriteBigEndian(outFile,(char*)&x1,sizeof(float),sizeof(float));
	WriteBigEndian(outFile,(char*)&y0,sizeof(float),sizeof(float));
	WriteBigEndian(outFile,(char*)&y1,sizeof(float),sizeof(float));
	WriteBigEndian(outFile,(char*)&z0,sizeof(float),sizeof(float));
	WriteBigEndian(outFile,(char*)&z1,sizeof(float),sizeof(float));
	
	for (s=0;s<fDim;s++)
	{
		if (fDim>1)
			cout << "writing frame " << s << endl;
		for (k=0;k<zDim;k++)
			for (j=0;j<yDim;j++)
				for (i=0;i<xDim;i++)
				{
					data = (*this)(i,j,k,s);
					WriteBigEndian(outFile,(char*)&data,sizeof(float),sizeof(float));
				}
	}
}

void dvdat::WriteFrame(ofstream* outFile,int64 theFrame)
{
	int64 i,j,k;
	
	float data;
	
	for (k=0;k<zDim;k++)
		for (j=0;j<yDim;j++)
			for (i=0;i<xDim;i++)
			{
				data = (*this)(i,j,k,theFrame);
				WriteBigEndian(outFile,(char*)&data,sizeof(float),sizeof(float));
			}
}

void dvdat::WriteVTKFormat(const string& filename)
{
	int64 i,j,k,s;
	std::string bareName;
	std::stringstream ss;
	
	float dx,dy,dz,data;
	dx = (x1 - x0)/float(xDim);
	dy = (y1 - y0)/float(yDim);
	dz = (z1 - z0)/float(zDim);
	
	ofstream *outFile;
	for (s=0;s<fDim;s++)
	{
		cout << "frame " << s << endl;
		ExtractNameWithoutExtension(&bareName,filename);
		ss.str("");
		ss << bareName << s << ".vtk";
		outFile = new ofstream(ss.str().c_str());
		*outFile << "# vtk DataFile Version 2.0" << endl;
		*outFile << "Converted DVDAT File" << endl;
		*outFile << "BINARY" << endl;
		*outFile << "DATASET STRUCTURED_POINTS" << endl;
		*outFile << "DIMENSIONS " << xDim << " " << yDim << " " << zDim << endl;
		*outFile << "ORIGIN " << x0 << " " << y0 << " " << z0 << endl;
		*outFile << "SPACING " << dx << " " << dy << " " << dz << endl;
		*outFile << "POINT_DATA " << xDim*yDim*zDim << endl;
		*outFile << "SCALARS volume_scalars float 1" << endl;
		*outFile << "LOOKUP_TABLE default" << endl;
		for (k=0;k<zDim;k++)
			for (j=0;j<yDim;j++)
				for (i=0;i<xDim;i++)
				{
					data = (*this)(i,j,k,s);
					WriteBigEndian(outFile,(char*)&data,sizeof(float),sizeof(float));
				}
		outFile->close();
		delete outFile;
	}
}

void dvdat::Allocate(dvdat *in)
{
	xDim = in->xDim;
	yDim = in->yDim;
	zDim = in->zDim;
	fDim = in->fDim;
	length = in->length;
	x0 = in->x0;
	y0 = in->y0;
	z0 = in->z0;
	x1 = in->x1;
	y1 = in->y1;
	z1 = in->z1;
	t0 = in->t0;
	t1 = in->t1;
	theData.resize(xDim*yDim*zDim*fDim);
}

void dvdat::Allocate(int64 xDim,int64 yDim,int64 zDim,int64 fDim,float x0,float x1,float y0,float y1,float z0,float z1,float t0,float t1)
{
	int64 headerSize = 16*sizeof(char) + 3*sizeof(int) + 6*sizeof(float);
	this->xDim = xDim;
	this->yDim = yDim;
	this->zDim = zDim;
	this->fDim = fDim;
	length = headerSize + sizeof(float)*xDim*yDim*zDim*fDim;
	this->x0 = x0;
	this->y0 = y0;
	this->z0 = z0;
	this->x1 = x1;
	this->y1 = y1;
	this->z1 = z1;
	this->t0 = t0;
	this->t1 = t1;
	theData.resize(xDim*yDim*zDim*fDim);
}

dvdat& dvdat::operator = (dvdat& in)
{
	xDim = in.xDim;
	yDim = in.yDim;
	zDim = in.zDim;
	fDim = in.fDim;
	length = in.length;
	x0 = in.x0;
	y0 = in.y0;
	z0 = in.z0;
	t0 = in.t0;
	x1 = in.x1;
	y1 = in.y1;
	z1 = in.z1;
	t1 = in.t1;
	theData.resize(xDim*yDim*zDim*fDim);
	theData = in.theData;
	return *this;
}

void dvdat::Extract(dvdat *src)
{
	int64 i,j,k,s;
	int64 i1,i2,j1,j2,k1,k2,s1,s2;
	float dt = (src->t1 - src->t0)/float(src->fDim);
	float dx = (src->x1 - src->x0)/float(src->xDim);
	float dy = (src->y1 - src->y0)/float(src->yDim);
	float dz = (src->z1 - src->z0)/float(src->zDim);
	cout << "Enter 8 integers indexed from 0 (x1,x2,y1,y2,z1,z2,f1,f2): " << endl;
	cin >> i1 >> i2 >> j1 >> j2 >> k1 >> k2 >> s1 >> s2;
	t0 = src->t0 + dt*s1;
	t1 = src->t0 + dt*s2;
	x0 = src->x0 + dx*i1;
	x1 = src->x0 + dx*i2;
	y0 = src->y0 + dy*j1;
	y1 = src->y0 + dy*j2;
	z0 = src->z0 + dz*k1;
	z1 = src->z0 + dz*k2;
	xDim = i2 - i1 + 1;
	yDim = j2 - j1 + 1;
	zDim = k2 - k1 + 1;
	fDim = s2 - s1 + 1;
	Allocate(xDim,yDim,zDim,fDim,x0,x1,y0,y1,z0,z1,t0,t1);
	for (s=0;s<fDim;s++)
		for (k=0;k<zDim;k++)
			for (j=0;j<yDim;j++)
				for (i=0;i<xDim;i++)
					(*this)(i,j,k,s) = (*src)(i1+i,j1+j,k1+k,s1+s);
}

void dvdat::Extract(dvdat *src,int64 i1,int64 i2,int64 j1,int64 j2,int64 k1,int64 k2,int64 s1,int64 s2)
{
	int64 i,j,k,s;
	float dt = (src->t1 - src->t0)/float(src->fDim);
	float dx = (src->x1 - src->x0)/float(src->xDim);
	float dy = (src->y1 - src->y0)/float(src->yDim);
	float dz = (src->z1 - src->z0)/float(src->zDim);
	t0 = src->t0 + dt*s1;
	t1 = src->t0 + dt*s2;
	x0 = src->x0 + dx*i1;
	x1 = src->x0 + dx*i2;
	y0 = src->y0 + dy*j1;
	y1 = src->y0 + dy*j2;
	z0 = src->z0 + dz*k1;
	z1 = src->z0 + dz*k2;
	xDim = i2 - i1 + 1;
	yDim = j2 - j1 + 1;
	zDim = k2 - k1 + 1;
	fDim = s2 - s1 + 1;
	Allocate(xDim,yDim,zDim,fDim,x0,x1,y0,y1,z0,z1,t0,t1);
	for (s=0;s<fDim;s++)
		for (k=0;k<zDim;k++)
			for (j=0;j<yDim;j++)
				for (i=0;i<xDim;i++)
					(*this)(i,j,k,s) = (*src)(i1+i,j1+j,k1+k,s1+s);
}

void dvdat::Downsample(dvdat *src)
{
	string ans;
	float summation;
	int64 i,j,k,s,ii,jj,kk,ss;
	int64 sx,sy,sz,sf;
	cout << "Enter 4 skip parameters (x,y,z,f): " << endl;
	cin >> sx >> sy >> sz >> sf;
	cout << "Average? (y/n):" << endl;
	cin >> ans;
	
	t0 = src->t0;
	t1 = src->t1;
	x0 = src->x0;
	x1 = src->x1;
	y0 = src->y0;
	y1 = src->y1;
	z0 = src->z0;
	z1 = src->z1;
	xDim = src->xDim/sx;
	yDim = src->yDim/sy;
	zDim = src->zDim/sz;
	fDim = src->fDim/sf;
	Allocate(xDim,yDim,zDim,fDim,x0,x1,y0,y1,z0,z1,t0,t1);
	for (s=0;s<src->fDim;s+=sf)
		for (k=0;k<src->zDim;k+=sz)
			for (j=0;j<src->yDim;j+=sy)
				for (i=0;i<src->xDim;i+=sx)
					if (ans=="y")
					{
						summation = 0.0;
						for (ss=0;ss<sf;ss++)
							for (kk=0;kk<sz;kk++)
								for (jj=0;jj<sy;jj++)
									for (ii=0;ii<sx;ii++)
										summation += (*src)(i+ii,j+jj,k+kk,s+ss);
						summation /= sf*sz*sy*sx;
						(*this)(i/sx,j/sy,k/sz,s/sf) = summation;
					}
					else
					{
						(*this)(i/sx,j/sy,k/sz,s/sf) = (*src)(i,j,k,s);
					}
}

void dvdat::Smooth()
{
	float last,temp;
	int64 i,j,k,s,p,passes;
	cout << "Passes: " << endl;
	cin >> passes;
	
	for (p=0;p<passes;p++)
	{
		for (s=0;s<fDim;s++)
		{
			if (xDim>2)
			{
				for (k=0;k<zDim;k++)
					for (j=0;j<yDim;j++)
					{
						last = (*this)(0,j,k,s);
						(*this)(0,j,k,s) = 0.75*(*this)(0,j,k,s) + 0.25*(*this)(1,j,k,s);
						for (i=1;i<xDim-1;i++)
						{
							temp = (*this)(i,j,k,s);
							(*this)(i,j,k,s) = 0.25*last + 0.5*(*this)(i,j,k,s) + 0.25*(*this)(i+1,j,k,s);
							last = temp;
						}
						(*this)(xDim-1,j,k,s) = 0.25*(*this)(xDim-2,j,k,s) + 0.75*(*this)(xDim-1,j,k,s);
					}
			}
			if (yDim>2)
			{
				for (k=0;k<zDim;k++)
					for (i=0;i<xDim;i++)
					{
						last = (*this)(i,0,k,s);
						(*this)(i,0,k,s) = 0.75*(*this)(i,0,k,s) + 0.25*(*this)(i,1,k,s);
						for (j=1;j<yDim-1;j++)
						{
							temp = (*this)(i,j,k,s);
							(*this)(i,j,k,s) = 0.25*last + 0.5*(*this)(i,j,k,s) + 0.25*(*this)(i,j+1,k,s);
							last = temp;
						}
						(*this)(i,yDim-1,k,s) = 0.25*(*this)(i,yDim-2,k,s) + 0.75*(*this)(i,yDim-1,k,s);
					}
			}
			if (zDim>2)
			{
				for (j=0;j<yDim;j++)
					for (i=0;i<xDim;i++)
					{
						last = (*this)(i,j,0,s);
						(*this)(i,j,0,s) = 0.75*(*this)(i,j,0,s) + 0.25*(*this)(i,j,1,s);
						for (k=1;k<zDim-1;k++)
						{
							temp = (*this)(i,j,k,s);
							(*this)(i,j,k,s) = 0.25*last + 0.5*(*this)(i,j,k,s) + 0.25*(*this)(i,j,k+1,s);
							last = temp;
						}
						(*this)(i,j,zDim-1,s) = 0.25*(*this)(i,j,zDim-2,s) + 0.75*(*this)(i,j,zDim-1,s);
					}
			}
		}
	}
}

void dvdat::Overlap()
{
	float R,r,z,dr,dz,ans,integrand;
	int64 i,j,k,s;
	dr = (x1-x0)/xDim;
	dz = (z1-z0)/zDim;
	
	for (s=0;s<fDim;s++)
	{
		ans = 0.0;
		for (k=0;k<zDim;k++)
			for (j=0;j<yDim;j++)
				for (i=0;i<xDim;i++)
				{
					r = x0 + dr*(i+0.5);
					z = z0 + dz*(k+0.5);
					R = sqrt(r*r+z*z);
					integrand = fabs(  (1.0/sqrt(32*pi))*R*exp(-0.5*R)*z/R  ); // p-state
					integrand *= sqrt((*this)(i,j,k,s)); // bohmian amplitude
					integrand *= dz*pi*(pow(r+0.5*dr,2.0) - pow(r-0.5*dr,2.0)); // volume element
					ans += integrand;
				}
		cout << ans*ans << endl;
	}
}

void dvdat::Add(dvdat *in1,dvdat *in2)
{
	int64 i,j,k,s;
	
	for (s=0;s<fDim;s++)
	{
		cout << "frame " << s << endl;
		for (k=0;k<zDim;k++)
			for (j=0;j<yDim;j++)
				for (i=0;i<xDim;i++)
					(*this)(i,j,k,s) = (*in1)(i,j,k,s) + (*in2)(i,j,k,s);
	}
}

void dvdat::Multiply(dvdat *in1,dvdat *in2)
{
	int64 i,j,k,s;
	
	for (s=0;s<fDim;s++)
	{
		cout << "frame " << s << endl;
		for (k=0;k<zDim;k++)
			for (j=0;j<yDim;j++)
				for (i=0;i<xDim;i++)
					(*this)(i,j,k,s) = (*in1)(i,j,k,s) * (*in2)(i,j,k,s);
	}
}

void dvdat::Intensity(dvdat *in1,dvdat *in2)
{
	int64 i,j,k,s;
	
	cout << "Intensity: I = X*X + Y*Y" << endl;
	
	for (s=0;s<fDim;s++)
	{
		cout << "frame " << s << endl;
		for (k=0;k<zDim;k++)
			for (j=0;j<yDim;j++)
				for (i=0;i<xDim;i++)
					(*this)(i,j,k,s) = (*in1)(i,j,k,s) * (*in1)(i,j,k,s) + (*in2)(i,j,k,s) * (*in2)(i,j,k,s);
	}
}

void dvdat::Azimuthal(dvdat *in1,dvdat *in2)
{
	int64 i,j,k,s;
	float dx,dy,x,y,phi;
	
	dx = (x1 - x0) / float(xDim);
	dy = (y1 - y0) / float(yDim);
	
	cout << "Azimuthal: Phi = -X sin(phi) + Y cos(phi)" << endl;
	
	for (s=0;s<fDim;s++)
	{
		cout << "frame " << s << endl;
		for (k=0;k<zDim;k++)
			for (j=0;j<yDim;j++)
				for (i=0;i<xDim;i++)
				{
					x = x0 + float(i)*dx + 0.5*dx;
					y = y0 + float(j)*dy + 0.5*dy;
					phi = atan2(y, x);
					(*this)(i,j,k,s) = -(*in1)(i,j,k,s)*sin(phi) + (*in2)(i,j,k,s)*cos(phi);
				}
	}
}

void dvdat::Radial(dvdat *in1,dvdat *in2)
{
	int64 i,j,k,s;
	float dx,dy,x,y,phi;
	
	dx = (x1 - x0) / float(xDim);
	dy = (y1 - y0) / float(yDim);
	
	cout << "Radial: R = X cos(phi) + Y sin(phi )" << endl;
	
	for (s=0;s<fDim;s++)
	{
		cout << "frame " << s << endl;
		for (k=0;k<zDim;k++)
			for (j=0;j<yDim;j++)
				for (i=0;i<xDim;i++)
				{
					x = x0 + float(i)*dx + 0.5*dx;
					y = y0 + float(j)*dy + 0.5*dy;
					phi = atan2(y, x);
					(*this)(i,j,k,s) = (*in1)(i,j,k,s)*cos(phi) + (*in2)(i,j,k,s)*sin(phi);
				}
	}
}

void dvdat::Divide(dvdat *in1,dvdat *in2)
{
	int64 i,j,k,s;
	
	for (s=0;s<fDim;s++)
	{
		cout << "frame " << s << endl;
		for (k=0;k<zDim;k++)
			for (j=0;j<yDim;j++)
				for (i=0;i<xDim;i++)
					(*this)(i,j,k,s) = (*in1)(i,j,k,s) / (*in2)(i,j,k,s);
	}
}

void dvdat::AngularSpectrum(dvdat *in)
{
	int64 i,j,k,s;
	float k2,kx,ky,kz;
	for (s=0;s<fDim;s++)
	{
		cout << "frame " << s << endl;
		for (k=0;k<zDim;k++)
			for (j=0;j<yDim;j++)
				for (i=0;i<xDim;i++)
				{
					kx = x0;
					ky = y0;
					kz = z0;
					k2 = kx*kx + ky*ky + kz*kz;
				}
	}
}

void dvdat::SliceEmittance(dvdat *out)
{
	int64 i,j,k,s;
	float e2,val,x,y,dx,dy,xc,yc;
	float x2Avg,y2Avg,xyAvg,numPars;
	
	dx = (x1 - x0) / float(xDim);
	dy = (y1 - y0) / float(yDim);
	
	for (s=0;s<fDim;s++)
	{
		numPars = xc = yc = 0.0;
		for (k=0;k<zDim;k++)
			for (j=0;j<yDim;j++)
				for (i=0;i<xDim;i++)
				{
					x = x0 + float(i)*dx + 0.5*dx;
					y = y0 + float(j)*dy + 0.5*dy;
					val = (*this)(i,j,k,s);
					xc += x*val;
					yc += y*val;
					numPars += val;
				}
		xc /= numPars;
		yc /= numPars;

		numPars = x2Avg = y2Avg = xyAvg = 0.0;
		for (k=0;k<zDim;k++)
			for (j=0;j<yDim;j++)
				for (i=0;i<xDim;i++)
				{
					x = x0 + float(i)*dx + 0.5*dx - xc;
					y = y0 + float(j)*dy + 0.5*dy - yc;
					val = (*this)(i,j,k,s);
					x2Avg += x*x*val;
					y2Avg += y*y*val;
					xyAvg += x*y*val;
					numPars += val;
				}
		x2Avg /= numPars;
		y2Avg /= numPars;
		xyAvg /= numPars;
		e2 = x2Avg*y2Avg - xyAvg*xyAvg;
		(*out)(s,0,0,0) = numPars==0 ? 0.0 : sqrt(e2);
	}
	
	cout << "slice emittance has been written to output slot (x axis is time)" << endl;
}

void dvdat::Emittance()
{
	int64 i,j,k,s;
	float e2,val,x,y,dx,dy,xc,yc;
	float x2Avg,y2Avg,xyAvg,numPars;
	
	dx = (x1 - x0) / float(xDim);
	dy = (y1 - y0) / float(yDim);
	
	cout << "Transverse: e^2 = <x^2><y^2> - <xy>^2" << endl;
	
	numPars = xc = yc = 0.0;
	for (s=0;s<fDim;s++)
		for (k=0;k<zDim;k++)
			for (j=0;j<yDim;j++)
				for (i=0;i<xDim;i++)
				{
					x = x0 + float(i)*dx + 0.5*dx;
					y = y0 + float(j)*dy + 0.5*dy;
					val = (*this)(i,j,k,s);
					xc += x*val;
					yc += y*val;
					numPars += val;
				}
	xc /= numPars;
	yc /= numPars;
	cout << "Centroid at ( " << xc << " , " << yc << " )" << endl;

	numPars = x2Avg = y2Avg = xyAvg = 0.0;
	for (s=0;s<fDim;s++)
		for (k=0;k<zDim;k++)
			for (j=0;j<yDim;j++)
				for (i=0;i<xDim;i++)
				{
					x = x0 + float(i)*dx + 0.5*dx - xc;
					y = y0 + float(j)*dy + 0.5*dy - yc;
					val = (*this)(i,j,k,s);
					x2Avg += x*x*val;
					y2Avg += y*y*val;
					xyAvg += x*y*val;
					numPars += val;
				}
	x2Avg /= numPars;
	y2Avg /= numPars;
	xyAvg /= numPars;
	e2 = x2Avg*y2Avg - xyAvg*xyAvg;
	cout << "sqrt(<x^2>) = " << sqrt(x2Avg) << endl;
	cout << "sqrt(<y^2>) = " << sqrt(y2Avg) << endl;
	cout << "e = " << sqrt(e2) << endl;
	
	dx = (y1 - y0) / float(yDim);
	dy = (t1 - t0) / float(fDim);
	
	cout << "Longitudinal: e^2 = <y^2><t^2> - <yt>^2" << endl;
	
	numPars = xc = yc = 0.0;
	for (s=0;s<fDim;s++)
		for (k=0;k<zDim;k++)
			for (j=0;j<yDim;j++)
				for (i=0;i<xDim;i++)
				{
					x = y0 + float(j)*dx + 0.5*dx;
					y = t0 + float(s)*dy + 0.5*dy;
					val = (*this)(i,j,k,s);
					xc += x*val;
					yc += y*val;
					numPars += val;
				}
	xc /= numPars;
	yc /= numPars;
	cout << "Centroid at ( " << xc << " , " << yc << " )" << endl;
	
	numPars = x2Avg = y2Avg = xyAvg = 0.0;
	for (s=0;s<fDim;s++)
		for (k=0;k<zDim;k++)
			for (j=0;j<yDim;j++)
				for (i=0;i<xDim;i++)
				{
					x = y0 + float(j)*dx + 0.5*dx - xc;
					y = t0 + float(s)*dy + 0.5*dy - yc;
					val = (*this)(i,j,k,s);
					x2Avg += x*x*val;
					y2Avg += y*y*val;
					xyAvg += x*y*val;
					numPars += val;
				}
	x2Avg /= numPars;
	y2Avg /= numPars;
	xyAvg /= numPars;
	e2 = x2Avg*y2Avg - xyAvg*xyAvg;
	cout << "sqrt(<y^2>) = " << sqrt(x2Avg) << endl;
	cout << "sqrt(<t^2>) = " << sqrt(y2Avg) << endl;
	cout << "e = " << sqrt(e2) << endl;
}

int main(int argc,char* argv[])
{
	dvdat *in1,*in2,*out;
	ifstream *inFile;
	ofstream *outFile;
	
	in1 = new dvdat;
	in2 = new dvdat;
	out = new dvdat;
	
	string cmd,filename,extension;
	bool quit = false;
	
	cout << "TurboWAVE Post Processing Tool" << endl;
	cout << "Version 2011-10-11" << endl;
	do
	{
		cout << "cmd>";
		getline(cin,cmd);
		if (cmd=="quit" || cmd=="stop" || cmd=="end" || cmd=="bye" || cmd=="exit")
			quit = true;
		if (cmd=="help" || cmd=="h" || cmd=="?")
			Usage();
		if (cmd[0]=='r')
		{
			dvdat *src = NULL;
			if (cmd[1]=='1')
				src = in1;
			if (cmd[1]=='2')
				src = in2;
			if (cmd[1]=='3')
				src = out;
			if (src!=NULL)
			{
				cout << "filename (include extension): ";
				getline(cin,filename);
				inFile = new ifstream;
				inFile->open(filename.c_str(),ios::binary);
				if (inFile->rdstate() & ios::failbit)
				{
					cout << "ERROR: couldn't read file" << endl;
				}
				else
				{
					src->ReadFile(inFile);
					inFile->close();
				}
				delete inFile;
			}
			else
				cout << "Syntax Error" << endl;
		}
		if (cmd[0]=='w')
		{
			dvdat *src = NULL;
			if (cmd[1]=='1')
				src = in1;
			if (cmd[1]=='2')
				src = in2;
			if (cmd[1]=='3')
				src = out;
			//Write files as vtk Vector Format
			if (cmd[1]=='v')
			{
				int64 i,j,k,s,xDim,yDim,zDim,fDim;
				std::string bareName;
				std::stringstream ss;
	
				float dx,dy,dz,vx,vy,vz,x0,x1,y0,y1,z0,z1;
				x0 = in1->dvdat::x0;
				x1 = in1->dvdat::x1;
				xDim = in1->dvdat::xDim;
				y0 = in1->dvdat::y0;
				y1 = in1->dvdat::y1;
				yDim = in1->dvdat::yDim;
				z0 = in1->dvdat::z0;
				z1 = in1->dvdat::z1;
				zDim = in1->dvdat::zDim;
				fDim = in1->dvdat::fDim;


				dx = (x1 - x0)/float(xDim);
				dy = (y1 - y0)/float(yDim);
				dz = (z1 - z0)/float(zDim);
	
				cout << "filename (must be vtk file format): ";
				getline(cin,filename);
				ExtractExtension(&extension,filename);
				ExtractNameWithoutExtension(&bareName,filename);
				for (s=0;s<fDim;s++)
				{
					ofstream *outFile;
					cout << "frame " << s << endl;
					ss.str("");
					ss << bareName << s << ".vtk";
					outFile = new ofstream(ss.str().c_str());
					*outFile << "# vtk DataFile Version 2.0" << endl;
					*outFile << "Converted DVDAT File" << endl;
					*outFile << "ASCII" << endl;
					*outFile << "DATASET STRUCTURED_POINTS" << endl;
					*outFile << "DIMENSIONS " << xDim << " " << yDim << " " << zDim << endl;
					*outFile << "ORIGIN " << x0 << " " << y0 << " " << z0 << endl;
					*outFile << "SPACING " << dx << " " << dy << " " << dz << endl;
					*outFile << "POINT_DATA " << xDim*yDim*zDim << endl;
					*outFile << "VECTORS Field_scalars float" << std::endl;
					for (k=0;k<zDim;k++)
						for (j=0;j<yDim;j++)
							for (i=0;i<xDim;i++)
							{
								vx = in1->dvdat::operator()(i,j,k,s);
								vy = in2->dvdat::operator()(i,j,k,s);
								vz = out->dvdat::operator()(i,j,k,s);
								*outFile << vx << "\t" << vy << "\t" << vz << std::endl;
							}
					outFile->close();
					delete outFile;
				}

			}
			if (src!=NULL&&cmd[1]!='v')
			{
				cout << "filename (extension determines format, dvdat, or vtk): ";
				getline(cin,filename);
				ExtractExtension(&extension,filename);
				if (extension=="dvdat")
				{
					outFile = new ofstream;
					outFile->open(filename.c_str(),ios::binary);
					src->WriteFile(outFile);
					outFile->close();
					delete outFile;
				}
				if (extension=="vtk")
					src->WriteVTKFormat(filename);
				if (extension!="dvdat" && extension!="vtk")
					cout << "Error: Unknown File Extension" << endl;
			}
			if (src==NULL&&cmd[1]!='v')
				cout << "Syntax Error" << endl;
		}
		if (cmd=="mv")
		{
			string str;
			dvdat *src = NULL;
			dvdat *dst = NULL;
			cout << "source slot (1,2,3): ";
			getline(cin,str);
			if (str=="1")
				src = in1;
			if (str=="2")
				src = in2;
			if (str=="3")
				src = out;
			cout << "destination slot (1,2,3): ";
			getline(cin,str);
			if (str=="1")
				dst = in1;
			if (str=="2")
				dst = in2;
			if (str=="3")
				dst = out;
			if (src!=NULL && dst!=NULL)
				*dst = *src;
		}
		if (cmd=="ex")
		{
			string str;
			dvdat *src = NULL;
			dvdat *dst = NULL;
			cout << "source slot (1,2,3): ";
			getline(cin,str);
			if (str=="1")
				src = in1;
			if (str=="2")
				src = in2;
			if (str=="3")
				src = out;
			cout << "destination slot (1,2,3): ";
			getline(cin,str);
			if (str=="1")
				dst = in1;
			if (str=="2")
				dst = in2;
			if (str=="3")
				dst = out;
			if (src!=NULL && dst!=NULL)
				dst->Extract(src);
		}
		if (cmd=="ds")
		{
			string str;
			dvdat *src = NULL;
			dvdat *dst = NULL;
			cout << "source slot (1,2,3): ";
			getline(cin,str);
			if (str=="1")
				src = in1;
			if (str=="2")
				src = in2;
			if (str=="3")
				src = out;
			cout << "destination slot (1,2,3): ";
			getline(cin,str);
			if (str=="1")
				dst = in1;
			if (str=="2")
				dst = in2;
			if (str=="3")
				dst = out;
			if (src!=NULL && dst!=NULL)
				dst->Downsample(src);
		}
		if (cmd=="sm")
		{
			string str;
			dvdat *in_out = NULL;
			cout << "smooth in-place slot (1,2,3): ";
			getline(cin,str);
			if (str=="1")
				in_out = in1;
			if (str=="2")
				in_out = in2;
			if (str=="3")
				in_out = out;
			if (in_out!=NULL)
				in_out->Smooth();
		}
		if (cmd=="ov")
		{
			string str;
			dvdat *in_out = NULL;
			cout << "overlap integral for slot (1,2,3): ";
			getline(cin,str);
			if (str=="1")
				in_out = in1;
			if (str=="2")
				in_out = in2;
			if (str=="3")
				in_out = out;
			if (in_out!=NULL)
				in_out->Overlap();
		}			
		if (cmd=="add")
		{
			out->Allocate(in1);
			out->Add(in1,in2);
		}
		if (cmd=="mul")
		{
			out->Allocate(in1);
			out->Multiply(in1,in2);
		}
		if (cmd=="int")
		{
			out->Allocate(in1);
			out->Intensity(in1,in2);
		}
		if (cmd=="azm")
		{
			out->Allocate(in1);
			out->Azimuthal(in1,in2);
		}
		if (cmd=="rad")
		{
			out->Allocate(in1);
			out->Radial(in1,in2);
		}
		if (cmd=="div")
		{
			out->Allocate(in1);
			out->Divide(in1,in2);
		}
		if (cmd=="as")
		{
			out->Allocate(in1);
			out->AngularSpectrum(in1);
		}
		
		if (cmd=="gp")
		{
			cout << "Grid Parameters in Slot 1" << endl;
			cout << "-------------------------" << endl;
			cout << in1->xDim << " x " << in1->yDim << " x " << in1->zDim << endl;
			cout << in1->fDim << " frames" << endl;
			cout << "x range = [ " << in1->x0 << " , " << in1->x1 << " ]" << endl;
			cout << "y range = [ " << in1->y0 << " , " << in1->y1 << " ]" << endl;
			cout << "z range = [ " << in1->z0 << " , " << in1->z1 << " ]" << endl;
			cout << "t range = [ " << in1->t0 << " , " << in1->t1 << " ]" << endl;
			cout << "Grid Parameters in Slot 2" << endl;
			cout << "-------------------------" << endl;
			cout << in2->xDim << " x " << in2->yDim << " x " << in2->zDim << endl;
			cout << in2->fDim << " frames" << endl;
			cout << "x range = [ " << in2->x0 << " , " << in2->x1 << " ]" << endl;
			cout << "y range = [ " << in2->y0 << " , " << in2->y1 << " ]" << endl;
			cout << "z range = [ " << in2->z0 << " , " << in2->z1 << " ]" << endl;
			cout << "t range = [ " << in2->t0 << " , " << in2->t1 << " ]" << endl;
			cout << "Grid Parameters in Slot 3" << endl;
			cout << "-------------------------" << endl;
			cout << out->xDim << " x " << out->yDim << " x " << out->zDim << endl;
			cout << out->fDim << " frames" << endl;
			cout << "x range = [ " << out->x0 << " , " << out->x1 << " ]" << endl;
			cout << "y range = [ " << out->y0 << " , " << out->y1 << " ]" << endl;
			cout << "z range = [ " << out->z0 << " , " << out->z1 << " ]" << endl;
			cout << "t range = [ " << out->t0 << " , " << out->t1 << " ]" << endl;
		}
		
		if (cmd=="lx")
		{
			float A,B;
			string str,theAxis;
			dvdat *theData = NULL;
			cout << "slot (1,2,3): ";
			getline(cin,str);
			if (str=="1")
				theData = in1;
			if (str=="2")
				theData = in2;
			if (str=="3")
				theData = out;
			cout << "axis (x,y,z,t): ";
			getline(cin,theAxis);
			cout << "transformation will be x' = Ax + B" << endl;
			cout << "A: ";
			cin >> A;
			cout << "B: ";
			cin >> B;
			if (theAxis=="x")
			{
				theData->x0 = A*theData->x0 + B;
				theData->x1 = A*theData->x1 + B;
			}
			if (theAxis=="y")
			{
				theData->y0 = A*theData->y0 + B;
				theData->y1 = A*theData->y1 + B;
			}
			if (theAxis=="z")
			{
				theData->z0 = A*theData->z0 + B;
				theData->z1 = A*theData->z1 + B;
			}
			if (theAxis=="t")
			{
				theData->t0 = A*theData->t0 + B;
				theData->t1 = A*theData->t1 + B;
			}
		}
		
		if (cmd=="em")
		{
			string str;
			dvdat *theData = NULL;
			cout << "slot (1,2,3): ";
			getline(cin,str);
			if (str=="1")
				theData = in1;
			if (str=="2")
				theData = in2;
			if (str=="3")
				theData = out;
			theData->Emittance();
		}
		if (cmd=="se")
		{
			string str;
			dvdat *theData = NULL;
			cout << "slot (1,2): ";
			getline(cin,str);
			if (str=="1")
				theData = in1;
			if (str=="2")
				theData = in2;
			out->Allocate(theData->fDim,1,1,1,theData->t0,theData->t1,0.0,1.0,0.0,1.0,0.0,1.0);
			theData->SliceEmittance(out);
		}
		if (cmd=="clip")
		{
			int64 i,i1,i2,j1,j2,k1,k2,s1,s2;
			dvdat *fullData = new dvdat;
			dvdat *clipData = new dvdat;
			cout << "filename (include extension): ";
			getline(cin,filename);
			inFile = new ifstream;
			inFile->open(filename.c_str(),ios::binary);
			if (inFile->rdstate() & ios::failbit)
			{
				cout << "ERROR: couldn't read file" << endl;
			}
			else
			{
				outFile = new ofstream;
				filename += "clip.dvdat";
				outFile->open(filename.c_str(),ios::binary);
				cout << "Enter 8 integers indexed from 0 (x1,x2,y1,y2,z1,z2,f1,f2): " << endl;
				cin >> i1 >> i2 >> j1 >> j2 >> k1 >> k2 >> s1 >> s2;
				for (i=s1;i<=s2;i++)
				{
					cout << "Processing Frame " << i << endl;
					fullData->ReadFrame(inFile,i);
					clipData->Extract(fullData,i1,i2,j1,j2,k1,k2,0,0);
					if (i==s1)
						clipData->WriteFile(outFile);
					else
						clipData->WriteFrame(outFile,0);
				}
				outFile->close();
				delete outFile;
				inFile->close();
			}
			delete inFile;
		}

	}
	while (!quit);
	return 0;
}