/*
 *  twpost.h
 *  twpost
 *
 *  Created by Daniel Gordon on 10/1/09.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */

void ExtractExtension(string* extension,const string& pathname);
void ExtractNameWithoutExtension(string* name,const string& pathname);

inline void ReverseBytes(char *bytes,int64 n)
{
	int64 i;
	char temp;
	for (i=0;i<n/2;i++)
	{
		temp = bytes[i];
		bytes[i] = bytes[n-1-i];
		bytes[n-1-i] = temp;
	}
}

inline void ReadBigEndian(ifstream *inFile,char *bytes,int64 len,int64 blockSize)
{
	int64 i;
	inFile->read(bytes,len);
	#ifdef THIS_MACHINE_LITTLE_ENDIAN
		for (i=0;i<len/blockSize;i++)
			ReverseBytes(&bytes[i*blockSize],blockSize);
	#endif
}

inline void WriteBigEndian(ofstream *outFile,char *bytes,int64 len,int64 blockSize)
{
	int64 i;
	#ifdef THIS_MACHINE_LITTLE_ENDIAN
		for (i=0;i<len/blockSize;i++)
			ReverseBytes(&bytes[i*blockSize],blockSize);
	#endif
	outFile->write(bytes,len);
	#ifdef THIS_MACHINE_LITTLE_ENDIAN
		for (i=0;i<len/blockSize;i++)
			ReverseBytes(&bytes[i*blockSize],blockSize);
	#endif
}

struct dvdat
{
	int64 xDim,yDim,zDim,fDim;
	int64 length;
	float x0,x1,y0,y1,z0,z1,t0,t1;
	valarray<float> theData;
	
	dvdat();
	~dvdat();
	float& operator () (int64 i,int64 j,int64 k,int64 s)
	{
		return theData[i + j*xDim + k*xDim*yDim + s*xDim*yDim*zDim];
	}
	dvdat& operator = (dvdat& src);
	void ReadFrame(ifstream* inFile, int64 f);
	void WriteFrame(ofstream* outFile, int64 f);
	void ReadFile(ifstream* inFile);
	void WriteFile(ofstream* outFile);
	void WriteVTKFormat(const string& filename);
	void Allocate(dvdat *in);
	void Allocate(int64 xDim,int64 yDim,int64 zDim,int64 fDim,float x0,float x1,float y0,float y1,float z0,float z1,float t0,float t1);
	void Extract(dvdat *src);
	void Extract(dvdat *src,int64 x1,int64 x2,int64 y1,int64 y2,int64 z1,int64 z2,int64 f1,int64 f2);
	void Downsample(dvdat *src);
	void Smooth();
	void Add(dvdat *in1,dvdat *in2);
	void Multiply(dvdat *in1,dvdat *in2);
	void Intensity(dvdat *in1,dvdat *in2);
	void Radial(dvdat *in1,dvdat *in2);
	void Azimuthal(dvdat *in1,dvdat *in2);
	void Divide(dvdat *in1,dvdat *in2);
	void AngularSpectrum(dvdat *in);
	void Emittance();
	void SliceEmittance(dvdat *out);
	void Overlap();
};

