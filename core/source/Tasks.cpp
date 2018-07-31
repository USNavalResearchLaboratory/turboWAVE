#include "definitions.h"
#include "tasks.h"
#include <algorithm>

////////////////////
//      TASK      //
////////////////////


Task::Task()
{
	tw::Int i;
	for (i=0;i<4;i++)
	{
		cornerCell[i] = 1;
		localCells[i] = 1;
		globalCells[i] = 1;
		domainIndex[i] = 0;
		domains[i] = 1;
		periodic[i] = 0;
	}
	deviceSearchString = "tesla";
	platformSearchString = "cuda";
}

Task::~Task()
{
	// WARNING: Assumes task was initialized
	#ifdef USE_OPENCL

	clReleaseCommandQueue(commandQueue);

	clReleaseProgram(fieldProgram);

	clReleaseKernel(k_fillVec4Field);
	clReleaseKernel(k_fieldToBoundary);
	clReleaseKernel(k_boundaryToField);
	clReleaseKernel(k_zeroGhostCells);
	clReleaseKernel(k_swapBuffers);
	clReleaseKernel(k_destructiveSum);
	clReleaseKernel(k_destructiveNorm1);
	clReleaseKernel(k_complexMod2);
	clReleaseKernel(k_destructiveComplexMod2);
	clReleaseKernel(k_MADD);
	clReleaseKernel(k_weightByVolume);

	clReleaseContext(context);
	#endif
}

void Task::Initialize(tw::Int *doms,tw::Int *gcells,tw::Int *cyclic)
{
	tw::Int i;
	tw::Int aperiodic[4] = { 0,0,0,0 };

	for (i=0;i<4;i++)
	{
		domains[i] = doms[i];
		globalCells[i] = gcells[i];
		periodic[i] = cyclic[i];
	}

	// Full 3D cartesian domain
	strip[0].Initialize3D(domains,periodic);
	strip[0].Get_coords(3,domainIndex);
	strip[0].Shift(1,1,&n0[1],&n1[1]);
	strip[0].Shift(2,1,&n0[2],&n1[2]);
	strip[0].Shift(3,1,&n0[3],&n1[3]);
	finiteStrip[0].Initialize3D(domains,aperiodic);

	// strip domains
	strip[1].InitializeStrip(strip[0],1);
	strip[2].InitializeStrip(strip[0],2);
	strip[3].InitializeStrip(strip[0],3);
	finiteStrip[1].InitializeStrip(finiteStrip[0],1);
	finiteStrip[2].InitializeStrip(finiteStrip[0],2);
	finiteStrip[3].InitializeStrip(finiteStrip[0],3);

	for (i=1;i<=3;i++)
	{
		localCells[i] = globalCells[i] / domains[i];
		cornerCell[i] = 1 + domainIndex[i]*localCells[i];
		if (localCells[i]==1)
			localCells2[i] = 1;
		else
			localCells2[i] = localCells[i]+2;
	}

	// Thread affinity for internal MPI threads

	#ifdef USE_TW_MPI
	if (affinityMask.size())
	{
		int mpi_rank,cpu_to_bind;
		tw::Thread *mpi_thread;
		MPI_Comm_rank(MPI_COMM_WORLD,&mpi_rank);
		cpu_to_bind = affinityMask[mpi_rank];
		mpi_thread = TW_MPI_ThreadRef(MPI_COMM_WORLD,mpi_rank);
		mpi_thread->SetAffinity(cpu_to_bind);
	}
	#endif

	#ifdef USE_OPENCL

	tw::Int whichDevice = -1;
	tw::Int whichPlatform = -1;
	std::vector<tw::Int> deviceMap;
	std::valarray<cl_platform_id> platforms;
	std::valarray<cl_device_id> devices;
	std::string name;
	std::stringstream messg;
	char buff[1024];
	size_t buffSize;
	cl_uint num;
	cl_int err;

	// Choose OpenCL Platform

	clGetPlatformIDs(0,NULL,&num); // max count needs to be zero for AMDAPP
	platforms.resize(num);
	clGetPlatformIDs(num,&platforms[0],NULL);

	for ( i=0 ; i<platforms.size() ; i++)
	{
		clGetPlatformInfo(platforms[i],CL_PLATFORM_NAME,sizeof(buff),&buff,&buffSize);
		name = std::string(buff,buffSize-1);
		std::transform(name.begin(),name.end(),name.begin(),::tolower);
		std::transform(platformSearchString.begin(),platformSearchString.end(),platformSearchString.begin(),::tolower);
		if (name.find(platformSearchString)!=std::string::npos)
			whichPlatform = i;
	}

	if (whichPlatform==-1)
	{
		messg << "WARNING: could not find platform search string '" << platformSearchString << "'" << std::endl;
		messg << "    (using first available)" << std::endl;
		whichPlatform = 0;
	}

	for ( i=0 ; i<platforms.size() ; i++)
	{
		messg << "Platform " << i+1 << ": ";
		clGetPlatformInfo(platforms[i],CL_PLATFORM_NAME,sizeof(buff),&buff,&buffSize);
		name = std::string(buff,buffSize-1);
		messg << name;
		if (whichPlatform==i)
			messg << " *";
		messg << std::endl;
	}

	// Choose OpenCL Device

	clGetDeviceIDs(platforms[whichPlatform],CL_DEVICE_TYPE_ALL,0,NULL,&num); // max count needs to be zero for AMDAPP
	devices.resize(num);
	clGetDeviceIDs(platforms[whichPlatform],CL_DEVICE_TYPE_ALL,num,&devices[0],NULL);

	// If device ID list is given, device map is directly supplied in input file
	if (deviceIDList.size()>0)
	{
		for ( i=0 ; i<deviceIDList.size() ; i++ )
			deviceMap.push_back(deviceIDList[i]);
	}
	// Otherwise, form a device map based on searching device names
	else
	{
		for ( i=0 ; i<devices.size() ; i++ )
		{
			clGetDeviceInfo(devices[i],CL_DEVICE_NAME,sizeof(buff),&buff,&buffSize);
			name = std::string(buff,buffSize-1);
			std::transform(name.begin(),name.end(),name.begin(),::tolower);
			std::transform(deviceSearchString.begin(),deviceSearchString.end(),deviceSearchString.begin(),::tolower);
			if (name.find(deviceSearchString)!=std::string::npos)
				deviceMap.push_back(i);
		}
	}

	if (deviceMap.size()==0)
	{
		messg << "WARNING: could not form a device map." << std::endl;
		messg << "Search string was '" << deviceSearchString << "'" << std::endl;
		messg << "Size of device ID list was " << deviceIDList.size() << std::endl;
		messg << "    (will try to use device number = MPI rank)" << std::endl;
		whichDevice = strip[0].Get_rank();
	}
	else
		whichDevice = deviceMap[strip[0].Get_rank() % deviceMap.size()];

	for ( i=0 ; i<devices.size() ; i++)
	{
		// Print all devices and indicate the one chosen for this Task
		messg << "Device " << i+1 << ": ";
		clGetDeviceInfo(devices[i],CL_DEVICE_NAME,sizeof(buff),&buff,&buffSize);
		name = std::string(buff,buffSize-1);
		messg << name;
		if (whichDevice==i)
			messg << " *";
		messg << std::endl;
	}

	// Set members gpu, context, commandQueue

	gpu = devices[whichDevice];
	context = clCreateContext(NULL,1,&gpu,NULL,NULL,&err);
	// Following is the updated function
	//commandQueue = clCreateCommandQueueWithProperties(context,gpu,NULL,&err);
	// Following is the deprecated function
	commandQueue = clCreateCommandQueue(context,gpu,0,&err);

	// Create kernels
	std::string buildLog;
	InitializeCLProgram(fieldProgram,"3dfields.cl",buildLog);
	if (buildLog.size()>4)
	{
		messg << "WARNING : Build log for 3dfields.cl is not empty:" << std::endl;
		messg << buildLog << std::endl;
	}

	k_fillVec4Field = clCreateKernel(fieldProgram,"FillVec4Field",&err);
	k_fieldToBoundary = clCreateKernel(fieldProgram,"FieldToBoundaryData",&err);
	k_boundaryToField = clCreateKernel(fieldProgram,"BoundaryDataToField",&err);
	k_ghostToField = clCreateKernel(fieldProgram,"GhostDataToField",&err);
	k_zeroGhostCells = clCreateKernel(fieldProgram,"ZeroGhostCells",&err);
	k_swapBuffers = clCreateKernel(fieldProgram,"SwapBuffers",&err);
	k_destructiveSum = clCreateKernel(fieldProgram,"DestructiveSum",&err);
	k_destructiveNorm1 = clCreateKernel(fieldProgram,"DestructiveNorm1",&err);
	k_complexMod2 = clCreateKernel(fieldProgram,"ComplexMod2",&err);
	k_destructiveComplexMod2 = clCreateKernel(fieldProgram,"DestructiveComplexMod2",&err);
	k_MADD = clCreateKernel(fieldProgram,"MADD",&err);
	k_weightByVolume = clCreateKernel(fieldProgram,"WeightByVolume",&err);

	initMessage = messg.str();

	#endif
}

#ifdef USE_OPENCL

void Task::InitializeCLProgram(cl_program& program,const std::string& fileName,std::string& buildLog)
{
	// assumes context and device already created

	cl_int err;
	size_t buffSize;
	std::valarray<char> buildLogBuff;
	std::valarray<char*> sourceList(1);
	std::valarray<char> sourceText;

	std::ifstream theFile ( fileName.c_str() );
	std::string	sourceString ( std::istreambuf_iterator<char>(theFile),( std::istreambuf_iterator<char>() ) );
	theFile.close();
	std::string header_and_source = CLDefinitions(gpu) + sourceString;
	sourceText.resize(header_and_source.size()+1);
	strcpy(&sourceText[0],header_and_source.c_str());
	sourceList[0] = &sourceText[0];

	program = clCreateProgramWithSource(context,1,(const char**)&sourceList[0],NULL,&err);
	clBuildProgram(program,1,&gpu,"-cl-strict-aliasing",NULL,NULL);
	clGetProgramBuildInfo(program,gpu,CL_PROGRAM_BUILD_LOG,0,NULL,&buffSize);
	buildLogBuff.resize(buffSize);
	clGetProgramBuildInfo(program,gpu,CL_PROGRAM_BUILD_LOG,buffSize,&buildLogBuff[0],NULL);
	buildLog = std::string(&buildLogBuff[0],buffSize);
}

#endif

void Task::ReadData(std::ifstream& inFile)
{
	inFile.read((char*)globalCells,sizeof(tw::Int)*4);
	inFile.read((char*)domains,sizeof(tw::Int)*4);
	inFile.read((char*)periodic,sizeof(tw::Int)*4);
}

void Task::WriteData(std::ofstream& outFile)
{
	outFile.write((char*)globalCells,sizeof(tw::Int)*4);
	outFile.write((char*)domains,sizeof(tw::Int)*4);
	outFile.write((char*)periodic,sizeof(tw::Int)*4);
}
