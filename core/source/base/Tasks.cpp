module;

#include "tw_includes.h"

export module tasks;
import base;

export namespace tw
{
	/// @brief C++ wrapper for MPI communicator.  Differences follow:
	/// 1. We are starting spatial axes at index 1, consistent with TW framework (but not MPI)
	/// 2. There is special handling to enable send to self using blocking send/recv
	/// 3. aliased buffers are treated specially
	/// 4. ALL data-count arguments are given in bytes (overloading distinguishes int/float where necessary)
	struct comm
	{
		MPI_Comm comm_cart;
		MPI_Status status;
		MPI_Datatype tw_float_type;
		MPI_Datatype tw_int_type;
		void *sendToSelfData;

		comm()
		{
			comm_cart = MPI_COMM_WORLD;
			sendToSelfData = NULL;
			if (sizeof(tw::Float)==sizeof(float))
				tw_float_type = MPI_FLOAT;
			if (sizeof(tw::Float)==sizeof(double))
				tw_float_type = MPI_DOUBLE;
			tw_int_type = MPI_SHORT;
			if (sizeof(tw::Int)==sizeof(int))
				tw_int_type = MPI_INT;
			if (sizeof(tw::Int)==sizeof(long))
				tw_int_type = MPI_LONG;
			if (sizeof(tw::Int)==sizeof(long long))
				tw_int_type = MPI_LONG_LONG_INT;
			if (tw_int_type==MPI_SHORT)
				throw tw::FatalError("Unsupported tw::Int type.");
		}
		~comm()
		{
			if (comm_cart!=MPI_COMM_WORLD)
				MPI_Comm_free(&comm_cart);
		}
		tw::Int Get_size() const
		{
			int ans;
			MPI_Comm_size(comm_cart,&ans);
			return ans;
		}
		tw::Int Get_rank() const
		{
			int ans;
			MPI_Comm_rank(comm_cart,&ans);
			return ans;
		}
		tw::Int Cart_rank(tw::Int i,tw::Int j,tw::Int k)
		{
			int ans;
			int coords[3] = {(int)i,(int)j,(int)k};
			MPI_Cart_rank(comm_cart,coords,&ans);
			return ans;
		}
		void Initialize3D(const node4& domains,const node4& periodic)
		{
			int d[3] = { (int)domains[1],(int)domains[2],(int)domains[3] };
			int p[3] = { (int)periodic[1],(int)periodic[2],(int)periodic[3] };
			MPI_Cart_create(MPI_COMM_WORLD,3,d,p,0,&comm_cart);
		}
		void InitializeStrip(comm& fullGrid,tw::Int axis)
		{
			int remain_dims[3] = { axis==1 , axis==2 , axis==3 };
			MPI_Cart_sub(fullGrid.comm_cart,remain_dims,&comm_cart);
		}
		void Clear()
		{
			if (comm_cart!=MPI_COMM_WORLD)
				MPI_Comm_free(&comm_cart);
			comm_cart = MPI_COMM_WORLD;
		}
		/// n.b. the first axis is 1, and for strip communicators the only axis is 1
		void Shift(tw::Int axis,tw::Int shift,tw::Int *src,tw::Int *dst) const
		{
			int srci,dsti;
			MPI_Cart_shift(comm_cart,axis-1,shift,&srci,&dsti);
			*src = srci;
			*dst = dsti;
		}
		/// get 4d domain indices for this rank, time will always be 0
		node4 Get_coords4()
		{
			return Get_coords4(Get_rank());
		}
		/// get 4d domain indices for other rank, time will always be 0
		node4 Get_coords4(tw::Int rank)
		{
			node4 ans;
			int temp[4];
			MPI_Cart_coords(comm_cart,rank,3,&temp[1]);
			ans[0] = 0;
			for (auto i=1; i<4; i++)
				ans[i] = temp[i];
			return ans;
		}
		void Get_periods(bool *periods)
		{
			int i,d[3],c[3],temp[3];
			MPI_Cart_get(comm_cart,3,d,temp,c);
			for (i=0;i<3;i++)
				periods[i+1] = temp[i];
		}
		void Send(void *buff,tw::Int buffSize,tw::Int dst)
		{
			if (dst==Get_rank())
				if (sendToSelfData!=NULL)
				{
					if (buff==sendToSelfData)
						std::cout << term::warning << ": memcpy in place during send to self." << std::endl;
					std::memcpy(sendToSelfData,buff,buffSize);
					sendToSelfData = NULL;
				}
				else
				{
					sendToSelfData = buff;
				}
			else
				MPI_Send(buff,buffSize,MPI_BYTE,dst,0,comm_cart);
		}
		void Recv(void *buff,tw::Int buffSize,tw::Int src)
		{
			if (src==Get_rank())
				if (sendToSelfData!=NULL)
				{
					if (buff==sendToSelfData)
						std::cout << term::warning << ": memcpy in place during recv from self." << std::endl;
					std::memcpy(buff,sendToSelfData,buffSize);
					sendToSelfData = NULL;
				}
				else
				{
					sendToSelfData = buff;
				}
			else
				MPI_Recv(buff,buffSize,MPI_BYTE,src,0,comm_cart,&status);
		}
		void Sum(void *sb,void *rb,tw::Int count,tw::Int root)
		{
			if (sb==rb && root==Get_rank())
				MPI_Reduce(MPI_IN_PLACE,rb,count/sizeof(tw::Float),tw_float_type,MPI_SUM,root,comm_cart);
			else
				MPI_Reduce(sb,rb,count/sizeof(tw::Float),tw_float_type,MPI_SUM,root,comm_cart);
		}
		void Bcast(void *buf,tw::Int count,tw::Int root)
		{
			MPI_Bcast(buf,count,MPI_BYTE,root,comm_cart);
		}
		void AllSum(void *sb,void *rb,tw::Int count,tw::Int root)
		{
			Sum(sb,rb,count,root);
			MPI_Bcast(rb,count,MPI_BYTE,root,comm_cart);
		}
		tw::Int GetMin(tw::Int test)
		{
			tw::Int ans;
			MPI_Reduce(&test,&ans,1,tw_int_type,MPI_MIN,0,comm_cart);
			MPI_Bcast(&ans,1,tw_int_type,0,comm_cart);
			return ans;
		}
		tw::Float GetMin(tw::Float test)
		{
			tw::Float ans;
			MPI_Reduce(&test,&ans,1,tw_float_type,MPI_MIN,0,comm_cart);
			MPI_Bcast(&ans,1,tw_float_type,0,comm_cart);
			return ans;
		}
		tw::Int GetMax(tw::Int test)
		{
			tw::Int ans;
			MPI_Reduce(&test,&ans,1,tw_int_type,MPI_MAX,0,comm_cart);
			MPI_Bcast(&ans,1,tw_int_type,0,comm_cart);
			return ans;
		}
		tw::Float GetMax(tw::Float test)
		{
			tw::Float ans;
			MPI_Reduce(&test,&ans,1,tw_float_type,MPI_MAX,0,comm_cart);
			MPI_Bcast(&ans,1,tw_float_type,0,comm_cart);
			return ans;
		}
		void Gather(void *sb,void *rb,tw::Int count,tw::Int root)
		{
			if ((char*)sb==(char*)rb+count*root && root==Get_rank())
				MPI_Gather(MPI_IN_PLACE,count,MPI_BYTE,rb,count,MPI_BYTE,root,comm_cart);
			else
				MPI_Gather(sb,count,MPI_BYTE,rb,count,MPI_BYTE,root,comm_cart);
		}
		void Scatter(void *sb,void *rb,tw::Int count,tw::Int root)
		{
			if ((char*)sb+count*root==(char*)rb && root==Get_rank())
				MPI_Scatter(sb,count,MPI_BYTE,MPI_IN_PLACE,count,MPI_BYTE,root,comm_cart);
			else
				MPI_Scatter(sb,count,MPI_BYTE,rb,count,MPI_BYTE,root,comm_cart);
		}
		void Barrier()
		{
			MPI_Barrier(comm_cart);
		}
	};
}

export struct Task
{
	tw::comm strip[4]; // all,x,y,z
	tw::comm finiteStrip[4]; // same as strip, but ignoring all periodicity
	tw::node4 domains,periodic,domainIndex; // blank,x,y,z
	tw::node4 n0,n1; // low and high neighbors : blank,x,y,z

	std::string unitTest,inputFileName,restartFileName; // keep with task so ComputeTool can access

	std::string deviceSearchString,platformSearchString,initMessage;

	UniformDeviate *uniformDeviate;
	GaussianDeviate *gaussianDeviate;

	#ifdef USE_OPENCL
	cl_program fieldProgram;
	cl_kernel k_fieldToBoundary,k_boundaryToField,k_ghostToField,k_zeroGhostCells;
	cl_kernel k_fillVec4Field,k_swapBuffers;
	cl_kernel k_destructiveSum,k_destructiveNorm1;
	cl_kernel k_MADD,k_weightByVolume;
	cl_kernel k_destructiveComplexMod2,k_complexMod2;

	cl_device_id gpu;
	cl_context context;
	cl_command_queue commandQueue;
	void InitializeCLProgram(cl_program& prog,const std::string& file,std::string& buildLog);
	#endif

	Task();
	virtual ~Task();
	void Initialize(const tw::node4& doms,const tw::node4& cyclic);
	tw::Int NumTasks() { return domains[1]*domains[2]*domains[3]; }
};

////////////////////
//      TASK      //
////////////////////


Task::Task()
{
	for (auto i=0;i<4;i++) {
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

void Task::Initialize(const tw::node4& doms,const tw::node4& cyclic)
{
	tw::node4 aperiodic { 0,0,0,0 };

	for (auto i=0;i<4;i++)
	{
		domains[i] = doms[i];
		periodic[i] = cyclic[i];
	}

	// Full 3D cartesian domain
	strip[0].Initialize3D(domains,periodic);
	domainIndex = strip[0].Get_coords4();
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
		messg << term::warning << ": could not find platform search string '" << platformSearchString << "'" << std::endl;
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

	// First try to interpret search string as a device number list
	std::string deviceList(deviceSearchString);
	do
	{
		if (deviceList.find(",")!=std::string::npos)
			deviceList.replace(deviceList.find(","),1," ");
	} while (deviceList.find(",")!=std::string::npos);
	std::stringstream deviceStream(deviceList);
	do
	{
		tw::Int deviceNum;
		if (deviceStream >> deviceNum)
		{
			if (deviceNum>=0 && deviceNum<devices.size())
				deviceMap.push_back(deviceNum);
			else
				deviceMap.push_back(-1);
		}
		else
			deviceMap.push_back(-1);
	} while(deviceStream.good());
	if (deviceMap.back()==-1)
		deviceMap.clear();

	// If it didn't work assume we are looking for a name
	if (deviceMap.size()==0)
		for ( i=0 ; i<devices.size() ; i++ )
		{
			clGetDeviceInfo(devices[i],CL_DEVICE_NAME,sizeof(buff),&buff,&buffSize);
			name = std::string(buff,buffSize-1);
			std::transform(name.begin(),name.end(),name.begin(),::tolower);
			std::transform(deviceSearchString.begin(),deviceSearchString.end(),deviceSearchString.begin(),::tolower);
			if (name.find(deviceSearchString)!=std::string::npos)
				deviceMap.push_back(i);
		}

	if (deviceMap.size()==0)
	{
		messg << term::warning << ": could not form a device map." << std::endl;
		messg << "Search string was '" << deviceSearchString << "'" << std::endl;
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
	InitializeCLProgram(fieldProgram,"fields.cl",buildLog);
	if (buildLog.size()>4)
	{
		messg << "WARNING : Build log for fields.cl is not empty:" << std::endl;
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
	// Substitute in turboWAVE-OpenCL protocols.  These are string pairs that match a key to a
	// frequently recurring code fragment.
	std::map<std::string,std::string> protocols = CLProtocols();
	for (auto p : protocols)
	{
		std::string::size_type f;
		do {
			f = sourceString.find(p.first);
			if (f!=std::string::npos)
				sourceString.replace(f,p.first.size(),p.second);
		} while(f!=std::string::npos);
	}
	std::string header_and_source = CLDefinitions(gpu) + sourceString;
	sourceText.resize(header_and_source.size()+1);
	// Following is the updated function
	//strcpy_s(&sourceText[0],sourceText.size(),header_and_source.c_str());
	// Following is the deprecated function
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
