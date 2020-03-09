// The following "tw::comm" object is, in essence, a C++ wrapper for the MPI communicator object
// Differences from MPI::Comm follow:
// 1. We are starting spatial axes at index 1, consistent with TW framework (but not MPI)
// 2. There is special handling to enable send to self using blocking send/recv
// 3. aliased buffers are treated specially
// 4. ALL data-count arguments are given in bytes (overloading distinguishes int/float where necessary)

namespace tw
{
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
		tw::Int Get_size()
		{
			int ans;
			MPI_Comm_size(comm_cart,&ans);
			return ans;
		}
		tw::Int Get_rank()
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
		void Initialize3D(tw::Int *domains,tw::Int *periodic)
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
		void Shift(tw::Int axis,tw::Int shift,tw::Int *src,tw::Int *dst)
		{
			int srci,dsti;
			MPI_Cart_shift(comm_cart,axis-1,shift,&srci,&dsti);
			*src = srci;
			*dst = dsti;
		}
		void Get_coords(tw::Int max_dims,tw::Int *xyz)
		{
			int i,temp[3];
			MPI_Cart_coords(comm_cart,Get_rank(),max_dims,temp);
			for (i=0;i<max_dims;i++)
				xyz[i+1] = temp[i];
		}
		void Get_coords(tw::Int max_dims,tw::Int rank,tw::Int *xyz)
		{
			int i,temp[3];
			MPI_Cart_coords(comm_cart,rank,max_dims,temp);
			for (i=0;i<max_dims;i++)
				xyz[i+1] = temp[i];
		}
		void Get_periods(tw::Int max_dims,bool *periods)
		{
			int i,d[3],c[3],temp[3];
			MPI_Cart_get(comm_cart,max_dims,d,temp,c);
			for (i=0;i<max_dims;i++)
				periods[i+1] = temp[i];
		}
		void Send(void *buff,tw::Int buffSize,tw::Int dst)
		{
			if (dst==Get_rank())
				if (sendToSelfData!=NULL)
				{
					if (buff==sendToSelfData)
						std::cout << "WARNING: memcpy in place during send to self." << std::endl;
					memcpy(sendToSelfData,buff,buffSize);
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
						std::cout << "WARNING: memcpy in place during recv from self." << std::endl;
					memcpy(buff,sendToSelfData,buffSize);
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

struct Task
{
	tw::comm strip[4]; // all,x,y,z
	tw::comm finiteStrip[4]; // same as strip, but ignoring all periodicity
	tw::Int domains[4],globalCells[4],periodic[4]; // blank,x,y,z
	tw::Int cornerCell[4],localCells[4],localCells2[4],domainIndex[4]; // blank,x,y,z
	tw::Int n0[4],n1[4]; // low and high neighbors : blank,x,y,z

	std::string deviceSearchString,platformSearchString,initMessage;
	std::valarray<tw::Int> deviceIDList; // accelerator devices such as GPGPU
	std::valarray<tw::Int> affinityMask; // IDs of logical processors for thread binding

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

	void Initialize(tw::Int *doms,tw::Int *gcells,tw::Int *cyclic);
	tw::Int LocalCellIndex(tw::Int i,tw::Int axis) { return i - cornerCell[axis] + 1; }
	tw::Int GlobalCellIndex(tw::Int i,tw::Int axis) { return i + cornerCell[axis] - 1; }

	tw::Int NumTasks()
	{
		return domains[1]*domains[2]*domains[3];
	}

	void ReadCheckpoint(std::ifstream& inFile);
	void WriteCheckpoint(std::ofstream& outFile);

	void Lock()
	{
		#ifdef USING_TW_MPI
		TW_MPI_Lock();
		#endif
	}
	void Unlock()
	{
		#ifdef USING_TW_MPI
		TW_MPI_Unlock();
		#endif
	}
	void Complete()
	{
		#ifdef USING_TW_MPI
		TW_MPI_ThreadRef(strip[0].comm_cart,strip[0].Get_rank())->Complete();
		#endif
	}
};
