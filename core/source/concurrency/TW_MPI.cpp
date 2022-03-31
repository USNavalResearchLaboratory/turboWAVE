#include "definitions.h"

// ------------------------------------------------------------------------------------
// TurboWAVE non-conforming internal implementation of MPI
// By no means is this a portable MPI implementation.  Departures include:
// * No need to launch using mpirun
// 	- Single startup process launches threads based on command line argument "-n [NUM]"
// * Always use "eager" protocol
// * Maximum communicator dimension = 3
// * Require contiguous data types
// * Tags are ignored
// * Minimal error checking , MPI_Status not updated
// * Subset of reduce operations (min , max , sum)
// * Subset of MPI functions
// * Gather/Scatter have to be homogeneous
// * NO PROMISE THAT THIS LIST IS EXHAUSTIVE
// * Over time we may try to relax assumptions
// ------------------------------------------------------------------------------------

// Here are the MPI_COMM_WORLD objects, one for each thread
// Note MPI_COMM_WORLD = NULL should be treated as a label

MPI_Comm_object comm_world[MAX_THREADS];
tw::mutex masterLock;

// Functions specific to TurboWAVE MPI

void TW_MPI_Lock()
{
	masterLock.Lock();
}

void TW_MPI_Unlock()
{
	masterLock.Unlock();
}

tw::Thread* TW_MPI_ThreadRef(MPI_Comm comm,int rank)
{
	int coord[4];
	int i,wrank;
	if (comm==MPI_COMM_WORLD)
		comm = TW_MPI_FindCommWorld();

	MPI_Cart_coords(comm,rank,comm->dims,&coord[1]);
	wrank = comm->wstride[0];
	for (i=1;i<=comm->dims;i++)
		wrank += comm->wstride[i]*coord[i];
	return comm_world[wrank].threadRef;
}

MPI_Comm TW_MPI_FindCommWorld()
{
	tw::Int i;
	for (i=0;i<comm_world[0].domains[1];i++)
		if (comm_world[i].threadRef->IsCallingThread())
			return &comm_world[i];
	std::cout << "ERROR: Could not find this thread's comm world." << std::endl;
	return MPI_COMM_WORLD;
}

void TW_MPI_Launch(std::vector<tw::Thread*>& threadList)
{
	tw::Int i;
	// Set up comm world for each thread
	for (i=0;i<threadList.size();i++)
	{
		comm_world[i].rank = i;
		comm_world[i].wrank = i;
		comm_world[i].dims = 1;
		comm_world[i].domains[1] = threadList.size();
		comm_world[i].periodic[1] = 0;
		comm_world[i].stride[0] = 0;
		comm_world[i].stride[1] = 1;
		comm_world[i].wstride[0] = 0;
		comm_world[i].wstride[1] = 1;
		comm_world[i].threadRef = threadList[i];
	}
	// Start all the threads
	for (i=0;i<threadList.size();i++)
		threadList[i]->Start();
}

// Standard MPI functions

int MPI_Init(int *argc, char ***argv)
{
	return 0;
}

int MPI_Finalize(void)
{
	return 0;
}

int MPI_Cart_rank(MPI_Comm comm,const int coords[],int *rank)
{
	int i;
	*rank = 0;
	if (comm==MPI_COMM_WORLD)
		comm = TW_MPI_FindCommWorld();

	for (i=1;i<=comm->dims;i++)
		*rank += coords[i-1]*comm->stride[i];
	return 0;
}

int MPI_Cart_coords(MPI_Comm comm,int rank,int maxdims,int *coords)
{
	if (comm==MPI_COMM_WORLD)
		comm = TW_MPI_FindCommWorld();

	if (maxdims==3)
	{
		coords[2] = rank/comm->stride[3];
		coords[1] = (rank-coords[2]*comm->stride[3])/comm->stride[2];
		coords[0] = (rank-coords[2]*comm->stride[3]-coords[1]*comm->stride[2])/comm->stride[1];
	}
	if (maxdims==2)
	{
		coords[1] = rank/comm->stride[2];
		coords[0] = (rank-coords[1]*comm->stride[2])/comm->stride[1];
	}
	if (maxdims==1)
	{
		coords[0] = rank/comm->stride[1];
	}
	return 0;
}

int MPI_Cart_get(MPI_Comm comm,int maxdims,int *dims,int *periods,int *coords)
{
	int i,r;
	if (comm==MPI_COMM_WORLD)
		comm = TW_MPI_FindCommWorld();
	for (i=0;i<maxdims;i++)
	{
		dims[i] = comm->domains[i+1];
		periods[i] = comm->periodic[i+1];
	}
	MPI_Comm_rank(comm,&r);
	MPI_Cart_coords(comm,r,maxdims,coords);
	return 0;
}

int MPI_Cart_create(MPI_Comm oldComm,int maxdims,int *domains,int *periodic,int reorder,MPI_Comm *newComm)
{
	int i;
	if (oldComm==MPI_COMM_WORLD)
		oldComm = TW_MPI_FindCommWorld();

	*newComm = (MPI_Comm)malloc(sizeof(MPI_Comm_object));
	(*newComm)->dims = maxdims;
	for (i=1;i<=maxdims;i++)
	{
		(*newComm)->domains[i] = domains[i-1];
		(*newComm)->periodic[i] = periodic[i-1];
	}
	(*newComm)->stride[0] = 0;
	(*newComm)->stride[1] = 1;
	(*newComm)->stride[2] = domains[0];
	(*newComm)->stride[3] = domains[0]*domains[1];
	(*newComm)->wstride[0] = 0;
	(*newComm)->wstride[1] = 1;
	(*newComm)->wstride[2] = domains[0];
	(*newComm)->wstride[3] = domains[0]*domains[1];
	(*newComm)->rank = oldComm->rank;
	(*newComm)->wrank = oldComm->wrank;
	(*newComm)->threadRef = oldComm->threadRef;
	return 0;
}

int MPI_Cart_sub(MPI_Comm oldComm,int *remain_dims,MPI_Comm *newComm)
{
	int i,d=0;
	int old_coords[4],new_coords[4];
	if (oldComm==MPI_COMM_WORLD)
		oldComm = TW_MPI_FindCommWorld();

	MPI_Cart_coords(oldComm,oldComm->rank,oldComm->dims,&old_coords[1]);

	*newComm = (MPI_Comm)malloc(sizeof(MPI_Comm_object));
	(*newComm)->wstride[0] = oldComm->wstride[0];
	(*newComm)->stride[0] = 0;
	for (i=1;i<=oldComm->dims;i++)
	{
		if (remain_dims[i-1])
		{
			d++;
			(*newComm)->periodic[d] = oldComm->periodic[i];
			(*newComm)->domains[d] = oldComm->domains[i];
			(*newComm)->wstride[d] = oldComm->wstride[i];
			new_coords[d] = old_coords[i];
			if (d==1)
				(*newComm)->stride[d] = 1;
			else
				(*newComm)->stride[d] = (*newComm)->stride[d-1] * (*newComm)->domains[d-1];
		}
		else
		{
			(*newComm)->wstride[0] += oldComm->wstride[i]*old_coords[i];
		}
	}
	(*newComm)->dims = d;
	(*newComm)->threadRef = oldComm->threadRef;

	(*newComm)->wrank = oldComm->wrank;
	MPI_Cart_rank(*newComm,&new_coords[1],&(*newComm)->rank);

	return 0;
}

int MPI_Cart_shift(MPI_Comm comm,int direction,int displ,int *src,int *dst)
{
	if (comm==MPI_COMM_WORLD)
		comm = TW_MPI_FindCommWorld();

	int src_coords[4],dst_coords[4];
	int a = direction+1;
	MPI_Cart_coords(comm,comm->rank,comm->dims,&src_coords[1]);
	MPI_Cart_coords(comm,comm->rank,comm->dims,&dst_coords[1]);

	src_coords[a] -= displ;
	dst_coords[a] += displ;

	if (dst_coords[a] >= comm->domains[a] || dst_coords[a]<0)
	{
		if (comm->periodic[a])
		{
			dst_coords[a] -= (displ/abs(displ))*comm->domains[a];
			MPI_Cart_rank(comm,&dst_coords[1],dst);
		}
		else
			*dst = MPI_PROC_NULL;
	}
	else
		MPI_Cart_rank(comm,&dst_coords[1],dst);

	if (src_coords[a] >= comm->domains[a] || src_coords[a]<0)
	{
		if (comm->periodic[a])
		{
			src_coords[a] += (displ/abs(displ))*comm->domains[a];
			MPI_Cart_rank(comm,&src_coords[1],src);
		}
		else
			*src = MPI_PROC_NULL;
	}
	else
		MPI_Cart_rank(comm,&src_coords[1],src);

	return 0;
}

int MPI_Send(void* buf, int count, MPI_Datatype datatype, int dest,int tag, MPI_Comm comm)
{
	if (comm==MPI_COMM_WORLD)
		comm = TW_MPI_FindCommWorld();

//	std::cout << comm->rank << " -> " << dest << std::endl;
	if (dest==MPI_PROC_NULL)
		return 0;

	if (dest==comm->rank)
	{
		comm->threadRef->self.Notify(buf,comm->threadRef);
		return 0;
	}

	tw::Thread* dst_thread = TW_MPI_ThreadRef(comm,dest);

	dst_thread->input.Notify(buf,comm->threadRef);
	dst_thread->output.Wait(NULL,comm->threadRef,-1);

	return 0;
}

int MPI_Recv(void* buf, int count, MPI_Datatype dt, int source,int tag, MPI_Comm comm, MPI_Status *status)
{
	if (comm==MPI_COMM_WORLD)
		comm = TW_MPI_FindCommWorld();

//	std::cout << comm->rank << " <- " << source << std::endl;
	void *srcData;
	int sizeofdata = TW_MPI_GetDataSize(dt);

	if (source==MPI_PROC_NULL)
		return 0;

	if (source==comm->rank)
	{
		comm->threadRef->self.Wait(&srcData,comm->threadRef,-1);
		memcpy(buf,srcData,count*sizeofdata);
		return 0;
	}

	tw::Thread* src_thread = TW_MPI_ThreadRef(comm,source);

	comm->threadRef->input.Wait(&srcData,src_thread,-1);
	memcpy(buf,srcData,count*sizeofdata);
	comm->threadRef->output.Notify(NULL,src_thread);

	return 0;
}

int MPI_Reduce(void *sb,void *rb,int count,MPI_Datatype dt,MPI_Op op,int root,MPI_Comm comm)
{
	if (comm==MPI_COMM_WORLD)
		comm = TW_MPI_FindCommWorld();

	int i,j,num;
	MPI_Status status;
	int sizeofdata = TW_MPI_GetDataSize(dt);
	MPI_Comm_size(comm,&num);
	if (comm->rank==root)
	{
		void *temp = malloc(count*sizeofdata);
		if (sb!=MPI_IN_PLACE)
			memcpy(rb,sb,count*sizeofdata);
		for (i=0;i<num;i++)
			if (i!=root)
			{
				MPI_Recv(temp,count,dt,i,0,comm,&status);
				switch (op)
				{
					case MPI_MIN:
						for (j=0;j<count;j++)
							TW_MPI_Minimize_Datatype(j,rb,temp,rb,dt);
						break;
					case MPI_MAX:
						for (j=0;j<count;j++)
							TW_MPI_Maximize_Datatype(j,rb,temp,rb,dt);
						break;
					case MPI_SUM:
						for (j=0;j<count;j++)
							TW_MPI_Add_Datatype(j,rb,temp,rb,dt);
						break;
				}
			}
		free(temp);
	}
	else
	{
		MPI_Send(sb,count,dt,root,0,comm);
	}
	return 0;
}

int MPI_Bcast(void *buf,int count,MPI_Datatype dt,int root,MPI_Comm comm)
{
	if (comm==MPI_COMM_WORLD)
		comm = TW_MPI_FindCommWorld();

	int i,num;
	MPI_Status status;
	MPI_Comm_size(comm,&num);
	if (comm->rank==root)
	{
		for (i=0;i<num;i++)
			if (i!=root)
				MPI_Send(buf,count,dt,i,0,comm);
	}
	else
	{
		MPI_Recv(buf,count,dt,root,0,comm,&status);
	}

	return 0;
}

int MPI_Gather(void *sb,int scount,MPI_Datatype sdt,void *rb,int rcount,MPI_Datatype rdt,int root,MPI_Comm comm)
{
	if (comm==MPI_COMM_WORLD)
		comm = TW_MPI_FindCommWorld();

	int i,j,num;
	MPI_Status status;
	int rsize = TW_MPI_GetDataSize(rdt);
	int ssize = TW_MPI_GetDataSize(sdt);
	MPI_Comm_size(comm,&num);
	if (comm->rank==root)
	{
		if (sb!=MPI_IN_PLACE)
			memcpy(TW_MPI_Offset_Pointer(root*rcount,rb,rdt),sb,scount*ssize);
		for (i=0;i<num;i++)
			if (i!=root)
				MPI_Recv(TW_MPI_Offset_Pointer(i*rcount,rb,rdt),rcount,rdt,i,0,comm,&status);
	}
	else
	{
		MPI_Send(sb,scount,sdt,root,0,comm);
	}
	return 0;
}

int MPI_Scatter(void *sb,int scount,MPI_Datatype sdt,void *rb,int rcount,MPI_Datatype rdt,int root,MPI_Comm comm)
{
	if (comm==MPI_COMM_WORLD)
		comm = TW_MPI_FindCommWorld();

	int i,j,num;
	MPI_Status status;
	int rsize = TW_MPI_GetDataSize(rdt);
	int ssize = TW_MPI_GetDataSize(sdt);
	MPI_Comm_size(comm,&num);
	if (comm->rank==root)
	{
		if (rb!=MPI_IN_PLACE)
			memcpy(rb,TW_MPI_Offset_Pointer(root*scount,sb,sdt),rcount*rsize);
		for (i=0;i<num;i++)
			if (i!=root)
				MPI_Send(TW_MPI_Offset_Pointer(i*scount,sb,sdt),scount,sdt,i,0,comm);
	}
	else
	{
		MPI_Recv(rb,rcount,rdt,root,0,comm,&status);
	}
	return 0;
}

int MPI_Barrier(MPI_Comm comm)
{
	if (comm==MPI_COMM_WORLD)
		comm = TW_MPI_FindCommWorld();
	char buf[255];
	MPI_Bcast(buf,1,MPI_BYTE,0,comm);
	return 0;
}
