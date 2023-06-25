// TurboWAVE non-conforming implementation of MPI

#define USE_TW_MPI

static const int MAX_THREADS = 256;
static const int MPI_PROC_NULL = -3;
static void * MPI_IN_PLACE = NULL;

enum MPI_Datatype { MPI_BYTE, MPI_SHORT, MPI_INT, MPI_LONG, MPI_LONG_LONG_INT, MPI_FLOAT, MPI_DOUBLE };
enum MPI_Op { MPI_MAX, MPI_MIN, MPI_SUM };

typedef struct {
   int source;
   int tag;
   int error;
   int len;
   int type;
} MPI_Status;

typedef struct {
	int rank,wrank;
	int dims;
	// spatial coordinate indices start at 1
	// index 0 is reserved for stride offset
	int domains[4];
	int periodic[4];
	int stride[4];
	int wstride[4];
	tw::Thread *threadRef;
} MPI_Comm_object;

typedef MPI_Comm_object * MPI_Comm;

static const MPI_Comm MPI_COMM_WORLD = NULL;

// TurboWAVE specific MPI helper functions (TW_MPI prefix)

inline void* TW_MPI_Offset_Pointer(int offset,void *ptr,MPI_Datatype dt)
{
	switch (dt)
	{
		case MPI_BYTE:
			return (char*)ptr+offset;
		case MPI_SHORT:
			return (short*)ptr+offset;
		case MPI_INT:
			return (int*)ptr+offset;
		case MPI_LONG:
			return (long*)ptr+offset;
		case MPI_LONG_LONG_INT:
			return (long long*)ptr+offset;
		case MPI_FLOAT:
			return (float*)ptr+offset;
		case MPI_DOUBLE:
			return (double*)ptr+offset;
	}
	return NULL;
}


inline int TW_MPI_GetDataSize(MPI_Datatype dt)
{
	switch (dt)
	{
		case MPI_BYTE:
			return 1;
		case MPI_SHORT:
			return sizeof(short);
		case MPI_INT:
			return sizeof(int);
		case MPI_LONG:
			return sizeof(long);
		case MPI_LONG_LONG_INT:
			return sizeof(long long);
		case MPI_FLOAT:
			return sizeof(float);
		case MPI_DOUBLE:
			return sizeof(double);
	}
	return 0;
}

inline void TW_MPI_Add_Datatype(int i,void *ans,void *x,void *y,MPI_Datatype dt)
{
	switch (dt)
	{
		case MPI_BYTE:
			*((char*)ans+i) = *((char*)x+i) + *((char*)y+i);
			break;
		case MPI_SHORT:
			*((short*)ans+i) = *((short*)x+i) + *((short*)y+i);
			break;
		case MPI_INT:
			*((int*)ans+i) = *((int*)x+i) + *((int*)y+i);
			break;
		case MPI_LONG:
			*((long*)ans+i) = *((long*)x+i) + *((long*)y+i);
			break;
		case MPI_LONG_LONG_INT:
			*((long long*)ans+i) = *((long long*)x+i) + *((long long*)y+i);
			break;
		case MPI_FLOAT:
			*((float*)ans+i) = *((float*)x+i) + *((float*)y+i);
			break;
		case MPI_DOUBLE:
			*((double*)ans+i) = *((double*)x+i) + *((double*)y+i);
			break;
	}
}

inline void TW_MPI_Minimize_Datatype(int i,void *ans,void *x,void *y,MPI_Datatype dt)
{
	switch (dt)
	{
		case MPI_BYTE:
			if (*((char*)x+i) < *((char*)y+i))
				*((char*)ans+i) = *((char*)x+i);
			else
				*((char*)ans+i) = *((char*)y+i);
			break;
		case MPI_SHORT:
			if (*((short*)x+i) < *((short*)y+i))
				*((short*)ans+i) = *((short*)x+i);
			else
				*((short*)ans+i) = *((short*)y+i);
			break;
		case MPI_INT:
			if (*((int*)x+i) < *((int*)y+i))
				*((int*)ans+i) = *((int*)x+i);
			else
				*((int*)ans+i) = *((int*)y+i);
			break;
		case MPI_LONG:
			if (*((long*)x+i) < *((long*)y+i))
				*((long*)ans+i) = *((long*)x+i);
			else
				*((long*)ans+i) = *((long*)y+i);
			break;
		case MPI_LONG_LONG_INT:
			if (*((long long*)x+i) < *((long long*)y+i))
				*((long long*)ans+i) = *((long long*)x+i);
			else
				*((long long*)ans+i) = *((long long*)y+i);
			break;
		case MPI_FLOAT:
			if (*((float*)x+i) < *((float*)y+i))
				*((float*)ans+i) = *((float*)x+i);
			else
				*((float*)ans+i) = *((float*)y+i);
			break;
		case MPI_DOUBLE:
			if (*((double*)x+i) < *((double*)y+i))
				*((double*)ans+i) = *((double*)x+i);
			else
				*((double*)ans+i) = *((double*)y+i);
			break;
	}
}

inline void TW_MPI_Maximize_Datatype(int i,void *ans,void *x,void *y,MPI_Datatype dt)
{
	switch (dt)
	{
		case MPI_BYTE:
			if (*((char*)x+i) > *((char*)y+i))
				*((char*)ans+i) = *((char*)x+i);
			else
				*((char*)ans+i) = *((char*)y+i);
			break;
		case MPI_SHORT:
			if (*((short*)x+i) > *((short*)y+i))
				*((short*)ans+i) = *((short*)x+i);
			else
				*((short*)ans+i) = *((short*)y+i);
			break;
		case MPI_INT:
			if (*((int*)x+i) > *((int*)y+i))
				*((int*)ans+i) = *((int*)x+i);
			else
				*((int*)ans+i) = *((int*)y+i);
			break;
		case MPI_LONG:
			if (*((long*)x+i) > *((long*)y+i))
				*((long*)ans+i) = *((long*)x+i);
			else
				*((long*)ans+i) = *((long*)y+i);
			break;
		case MPI_LONG_LONG_INT:
			if (*((long long*)x+i) > *((long long*)y+i))
				*((long long*)ans+i) = *((long long*)x+i);
			else
				*((long long*)ans+i) = *((long long*)y+i);
			break;
		case MPI_FLOAT:
			if (*((float*)x+i) > *((float*)y+i))
				*((float*)ans+i) = *((float*)x+i);
			else
				*((float*)ans+i) = *((float*)y+i);
			break;
		case MPI_DOUBLE:
			if (*((double*)x+i) > *((double*)y+i))
				*((double*)ans+i) = *((double*)x+i);
			else
				*((double*)ans+i) = *((double*)y+i);
			break;
	}
}

void TW_MPI_Lock();
void TW_MPI_Unlock();
MPI_Comm TW_MPI_FindCommWorld();
tw::Thread* TW_MPI_ThreadRef(MPI_Comm comm,int rank);
void TW_MPI_Launch(std::vector<tw::Thread*>& threadList);

// Standard MPI

int MPI_Init(int *argc, char ***argv);
int MPI_Finalize(void);

inline int MPI_Comm_size(MPI_Comm comm, int *size)
{
	int i;
	*size = 1;

	if (comm==MPI_COMM_WORLD)
		comm = TW_MPI_FindCommWorld();
	for (i=1;i<=comm->dims;i++)
		*size *= comm->domains[i];
	return 0;
}

inline int MPI_Comm_rank(MPI_Comm comm, int *rank)
{
	if (comm==MPI_COMM_WORLD)
		comm = TW_MPI_FindCommWorld();
	*rank = comm->rank;
	return 0;
}

inline int MPI_Comm_free(MPI_Comm *comm)
{
	if (*comm!=MPI_COMM_WORLD)
		free(*comm);
	return 0;
}

int MPI_Cart_create(MPI_Comm oldComm,int maxdims,int *domains,int *periodic,int reorder,MPI_Comm *newComm);
int MPI_Cart_sub(MPI_Comm oldComm,int *remain_dims,MPI_Comm *newComm);
int MPI_Cart_rank(MPI_Comm comm,const int coords[],int *rank);
int MPI_Cart_shift(MPI_Comm comm,int direction,int displ,int *src,int *dst);
int MPI_Cart_coords(MPI_Comm comm,int rank,int maxdims,int *coords);
int MPI_Cart_get(MPI_Comm comm,int maxdims,int *dims,int *periods,int *coords);

int MPI_Send(void* buf, int count, MPI_Datatype datatype, int dest,int tag, MPI_Comm comm);
int MPI_Recv(void* buf, int count, MPI_Datatype datatype, int source,int tag, MPI_Comm comm, MPI_Status *status);
int MPI_Reduce(void *sb,void *rb,int count,MPI_Datatype dt,MPI_Op op,int root,MPI_Comm comm);
int MPI_Bcast(void *buf,int count,MPI_Datatype dt,int root,MPI_Comm comm);
int MPI_Gather(void *sb,int scount,MPI_Datatype sdt,void *rb,int rcount,MPI_Datatype rdt,int root,MPI_Comm comm);
int MPI_Scatter(void *sb,int scount,MPI_Datatype sdt,void *rb,int rcount,MPI_Datatype rdt,int root,MPI_Comm comm);
int MPI_Barrier(MPI_Comm comm);
