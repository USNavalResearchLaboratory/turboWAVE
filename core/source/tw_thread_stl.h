// tw::Thread is a wrapper for std::thread with additional structure
// for a message passing paradigm, i.e., tw::queue.  TW_MPI is built on this.
namespace tw
{
	struct packet
	{
		void *data,*sender;
	};

	struct mutex
	{
		std::mutex theMutex;

		void Lock()
		{
			theMutex.lock();
		}
		void Unlock()
		{
			theMutex.unlock();
		}
	};

	struct queue
	{
		std::vector<tw::packet> theQueue;
		std::mutex protectQueue;
		std::condition_variable_any theCondition;

		void Notify(void* data,void* sender)
		{
			protectQueue.lock();
			tw::packet thePacket;
			thePacket.data = data;
			thePacket.sender = sender;
			theQueue.push_back(thePacket);
			theCondition.notify_all();
			protectQueue.unlock();
		}

		tw::Int Wait(void** data,void* sender,tw::Int timeout)
		{
			// wait for data from sender
			// if sender==NULL then accept data from anyone
			// timeout is in secs; if -1 wait forever

			tw::Int i;
			tw::Int pos = -1;

			protectQueue.lock();

			for (i=0;i<theQueue.size();i++)
				if (sender==theQueue[i].sender || sender==NULL)
				{
					pos = i;
					if (data!=NULL)
						*data = theQueue[i].data;
				}

			while (pos==-1)
			{
				if (timeout==-1)
				{
					theCondition.wait(protectQueue);
				}
				else
				{
				    if (theCondition.wait_for(protectQueue,std::chrono::seconds(timeout))==std::cv_status::timeout)
						return 1; // timeout
				}

				for (i=0;i<theQueue.size();i++)
					if (sender==theQueue[i].sender || sender==NULL)
					{
						pos = i;
						if (data!=NULL)
							*data = theQueue[i].data;
					}
			}

			theQueue.erase(theQueue.begin() + pos);

			protectQueue.unlock();

			return 0;
		}

		tw::Int WaitAny(void** data,void** sender,tw::Int timeout)
		{
			// wait for data from any sender

			tw::Int pos = -1;

			protectQueue.lock();

			if (theQueue.size()>0)
			{
				pos = 0;
				if (data!=NULL)
					*data = theQueue[0].data;
				if (sender!=NULL)
					*sender = theQueue[0].sender;
			}

			while (pos==-1)
			{
				if (timeout==-1)
				{
					theCondition.wait(protectQueue);
				}
				else
				{
				    if (theCondition.wait_for(protectQueue,std::chrono::seconds(timeout))==std::cv_status::timeout)
						return 1; // timeout
				}

				if (theQueue.size()>0)
				{
					pos = 0;
					if (data!=NULL)
						*data = theQueue[0].data;
					if (sender!=NULL)
						*sender = theQueue[0].sender;
				}
			}

			theQueue.erase(theQueue.begin() + pos);

			protectQueue.unlock();

			return 0;
		}
	};



	//////////////
	//  THREAD  //
	//////////////


	struct Thread
	{
		private:

		tw::Int refNum;
		std::thread sysThread;

		public:

		tw::queue input,output,self;

		Thread(tw::Int n)
		{
			refNum = n;
		}
		bool IsCallingThread()
		{
			return sysThread.get_id()==std::this_thread::get_id();
		}
		tw::Int GetRefNum()
		{
			return refNum;
		}
		void Dispatch()
		{
			Run();
		}
		void Start()
		{
			sysThread = std::thread(&Thread::Dispatch,this);
		}
		virtual void Run()
		{
		}
		virtual void Complete()
		{
			// must be called from parent thread, not this thread
			sysThread.join();
		}
		void SetAffinity(tw::Int cpuid)
		{
			// pthreads affinity; so far not found useful.
			// must rewrite for std::thread at any rate.
	// 		int err;
	// 		cpu_set_t cpuset;
	// 		CPU_ZERO(&cpuset);
	// 		CPU_SET(cpuid,&cpuset);
	// 		err = pthread_setaffinity_np(sysThread,sizeof(cpu_set_t),&cpuset);
	// 		if (err)
	// 			throw tw::FatalError("Could not bind thread.  Repair or delete affinity setting.");
		}
	};
}
