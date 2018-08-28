#include "sim.h"

////////////////////
// ERROR HANDLING //
////////////////////

void out_of_store();

void out_of_store()
{
	std::cout << "operator new failed: out of memory" << std::endl;
	exit(1);
}

//////////////////////////////
// Compute Thread Launcher  //
// and interactive thread   //
// (used with internal MPI) //
//////////////////////////////

#ifdef USING_TW_MPI

struct Launcher : tw::Thread
{
	Grid *tw;
	tw::Int numOMPThreads;
	Launcher(tw::Int rank,tw::Int c) : tw::Thread(rank)
	{
		numOMPThreads=c;
		tw = new Grid;
	}
	~Launcher()
	{
		delete tw; // delete grid must go before MPI_Finalize
		MPI_Finalize();
	}
	virtual void Run();
};

struct TW_Interactive : tw::Thread
{
	Grid *oneGrid;
	TW_Interactive(Grid *theGrid) : tw::Thread(0) { oneGrid = theGrid; }
	virtual void Run();
};


void Launcher::Run()
{
	#ifdef USE_OPENMP
	omp_set_num_threads(numOMPThreads);
	#endif
	MPI_Init(NULL,NULL);
	std::string messg = tw->InputFileFirstPass();
	*tw->tw_out << messg;
	tw->Run();
}

void TW_Interactive::Run()
{
	std::string cmd;
	do
	{
		std::getline(std::cin,cmd);
		if (!oneGrid->Completed())
			oneGrid->InteractiveCommand(cmd,&std::cout);
	} while (!oneGrid->Completed());
}

#endif

/////////////////////////
// PROGRAM ENTRY POINT //
// Standard MPI Launch //
/////////////////////////

#ifndef USING_TW_MPI
int main(int argc,char *argv[])
{
	std::string initMessage;
	std::set_new_handler(&out_of_store);
	tw::Int bitsPerFloat = sizeof(tw::Float)*8;

	MPI_Init(&argc,&argv);

	Grid *tw = new Grid;
	initMessage = tw->InputFileFirstPass();
	*tw->tw_out << std::endl << "*** Starting turboWAVE Session ***" << std::endl;
	*tw->tw_out << "Floating point precision = " << bitsPerFloat << " bits" << std::endl;
	#ifdef USE_OPENMP
	*tw->tw_out << "Maximum OpenMP threads = " << omp_get_max_threads() << std::endl;
	#endif
	*tw->tw_out << std::endl;
	*tw->tw_out << initMessage;

	tw->Run();

	delete tw; // since this might call MPI release functions, has to precede MPI_Finalize()

	*tw->tw_out << std::endl << "*** TurboWAVE Session Completed ***" << std::endl;

	MPI_Finalize();

	return 0;
}
#endif

/////////////////////////
// PROGRAM ENTRY POINT //
// Internal MPI Launch //
/////////////////////////

#ifdef USING_TW_MPI
int main(int argc,char *argv[])
{
	int numMPIThreads=1,numOMPThreads=1;
	bool interactive = true;
	tw::Int i,numCompleted=0;
	std::string arg;
	std::set_new_handler(&out_of_store);
	tw::Int bitsPerFloat = sizeof(tw::Float)*8;

	std::cout << std::endl << "*** Starting turboWAVE Session ***" << std::endl;
	std::cout << "Floating point precision = " << bitsPerFloat << " bits" << std::endl;

	try
	{
		#ifdef USE_OPENMP
		if (argc==1)
		{
			numOMPThreads = omp_get_max_threads();
			std::cout << "INFO: Defaulting to single MPI task with OpenMP parameters based on shell environment." << std::endl;
			std::cout << "INFO: Better performance may be obtained with hybrid MPI/OpenMP." << std::endl;
		}
		#endif

		tw::Int idx = 1;
		while (idx<argc)
		{
			arg = std::string(argv[idx]);

			if (arg!="-n" && arg!="-c" && arg!="-nointeractive" && arg!="--version" && arg!="--help")
				throw tw::FatalError("Unrecognized argument");

			if (arg=="--version")
			{
				std::cout << "turboWAVE version 3.1b" << std::endl;
				if (argc==2)
					exit(0);
			}

			if (arg=="--help")
			{
				std::cout << "This is turboWAVE, a PIC/hydro/quantum simulation code." << std::endl;
				std::cout << "Please see the documentation at http://readthedocs.org." << std::endl;
				if (argc==2)
					exit(0);
			}

			if (arg=="-n")
			{
				idx++;
				if (idx<argc)
					numMPIThreads = atoi(argv[idx]);
				else
					throw tw::FatalError("Incomplete arguments");
			}

			if (arg=="-c")
			{
				idx++;
				if (idx<argc)
					numOMPThreads = atoi(argv[idx]);
				else
					throw tw::FatalError("Incomplete arguments");
			}

			if (arg=="-nointeractive")
				interactive = false;

			idx++;
		}
		if (numMPIThreads<1)
			throw tw::FatalError("Number of MPI threads < 1");
		if (numOMPThreads<1)
			throw tw::FatalError("Number of OpenMP threads < 1");
		#ifdef USE_OPENMP
		if (numOMPThreads==1)
			std::cout << "INFO: Use '-c' option to set OpenMP threads per MPI task." << std::endl;
		#endif
		#ifndef USE_OPENMP
		if (numOMPThreads>1)
			throw tw::FatalError("Requested OpenMP threads but this executable not OpenMP enabled");
		#endif
	}
	catch (tw::FatalError& e)
	{
		std::cout << "COMMAND LINE ERROR: " << e.what() << std::endl;
		exit(1);
	}

	#ifdef USE_OPENMP
	// note: cannot set OpenMP threads here due to scope (must be in MPI thread)
	std::cout << "Requested OpenMP threads = " << numOMPThreads << std::endl;
	#endif
	std::cout << std::endl;

	TW_MPI_Lock();

	std::valarray<tw::Thread*> launcher(numMPIThreads);
	for (i=0;i<numMPIThreads;i++)
		launcher[i] = new Launcher(i,numOMPThreads);
	TW_MPI_Launch(numMPIThreads,&launcher[0]);

	std::cout << "Internal MPI Startup Complete" << std::endl;

	std::cout << "Launch Interactive Thread..." << std::endl;

	TW_Interactive interactiveThread(((Launcher*)(launcher[0]))->tw);
	if (interactive)
		interactiveThread.Start();

	TW_MPI_Unlock();

	for (i=0;i<numMPIThreads;i++)
		launcher[i]->Complete();

	std::cout << std::endl << "*** TurboWAVE Session Completed ***" << std::endl;

	if (interactive)
	{
		std::cout << std::endl << "Press Enter to terminate interactive thread." << std::endl;
		interactiveThread.Complete();
	}

	for (i=0;i<numMPIThreads;i++)
		delete launcher[i];

	return 0;
}
#endif
