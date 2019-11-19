#include "simulation.h"
#define TW_VERSION_STRING "3.4.1"

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
	Simulation *tw;
	tw::Int numOMPThreads;
	Launcher(tw::Int rank,tw::Int c,const std::string& fileName) : tw::Thread(rank)
	{
		numOMPThreads=c;
		tw = new Simulation(fileName);
	}
	virtual ~Launcher()
	{
		delete tw; // delete grid must go before MPI_Finalize
		MPI_Finalize();
	}
	virtual void Run();
};

struct TW_Interactive : tw::Thread
{
	Simulation *oneSim;
	TW_Interactive(Simulation *sim) : tw::Thread(0) { oneSim = sim; }
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
		if (!oneSim->Completed())
			oneSim->InteractiveCommand(cmd,&std::cout);
	} while (!oneSim->Completed());
}

#endif

/////////////////////////
// PROGRAM ENTRY POINT //
// Standard MPI Launch //
/////////////////////////

#ifndef USING_TW_MPI
int main(int argc,char *argv[])
{
	int numOMPThreads=0; // indicates -c argument was not given
	std::string initMessage,arg,inputFileName("stdin");
	std::set_new_handler(&out_of_store);
	tw::Int bitsPerFloat = sizeof(tw::Float)*8;

	MPI_Init(&argc,&argv);

	try
	{
		tw::Int idx = 1;
		while (idx<argc)
		{
			arg = std::string(argv[idx]);

			if (arg!="--input-file" && arg!="--no-interactive" && arg!="--version" && arg!="--help" && arg!="-c")
				throw tw::FatalError("Unrecognized argument");

			if (arg=="--version")
			{
				std::cout << "turboWAVE version " << TW_VERSION_STRING << std::endl;
				if (argc==2)
					exit(0);
			}

			if (arg=="--help")
			{
				std::cout << "This is turboWAVE, a PIC/hydro/quantum simulation code." << std::endl;
				std::cout << "Usage: <launcher> -np <procs> tw3d [-c <threads>] [--input-file <path>] [--version] [--help] [--no-interactive]" << std::endl;
				std::cout << "<launcher>             MPI launcher such as mpirun, mpiexec, etc.." << std::endl;
				std::cout << "-np <procs>            Launch <procs> MPI processes (flag may vary with launcher program)." << std::endl;
				std::cout << "-c <threads>           Fork <threads> OpenMP threads per MPI process." << std::endl;
				std::cout << "--input-file <path>    Use <path> as the input file." << std::endl;
				std::cout << "--version              Display the version number (no simulation if only argument)." << std::endl;
				std::cout << "--help                 Display this message (no simulation if only argument)." << std::endl;
				std::cout << "--no-interactive       Suppress the interactive thread." << std::endl;
				std::cout << "Full documentation can be found at https://turbowave.readthedocs.io" << std::endl;
				if (argc==2)
					exit(0);
			}
			if (arg=="--input-file")
			{
				idx++;
				if (idx<argc)
					inputFileName = std::string(argv[idx]);
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
				if (numOMPThreads<1)
					throw tw::FatalError("Number of OpenMP threads < 1");
			}
			idx++;
		}
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

	Simulation *tw = new Simulation(inputFileName);
	initMessage = tw->InputFileFirstPass();
	*tw->tw_out << std::endl << "*** Starting turboWAVE Session ***" << std::endl;
	*tw->tw_out << "Floating point precision = " << bitsPerFloat << " bits" << std::endl;
	#ifdef USE_OPENMP
	if (numOMPThreads>0)
		omp_set_num_threads(numOMPThreads);
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
	std::string arg,inputFileName("stdin");
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

			if (arg!="-n" && arg!="-c" && arg!="--no-interactive" && arg!="--version" && arg!="--help" && arg!="--input-file")
				throw tw::FatalError("Unrecognized argument");

			if (arg=="--input-file")
			{
				idx++;
				if (idx<argc)
					inputFileName = std::string(argv[idx]);
				else
					throw tw::FatalError("Incomplete arguments");
			}

			if (arg=="--version")
			{
				std::cout << "turboWAVE version " << TW_VERSION_STRING << std::endl;
				if (argc==2)
					exit(0);
			}

			if (arg=="--help")
			{
				std::cout << "This is turboWAVE, a PIC/hydro/quantum simulation code." << std::endl;
				std::cout << "Usage: tw3d [-n <threads>] [-c <threads>] [--input-file <path>] [--version] [--help] [--no-interactive]" << std::endl;
				std::cout << "-n <threads>           Start <threads> MPI threads." << std::endl;
				std::cout << "-c <threads>           Fork <threads> OpenMP threads per MPI thread." << std::endl;
				std::cout << "--input-file <path>    Use <path> as the input file." << std::endl;
				std::cout << "--version              Display the version number (no simulation if only argument)." << std::endl;
				std::cout << "--help                 Display this message (no simulation if only argument)." << std::endl;
				std::cout << "--no-interactive       Suppress the interactive thread." << std::endl;
				std::cout << "Full documentation can be found at https://turbowave.readthedocs.io" << std::endl;
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

			if (arg=="--no-interactive")
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

	std::vector<tw::Thread*> launcher(numMPIThreads);
	for (i=0;i<numMPIThreads;i++)
		launcher[i] = new Launcher(i,numOMPThreads,inputFileName);
	TW_MPI_Launch(launcher);

	std::cout << "Internal MPI Startup Complete" << std::endl;

	TW_Interactive interactiveThread(((Launcher*)(launcher[0]))->tw);
	if (interactive)
	{
		std::cout << "Launch Interactive Thread..." << std::endl;
		interactiveThread.Start();
	}

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
