#include "simulation.h"

////////////////////
// ERROR HANDLING //
////////////////////

void out_of_store();

void out_of_store()
{
	std::cout << "operator new failed: out of memory" << std::endl;
	exit(1);
}

////////////////////////////
// COMMAND LINE INTERFACE //
////////////////////////////

/// Representation of a command line argument
///
/// All arguments are in the form of a flag followed by values.
/// Positional arguments are not allowed.
struct Arg
{
	std::string short_flag; ///< short form of the flag, can be empty if long form is not
	std::string long_flag; ///< long form of the flag, can be empty if short form is not
	std::string value_label; ///< label for the value(s) following the flag, can be empty
	std::string help; ///< help string for the argument
	std::vector<std::string> values; /// values following the flag, can be empty
	int num; /// 0: no values, 1: single value, 2: multiple values
	Arg() { num = 0; }
	Arg(const std::string& sflag,const std::string& lflag,const std::string& lab,const std::string& h,int n)
	{
		short_flag = sflag;
		long_flag = lflag;
		value_label = lab;
		help = h;
		num = n;
	}
	void AddValue(const std::string& val)
	{
		values.push_back(val);
	}
	bool operator == (const std::string& s) const
	{
		return s==short_flag || s==long_flag;
	}
};

/// Command line argument handler
///
/// The handler is responsible for parsing the command line, and
/// printing the help message.
/// Particulars are defined in the constructor and `Evaluate`.
class CommandHandler
{
	bool tw_mpi;
	std::map<std::string,Arg> arg_map; // searchable
	std::vector<Arg> option_vec; // listable
	std::vector<std::string> cmd_line;
	std::vector<std::string>::iterator cmd_iter;
	public:
	CommandHandler(int argc,char *argv[],bool tw_mpi);
	void AddSomething(const std::string& sflag,const std::string& lflag,const std::string& lab,const std::string& help,int n)
	{
		auto arg = Arg(sflag,lflag,lab,help,n);
		if (arg_map.find(sflag)!=arg_map.end())
			throw tw::FatalError("attempt to redefine command line key");
		if (sflag!="")
			arg_map[sflag] = arg;
		if (arg_map.find(lflag)!=arg_map.end())
			throw tw::FatalError("attempt to redefine command line key");
		if (lflag!="")
			arg_map[lflag] = arg;
		option_vec.push_back(arg);
	}
	void AddFlag(const std::string& sflag,const std::string& lflag,const std::string& help)
	{
		AddSomething(sflag,lflag,"",help,0);
	}
	void AddOption(const std::string& sflag,const std::string& lflag,const std::string& lab,const std::string& help)
	{
		AddSomething(sflag,lflag,lab,help,1);
	}
	void AddOptionMultiValued(const std::string& sflag,const std::string& lflag,const std::string& lab,const std::string& help)
	{
		AddSomething(sflag,lflag,lab,help,2);
	}
	void Usage();
	bool Next(Arg& arg);
	void Evaluate(const Arg& arg,
		int& mpi_procs,
		int& omp_threads,
		std::string& ipath,
		std::string& rpath,
		std::string& platform,
		std::string& device,
		std::string& unitTest,
		bool& interactive,
		tw::Int& outputLevel);
};

CommandHandler::CommandHandler(int argc,char *argv[],bool tw_mpi)
{
	this->tw_mpi = tw_mpi;
	for (int i=0;i<argc;i++)
		cmd_line.push_back(argv[i]);
	cmd_iter = cmd_line.begin();
	if (tw_mpi)
	{
		AddOption("-n","","<threads>","number of MPI threads to start");
		AddFlag("","--no-interactive","suppress the interactive thread");
	}
	AddOption("-c","","<threads>","number of OpenMP threads per MPI thread");
	AddOption("-i","--input-file","<path>","path to the input file");
	AddFlag("-v","--version","display the version number (can be only argument)");
	AddFlag("-h","--help","display this message (can be only argument)");
	AddFlag("","--restart","load checkpoint data");
	AddOption("","--platform","<plat_str>","OpenCL platform search string");
	AddOption("","--device","<dev_str>","OpenCL device search string");
	AddOption("","--output-level","<level>","0,1,...");
	AddOptionMultiValued("","--unit-test","<test>","name of the test");
}

void CommandHandler::Usage()
{
	std::cout << "turboWAVE, https://turbowave.readthedocs.io" << std::endl;
	if (tw_mpi)
	{
		std::cout << "Usage: tw3d [optional arguments...]" << std::endl;
	}
	else
	{
		std::cout << "Usage: <launcher> -np <procs> tw3d [optional arguments...]" << std::endl;
		std::cout << "  <launcher> : MPI launcher such as mpirun, mpiexec, etc.." << std::endl;
		std::cout << "  -np <procs> : Launch <procs> MPI processes (flag may vary with launcher program)." << std::endl;
	}
	std::string line;
	size_t tab1=0,tab2=0;
	// set tab stops
	for (auto it=option_vec.begin();it!=option_vec.end();++it)
	{
		line = it->short_flag + ", " + it->long_flag;
		if (line.size()>tab1)
			tab1 = line.size();
		line = it->value_label;
		if (line.size()>tab2)
			tab2 = line.size();
	}
	tab1 += 2;
	tab2 += 4;
	// print
	for (auto it=option_vec.begin();it!=option_vec.end();++it)
	{
		line = it->short_flag;
		if (it->short_flag!="" && it->long_flag!="")
			line += ", ";
		line += it->long_flag;
		line.append(tab1-line.size(),' ');
		line += it->value_label;
		line.append(tab1+tab2-line.size(),' ');
		line += it->help;
		std::cout << line << std::endl;
	}
}

bool CommandHandler::Next(Arg& arg)
{
	++cmd_iter;
	if (cmd_iter!=cmd_line.end())
	{
		auto it = arg_map.find(*cmd_iter);
		if (it!=arg_map.end())
		{
			arg = it->second;
			if (it->second.num)
			{
				do {
					++cmd_iter;
					if (cmd_iter!=cmd_line.end() && arg_map.find(*cmd_iter)==arg_map.end())
						arg.AddValue(*cmd_iter);
					else
					{
						--cmd_iter;
						break;
					}
				} while (true);
				if (it->second.num==1 && arg.values.size()!=1)
					throw tw::FatalError("expected one value for " + it->first);
				if (it->second.num>1 && arg.values.size()<1)
					throw tw::FatalError("expected at least one value for " + it->first);
			}
			return true;
		}
		throw tw::FatalError("unrecognized argument <" + *cmd_iter + ">");
	}
	return false;
}

void CommandHandler::Evaluate(const Arg& arg,
	int& mpi_procs,
	int& omp_threads,
	std::string& ipath,
	std::string& rpath,
	std::string& platform,
	std::string& device,
	std::string& unitTest,
	bool& interactive,
	tw::Int& outputLevel)
{
	if (arg=="--version")
	{
		std::cout << "turboWAVE version " << TW_VERSION_STRING << std::endl;
		if (cmd_line.size()==2)
			exit(0);
	}
	if (arg=="--help")
	{
		Usage();
		if (cmd_line.size()==2)
			exit(0);
	}
	if (arg=="-n")
	{
		mpi_procs = std::stoi(arg.values[0]);
		if (mpi_procs<1)
			throw tw::FatalError("Number of MPI threads < 1");
	}
	if (arg=="-c")
	{
		omp_threads = std::stoi(arg.values[0]);
		if (omp_threads<1)
			throw tw::FatalError("Number of OpenMP threads < 1");
	}
	if (arg=="--input-file")
		ipath = arg.values[0];
	if (arg=="--restart")
		rpath = "dump";
	if (arg=="--platform")
		platform = arg.values[0];
	if (arg=="--device")
		device = arg.values[0];
	if (arg=="--unit-test")
	{
		unitTest = "";
		for (auto val : arg.values)
			unitTest += val + " ";
		unitTest.pop_back();
	}
	if (arg=="--no-interactive")
		interactive = false;
	if (arg=="--output-level")
	{
		outputLevel = std::stoi(arg.values[0]);
		if (outputLevel<0 || outputLevel>1)
			throw tw::FatalError("Illegal output level");
	}
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
	Launcher(tw::Int rank,tw::Int c,
		const std::string& unitTest,
		const std::string& inputFileName,
		const std::string& restartFileName,
		const std::string& platform,
		const std::string& device,
		const tw::Int& outputLevel,
		const tw::Int& errorCheckingLevel) : tw::Thread(rank)
	{
		numOMPThreads=c;
		tw = new Simulation(unitTest,inputFileName,restartFileName,platform,device,outputLevel,errorCheckingLevel);
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
	// OpenMP threads must be set up within MPI thread
	omp_set_num_threads(numOMPThreads);
	#endif
	MPI_Init(NULL,NULL);
	if (tw->unitTest=="tw::none")
		tw->Run();
	else
		tw->Test();
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
	tw::Int outputLevel = 0, errorCheckingLevel = 0;
	std::string initMessage,arg,inputFileName("stdin"),restartFileName("tw::none"),platform("cuda"),device("tesla"),unitTest("tw::none");
	std::set_new_handler(&out_of_store);

	MPI_Init(&argc,&argv);

	try
	{
		int dummyProcs;
		bool dummyInteractive;
		Arg arg;
		auto cli = CommandHandler(argc,argv,false);
		while (cli.Next(arg))
			cli.Evaluate(arg,dummyProcs,numOMPThreads,inputFileName,restartFileName,platform,device,unitTest,dummyInteractive,outputLevel);

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

	Simulation *tw = new Simulation(unitTest,inputFileName,restartFileName,platform,device,outputLevel,errorCheckingLevel);
	#ifdef USE_OPENMP
	if (numOMPThreads>0) // if not given assume environment determines
		omp_set_num_threads(numOMPThreads);
	#endif

	if (tw->unitTest=="tw::none")
		tw->Run();
	else
		tw->Test();

	delete tw; // since this might call MPI release functions, has to precede MPI_Finalize()

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
	tw::Int i,outputLevel=0,errorCheckingLevel=0,numCompleted=0;
	std::string arg,inputFileName("stdin"),restartFileName("tw::none"),platform("cuda"),device("tesla"),unitTest("tw::none");
	std::set_new_handler(&out_of_store);
	tw::Int bitsPerFloat = sizeof(tw::Float)*8;

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

		Arg arg;
		auto cli = CommandHandler(argc,argv,true);
		while (cli.Next(arg))
			cli.Evaluate(arg,numMPIThreads,numOMPThreads,inputFileName,restartFileName,platform,device,unitTest,interactive,outputLevel);
		if (unitTest!="tw::none")
			interactive = false;

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

	TW_MPI_Lock();

	std::vector<tw::Thread*> launcher(numMPIThreads);
	for (i=0;i<numMPIThreads;i++)
		launcher[i] = new Launcher(i,numOMPThreads,unitTest,inputFileName,restartFileName,platform,device,outputLevel,errorCheckingLevel);
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
