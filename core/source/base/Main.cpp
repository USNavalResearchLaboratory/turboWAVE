#include "mpi.h"
#include "tw_includes.h"
#include <thread>
import base;
import twmodule;
import logger;

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
	std::map<std::string,Arg> arg_map; // searchable
	std::vector<Arg> option_vec; // listable
	std::vector<std::string> cmd_line;
	std::vector<std::string>::iterator cmd_iter;
	public:
	CommandHandler(int argc,char *argv[]);
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
		int& omp_threads,
		std::string& ipath,
		std::string& rpath,
		std::string& platform,
		std::string& device,
		std::string& unitTest,
		bool& interactive,
		tw::Int& outputLevel);
};

CommandHandler::CommandHandler(int argc,char *argv[])
{
	for (int i=0;i<argc;i++)
		cmd_line.push_back(argv[i]);
	cmd_iter = cmd_line.begin();
	AddFlag("","--interactive","MPI rank 0 will spawn an interactive thread");
	AddOption("-c","","<threads>","number of OpenMP threads per MPI process");
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
	std::cout << "Usage: <launcher> -np <procs> tw3d [optional arguments...]" << std::endl;
	std::cout << "  <launcher> : MPI launcher such as mpirun, mpiexec, etc.." << std::endl;
	std::cout << "  -np <procs> : Launch <procs> MPI processes (flag may vary with launcher program)." << std::endl;
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
	std::cout << "Logging level is controlled by TW_LOG environment variable" << std::endl;
	std::cout << "Valid levels are error, warn, info, debug, trace." << std::endl;
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
	if (arg=="--interactive")
		interactive = true;
	if (arg=="--output-level")
	{
		outputLevel = std::stoi(arg.values[0]);
		if (outputLevel<0 || outputLevel>1)
			throw tw::FatalError("Illegal output level");
	}
}

/////////////////////////
// INTERACTIVE THREAD  //
/////////////////////////

void start_interactive(Simulation *tw)
{
	std::println(std::cout,"interactive thread is running (type help)");
	std::flush(std::cout);
	std::string cmd;
	auto running = [&] {
		return tw->StepNow() < tw->StepsToTake() && tw->WindowPos(0) < tw->MaxWindowPos(0);
	};
	do {
		std::getline(std::cin,cmd);
		if (running())
			tw->InteractiveCommand(cmd,&std::cout);
	} while (running());
}

/////////////////////////
// PROGRAM ENTRY POINT //
/////////////////////////

int main(int argc,char *argv[])
{
	std::thread ithread;
	tw::Int failure_count = 0;
	bool interactive = false;
	int numOMPThreads=1;
	tw::Int outputLevel = 0, errorCheckingLevel = 0;
	std::string initMessage,arg,inputFileName("stdin"),restartFileName("tw::none"),platform("cuda"),device("tesla"),unitTest("tw::none");

	MPI_Init(&argc,&argv);
	logger::init();

	try
	{
		Arg arg;
		auto cli = CommandHandler(argc,argv);
		while (cli.Next(arg))
			cli.Evaluate(arg,numOMPThreads,inputFileName,restartFileName,platform,device,unitTest,interactive,outputLevel);

	}
	catch (tw::FatalError& e)
	{
		std::cout << "COMMAND LINE ERROR: " << e.what() << std::endl;
		exit(1);
	}

	Simulation *tw = new Simulation(interactive,unitTest,inputFileName,restartFileName,platform,device,outputLevel,errorCheckingLevel);
	omp_set_num_threads(numOMPThreads);
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	if (interactive && rank==0) {
		ithread = std::thread(start_interactive,tw);
	}

	if (tw->unitTest=="tw::none") {
		tw->Run();
	} else {
		tw->Test();
		failure_count = tw->failure_count;
	}

	std::flush(std::cout);
	if (interactive && rank==0) {
		ithread.join();
	}
	
	delete tw; // since this might call MPI release functions, has to precede MPI_Finalize()

	MPI_Finalize();

	return failure_count;
}
