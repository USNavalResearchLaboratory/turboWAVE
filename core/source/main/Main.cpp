#include "mpi.h"
#include "tw_includes.h"
#ifndef USE_STD_MODULE
#include <thread>
#endif
#ifdef _WIN64
#include <Windows.h>
#endif
import base;
import driver;
import tasks;
import metric_space;
import simulation;
import commands;
import logger;

/////////////////////////
// INTERACTIVE THREAD  //
/////////////////////////

void start_interactive(Simulation *tw)
{
	std::println(std::cout,"interactive thread is running (type help)");
	std::flush(std::cout);
	std::string cmd;
	auto running = [&] {
		return tw->space->StepNow() < tw->space->StepsToTake() && tw->space->WindowPos(0) < tw->space->MaxWindowPos(0);
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
	#ifdef _WIN64
	SetConsoleCP(CP_UTF8);
	SetConsoleOutputCP(CP_UTF8);
	#endif
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

	Task *universalTask = new Task();
	MetricSpace *universalSpace = new MetricSpace();
	Simulation *tw = new Simulation(interactive,unitTest,inputFileName,restartFileName,
		platform,device,outputLevel,errorCheckingLevel,
		universalSpace,
		universalTask
	);
	omp_set_num_threads(numOMPThreads);
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	if (interactive && rank==0) {
		ithread = std::thread(start_interactive,tw);
	}

	try {
		if (universalTask->unitTest=="tw::none") {
			tw->Run();
		} else {
			tw->Test();
			failure_count = tw->failure_count;
		}
	} catch (tw::FatalError& e) {
		std::println("TOP LEVEL ERROR {}",e.what());
	}

	std::flush(std::cout);
	if (interactive && rank==0) {
		ithread.join();
	}
	
	// release anything that might call MPI release functions, must precede MPI_Finalize()
	delete tw;
	delete universalTask;
	delete universalSpace;

	MPI_Finalize();

	return failure_count;
}
