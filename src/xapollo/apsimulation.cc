// includes
#include <apsimulation.h>

// include getopt or mygetopt
#ifdef _NO_GETOPT_
# include <mygetopt.h>
#else
# include <getopt.h>
#endif

// lib includes
#include <wxobject.h>
#include <wxlogger.h>
#include <wxlogstream.h>

// std includes
#include <ctime>
#include <iostream>

template <typename REAL>
ApSimulation<REAL>::ApSimulation(int argc, char **argv)
  : solver(0), inpFileName("apollo.inp"),
    outPrefixSet(false), restartSim(false), realType("double")
{
  // parse command line parameters
  parseCmdLine(argc, argv);
//  PetscInitialize(&argc, &argv, PETSC_NULL, PETSC_NULL);
}

template <typename REAL>
ApSimulation<REAL>::~ApSimulation()
{
  // delete top level solver object
  delete solver;
//  PetscFinalize();
  // delete message object
  // WxLogger::cleanUp();
}

template <typename REAL>
std::string
ApSimulation<REAL>::getInpFileName() const
{
  return inpFileName;
}

template <typename REAL>
std::string
ApSimulation<REAL>::getRunName() const
{
  return runName;
}

template <typename REAL>
bool
ApSimulation<REAL>::isRestarting() const
{
  return restartSim;
}

template <typename REAL>
std::string
ApSimulation<REAL>::getRestartFile() const
{
  return restartFile;
}

template <typename REAL>
void
ApSimulation<REAL>::init()
{
  // initialize top level solver
  solver->init();
}

//template <typename REAL>
//void
//ApSimulation<REAL>::load(WxIoBase& io, const WxIoNodeType& grpNode)
//{
//  // open timeData group
//  WxIoNodeType timeGrp = io.openGroup(grpNode, "timeData");
//  // get time at which we are restarting
//  REAL tstart;
//  io.readAttribute<REAL>(timeGrp, "time", tstart);
//  solver->setCurrentTime(tstart);
//  // get frame from which we are starting
//  unsigned frame;
//  io.readAttribute<unsigned>(timeGrp, "step", frame);
//  solver->setStartFrame(frame);

//  // open solver group
//  WxIoNodeType solverGrp = io.openGroup(grpNode, solver->getSolverName());
//  // load solver into memory
//  solver->load(io, solverGrp);
//}

template <typename REAL>
void
ApSimulation<REAL>::setup(const WxCryptSet& wxc)
{
    // determine run name
    if (outPrefixSet)
        runName = outPrefix;
    else
        runName = stripName(inpFileName);

    // create Apollo root logger
    WxLogger *wr = WxLogger::get("apollo-root");

    // set root logger's verbosity level
    std::string level;
    if (wxc.has("GlobalVerbosity"))
      level = wxc.get<std::string>("GlobalVerbosity");
    else
      level = "debug";
    wr->setLevel(level);

    // add file handler to root logger
    std::ostringstream fn;
    fn << runName << "_" << this->getMsg().rank() << ".log";
    WxLogRecordHandler *wrfhndl = new WxFileHandler(
        fn.str());
    wr->addHandler(wrfhndl);

    // create Apollo console logger
    WxLogger *wrc = WxLogger::get("apollo-root.console");

    // set console logger's verbosity level
    if (wxc.has("Verbosity"))
      level = wxc.get<std::string>("Verbosity");
    else
      level = "debug";
    wrc->setLevel(level);

    // add stream handler to console
    WxLogRecordHandler *wrshndl = new WxStreamHandler();
    if (this->getMsg().rank() == 0)
      // add console stream only on rank 0
      wrc->addHandler(wrshndl);

    // get output streams from newly created loggers
    WxLogStream infStrm = wrc->getInfoStream();
//    WxLogStream errStrm = wrc->getErrorStream();
//    WxLogStream wrnStrm = wrc->getWarningStream();

    // now setup top level solver
    infStrm << "Setting up Apollo simulation..." << std::endl;
    std::string simName;
    // name of simulation to run
    if (wxc.has("Simulation"))
      simName = wxc.get<std::string>("Simulation");
    else
    { // no simulation specified, so throw exception
      WxExcept wxe("No Simulation key found in ");
      wxe << inpFileName << std::endl;
      throw wxe;
    }

    if (!wxc.hasSet(simName))
    { // solver not found
      WxExcept wxe("ERROR: Solver set ");
      wxe << simName << " not found" << std::endl;
      throw wxe;
    }
    // get hold of solver's cryptset
    const WxCryptSet& solverCrypt = wxc.getSet(simName);
    infStrm << "Simulation name is " << simName << std::endl;

    // create new solver
    solver = new ApSolver<REAL>(simName);

    // set run name
    solver->setRunName(runName);
    // set name of solver
    solver->setSolverName(simName);

    // set I/O for use in solver
    //solver->setIo(this->getIo());
    // set msg for use in solver
    //solver->setMsg(this->getMsg());
    // setup solver
    solver->setup(solverCrypt);
}

template <typename REAL>
void
ApSimulation<REAL>::simulate()
{
    WxLogStream debStrm = WxLogger::get("apollo-root.console")->getDebugStream();
    WxLogStream infStrm = WxLogger::get("apollo-root.console")->getInfoStream();
    debStrm << "Running simulation...\n";

    time_t start = time(0); // time at start of main loop
    struct tm * timeinfo;
    timeinfo = localtime ( &start );
    infStrm << "Simulation started at time " << asctime(timeinfo) << std::endl;

    // run simulation
    solver->solve();
    time_t end = time(0); // time at end of main loop
    timeinfo = localtime ( &end );
    infStrm << "Simulation finished at time " << asctime(timeinfo) << std::endl;
}

template <typename REAL>
void
ApSimulation<REAL>::parseCmdLine(int argc, char **argv)
{
  // options description for use in cmd line parser
  struct option longopts[] = {
    { "input-file", required_argument, NULL, 'i'},
    { "output-prefix", required_argument, NULL, 'o'},
    { "restart", required_argument, NULL, 'r'},
    { "real-type", required_argument, NULL, 0},
    { "help", no_argument, NULL, 1},
    { NULL, 0, NULL, 0}
  };

  // getopt_long returns int, and signals "no more options" with -1. Storing that
  // in a char is only correct where char happens to be signed: on ARM and
  // PowerPC, where it is unsigned, -1 becomes 255 and the loop never ends.
  int ch;
  // parse command line parameters. 'h' is in the option string because usage()
  // advertises -h; without it getopt rejected the documented flag.
  while ((ch = getopt_long(argc, argv, "hi:r:o:", longopts, NULL)) != -1)
  {
    switch (ch)
    {
      case 'i':
        // input file name
        inpFileName = optarg;
        break;

      case 'o':
        // output prefix
        outPrefix = optarg;
        outPrefixSet = true;
        break;

      case 'r':
        // restart simulation
        restartFile = optarg; // name of class to describe
        restartSim = true;
        break;

      case 0:
        // real number type to use
        realType = optarg;
        break;

      case 'h':
      case 1:
        // help
        usage();
        exit(0);

      default:
        // getopt_long has already written its own diagnostic; do not continue
        // with a half-parsed command line and a default input file.
        usage();
        exit(2);
    }
  }
  // Consume the options; whatever is left is a positional argument.
  argc -= optind;
  argv += optind;

  if (argc > 0)
  { // Apollo takes its input file with -i. A bare filename is not the input
    // file, and silently falling back to apollo.inp made that look like a
    // missing-file error rather than the usage error it is.
    std::cerr << "Apollo: unexpected argument '" << argv[0]
              << "'.\n  The input file is given with -i, for example:\n"
                 "    apollo -i " << argv[0] << std::endl;
    exit(2);
  }

  if (realType == "float")
    std::cerr << "Apollo: warning: --real-type=float is accepted but has no "
                 "effect; the solver is instantiated as double." << std::endl;
}

template <typename REAL>
void
ApSimulation<REAL>::usage()
{
  // print help
  std::cout << "*** Welcome to Apollo ***" << std::endl;
  std::cout << "Version XXX from svn revision YYY" << std::endl;
  std::cout << "Apollo accepts the following command line options" << std::endl;

  // for help message
  std::cout << " -h\n"
            << " --help\n"
            << "    Print this help message.\n" << std::endl;
  // for input file
  std::cout << " -i <file-name>\n"
            << " --input-file=<file-name>\n"
            << "    Read input from <file-name>. Defaults to apollo.inp.\n" << std::endl;
  // for output prefix
  std::cout << " -o <output-prefix>\n"
            << " --output-prefix=<output-prefix>\n"
            << "    Use this as output prefix. Defaults to input file name without extension.\n" << std::endl;
  // for restart file
  std::cout << "-r <file-name>\n" 
            << " --restart=<file-name>\n"
            << "    Restart simulation from file <file-name>.\n" << std::endl;
  // for specifing real number type to use
  std::cout << " --real-type=[float|double]\n"
            << "    Real number type to use. Defaults to double.\n" << std::endl;

  std::cout << std::endl;
}

template <typename REAL>
std::string
ApSimulation<REAL>::stripName(const std::string& nm)
{
  std::string snm = nm;
  // find_last_of returns npos when there is no '.', and npos is not 0, so the
  // old `if (trunc > 0)` test passed and erase(npos, ...) threw std::out_of_range
  // for any input file whose name has no extension.
  std::string::size_type trunc = snm.find_last_of('.');
  if (trunc != std::string::npos && trunc > 0)
    snm.erase(trunc);
  return snm;
}

// instantiations
//template class ApSimulation<float>;
template class ApSimulation<double>;
