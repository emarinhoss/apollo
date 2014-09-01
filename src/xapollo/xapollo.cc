//includes
#include <apsimulation.h>
#include <wxexcept.h>
#include <wxmpimsg.h>
#include <wxlogger.h>
#include <wxlogstream.h>

// std includes
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>

template <typename REAL>
void
apolloMain(int argc, char **argv)
{
  // initialize message passing object
  // parallel simulation
  WxMpiMsg msg;

  // create new simulation
  ApSimulation<REAL> sim(argc, argv);

  // set its I/O and messaging objects
  //sim.setIo(io);
  sim.setMsg(msg);

  // open input file for reading
  std::string inpFileName = sim.getInpFileName();
  std::ifstream inp(inpFileName.c_str());
  if (!inp)
  { // we need to print to cerr as loggers have not been initialized
       std::cerr << "Input filename " << inpFileName << " not found.  Exiting." << std::endl;
       exit(1);
  }

  try
  {
    // create cryptset for complete simulation
    WxCryptSet inputSet(inp);

    //step 1: setup the class using its cryptset
    sim.setup(inputSet);

    //step 2: run init() or load().
//    if (sim.isRestarting())
//    { // restarting old simulation
//      std::string rf = sim.getRestartFile();
//      // open H5 file
//      WxIoNodeType h5f = io.openFile(rf, "r");
//      WxIoNodeType rootGrp = io.openGroup(h5f, "/");
//      sim.load(io, rootGrp);
//    }
//    else
//    { // starting new simulation
      sim.init();
//    }

    // simulation is now fully constructed: run it!
//    sim.simulate();
  }
  catch (const WxExcept&e )
  { // an error has occured: print message
    if (msg.rank() == 0)
    {
      std::cerr << "Exception caught...." << std::endl;
      std::cerr << e.what() << std::endl;
      exit(1);
    }
  }
}

// 
// Main entry point into system
//
int
main(int argc, char **argv)
{
  MPI_Init(&argc, &argv);

  // run top level main
  apolloMain<double>(argc, argv);

  MPI_Finalize();
}
