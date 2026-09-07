//includes
#include <apsimulation.h>
#include <wxexcept.h>
#include <wxmpimsg.h>
#include <wxlogger.h>
#include <wxlogstream.h>

// std includes
#include <mpi.h>
#include <cstdlib>
#include <stdexcept>
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
       if (msg.rank() == 0)
       {
         std::cerr << "Apollo: input file '" << inpFileName
                   << "' could not be opened." << std::endl;
         if (inpFileName.size() > 4 &&
             inpFileName.compare(inpFileName.size()-4, 4, ".pin") == 0)
           std::cerr << "  Decks under examples/ are .pin templates. Expand one "
                        "into the .inp the solver reads:\n"
                        "    python3 scripts/wxinpparse.py -i "
                     << inpFileName << std::endl;
       }
       // Every rank reads the same file, so every rank fails here together; even
       // so, abort through MPI rather than exit() so no rank is left in a
       // collective waiting for one that has already gone.
       MPI_Abort(MPI_COMM_WORLD, 1);
  }

  if (sim.isRestarting())
  { // The restart machinery was never finished: the load() branch below is
    // commented out, so -r is accepted and then ignored, and the run silently
    // starts from the initial condition instead of the checkpoint. For a job
    // being resumed that is worse than refusing.
    if (msg.rank() == 0)
      std::cerr << "Apollo: --restart is not implemented. The option is parsed "
                   "but no checkpoint is ever read, so the run would silently "
                   "start from the initial condition." << std::endl;
    MPI_Abort(MPI_COMM_WORLD, 1);
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
       sim.simulate();
  }
  catch (const WxExcept& e)
  { // an error has occurred: report it and take the whole job down.
    //
    // Reporting only from rank 0 loses the message whenever the failure is on
    // another rank - a mesh partition that will not load, a boundary condition
    // that only exists on part of the domain. Worse, the ranks that were not
    // rank 0 used to fall out of this catch and return normally, so they reached
    // PetscFinalize() and exited 0 while their peers were already gone: the job
    // either hung in a collective or reported success after failing.
    std::cerr << "Apollo: error on MPI rank " << msg.rank() << " of "
              << msg.numProcs() << ":\n  " << e.what() << std::endl;
    std::cerr.flush();
    MPI_Abort(MPI_COMM_WORLD, 1);
  }
  catch (const std::exception& e)
  {
    std::cerr << "Apollo: unexpected error on MPI rank " << msg.rank()
              << ": " << e.what() << std::endl;
    std::cerr.flush();
    MPI_Abort(MPI_COMM_WORLD, 1);
  }
}

// 
// Main entry point into system
//
int
main(int argc, char **argv)
{
    PetscInitialize(&argc, &argv, PETSC_NULL, PETSC_NULL); //MPI_Init(&argc, &argv);

    // run top level main
    apolloMain<double>(argc, argv);

    PetscFinalize(); //MPI_Finalize();

    // main() without a return statement yields 0, which is right here but only
    // by accident; say so.
    return 0;
}
