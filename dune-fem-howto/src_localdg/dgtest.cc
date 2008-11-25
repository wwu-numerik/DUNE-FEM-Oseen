/**
 *	\file dgtest.cc
 *	\brief dgtest.cc
 **/

// Dune includes
#include <config.h>

#include <dune/common/fvector.hh>
#include <dune/fem/misc/utility.hh>
#include <dune/fem/gridpart/gridpart.hh>

#include <dune/common/misc.hh>
#include <dune/grid/common/quadraturerules.hh>

#include <dune/fem/misc/l2error.hh>

#include <iostream>
#include <string>

#include <dune/fem/io/parameter.hh>
#include <dune/fem/io/file/datawriter.hh>
#include <dune/common/timer.hh>
#include <dune/common/mpihelper.hh>

#include "models.hh"
// Include helper function for EOC LaTeX output and L2Projection
#include "stuff.cc"

using namespace Dune;

/**
 * @brief main function for the LocalDG Advection-Diffusion application
 *
 * \c main starts the Simulation of an advection-diffusion pde with
 * the localdg method with EOC analysis and visual output to grape, paraview or
 * gnuplot.  
 * \attention The localdg implementation uses the \c Dune::Pass
 * concept.
 *
 * @param argc number of arguments from command line
 * @param argv array of arguments from command line
 * @param envp array of environmental variables
 * @return 0 we don't program bugs. :)
 */
int main(int argc, char ** argv, char ** envp) {

  /* Parallelization is not implemented */
  MPIHelper::instance(argc,argv);

  try {

  // *** Initialization
  Parameter::append(argc,argv);      /*@\label{dg:param0}@*/
  if (argc == 2) {
    Parameter::append(argv[1]);
  } else {
    Parameter::append("parameter");  /*@\label{dg:paramfile}@*/
  }

  // initialize grid
  std::string filename;
  Parameter::get("femhowto.localdg.gridFile", filename);
  GridPtr<GridType> grid(filename); // ,MPI_COMM_WORLD);

  // ----- read in runtime parameters ------
  int eocSteps   = Parameter::getValue<int>("femhowto.localdg.eocSteps", 1);
  int startLevel = Parameter::getValue<int>("femhowto.localdg.startLevel", 0);
  int printCount = Parameter::getValue<int>("femhowto.localdg.printCount", -1);
  int verbose    = Parameter::getValue<int>("femhowto.localdg.verbose", 0);
  const double maxTimeStep = 
      Parameter::getValue("femhowto.localdg.maxTimeStep", std::numeric_limits<double>::max());
  std::string eocOutPath = Parameter::getValue<std::string>("femhowto.localdg.eocOutputPath",  /*@\label{dg:param1}@*/
                                                   std::string("."));            
  EocOutput eocOutput(eocOutPath + std::string("/eoc.tex"));
        
  // ----- read in model parameters ------
  double cfl; // CFL coefficient for SSP ODE Solver     /*@\label{dg:cfl0}@*/
  switch (order) 
  {
    case 0: cfl=0.9;  break;
    case 1: cfl=0.2;  break;
    case 2: cfl=0.15; break;
    case 3: cfl=0.05; break;
    case 4: cfl=0.09; break;
  }
  Parameter::get("femhowto.localdg.cfl", cfl, cfl);     /*@\label{dg:cfl1}@*/
  double startTime   = Parameter::getValue<double>("femhowto.localdg.startTime", 0.0);
  double endTime     = Parameter::getValue<double>("femhowto.localdg.endTime", 0.9);
        
  std::cout << " CFL : " << cfl << std::endl;
  
  
  // InitialDataType is a Dune::Operator that evaluates to $u_0$ and also has a
  // method that gives you the exact solution.
  InitialDataType problem;
  // Initialize the model      
  ModelType model(problem);
  // Initialize flux for advection discretization (UpwindFlux)
  FluxType convectionFlux(model);

  // Initialize Timer for CPU time measurements
  Timer timer;
  // Initialize variables needed for EOC computation
  double prevTime  = 0.;
  double error     = 0.;
  double prevError = 0.;
  int    level     = 0;
  double maxdt     = 0.;
  double mindt     = 1.e10;
  double averagedt = 0.;

  // Initialize L2Error for computation of error between discretized and exact
  // solution
  L2Error<DgType::DestinationType> L2err;
  FieldVector<double,ModelType::dimRange> err;


  // Refine the grid until the startLevel is reached
  for(int level=0; level < startLevel ; ++level)
    grid->globalRefine(DGFGridInfo<GridType>::refineStepsForHalf());

  /*********************************************
   * EOC Loop                                  *
   ********************************************/
  for(int eocloop=0; eocloop < eocSteps; ++eocloop) 
  {
    // Initialize DG Operator
    DgType dg(*grid, convectionFlux);                        /*@\label{dg:dgop}@*/
    // Initialize TimeProvider
    GridTimeProvider< GridType > tp(startTime, cfl, *grid);  /*@\label{dg:timeprovider}@*/
    // Initialize ODE Solver needed for the time discretization (Runge-Kutta)
    ODEType ode(dg, tp, rkSteps, verbose); 

    // Initialize the discrete Function $u$
    DgType::DestinationType U("U", dg.space());
    initialize(problem,U);

    // Print problem info to the eocOutput file
    if (eocloop==0) {
      eocOutput.printInput(problem, *grid, ode, filename);
    }
    
    // Initialize the DataWriter that writes the solution on the harddisk in a
    // format readable by Paraview e.g. 
    // IOTupleType is a Tuple of discrete functions
    IOTupleType dataTup ( &U );                                                    /*@\label{dg:datawriterdef0}@*/
    typedef DataWriter<GridType, IOTupleType> DataWriterType;
    DataWriterType dataWriter( *grid, filename, dataTup, startTime, endTime );     /*@\label{dg:datawriterdef1}@*/

    dataWriter.write(0.0, 0);
    
    FieldVector<double,ModelType::dimRange> projectionError =
      L2err.norm(problem, U, startTime);  
    std::cout << "Projection error " << problem.myName << ": " << projectionError <<
      std::endl;
        
    /**********************************************
     * Time Loop                                  *
     *********************************************/
    tp.provideTimeStepEstimate(maxTimeStep);              /*@\label{dg:maxtimestep1}@*/
    // ode.initialize applies the DG Operator once in order to get an estimate for dt.
    ode.initialize(U);                                    /*@\label{dg:odeinit}@*/
    for( tp.init() ; tp.time() < endTime ; tp.next() )
    {
      tp.provideTimeStepEstimate(maxTimeStep);            /*@\label{dg:maxtimestep2}@*/
      const double tnow  = tp.time();
      const double ldt   = tp.deltaT();
      const int counter  = tp.timeStep();

      /************************************************
       * Compute an ODE timestep                      *
       ***********************************************/
      ode.solve(U);                                      /*@\label{dg:odesolve}@*/

      if (!U.dofsValid()) {
        std::cout << "Invalid DOFs" << std::endl;
        if(eocloop == eocSteps-1) {
          dataWriter.write(1e10, counter+1);
        }
        abort();
      }

      if(verbose > 1 && printCount > 0 && counter % printCount == 0) {
        std::cout << "step: " << counter << "  time = " << tnow << ", dt = " << ldt << "\n";
      }

      if( eocloop == eocSteps -1 ) {
        dataWriter.write(tnow, counter);                 /*@\label{dg:datawrite1}@*/
      }
      
      // some statistics
      mindt = (ldt<mindt)?ldt:mindt;
      maxdt = (ldt>maxdt)?ldt:maxdt;
      averagedt += ldt;

      // Abort if the ODE solver does not converge
      if(counter%100 == 0) 
      {
        err = L2err.norm(problem, U, tp.time());
        if(err.one_norm() > 1e5 || ldt < 1e-10) 
        {
          averagedt /= double(counter);
          std::cout << "Solution doing nasty things!" << std::endl;
          std::cout << tnow << std::endl;
          eocOutput.printTexAddError(err[0], prevError, -1, grid->size(0), counter, averagedt);
          eocOutput.printTexEnd(timer.elapsed());
          exit(EXIT_FAILURE);
        }
      }

    } /****** END of time loop *****/

    averagedt /= double(tp.timeStep());
    if(verbose > 3)
    {
      std::cout << "Minimum dt: " << mindt 
         << "\nMaximum dt: " << maxdt 
         << "\nAverage dt: " << averagedt << endl;
    }

    // Write solution to hd                            
    if( eocloop == eocSteps-1 )
    {
      dataWriter.write(tp.time(), tp.timeStep());      /*@\label{dg:datawrite2}@*/
    }

    // Compute L2 error of discretized solution ...
    err = L2err.norm(problem, U, tp.time());
    std::cout << "Error " << problem.myName << ": " << err << endl;
    error       = err.two_norm();            
    double time = timer.elapsed() - prevTime;

    // ... and print the statistics out to the eocOutputPath file
    eocOutput.printTexAddError(error, prevError, time, grid->size(0), tp.timeStep(), averagedt);
    prevTime  = time;
    prevError = error;
   
    // Stop if too much time passed by (energy prices are soo high!)
    if(time > 3000.)
      break;
    
    // Refine the grid for the next EOC Step.
    if(eocloop < eocSteps-1) {
      grid->globalRefine(DGFGridInfo<GridType>::refineStepsForHalf());
      ++level;
    }
  } /***** END of EOC Loop *****/
  
  // Close the eocOutputPath file
  eocOutput.printTexEnd(timer.elapsed());

  Parameter::write("parameter.log");

  }
  catch (Dune::Exception &e) {
    std::cerr << e << std::endl;
    return 1;
  } catch (...) {
    std::cerr << "Generic exception!" << std::endl;
    return 2;
  }

  return 0;  
} 

