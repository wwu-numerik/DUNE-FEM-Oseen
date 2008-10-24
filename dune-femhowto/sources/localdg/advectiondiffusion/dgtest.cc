// always include this first 
#include <config.h>

//- system includes 
#include <iostream>
#include <string>

//- dune-common includes 
#include <dune/common/misc.hh>
#include <dune/fem/misc/utility.hh>
#include <dune/common/timer.hh>
#include <dune/common/mpihelper.hh>

//- dune-grid includes 
#include <dune/grid/common/gridpart.hh>
#if HAVE_GRAPE
#include <dune/grid/io/visual/grapedatadisplay.hh>
#endif

//- dune-fem includes 
#include <dune/fem/io/file/grapedataio.hh>

// l2 projection 
#include <dune/fem/operator/projection/l2projection.hh>

// l2 error class 
#include <dune/fem/misc/l2error.hh>

//- local includes 
// in models the type of application is chosen 
#include "models.hh"
// some helper functions 
#include "stuff.cc"

#include "artdiff.cc"

#include <dune/fem/operator/projection/vtxprojection.hh> 
#include <dune/fem/operator/projection/p1projection.hh> 

int main(int argc, char ** argv, char ** envp) 
{
  // initialize MPI if available 
  MPIHelper::instance(argc,argv);

  // try and catch exceptions 
  try {
    
    // *** Initialization
    if (argc<2) {
      cout << "Call: dgtest gridfilename [ref-steps=1] [start-level=0] [epsilon=0.01] [use-grape=0]" << endl;
      exit(EXIT_FAILURE);  
    }
    
    // Polynomial and ODE order
    // Grid:
    GridPtr<GridType> grid(argv[1]); // ,MPI_COMM_WORLD);
    
    int repeats = 1;
    if (argc>2)
      repeats=atoi(argv[2]);
    
    int startlevel = 0;
    if (argc>3)
      startlevel=atoi(argv[3]);
    
    double epsilon = 0.01;
    if (argc>4)
      epsilon=atof(argv[4]);
    
    int graped = 0;
    if (argc>5)
      graped=atoi(argv[5]);
      
    // CFL: default value
    double cfl = 0.45 / (2.0 * order+1);;

    if (argc>6)
      cfl=atof(argv[6]);
    
    InitialDataType problem(epsilon,true);
    
    string myoutput = "eoc.tex";
    EocOutput eocoutput(myoutput);
    // *** Models
    ModelType model(*grid,problem);
    // *** Fluxes 
    FluxType eulerflux(model);

    Timer timer;
    int maxit=0;
    double zeit=0.;
    double prevzeit=0.;
    double fehler=0.;
    double prevfehler=0.;
    
    L2Error<DgType::DestinationType> L2err;
    FieldVector<double,ModelType::dimRange> err;

    // printSGrid(0,0,dg.space(),U);
    
    //int n=0;
    //double nextsave=0.;
    //double savestep=0.05;
    double maxtime = problem.endtime();
    int level=0;
    
    FieldVector<double,dimworld> upwind(1.);
    upwind[0] *= 0.37;
    
    for(level=0; level < startlevel ; ++level)
      grid->globalRefine(DGFGridInfo<GridType>::refineStepsForHalf());

    std::vector<double> errors (repeats);

    for(int eocloop=0;eocloop < repeats; ++eocloop) {
      // *** Operator typedefs
      DgType dg(*grid,eulerflux,upwind);
      ODEType ode(dg,rksteps,cfl); 

      // *** solution data
      DgType::DestinationType U("U", dg.space());

      GridPartType gridPart( *grid );
      ArtDiffSpace artSpace ( gridPart ); 
      ArtificialDiffusionType artDiff( "artdiff", artSpace );
      artDiff.clear();

      // initialize data 
      initialize(problem,U);

      if (eocloop==0) 
        eocoutput.printInput(problem,*grid,ode,argv[1]);
      
      /*
      LagrangeSpaceType lSpace ( gridPart ); 
      P1FunctionType p1Function( "p1-func", lSpace );
      //WeightDefault<GridPartType> weight;
      P1ProjectionImpl::project( U , p1Function );
      p1Function.print(std::cout);
      */

#if HAVE_GRAPE 
      if(graped > 0) {
        GrapeDataDisplay< GridType > grape(U.space().gridPart());
        //grape.addData( p1Function );
        grape.dataDisplay(U);
      }
#endif 
      
      double t=0.0;
      int counter=0;
      FieldVector<double,ModelType::dimRange> projectionError = 
        L2err.norm(problem,U,t);  
      cout << "Projection error " << problem.myName << ": " << projectionError << endl;
    
      double maxdt=0.,mindt=1.e10,averagedt=0.;
      // *** Time loop
      while (t<maxtime) 
      {
        //ArtificialDiffusion :: calculate ( U, artDiff );
        double ldt = -t;
        t = ode.solve(U);
        if( U.space().order() > 0 )
        {
          DgType::DestinationType tmp( U ); 
          dg.limit(tmp,U);
        }

        // switch diffusion vector
        dg.switchupwind();

#if HAVE_GRAPE 
        if(graped > 0)
        {
          if(counter%graped == 0 && counter > 0) 
          {
            GrapeDataDisplay< GridType > grape(U.space().gridPart());
            grape.dataDisplay(U);
          }
        }
#endif 
        ldt += t;
        std::cout << "Current time = " << t << "\n";
        mindt = (ldt<mindt)?ldt:mindt;
        maxdt = (ldt>maxdt)?ldt:maxdt;
        averagedt += ldt;
        
        if(0 && counter%100 == 0) 
        {
          err = L2err.norm(problem,U,t);
          if(err.one_norm() > 1e5 || ldt < 1e-10) 
          {
            averagedt /= double(counter);
            cout << "Solution doing nasty things!" << std::endl;
            cout << t << endl;
            eocoutput.printTexAddError(err[0],prevfehler,-1,grid->size(0),counter,averagedt);
            eocoutput.printTexEnd(timer.elapsed());
            exit(EXIT_FAILURE);
          }
        }
        ++counter;
      }
      
      averagedt /= double(counter);

      err = L2err.norm(problem,U,t);
      errors[eocloop] = err;
      cout << "Error " << problem.myName << ": " << err << endl;
      if(eocloop > 0)
      {
        double eoc = log(errors[eocloop-1]/errors[eocloop])/M_LN2;
        cout << "EOC " << problem.myName << ": " << eoc << endl;
      }
      
      fehler = err.two_norm();		
      zeit = timer.elapsed()-prevzeit;
      eocoutput.printTexAddError(fehler,prevfehler,zeit,grid->size(0),counter,averagedt);
     
#if HAVE_GRAPE
      if(graped>0)
      {
        GrapeDataDisplay< GridType > grape(*grid);
        grape.dataDisplay(U);
      }
#endif
      
      //if(zeit > 3000.)
      //  break;
      
      if(eocloop < repeats-1) {
        grid->globalRefine(DGFGridInfo<GridType>::refineStepsForHalf());
        ++level;
      }
      prevzeit = zeit;
      prevfehler = fehler;
      ++maxit;
    }
    
    eocoutput.printTexEnd(timer.elapsed());

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
