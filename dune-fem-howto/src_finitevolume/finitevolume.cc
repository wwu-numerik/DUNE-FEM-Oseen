#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <iostream>               // for input/output to shell
#include <sstream>
#include <fstream>                // for input/output to files
#include <vector>                 // STL vector class


#include <dune/common/mpihelper.hh> // An initializer of MPI
#include <dune/common/exceptions.hh> // We use exceptions
#include <dune/grid/io/file/dgfparser/dgfgridtype.hh> // for the grid

#include <fem/space/fvspace.hh>             //FiniteVolumeSpace
//#include <fem/space/common/adaptiveleafgridpart.hh>  //for adaptiveleafgridpart
//#include <fem/space/dgspace/dgadaptiveleafgridpart.hh>  //for DGAdaptiveLeafGridPart
#include <dune/fem/function/adaptivefunction.hh>  //Adaptive Function
#include <dune/fem/function/blockvectorfunction.hh>  //Blockvector Function
// checks for defined gridtype and includes appropriate dgfparser implementation 
#include <dune/common/timer.hh>   // timer 

#include <dune/fem/io/file/vtkio.hh>

//#include "vtkout.hh"
//#include "grapeviz.hh"
#include "transportproblem.hh"
#include "finitevolumescheme.hh"
#include "adaptation.hh"
#include "norm.hh"

#ifndef WANT_GRAPE
  #define WANT_GRAPE 0
#endif

#ifndef POLORDER
  #define POLORDER 0
#endif

#if ENABLE_MPI
        typedef Dune::CollectiveCommunication< MPI_Comm > CollectiveCommunication;
#else
        typedef Dune::CollectiveCommunication< double > CollectiveCommunication;
#endif

// use namespace Dune 
using namespace Dune;

//===============================================================
// the time loop function working for all types of grids
//===============================================================
template<bool adaptive, class GridType>
//void timeloop (GridType& grid, int minLevel,int maxLevel, const double endTime, bool display)
void timeloop (GridType& grid, int minLevel,int maxLevel, const double endTime)
{
  /*
  // in case of structured grid do nothing
  if( Capabilities::IsUnstructured<GridType> :: v == false || ! adaptive)
  {
    // in case of structured grids choose other min level 
    minLevel += 1 * DGFGridInfo<GridType>::refineStepsForHalf(); 
  }
  */

  // *** refine grid until upper limit of level 
  if (adaptive)
    grid.globalRefine(minLevel);
  else
    grid.globalRefine(maxLevel);


  // *** some important typedefs:
//  typedef DGAdaptiveLeafGridPart<GridType> GridPartType;  // use leaf grid
// typedef Dune::LeafGridPart<GridType> GridPartType;
  typedef HierarchicGridPart<GridType> GridPartType;
  GridPartType gridPart(grid);

  //problem type
  typedef ProblemImp<GridType> Problem;

  //extract ProblemType (see transportproblem.hh)
  typedef typename Problem::FunctionSpaceType FunctionSpaceType;
  

  // use linear upwind flux (See transportproblem.hh)
  typedef LinearUpwindFlux< FieldVector<double,GridType::dimension> , double> NumericalFluxType;  
  NumericalFluxType flux;
  
  // define the discrete function space our unkown belongs to, here we use a Finite Volume Space
  typedef FiniteVolumeSpace < FunctionSpaceType, GridPartType, POLORDER >  DiscreteFunctionSpaceType;
  
  typedef AdaptiveDiscreteFunction < DiscreteFunctionSpaceType > DiscreteFunctionType;
  // typedef BlockVectorDiscreteFunction < DiscreteFunctionSpaceType > DiscreteFunctionType;
  
  // type of adaptation class (See adaption.hh)
  typedef Adaptation< DiscreteFunctionType , Capabilities::IsUnstructured<GridType> :: v && adaptive> AdaptationType;
  //types needed for adaption
  //typedef RestrictProlongDiscontinuousSpace< DiscreteFunctionType, POLORDER >  RestrictProlongType; //Polorder = 0
  typedef RestrictProlongPieceWiseConstantData< DiscreteFunctionType >  RestrictProlongType;
  //typedef RestrictProlongDefaultImplementation< DiscreteFunctionType, DiscreteFunctionSpaceType >  RestrictProlongType;
  typedef AdaptationManager < GridType, RestrictProlongType > ADOperatorType;
  
  // type of finite volume scheme (See finitevolumescheme.hh)
  typedef FiniteVolumeScheme<DiscreteFunctionType> FiniteVolumeSchemeType;
  
  // *** some initialation stuff

  //Test initialization for discrete function of FiniteVolumeSpace
  DiscreteFunctionSpaceType discreteFunctionSpace ( gridPart );
  DiscreteFunctionType solution ( "solution" , discreteFunctionSpace );
  FiniteVolumeSchemeType::initialize(solution);  // initialize concentration with initial values
  
  /******************************************************************************************  Wird eh unten uberschrieben
  std::string name("solution");
  VTKIO<GridPartType> vtkio( gridPart );
  vtkio.addCellData( solution ); 
  vtkio.write(name.c_str(),Dune::VTKOptions::ascii);
  ***************************************************************************************************************************/  

  
  //create  Adaption-Types and initialize adaption parameters
  RestrictProlongType rp(solution);                                        //eigentlich nur noetig, wenn adaptiv gerechnet werden soll!!!!!!!!!
  ADOperatorType adop(grid,rp);
  
  // tolerance value for refinement strategy
  const double refineTol  = 0.1;
  // tolerance value for coarsening strategy
  const double coarsenTol = 0.1 * refineTol;

   // make initial adaptation if adaptive
  if (adaptive) {
    for(int i=minLevel; i<maxLevel; ++i)
      { 
	// adapt grid to avoid initial error 
//	std::cout<< "Round "<<i-minLevel+1<<"..."<<std::endl;                            <-to be deleted
//	std::cout.flush();
	// AdaptationType::adaptFem(gridPart,solutionFem,minLevel,maxLevel);
	
	// mark entities for coarsening and refinement 
//	std::cout<<"  Marking entities...";                                                 <-to be deleted
	// mark entities for coarsening and refinement 
	AdaptationType::markEntities(gridPart,solution,minLevel,maxLevel,refineTol,coarsenTol);
//	std::cout<<"  done"<<std::endl; <-to be deleted
	
	//call adaptation method
	adop.adapt();
    
	// initialize concentration with initial values
//	std::cout<< "Round "<<i-minLevel+1<<"done"<<std::endl;                                  <-to be deleted
//	std::cout.flush();
	solution.clear();
	FiniteVolumeSchemeType::initialize(solution);
      }
  }
  //if (0)
  {  
    const std::string name("solution");
    VTKIO<GridPartType> vtkio( gridPart );
    vtkio.addCellData( solution ); 
    vtkio.write(name.c_str(),Dune::VTKOptions::ascii);
  }
  #if WANT_GRAPE
  GrapeDataDisplay< GridType > grape( grid );
//  grape.dataDisplay( solution );
  bool isVectorial = false;
  grape.addData(solution, POLORDER, isVectorial);
  grape.display();
  #endif
//   std::cout<< "nach schreiben der initialwerte"<<std::endl;                                 <-to be deleted
//   std::cout.flush(); 
 

  // now do the time steps
  double time = 0.0;
  int k=0;
  const double saveInterval = 0.01; 
  double saveStep = saveInterval;
  int counter = 1;
  //update Function for the update each time step
  DiscreteFunctionType update ( "update" , discreteFunctionSpace );
  
     
  // ********** time loop **************
  while (time < endTime)
  {
    ++k;    // augment time step counter 
    // ******* apply finite volume scheme 
    double dt = FiniteVolumeSchemeType::
      finiteVolumeScheme(flux,solution,update,time,0.9); 
    
    // ******* update solution 
  //////  FiniteVolumeSchemeType::printFunctionValues(solution);
    update *= dt;
    solution += update;

    
    // get global min of dt (parallel stuff)
//    dt = gridPart.grid().comm().min( dt );                             
/*********
   

    // ******* increase time
**************/
    time += dt;
    // ******* check if data should be written
    if (time >= saveStep)
    {
     
      // output of Fem-Data 
      #if 1
      {
	std::stringstream nameStream;
	nameStream <<"solution";
	if (counter<10) nameStream << "0";
	nameStream<<counter;
	//	char fname[128];
	//      sprintf(fname,"%s-%05d",name.c_str(),counter);
	std::string name = nameStream.str();
	VTKIO<GridPartType> vtkio( gridPart );
	vtkio.addCellData( solution ); 
	vtkio.write(name.c_str(),Dune::VTKOptions::ascii);
      }
      #endif
      #if WANT_GRAPE
      if (counter==30) {
	GrapeDataDisplay< GridType > grape( grid );
	bool isVectorial = false;
	grape.addData(solution, POLORDER, isVectorial);
	grape.display();
	//grape.dataDisplay( solution );
      }
      #endif

      //set saveStep for next intervall
      //saveStep += saveInterval;
      //++counter;
    }

    // ********** call adaptation algorithm
    if (adaptive)
      {
	// mark entities for coarsening and refinement 
	AdaptationType::markEntities(gridPart,solution,minLevel,maxLevel,refineTol,coarsenTol);
    
	//call adaptation method
	adop.adapt();          //this method manages the adaptionprocess: gridrefinement,restriction and prolongation of the solution function
      }
    
    //write out solution file after adaptation process
    if (time >= saveStep)
      {
	if (adaptive) {
	  std::stringstream nameStream;
	  nameStream <<"solution";
	  if (counter<10) nameStream << "0";
	  nameStream<<counter<<"adapted";
	  //	char fname[128];
	  //      sprintf(fname,"%s-%05d",name.c_str(),counter);
	  std::string name = nameStream.str();
	  VTKIO<GridPartType> vtkio( gridPart );
	  vtkio.addCellData( solution ); 
	  vtkio.write(name.c_str(),Dune::VTKOptions::ascii); 	
	}	
	saveStep += saveInterval;
	++counter;
      }
    

/*****************
    //  communicate data (parallel stuff)
    {
      typedef typename GridPartType :: IndexSetType IndexSetType;
      // create data handle 
      FVDataHandle<GridType,IndexSetType,std::vector<double> > dataHandle(gridPart.indexSet(),solution);  //FEM

      // communicate data 
      gridPart.communicate( dataHandle, InteriorBorder_All_Interface, ForwardCommunication );    //FEM
    }
**********************/
    // print info about time, timestep size and counter  
    std::cout << "elements = " << gridPart.indexSet().size(0);
    std::cout << "   step = " << k << "   time = " << time << "   dt = " << dt << std::endl;
      
  }
  // ********** end of the time loop ****************

//print exact solution on the same gridPart
  DiscreteFunctionType exactSolution ( "ExactSolution" , discreteFunctionSpace );
  FiniteVolumeSchemeType::calcExactSolution( exactSolution,time);  // calculate exact solution
  
  std::string  name = "ExactSolution"; 
  
  VTKIO<GridPartType> vtkio( gridPart );
  vtkio.clear();
  vtkio.addCellData(exactSolution); 
  vtkio.write(name.c_str(),Dune::VTKOptions::ascii);

  DiscreteFunctionType errorFunction ( "ErrorFunction" , discreteFunctionSpace );
  errorFunction.clear();
  errorFunction += exactSolution;
  errorFunction -= solution;
  double l1error = Norm::L1Norm(errorFunction);
  std::cout << "L1-Fehler = " << l1error << std::endl;
}

//===============================================================
// The main function creates objects and does the time loop
//===============================================================
int main (int argc , char ** argv)
{ 
  // initialize MPI, finalize is done automatically on exit (parallel stuff)
  MPIHelper::instance(argc,argv);

  // start try/catch block to get error messages from dune
  try {
    // read macrogrid (unitcube) from file
    std::stringstream dgfFileName;
    dgfFileName << "../macrogrids/unitcube" << GridType :: dimension << ".dgf";
    GridPtr<GridType> gridPtr( dgfFileName.str() );    // create grid pointer, GridType is defined by gridtype.hh
    GridType& grid = *gridPtr;    // grid reference 
    
      // balance grid (parallel stuff) 
    grid.loadBalance();

    // Parameters for min- and max-Level Refinement
 //    int minLevel = 1 * DGFGridInfo<GridType>::refineStepsForHalf();    // half grid 1 time (initial refinement) 
    int minLevel = 2 * DGFGridInfo<GridType>::refineStepsForHalf();    // half grid 1 time (initial refinement)
 //   int maxLevel = minLevel + 5 * DGFGridInfo<GridType>::refineStepsForHalf();    // allow maximum level of refinement of minLevel + 5 * half factor
    int maxLevel = minLevel + 4* DGFGridInfo<GridType>::refineStepsForHalf();    // allow maximum level of refinement of minLevel + 5 * half factor


    // define some more parameters
    const bool adaptation = true;
    const double endTime = 0.5;
  //  bool grape = false; // toggle grape display  

    // start timer
    Timer timer;

    // ****** do time loop until end time ********
    //timeloop<adaptation>(grid,minLevel,maxLevel,endTime,grape);
    timeloop<adaptation>(grid,minLevel,maxLevel,endTime);
    // ******

    std::cout << "Simulation finished: cputime = " << timer.elapsed() << " sec!" << std::endl;
  }
  catch (std::exception & e) {
    std::cout << "STL ERROR: " << e.what() << std::endl;
    return 1;
  }
  catch (Dune::Exception & e) {
    std::cout << "DUNE ERROR: " << e.what() << std::endl;
    return 1;
  }
  catch (...) {
    std::cout << "Unknown ERROR" << std::endl;
    return 1;
  }

  // done
  return 0;
}
