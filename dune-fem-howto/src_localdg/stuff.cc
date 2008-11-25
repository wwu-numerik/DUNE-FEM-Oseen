/**
 *	\file	stuff.cc
 *	\brief	stuff.cc
 **/

#include  <fstream>
/* Helper class for projecting an analytic function onto a discrete function space */
template <class DiscreteFunctionType, class FunctionType, int polOrd>
class L2ProjectionLocal
{
  typedef typename DiscreteFunctionType::FunctionSpaceType FunctionSpaceType;
  
 public:
  static void project (const FunctionType &f, DiscreteFunctionType &discFunc) 
  {
    typedef typename DiscreteFunctionType::Traits::DiscreteFunctionSpaceType FunctionSpaceType;
    typedef typename FunctionSpaceType::Traits::GridPartType GridPartType;
    typedef typename FunctionSpaceType::Traits::GridType GridType;
    typedef typename GridType :: template Codim<0> :: Geometry Geometry;
    typedef typename FunctionSpaceType::Traits::IteratorType Iterator;
    
    const FunctionSpaceType& space =  discFunc.space();
    
    discFunc.clear();
    
    typedef typename DiscreteFunctionType::LocalFunctionType LocalFuncType;
    
    typename FunctionSpaceType::RangeType ret (0.0);
    typename FunctionSpaceType::RangeType phi (0.0);
    
    Iterator it = space.begin();
    Iterator endit = space.end();

    // if empty grid, do nothing 
    if( it == endit ) return ;
    
    // Get quadrature rule
    CachingQuadrature<GridPartType,0> quad(*it, 2*polOrd+1);
    
    const int quadNop = quad.nop();
    for( ; it != endit ; ++it) 
    {
      LocalFuncType lf = discFunc.localFunction(*it);

      const Geometry& geo = (*it).geometry();
      for(int qP = 0; qP < quadNop; ++qP) 
      {
        f.evaluate(geo.global(quad.point(qP)), ret);
        ret *=  quad.weight(qP);

        lf.axpy( quad[ qP ], ret );
      }
    }
  }
};

/**
 * @brief projects analytical function onto discrete function space
 *
 * Projects analytical function \c f onto discrete function space and stores
 * result in \c df
 *
 */
template <class StupidFunction,class DFType>
void initialize(const StupidFunction& f,DFType& df)
{
  //- Typedefs and enums
  typedef typename DFType::Traits::DiscreteFunctionSpaceType SpaceType;
  typedef typename SpaceType::Traits::IteratorType Iterator;
  typedef typename DFType::Traits::LocalFunctionType LocalFunction;
  typedef typename SpaceType::Traits::GridType Grid;
  typedef typename Grid::template Codim<0>::Entity::Geometry Geometry;

  typedef typename SpaceType::Traits::RangeType RangeType;
  typedef typename SpaceType::Traits::DomainType DomainType;

  enum { dim = Grid::dimension };

  typedef FieldVector<double, dim> Coordinate;

  //- Actual method
  L2ProjectionLocal<DFType, StupidFunction, 2>::project(f, df);

  typedef typename DFType::DofIteratorType DofIterator;
  /*for (DofIterator it = df.dbegin(); it != df.dend(); ++it) {
    std::cout << *it << std::endl;
    }*/
}

/**
 * @brief EocOutput is a helper class for making a LaTeX output file with an EOC error table.
 *
 * \c EocOutput writes a LaTeX file whose name is given in the constructor and
 * fills it with EOC data given to an instance of \c EocOutput by the function
 * \c printTexAddError.
 *
 */
class EocOutput {

  std::string outputFile;
	
 public:
		
  /**
   * @brief Constructor of EocOutput
   *
   * Opens the file given by \c name.
   *
   * @param name name of file to write to.
   */
  EocOutput(std::string name)
  {
    outputFile = name;
    std::cout << outputFile << std::endl;
    mode_t mode =   S_IRUSR | S_IWUSR | S_IXUSR |
                    S_IRGRP | S_IWGRP | S_IXGRP |
                    S_IROTH | S_IXOTH;
    // get dirname of outputFile
    std::string dir = outputFile.substr(0,outputFile.find_last_of("/"));
    // test if dir exists
    std::ifstream fin(dir.c_str());
    if( fin == NULL ) {
      std::cerr << "EocOutput path not existent. Creating..." << std::endl;
    } else {
      fin.close();
    }

    //the dir to be created is what stands after 
    //the last "/" in given outputpath
    std::string sub_dir = outputFile.substr(0, outputFile.find_last_of("/"));
    if ( 0 != mkdir( sub_dir.c_str(), mode ) ) {
      std::cerr << "Failed creating dir: " << sub_dir << std::endl;
    }

    std::ostringstream filestream;
    filestream << outputFile;

    std::ofstream ofs(filestream.str().c_str(), std::ios::out);

    ofs << "\\documentclass[12pt,english]{article}\n"
	<< "\\usepackage[T1]{fontenc}\n"
	<< "\\usepackage[latin1]{inputenc}\n"
	<< "\\usepackage{setspace}\n"
	<< "\\onehalfspacing\n"
	<< "\\makeatletter\n"
	<< "\\providecommand{\\boldsymbol}[1]{\\mbox{\\boldmath $#1$}}\n"
	<< "\\providecommand{\\tabularnewline}{\\\\}\n"
	<< "\\usepackage{babel}\n"
	<< "\\makeatother\n"
	<< "\\begin{document}\n";
				 
    ofs.close();	
  }
	
  /**
   * @brief closes the file
   *
   * Closes the file and adds statistics about the total cpu-time.
   *
   * @param totaltime total cpu time elapsed
   */
  void printTexEnd(double totaltime)
  {
    std::ostringstream filestream;
    filestream << outputFile;

    std::ofstream ofs(filestream.str().c_str(), std::ios::app);
		
    ofs  << "\\end{tabular}\\\\\n\n"
	 << "Total time: " << totaltime << "\n"
	 << "\\end{document}\n" << std::endl;
		
    ofs.close();
  }
	
  /**
   * @brief adds a line in the eoc error table
   *
   * Adds a line in the eoc error table with data about the error and the
   * average time discretization step at the current refinement level, and the
   * EOC by taking into account the previous discretization error.
   *
   * @param error error at current refinement level
   * @param prevError previous error
   * @param time cpu-time elapsed for this solution
   * @param level current refinement level
   * @param counter number of timesteps
   * @param averagedt average time discretization step length
   */
  void printTexAddError(double error, double prevError, double time, int level, int counter,double averagedt)
  {
    std::ostringstream filestream;
    filestream << outputFile;

    std::ofstream ofs(filestream.str().c_str(), std::ios::app);
		
    if(prevError > 0.0)
      {	       
	ofs <<  "\\hline \n"
	    << level << " & " << error << " & " << log(prevError/error)/M_LN2 << " & " << time << " & " << counter << " & " << averagedt << "\n"
	    << "\\tabularnewline\n"
	    << "\\hline \n";
      }
    else
      {	       
	ofs << "\\begin{tabular}{|c|c|c|c|c|c|}\n"
	    << "\\hline \n"
	    << "Size & $\\left\\Vert u-u_{h}\\right\\Vert _{L_{2}}$ & EOC & CPU & \\#Iterations & a-dt\n"
	    << "\\tabularnewline\n"
	    << "\\hline\n"
	    << "\\hline\n"
	    << level << " & " << error << " & " << "---" << " & " << time << " & " << counter << " & " << averagedt << "\n"
	    << "\\tabularnewline\n"
	    << "\\hline \n";
      }
		
    ofs.close();
  }

		
  /**
   * @brief prints a header into the LaTeX file
   *
   * Prints out information about the Model that was computed.
   *
   * @param u0 an instance of the computed problem
   * @param grid grid instance
   * @param ode ode solver
   * @param arg file name of the dgf file used
   */
  void printInput(InitialDataType& u0, GridType& grid, ODEType& ode, 
		  std::string arg)
  {
    std::ostringstream filestream;
    filestream << outputFile;

    std::ofstream ofs(filestream.str().c_str(), std::ios::app);

    ofs  << "Grid: " << grid.name() << "\n\n"
	 << "Macrogrid: " << arg << "\\\\\n\n";
		
    ofs.close();
    u0.printmyInfo(outputFile);		
  }
};
