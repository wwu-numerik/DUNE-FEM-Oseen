/**
 *  \brief  runmanager.hh
 *
 *  \todo   doc
 **/

#ifndef RUNMANAGER_HH_INCLUDED
#define RUNMANAGER_HH_INCLUDED

#ifndef ENABLE_ALUGRID
    #define ENABLE_ALUGRID
#endif
#ifndef HAVE_ALUGRID
    #define HAVE_ALUGRID
#endif
#ifndef ALUGRID_SIMPLEX
    #define ALUGRID_SIMPLEX
#endif
#include <dune/grid/io/file/dgfparser/dgfalu.hh>

#include <dune/fem/gridpart/adaptiveleafgridpart.hh>
#include <dune/fem/io/file/vtkio.hh>
#include <dune/fem/space/dgspace.hh>

#include <dune/stokes/discretestokesmodelinterface.hh>
#include <dune/stokes/stokespass.hh>

#include <dune/stuff/misc.hh>
#include <dune/stuff/logging.hh>
#include <dune/stuff/grid.hh>
#include <dune/stuff/functions.hh>

#include "analyticaldarcydata.hh"

#include <fstream>


//template< class GridImp >
class RunManager
{
    public:

        typedef Dune::ALUSimplexGrid< gridDim, gridDim >
            GridType;

        typedef Dune::AdaptiveLeafGridPart< GridType >
            GridPartType;

    private:

        typedef Dune::GridPtr< GridType >
            GridPointerType;

        typedef Dune::DiscreteStokesModelDefaultTraits<
                        GridPartType,
                        Darcy::ConstantFunction,
                        Darcy::ConstantFunction,
                        gridDim,
                        MICRO_POLORDER,
                        MICRO_POLORDER,
                        MICRO_POLORDER >
            MicroStokesModelTraits;

        typedef Dune::DiscreteStokesModelDefault< MicroStokesModelTraits >
            MicroStokesModelType;

        typedef MicroStokesModelType::DiscreteStokesFunctionSpaceWrapperType
            MicroDiscreteStokesFunctionSpaceWrapperType;

        typedef MicroStokesModelType::DiscreteStokesFunctionWrapperType
            MicroDiscreteStokesFunctionWrapperType;

    public:

        typedef MicroDiscreteStokesFunctionWrapperType::DiscreteVelocityFunctionType
            DiscreteFunctionType;

    private:

        typedef MicroStokesModelTraits::AnalyticalForceType
            MicroAnalyticalForceType;

        typedef MicroStokesModelTraits::AnalyticalDirichletDataType
            MicroAnalyticalDirichletDataType;

        typedef Dune::StartPass< MicroDiscreteStokesFunctionWrapperType, -1 >
            MicroStartPassType;
        MicroStartPassType microStartPass;

        typedef Dune::StokesPass< MicroStokesModelType, MicroStartPassType, 0 >
            MicroStokesPassType;

    public:

        RunManager( const int verbose = 1, const std::string outputPrefix = "" )
            : hasReferenceSolution_( false ),
            verbose_( 0 ),
            outputPrefix_( "" ),
            dgfFileName_ (Dune::Parameter::getValue( "micro_grid_2d", std::string("micro_grid_2d.dgf") ) )
        {
            setVerboseLevel( verbose );
            setOutputPrefix( outputPrefix );
        }

        ~RunManager()
        {}

        void generateReferenceSolution( const int refineLevel )
        {
            // logging
            Logging::LogStream& info = Logger().Info();
            Logging::LogStream& debug = Logger().Dbg();
            Logging::LogStream& error = Logger().Err();

            debug << outputPrefix_ << "Computing reference Solution..." << std::endl;

            // grid
            debug << outputPrefix_ << "\tInitialising micro grid..." << std::endl;

//            const int refine_level = Dune::Parameter::getValue( "micro_reference_solution_refine", 0 ) * Dune::DGFGridInfo< GridType >::refineStepsForHalf();
            const int refine_level = refineLevel * Dune::DGFGridInfo< GridType >::refineStepsForHalf();
            referenceSolutionRefineLevel_ = refine_level;
            GridPointerType gridPointer( dgfFileName_ );
            gridPointer->globalRefine( refine_level );
            GridPartType microGridPart( *gridPointer );

            info << outputPrefix_ << "\tInitialised micro grid (with " << microGridPart.grid().size( 0 ) << " elements)." << std::endl;

            // problem

            debug << outputPrefix_ << "\tInitialising micro problem..." << std::endl;

            const double microViscosity = Dune::Parameter::getValue( "micro_viscosity", 1.0 );

            // function wrapper for the solutions
            MicroDiscreteStokesFunctionSpaceWrapperType
                microDiscreteStokesFunctionSpaceWrapper( microGridPart );

            MicroDiscreteStokesFunctionWrapperType
                microSolutionsX( "microX", microDiscreteStokesFunctionSpaceWrapper );

            MicroDiscreteStokesFunctionWrapperType
                microSolutionsY( "microY", microDiscreteStokesFunctionSpaceWrapper );

            MicroDiscreteStokesFunctionWrapperType
                dummy( "dummy", microDiscreteStokesFunctionSpaceWrapper );

            // analytical data and model
            // for x dimension
            MicroAnalyticalForceType::RangeType unitVectorX( 0.0 );
            unitVectorX[ 0 ] = 1.0;

            MicroAnalyticalForceType microAnalyticalForceX( microDiscreteStokesFunctionSpaceWrapper.discreteVelocitySpace(), unitVectorX );

            MicroAnalyticalDirichletDataType microAnalyticalDirichletDataX( microDiscreteStokesFunctionSpaceWrapper.discreteVelocitySpace(), unitVectorX );

            MicroStokesModelType microStokesModelX(  Dune::StabilizationCoefficients::getDefaultStabilizationCoefficients(),
                                                        microAnalyticalForceX,
                                                        microAnalyticalDirichletDataX,
                                                        microViscosity );

            // for y dimension
            MicroAnalyticalForceType::RangeType unitVectorY( 0.0 );
            unitVectorX[ 1 ] = 1.0;

            MicroAnalyticalForceType microAnalyticalForceY( microDiscreteStokesFunctionSpaceWrapper.discreteVelocitySpace(), unitVectorY );

            MicroAnalyticalDirichletDataType microAnalyticalDirichletDataY( microDiscreteStokesFunctionSpaceWrapper.discreteVelocitySpace(), unitVectorY );

            MicroStokesModelType microStokesModelY(  Dune::StabilizationCoefficients::getDefaultStabilizationCoefficients(),
                                                        microAnalyticalForceY,
                                                        microAnalyticalDirichletDataY,
                                                        microViscosity );

            // passes
            MicroStokesPassType microStokesPassX(   microStartPass,
                                                    microStokesModelX,
                                                    microGridPart,
                                                    microDiscreteStokesFunctionSpaceWrapper );

            MicroStokesPassType microStokesPassY(   microStartPass,
                                                    microStokesModelY,
                                                    microGridPart,
                                                    microDiscreteStokesFunctionSpaceWrapper );

            info << outputPrefix_ << "\tInitialised micro problem." << std::endl;

            debug << outputPrefix_ << "\tSolving micro system..." << std::endl;

            dummy.clear();
            microSolutionsX.clear();
            microSolutionsY.clear();
            microStokesPassX.apply( dummy, microSolutionsX );
            microStokesPassY.apply( dummy, microSolutionsX );

            info << outputPrefix_ << "\tMicro system solved." << std::endl;

            debug << outputPrefix_ << "\tSaving micro reference solution..." << std::endl;

            saveReferenceSolution( microSolutionsX.discreteVelocity(), std::string( Dune::Parameter::getValue( "micro_reference_solution_filename", std::string( "micro_reference_velocity" ) ) + "_X" + "_ref_" + Stuff::toString( referenceSolutionRefineLevel_ ) ) );
            saveReferenceSolution( microSolutionsY.discreteVelocity(), std::string( Dune::Parameter::getValue( "micro_reference_solution_filename", std::string( "micro_reference_velocity" ) ) + "_Y" + "_ref_" + Stuff::toString( referenceSolutionRefineLevel_ ) ) );

            info << outputPrefix_ << "\tMicro reference solution saved." << std::endl;

        }

        int setVerboseLevel( const int verbose )
        {
            Logging::LogStream& info = Logger().Info();
            Logging::LogStream& debug = Logger().Dbg();

            if ( verbose == 0) {
                info.Suspend();
                debug.Suspend();
                verbose_ = verbose;
                return 0;
            }
            else if ( verbose == 1 ) {
                info.Resume();
                debug.Suspend();
                verbose_ = verbose;
                return 0;
            }
            else if ( verbose == 2 ) {
                info.Resume();
                debug.Resume();
                verbose_ = verbose;
                return 0;
            }
            else {
                return 1;
            }
        }

        void setOutputPrefix( const std::string outputPrefix )
        {
            outputPrefix_ = outputPrefix;
        }

        int saveReferenceSolution( const DiscreteFunctionType& discreteFunction, const std::string filename )
        {
            Logging::LogStream& error = Logger().Err();

            int errorState = discreteFunction.write_ascii( filename + std::string( ".ddf" ) );

            // save dgf file with information
            std::ifstream readFile( dgfFileName_.c_str() );
            if ( readFile.is_open() ) {
                std::ofstream writeFile( std::string( filename + std::string( ".dgf" ) ).c_str() );
                if ( writeFile.is_open() ) {
                    // shift until "DGF"
                    while( !( readFile.eof() ) ) {
                        std::string line( "" );
                        std::getline( readFile, line );
                        if ( line.size() ) {
                            // remove lines with redundant information == "# ddf_"
                            if( !( ( line[0] == '#' ) && ( line[1] == ' ' ) && ( line[2] == 'd' ) && ( line[3] == 'd' ) && ( line[4] == 'f' ) && ( line[5] == '_' ) ) ) {
                                writeFile << line << std::endl;
                            }
                            if( ( line[0] == 'D' ) && ( line[1] == 'G' ) && ( line[2] == 'F' ) ) {
                                break;
                            }
                        }
                    }
                    // insert informations
                    writeFile << "# ddf_gridtype: " << "ALUGRID_SIMPLEX" << std::endl;
                    writeFile << "# ddf_refine_level: " << referenceSolutionRefineLevel_ << std::endl;
                    // rest of dgf file
                    while( !( readFile.eof() ) ) {
                        std::string line( "" );
                        std::getline( readFile, line );
                        if ( line.size() ) {
                            // remove lines with redundant information == "# ddf_"
                            if( !( ( line[0] == '#' ) && ( line[1] == ' ' ) && ( line[2] == 'd' ) && ( line[3] == 'd' ) && ( line[4] == 'f' ) && ( line[5] == '_' ) ) ) {
                                writeFile << line << std::endl;
                            }
                        }
                    }
                    writeFile.flush();
                    writeFile.close();
                    readFile.close();
                }
                else {
                    ++errorState;
                    error << "Error: could not open " << filename + std::string( ".dgf" ) << "!" << std::endl;
                }
            }
            else {
                ++errorState;
                error << "Error: could not open " << dgfFileName_ << "!" << std::endl;
            }
            return errorState;
        }

        int loadReferenceSolution( const std::string filename )//, DiscreteFunctionType& referenceSolution )
        {
            Logging::LogStream& error = Logger().Err();

            int errorState( 0 );

            std::string gridType( "" );
            int refineLevel( 0 );

            bool gridTypeRead( false );
            bool refineLevelRead( false );

            std::ifstream readFile( std::string( filename + std::string( ".dgf" ) ).c_str() );
            if ( readFile.is_open() ) {
                while( !( readFile.eof() ) ) {
                    if ( !( gridTypeRead || refineLevelRead ) ) {
                        std::string line( "" );
                        std::getline( readFile, line );
                        if ( line.size() ) {
                            if ( line.substr( 0, 6 ) == "# ddf_" ) {
                                unsigned int key_start = 6;
                                unsigned int key_end = key_start;
                                for ( ; key_end < line.size(); ++key_end ) {
                                    const char &c = line[key_end];
                                    if( (c == ' ') || (c == '\t') || (c == ':') ) {
                                        break;
                                    }
                                }
                                std::string key = line.substr( key_start, key_end - key_start );
                                unsigned int value_start = key_end;
                                for( ; value_start < line.size() ; ++value_start ) {
                                    if( line[value_start] == ':' ) {
                                        break;
                                    }
                                }
                                ++value_start;
                                for( ; value_start < line.size(); ++value_start ) {
                                    if( ( line[value_start] != ' ' ) && ( line[value_start] != '\t' ) ) {
                                        break;
                                    }
                                }
                                if( value_start >= line.size() ) {
                                    ++errorState;
                                }
                                std::string value = line.substr( value_start, line.size() - value_start );
                                if ( key == "gridtype" ) {
                                    gridType = value;
                                    gridTypeRead = true;
                                }
                                else if ( key == "refine_level" ) {
                                    refineLevel = Stuff::fromString< int >( value );
                                    refineLevelRead = true;
                                }
                            }
                        }
                    }
                    else {
                        std::cout << "gridtype is " << gridType << std::endl;
                        std::cout << "refine level is " << refineLevel << std::endl;
                        break;
                    }
                }
            }
            else {
                ++errorState;
                error << "Error: could not open " << std::string( filename + std::string( ".dgf" ) ) << "!" << std::endl;
            }
            return errorState;
        }

    private:

        int referenceSolutionRefineLevel_;
        bool hasReferenceSolution_;
        int verbose_;
        std::string outputPrefix_;
        const std::string dgfFileName_;
};

#endif
