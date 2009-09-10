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
#include <dune/fem/misc/l2norm.hh>
#include <dune/fem/operator/projection/l2projection.hh>
#include <dune/fem/space/dgspace.hh>

#include <dune/stokes/discretestokesmodelinterface.hh>
#include <dune/stokes/stokespass.hh>

#include <dune/stuff/misc.hh>
#include <dune/stuff/logging.hh>
#include <dune/stuff/grid.hh>
#include <dune/stuff/functions.hh>

//#ifndef MICRO_PROBLEM
//    #define MICRO_PROBLEM MICRO_COCKBURN_PROBLEM
//    #warning #warning ("MICRO_PROBLEM undefined, defaulting to MICRO_COCKBURN_PROBLEM")
//#endif

//#if MICRO_PROBLEM==MICRO_COCKBURN_PROBLEM
//    #define COCKBURN_PROBLEM
//    #include "../src/analyticaldata.hh"
//#else
    #include "analyticaldarcydata.hh"
//#endif

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

//#if MICRO_PROBLEM==MICRO_COCKBURN_PROBLEM
//        typedef Dune::DiscreteStokesModelDefaultTraits<
//                        GridPartType,
//                        Force,
//                        DirichletData,
//                        gridDim,
//                        MICRO_POLORDER,
//                        MICRO_POLORDER,
//                        MICRO_POLORDER >
//            MicroStokesModelTraits;
//#else
        typedef Dune::DiscreteStokesModelDefaultTraits<
                        GridPartType,
                        Darcy::ConstantFunction,
                        Darcy::ConstantFunction,
                        gridDim,
                        MICRO_POLORDER,
                        MICRO_POLORDER,
                        MICRO_POLORDER >
            MicroStokesModelTraits;
//#endif


        typedef Dune::DiscreteStokesModelDefault< MicroStokesModelTraits >
            MicroStokesModelType;

        typedef MicroStokesModelType::DiscreteStokesFunctionSpaceWrapperType
            MicroDiscreteStokesFunctionSpaceWrapperType;

        typedef MicroStokesModelType::DiscreteStokesFunctionWrapperType
            MicroDiscreteStokesFunctionWrapperType;

    public:

        typedef MicroDiscreteStokesFunctionSpaceWrapperType::DiscreteVelocityFunctionSpaceType
            DiscreteFunctionSpaceType;

        typedef MicroDiscreteStokesFunctionWrapperType::DiscreteVelocityFunctionType
            DiscreteFunctionType;

    private:

        typedef MicroStokesModelTraits::AnalyticalForceType
            MicroAnalyticalForceType;

        typedef MicroStokesModelTraits::AnalyticalDirichletDataType
            MicroAnalyticalDirichletDataType;

        typedef Dune::StartPass< MicroDiscreteStokesFunctionWrapperType, -1 >
            MicroStartPassType;

        MicroStartPassType microStartPass_;

        typedef Dune::StokesPass< MicroStokesModelType, MicroStartPassType, 0 >
            MicroStokesPassType;

        typedef Dune::VTKIO< GridPartType >
            VTKWriterType;
    public:

        RunManager( const int verbose = 1, const std::string outputPrefix = "" )
            : hasReferenceSolution_( false ),
            verbose_( 0 ),
            outputPrefix_( "" ),
            dgfFileName_ (Dune::Parameter::getValue( "micro_grid_2d", std::string("micro_grid_2d.dgf") ) ),
            referenceSolution_( NULL )
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
            info.Flush();
            debug.Flush();
            error.Flush();

//#if MICRO_PROBLEM==MICRO_COCKBURN_PROBLEM
//            debug << outputPrefix_ << "Computing COCKBURN reference Solution..." << std::endl;
//#else
            debug << outputPrefix_ << "Computing reference Solution..." << std::endl;
//#endif

            // grid
            debug << outputPrefix_ << "Initialising micro grid..." << std::endl;

            referenceSolutionRefineLevel_ = refineLevel;
            GridPointerType gridPointer( dgfFileName_ );
            gridPointer->globalRefine( refineLevel );
            GridPartType microGridPart( *gridPointer );

            info << outputPrefix_ << "Initialised micro grid (with " << microGridPart.grid().size( 0 ) << " elements)." << std::endl;

            // problem

            debug << outputPrefix_ << "Initialising micro problem..." << std::endl;

            const double microViscosity = Dune::Parameter::getValue( "micro_viscosity", 1.0 );

            // function wrapper for the solutions
            MicroDiscreteStokesFunctionSpaceWrapperType
                microDiscreteStokesFunctionSpaceWrapper( microGridPart );

//#if MICRO_PROBLEM==MICRO_COCKBURN_PROBLEM
//            MicroDiscreteStokesFunctionWrapperType
//                microSolutions( "micro", microDiscreteStokesFunctionSpaceWrapper );
//
//            MicroDiscreteStokesFunctionWrapperType
//                dummy( "dummy", microDiscreteStokesFunctionSpaceWrapper );
//
//            MicroAnalyticalForceType microAnalyticalForce( microViscosity , microDiscreteStokesFunctionSpaceWrapper.discreteVelocitySpace() );
//
//            MicroAnalyticalDirichletDataType microAnalyticalDirichletData( microDiscreteStokesFunctionSpaceWrapper.discreteVelocitySpace() );
//
//            MicroStokesModelType microStokesModel(  Dune::StabilizationCoefficients::getDefaultStabilizationCoefficients(),
//                                                        microAnalyticalForce,
//                                                        microAnalyticalDirichletData,
//                                                        microViscosity );
//#else
            MicroDiscreteStokesFunctionWrapperType
                microSolutionsX( "microX", microDiscreteStokesFunctionSpaceWrapper );

            MicroDiscreteStokesFunctionWrapperType
                microSolutionsY( "microY", microDiscreteStokesFunctionSpaceWrapper );

            MicroDiscreteStokesFunctionWrapperType
                dummy( "dummy", microDiscreteStokesFunctionSpaceWrapper );
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
//#endif

//#if MICRO_PROBLEM==MICRO_COCKBURN_PROBLEM
//            // passes
//            MicroStokesPassType microStokesPass(   microStartPass_,
//                                                    microStokesModel,
//                                                    microGridPart,
//                                                    microDiscreteStokesFunctionSpaceWrapper );
//
//#else
            MicroStokesPassType microStokesPassX(   microStartPass_,
                                                    microStokesModelX,
                                                    microGridPart,
                                                    microDiscreteStokesFunctionSpaceWrapper );

            MicroStokesPassType microStokesPassY(   microStartPass_,
                                                    microStokesModelY,
                                                    microGridPart,
                                                    microDiscreteStokesFunctionSpaceWrapper );
//#endif

            info << outputPrefix_ << "Initialised micro problem." << std::endl;

            debug << outputPrefix_ << "Solving micro system..." << std::endl;

            dummy.clear();

//#if MICRO_PROBLEM==MICRO_COCKBURN_PROBLEM
//            microSolutions.clear();
//            microStokesPass.apply( dummy, microSolutions );
//#else
            microSolutionsX.clear();
            microSolutionsY.clear();
            microStokesPassX.apply( dummy, microSolutionsX );
//            microStokesPassY.apply( dummy, microSolutionsX );
//#endif

            info << outputPrefix_ << "Micro system solved." << std::endl;

            debug << outputPrefix_ << "Saving micro reference solution..." << std::endl;

//#if MICRO_PROBLEM==MICRO_COCKBURN_PROBLEM
//            saveReferenceSolution(
//                microSolutions.discreteVelocity(),
//                std::string(
//                    Dune::Parameter::getValue( "micro_reference_solution_filename", std::string( "micro_reference_velocity" ) )
//                        + "_cockburn"
//                        + "_ref_"
//                        + Stuff::toString( referenceSolutionRefineLevel_ ) ) );
//#else
            saveReferenceSolution(
                microSolutionsX.discreteVelocity(),
                std::string(
                    Dune::Parameter::getValue( "micro_reference_solution_filename", std::string( "micro_reference_velocity" ) )
                        + "_X"
                        + "_ref_"
                        + Stuff::toString( referenceSolutionRefineLevel_ ) ) );

//            saveReferenceSolution(
//                microSolutionsY.discreteVelocity(),
//                std::string( Dune::Parameter::getValue( "micro_reference_solution_filename", std::string( "micro_reference_velocity" ) )
//                    + "_Y"
//                    + "_ref_"
//                    + Stuff::toString( referenceSolutionRefineLevel_ ) ) );
//#endif

            info << outputPrefix_ << "Micro reference solution saved." << std::endl;

            debug << outputPrefix_ << "Writing micro output..." << std::endl;

            VTKWriterType microVtkWriter( microGridPart );
            std::string outputFilename = "";

//#if MICRO_PROBLEM==MICRO_COCKBURN_PROBLEM
//            outputFilename = std::string( "data/micro_reference_velocity_cockburn_ref_" ) + Stuff::toString( referenceSolutionRefineLevel_ );
//            microVtkWriter.addVertexData( microSolutions.discreteVelocity() );
//            microVtkWriter.write( outputFilename.c_str() );
//            microVtkWriter.clear();
//
//            microVtkWriter.addVertexData( microSolutions.discretePressure() );
//            outputFilename = std::string( "data/micro_reference_pressure_cockburn_ref_" ) + Stuff::toString( referenceSolutionRefineLevel_ );
//            microVtkWriter.write( outputFilename.c_str() );
//            microVtkWriter.clear();
//#else
            outputFilename = std::string( "data/micro_reference_velocity_X_ref_" ) + Stuff::toString( referenceSolutionRefineLevel_ );
            microVtkWriter.addVertexData( microSolutionsX.discreteVelocity() );
            microVtkWriter.write( outputFilename.c_str() );
            microVtkWriter.clear();

            microVtkWriter.addVertexData( microSolutionsX.discretePressure() );
            outputFilename = std::string( "data/micro_reference_pressure_X_ref_" ) + Stuff::toString( referenceSolutionRefineLevel_ );
            microVtkWriter.write( outputFilename.c_str() );
            microVtkWriter.clear();

//            microVtkWriter.addVertexData( microSolutionsY.discreteVelocity() );
//            outputFilename = std::string( "data/micro_reference_velocity_Y_ref_" ) + Stuff::toString( referenceSolutionRefineLevel_ );
//            microVtkWriter.write( outputFilename.c_str() );
//            microVtkWriter.clear();

//            microVtkWriter.addVertexData( microSolutionsY.discretePressure() );
//            outputFilename = std::string( "data/micro_reference_pressure_Y_ref_" ) + Stuff::toString( referenceSolutionRefineLevel_ );
//            microVtkWriter.write( outputFilename.c_str() );
//            microVtkWriter.clear();
//#endif

            info << outputPrefix_ << "Micro Output written." << std::endl;

//#if MICRO_PROBLEM==MICRO_COCKBURN_PROBLEM
//            info << outputPrefix_ << "COCKBURN reference solution computed." << std::endl;
//#else
            info << outputPrefix_ << "Reference solution computed." << std::endl;
//#endif

        }

        int setVerboseLevel( const int verbose )
        {
            Logging::LogStream& info = Logger().Info();
            Logging::LogStream& debug = Logger().Dbg();
            info.Flush();
            debug.Flush();

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

private:

        int saveReferenceSolution( const DiscreteFunctionType& discreteFunction, const std::string filename )
        {
            Logging::LogStream& error = Logger().Err();
            error.Flush();

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

public:

        int loadReferenceSolution( const std::string filename )
        {
            Logging::LogStream& info = Logger().Info();
            Logging::LogStream& debug = Logger().Dbg();
            Logging::LogStream& error = Logger().Err();
            info.Flush();
            debug.Flush();
            error.Flush();

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
                        break;
                    }
                }
            }
            else {
                ++errorState;
                error << "Error: could not open " << std::string( filename + std::string( ".dgf" ) ) << "!" << std::endl;
                return errorState;
            }

            debug << outputPrefix_ << "RunManager: Loading reference solution..." << std::endl;

            GridPointerType referenceGridPointer( std::string( filename + std::string( ".dgf" ) ) );
            referenceGridPointer->globalRefine( refineLevel );
            GridPartType referenceGridPart( *referenceGridPointer );

            DiscreteFunctionSpaceType referenceSpace( referenceGridPart );

            DiscreteFunctionType referenceSolution( "referenceSolution", referenceSpace );
            const bool readFromAscii = referenceSolution.read_ascii( std::string( filename + std::string( ".ddf" ) ) );
            if ( readFromAscii ) {
                referenceSolution_ = &referenceSolution;
                hasReferenceSolution_ = true;

                VTKWriterType referenceVtkWriter( referenceGridPart );

                std::string outputFilename = std::string( "data/loaded_reference_solution" );
                referenceVtkWriter.addVertexData( referenceSolution );
                referenceVtkWriter.write( outputFilename.c_str() );
                referenceVtkWriter.clear();

                info << outputPrefix_ << "RunManager: Reference solution loaded." << std::endl;
            }
            else {
                ++errorState;
                error << "Error: could not read discrete function from  " << std::string( filename + std::string( ".ddf" ) ) << "!" << std::endl;
                return errorState;
            }
            return errorState;
        }

        double computeDifferenceToReferenceSolution( const DiscreteFunctionType& computedFunction )
        {
            if ( hasReferenceSolution_ && referenceSolution_ ) {

                const DiscreteFunctionType referenceSolution = *referenceSolution_;

                DiscreteFunctionType refinedComputedFunction( "refined_computed_function", referenceSolution.space() );
                refinedComputedFunction.clear();

                typedef Dune::L2Projection< double, double, DiscreteFunctionType, DiscreteFunctionType >
                    ProjectionType;
                ProjectionType projection;
                projection( computedFunction, refinedComputedFunction );

                DiscreteFunctionType error( referenceSolution );
                error -= refinedComputedFunction;

                Dune::L2Norm< GridPartType > l2_Error( referenceSolution.space().gridPart() );
                return l2_Error.norm( error );
            }
            else {
                return -1.0;
            }
        }

    private:

        int referenceSolutionRefineLevel_;
        bool hasReferenceSolution_;
        int verbose_;
        std::string outputPrefix_;
        const std::string dgfFileName_;
        DiscreteFunctionType* referenceSolution_;
};

#endif
