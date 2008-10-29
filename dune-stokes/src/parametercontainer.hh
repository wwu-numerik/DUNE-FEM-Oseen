/** \file parametercontainer.hh
    \brief  containing class ParameterContainer
 **/

#ifndef PARAMETERHANDLER_HH
#define PARAMETERHANDLER_HH

#include "dune/fem/io/parameter.hh"
//#include <fstream>
//#include <ostream>
//#include <map>
//#include "stuff.hh"
//#include "logging.hh"

///**
// *  \brief class processing parameter file
// *
// *  \c ParameterHandler reads a parameter file once and stores all found values internally
// **/
//class ParameterHandler
//{
//    private:
//        typedef std::map< std::string, std::string > MapType;
//        MapType parameter_map_;
//        bool status_;
//
//    public:
//        /** \ Please doc me!
//         *
//         **/
//        ParameterHandler( const std::string filename )
//            :status_( false )
//        {
//            ParseParamFile( filename );
//        }
//
//        ParameterHandler(  )
//            :status_( false )
//        {
//        }
//
//        /** \brief function used for parametrized ctor and two-step creation
//        **/
//        bool ParseParamFile( const std::string filename )
//        {
//            std::ifstream parameter_file( filename.c_str() );
//            if( parameter_file.is_open() )
//            {
//                while( !parameter_file.eof() )
//                {
//                    std :: string line;
//                    std :: getline( parameter_file, line );
//
//                    if( line.size() == 0 )
//                        continue;
//
//                    if( (line[ 0 ] != '%') && (line[ 0 ] != '#') )
//                    {
//                        std::string param_name = line.substr( 0, line.find_first_of(":") );
//                        std::string value = line.substr( line.find_first_of(":")+1, line.length() );
//                        parameter_map_[param_name] = value;
//                    }
//
//                }
//                status_ = true;
//                parameter_file.close();
//            }
//            else {
//                //LOGERROR
//                status_ = false;
//                std::cerr << "ERROR: file " << filename << " not found!\n";
//            }
//            return Ok();
//        }
//
//        /**
//         *  \brief checks, if a parameter is found in the parameterfile
//         *  \arg const std::string name name of the parameter to be found
//         *  \return returns true, if the parameter is found
//         **/
//         bool ParameterExists( const std::string name ) const
//         {
//             if ( !( status_ ) ) {
//                 return false;
//             }
//             else {
//                 MapType::const_iterator it = parameter_map_.find( name );
//                 if ( it != parameter_map_.end() ) {
//                     return true;
//                 }
//                 else {
//                     return false;
//                 }
//
//             }
//         }
//
//        /** \todo Please doc me! */
//        template < class ReturnType >
//        ReturnType GetParameter( const std::string name) const
//        {
//            if ( !( status_ ) ) {
//                return (ReturnType) 0;
//            }
//            else {
//                MapType::const_iterator it = parameter_map_.find( name ) ;
//                if ( it != parameter_map_.end() ){
//                    return Stuff::fromString<ReturnType>( it->second );
//                }
//                else {
//                    //LogError
//                    return (ReturnType) 0;
//                }
//            }
//        }
//
//        /** \todo Please doc me! */
//        void Print( LogStream &out ) const
//        {
//            assert( status_ );
//            for (MapType::const_iterator it = parameter_map_.begin(); parameter_map_.end() != it; ++it){
//                out << it->first << ":" << it->second << "\n" ;
//            }
//            //out << std::endl;
//        }
//
//        /** \todo Please doc me! */
//        bool Ok()
//        {
//            return status_;
//        }
//
//        /** \todo Please doc me! */
//        ~ParameterHandler(){}
//
//};

///** \brief global singelton for paramhandler
//**/
//ParameterHandler& params()
//{
//    static ParameterHandler param;
//    return param;
//}

/**
 *  \brief class containing global parameters
 *
 *  \c ParameterContainer contains all the needed global parameters getting them via Dune::Parameter
 **/
class ParameterContainer
{
    public:
        /**
         *  \brief  constuctor
         *  \arg    int argc   number of command line arguments
         *  \arg    char** argv array of command line arguments
         **/
        ParameterContainer( int argc, char** argv )
        {
            all_set_up_ = false;
            argc_ = argc;
            argv_ = argv;
            eps_ = 1.0e20;
            all_set_up_ = false;
        }

        /**
         *  \brief  destuctor
         *  \todo   implement + doc me
         **/
        ~ParameterContainer()
        {
        }

        /**
         *  \brief  prints all parameters
         *  \todo   implement me
         *  \arg    std::ostream& out stream to print out
         **/
        void Print( std::ostream& out ) const
        {
            out << "\nthis is the ParameterContainer.Print() function" << std::endl;
        }

         /**
         *  \brief  checks command line parameters
         *  \return bool returns true, if comman line arguments are valid
         **/
       bool ReadCommandLine()
        {
            if ( argc_ == 2 )
            {
                parameter_filename_ = argv_[1];
                return true;
            }
            else
            {
                std::cerr << "\nUsage: " << argv_[0] << " parameterfile\n" << std::endl;
                return false;
            }
        }

        /**
         *  \brief  sets all needed parameters
         *  \return bool returns true, if all needed parameters are set up
         **/
        bool SetUp()
        {
            Dune::Parameter::append( parameterFilename() );
            Dune::Parameter::append( argc_, argv_ );
            bool has_not_worked = false;
            if ( !( Dune::Parameter::exists( "dimension" ) ) )
            {
                has_not_worked = ( has_not_worked || true );
            }
            else {
                Dune::Parameter::get( "dimension", dimension_ );
            }
            if ( !( Dune::Parameter::exists( "polynomial_order" ) ) )
            {
                has_not_worked = ( has_not_worked || true );
            }
            else {
                Dune::Parameter::get( "polynomial_order", pol_order_ );
            }
            if ( has_not_worked ) {
            std::cerr << "\nError: not all parameters found in " << parameterFilename() << "!";
            PrintParameterSpecs( std::cerr );
            }
            return !( has_not_worked );
        }

        /**
         *  \brief prints on a given stream, how a parameterfile schould look like
         *  \arg std::ostream& out where to print
         **/
        void PrintParameterSpecs( std::ostream& out ) const
        {
            out << "\na valid parameterfile should at least specify the following parameters:" << std::endl;
            out << "dimension: " << std::endl;
            out << "polynomial_order: " << std::endl;
            out << std::endl;
        }

        /**
         *  \brief  returns the filename of the parameterfile
         *  \return std::string filename of the parameterfile
         **/
        std::string parameterFilename() const
        {
            return parameter_filename_;
        }

    private:
        int dimension_;
        int pol_order_;
        double eps_;
        int argc_;
        char** argv_;
        bool all_set_up_;
        std::string parameter_filename_;
};

#endif // end of PARAMETERHANDLER.HH
