/** \file parameterhandler.hh
    \brief  containing class ParameterHandler
 **/

#ifndef PARAMETERHANDLER_HH
#define PARAMETERHANDLER_HH

#include <fstream>
#include <map>
#include "stuff.hh"

/**
 *  \brief class processing parameter file
 *
 *  \c ParameterHandler reads a parameter file once and stores all found values internally
 **/
class ParameterHandler
{
    private:
        typedef std::map< std::string, std::string > MapType;
        MapType parameter_map_;
        bool status_;

    public:
        /** \ Please doc me!
         *
         **/
        ParameterHandler( const std::string filename )
            :status_( false )
        {
            std::ifstream parameter_file( filename.c_str() );
            if( parameter_file.is_open() )
            {
                while( !parameter_file.eof() )
                {
                    std :: string line;
                    std :: getline( parameter_file, line );

                    if( line.size() == 0 )
                        continue;

                    if( (line[ 0 ] != '%') && (line[ 0 ] != '#') )
                    {
                        std::string param_name = line.substr( 0, line.find_first_of(":") );
                        std::string value = line.substr( line.find_first_of(":")+1, line.length() );
                        parameter_map_[param_name] = value;
                    }

                }
                status_ = true;
                parameter_file.close();
            }
            else {
                //LOGERROR
                status_ = false;
                std::cerr << "ERROR: file " << filename << " not found!\n";
            }
        }

        /** \todo Please doc me! */
        template < class ReturnType >
        ReturnType GetParameter(const std::string name) const{
            MapType::const_iterator it = parameter_map_.find( name ) ;
            if ( it != parameter_map_.end() ){
                return Stuff::fromString<ReturnType>( it->second );
            }
            else {
                //LogError
                return (ReturnType) 0;
            }

        }

        /** \todo Please doc me! */
        void Print( std::ostream &out ) const
        {
            for (MapType::const_iterator it = parameter_map_.begin(); parameter_map_.end() != it; ++it){
                out << it->first << ":" << it->second << "\n" ;
            }
            out << std::endl;
        }

        /** \todo Please doc me! */
        bool Ok()
        {
            return status_;
        }

        /** \todo Please doc me! */
        ~ParameterHandler(){}

};

/**
 *  \brief class containing global parameters
 *
 *  \c ParameterContainer contains all the needed global parameters
 **/
class ParameterContainer
{
    private:

    public:
        /**
         *  \brief std constuctor
         *  \todo   implement + doc me
         **/
        ParameterContainer(){}

        /**
         *  \brief std destuctor
         *  \todo   implement + doc me
         **/
        ~ParameterContainer(){}
};

#endif // end of PARAMETERHANDLER.HH
