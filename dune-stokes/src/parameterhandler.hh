/** \file parameterhandler.hh
    \brief  brief
 **/

#ifndef PARAMETERHANDLER_HH
#define PARAMETERHANDLER_HH

#include <fstream>
#include <map>
#include "stuff.hh"

/** \todo Please doc me! */
class ParameterHandler
{
    private:
        typedef std::map< std::string, std::string > MapType;
        MapType parameter_map_;
        bool status_;

    public:
        /** \todo Please doc me!
         *  private function has_worked() fehlt
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
        void Print( ) const
        {
            for (MapType::const_iterator it = parameter_map_.begin(); parameter_map_.end() != it; ++it){
                std::cout << it->first << ":" << it->second << "\n" ;
            }
            std::cout << std::endl;
        }

        /** \todo Please doc me! */
        bool Ok()
        {
            return status_;
        }

        /** \todo Please doc me! */
        ~ParameterHandler(){}

};

#endif // end of PARAMETERHANDLER.HH
