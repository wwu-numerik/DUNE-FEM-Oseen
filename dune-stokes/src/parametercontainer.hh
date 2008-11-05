/** \file parametercontainer.hh
    \brief  containing class ParameterContainer
 **/

#ifndef PARAMETERCONTAINER_HH_INCLUDED
#define PARAMETERCONTAINER_HH_INCLUDED

#include "dune/fem/io/parameter.hh"

#include "stuff.hh"
#include "logging.hh"


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
         *  \attention  call ReadCommandLine to set up parameterContainer
         **/
        ParameterContainer()
        {
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
       bool ReadCommandLine( int argc, char** argv )
        {
            if ( argc_ == 2 )
            {
                argc_ = argc;
                argv_ = argv;
                eps_ = 1.0e20;
                parameter_filename_ = argv_[1];
                return true;
            }
            else
            {
                std::cerr << "\nUsage: " << argv_[0] << " parameterfile" << std::endl;
                PrintParameterSpecs( std::cerr );
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
            if ( !( Dune::Parameter::exists( "grid_dimension" ) ) ) {
                std::cerr << "\nError: not all parameters found in " << parameterFilename() << "!";
                PrintParameterSpecs( std::cerr );
                std::cerr << "missing parameters are: grid_dimension" << std::endl;
                has_not_worked = true;

            }
            else {
                Dune::Parameter::get( "grid_dimension", grid_dimension_ );
            }
            if ( !( Dune::Parameter::exists( "polynomial_order" ) ) ) {
                if ( !( has_not_worked ) ) {
                    std::cerr << "\nError: not all parameters found in " << parameterFilename() << "!";
                    PrintParameterSpecs( std::cerr );
                    std::cerr << "missing parameters are: polynomial_order" << std::endl;
                }
                else {
                    std::cerr << "                        polynomial_order" << std::endl;
                }
                has_not_worked = true;
            }
            else {
                Dune::Parameter::get( "polynomial_order", pol_order_ );
            }
            if ( has_not_worked ) {
                std::cerr << "\nError: not all parameters found in " << parameterFilename() << "!";
                PrintParameterSpecs( std::cerr );
                std::cerr << std::endl;
            }
            return !( has_not_worked );
        }

        /**
         *  \brief prints on a given stream, how a parameterfile schould look like
         *  \arg std::ostream& out where to print
         **/
        void PrintParameterSpecs( std::ostream& out )
        {
            out << "\na valid parameterfile should at least specify the following parameters:";
            out << "\n(copy this into your parameterfile)" << std::endl;
            //out << "grid_dimension: " << std::endl;
            //out << "polynomial_order: " << std::endl;
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

        /**
         *  \brief  sets the filename of the parameterfile
         **/
        void SetParameterFilename( const std::string param_filename )
        {
            parameter_filename_ = param_filename;
        }

        /**
         *  \brief  sets the dimension
         **/
        void SetGridDimension( const int grid_dim )
        {
            grid_dimension_ = grid_dim;
        }

        /**
         *  \brief  returns the dimension
         *  \return int dimension
         **/
        int gridDimension() const
        {
            return grid_dimension_;
        }

        /**
         *  \brief  returns the polynomial order
         *  \return int polynomial order
         **/
        int polOrder() const
        {
            return pol_order_;
        }

        /**
         *  \brief  sets the polynomial order
         **/
        void SetPolOrder( const int pol_order )
        {
            pol_order_ = pol_order;
        }

    private:
        int grid_dimension_;
        int pol_order_;
        double eps_;
        int argc_;
        char** argv_;
        bool all_set_up_;
        std::string parameter_filename_;
};

///global ParameterContainer instance
ParameterContainer& Parameters()
{
    static ParameterContainer parameters;
    return parameters;
}
#endif // end of PARAMETERHANDLER.HH
