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
 *
 *  \todo needs reorganizing!
 **/
class ParameterContainer
{
    public:
        /**
         *  \brief  constuctor
         *  \attention  call ReadCommandLine to set up parameterContainer
         **/
        ParameterContainer()
            : all_set_up_( false )
        {
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
         *  \arg    std::ostream& out stream to print to
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
            if ( argc == 2 )
            {
                argc_ = argc;
                argv_ = argv;
                eps_ = 1.0e20;
                parameter_filename_ = argv[1];
                return true;
            }
            else
            {
                std::cerr << "\nUsage: " << argv[0] << " parameterfile" << std::endl;
                PrintParameterSpecs( std::cerr );
                std::cerr << std::endl;
                return false;
            }
        }

        /**
         *  \brief  sets all needed parameters
         *  \return bool returns true, if all needed parameters are set up
         *  \attention  SetGridDimension should be called next
         **/
        bool SetUp()
        {
            Dune::Parameter::append( ParameterFilename() );
            bool has_not_worked = false;
            if ( !( Dune::Parameter::exists( "dgf_file_2d" ) ) ) {
                std::cerr << "\nError: not all parameters found in " << ParameterFilename() << "!";
                PrintParameterSpecs( std::cerr );
                std::cerr << "\nmissing parameters are: dgf_file_2d" << std::endl;
                has_not_worked = true;

            }
            else {
                Dune::Parameter::get( "dgf_file_2d", dgf_filenames_[1] );
            }
            if ( !( Dune::Parameter::exists( "dgf_file_3d" ) ) ) {
                if ( !( has_not_worked ) ) {
                    std::cerr << "\nError: not all parameters found in " << ParameterFilename() << "!";
                    PrintParameterSpecs( std::cerr );
                    std::cerr << "\nmissing parameters are: dgf_file_3d" << std::endl;
                }
                else {
                    std::cerr << "                        dgf_file_3d\n" << std::endl;
                }
                has_not_worked = true;
            }
            else {
                Dune::Parameter::get( "dgf_file_3d", dgf_filenames_[2] );
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
            out << "\nRemark: the correpondig files have to exist!" << std::endl;
            out << "\n(copy this into your parameterfile)" << std::endl;
            out << "dgf_file_2d: " << std::endl;
            out << "dgf_file_3d: " << std::endl;
        }

        /**
         *  \brief  returns the filename of the parameterfile
         *  \return std::string filename of the parameterfile
         **/
        std::string ParameterFilename() const
        {
            return parameter_filename_;
        }

        /**
         *  \brief  returns the filename of the dgf file
         *  \return std::string filename of the dgf file
         **/
        std::string DgfFilename() const
        {
            return dgf_filenames_[ grid_dimension_ - 1 ];
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
        int GridDimension() const
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
        std::string dgf_filenames_[3];
};

///global ParameterContainer instance
ParameterContainer& Parameters()
{
    static ParameterContainer parameters;
    return parameters;
}
#endif // end of PARAMETERHANDLER.HH
