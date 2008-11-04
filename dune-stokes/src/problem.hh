/** \file problem.hh
    \brief problem.hh
 **/

#ifndef PROBLEM_HH
#define PROBLEM_HH

#include <cmath>
#include "dune/common/fvector.hh"

/**
 *  \brief describes the velocity
 *  \f[u(x_{1},x_{2}):=-e^{x_{1}}\Big( x_{2} cos(x_{2}) + sin(x_{2}) \Big)\f]
 *
 *  and the gradient of the velocity
 *
 *  as the solution of the stokes test problem
 *  \todo doc test problem
 **/
template < int grid_dim >    // class C should be Dune::FieldVector
class Velocity
{
    public:
        enum {
            grid_dim = GridType::dimensionworld
        };
        typedef Dune::FieldVector< double, grid_dim > DomainType;
        typedef Dune::FieldVector< double, grid_dim > RangeType;
        typedef Dune::FieldVector< RangeType, grid_dim > GradientRangeType;
        typedef Dune::FieldVector< double, 1 > DivergenceRangeType;
        /**
         *  \brief constructor
         **/
        Velocity()
        {
        }

        /**
         *  \brief  destructor
         **/
        ~Velocity()
        {
        }

        /**
         *  \brief evaluates the velocity
         *  \arg const C& arg
         *  \arg C& ret returns velocity at arg
         *  \return true if has worked
         **/
        bool evaluate( const DomainType& arg, DomainType& ret )
        {
            if ( grid_dim == 2 ) {
                if ( arg.dim() == grid_dim ) {
                    // some computations
                    double x1 = arg[0];
                    double x2 = arg[1];
                    double exp_of_x1 = std::exp( x1 );
                    double sin_of_x2 = std::sin( x2 );
                    // return
                    ret[0] = -1.0 * exp_of_x1 * ( x2 * std::cos( x2 ) + sin_of_x2 );
                    ret[1] = exp_of_x1 * x2 * sin_of_x2;
                    return true;
                }
                else {
                    std::cerr << "\nError in Velocity.evaluate(): "
                        << "dimension of argument does not match dimension "
                        << "of grid!" << std::endl;
                    return false;
                }
            }
            else {
                std::cerr << "\nError in Velocity.gradient(): "
                    << "dimension of grid is not yet implemented!" << std::endl;
                return false;
            }
        }

        /**
         *  \brief evaluates the gradient of the velocity
         **/
        bool gradient( const DomainType& arg, GradientRangeType& ret )
        {
            if ( grid_dim == 2 ) {
                if ( arg.dim() == grid_dim ) {
                    // some computations
                    double x1 = arg[0];
                    double x2 = arg[1];
                    double exp_of_x1 = std::exp( x1 );
                    double sin_of_x2 = std::sin( x2 );
                    double cos_of_x2 = std::cos( x2 );
                    // gradient of u_{1}
                    double du1_over_dx1 = -1.0 * exp_of_x1 *
                        ( x2 * cos_of_x2 + sin_of_x2 );
                    double du1_over_dx2 = exp_of_x1 *
                        ( x2 * sin_of_x2 - 2.0 * cos_of_x2 );
                    RangeType grad_u1;
                    grad_u1[0] = du1_over_dx1;
                    grad_u1[1] = du1_over_dx2;
                    // gradient of u_{2}
                    double du2_over_dx1 = exp_of_x1 * x2 * sin_of_x2;
                    double du2_over_dx2 = exp_of_x1 *
                        ( x2 * cos_of_x2 + sin_of_x2 );
                    RangeType grad_u2;
                    grad_u2[0] = du2_over_dx1;
                    grad_u2[1] = du2_over_dx2;
                    // return
                    ret[0] = grad_u1;
                    ret[1] = grad_u2;
                    return true;
                }
                else {
                    std::cerr << "\nError in Velocity.gradient(): "
                        << "dimension of argument "
                        << arg.dim() << " does not match dimension "
                        << "of grid " << grid_dim << "!"
                        << std::endl;
                    return false;
                }

                return true;
            }
            else {
                std::cerr << "\nError in Velocity.gradient(): "
                    << "dimension "
                    << grid_dim << " is not yet implemented!" << std::endl;
                return false;
            }
        }

    private:



};

#endif  // end of problem.hh

