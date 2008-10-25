/**************************************************************************
**       Title: examples of elliptic models 
**    $RCSfile$
**   $Revision: 1867 $$Name$
**       $Date: 2007-07-19 11:05:24 +0200 (Thu, 19 Jul 2007) $
**   Copyright: GPL $Author: nolte $
** Description: implementation of default elliptic models for 
**              use in solvers such as FEOp. A model is a class providing
**              all data functions of a problem. The model is assumed to
**              be instantiated once and passed to solvers. It is depending
**              on a traits class providing typedefs.
**
**************************************************************************/
#ifndef DUNE_ELLIPTICMODEL_HH
#define DUNE_ELLIPTICMODEL_HH

/** \file
    \brief mystokesmodell.hh implementation of default elliptic models
 */

#include <sstream>
#include <dune/fem/operator/model/linearellipticmodel.hh>
#include <dune/fem/space/common/functionspace.hh>
#include <dune/fem/function/common/function.hh>
#include <dune/fem/misc/boundaryidentifier.hh>
namespace Dune
{
  
  /** \class PoissonModel
   *  \brief The PoissonModel class provides a default model for an elliptic
   *  problem to be handled by FEOp
   *
   *  The model problem simply is
   *  \f[ - \mathrm{div} \nabla u = n \pi^2 \prod_{i=1}^n \sin( \pi x_i ). \f]
   *  Using homogeneous Dirichlet boundary values, the exact solution on the unit
   *  square is \f$u( x ) = \prod_{i=1}^n \sin( \pi x_i ) \f$
   *
   *  All types are extracted from the TraisImp, which defaults to
   *  DefaultElementMatrixTraits.
   */
  
  //  static double nu=1.;
  static double shift=1.001;

template <typename Field, class GridImp, int dimR>
struct ModelParam
{
  enum { dim      = GridImp::dimension };
  enum { dimworld = GridImp::dimensionworld };
  
  enum { dimRange  = dimR };
  enum { dimDomain = dimworld };

  typedef Field FieldType;
  
  typedef GridImp GridType;
};
  /*
template <class Grid>
class ModelTraits {
  public:
    enum { dimRange2 = 1 };
    enum { dimRange1= dimRange2*Grid::dimensionworld };
    typedef Grid GridType;
    enum { dimDomain = GridType::dimensionworld };
    enum { dimRange = dimRange2 }; 
    enum { dimGradRange = dimRange1 };
    typedef FieldVector<double, dimDomain> DomainType;
    typedef FieldVector<double, dimDomain-1> FaceDomainType;
    typedef FieldVector<double,dimRange> RangeType;
    typedef FieldVector<double,dimGradRange> GradientType;
    typedef FieldMatrix<double,dimDomain+1,dimDomain> FluxRangeType;
    typedef FieldMatrix<double,dimGradRange,dimDomain> DiffusionRangeType;
    typedef typename GridType::template Codim<0>::Entity EntityType;
   };
  */
  
  template< class ModelParamType,class FunctionSpaceImp ,class PressureSpaceImp>
  class StokesModel
    : public LinearEllipticModelDefault< FunctionSpaceImp, StokesModel< ModelParamType,FunctionSpaceImp,PressureSpaceImp > >
  {
  public:  
    typedef FunctionSpaceImp FunctionSpaceType;
    typedef FunctionSpaceType FuncSpaceType;
    typedef PressureSpaceImp PressureSpaceType;
  public:
    enum {dim      = ModelParamType :: dim };
    enum {dimworld = ModelParamType :: dimworld };
    
    enum {dimRange  = ModelParamType :: dimRange };
    enum {dimDomain = ModelParamType :: dimDomain };
    enum {dimGradRange = dimRange * dimDomain }; 

    //! boundary types that may be used by the model
    enum BndType { Dirichlet, Neumann, NoFlow, OutFlow };
    
    typedef typename ModelParamType :: GridType          GridType; 
    typedef typename ModelParamType :: FieldType         FieldType; 
    
    // typedef FunctionSpace < FieldType , FieldType, dim , dimRange > FuncSpaceType;
    //   typedef FunctionSpace < FieldType , FieldType, dim , dimRange > FunctionSpaceType;

  
    typedef MatrixFunctionSpace < FieldType , FieldType, dimDomain, dimRange,dimDomain > GradFuncSpaceType;

    typedef typename FuncSpaceType :: RangeFieldType RangeFieldType;
  
    typedef typename FuncSpaceType::RangeType                RangeType;
    typedef typename FuncSpaceType::DomainType               DomainType;
    typedef typename FuncSpaceType::JacobianRangeType        JacobianRangeType;
  
    typedef typename GradFuncSpaceType::RangeType            GradRangeType;
    typedef typename GradFuncSpaceType::DomainType           GradDomainType;
    typedef typename GradFuncSpaceType::JacobianRangeType    GradJacobianRangeType;
  
    typedef typename PressureSpaceType::RangeType PressureRangeType;
    
    typedef FieldMatrix<double,dimGradRange,dimGradRange> DiffusionRangeType;
    typedef BoundaryIdentifier BoundaryIdentifierType ;
    typedef typename GridType::template Codim<0>::Entity Entity;
  public:
    struct Traits
    {
      typedef typename ModelParamType :: GridType GridType;
      enum {dimRange = ModelParamType :: dimRange };
      enum {dimDomain = dim };
      enum {dimGradRange = dimDomain * dimRange };
      typedef FieldVector< FieldType, dimDomain-1 > FaceDomainType;
      typedef typename GridType::template Codim<0>::Entity Entity;
    
      typedef FunctionSpace < FieldType , FieldType, dim , dimRange > FuncSpaceType;
      
      typedef MatrixFunctionSpace < FieldType , FieldType, dimDomain, dimRange,dimDomain > GradFuncSpaceType;
      
      typedef typename FuncSpaceType::RangeType          RangeType;
      typedef typename FuncSpaceType::DomainType         DomainType;
      typedef typename FuncSpaceType::JacobianRangeType  JacobianRangeType;
      typedef FieldMatrix<double,dimDomain+1,dimDomain> FluxRangeType;
    };



  private:
    typedef StokesModel< ModelParamType,FunctionSpaceType,PressureSpaceType > ThisType;
    typedef LinearEllipticModelDefault< FunctionSpaceType, ThisType > BaseType;
    
  public:
    typedef typename BaseType :: BoundaryType BoundaryType;

 //    typedef typename FunctionSpaceType :: DomainType DomainType;
    // typedef typename FunctionSpaceType :: RangeType RangeType;
//     typedef typename FunctionSpaceType :: JacobianRangeType JacobianRangeType;

//     typedef typename FunctionSpaceType :: DomainFieldType DomainFieldType;
//     typedef typename FunctionSpaceType :: RangeFieldType RangeFieldType;

  public:
    //! return boundary type of a boundary point p used in a quadrature
    template< class IntersectionIteratorType >
    inline BoundaryType boundaryType( const IntersectionIteratorType &intersection ) const
    {
        return BaseType :: Dirichlet;
    }

    //! determine dirichlet value in a boundary point used in a quadrature
    template< class IntersectionIteratorType, class QuadratureType >
    inline void dirichletValues( const IntersectionIteratorType &intersection,
                                 const QuadratureType &quadrature,
                                 int p,
                                 RangeType &ret ) const
    {
      typedef typename IntersectionIteratorType :: Entity EntityType;
      
      const int dimension = DomainType :: dimension;
      
      const DomainType &x = intersection.inside()->geometry().global( quadrature.point( p ) );
      //u1
      ret[0]=cos(x[1]);
      ret[0]*=x[1];
      ret[0]+=sin(x[1]);
      ret[0]*=exp(x[0]);
      ret[0]*=-1.0;
     
      //u2
      ret[1]=sin(x[1]);
      ret[1]*=x[1];
      ret[1]*=exp(x[0]);

    }

    //! determine neumann value in a boundary point used in a quadrature
    template< class IntersectionIteratorType, class QuadratureType >
    inline void neumannValues( const IntersectionIteratorType &intersection,
                               const QuadratureType &quadrature,
                               int p,
                               RangeType &ret ) const
    {
      std :: cout << "Neumann boundary values are not implemented." << std :: endl;
      assert( false );

      //const DomainType &x = entity.geometry().global( quadrature.point( p ) ); 
            
      ret[ 0 ] = 0.0;
    }

    //! determine robin value in a boundary point used in a quadrature
    template< class IntersectionIteratorType, class QuadratureType >
    inline void robinValues( const IntersectionIteratorType &intersection,
                             const QuadratureType &quadrature,
                             int p, 
                             RangeType &ret ) const
    {
      std :: cout << "Robin boundary values are not implemented." << std :: endl;
      assert( false );

      //const DomainType &x = entity.geometry().global( quadrature.point( p ) ); 
            
      ret[ 0 ] = 0.0;
    }

    //! determine mass (i.e. value of the function c) in a quadrature point
    template< class EntityType, class QuadratureType >
    inline void mass( const EntityType &entity,
                      const QuadratureType &quadrature,
                      int p,
                      RangeType &ret ) const
    {

      assert(false);
      //const DomainType &x = entity.geometry().global( quadrature.point( p ) ); 

      ret[ 0 ] = 0.0;
    }

    //! Determine source (i.e. value of the function f) in a quadrature point
    template< class EntityType, class QuadratureType >
    inline void source( const EntityType &entity,
                        const QuadratureType &quadrature,
                        int p,
                        RangeType &ret ) const
    {
      const int dimension = DomainType :: dimension;
      
      const DomainType &x = entity.geometry().global( quadrature.point( p ) );
      
      double factor=1-nu_;
      
      ret[0]=sin(x[1]);
      ret[0]*=2;
      ret[0]*=exp(x[0]);
      ret[0]*=factor;
      ret[1]=cos(x[1]);
      ret[1]*=2;
      ret[1]*=exp(x[0]);
      ret[1]*=factor;
      //    return ret; 

    
    }

    //! no direct access to stiffness and velocity, but whole flux, i.e.
    //! diffflux = stiffness * grad( phi )
    template< class EntityType, class QuadratureType >
    inline void diffusiveFlux( const EntityType &entity,
                               const QuadratureType &quadrature,
                               int p,
                               const JacobianRangeType &gradphi, 
                               JacobianRangeType &ret ) const
    {
      assert(false);
      // - laplace phi = -div (1 * grad phi - 0 * phi)
      ret = gradphi;          
    }

    //! no direct access to stiffness and velocity, but whole flux, i.e.
    //! convectiveFlux =  - velocity * phi 
    template <class EntityType, class QuadratureType>  
    inline void convectiveFlux( const EntityType &entity,
                                const QuadratureType &quadrature,
                                int p,
                                const RangeType &phi, 
                                DomainType &ret ) const
    {
      assert(false);
      // - laplace phi = -div (1 * grad phi - 0 * phi)
      ret = 0.0;
    }
    
    template <class EntityType, class QuadratureType>  
    inline void convectiveFlux( const EntityType &entity,
                                const QuadratureType &quadrature,
                                int p,
                                const RangeType &phi, 
                                JacobianRangeType &ret ) const
    {
      // - laplace phi = -div (1 * grad phi - 0 * phi)
      ret = 0.0;
    }





    //! the coefficient for robin boundary condition
    template< class IntersectionIteratorType, class QuadratureType >
    inline double robinAlpha( const IntersectionIteratorType &intersection,
                              const QuadratureType &quadrature,
                              int p ) const
    { 
      return 1.0;
    }

    void setnu(double nu){nu_=nu;}
    
    double nu_;
  };  // end of PoissonModel class
  /*======================================================================*/
  /*!
   *  \class PoissonExactSolution
   *  \brief The class provides the exact solution for the model given by 
   *         the PoissonModel class
   *
   *  The function represents u = x(1-x)y(1-y), which is the solution of
   *  the model problem  - div grad u = 2(x(1-x)+ y(1-y)) on
   *  the unit square with homogeneous Dirichlet boundary values. 
   *  The function can be used for EOC calculation
   */
  /*======================================================================*/

  template< class FunctionSpaceImp >
  class StokesExactSolution
    : public Function< FunctionSpaceImp, StokesExactSolution< FunctionSpaceImp > >
  {
  public:
    typedef FunctionSpaceImp FunctionSpaceType;

  private:
    typedef StokesExactSolution< FunctionSpaceType > ThisType;
    typedef Function< FunctionSpaceType, ThisType > BaseType;

  public:
    typedef typename FunctionSpaceType :: DomainType DomainType;
    typedef typename FunctionSpaceType :: RangeType RangeType;
    typedef typename FunctionSpaceType :: JacobianRangeType JacobianRangeType;

    typedef typename FunctionSpaceType :: DomainFieldType DomainFieldType;
    typedef typename FunctionSpaceType :: RangeFieldType RangeFieldType;

  public:
    inline StokesExactSolution ( const FunctionSpaceType &functionSpace )
    : BaseType( functionSpace )
    {
    }

 
    //! u( x ) = sin( pi x_1 ) * ... * sin( pi x_n )
    inline void evaluate ( const DomainType &x, RangeType &ret ) const
    {
      enum { dimension = DomainType :: dimension };
    
    
      //u1
      ret[0]=cos(x[1]);
      ret[0]*=x[1];
      ret[0]+=sin(x[1]);
      ret[0]*=exp(x[0]);
      ret[0]*=-1.0;
     
      //u2
      ret[1]=sin(x[1]);
      ret[1]*=x[1];
      ret[1]*=exp(x[0]);

      
      // ret[0]+=1.0;
    }

    inline void evaluate( const DomainType &x, RangeFieldType t, RangeType &ret ) const
    {
      evaluate( x , ret );
    }
  };
/*======================================================================*/
  /*!
   *  \class PoissonExactGradient
   */
  /*======================================================================*/
 template< class FunctionSpaceImp >
  class StokesExactPressure
    : public Function< FunctionSpaceImp, StokesExactPressure< FunctionSpaceImp > >
  {
  public:
    typedef FunctionSpaceImp FunctionSpaceType;

  private:
    typedef StokesExactPressure< FunctionSpaceType > ThisType;
    typedef Function< FunctionSpaceType, ThisType > BaseType;

  public:
    typedef typename FunctionSpaceType :: DomainType DomainType;
    typedef typename FunctionSpaceType :: RangeType RangeType;
    typedef typename FunctionSpaceType :: JacobianRangeType JacobianRangeType;

    typedef typename FunctionSpaceType :: DomainFieldType DomainFieldType;
    typedef typename FunctionSpaceType :: RangeFieldType RangeFieldType;

  public:
    inline StokesExactPressure ( const FunctionSpaceType &functionSpace )
    : BaseType( functionSpace )
    {
    }

 
    //! u( x ) = sin( pi x_1 ) * ... * sin( pi x_n )
    inline void evaluate ( const DomainType &x, RangeType &ret ) const
    {
      enum { dimension = DomainType :: dimension };
       
 
      ret[0] = sin(x[1]);
      ret[0] *=exp(x[0]);
      ret[0] *=2.0;

    }

    inline void evaluate( const DomainType &x, RangeFieldType t, RangeType &ret ) const
    {
      evaluate( x , ret );
    }
  };

 

 

  template< class ModelParamType,class FunctionSpaceImp ,class PressureSpaceImp>
  class StokesSingModel
    : public LinearEllipticModelDefault< FunctionSpaceImp, StokesSingModel< ModelParamType,FunctionSpaceImp,PressureSpaceImp > >
  {
  public:  
    typedef FunctionSpaceImp FunctionSpaceType;
    typedef FunctionSpaceType FuncSpaceType;
    typedef PressureSpaceImp PressureSpaceType;
  public:
    enum {dim      = ModelParamType :: dim };
    enum {dimworld = ModelParamType :: dimworld };
    
    enum {dimRange  = ModelParamType :: dimRange };
    enum {dimDomain = ModelParamType :: dimDomain };
    enum {dimGradRange = dimRange * dimDomain }; 

    //! boundary types that may be used by the model
    enum BndType { Dirichlet, Neumann, NoFlow, OutFlow };
    
    typedef typename ModelParamType :: GridType          GridType; 
    typedef typename ModelParamType :: FieldType         FieldType; 
    
    // typedef FunctionSpace < FieldType , FieldType, dim , dimRange > FuncSpaceType;
    //   typedef FunctionSpace < FieldType , FieldType, dim , dimRange > FunctionSpaceType;

  
    typedef MatrixFunctionSpace < FieldType , FieldType, dimDomain, dimRange,dimDomain > GradFuncSpaceType;

    typedef typename FuncSpaceType :: RangeFieldType RangeFieldType;
  
    typedef typename FuncSpaceType::RangeType                RangeType;
    typedef typename FuncSpaceType::DomainType               DomainType;
    typedef typename FuncSpaceType::JacobianRangeType        JacobianRangeType;
  
    typedef typename GradFuncSpaceType::RangeType            GradRangeType;
    typedef typename GradFuncSpaceType::DomainType           GradDomainType;
    typedef typename GradFuncSpaceType::JacobianRangeType    GradJacobianRangeType;
  
    typedef typename PressureSpaceType::RangeType PressureRangeType;
    
    typedef FieldMatrix<double,dimGradRange,dimGradRange> DiffusionRangeType;
    typedef BoundaryIdentifier BoundaryIdentifierType ;
    typedef typename GridType::template Codim<0>::Entity Entity;
  public:
    struct Traits
    {
      typedef typename ModelParamType :: GridType GridType;
      enum {dimRange = ModelParamType :: dimRange };
      enum {dimDomain = dim };
      enum {dimGradRange = dimDomain * dimRange };
      typedef FieldVector< FieldType, dimDomain-1 > FaceDomainType;
      typedef typename GridType::template Codim<0>::Entity Entity;
    
      typedef FunctionSpace < FieldType , FieldType, dim , dimRange > FuncSpaceType;
      
      typedef MatrixFunctionSpace < FieldType , FieldType, dimDomain, dimRange,dimDomain > GradFuncSpaceType;
      
      typedef typename FuncSpaceType::RangeType          RangeType;
      typedef typename FuncSpaceType::DomainType         DomainType;
      typedef typename FuncSpaceType::JacobianRangeType  JacobianRangeType;
      typedef FieldMatrix<double,dimDomain+1,dimDomain> FluxRangeType;
    };



  private:
    typedef StokesSingModel< ModelParamType,FunctionSpaceType,PressureSpaceType > ThisType;
    typedef LinearEllipticModelDefault< FunctionSpaceType, ThisType > BaseType;
    
  public:
    typedef typename BaseType :: BoundaryType BoundaryType;

 //    typedef typename FunctionSpaceType :: DomainType DomainType;
    // typedef typename FunctionSpaceType :: RangeType RangeType;
//     typedef typename FunctionSpaceType :: JacobianRangeType JacobianRangeType;

//     typedef typename FunctionSpaceType :: DomainFieldType DomainFieldType;
//     typedef typename FunctionSpaceType :: RangeFieldType RangeFieldType;

  public:
    //! return boundary type of a boundary point p used in a quadrature
    template< class IntersectionIteratorType >
    inline BoundaryType boundaryType( const IntersectionIteratorType &intersection ) const
    {
        return BaseType :: Dirichlet;
    }

    //! determine dirichlet value in a boundary point used in a quadrature
    template< class IntersectionIteratorType, class QuadratureType >
    inline void dirichletValues( const IntersectionIteratorType &intersection,
                                 const QuadratureType &quadrature,
                                 int p,
                                 RangeType &ret ) const
    {
      typedef typename IntersectionIteratorType :: Entity EntityType;
      
      const int dimension = DomainType :: dimension;
      
      const DomainType &x = intersection.inside()->geometry().global( quadrature.point( p ) );
   
      DomainType tmp=x;
      tmp[0]-=shift;
      tmp[1]-=shift;
      double abs=tmp*tmp  ;
      
      ret[0]=tmp[0];
      ret[0]/=abs;
      ret[1]=tmp[1];
      ret[1]/=abs;    

  

    }

    //! determine neumann value in a boundary point used in a quadrature
    template< class IntersectionIteratorType, class QuadratureType >
    inline void neumannValues( const IntersectionIteratorType &intersection,
                               const QuadratureType &quadrature,
                               int p,
                               RangeType &ret ) const
    {
      std :: cout << "Neumann boundary values are not implemented." << std :: endl;
      assert( false );

      //const DomainType &x = entity.geometry().global( quadrature.point( p ) ); 
            
      ret[ 0 ] = 0.0;
    }

    //! determine robin value in a boundary point used in a quadrature
    template< class IntersectionIteratorType, class QuadratureType >
    inline void robinValues( const IntersectionIteratorType &intersection,
                             const QuadratureType &quadrature,
                             int p, 
                             RangeType &ret ) const
    {
      std :: cout << "Robin boundary values are not implemented." << std :: endl;
      assert( false );

      //const DomainType &x = entity.geometry().global( quadrature.point( p ) ); 
            
      ret[ 0 ] = 0.0;
    }

    //! determine mass (i.e. value of the function c) in a quadrature point
    template< class EntityType, class QuadratureType >
    inline void mass( const EntityType &entity,
                      const QuadratureType &quadrature,
                      int p,
                      RangeType &ret ) const
    {

      assert(false);
      //const DomainType &x = entity.geometry().global( quadrature.point( p ) ); 

      ret[ 0 ] = 0.0;
    }

    //! Determine source (i.e. value of the function f) in a quadrature point
    template< class EntityType, class QuadratureType >
    inline void source( const EntityType &entity,
                        const QuadratureType &quadrature,
                        int p,
                        RangeType &ret ) const
    {
      const int dimension = DomainType :: dimension;
      
      const DomainType &x = entity.geometry().global( quadrature.point( p ) );
      
      ret[0]=sin(x[1]);
      ret[0]*=2;
      ret[0]*=exp(x[0]);
   
      
      ret[1]=cos(x[1]);
      ret[1]*=2;
      ret[1]*=exp(x[0]);
   
      //    return ret; 

    
    }

    //! no direct access to stiffness and velocity, but whole flux, i.e.
    //! diffflux = stiffness * grad( phi )
    template< class EntityType, class QuadratureType >
    inline void diffusiveFlux( const EntityType &entity,
                               const QuadratureType &quadrature,
                               int p,
                               const JacobianRangeType &gradphi, 
                               JacobianRangeType &ret ) const
    {
      assert(false);
      // - laplace phi = -div (1 * grad phi - 0 * phi)
      ret = gradphi;          
    }

    //! no direct access to stiffness and velocity, but whole flux, i.e.
    //! convectiveFlux =  - velocity * phi 
    template <class EntityType, class QuadratureType>  
    inline void convectiveFlux( const EntityType &entity,
                                const QuadratureType &quadrature,
                                int p,
                                const RangeType &phi, 
                                DomainType &ret ) const
    {
      assert(false);
      // - laplace phi = -div (1 * grad phi - 0 * phi)
      ret = 0.0;
    }
    
    template <class EntityType, class QuadratureType>  
    inline void convectiveFlux( const EntityType &entity,
                                const QuadratureType &quadrature,
                                int p,
                                const RangeType &phi, 
                                JacobianRangeType &ret ) const
    {
      // - laplace phi = -div (1 * grad phi - 0 * phi)
      ret = 0.0;
    }





    //! the coefficient for robin boundary condition
    template< class IntersectionIteratorType, class QuadratureType >
    inline double robinAlpha( const IntersectionIteratorType &intersection,
                              const QuadratureType &quadrature,
                              int p ) const
    { 
      return 1.0;
    }
  };  // end of PoissonModel class


  /*======================================================================*/
  /*!
   *  \class PoissonExactSolution
   *  \brief The class provides the exact solution for the model given by 
   *         the PoissonModel class
   *
   *  The function represents u = x(1-x)y(1-y), which is the solution of
   *  the model problem  - div grad u = 2(x(1-x)+ y(1-y)) on
   *  the unit square with homogeneous Dirichlet boundary values. 
   *  The function can be used for EOC calculation
   */
  /*======================================================================*/

  template< class FunctionSpaceImp >
  class StokesSingExactSolution
    : public Function< FunctionSpaceImp, StokesSingExactSolution< FunctionSpaceImp > >
  {
  public:
    typedef FunctionSpaceImp FunctionSpaceType;

  private:
    typedef StokesSingExactSolution< FunctionSpaceType > ThisType;
    typedef Function< FunctionSpaceType, ThisType > BaseType;

  public:
    typedef typename FunctionSpaceType :: DomainType DomainType;
    typedef typename FunctionSpaceType :: RangeType RangeType;
    typedef typename FunctionSpaceType :: JacobianRangeType JacobianRangeType;

    typedef typename FunctionSpaceType :: DomainFieldType DomainFieldType;
    typedef typename FunctionSpaceType :: RangeFieldType RangeFieldType;

  public:
    inline StokesSingExactSolution ( const FunctionSpaceType &functionSpace )
    : BaseType( functionSpace )
    {
    }

 
    //! u( x ) = sin( pi x_1 ) * ... * sin( pi x_n )
    inline void evaluate ( const DomainType &x, RangeType &ret ) const
    {
      enum { dimension = DomainType :: dimension };
      
      DomainType tmp=x;
      tmp[0]-=shift;
      tmp[1]-=shift;
      double abs=tmp*tmp  ;
      
      ret[0]=tmp[0];
      ret[0]/=abs;
      ret[1]=tmp[1];
      ret[1]/=abs;    

     

      
      // ret[0]+=1.0;
    }

    inline void evaluate( const DomainType &x, RangeFieldType t, RangeType &ret ) const
    {
      evaluate( x , ret );
    }
  };
/*======================================================================*/
  /*!
   *  \class PoissonExactGradient
   */
  /*======================================================================*/
 template< class FunctionSpaceImp >
  class StokesSingExactPressure
    : public Function< FunctionSpaceImp, StokesSingExactPressure< FunctionSpaceImp > >
  {
  public:
    typedef FunctionSpaceImp FunctionSpaceType;

  private:
    typedef StokesSingExactPressure< FunctionSpaceType > ThisType;
    typedef Function< FunctionSpaceType, ThisType > BaseType;

  public:
    typedef typename FunctionSpaceType :: DomainType DomainType;
    typedef typename FunctionSpaceType :: RangeType RangeType;
    typedef typename FunctionSpaceType :: JacobianRangeType JacobianRangeType;

    typedef typename FunctionSpaceType :: DomainFieldType DomainFieldType;
    typedef typename FunctionSpaceType :: RangeFieldType RangeFieldType;

  public:
    inline StokesSingExactPressure ( const FunctionSpaceType &functionSpace )
    : BaseType( functionSpace )
    {
    }

 
    //! u( x ) = sin( pi x_1 ) * ... * sin( pi x_n )
    inline void evaluate ( const DomainType &x, RangeType &ret ) const
    {
      enum { dimension = DomainType :: dimension };
       
 
      ret[0] = sin(x[1]);
      ret[0] *=exp(x[0]);
      ret[0] *=2.0;

    }

    inline void evaluate( const DomainType &x, RangeFieldType t, RangeType &ret ) const
    {
      evaluate( x , ret );
    }
  };





 






} // end namespace Dune

#endif
