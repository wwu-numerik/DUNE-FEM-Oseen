#ifndef DUNE_MATRIX_OPERATOR_HH
#define DUNE_MATRIX_OPERATOR_HH

/** \file
    \brief matrixoperator.hh
 */

#include <dune/fem/function/common/discretefunction.hh>
//#include "../discretefunction/common/discretefunction.hh"
#include <dune/fem/operator/common/operator.hh>
//#include "common/operator.hh"

//#include <dune/fem/operator/feop/spmatrix.hh>

namespace Dune {
  
  //Matrix-Vector-Multiplikation for DiscreteFunctions
  
  template<class MatrixType,class ArgumentType,class DestinationType>
  class MatrixOperator:public Operator<
    typename ArgumentType::DomainFieldType,
    typename ArgumentType::RangeFieldType,
    ArgumentType,DestinationType>
  {
    
    
  public:
    MatrixOperator():myMatrix_(){};

    MatrixOperator (MatrixType& A):myMatrix_(A){/*Dimension Checken*/}
    
    virtual void operator()(const ArgumentType& arg,DestinationType& dest) const
    {
      const double* argpt=arg.leakPointer();
      dest.clear();
      double* destpt=dest.leakPointer();
      
      myMatrix_.multOEM(argpt,destpt);
    }
    

    


  protected:
    mutable MatrixType&  myMatrix_;

  };

}
#endif
