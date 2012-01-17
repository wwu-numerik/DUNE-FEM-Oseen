#ifndef BCRSTRAITS_HH
#define BCRSTRAITS_HH

#include <dune/fem/operator/2order/dgmatrixsetup.hh>
#include <dune/oseen/assembler/mod_istlmatrix.hh>
#include <dune/stuff/matrix_patterns.hh>

namespace Dune {
namespace Stokes {
namespace Integrators {



template <class RowFunctionImp, class ColFunctionImp >
struct ModifiedDGMatrixTraits
{
    typedef typename ColFunctionImp::DiscreteFunctionSpaceType
        ColSpaceImp;
    typedef typename RowFunctionImp::DiscreteFunctionSpaceType
        RowSpaceImp;
  typedef RowSpaceImp RowSpaceType;
  typedef ColSpaceImp ColumnSpaceType;
  typedef Stuff::Matrix::ElementNeighborStencil<RowSpaceType,ColumnSpaceType> StencilType;
  typedef Dune::ParallelScalarProduct< ColumnSpaceType > ParallelScalarProductType;

  template< class M >
  struct Adapter
  {
    typedef Dune::DGParallelMatrixAdapter< M > MatrixAdapterType;
  };

};


} //namespace Dune
} //namespace Stokes
} //namespace Integrators


#endif // BCRSTRAITS_HH
