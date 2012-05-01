/**************************************************************************
**       Title: Example of a ElementIntegratorTraits class implementation
**    $RCSfile$
**   $Revision: 2143 $$Name$
**       $Date: 2007-08-16 10:00:52 +0200 (Thu, 16 Aug 2007) $
**   Copyright: GPL $Author: nolte $
** Description: the file contains a default implementation of a
**              ElementIntegratorTraits class to be used with an
**              appropriate Model in an FEOp or RhsAssembler for
**              solving a general elliptic problem.
**
**************************************************************************/

#ifndef DUNE_ELEMENTINTEGRATORTRAITS_HH
#define DUNE_ELEMENTINTEGRATORTRAITS_HH

#include <dune/fem/space/lagrangespace.hh>
#include <dune/fem/function/adaptivefunction.hh>
#include <dune/fem/quadrature/cachequad.hh>
#include <dune/fem/operator/matrixadapter.hh>

namespace Dune
{

    /*======================================================================*/
    /*!
     *  \class DefaultElementIntegratorTraits
     *  \brief The DefaultElementIntegratorTraits provides Type-Information
     *         for the ElementMatrices and FEOp operator.
     *
     *  default implementation of a ElementIntegratorTraits class to be used
     *  with an appropriate Model in an FEOp for solving a general elliptic
     *  problem.
     *
     *  It is only
     *  considered to yield information by its types, no member variables or
     *  methods are provided, neither is it instantiated at any time.
     *
     *  Currently scalar functions and Lagrange-Basis of degree 1 are used,
     *  elementquadratures are chosen for any quadrature in the FEOp and
     *  ElementIntegrators.
     *
     *  All types are derived from GridType and dimworld obtained by
     *  inclusion of gridtype.hh
     *
     *  Essential Datatypes without explicit interface:
     *
     *  required for ElementQuadratureTypes:
     *     constructor with arguments (entity, order)
     *     nop, weight, point methods
     *
     *  required for IntersectionQuadratureTypes:
     *     enum INSIDE
     *     constructor with arguments (entity, order, INSIDE)
     *     nop, weight, point methods
     *     localpoint, geometry methods
     *
     *  required for LocalMatrixType (e.g. satisfied by
     *  FieldMatrixAdapter<Fieldmatrix<...>>):
     *     constructor without arguments
     *     rows(), cols() methods
     *     add(rown, coln, value) allows writable additive access
     *     to ij-th component.
     *
     *  The class can be taken as an example for own implementations.
     */
    /*======================================================================*/
    template< class GridImp, int polOrder >
    struct EllipticElementIntegratorTraits
    {
        typedef GridImp GridType;

        typedef LeafGridPart< GridType > GridPartType;

        enum { dimworld = GridType :: dimensionworld };
        enum { dim = GridType :: dimension };

        typedef FunctionSpace< double, double, dimworld, 1 > FunctionSpaceType;
        typedef LagrangeDiscreteFunctionSpace< FunctionSpaceType, GridPartType, polOrder >
        DiscreteFunctionSpaceType;
        typedef AdaptiveDiscreteFunction< DiscreteFunctionSpaceType >
        DiscreteFunctionType;

        enum { elementMatrixSize = 100 };
        typedef FieldMatrixAdapter
        < FieldMatrix< double, elementMatrixSize,elementMatrixSize > >
        ElementMatrixType;

        typedef CachingQuadrature< GridPartType, 0 > ElementQuadratureType;
        typedef CachingQuadrature< GridPartType, 1 > IntersectionQuadratureType;
        enum { quadDegree = 2*polOrder+1 }; //<! degree of quadrature

        //! derived types
        typedef typename GridPartType :: IntersectionIteratorType
        IntersectionIteratorType;
        typedef typename GridType :: template Codim< 0 > :: Entity EntityType;
        typedef typename GridType :: template Codim< 0 > :: EntityPointer
        EntityPointerType;
        typedef typename EntityType :: ctype CoordType;

        typedef typename DiscreteFunctionSpaceType :: DomainFieldType
        DomainFieldType;
        typedef typename DiscreteFunctionSpaceType :: RangeFieldType
        RangeFieldType;
        typedef typename DiscreteFunctionSpaceType :: DomainType DomainType;
        typedef typename DiscreteFunctionSpaceType :: RangeType RangeType;
        typedef typename DiscreteFunctionSpaceType :: JacobianRangeType
        JacobianRangeType;

        typedef typename DiscreteFunctionSpaceType :: BaseFunctionSetType
        BaseFunctionSetType;
    };

}

#endif

/** Copyright (c) 2012, Felix Albrecht, Rene Milk      
 * All rights reserved.
 * 
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met: 
 * 
 * 1. Redistributions of source code must retain the above copyright notice, this
 *    list of conditions and the following disclaimer. 
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution. 
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
 * ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 * 
 * The views and conclusions contained in the software and documentation are those
 * of the authors and should not be interpreted as representing official policies, 
 * either expressed or implied, of the FreeBSD Project.
**/

