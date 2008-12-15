/** \file stokespass.hh
    \brief  stokespass.hh
 **/

#ifndef STOKESPASS_HH
#define STOKESPASS_HH

#include <dune/fem/pass/pass.hh>
#include <dune/fem/operator/matrix/spmatrix.hh>
#include <dune/fem/space/dgspace.hh>

#include "../src/stuff.hh" // should be removed in the end

namespace Dune
{
template <  class DiscreteModelImp,
            class PreviousPassImp,
            int PassID = 0 >
class StokesPass
    : public Pass < DiscreteModelImp, PreviousPassImp, PassID >
{
    public:

        //! base type
        typedef Pass < DiscreteModelImp, PreviousPassImp, PassID >
            BaseType;

        //! previous pass type
        typedef PreviousPassImp
            PreviousPassType;

        //! discrete model type
        typedef DiscreteModelImp
            DiscreteModelType;

        //! volume quadrature type
        typedef typename DiscreteModelType::VolumeQuadratureType
            VolumeQuadratureType;

        //! face quadrature type
        typedef typename DiscreteModelType::FaceQuadratureType
            FaceQuadratureType;

        //! type of discrete function space wrapper
        typedef typename DiscreteModelType::DiscreteStokesFunctionSpaceWrapperType
            DiscreteStokesFunctionSpaceWrapperType;

        //! discrete function wrapper type
        typedef typename DiscreteModelType::DiscreteStokesFunctionWrapperType
            DiscreteStokesFunctionWrapperType;

        //! discrete function type for the velocity
        typedef typename DiscreteStokesFunctionWrapperType::DiscreteVelocityFunctionType
            DiscreteVelocityFunctionType;

        //! discrete function space type for the velocity
        typedef typename DiscreteVelocityFunctionType::DiscreteFunctionSpaceType
            DiscreteVelocityFunctionSpaceType;

        //! discrete function type for sigma
        typedef typename DiscreteModelType::DiscreteSigmaFunctionType
            DiscreteSigmaFunctionType;

        //! discrete function space type for sigma
        typedef typename DiscreteSigmaFunctionType::DiscreteFunctionSpaceType
            DiscreteSigmaFunctionSpaceType;

        //! discrete fucntion type for the pressure
        typedef typename DiscreteStokesFunctionWrapperType::DiscretePressureFunctionType
            DiscretePressureFunctionType;

        //! discrete function space type for the pressure
        typedef typename DiscretePressureFunctionType::DiscreteFunctionSpaceType
            DiscretePressureFunctionSpaceType;

        //! Coordinate type (world coordinates)
        typedef typename DiscreteVelocityFunctionSpaceType::DomainType
            WorldCoordinateType;

        //! Vector type of the velocity's discrete function space's range
        typedef typename DiscreteVelocityFunctionSpaceType::RangeType
            VelocityRangeType;

        //! vector type of sigmas' discrete functions space's range
        typedef typename DiscreteSigmaFunctionSpaceType::RangeType
            SigmaRangeType;

        //! Vector type of the pressure's discrete function space's range
        typedef typename DiscretePressureFunctionSpaceType::RangeType
            PressureRangeType;

        //! Type of GridPart
        typedef typename DiscreteVelocityFunctionSpaceType::GridPartType
            GridPartType;

        //! Intersection iterator of the gridpart
        typedef typename GridPartType::IntersectionIteratorType
            IntersectionIteratorType;

        //! entity iterator of the gridpart
        typedef typename GridPartType::template Codim< 0 >::IteratorType
            EntityIteratorType;

        /**
         *  \name typedefs for interface compliance
         *  \{
         **/
        typedef typename BaseType::DestinationType
            DestinationType;

        typedef typename BaseType::DomainType
            DomainType;

        typedef typename BaseType::RangeType
            RangeType;

        typedef typename BaseType::TotalArgumentType
            TotalArgumentType;
        /**
         *  \}
         **/


        /**
         *  \brief  constructor
         *  \todo   doc
         **/
        StokesPass( PreviousPassType& prevPass,
                    DiscreteModelType& discreteModel,
                    GridPartType& gridPart,
                    DiscreteStokesFunctionSpaceWrapperType& spaceWrapper )
                : BaseType( prevPass ),
                discreteModel_( discreteModel ),
                gridPart_( gridPart ),
                spaceWrapper_( spaceWrapper ),
                velocitySpace_( spaceWrapper.discreteVelocitySpace() ),
                pressureSpace_( spaceWrapper.discretePressureSpace() ),
                sigmaSpace_( gridPart )
        {}

        /**
         *  \brief  empty constructor
         **/
        StokesPass()
        {}

        virtual void apply( const DomainType &arg, RangeType &dest) const
        {
            // walk the grid
            EntityIteratorType entityItEnd = velocitySpace_.end();
            for ( EntityIteratorType entityIt = velocitySpace_.begin(); entityIt != entityIt; ++entityIt ) {
            } // done walking the grid
        }

        virtual void compute( const TotalArgumentType &arg, DestinationType &dest) const
        {}

        virtual void allocateLocalMemory()
        {}

    private:
        DiscreteModelType& discreteModel_;
        GridPartType& gridPart_;
        DiscreteStokesFunctionSpaceWrapperType& spaceWrapper_;
        DiscreteVelocityFunctionSpaceType& velocitySpace_;
        DiscretePressureFunctionSpaceType& pressureSpace_;
        DiscreteSigmaFunctionSpaceType sigmaSpace_;

};

}
#endif  // end of stokespass.hh




//            typedef SparseRowMatrix<double> MatrixType;
//            typedef MatrixType B_OperatorType;
//            typedef MatrixType B_Transposed_OperatorType;
//            typedef MatrixType C_OperatorType;
//

//        private:
//            enum { VelocityDimDomain = VelocityDiscreteFunctionSpaceType::DimDomain };
//            enum { VelocityDimRange = VelocityDiscreteFunctionSpaceType::DimRange };
//            enum { VelocityPolOrder = VelocityDiscreteFunctionSpaceType::polOrd };
//
//            //no need for diff velo/press typedefs, these will always be the same
//            typedef typename VelocityDiscreteFunctionSpaceType::GridType GridType;
//            typedef typename VelocityDiscreteFunctionSpaceType::GridPartType GridPartType;
//
//            typedef typename VelocityDiscreteFunctionSpaceType::DomainFieldType VelocityDomainFieldType;
//            typedef typename VelocityDiscreteFunctionSpaceType::RangeFieldType VelocityRangeFieldType;
//
//            typedef typename DiscreteModelType::VolumeQuadratureType VolumeQuadratureType;
//            typedef typename DiscreteModelType::FaceQuadratureType FaceQuadratureType;
//
//            //afaics these are ony ever used internally
//            typedef MatrixFunctionSpace<VelocityDomainFieldType, VelocityRangeFieldType, VelocityDimDomain, VelocityDimRange, VelocityDimRange > GradientSpaceType;
//            typedef DiscontinuousGalerkinSpace<GradientSpaceType, GridPartType, VelocityPolOrder> DiscreteGradientSpaceType;
//
//            typedef SparseRowMatrixObject<DiscreteGradientSpaceType,VelocityDiscreteFunctionSpaceType> VelocityGradMatType;
//            typedef typename VelocityGradMatType::LocalMatrixType LocalVelocityGradMatType;
//
//            typedef SparseRowMatrixObject<VelocityDiscreteFunctionSpaceType,DiscreteGradientSpaceType> VelocityDivMatType;
//            typedef typename VelocityDivMatType::LocalMatrixType LocalVelocityDivMatType;
//
//            typedef SparseRowMatrixObject<VelocityDiscreteFunctionSpaceType,VelocityDiscreteFunctionSpaceType> VelocityStabMatType;
//            typedef typename VelocityStabMatType::LocalMatrixType LocalVelocityStabMatType;
//
//            typedef SparseRowMatrixObject<VelocityDiscreteFunctionSpaceType,PressureDiscreteFunctionSpaceType> PressureGradMatType;
//            typedef typename PressureGradMatType::LocalMatrixType LocalPressureGradMatType;
//
//            typedef SparseRowMatrixObject<PressureDiscreteFunctionSpaceType,VelocityDiscreteFunctionSpaceType> PressureDivMatType;
//            typedef typename PressureDivMatType::LocalMatrixType LocalPressureDivMatType;
//
//            typedef SparseRowMatrixObject<PressureDiscreteFunctionSpaceType,PressureDiscreteFunctionSpaceType> PressureStabMatType;
//            typedef typename PressureStabMatType::LocalMatrixType LocalPressureStabMatType;
//
//
//        public:
//
//            StokesPass( PreviousPassType& prev_pass,
//                        const DiscreteFunctionSpaceType& velo_space,
//                        const PressureDiscreteFunctionSpaceType& press_space,
//                        const DiscreteModelType& disc_model)
//                : BaseType( prev_pass, velo_space ),
//                velo_space_( velo_space ),
//                press_space_( press_space ),
//                grad_space_( press_space.gridPart() ),
//                disc_model_( disc_model )
//            {}
//
//            const B_OperatorType& Get_B_Operator() const { return b_op_; }
//            const B_Transposed_OperatorType& Get_B_Transposed_Operator() const { return b_transp_op_; }
//            const C_OperatorType& Get_C_Operator() const { return c_op_; }
//
//            MatrixType& systemMatrix()
//            {}
//
//            const VelocityDiscreteFunctionType& rhs1() const
//            {}
//
//
//            virtual void applyLocal(EntityType& entity) const
//            {
//                    //- typedefs
//                typedef typename VelocityDiscreteFunctionSpaceType::IndexSetType IndexSetType;
//                typedef typename VelocityDiscreteFunctionSpaceType::BaseFunctionSetType BaseFunctionSetType;
//
//                typedef typename DiscreteGradientSpaceType::IndexSetType GradientIndexSetType;
//                typedef typename DiscreteGradientSpaceType::BaseFunctionSetType GradientBaseFunctionSetType;
//
//                typedef typename PressureDiscreteFunctionSpaceType::IndexSetType PressureIndexSetType;
//                typedef typename PressureDiscreteFunctionSpaceType::BaseFunctionSetType PressureBaseFunctionSetType;
//                //- statements
//
//                typedef typename DestinationType::LocalFunctionType LocalFuncType;
//                VolumeQuadratureType volQuad(entity, 2); //! \todo get order from model too
//
//
//                double massVolElInv;
//                double vol = volumeElement(entity, volQuad,massVolElInv);
//
//                const BaseFunctionSetType& bsetEn = velo_space_.baseFunctionSet(entity);
//                const int numDofs = bsetEn.numBaseFunctions();
//
//                const GradientBaseFunctionSetType& gradbsetEn = grad_space_.baseFunctionSet(entity);
//                const int gradientNumDofs = gradbsetEn.numBaseFunctions();
//
//                const PressureBaseFunctionSetType& pressurebsetEn = press_space_.baseFunctionSet(entity);
//                const int pressureNumDofs = pressurebsetEn.numBaseFunctions();
//                LocalVelocityGradMatType en_grad = gradMatrix_.localMatrix(entity,entity);
//                LocalVelocityDivMatType en_div= divMatrix_.localMatrix(entity,entity);
//                LocalVelocityStabMatType en_stab=stabMatrix_.localMatrix(entity,entity);
//                LocalPressureGradMatType en_pressGrad=pressureGradMatrix_.localMatrix(entity,entity);
//                LocalPressureDivMatType en_pressDiv=pressureDivMatrix_.localMatrix(entity,entity);
//                LocalPressureStabMatType en_pressStab=pressureStabMatrix_.localMatrix(entity,entity);


                //compute volume integral contributions
//                for (int l = 0; l < quadNop ; ++l)
//                {
//                } // end volume integral contributions


                //compute surface integral contributions
//                IntersectionIterator endnit = entity.ileafend();
//                double dtLocal = 0.0;
//                double minvol = vol;
//                for (IntersectionIterator nit = entity.ileafbegin(); nit != endnit; ++nit)
//                {
//                    //	  int twistSelf = twistUtil_.twistInSelf(nit);
//                    FaceQuadratureType faceQuadInner(spc_.gridPart(),nit, faceQuadOrd_,
//                    FaceQuadratureType::INSIDE);
//
//                    if (nit.neighbor())
//                    {
//                    } //end if inner
//
//                    if (nit.boundary())
//                    {
//                    } // end if boundary
//                }
                //end compute surface integral contributions
            //}//end void applyLocal(EntityType& entity) const


        //private:
            //will prolly be generated on-the-fly, not stored as members
//            B_OperatorType b_op_;
//            B_Transposed_OperatorType b_transp_op_;
//            C_OperatorType c_op_;
//
//            //more dummies
//            // VelocityDiscreteFunctionType& rhs;
//
//
//            VelocityGradMatType gradMatrix_;
//            VelocityDivMatType divMatrix_;
//            VelocityStabMatType stabMatrix_;
//            PressureGradMatType pressureGradMatrix_;
//            PressureDivMatType pressureDivMatrix_;
//            PressureStabMatType pressureStabMatrix_;
//
