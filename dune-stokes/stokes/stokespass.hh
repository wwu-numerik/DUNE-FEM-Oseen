/** \file stokespass.hh
    \brief  stokespass.hh
 **/

#ifndef STOKESPASS_HH
#define STOKESPASS_HH

#include <dune/fem/pass/pass.hh>
#include <dune/fem/operator/matrix/spmatrix.hh>
#include <dune/fem/space/dgspace.hh>

namespace Dune
{
template <class VelocityDiscreteFunctionImp,
        class PressureDiscreteFunctionImp, class DiscreteModelImp, class PreviousPassImp, int PassID = 0 >
class StokesPass : public LocalPass < DiscreteModelImp, PreviousPassImp, PassID >
{
    public:

        //! base type
        typedef LocalPass < DiscreteModelImp, PreviousPassImp, PassID >
            BaseType;

        //! previous pass type
        typedef PreviousPassImp
            PreviousPassType;

        //! discrete model type
        typedef DiscreteModelImp
            DiscreteModelType;

        //! discrete function type for the velocity
//        typedef  VelocityDiscreteFunctionType;
//            typedef PressureDiscreteFunctionImp PressureDiscreteFunctionType;
//            typedef typename PressureDiscreteFunctionType::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;
//            typedef typename PressureDiscreteFunctionType::DiscreteFunctionSpaceType PressureDiscreteFunctionSpaceType;
//            typedef typename VelocityDiscreteFunctionType::DiscreteFunctionSpaceType VelocityDiscreteFunctionSpaceType;
//
//            typedef SparseRowMatrix<double> MatrixType;
//            typedef MatrixType B_OperatorType;
//            typedef MatrixType B_Transposed_OperatorType;
//            typedef MatrixType C_OperatorType;
//
//            //typedefs for interface compliance, not definitive
//            typedef typename BaseType::ArgumentType ArgumentType;
//            typedef typename BaseType::DestinationType DestinationType;
//            typedef typename BaseType::Entity EntityType;
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
//            virtual void prepare( const ArgumentType& arg,
//                                    DestinationType& dest) const
//            {}
//
//            virtual void finalize(const ArgumentType& arg,
//                                    DestinationType& dest) const
//            {}
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
//            VelocityDiscreteFunctionSpaceType& velo_space_;
//            PressureDiscreteFunctionSpaceType& press_space_;
//            DiscreteGradientSpaceType& grad_space_;
//
//            VelocityGradMatType gradMatrix_;
//            VelocityDivMatType divMatrix_;
//            VelocityStabMatType stabMatrix_;
//            PressureGradMatType pressureGradMatrix_;
//            PressureDivMatType pressureDivMatrix_;
//            PressureStabMatType pressureStabMatrix_;
//
//            const DiscreteModelType& disc_model_;

};

}
#endif  // end of stokespass.hh
