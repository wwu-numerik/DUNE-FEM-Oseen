/** \file stokespass.hh
    \brief  stokespass.hh
 **/

#ifndef STOKESPASS_HH
#define STOKESPASS_HH

#include <dune/fem/pass/pass.hh>
#include<dune/fem/operator/matrix/spmatrix.hh>
#include <dune/fem/space/dgspace.hh>

namespace Dune
{
    template <class VelocityDiscreteFunctionImp,
            class PressureDiscreteFunctionImp, class DiscreteModelImp, class PreviousPassImp, int PassID = 0 >
    class StokesPass : public LocalPass < DiscreteModelImp, PreviousPassImp, PassID >
    {

        public:

            //! template repetions etc
            typedef LocalPass < DiscreteModelImp, PreviousPassImp, PassID > BaseType;
            typedef PreviousPassImp PreviousPassType;
            typedef VelocityDiscreteFunctionImp VelocityDiscreteFunctionType;
            typedef PressureDiscreteFunctionImp PressureDiscreteFunctionType;
            typedef typename PressureDiscreteFunctionType::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;
            typedef typename PressureDiscreteFunctionType::DiscreteFunctionSpaceType PressureDiscreteFunctionSpaceType;
            typedef typename VelocityDiscreteFunctionType::DiscreteFunctionSpaceType VelocityDiscreteFunctionSpaceType;

            typedef SparseRowMatrix<double> MatrixType;
            typedef MatrixType B_OperatorType;
            typedef MatrixType B_Transposed_OperatorType;
            typedef MatrixType C_OperatorType;

            //typedefs for interface compliance, not definitive
            typedef typename BaseType::ArgumentType ArgumentType;
            typedef typename BaseType::DestinationType DestinationType;
            typedef typename BaseType::Entity EntityType;

        private:
            enum { VelocityDimDomain = VelocityDiscreteFunctionSpaceType::DimDomain };
            enum { VelocityDimRange = VelocityDiscreteFunctionSpaceType::DimRange };
            enum { VelocityPolOrder = VelocityDiscreteFunctionSpaceType::polOrd };

            //no need for diff velo/press typedefs, these will always be the same
            typedef typename VelocityDiscreteFunctionSpaceType::GridType GridType;
            typedef typename VelocityDiscreteFunctionSpaceType::GridPartType GridPartType;

            typedef typename VelocityDiscreteFunctionSpaceType::DomainFieldType VelocityDomainFieldType;
            typedef typename VelocityDiscreteFunctionSpaceType::RangeFieldType VelocityRangeFieldType;

            //afaics these are ony ever used internally
            typedef MatrixFunctionSpace<VelocityDomainFieldType, VelocityRangeFieldType, VelocityDimDomain, VelocityDimRange, VelocityDimRange > GradientSpaceType;
            typedef DiscontinuousGalerkinSpace<GradientSpaceType, GridPartType, VelocityPolOrder> DiscreteGradientSpaceType;



        public:

            StokesPass(PreviousPassType& prev_pass, const DiscreteFunctionSpaceType& disc_space)
                : BaseType( prev_pass, disc_space )
                //rhs( "pass_rhs"
            {}

            const B_OperatorType& Get_B_Operator() const { return b_op_; }
            const B_Transposed_OperatorType& Get_B_Transposed_Operator() const { return b_transp_op_; }
            const C_OperatorType& Get_C_Operator() const { return c_op_; }

            MatrixType& systemMatrix()
            {}

            const VelocityDiscreteFunctionType& rhs1() const
            {}

            virtual void prepare( const ArgumentType& arg,
                                    DestinationType& dest) const
            {}

            virtual void finalize(const ArgumentType& arg,
                                    DestinationType& dest) const
            {}

            virtual void applyLocal(EntityType& entity) const
            {
                    //- typedefs
                typedef typename VelocityDiscreteFunctionSpaceType::IndexSetType IndexSetType;
                typedef typename VelocityDiscreteFunctionSpaceType::BaseFunctionSetType BaseFunctionSetType;

                typedef typename DiscreteGradientSpaceType::IndexSetType GradientIndexSetType;
                typedef typename DiscreteGradientSpaceType::BaseFunctionSetType GradientBaseFunctionSetType;

                typedef typename PressureDiscreteFunctionSpaceType::IndexSetType PressureIndexSetType;
                typedef typename PressureDiscreteFunctionSpaceType::BaseFunctionSetType PressureBaseFunctionSetType;
                //- statements

                typedef typename DestinationType::LocalFunctionType LocalFuncType;
                VolumeQuadratureType volQuad(entity, volumeQuadOrd_);


                double massVolElInv;
                double vol = volumeElement(entity, volQuad,massVolElInv);

                const BaseFunctionSetType& bsetEn = spc_.baseFunctionSet(entity);
                const int numDofs = bsetEn.numBaseFunctions();

                const GradientBaseFunctionSetType& gradbsetEn = gradientspc_.baseFunctionSet(entity);
                const int gradientNumDofs = gradbsetEn.numBaseFunctions();

                const PressureBaseFunctionSetType& pressurebsetEn = pressurespc_.baseFunctionSet(entity);
                const int pressureNumDofs = pressurebsetEn.numBaseFunctions();
                LocalGradMatType en_grad = gradMatrix_.localMatrix(entity,entity);
                LocalDivMatType en_div= divMatrix_.localMatrix(entity,entity);
                LocalStabMatType en_stab=stabMatrix_.localMatrix(entity,entity);
                LocalPressureGradMatType en_pressGrad=pressureGradMatrix_.localMatrix(entity,entity);
                LocalPressureDivMatType en_pressDiv=pressureDivMatrix_.localMatrix(entity,entity);
                LocalPressureStabMatType en_pressStab=pressureStabMatrix_.localMatrix(entity,entity);
                double nu=problem_.get_nu();


                //compute volume integral contributions
                for (int l = 0; l < quadNop ; ++l)
                {
                } // end volume integral contributions


                //compute surface integral contributions
                IntersectionIterator endnit = entity.ileafend();
                double dtLocal = 0.0;
                double minvol = vol;
                for (IntersectionIterator nit = entity.ileafbegin(); nit != endnit; ++nit)
                {
                    //	  int twistSelf = twistUtil_.twistInSelf(nit);
                    FaceQuadratureType faceQuadInner(spc_.gridPart(),nit, faceQuadOrd_,
                    FaceQuadratureType::INSIDE);

                    if (nit.neighbor())
                    {
                    } //end if inner

                    if (nit.boundary())
                    {
                    } // end if boundary
                }
                //end compute surface integral contributions
            }//end void applyLocal(EntityType& entity) const


        private:
            //will prolly be generated on-the-fly, not stored as members
            B_OperatorType b_op_;
            B_Transposed_OperatorType b_transp_op_;
            C_OperatorType c_op_;

            //more dummies
           // VelocityDiscreteFunctionType& rhs;

    };
}
#endif  // end of stokespass.hh
