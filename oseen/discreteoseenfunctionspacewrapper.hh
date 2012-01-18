/**
 *  \file   discretestokesfunctionspacewrapper.hh
 *  \todo   doc
 **/

#ifndef DISCRETESTOKESFUNCTIONSPACEWRAPPER_HH_INCLUDED
#define DISCRETESTOKESFUNCTIONSPACEWRAPPER_HH_INCLUDED

#include <dune/fem/space/common/discretefunctionspace.hh>
#include <dune/fem/space/common/restrictprolonginterface.hh>
#include <dune/fem/io/file/vtkio.hh>
#include <dune/fem/operator/projection/l2projection.hh>
#include <dune/stuff/adaption.hh>

namespace Dune
{

// forward
template < class DiscreteOseenFunctionSpaceWrapperTraitsImp >
class DiscreteOseenFunctionSpaceWrapper;

/**
 *  \todo   doc
 **/
template < class DiscreteVelocitySpaceImp, class DiscretePressureSpaceImp >
class DiscreteOseenFunctionSpaceWrapperTraits
{
    public:

        //! type of discrete velocity function space
        typedef DiscreteVelocitySpaceImp
            DiscreteVelocityFunctionSpaceType;

        //! type of discrete pressure function space
        typedef DiscretePressureSpaceImp
            DiscretePressureFunctionSpaceType;

        /**
         *  \name for interface compliance
         *  \{
         **/
        //! own type (CRTP) (for interface compliance)
        typedef DiscreteOseenFunctionSpaceWrapper<
                    DiscreteOseenFunctionSpaceWrapperTraits<   DiscreteVelocityFunctionSpaceType,
                                                                DiscretePressureFunctionSpaceType > >
            DiscreteFunctionSpaceType;

        //! type of function space
        typedef typename DiscreteVelocityFunctionSpaceType::FunctionSpaceType
            FunctionSpaceType;

        //! type of base function set
        typedef typename DiscreteVelocityFunctionSpaceType::BaseFunctionSetType
            BaseFunctionSetType;

        //! type of DoF mapper
        typedef typename DiscreteVelocityFunctionSpaceType::MapperType
            MapperType;

        //! type of block mapper
        typedef typename DiscreteVelocityFunctionSpaceType::BlockMapperType
            BlockMapperType;

        //! type of underlying grid part
        typedef typename DiscreteVelocityFunctionSpaceType::GridPartType
            GridPartType;

		template< class DiscreteFunction,
				  class Operation =  // get default type from traits
				  typename DiscreteVelocityFunctionSpaceType::Traits:: template CommDataHandle< DiscreteFunction > :: OperationType
				>
		//! TODO this probably isn't too sane to do...
		struct CommDataHandle
		{
		  //! type of communication data handle
		  typename DiscreteVelocityFunctionSpaceType::Traits
			:: template CommDataHandle< DiscreteFunction, Operation > :: Type
			Type;

		  //! type of operation to perform on scatter
		  typename DiscreteVelocityFunctionSpaceType::Traits
			:: template CommDataHandle< DiscreteFunction, Operation > :: OperationType
			OperationType;
		};

		//! type of communication manager
		typedef CommunicationManager< DiscreteFunctionSpaceType >
				CommunicationManagerType;
        /**
         *  \}
         **/
}; // end of DiscreteOseenFunctionSpaceWrapperTraits

/**
 *  \todo   doc
 **/
template < class DiscreteOseenFunctionSpaceWrapperTraitsImp >
class DiscreteOseenFunctionSpaceWrapper
    : public DiscreteFunctionSpaceDefault< DiscreteOseenFunctionSpaceWrapperTraitsImp >
{
    public:

        typedef DiscreteOseenFunctionSpaceWrapperTraitsImp
            Traits;

        //! base type
        typedef DiscreteFunctionSpaceDefault< DiscreteOseenFunctionSpaceWrapperTraitsImp >
            BaseType;

        //! type of discrete velocity function space
        typedef typename Traits::DiscreteVelocityFunctionSpaceType
            DiscreteVelocityFunctionSpaceType;

        //! type of discrete pressure function space
        typedef typename Traits::DiscretePressureFunctionSpaceType
            DiscretePressureFunctionSpaceType;

        //! own type (CRTP)
        typedef typename BaseType::DiscreteFunctionSpaceType
            DiscreteFunctionSpaceType;

        //! type of function space (of the discrete velocity function space)
        typedef typename BaseType::FunctionSpaceType
            FunctionSpaceType;

        //! type of base function set (of the discrete velocity function space)
        typedef typename BaseType::BaseFunctionSetType
            BaseFunctionSetType;

        //! type of DoF mapper (of the discrete velocity function space)
        typedef typename BaseType::MapperType
            MapperType;

        //! type of block mapper (of the discrete velocity function space)
        typedef typename BaseType::BlockMapperType
            BlockMapperType;

        //! type of underlying grid part (of the discrete velocity function space)
        typedef typename BaseType::GridPartType
            GridPartType;

        //! type of underlying grid (of the discrete velocity function space)
        typedef typename GridPartType::GridType
            GridType;

        //! type of used index set (of the discrete velocity function space)
        typedef typename GridPartType::IndexSetType
            IndexSetType;

        //! type of iterator of codim 0 entities for grid traversal (of the discrete velocity function space)
        typedef typename GridPartType::template Codim< 0 >::IteratorType
            IteratorType;

        //! type of codim 0 entity (of the discrete velocity function space)
        typedef typename IteratorType::Entity
            EntityType;

        /**
         *  \brief  constructor
         *  \todo   doc
         **/
        DiscreteOseenFunctionSpaceWrapper( GridPartType& gridPart )
            : BaseType( gridPart ),
            velocitySpace_( gridPart ),
            pressureSpace_( gridPart )
        {}

        /**
         *  \brief  destructor
         **/
		virtual ~DiscreteOseenFunctionSpaceWrapper()
        {}

        /**
         *  \brief  returns the discrete velocity space
         *  \todo   doc
         **/
        DiscreteVelocityFunctionSpaceType& discreteVelocitySpace()
        {
            return velocitySpace_;
        }
        const DiscreteVelocityFunctionSpaceType& discreteVelocitySpace() const
        {
            return velocitySpace_;
        }

        /**
         *  \brief  return the discrete pressure space
         *  \todo   doc
         **/
        DiscretePressureFunctionSpaceType& discretePressureSpace()
        {
            return pressureSpace_;
        }
        const DiscretePressureFunctionSpaceType& discretePressureSpace() const
        {
            return pressureSpace_;
        }

        /**
         *  \brief  return type identifier of the discrete velocity function space
         **/
        DFSpaceIdentifier type() const
        {
            return velocitySpace_.type();
        }

        /**
         *  \brief  get base function set of the discrete velocity function space for given entity
         *  \todo   doc
         **/
        inline const BaseFunctionSetType baseFunctionSet( const EntityType& entity ) const
        {
            return velocitySpace_.baseFunctionSet( entity );
        }

        /**
         *  \brief  returns true if the discrete velocity function space contains DoFs of given codimension
         *  \todo   doc
         **/
        inline bool contains( const int codim ) const
        {
            return velocitySpace_.contains( codim );
        }

        /**
         *  \brief  returns true if the discrete velocity function space contains only globally continuous functions
         *  \todo   doc
         **/
        inline bool continuous() const
        {
            return velocitySpace_.continuous();
        }

        /**
         *  \brief  get index of the sequence in the discrete velocity function spaces grid sequences
         *  \todo   doc
         **/
        inline int sequence() const
        {
            return velocitySpace_.sequence();
        }

        /**
         *  \brief  get global order of the discrete velocity function space
         *  \todo doc
         **/
        inline int order() const
        {
            return velocitySpace_.order();
        }

        /**
         *  \brief  get a reference to the discrete velocity function spacees DoF mapper
         *  \todo   doc
         **/
        inline MapperType& mapper() const
        {
            return velocitySpace_.mapper();
        }

        /**
         *  \brief  get a reference to the discrete velocity function spaces block mapper
         *  \todo   doc
         **/
        inline BlockMapperType& blockMapper() const
        {
            return velocitySpace_.blockMapper();
        }

        /**
         *  \brief  get reference to grid the discrete velocity function space belongs to
         *  \todo   doc
         **/
        inline const GridType& grid() const
        {
            return velocitySpace_.grid();
        }

        /**
         *  \brief  get reference to grid the discrete velocity function space belongs to
         *  \todo   doc
         **/
        inline GridType& grid()
        {
            return velocitySpace_.grid();
        }

        /**
         *  \brief  get a reference to the discrete velocity function space grid partition
         *  \todo   doc
         **/
        inline const GridPartType& gridPart() const
        {
            return velocitySpace_.gridPart();
        }

        /**
         *  \brief  get a reference to the discrete velocity function space grid partition
         *  \todo   doc
         **/
        inline GridPartType& gridPart()
        {
            return velocitySpace_.gridPart();
        }

        /**
         *  \brief  get a reference to the discrete velocity function space index set
         *  \todo   doc
         **/
        inline const IndexSetType& indexSet() const
        {
            return velocitySpace_.indexSet();
        }

        /**
         *  \brief  get number of DoFs for the discrete velocity function space
         *  \todo   doc
         **/
        inline int size() const
        {
            return velocitySpace_.size();
        }

        /**
         *  \brief  get discrete velocity function space iterator pointing to the first entity of the associated grid partition
         *  \todo   doc
         **/
        inline IteratorType begin() const
        {
            return velocitySpace_.begin();
        }

        /**
         *  \brief  get discrete velocity function space iterator pointing behind the last entity of the associated grid partition
         *  \todo   doc
         **/
        inline IteratorType end() const
        {
            return velocitySpace_.end();
        }

        /**
         *  \brief  apply discrete velocity function space functor to each entity in the associated grid partition
         *  \todo   doc
         **/
        template< class FunctorType >
        inline void forEach( FunctorType& f ) const
        {
            return velocitySpace_.forEach( f );
        }

        /**
         *  \brief  returns true if the discrete velocity function space grid has more than one geometry type
         *  \todo   doc
         **/
        inline bool multipleGeometryTypes() const
        {
            return velocitySpace_.multipleGeometryTypes();
        }

        /**
         *  \brief  returns true if discrete velocity fucntion space base function sets depend on the entity
         *  \todo   doc
         **/
        inline bool multipleBaseFunctionSets() const
        {
            return velocitySpace_.multipleBaseFunctionSets();
        }

        /**
         *  \brief  map local DoF number to global DoF number (discrete velocity function space)
         *  \todo   doc
         **/
        inline int mapToGlobal( const EntityType& entity,
                                const int localDof ) const
        {
            return velocitySpace_.mapToGlobal( entity, localDof );
        }

        /**
         *  \brief  return maximal number of local DoFs in discrete velocity funtion space
         *  \todo   doc
         **/
        inline int maxNumLocalDofs() const
        {
            return velocitySpace_.maxNumLocalDofs();
        }

        /**
         *  \brief  creates DataHandle for given discrete function (from dicrete velocity function space)
         *  \todo   doc
         **/
        template< class DiscreteFunction, class Operation >
        inline typename BaseType::template CommDataHandle< DiscreteFunction, Operation >::Type createDataHandle(    DiscreteFunction& discreteFunction,
                                                                                                                    const Operation* operation ) const
        {
            return velocitySpace_.createDataHandle( discreteFunction,
                                                    operation );
        }

    private:

        DiscreteVelocityFunctionSpaceType velocitySpace_;
        DiscretePressureFunctionSpaceType pressureSpace_;


}; // end of DiscreteOseenFunctionSpaceWrapper

//! forward
template < class DiscreteOseenFunctionWrapperTraitsImp >
class DiscreteOseenFunctionWrapper;

/**
 *  \todo   doc
 **/
template <  class DiscreteOseenFunctionSpaceWrapperImp,
            class DiscreteVelocityFunctionImp,
            class DiscretePressureFunctionImp >
class DiscreteOseenFunctionWrapperTraits
{
    public:

        //! own type (CRTP)
        typedef DiscreteOseenFunctionWrapper<
                    DiscreteOseenFunctionWrapperTraits<
                        DiscreteOseenFunctionSpaceWrapperImp,
                        DiscreteVelocityFunctionImp,
                        DiscretePressureFunctionImp > >
            DiscreteFunctionType;

        //! type of associated discrete function space
        typedef DiscreteOseenFunctionSpaceWrapperImp
            DiscreteFunctionSpaceType;

        //! type of discrete velocity function
        typedef DiscreteVelocityFunctionImp
            DiscreteVelocityFunctionType;

        //! type of discrete pressure function
        typedef DiscretePressureFunctionImp
            DiscretePressureFunctionType;

		typedef VTKIO < typename DiscretePressureFunctionType::DiscreteFunctionSpaceType::GridPartType >
			VtkWriterType;

        typedef Dune::tuple<const DiscreteVelocityFunctionType*,const DiscretePressureFunctionType*>
			FunctionTupleType;

}; // end of DiscreteOseenFunctionWrapperTraits

/**
 *  \todo   doc,
 *          should comply with the DiscreteFunctionInterface some time
 **/
template < class DiscreteOseenFunctionWrapperTraitsImp >
class DiscreteOseenFunctionWrapper
//    : public DiscreteFunctionInterface< DiscreteOseenFunctionWrapperTraitsImp >
{
    public:

        typedef DiscreteOseenFunctionWrapperTraitsImp
            Traits;
		typedef DiscreteOseenFunctionWrapper<Traits>
			ThisType;

        typedef typename Traits::DiscreteFunctionType
            DiscreteFunctionType;

        //! DiscreteOseenFunctionSpaceWrapper
        typedef typename Traits::DiscreteFunctionSpaceType
            DiscreteFunctionSpaceType;

        //! type of discrete velocity function
        typedef typename Traits::DiscreteVelocityFunctionType
            DiscreteVelocityFunctionType;

        //! type of discrete pressure function
        typedef typename Traits::DiscretePressureFunctionType
            DiscretePressureFunctionType;

        //! type of range field
        typedef typename DiscreteVelocityFunctionType::RangeFieldType
            RangeFieldType;
		typedef typename DiscreteFunctionSpaceType::GridPartType
			GridPartType;
		typedef typename GridPartType::GridType
			GridType;

//the adaption manager, etc. is not really tested for all gridparts/spaces so we need to be able to disable it
#if ENABLE_ADAPTIVE
    protected:

		typedef DiscreteOseenFunctionWrapperAdaptionManager< ThisType >
            AdaptionManagerType;

    public:
#endif


        /**
         *  \brief  constructor
         *  \todo   doc
         **/
        DiscreteOseenFunctionWrapper(  const std::string name,
                                        DiscreteFunctionSpaceType& space,
										GridPartType& gridPart )
            : space_( space ),
            velocity_( name + std::string("_velocity"), space.discreteVelocitySpace() ),
			pressure_( name + std::string("_pressure"), space.discretePressureSpace() ),
            #if ENABLE_ADAPTIVE
				adaptionManager_ ( gridPart.grid(), *this ),
            #endif
			vtkWriter_( gridPart )
        {}

        /**
         *  \brief      constructor
         *
         *              wraps existing velocity and pressure
         *  \attention  by copying
         **/
		DiscreteOseenFunctionWrapper(  DiscreteFunctionSpaceType& space,
                                        DiscreteVelocityFunctionType& velocity,
                                        DiscretePressureFunctionType& pressure )
            : space_( space ),
            velocity_( velocity ),
			pressure_( pressure ),
            #if ENABLE_ADAPTIVE
				adaptionManager_ ( pressure.space().gridPart().grid(), *this ),
            #endif
			vtkWriter_( pressure.space().gridPart() )
        {}

        /**
         *  \brief  destructor
         *  \todo   doc
         **/
		virtual ~DiscreteOseenFunctionWrapper()
        {}

        /**
         *  \todo   doc
         **/
        const DiscreteVelocityFunctionType& discreteVelocity() const
        {
            return velocity_;
        }

        /**
         *  \todo   doc
         **/
        DiscreteVelocityFunctionType& discreteVelocity()
        {
            return velocity_;
        }

        /**
         *  \todo   doc
         **/
        const DiscretePressureFunctionType& discretePressure() const
        {
            return pressure_;
        }

        /**
         *  \todo   doc
         **/
        DiscretePressureFunctionType& discretePressure()
        {
            return pressure_;
        }

        /**
         *  \brief  obtain the name of the discrete function
         *  \todo   doc
         **/
        inline const std::string& name() const
        {
            return velocity_.name();
        }

//        //! obtain a local function for an entity (read-only)
//        template< class EntityType >
//        inline const LocalFunctionType localFunction( const EntityType& entity ) const
//        {
//            return velocity_.localFunction( entity );
//        }
//
//        //! obtain a local function for an entity
//        template< class EntityType >
//        inline LocalFunctionType localFunction( const EntityType& entity )
//        {
//            return velocity_.localFunction( entity );
//        }

		//! get the grid from velo space (we're assuming same grid for both functions everywhere anyways)
		const typename DiscreteFunctionSpaceType::GridType& grid() const { return velocity_.space().grid(); }

        /**
         *  \brief  set all degrees of freedom to zero
         *  \todo   doc
         **/
        inline void clear()
        {
            velocity_.clear();
        }

        /**
         *  \brief  obtain total number of DoFs
         *  \todo   doc
         **/
        inline int size() const
        {
            return velocity_.size();
        }

        /**
         *  \brief  add another discrete function to this one
         *  \todo   doc
         **/
        inline DiscreteFunctionType& operator+= ( const DiscreteFunctionType& arg )
        {
            velocity_ += arg.discreteVelocity();
			pressure_ += arg.discretePressure();
            return *this;
        }

        /**
         *  \brief  substract all degrees of freedom from given discrete function using the dof iterators
         *  \todo   doc
         **/
        template < class DFType >
        DiscreteFunctionType& operator-=( const DFType& arg )
        {
            velocity_ -= arg.discreteVelocity();
			pressure_ -= arg.discretePressure();
            return *this;
        }

        /**
         *  \brief  multiply all DoFs by a scalar factor
         *  \todo   doc
         **/
        inline DiscreteFunctionType& operator*=( const RangeFieldType& scalar )
        {
            velocity_ *= scalar;
			pressure_ *= scalar;
            return *this;
        }

        /**
         *  \brief  devide all DoFs by a scalar factor
         *  \todo   doc
         **/
        inline DiscreteFunctionType& operator/=( const RangeFieldType& scalar )
        {
            velocity_ /= scalar;
			pressure_ /= scalar;
            return *this;
        }

		/**
         *  \brief  discreteFunction/grid adaption
         *  this will both refine/coarsen all marked grid entities and adapt both member functions to that new grid layout
         **/
        void adapt()
        {
            #if ENABLE_ADAPTIVE
                adaptionManager_.adapt();
            #else
                //output warning
            #endif
        }

		void writeVTK( const std::string& path, const int number_postfix )
		{
			std::stringstream s;
			s << std::setfill('0') << std::setw(6) << number_postfix;
			writeVTK( path, s.str() );
		}

		//! write both wrapped functions to "path/{pressure,velocity}.name()+postfix+.vtk"
		void writeVTK( const std::string& path, const std::string postfix = std::string() )
		{
			if ( DiscreteVelocityFunctionType::FunctionSpaceType::DimRange > 1 ){
				vtkWriter_.addVectorVertexData( velocity_ );
				vtkWriter_.addVectorCellData( velocity_ );
			}
			else {
				vtkWriter_.addVertexData( velocity_ );
				vtkWriter_.addCellData( velocity_ );
			}

			vtkWriter_.write( getPath( velocity_, path, postfix ) );
			vtkWriter_.clear();

			if ( DiscretePressureFunctionType::FunctionSpaceType::DimRange > 1 ){
				vtkWriter_.addVectorVertexData( pressure_ );
				vtkWriter_.addVectorCellData( pressure_ );
			}
			else {
				vtkWriter_.addVertexData( pressure_ );
				vtkWriter_.addCellData( pressure_ );
			}
			vtkWriter_.write( getPath( pressure_, path, postfix ) );
			vtkWriter_.clear();
		}

		typename Traits::FunctionTupleType& functionTuple() const
		{
			static typename Traits::FunctionTupleType tuple( &velocity_, &pressure_ );
			return tuple;
		}

		template < class ContinuousVelocityType, class ContinuousPressureType >
		void projectInto( const ContinuousVelocityType& continuous_velocity, const ContinuousPressureType& continuous_pressure )
		{
			typedef Dune::L2Projection< double, double, ContinuousPressureType, DiscretePressureFunctionType >
				PressureProjection;
			PressureProjection().operator()( continuous_pressure, pressure_ );

			typedef Dune::L2Projection< double, double, ContinuousVelocityType, DiscreteVelocityFunctionType >
				VelocityProjection;
			VelocityProjection().operator()( continuous_velocity, velocity_ );
		}

		void assign( const ThisType& other )
		{
			velocity_.assign( other.discreteVelocity() );
			pressure_.assign( other.discretePressure() );
		}

		const DiscreteFunctionSpaceType& space() const
		{
			return space_;
		}

		DiscreteFunctionSpaceType& space()
		{
			return space_;
		}

    private:
		template <class FunctionType>
		inline const char* getPath( const FunctionType& f, const std::string& base_path, const std::string& postfix ) const
		{
			std::stringstream ss;
			ss << base_path;
			if ( base_path.at(base_path.size()-1) != '/' )
				ss << '/';
			ss << f.name() << postfix;
			return ss.str().c_str();
		}

		DiscreteFunctionSpaceType& space_;
        DiscreteVelocityFunctionType velocity_;
        DiscretePressureFunctionType pressure_;

        //declaration order is important here, do not change
    #if ENABLE_ADAPTIVE
		AdaptionManagerType adaptionManager_;
    #endif
		typename Traits::VtkWriterType vtkWriter_;

		// we can uncomment this if the adpation manager copy-problem is resolved
		//DiscreteOseenFunctionWrapper( const DiscreteOseenFunctionWrapper& );

}; // end of DiscreteOseenFunctionWrapper

/** \brief a minimal function wrapper making a function pair usable in OEM solver size() and grid() methods
  */
template < class DiscreteFunctionWrapper >
class CombinedDiscreteFunctionSpace
{
	public:
		CombinedDiscreteFunctionSpace( const DiscreteFunctionWrapper& function_wrapper )
			:function_wrapper_( function_wrapper )
		{}

		const typename DiscreteFunctionWrapper::DiscreteFunctionSpaceType::GridType& grid() const
		{
			return function_wrapper_.grid();
		}

		size_t size() const
		{
			return function_wrapper_.discreteVelocity().space().size()
				   + function_wrapper_.discretePressure().space().size();
		}

	private:
		const DiscreteFunctionWrapper& function_wrapper_;
};

/** \brief a minimal function wrapper making a function pair usable in OEM solver requiring continous access to dof arrays
  *
  **/
template < class DiscreteFunctionWrapper >
class CombinedDiscreteFunction
{
public:
	typedef CombinedDiscreteFunction< DiscreteFunctionWrapper >
		ThisType;
	typedef typename DiscreteFunctionWrapper::DiscreteVelocityFunctionType
		DiscreteFunctionTypeA;
	typedef typename DiscreteFunctionTypeA::DiscreteFunctionSpaceType
		DiscreteFunctionSpaceTypeA;
	typedef typename DiscreteFunctionTypeA::DofIteratorType
		DofIteratorTypeA;
	typedef typename DiscreteFunctionTypeA::ConstDofIteratorType
		ConstDofIteratorTypeA;

	typedef typename DiscreteFunctionWrapper::DiscretePressureFunctionType
		DiscreteFunctionTypeB;
	typedef typename DiscreteFunctionTypeB::DiscreteFunctionSpaceType
		DiscreteFunctionSpaceTypeB;
	typedef typename DiscreteFunctionTypeB::DofIteratorType
		DofIteratorTypeB;
	typedef typename DiscreteFunctionTypeB::ConstDofIteratorType
		ConstDofIteratorTypeB;

	typedef typename DiscreteFunctionSpaceTypeA::RangeFieldType
		RangeFieldType;
	typedef typename DiscreteFunctionSpaceTypeA::DomainFieldType
		DomainFieldType;
	typedef MutableArray< RangeFieldType >
		DofStorageType;
	typedef CombinedDiscreteFunctionSpace< DiscreteFunctionWrapper >
		FunctionSpaceType;

	CombinedDiscreteFunction( const CombinedDiscreteFunction& other )
		:combined_space_( other.space() ),
		space_A_( other.space_A() ),
		space_B_( other.space_B() ),
		numDofs_( space_A_.size() + space_B_.size() ),
		dofVec_( numDofs_ )
	{
		for ( size_t i = 0; i < numDofs_; ++i )
			dofVec_[i] += other.leakPointer()[i];
	}

	CombinedDiscreteFunction( const DiscreteFunctionWrapper& wrapper )
		:combined_space_( wrapper ),
		space_A_( wrapper.discreteVelocity().space() ),
		space_B_( wrapper.discretePressure().space() ),
		numDofs_( space_A_.size() + space_B_.size() ),
		dofVec_( numDofs_ )
	{
		size_t i = 0;
		for (	ConstDofIteratorTypeA it = wrapper.discreteVelocity().dbegin();
				int(i) < space_A_.size() && wrapper.discreteVelocity().dend() != it;
				++it,++i)
		{
			dofVec_[i] = *it;
		}

		i = space_A_.size();
		for (	ConstDofIteratorTypeB it = wrapper.discretePressure().dbegin();
				int(i) < space_A_.size() + space_B_.size() && wrapper.discretePressure().dend() != it;
				++it,++i )
		{
			dofVec_[i] = *it;
		}
	}


	//! copy dofs back to respective functions
	void copyBack( DiscreteFunctionWrapper& wrapper ) const
	{
		size_t i = 0;
		for (	DofIteratorTypeA it = wrapper.discreteVelocity().dbegin();
				int(i) < space_A_.size() && wrapper.discreteVelocity().dend() != it;
				++it,++i)
		{
			*it = dofVec_[i];
		}

		i = space_A_.size();
		for (	DofIteratorTypeB it = wrapper.discretePressure().dbegin();
				int(i) < space_A_.size() + space_B_.size() && wrapper.discretePressure().dend() != it;
				++it,++i )
		{
			*it = dofVec_[i];
		}
	}

	//! return pointer to internal array for use of BLAS routines
	RangeFieldType * leakPointer () { return dofVec_.leakPointer();  }
	//! return pointer to internal array for use of BLAS routines
	const RangeFieldType * leakPointer () const { return dofVec_.leakPointer(); }

	RangeFieldType * leakPointerB () { return dofVec_.leakPointer() + space_A().size();  }
	//! return pointer to internal array for use of BLAS routines
	const RangeFieldType * leakPointerB () const { return dofVec_.leakPointer() + space_A().size(); }

	const FunctionSpaceType& space() const
	{
		return combined_space_;
	}

	ThisType& operator+= (const ThisType &g)
	{
		for ( size_t i = 0; i < numDofs_; ++i )
			dofVec_[i] += g.leakPointer()[i];
		return *this;
	}

	ThisType& operator-= (const ThisType &g)
	{
		for ( size_t i = 0; i < numDofs_; ++i )
			dofVec_[i] -= g.leakPointer()[i];
		return *this;
	}

	ThisType& operator*= (const RangeFieldType &scalar)
	{
		for ( size_t i = 0; i < numDofs_; ++i )
			dofVec_[i] *= scalar;
		return *this;
	}

	ThisType& operator/= (const RangeFieldType &scalar)
	{
		for ( size_t i = 0; i < numDofs_; ++i )
			dofVec_[i] /= scalar;
		return *this;
	}

private:
	FunctionSpaceType combined_space_;
	const DiscreteFunctionSpaceTypeA& space_A_;
	const DiscreteFunctionSpaceTypeB& space_B_;
	const size_t numDofs_;
	mutable DofStorageType dofVec_;

	const DiscreteFunctionSpaceTypeA& space_A() const { return space_A_; }
	const DiscreteFunctionSpaceTypeB& space_B() const { return space_B_; }

};

} // end of namespace Dune

#endif // end of discretestokesfunctionspacewrapper.hh
