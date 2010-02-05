/**
 *  \file   analyticaldata.hh
 *  \brief  contains classes representing analytical data
 **/

#ifndef ANALYTICALDATA_HH
#define ANALYTICALDATA_HH

#include <cmath>

#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>
#include <dune/stuff/parametercontainer.hh>
#include <dune/stokes/boundaryinfo.hh>

/**
 *  \todo   texdoc
 **/
template < class FunctionSpaceImp >
class Force : public Dune::Function < FunctionSpaceImp , Force < FunctionSpaceImp > >
{
    public:
        typedef Force< FunctionSpaceImp >
            ThisType;
        typedef Dune::Function< FunctionSpaceImp, ThisType >
            BaseType;
        typedef typename BaseType::DomainType
            DomainType;
        typedef typename BaseType::RangeType
            RangeType;

        /**
         *  \brief  constructor
         *  \param  viscosity   viscosity \f$\mu\f$ of the fluid
         **/
        Force( const double viscosity, const FunctionSpaceImp& space, const double alpha = 0.0 )
            : BaseType ( space ),
              viscosity_( viscosity ),
              alpha_( alpha ),
              dim_( FunctionSpaceImp::dimDomain )
        {}

        /**
         *  \brief  destructor
         *  doing nothing
         **/
        ~Force()
        {}

        /**
         *  \brief  evaluates the force
         *  \param  arg
         *          point to evaluate at
         *  \param  ret
         *          value of force at given point
         **/
        inline void evaluate( const DomainType& arg, RangeType& ret ) const
        {
            if ( dim_ == 1 ) {
                assert( !"force not implemented in 1D!" );
            }
            else if ( dim_ == 2 ) {
                const double x1 = arg[0];
                const double x2 = arg[1];
#ifdef SIMPLE_PROBLEM
                ret[0] = 0.0;//arg[1];
                ret[1] = 0.0;//-1.0;//arg[0];
#elif defined(CONSTANT_PROBLEM)
                ret[0] = 0.0;//arg[1];
                ret[1] = -1.0;//arg[0];
#elif defined(ROTATE_PROBLEM)
                ret[0] = arg[1];
                ret[1] = -1.0 * arg[0];
#elif defined(POROSITY_PROBLEM)
                ret[ 0 ] = 0.0;
                ret[ 1 ] = 1.0;
#elif defined(POROSITY_PROBLEM_WOIDS)
                ret[ 0 ] = 0.0;
                ret[ 1 ] = 0.0;
#elif defined(GENRALIZED_STOKES_PROBLEM)
                const double x = arg[0];
                const double y = arg[1];
                const double tmp = alpha_ * std::cos( M_PI_2 * (x+y) ) + M_PI_2 * M_PI * std::cos( M_PI_2 * ( x + y ) ) + M_PI_2 * std::cos( M_PI_2 * ( x - y ) ) ;
                ret[0]  =      tmp;
                ret[1] =    -  tmp;

#elif defined(DARCY_PROBLEM)
                // im verhältnis zu [-1,1]^2
                double scaleX = Parameters().getParam( "domain_scale_x", 2.0 );
                double scaleY = Parameters().getParam( "domain_scale_y", 2.0 );
                double shiftX = Parameters().getParam( "domain_shift_x", -1.0 );
                double shiftY = Parameters().getParam( "domain_shift_y", -1.0 );
                ret[0] = ( arg[1] * scaleY ) + shiftY;
                ret[1] = -1.0 * ( ( arg[0] * scaleX ) + shiftX );
#elif defined(MICRO_PROBLEM_X)
                ret[0] = 1.0;
                ret[1] = 0.0;
#elif defined(MICRO_PROBLEM_Y)
                ret[0] = 0.0;
                ret[1] = 1.0;
#else
                ret[0] = 0.0;//arg[1];
                ret[1] = 0.0;//arg[0];
#endif
            }
            else if ( dim_ == 3 ) {
//                const double x1 = arg[0];
//                const double x2 = arg[1];
//                const double x3 = arg[2];
#ifdef SIMPLE_PROBLEM
                ret[0] = 0.0;//arg[1];
                ret[1] = 0.0;//-1.0;//arg[0];
                ret[2] = 0.0;
#elif defined(CONSTANT_PROBLEM)
                ret[0] = 0.0;//arg[1];
                ret[1] = 0.0;//arg[0];
                ret[2] = -1.0;//arg[0];
#elif defined(ROTATE_PROBLEM)
                assert( !"ROTATE_PROBLEM not implemented in 3D!" );
#elif defined(POROSITY_PROBLEM)
                assert( !"POROSITY_PROBLEM not implemented in 3D!" );
#elif defined(POROSITY_PROBLEM_WOIDS)
                assert( !"POROSITY_PROBLEM_WOIDS not implemented in 3D!" );
#elif defined(GENRALIZED_STOKES_PROBLEM)
                assert( !"GENRALIZED_STOKES_PROBLEM not implemented in 3D!" );
#elif defined(DARCY_PROBLEM)
                assert( !"DARCY_PROBLEM not implemented in 3D!" );
#elif defined(AORTA_PROBLEM)
                ret[0] = 0.0;//arg[1];
                ret[1] = 0.0;//-1.0;//arg[0];
                ret[2] = 0.0;
#else
                assert( !"force not implemented in 3D!" );
#endif
            }
            else {
                assert( !"force not implemented for more then 3 dimensions!" );
            }
        }

    private:
        const double viscosity_;
        const double alpha_;
        const int dim_;
};

/**
 *  \brief  describes the dirichlet boundary data
 *
 *  \tparam DirichletTraitsImp
 *          types like functionspace, range type, etc
 *
 *  \todo   extensive docu with latex
 **/
template < class FunctionSpaceImp >
class DirichletData : public Dune::Function < FunctionSpaceImp, DirichletData < FunctionSpaceImp > >
{
    public:
        typedef DirichletData< FunctionSpaceImp >
            ThisType;
        typedef Dune::Function< FunctionSpaceImp, ThisType >
            BaseType;
        typedef typename BaseType::DomainType
            DomainType;
        typedef typename BaseType::RangeType
            RangeType;

        /**
         *  \brief  constructor
         *
         *  doing nothing besides Base init
         **/
        DirichletData( const FunctionSpaceImp& space )
            : BaseType( space ),
              dim_( FunctionSpaceImp::dimDomain )
        {}

        /**
         *  \brief  destructor
         *
         *  doing nothing
         **/
         ~DirichletData()
         {}

        template < class IntersectionIteratorType >
        void evaluate( const DomainType& arg, RangeType& ret, const IntersectionIteratorType& faceIter ) const
        {
            const int id = faceIter.boundaryId();
            if ( dim_ == 1 ) {
                assert( !"dirichlet data not implemented in 1D!" );
            }
            else if ( dim_ == 2 ) {
                // some computations
                const double x1 = arg[0];
                const double x2 = arg[1];
#ifdef SIMPLE_PROBLEM
                ret[0] = 1.0;
                ret[1] = 0.0;
#elif defined(CONSTANT_PROBLEM)
                ret[0] = 0.0;
                ret[1] = 0.0;
#elif defined(ROTATE_PROBLEM)
                ret[0] = 0.0;
                ret[1] = 0.0;
#elif defined(POROSITY_PROBLEM)
                if ( id == 2 ) { // faces on inner hole
                    ret[ 0 ] = 0.0;
                    ret[ 1 ] = 0.0;
                }
                else if ( id == 3 ) { // bottom faces
                    ret[ 0 ] = 0.0;
                    ret[ 1 ] = 1.0;
                }
                else if ( id == 4 ) { // right faces
                    ret[ 0 ] = 0.0;
                    ret[ 1 ] = 0.0;
                }
                else if ( id == 5 ) { // top faces
                    ret[ 0 ] = 0.0;
                    ret[ 1 ] = 1.0;
                }
                else if ( id == 6 ) { // left faces
                    ret[ 0 ] = 0.0;
                    ret[ 1 ] = 0.0;
                }
#elif defined(POROSITY_PROBLEM_WOIDS)
                const double x1 = arg[0];
                const double x2 = arg[1];
                if ( !( x2 > 0.0 ) ) { // bottom faces
                    ret[ 0 ] = 0.0;
                    ret[ 1 ] = 0.0;
                }
                else if ( !( x1 < 1.0 ) ) { // right faces
                    ret[ 0 ] = 1.0;
                    ret[ 1 ] = 0.0;
                }
                else if ( !( x2 < 1.0 ) ) { // top faces
                    ret[ 0 ] = 0.0;
                    ret[ 1 ] = 0.0;
                }
                else if ( !( x1 > 0.0 ) ) { // left faces
                    ret[ 0 ] = 1.0;
                    ret[ 1 ] = 0.0;
                }
                else {
                    ret[ 0 ] = 0.0;
                    ret[ 1 ] = 0.0;
                }
#elif defined(GENRALIZED_STOKES_PROBLEM)
                const double tmp = std::cos( ( M_PI_2 ) * ( x1 + x2 ) );
                ret[0] = tmp;
                ret[1] = -1.0 * tmp;
#elif defined(DARCY_PROBLEM)
                ret[0] = 0.0;
                ret[1] = 0.0;
#elif defined(MICRO_PROBLEM_X)
                if ( id == 2 ) { // faces on inner hole
                    ret[ 0 ] = 0.0;
                    ret[ 1 ] = 0.0;
                }
                else if ( id == 3 ) { // bottom faces
                    ret[ 0 ] = 0.0;
                    ret[ 1 ] = 0.0;
                }
                else if ( id == 4 ) { // right faces
                    ret[ 0 ] = 1.0;
                    ret[ 1 ] = 0.0;
                }
                else if ( id == 5 ) { // top faces
                    ret[ 0 ] = 0.0;
                    ret[ 1 ] = 0.0;
                }
                else if ( id == 6 ) { // left faces
                    ret[ 0 ] = 1.0;
                    ret[ 1 ] = 0.0;
                }
#elif defined(MICRO_PROBLEM_Y)
                if ( id == 2 ) { // faces on inner hole
                    ret[ 0 ] = 0.0;
                    ret[ 1 ] = 0.0;
                }
                else if ( id == 3 ) { // bottom faces
                    ret[ 0 ] = 0.0;
                    ret[ 1 ] = 1.0;
                }
                else if ( id == 4 ) { // right faces
                    ret[ 0 ] = 0.0;
                    ret[ 1 ] = 0.0;
                }
                else if ( id == 5 ) { // top faces
                    ret[ 0 ] = 0.0;
                    ret[ 1 ] = 1.0;
                }
                else if ( id == 6 ) { // left faces
                    ret[ 0 ] = 0.0;
                    ret[ 1 ] = 0.0;
                }
#else
                if ( id == 2 ) { // faces on inner hole
                    ret[ 0 ] = 0.0;
                    ret[ 1 ] = 0.0;
                }
                else {
                    double exp_of_x1 = std::exp( x1 );
                    double sin_of_x2 = std::sin( x2 );
                    double cos_of_x2 = std::cos( x2 );
                    //return
                    ret[0] = x2 * cos_of_x2;
                    ret[0] += sin_of_x2;
                    ret[0] *= -1.0 * exp_of_x1;
                    ret[1] = exp_of_x1 * x2 * sin_of_x2;
                }
#endif
            }
            else if ( dim_ == 3 ) {
//                double x1 = arg[0];
//                double x2 = arg[1];
//                double x3 = arg[2];
                // some computations
#ifdef SIMPLE_PROBLEM
                ret[0] = 1.0;
                ret[1] = 0.0;
                ret[2] = 0.0;
#elif defined(CONSTANT_PROBLEM)
                ret[0] = 0.0;
                ret[1] = 0.0;
                ret[2] = 0.0;
#elif defined(ROTATE_PROBLEM)
                assert( !"ROTATE_PROBLEM not implemented in 3D!" );
#elif defined(POROSITY_PROBLEM)
                assert( !"POROSITY_PROBLEM not implemented in 3D!" );
#elif defined(POROSITY_PROBLEM_WOIDS)
                assert( !"POROSITY_PROBLEM_WOIDS not implemented in 3D!" );
#elif defined(GENRALIZED_STOKES_PROBLEM)
                assert( !"GENRALIZED_STOKES_PROBLEM not implemented in 3D!" );
#elif defined(DARCY_PROBLEM)
                assert( !"DARCY_PROBLEM not implemented in 3D!" );
#elif defined(AORTA_PROBLEM)

                switch ( id ) {
                    case 1: {
                        ret[0] = 0.0;//arg[1];
                        ret[1] = 0.0;//-1.0;//arg[0];
                        ret[2] = 0.0;
                        return;
                    }
                    case 2: {
                        ret[0] = 1000.0;//arg[1];
                        ret[1] = 1000.0;//-1.0;//arg[0];
                        ret[2] = 1000.0;
                        return;
                    }
                    case 6:
                    case 5:
                    case 4:
                    case 3: {
                        ret[0] = 1000.0;//arg[1];
                        ret[1] = 1000.0;//-1.0;//arg[0];
                        ret[2] = 1000.0;
                        return;
                    }
                    default:
                        assert( false );
                        return;
                }
#else
                assert( !"dirichlet data not implemented in 3D!" );
#endif
            }
            else {
                assert( !"dirichlet data not implemented for more than 3 dimensions!" );
            }
        }

         /**
          * \brief  evaluates the dirichlet data
          * \param  arg
          *         point to evaluate at
          * \param  ret
          *         value of dirichlet boundary data at given point
          **/
        inline void evaluate( const DomainType& arg, RangeType& ret ) const {}

    private:
        const int dim_;

};

struct SimpleDirichletDataTraits {

	template < class FunctionSpaceImp, class GridPartImp >
	struct Implementation {
		typedef DirichletData< FunctionSpaceImp>
			AnalyticalDirichletDataType;

		template <class DiscreteStokesFunctionWrapper >
		static AnalyticalDirichletDataType getInstance( const DiscreteStokesFunctionWrapper& wrapper ) {
			return 	AnalyticalDirichletDataType( wrapper.discreteVelocitySpace() );
		}
	};
};


template < class FunctionSpaceImp >
class InOutFluxDirichletData : public Dune::Function < FunctionSpaceImp, InOutFluxDirichletData < FunctionSpaceImp > >
{
	public:
		enum BoundaryType {
			zeroBoundary	= 0,
			influxBoundary	= 1,
			outfluxBoundary	= 2
		};

		typedef InOutFluxDirichletData< FunctionSpaceImp >
			ThisType;
		typedef Dune::Function< FunctionSpaceImp, ThisType >
			BaseType;
		typedef typename BaseType::DomainType
			DomainType;
		typedef typename BaseType::RangeType
			RangeType;
		typedef std::map< int, BoundaryType >
			BoundaryIdTypeMapType;
		typedef typename BoundaryIdTypeMapType::const_iterator
			BoundaryIdTypeMapTypeConstIterator;

		/**
			*  \brief  constructor
			*
			*
			**/
		InOutFluxDirichletData( const FunctionSpaceImp& space )
			: BaseType( space ),
				dim_( FunctionSpaceImp::dimDomain )
		{
			zeroBoundaryIds_ 	= Parameters().getList( "zeroBoundaryIds" , 1 );
			influxBoundaryIds_	= Parameters().getList( "influxBoundaryIds" , 2 );
			outfluxBoundaryIds_	= Parameters().getList( "outfluxBoundaryIds" , 3 );
			setupBoundaryIdTypeMap_();
		}

		/**
			*  \brief  destructor
			*
			*  doing nothing
			**/
		~InOutFluxDirichletData()
		{}

		template < class IntersectionIteratorType >
		void evaluate( const DomainType& arg, RangeType& ret, const IntersectionIteratorType& faceIter ) const
		{
			const int id = faceIter.boundaryId();
		#if defined(AORTA_PROBLEM)


			typedef Dune::FieldVector< typename IntersectionIteratorType::ctype, IntersectionIteratorType::dimension - 1 >
				LocalVectorType;

			LocalVectorType center = Stuff::getBarycenterLocal( faceIter.intersectionSelfLocal() );
			RangeType normal = faceIter.unitOuterNormal( center );
			static const double gd_factor = Parameters().getParam( "gd_factor", 1.0 );
			ret = normal;
			BoundaryIdTypeMapTypeConstIterator id_it = boundaryIdTypeMap_.find( id );
			assert ( id_it != boundaryIdTypeMap_.end() );

			switch ( id_it->second ) {
				case zeroBoundary: {
					ret *= 0;
					return;
				}
				case influxBoundary:
					ret *= -1;
				case outfluxBoundary: {
					ret *= gd_factor;
					return;
				}
				default:
					assert( false );
					return;
			}
		#else
			ASSERT_EXCEPTION( false, "only valid in aorta problem" );
		#endif
		}

		inline void evaluate( const DomainType& arg, RangeType& ret ) const { assert(false); }

	protected:
		void setupBoundaryIdTypeMap_() {
			Logger().Info() << "\t- Using ids ";
			for( std::vector< int >::const_iterator it = zeroBoundaryIds_.begin(); it != zeroBoundaryIds_.end(); ++it ) {
				boundaryIdTypeMap_[*it] = zeroBoundary;
				Logger().Info() << *it << " ";
			}
			Logger().Info() << " \t\tfor g_d = 0 \n\t        ids ";
			for( std::vector< int >::const_iterator it = influxBoundaryIds_.begin(); it != influxBoundaryIds_.end(); ++it ) {
				boundaryIdTypeMap_[*it] = influxBoundary;
				Logger().Info() << *it << " ";
			}
			Logger().Info() << " \tfor g_d = - n \n\t        ids ";
			for( std::vector< int >::const_iterator it = outfluxBoundaryIds_.begin(); it != outfluxBoundaryIds_.end(); ++it ) {
				boundaryIdTypeMap_[*it] = outfluxBoundary;
				Logger().Info() << *it << " ";
			}
			Logger().Info() << " \tfor g_d = n " << std::endl;
		}

		std::vector< int > zeroBoundaryIds_;
		std::vector< int > influxBoundaryIds_;
		std::vector< int > outfluxBoundaryIds_;
		BoundaryIdTypeMapType boundaryIdTypeMap_;

		const int dim_;
};


template < class FunctionSpaceImp, class GridPartType >
class FirstOrderBoundaryShapeFunction : public Dune::BoundaryShapeFunctionBase< FunctionSpaceImp, GridPartType >
{
	typedef Dune::BoundaryShapeFunctionBase < FunctionSpaceImp, GridPartType >
        ParentType;

    public:
        using ParentType::p1;
        using ParentType::p2;
        using ParentType::m;
        using ParentType::scale_factor_;
        using ParentType::direction_;
		using ParentType::edge_distance_;


		FirstOrderBoundaryShapeFunction( const FunctionSpaceImp& space, typename ParentType::PointInfo pf, typename ParentType::RangeType direction, double scale_factor = 1.0 )
            : ParentType ( space, pf, direction, scale_factor )
        {}

        virtual void evaluate( const typename ParentType::DomainType& arg, typename ParentType::RangeType& ret ) const
        {
            typename ParentType::RangeType tmp = direction_;
			const double dist  = ( arg - m ).two_norm();
			tmp *= scale_factor_ * std::abs( edge_distance_ - dist ) / edge_distance_;
            ret = tmp;
        }
};

template < class FunctionSpaceImp, class GridPartType >
class SecondOrderBoundaryShapeFunction : public Dune::BoundaryShapeFunctionBase< FunctionSpaceImp, GridPartType >
{
	typedef Dune::BoundaryShapeFunctionBase < FunctionSpaceImp, GridPartType >
		ParentType;

	public:
		using ParentType::p1;
		using ParentType::p2;
		using ParentType::m;
		using ParentType::scale_factor_;
		using ParentType::direction_;
		using ParentType::edge_distance_;


		SecondOrderBoundaryShapeFunction( const FunctionSpaceImp& space, typename ParentType::PointInfo pf, typename ParentType::RangeType direction, double scale_factor = 1.0 )
			: ParentType ( space, pf, direction, scale_factor ),
			parabolic_stretch_( Parameters().getParam("parabolic_stretch", 1.0) )
		{
		}

		virtual void evaluate( const typename ParentType::DomainType& arg, typename ParentType::RangeType& ret ) const
		{
			typename ParentType::RangeType tmp = direction_;
			const double fac = parabolic_stretch_ * scale_factor_ * ( (arg - p1).two_norm() * (arg - p2).two_norm() );
			tmp *= fac;
			ret = tmp;
		}

	protected:
		const double parabolic_stretch_;
};

template < class FunctionSpaceImp, class GridPartType, template <class ,class> class BoundaryFunctionImp =  FirstOrderBoundaryShapeFunction >
class VariableDirichletData : public Dune::Function < FunctionSpaceImp, VariableDirichletData < FunctionSpaceImp, GridPartType > >
{
	public:
		enum BoundaryType {
			zeroBoundary	= 0,
			influxBoundary	= 1,
			outfluxBoundary	= 2
		};

		typedef VariableDirichletData< FunctionSpaceImp, GridPartType >
			ThisType;
		typedef Dune::Function< FunctionSpaceImp, ThisType >
			BaseType;
		typedef typename BaseType::DomainType
			DomainType;
		typedef typename BaseType::RangeType
			RangeType;
		typedef std::map< int, BoundaryType >
			BoundaryIdTypeMapType;
		typedef typename BoundaryIdTypeMapType::const_iterator
			BoundaryIdTypeMapTypeConstIterator;
        typedef std::map< int, RangeType>
            ID_ValueMapType;
        typedef Dune::BoundaryInfo< GridPartType >
            BoundaryInfoType;
		typedef BoundaryFunctionImp< FunctionSpaceImp, GridPartType >
            BoundaryFunctionType;
        typedef std::vector<BoundaryFunctionType*>
            BoundaryFunctionListType;
		/**
			*  \brief  constructor
			*
			*
			**/
		VariableDirichletData( const FunctionSpaceImp& space, const GridPartType& gridpart )
			: BaseType( space ),
            gridpart_( gridpart ),
            dim_( FunctionSpaceImp::dimDomain ),
            boundaryInfo_( gridpart )
		{
			zeroBoundaryIds_ 	= Parameters().getList( "zeroBoundaryIds" , 1 );
			influxBoundaryIds_	= Parameters().getList( "influxBoundaryIds" , 2 );
			outfluxBoundaryIds_	= Parameters().getList( "outfluxBoundaryIds" , 3 );
			setupBoundaryIdTypeMap_();
		}

		/**
			*  \brief  destructor
			*
			*  doing nothing
			**/
		~VariableDirichletData()
		{}

		template < class IntersectionIteratorType >
		void evaluate( const DomainType& arg, RangeType& ret, const IntersectionIteratorType& faceIter ) const
		{
		#if defined(AORTA_PROBLEM)
			const int id = faceIter.boundaryId();
//			assert( boundaryFunctionList_.size() >= id );
			const BoundaryFunctionType* boundaryFunction = boundaryFunctionList_[id];
			assert ( boundaryFunction );
            boundaryFunction->evaluate( arg, ret );
//			if ( zeroBoundaryIds_.find(id) != zeroBoundaryIds_.end() )
//				ret *=0;
		#else
			ASSERT_EXCEPTION( false, "only valid in aorta problem" );
		#endif
		}

		inline void evaluate( const DomainType& arg, RangeType& ret ) const { assert(false); }

	protected:
		void setupBoundaryIdTypeMap_() {
			Logger().Info() << "\t- Using ids ";
			for( std::vector< int >::const_iterator it = zeroBoundaryIds_.begin(); it != zeroBoundaryIds_.end(); ++it ) {
				boundaryIdTypeMap_[*it] = zeroBoundary;
				Logger().Info() << *it << " ";
			}
			Logger().Info() << " \t\tfor g_d = 0 \n\t        ids ";
			for( std::vector< int >::const_iterator it = influxBoundaryIds_.begin(); it != influxBoundaryIds_.end(); ++it ) {
				boundaryIdTypeMap_[*it] = influxBoundary;
				Logger().Info() << *it << " ";
			}
			Logger().Info() << " \tfor g_d = - n \n\t        ids ";
			for( std::vector< int >::const_iterator it = outfluxBoundaryIds_.begin(); it != outfluxBoundaryIds_.end(); ++it ) {
				boundaryIdTypeMap_[*it] = outfluxBoundary;
				Logger().Info() << *it << " ";
			}
			Logger().Info() << " \tfor g_d = n " << std::endl;
			const double gd_factor = Parameters().getParam( "gd_factor", 1.0 );
			boundaryFunctionList_.reserve( boundaryIdTypeMap_.size() + 1 ); //+1 since BIDs are not 0 indexed
			for ( BoundaryIdTypeMapTypeConstIterator it = boundaryIdTypeMap_.begin(); it != boundaryIdTypeMap_.end(); ++it ) {
                const int id = it->first;
                std::string paramname = std::string( "gd_" ) + Stuff::toString( id );
                std::vector< double > components = Parameters().getList( paramname, double(id) ); //senseless def val here...
                RangeType value;
                assert( components.size() == value.dim() );
                for ( int i = 0; i < value.dim(); ++i )
                    value[i] = components[i];
                id_value_map_[ id ] = value;
                boundaryFunctionList_[id] = new BoundaryFunctionType(BaseType::space(), boundaryInfo_.GetPointInfo( id ),  value, gd_factor );
			}
		}

		std::vector< int > zeroBoundaryIds_;
		std::vector< int > influxBoundaryIds_;
		std::vector< int > outfluxBoundaryIds_;
		BoundaryIdTypeMapType boundaryIdTypeMap_;
        ID_ValueMapType	id_value_map_;
        const GridPartType& gridpart_;

		const int dim_;
		BoundaryInfoType boundaryInfo_;
		BoundaryFunctionListType boundaryFunctionList_;
};

#endif // end of analyticaldata.hh
