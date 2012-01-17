#ifndef DUNE_OSEEN_STAB_COEFF_HH
#define DUNE_OSEEN_STAB_COEFF_HH

#include <dune/stuff/parametercontainer.hh>
#include <map>

namespace Dune {
class StabilizationCoefficients {
	public:
		typedef int
			PowerType;
		typedef double
			FactorType;

		typedef std::pair<PowerType,FactorType>
			ValueType;

		typedef std::map< std::string, ValueType >
			CoefficientMap;

		static const PowerType invalid_power;
		static const FactorType invalid_factor;

	protected:
		CoefficientMap map_;

	public:
		StabilizationCoefficients(  const PowerType pow,
									const FactorType fac )
		{
			*this = StabilizationCoefficients( pow,pow,pow,pow,fac,fac,fac,fac );
		}

		StabilizationCoefficients(  const PowerType C11_pow,
									const PowerType C12_pow,
									const PowerType D11_pow,
									const PowerType D12_pow,
									const FactorType C11_fac,
									const FactorType C12_fac,
									const FactorType D11_fac,
									const FactorType D12_fac )
		{
			map_["C11"] = ValueType( C11_pow, C11_fac );
			map_["C12"] = ValueType( C12_pow, C12_fac );
			map_["D11"] = ValueType( D11_pow, D11_fac );
			map_["D12"] = ValueType( D12_pow, D12_fac );
		}

		void Add( const std::string& name, FactorType factor, PowerType power = invalid_power )
		{
			if ( map_.find(name) == map_.end() )
				map_[name] = ValueType( power, factor );
		}

		void Add( const std::string& name )
		{
			Add( name, Parameters().getParam( name, FactorType() ), invalid_power );
		}


		FactorType Factor( const std::string& name ) const {
			return getValue(name).second;
		}

		void Factor( const std::string& name, FactorType new_factor ) {
			getValue(name).second = new_factor;
		}

		void FactorFromParams( const std::string& name, const FactorType default_value = FactorType() ) {
			getValue(name).second = Parameters().getParam( name, default_value );
		}

		PowerType Power( const std::string& name ) const {
			return getValue(name).first;
		}

		void Power( const std::string& name, PowerType new_power ) {
			getValue(name).first = new_power;
		}

		static StabilizationCoefficients getDefaultStabilizationCoefficients() {
			return StabilizationCoefficients( -1,   invalid_power,    1,    invalid_power,
											 1.0,               0,  1.0,                0 );
		}

		template < class Stream >
		void print( Stream& stream ) {
			if ( this->Equals( getDefaultStabilizationCoefficients() ) )
				stream << "default stabilisation coefficients used " ;
			else {
				CoefficientMap::const_iterator it = map_.begin();
				for ( ; it != map_.end(); ++it ) {
					stream  << std::endl << (*it).first << " (factor/power): " << (*it).second.second
							<< " / " << (*it).second.first ;
				}
			}
			stream << std::endl;
		}

		std::vector<std::string> getCoefficientNames() const {
			std::vector<std::string> ret;
			for( CoefficientMap::const_iterator it = map_.begin(); it != map_.end(); ++it )
				ret.push_back( it->first );
			return ret;
		}

		bool Equals( const StabilizationCoefficients& other ) {
			if ( map_.size() != other.map_.size() )
				return false;
			return std::equal( map_.begin(), map_.end(), other.map_.begin() );
		}
	private:
		ValueType& getValue( const std::string name )
		{
			CoefficientMap::iterator it = map_.find( name );
			if ( it == map_.end() )
				DUNE_THROW( ParameterInvalid, "Stabilization Parameter '" << name << "' missing. (Did you forget to ::Add(name) it?" );
			return it->second;
		}

		const ValueType& getValue( const std::string name ) const
		{
			CoefficientMap::const_iterator it = map_.find( name );
			if ( it == map_.end() )
				DUNE_THROW( ParameterInvalid, "Stabilization Parameter '" << name << "' missing. (Did you forget to ::Add(name) it?" );
			return it->second;
		}
	public:
		//! gives a vector c such that c * normal = signum( v * n ) / 2
		template < class FieldVectorType >
		class C12 : public FieldVectorType {
		public:
			C12 (	const FieldVectorType& normal,
					const StabilizationCoefficients& coeff,
					const FieldVectorType v = FieldVectorType(1) )
				: FieldVectorType( 0.0 )
			{
				const double v_normal = v * normal;
				if ( v_normal != 0.0 ) {
					const double a = copysign(1.0, v * normal ) / 2.0;
					if ( normal[1] != 0 ) {
						(*this)[0] = 1;
						(*this)[1] = ( a - normal[0])/normal[1];
					}
					else {
						(*this)[1] = 1;
						(*this)[0] = ( a - normal[1])/normal[0];
					}
				}
				*this *= coeff.Factor("C12");
			}
		};
};
const StabilizationCoefficients::PowerType StabilizationCoefficients::invalid_power = -9;
const StabilizationCoefficients::FactorType StabilizationCoefficients::invalid_factor = -9.0;

} //namespace

#endif // DUNE_OSEEN_STAB_COEFF_HH
