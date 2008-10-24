#ifndef DUNE_DGSPACEDATAHANDLE_HH
#define DUNE_DGSPACEDATAHANDLE_HH

#include <dune/fem/misc/utility.hh>
//- Dune includes 
#include <dune/grid/common/datahandleif.hh>

#include <dune/fem/space/common/commoperations.hh>

namespace Dune {

/** @addtogroup DFComm  
    @{
**/

  /** \brief Communication data handle for DiscreteFunctions based on
      discontinuous spaces such as DG or FV spaces.
      \param DiscreteFunctionImp type of discrete function to be
      communicated 
  */
  template <class DiscreteFunctionImp, class OperationImp>
  class DGCommunicationHandler 
   : public CommDataHandleIF<DGCommunicationHandler <DiscreteFunctionImp ,OperationImp> ,
                             typename DiscreteFunctionImp :: RangeFieldType >
  {
    //! empty for higher codims 
    template <class MessageBufferImp, class EntityType, int codim> 
    struct HandleData
    {
      static void gather (DiscreteFunctionImp& discreteFunction,
                          MessageBufferImp& buff, const EntityType& en)
      {
      }
      static void scatter (DiscreteFunctionImp& discreteFunction,
          MessageBufferImp& buff, const EntityType& en, size_t n)
      {
      }
      static size_t size (DiscreteFunctionImp& discreteFunction,const EntityType& en)
      {
        return 0;
      }
    };  
    
    template <class MessageBufferImp, class EntityType> 
    struct HandleData<MessageBufferImp,EntityType,0>
    {
      typedef typename DiscreteFunctionImp :: LocalFunctionType LocalFunctionType; 

      //! gather data 
      static void gather (DiscreteFunctionImp& discreteFunction,
                          MessageBufferImp& buff, const EntityType& en)
      {
        // get local function 
        LocalFunctionType lf = discreteFunction.localFunction(en);
        const int numDofs = lf.numDofs(); 
        // for all local dofs, write data to buffer 
        for(int i=0; i<numDofs; ++i) 
        {
          buff.write( lf[i] );
        }
      }
      
      //! scatter data 
      static void scatter (DiscreteFunctionImp& discreteFunction,
          MessageBufferImp& buff, const EntityType& en, size_t n)
      {
        LocalFunctionType lf = discreteFunction.localFunction(en);
        const int numDofs = lf.numDofs(); 
        DataType val; 
        // for all local dofs, read data from buffer 
        // and apply operation 
        for(int i=0; i<numDofs; ++i) 
        {
          buff.read( val );

          // apply given operation  
          OperationImp::apply(val , lf[i]);
        }
      }
      
      //! return local dof size to be communicated 
      static size_t size (DiscreteFunctionImp& discreteFunction,const EntityType& en)
      {
        // return size of local function 
        LocalFunctionType lf = discreteFunction.localFunction(en);
        return lf.numDofs(); 
      }
    };
    
  public:  
    typedef DiscreteFunctionImp DiscreteFunctionType;
    typedef typename DiscreteFunctionType::LocalFunctionType LocalFunctionType;
    typedef typename DiscreteFunctionType::RangeFieldType DataType;
  private:  
    //! cannot be implemented because of the reference
    DGCommunicationHandler & operator = (const DGCommunicationHandler & org);

    // discrete function to communicate 
    mutable DiscreteFunctionType & discreteFunction_; 

    const bool fixedSize_;
  public:
    DGCommunicationHandler(DiscreteFunctionType & df) 
      : discreteFunction_(df) 
      , fixedSize_ (! df.space().multipleGeometryTypes())
    {
      // if space is continuous, check contained codim again 
      assert( ! df.space().continuous() );
      //std::cout << fixedSize_ << "\n";
    }
    
    DGCommunicationHandler(const DGCommunicationHandler & org)
      : discreteFunction_(org.discreteFunction_) 
      , fixedSize_(org.fixedSize_)
    {
    }

    bool contains (int dim, int codim) const
    {
      return discreteFunction_.space().contains(codim);
    }

    bool fixedsize (int dim, int codim) const
    {
      return fixedSize_;
    }

    //! read buffer and apply operation 
    template<class MessageBufferImp, class EntityType>
    void gather (MessageBufferImp& buff, const EntityType& en) const
    {
      enum { codim = EntityType :: codimension };
      HandleData<MessageBufferImp,EntityType,codim>::gather(discreteFunction_,buff,en);
    }

    //! read buffer and apply operation 
    template<class MessageBufferImp, class EntityType>
    void scatter (MessageBufferImp& buff, const EntityType& en, size_t n)
    {
      enum { codim = EntityType :: codimension };
      HandleData<MessageBufferImp,EntityType,codim>::scatter(discreteFunction_,buff,en,n);
    }

    //! return local dof size to be communicated 
    template <class EntityType>
    size_t size (const EntityType& en) const
    {
      enum { codim = EntityType :: codimension };
      return HandleData<EntityType,EntityType,codim>::size(discreteFunction_,en);
    }
  };
  
//@} 
} // end namespace Dune
#endif