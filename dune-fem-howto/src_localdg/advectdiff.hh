/**
 *	\file advectdiff.hh
 *	\brief advectdiff.hh
 **/

#ifndef DUNE_DGOPERATORS_HH
#define DUNE_DGOPERATORS_HH

//- system includes 
#include <string>

// Dune includes
#include <dune/fem/gridpart/gridpart.hh>

#include <dune/fem/pass/dgpass.hh>
#include <dune/fem/pass/discretemodel.hh>
#include <dune/fem/pass/selection.hh>
#include <dune/fem/space/dgspace.hh>
#include <dune/fem/operator/common/spaceoperatorif.hh>
#include "discretemodels.hh"

/****************************************
 * DG Operator                          *
 ***************************************/

namespace Dune {  

/**
 * @brief LDG Operator that is evaluated at every time step by the ODE Solver.
 *
 * LDG Operator \f$\mathcal{L}\f$ that implements the two Passes for the
 * advection-diffusion problem
 */
template <class Model,template<class M> class NumFlux,int polOrd >
class DGAdvectionDiffusionOperator : 
  public SpaceOperatorInterface<
  typename PassTraits<Model, Model::Traits::dimRange, polOrd>::DestinationType> 
  {
public:
  // Id's for the three Passes (including StartPass)
  enum PassIdType{ u, pass1, pass2 };    /*@\label{ad:passids}@*/

  enum { dimRange = Model::dimRange };
  enum { dimDomain = Model::Traits::dimDomain };

  typedef NumFlux<Model> NumFluxType;
  typedef typename Model::Traits::GridType GridType;

  // Pass 1 Model
  typedef AdvDiffDModel1<Model, NumFluxType, polOrd, u>                /*@\label{ad:typedefmodel1}@*/
    DiscreteModel1Type;
  // Pass 2 Model
  typedef AdvDiffDModel2<Model, NumFluxType, polOrd, u, pass1>         /*@\label{ad:typedefmodel2}@*/
    DiscreteModel2Type;

  typedef typename DiscreteModel1Type::Traits Traits1;
  typedef typename DiscreteModel2Type::Traits Traits2;

  typedef typename Traits2::DomainType DomainType;
  typedef typename Traits2::DiscreteFunctionType DiscreteFunction2Type;

  /**************************************
   * Join the Passes 0-2                *
   *************************************/
  typedef StartPass<DiscreteFunction2Type, u> Pass0Type;
  typedef LocalDGPass<DiscreteModel1Type, Pass0Type, pass1> Pass1Type; /*@\label{ad:typedefpass1}@*/
  typedef LocalDGPass<DiscreteModel2Type, Pass1Type, pass2> Pass2Type; /*@\label{ad:typedefpass2}@*/

  typedef typename Traits1::DiscreteFunctionSpaceType 
    Space1Type;
  typedef typename Traits2::DiscreteFunctionSpaceType  
    Space2Type;
  typedef typename Traits1::DestinationType Destination1Type;
  typedef typename Traits2::DestinationType Destination2Type;
  typedef Destination2Type DestinationType;
  typedef Space2Type SpaceType;

  typedef typename Traits1::GridPartType GridPartType;

public:
  /**
   * @brief Constructor
   *
   * initializes the LDGPasses
   *
   * @param grid underlying grid instance
   * @param numf instance of a numerical flux for the advection term
   */
  DGAdvectionDiffusionOperator(GridType& grid,
                               const NumFluxType& numf) :
    grid_(grid),
    model_(numf.model()),
    numflux_(numf),
    gridPart_(grid_),
    space1_(gridPart_),
    space2_(gridPart_),
    problem1_(model_,numflux_),
    problem2_(model_,numflux_),
    pass1_(problem1_, pass0_, space1_),    /*@\label{ad:initialisepass1}@*/
    pass2_(problem2_, pass1_, space2_)     /*@\label{ad:initialisepass2}@*/ 
  {}

  /**
   * @brief This methods gets called by the TimeProvider in order to compute
   * the next time step.
   */
  double timeStepEstimate() const
  {
    return pass2_.timeStepEstimate();
  }

  /**
   * @brief This method is called by the TimeProvider in order to inform the
   * Operator about the current time.
   */
  void setTime(double time) {
    pass2_.setTime(time);
  }

  void operator()(const DestinationType& arg, DestinationType& dest) const {
    pass2_(arg,dest);             /*@\label{ad:callpass2}@*/
  }

  const SpaceType& space() const {
    return space2_;
  }
  
  /**
   * @brief LaTeX Information printed by EocOutput
   */
  void printmyInfo(std::string filename) const {
    std::ostringstream filestream;
    filestream << filename;
    std::ofstream ofs(filestream.str().c_str(), std::ios::app);
    ofs << "Advection-Diffusion Op., polynomial order: " << polOrd << "\\\\\n\n";
    ofs.close();
  }
private:
  GridType& grid_;
  const Model& model_;
  const NumFluxType& numflux_;
  GridPartType gridPart_;
  Space1Type space1_;
  Space2Type space2_;
  DiscreteModel1Type problem1_;
  DiscreteModel2Type problem2_;
  Pass0Type pass0_;
  Pass1Type pass1_;
  Pass2Type pass2_;
};
}
#endif
