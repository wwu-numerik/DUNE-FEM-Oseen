#ifndef UNITCUBE_ONEDGRID_HH
#define UNITCUBE_ONEDGRID_HH

#include<dune/grid/onedgrid.hh>

// OneDGrid specialization
template<>
class UnitCube<Dune::OneDGrid,1> 
{
public:
  typedef Dune::OneDGrid GridType;

  UnitCube () : grid_(1,0.0,1.0)
  {}

  Dune::OneDGrid& grid ()
  {
	return grid_;
  }

private:  
  Dune::OneDGrid grid_;
};

#endif
