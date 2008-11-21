#ifndef UNITCUBE_YASPGRID_HH
#define UNITCUBE_YASPGRID_HH

#include <dune/grid/yaspgrid.hh>

// YaspGrid specialization
template<int dim, int size>
class UnitCube<Dune::YaspGrid<dim,dim>,size> 
{
public:
  typedef Dune::YaspGrid<dim,dim> GridType;

  UnitCube () : Len(1.0), s(size), p(false),
#if HAVE_MPI
  grid_(MPI_COMM_WORLD,Len,s,p,1)
#else
  grid_(Len,s,p,1)
#endif
  {  }

  Dune::YaspGrid<dim,dim>& grid ()
  {
	return grid_;
  }

private:  
  Dune::FieldVector<double,dim> Len; 
  Dune::FieldVector<int,dim> s;
  Dune::FieldVector<bool,dim> p;
  Dune::YaspGrid<dim,dim> grid_;
};

#endif
