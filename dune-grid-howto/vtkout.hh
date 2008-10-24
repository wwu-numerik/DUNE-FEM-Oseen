#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <stdio.h>

template<class G, class V>
void vtkout (const G& grid, const V& c, char* name, int k)
{
  Dune::VTKWriter<G> vtkwriter(grid);
  char fname[128];
  sprintf(fname,"%s-%05d",name,k);
  vtkwriter.addCellData(c,"celldata");
  vtkwriter.write(fname,Dune::VTKOptions::ascii);
}
