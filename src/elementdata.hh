#ifndef DUNE_OSEEN_ELEMENTDATA_HH
#define DUNE_OSEEN_ELEMENTDATA_HH

#include <dune/stuff/common/filesystem.hh>
#include <dune/stuff/aliases.hh>

#include<dune/geometry/genericgeometry/referenceelements.hh>
#include<dune/grid/common/mcmgmapper.hh>
#include<dune/grid/io/file/vtk/vtkwriter.hh>
#if HAVE_GRAPE
#include <dune/grid/io/visual/grapedatadisplay.hh>
#endif

//! Parameter for mapper class
template<int dim>
struct P0Layout
{
  bool contains (Dune::GeometryType gt)
  {
	if (gt.dim()==dim) return true;
	return false;
  }
};

// demonstrate attaching data to elements
template<class G, class F>
void elementdata (const G& grid, const F& f)
{
  // the usual stuff
  typedef typename G::ctype ct;
  typedef typename G::LeafGridView GridView;
  typedef typename GridView::template Codim<0>::Iterator ElementLeafIterator;

  // get grid view on leaf part
  GridView gridView = grid.leafView();

  // make a mapper for codim 0 entities in the leaf grid
  Dune::LeafMultipleCodimMultipleGeomTypeMapper<G,P0Layout>
	  mapper(grid);                                    /*@\label{edh:mapper}@*/

  // allocate a vector for the data
  std::vector<double> c(mapper.size());                /*@\label{edh:c}@*/

  // iterate through all entities of codim 0 at the leafs
  for (ElementLeafIterator it = gridView.template begin<0>(); /*@\label{edh:loop0}@*/
	   it!=gridView.template end<0>(); ++it)
	{
	  // evaluate functor and store value
	  c[mapper.map(*it)] = f(*it);	       /*@\label{edh:feval}@*/
	}                                              /*@\label{edh:loop1}@*/

  // generate a VTK file
  // Dune::LeafP0Function<G,double> cc(grid,c);
  Dune::VTKWriter<typename G::LeafGridView> vtkwriter(gridView); /*@\label{edh:vtk0}@*/
  vtkwriter.addCellData(c,"data");
  DSC::testCreateDirectory( f.filename() );
  vtkwriter.write( f.filename().c_str(), DSC_CONFIG_GET( "binary_vtk", true ) ? Dune::VTK::base64 : Dune::VTK::ascii );

  // online visualization with Grape
#if HAVE_GRAPE                                         /*@\label{edh:grape0}@*/
  if ( DSC_CONFIG_GET("use_grape",false) )
  {
    const int polynomialOrder = 0; // we piecewise constant data
    const int dimRange = 1; // we have scalar data here
    // create instance of data display
    Dune::GrapeDataDisplay<G> grape(grid);
    // display data
    grape.displayVector("concentration", // name of data that appears in grape
                        c,  // data vector
                        grid.leafIndexSet(), // used index set
                        polynomialOrder, // polynomial order of data
                        dimRange); // dimRange of data
  }
#endif                                                 /*@\label{edh:grape1}@*/
}

#endif //DUNE_OSEEN_ELEMENTDATA_HH

/** Copyright (c) 2012, Rene Milk 
 * All rights reserved.
 * 
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met: 
 * 
 * 1. Redistributions of source code must retain the above copyright notice, this
 *    list of conditions and the following disclaimer. 
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution. 
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
 * ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 * 
 * The views and conclusions contained in the software and documentation are those
 * of the authors and should not be interpreted as representing official policies, 
 * either expressed or implied, of the FreeBSD Project.
**/

