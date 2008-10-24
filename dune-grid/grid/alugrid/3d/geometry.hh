#ifndef DUNE_ALU3DGRIDGEOMETRY_HH
#define DUNE_ALU3DGRIDGEOMETRY_HH

// System includes

// Dune includes
#include <dune/common/misc.hh>
#include <dune/grid/common/grid.hh>

// Local includes
#include "alu3dinclude.hh"
#include "topology.hh"

namespace Dune {
  // Forward declarations
  template<int cd, int dim, class GridImp> 
  class ALU3dGridEntity;
  template<int cd, class GridImp >
  class ALU3dGridEntityPointer;
  template<int mydim, int coorddim, class GridImp>
  class ALU3dGridGeometry;
  template<int dim, int dimworld, ALU3dGridElementType elType> 
  class ALU3dGrid;
  class BilinearSurfaceMapping;
  class TrilinearMapping;

  //**********************************************************************
  //
  // --ALU3dGridGeometry
  // --Geometry
  /*!
    Defines the geometry part of a mesh entity. Works for all dimensions, element types and dimensions
    of world. Provides reference element and mapping between local and global coordinates.
    The element may have different implementations because the mapping can be
    done more efficient for structured meshes than for unstructured meshes.

    dim: An element is a polygonal in a hyperplane of dimension dim. 0 <= dim <= 3 is typically 
    dim=0 is a point.

    dimworld: Each corner is a point with dimworld coordinates.
  */
  //! ALU3dGridGeometry 
  //! Empty definition, needs to be specialized for element type
  template <int mydim, int cdim, class GridImp>
  class ALU3dGridGeometry : 
    public GeometryDefaultImplementation <mydim,cdim,GridImp,ALU3dGridGeometry> {};

  //! Specialisation for tetrahedra
  template <int mydim, int cdim>
  class ALU3dGridGeometry<mydim, cdim, const ALU3dGrid<3, 3, tetra> > : 
    public GeometryDefaultImplementation<mydim, cdim, const ALU3dGrid<3, 3, tetra>, 
                           ALU3dGridGeometry> {
    typedef const ALU3dGrid<3, 3, tetra> GridImp;
  
    typedef typename ALU3dImplTraits<tetra>::IMPLElementType IMPLElementType;
    typedef typename ALU3dImplTraits<tetra>::PLLBndFaceType PLLBndFaceType;
    typedef typename ALU3dImplTraits<tetra>::GEOFaceType GEOFaceType;
    typedef typename ALU3dImplTraits<tetra>::GEOEdgeType GEOEdgeType;
    typedef typename ALU3dImplTraits<tetra>::GEOVertexType GEOVertexType;

    typedef ElementTopologyMapping<tetra> ElementTopo;
    typedef FaceTopologyMapping<tetra> FaceTopo;
    //! know number of vertices 
    enum { corners_ = mydim+1 };
  public:
    typedef FieldMatrix<alu3d_ctype, 3, 3> FaceCoordinatesType;
    //! for makeRefGeometry == true a Geometry with the coordinates of the 
    //! reference element is made 
    inline ALU3dGridGeometry();

    //! return the element type identifier
    //! line , triangle or tetrahedron, depends on dim 
    inline const GeometryType & type () const;

    //! return the number of corners of this element. Corners are numbered 0...n-1
    inline int corners () const;

    //! access to coordinates of corners. Index is the number of the corner 
    inline const FieldVector<alu3d_ctype, cdim>& operator[] (int i) const;

    //! maps a local coordinate within reference element to 
    //! global coordinate in element 
    inline FieldVector<alu3d_ctype, cdim> global (const FieldVector<alu3d_ctype, mydim>& local) const;
  
    //! maps a global coordinate within the element to a 
    //! local coordinate in its reference element
    inline FieldVector<alu3d_ctype,  mydim> local (const FieldVector<alu3d_ctype, cdim>& global) const;
  
    //! returns true if the point in local coordinates is inside reference element
    inline bool checkInside(const FieldVector<alu3d_ctype, mydim>& local) const;

    //! A(l) , see grid.hh 
    inline alu3d_ctype integrationElement (const FieldVector<alu3d_ctype, mydim>& local) const;

    //! can only be called for dim=dimworld!
    inline const FieldMatrix<alu3d_ctype,mydim,mydim>& jacobianInverseTransposed (const FieldVector<alu3d_ctype, mydim>& local) const;

    //! returns volume of geometry 
    inline alu3d_ctype volume () const;
    //***********************************************************************
    //!  Methods that not belong to the Interface, but have to be public
    //***********************************************************************
    //! generate the geometry for out of given ALU3dGridElement  
    inline bool buildGeom(const IMPLElementType & item);
    inline bool buildGeom(const ALU3DSPACE HFaceType & item, int twist, int face );
    inline bool buildGeom(const FaceCoordinatesType& coords);
    inline bool buildGeom(const ALU3DSPACE HEdgeType & item, int twist, int );
    inline bool buildGeom(const ALU3DSPACE VertexType & item, int twist, int);
 
    //! build geometry of local coordinates relative to father 
    template <class GeometryImpType>
    inline bool buildGeomInFather(const GeometryImpType &fatherGeom , const GeometryImpType & myGeom);
        
    //! print internal data
    //! no interface method
    inline void print (std::ostream& ss) const;

  private:
    // calculates determinant and volume 
    void calculateDeterminant () const ;
    
    //! calculates the vertex index in the reference element out of a face index
    //! and a local vertex index
    inline int faceIndex(int faceIdx, int vtxIdx) const;

    // generate transposed Jacobian Inverse and calculate integration_element 
    inline void buildJacobianInverseTransposed() const; 

    // calculates the element matrix for calculation of the jacobian inverse
    inline void calcElMatrix () const;

    // copies the values of point to the values of coord 
    inline void copyCoordVec (const alu3d_ctype (& point)[cdim] ,
            FieldVector<alu3d_ctype,cdim> & coord ) const;

    //! the vertex coordinates 
    mutable FieldMatrix<alu3d_ctype, corners_, cdim> coord_; 
  
    mutable FieldMatrix<alu3d_ctype,mydim,mydim> Jinv_; //!< storage for inverse of jacobian
    mutable alu3d_ctype detDF_;              //!< storage of integration_element
    mutable alu3d_ctype volume_;             //!< storage of volume 
    mutable FieldMatrix<alu3d_ctype, cdim , mydim> A_;    //!< transformation matrix (transposed) 
    mutable FieldMatrix<alu3d_ctype, mydim, mydim> AT_A_;    //!< transformation matrix (transposed) 

    mutable FieldVector<alu3d_ctype, mydim> AT_x_;
    mutable FieldVector<alu3d_ctype, mydim> localCoord_;
    mutable FieldVector<alu3d_ctype, cdim>  globalCoord_;

    mutable FieldVector<alu3d_ctype,cdim> tmpV_; //! temporary memory 
    mutable FieldVector<alu3d_ctype,cdim> tmpU_; //! temporary memory

    const GeometryType myGeomType_;

    //! is true if Jinv_, A and detDF_ is calced
    mutable bool builtinverse_;
    mutable bool builtA_;
    mutable bool builtDetDF_;
  };

  //! Specialisation for hexahedra
  template <int mydim, int cdim>
  class ALU3dGridGeometry<mydim, cdim, const ALU3dGrid<3, 3, hexa> > : 
    public GeometryDefaultImplementation<mydim, cdim, const ALU3dGrid<3, 3, hexa>,
                           ALU3dGridGeometry> {
    typedef const ALU3dGrid<3, 3, hexa> GridImp;
    friend class ALU3dGridIntersectionIterator<GridImp>;

    typedef typename ALU3dImplTraits<hexa>::IMPLElementType IMPLElementType;
    typedef typename ALU3dImplTraits<hexa>::PLLBndFaceType PLLBndFaceType;
    typedef typename ALU3dImplTraits<hexa>::GEOFaceType GEOFaceType;
    typedef typename ALU3dImplTraits<hexa>::GEOEdgeType GEOEdgeType;
    typedef typename ALU3dImplTraits<hexa>::GEOVertexType GEOVertexType;

    typedef ElementTopologyMapping<hexa> ElementTopo;
    typedef FaceTopologyMapping<hexa> FaceTopo;

    template <int dummy, int dim>
    struct MappingChooser
    {
      typedef bool TrilinearMappingType;
      typedef bool BilinearSurfaceMappingType; 
    };

    template <int dummy>
    struct MappingChooser<dummy,2>
    {
      typedef bool TrilinearMappingType;
      typedef BilinearSurfaceMapping  BilinearSurfaceMappingType; 
    };

    template <int dummy>
    struct MappingChooser<dummy,3>
    {
      typedef TrilinearMapping TrilinearMappingType;
      typedef bool BilinearSurfaceMappingType; 
    };

    typedef typename MappingChooser<0,mydim> :: TrilinearMappingType TrilinearMappingType;
    typedef typename MappingChooser<0,mydim> :: BilinearSurfaceMappingType BilinearSurfaceMappingType;

    enum { corners_ = Power_m_p<2,mydim>::power };
  public:
    typedef FieldMatrix<alu3d_ctype, 4, 3> FaceCoordinatesType;

    //! for makeRefGeometry == true a Geometry with the coordinates of the 
    //! reference element is made 
    ALU3dGridGeometry();

    //! Destructor
    ~ALU3dGridGeometry();

    //! return the element type identifier
    //! line , triangle or tetrahedron, depends on dim 
    const GeometryType & type () const;

    //! return the number of corners of this element. Corners are numbered 0..n-1
    int corners () const;
  
    //! access to coordinates of corners. Index is the number of the corner 
    const FieldVector<alu3d_ctype, cdim>& operator[] (int i) const;

    //! maps a local coordinate within reference element to 
    //! global coordinate in element 
    FieldVector<alu3d_ctype, cdim> global (const FieldVector<alu3d_ctype, mydim>& local) const;
  
    //! maps a global coordinate within the element to a 
    //! local coordinate in its reference element
    FieldVector<alu3d_ctype,  mydim> local (const FieldVector<alu3d_ctype, cdim>& global) const;
  
    //! returns true if the point in local coordinates is inside reference 
    //! element
    bool checkInside(const FieldVector<alu3d_ctype, mydim>& local) const;

    //! A(l) , see grid.hh 
    alu3d_ctype integrationElement (const FieldVector<alu3d_ctype, mydim>& local) const;

    //! can only be called for dim=dimworld! (Trivially true, since there is no
    //! other specialization...)
    const FieldMatrix<alu3d_ctype,mydim,mydim>& jacobianInverseTransposed (const FieldVector<alu3d_ctype, mydim>& local) const;

    //! returns volume of geometry 
    alu3d_ctype volume () const;

    //***********************************************************************
    //!  Methods that not belong to the Interface, but have to be public
    //***********************************************************************
    //! generate the geometry out of a given ALU3dGridElement  
    bool buildGeom(const IMPLElementType & item);
    bool buildGeom(const ALU3DSPACE HFaceType & item, int twist, int faceNum);
    bool buildGeom(const FaceCoordinatesType& coords);
    bool buildGeom(const ALU3DSPACE HEdgeType & item, int twist, int);
    bool buildGeom(const ALU3DSPACE VertexType & item, int twist, int);
 
    //! build ghost out of internal boundary segment 
    bool buildGhost(const PLLBndFaceType & ghost);
 
    //! build geometry of local coordinates relative to father 
    template <class GeometryType>
    bool buildGeomInFather(const GeometryType &fatherGeom , const GeometryType & myGeom);
        
    //! print internal data
    //! no interface method
    void print (std::ostream& ss) const;

  private:
    // create triMap from coordinates , deletes old mapping 
    void buildMapping() const;
    
    // create biMap_ from coordinates , deletes old mapping 
    void buildBilinearMapping() const;
    
    //! the vertex coordinates 
    mutable FieldMatrix<alu3d_ctype, corners_ , cdim> coord_;
    
    //! the vertex coordinates 
    typedef double CoordPtrType[cdim];
    mutable FieldVector<const CoordPtrType*, corners_ > coordPtr_;
    
    mutable FieldVector<alu3d_ctype, mydim> tmp1_;
    mutable FieldVector<alu3d_ctype, cdim> tmp2_;

    const GeometryType myGeomType_;

    mutable TrilinearMappingType triMap_;
    mutable BilinearSurfaceMappingType biMap_;

    mutable bool buildTriMap_;
    mutable bool buildBiMap_;

    mutable FieldVector<alu3d_ctype, mydim> localBaryCenter_;
    mutable alu3d_ctype volume_;
  };

} // end namespace Dune

#include "geometry_imp.cc"

#endif