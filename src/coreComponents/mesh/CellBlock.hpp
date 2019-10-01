/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2019 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2019 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2019 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All right reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file CellBlock.hpp
 */

#ifndef CELLBLOCK_HPP_
#define CELLBLOCK_HPP_

#include "ElementSubRegionBase.hpp"
#include "FaceManager.hpp"
#include "meshUtilities/ComputationalGeometry.hpp"
#include "rajaInterface/GEOS_RAJA_Interface.hpp"

#define USE_ELEM_PATCHES 1

#if USE_ELEM_PATCHES
  #define ELEM_PATCH_VIZ 1
#endif

class StableTimeStep;

namespace geosx
{

/**
 * Class to manage the data stored at the element level.
 */
class CellBlock : public ElementSubRegionBase
{
public:

  using NodeMapType=FixedOneToManyRelation;
  using FaceMapType=FixedOneToManyRelation;

  /**
   * @name Static Factory Catalog Functions
   */
  ///@{

  /**
   *
   * @return the name of this type in the catalog
   */
  static const string CatalogName()
  { return "CellBlock"; }

  /**
   *
   * @return the name of this type in the catalog
   */
  virtual const string getCatalogName() const override final
  { return CellBlock::CatalogName(); }


  ///@}


  /// deleted default constructor
  CellBlock() = delete;

  /**
   * @brief constructor
   * @param name the name of the object in the data repository
   * @param parent the parent object of this object in the data repository
   */
  CellBlock( string const & name, Group * const parent );

  /**
   * @brief copy constructor
   * @param init the source to copy
   */
  CellBlock(const CellBlock& init);


  virtual ~CellBlock() override;

  virtual void SetElementType( string const & elementType ) override;

  localIndex GetNumFaceNodes( localIndex const elementIndex,
                              localIndex const localFaceIndex) const;

  localIndex GetFaceNodes( localIndex const elementIndex,
                           localIndex const localFaceIndex,
                           localIndex * const nodeIndicies) const;

  /**
   * @brief function to return the localIndices of the nodes in a face of the element
   * @param elementIndex The localIndex of the target element
   * @param localFaceIndex the element local localIndex of the target face (this will be 0-numFacesInElement
   * @param nodeIndicies the node indices of the face
   */
  void GetFaceNodes( const localIndex elementIndex,
                     const localIndex localFaceIndex,
                     localIndex_array& nodeIndicies) const;

  /**
   * @brief function to return element center. this should be depricated.
   * @param k
   * @param nodeManager
   * @param useReferencePos
   * @return
   */
  R1Tensor const & calculateElementCenter( localIndex k,
                                           const NodeManager& nodeManager,
                                           const bool useReferencePos = true) const override;

  void calculateElementCenters( arrayView1d<R1Tensor const> const & X ) const
  {
    arrayView1d<R1Tensor> const & elementCenters = m_elementCenter;
    localIndex nNodes = numNodesPerElement();

    if (!m_elementTypeString.compare(0, 4, "C3D6"))
    {
      nNodes -= 2;
    }

    forall_in_range<parallelHostPolicy>( 0, size(), GEOSX_LAMBDA( localIndex const k )
    {
      elementCenters[k] = 0;
      for ( localIndex a = 0 ; a < nNodes ; ++a)
      {
        const localIndex b = m_toNodesRelation[k][a];
        elementCenters[k] += X[b];
      }
      elementCenters[k] /= nNodes;
    });
  }

  virtual void CalculateElementGeometricQuantities( NodeManager const & nodeManager,
                                                    FaceManager const & facemanager ) override;

  inline void CalculateCellVolumesKernel( localIndex const k,
                                          array1d<R1Tensor> const & X ) const
  {
    R1Tensor & center = m_elementCenter[k];
    center = 0.0;

    R1Tensor Xlocal[10];

    for (localIndex a = 0; a < m_numNodesPerElement; ++a)
    {
      Xlocal[a] = X[m_toNodesRelation[k][a]];
      center += Xlocal[a];
    }
    center /= m_numNodesPerElement;

    if( m_numNodesPerElement == 8 )
    {
      m_elementVolume[k] = computationalGeometry::HexVolume(Xlocal);
    }
    else if( m_numNodesPerElement == 4)
    {
      m_elementVolume[k] = computationalGeometry::TetVolume(Xlocal);
    }
    else if( m_numNodesPerElement == 6)
    {
      m_elementVolume[k] = computationalGeometry::WedgeVolume(Xlocal);
    }
    else if ( m_numNodesPerElement == 5)
    {
      m_elementVolume[k] = computationalGeometry::PyramidVolume(Xlocal);
    }
    else
    {
        GEOS_ERROR("GEOX does not support cells with " << m_numNodesPerElement << " nodes");
    }
  }


  virtual void setupRelatedObjectsInRelations( MeshLevel const * const mesh ) override;


  virtual arraySlice1dRval<localIndex const> nodeList( localIndex const k ) const override
  {
    return m_toNodesRelation[k];
  }

  virtual arraySlice1dRval<localIndex> nodeList( localIndex const k ) override
  {
    return m_toNodesRelation[k];
  }


  /**
   * @return the element to node map
   */
  FixedOneToManyRelation & nodeList()                    { return m_toNodesRelation; }

  /**
   * @return the element to node map
   */
  FixedOneToManyRelation const & nodeList() const        { return m_toNodesRelation; }

  /**
   * @return the element to node map
   */
  localIndex & nodeList( localIndex const k, localIndex a ) { return m_toNodesRelation[k][a]; }

  /**
   * @return the element to node map
   */
  localIndex const & nodeList( localIndex const k, localIndex a ) const { return m_toNodesRelation[k][a]; }

  /**
   * @return the element to face map
   */
  FixedOneToManyRelation       & faceList()       { return m_toFacesRelation; }

  /**
   * @return the element to face map
   */
  FixedOneToManyRelation const & faceList() const { return m_toFacesRelation; }

  /**
   * @brief Add a property on the CellBlock
   * @param[in] propertyName the name of the property
   * @return a non-const reference to the property
   */
  template<typename T>
  T & AddProperty( string const & propertyName )
  {
    m_externalPropertyNames.push_back( propertyName );
    return this->registerWrapper< T >( propertyName )->reference();
  }

  template< typename LAMBDA >
  void forExternalProperties( LAMBDA && lambda ) const
  {
    for( auto & externalPropertyName : m_externalPropertyNames )
    {
      const dataRepository::WrapperBase * wrapper = this->getWrapperBase( externalPropertyName );
      lambda( wrapper );
    }
  }

#if USE_ELEM_PATCHES
  NodeMapType & patchNodeList()
  { return m_patchToNodesRelation; }

  NodeMapType const & patchNodeList() const
  { return m_patchToNodesRelation; }

  LvArray::ArrayOfArrays< localIndex, localIndex > & patchNodes()
  { return m_patchNodes; }

  LvArray::ArrayOfArrays< localIndex, localIndex > const & patchNodes() const
  { return m_patchNodes; }

  array1d< localIndex > & patchOffsets()
  { return m_patchOffsets; }

  array1d< localIndex > const & patchOffsets() const
  { return m_patchOffsets; }

#if ELEM_PATCH_VIZ
  array1d< localIndex > & patchIndex()
  { return m_patchIndex; }

  array1d< localIndex > const & patchIndex() const
  { return m_patchIndex; }

  array1d< localIndex > & elemIndex()
  { return m_elemIndex; }

  array1d< localIndex > const & elemIndex() const
  { return m_elemIndex; }
#endif

#endif

protected:


  /// The elements to nodes relation
  NodeMapType  m_toNodesRelation;

  /// The elements to faces relation
  FaceMapType  m_toFacesRelation;

#if USE_ELEM_PATCHES
  /// Offset of each patch (index of first element)
  array1d<localIndex> m_patchOffsets;

  /// List of nodes in each patch
  LvArray::ArrayOfArrays<localIndex, localIndex> m_patchNodes;

  /// Elem-to-node map with patch-local node indices
  NodeMapType m_patchToNodesRelation;

#if ELEM_PATCH_VIZ
  /// Patch index for each element (for visualization/debugging)
  array1d<localIndex> m_patchIndex;

  /// Patch-local index of each element (for visualization/debugging)
  array1d<localIndex> m_elemIndex;
#endif
#endif

private:
  /// Name of the properties register from an external mesh
  string_array m_externalPropertyNames;

};



///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////



}



#endif /* ELEMENTOBJECTT_H_ */
