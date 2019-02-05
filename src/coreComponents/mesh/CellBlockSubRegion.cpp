/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2019, Lawrence Livermore National Security, LLC.
 *
 * Produced at the Lawrence Livermore National Laboratory
 *
 * LLNL-CODE-746361
 *
 * All rights reserved. See COPYRIGHT for details.
 *
 * This file is part of the GEOSX Simulation Framework.
 *
 * GEOSX is a free software; you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License (as published by the
 * Free Software Foundation) version 2.1 dated February 1999.
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */


#include "CellBlockSubRegion.hpp"
#include "constitutive/ConstitutiveManager.hpp"

namespace geosx
{
using namespace dataRepository;
using namespace constitutive;

CellBlockSubRegion::CellBlockSubRegion( string const & name, ManagedGroup * const parent ):
  CellBlock( name, parent )
{
  RegisterViewWrapper( viewKeyStruct::constitutiveGroupingString, &m_constitutiveGrouping, 0)->
    setSizedFromParent(0);

  RegisterViewWrapper( viewKeyStruct::constitutiveMapString,
                       &m_constitutiveMapView, 0);

  RegisterViewWrapper( viewKeyStruct::dNdXString, &m_dNdX, 0);


  RegisterViewWrapper( viewKeyStruct::constitutivePointVolumeFraction,
                       &m_constitutivePointVolumeFraction, 0);

  RegisterViewWrapper( viewKeyStruct::dNdXString, &m_dNdX, 0)->setSizedFromParent(1);

}

CellBlockSubRegion::~CellBlockSubRegion()
{
  // TODO Auto-generated destructor stub
}

void CellBlockSubRegion::CopyFromCellBlock( CellBlock const * source )
{
  this->SetElementType(source->GetElementType());
  this->numNodesPerElement() = source->numNodesPerElement();
  this->numFacesPerElement() = source->numFacesPerElement();
  this->resize(source->size());
  this->nodeList() = source->nodeList();
  this->m_localToGlobalMap = source->m_localToGlobalMap;
  this->ConstructGlobalToLocalMap();
}

void CellBlockSubRegion::ConstructSubRegionFromFaceSet( FaceManager const * const faceManager,
                                                        string const & setName )
{
  set<localIndex> const & targetSet = faceManager->sets()->getReference<set<localIndex> >(setName);
  m_toFacesRelation.resize(0,2);
  this->resize( targetSet.size() );


}

void CellBlockSubRegion::MaterialPassThru( string const & matName,
                                           string const & setName,
                                           set<localIndex> & materialSet,
                                           ManagedGroup * material )
{}





void CellBlockSubRegion::ViewPackingExclusionList( set<localIndex> & exclusionList ) const
{
  ObjectManagerBase::ViewPackingExclusionList(exclusionList);
  exclusionList.insert(this->getWrapperIndex(viewKeyStruct::nodeListString));
//  exclusionList.insert(this->getWrapperIndex(this->viewKeys.edgeListString));
  exclusionList.insert(this->getWrapperIndex(viewKeyStruct::faceListString));
}


localIndex CellBlockSubRegion::PackUpDownMapsSize( arrayView1d<localIndex const> const & packList ) const
{
  buffer_unit_type * junk = nullptr;
  return PackUpDownMapsPrivate<false>( junk, packList );
}


localIndex CellBlockSubRegion::PackUpDownMaps( buffer_unit_type * & buffer,
                                               arrayView1d<localIndex const> const & packList ) const
{
  return PackUpDownMapsPrivate<true>( buffer, packList );
}

template< bool DOPACK >
localIndex CellBlockSubRegion::PackUpDownMapsPrivate( buffer_unit_type * & buffer,
                                                      arrayView1d<localIndex const> const & packList ) const
{
  localIndex packedSize = 0;

  packedSize += bufferOps::Pack<DOPACK>( buffer,
                                         nodeList().Base(),
                                         m_unmappedGlobalIndicesInNodelist,
                                         packList,
                                         this->m_localToGlobalMap,
                                         nodeList().RelatedObjectLocalToGlobal() );

  packedSize += bufferOps::Pack<DOPACK>( buffer,
                                         faceList().Base(),
                                         m_unmappedGlobalIndicesInFacelist,
                                         packList,
                                         this->m_localToGlobalMap,
                                         faceList().RelatedObjectLocalToGlobal() );

  return packedSize;
}


localIndex CellBlockSubRegion::UnpackUpDownMaps( buffer_unit_type const * & buffer,
                                                 localIndex_array & packList )
{
  localIndex unPackedSize = 0;



  unPackedSize += bufferOps::Unpack( buffer,
                                     nodeList().Base(),
                                     packList,
                                     m_unmappedGlobalIndicesInNodelist,
                                     this->m_globalToLocalMap,
                                     nodeList().RelatedObjectGlobalToLocal() );

  unPackedSize += bufferOps::Unpack( buffer,
                                     faceList().Base(),
                                     packList,
                                     m_unmappedGlobalIndicesInFacelist,
                                     this->m_globalToLocalMap,
                                     faceList().RelatedObjectGlobalToLocal() );

  return unPackedSize;
}

void CellBlockSubRegion::FixUpDownMaps( bool const clearIfUnmapped )
{
  ObjectManagerBase::FixUpDownMaps( nodeList(),
                                    m_unmappedGlobalIndicesInNodelist,
                                    clearIfUnmapped );

  ObjectManagerBase::FixUpDownMaps( faceList(),
                                    m_unmappedGlobalIndicesInFacelist,
                                    clearIfUnmapped );
}


} /* namespace geosx */
