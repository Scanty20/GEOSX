/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2018, Lawrence Livermore National Security, LLC.
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

/*
 * @file FiniteVolumeManager.cpp
 *
 */

#include "FiniteVolumeManager.hpp"

#include "finiteVolume/FluxApproximationBase.hpp"
#include "managers/DomainPartition.hpp"
#include "mesh/MeshForLoopInterface.hpp"
#include "meshUtilities/ComputationalGeometry.hpp"

namespace geosx
{

using namespace dataRepository;


FiniteVolumeManager::FiniteVolumeManager(string const &name, ManagedGroup *const parent)
  : ManagedGroup(name, parent)
{
  setSchemaFlags(SchemaFlags::UNIQUE_NODE);
}

FiniteVolumeManager::~FiniteVolumeManager()
{

}

ManagedGroup * FiniteVolumeManager::CreateChild(string const &childKey, string const &childName)
{
  std::unique_ptr<FluxApproximationBase> approx = FluxApproximationBase::CatalogInterface::Factory(childKey, childName, this);
  return this->RegisterGroup<FluxApproximationBase>(childName, std::move(approx));
}


void FiniteVolumeManager::ExpandObjectCatalogs()
{
  // During schema generation, register one of each type derived from FluxApproximationBase here
  for (auto& catalogIter: FluxApproximationBase::GetCatalog())
  {
    CreateChild( catalogIter.first, catalogIter.first );
  }
}


FluxApproximationBase const * FiniteVolumeManager::getFluxApproximation(std::string const &name) const
{
  return this->GetGroup<FluxApproximationBase>(name);
}

void FiniteVolumeManager::IntermediateInitializationPreSubGroups(ManagedGroup * const rootGroup)
{
  DomainPartition * domain = rootGroup->GetGroup<DomainPartition>(keys::domain);
  precomputeFiniteVolumeData(domain);
}

void FiniteVolumeManager::precomputeFiniteVolumeData(DomainPartition * const domain)
{
  MeshLevel * const mesh = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);
  NodeManager * const nodeManager = mesh->getNodeManager();
  ElementRegionManager * const elemManager = mesh->getElemManager();
  FaceManager * const faceManager = mesh->getFaceManager();

  r1_array const & X = nodeManager->referencePosition();

  ElementRegionManager::ElementViewAccessor<arrayView1d<R1Tensor>> elemCenter = 
    elemManager->ConstructViewAccessor<array1d<R1Tensor>, arrayView1d<R1Tensor>>(
                                        CellBlock::viewKeyStruct::elementCenterString);

  ElementRegionManager::ElementViewAccessor<arrayView1d<real64>> elemVolume = 
    elemManager->ConstructViewAccessor<array1d<real64>, arrayView1d<real64>>(
                                        CellBlock::viewKeyStruct::elementVolumeString);

  ElementRegionManager::ElementViewAccessor<arrayView2d<localIndex>> const elemsToNodes =
    elemManager->ConstructViewAccessor<FixedOneToManyRelation, arrayView2d<localIndex>>(
                                    CellBlockSubRegion::viewKeyStruct::nodeListString);


  // Loop over all the elements and calculate element centers, and element volumes
  forAllElemsInMesh( mesh, [&]( localIndex const er,
                                localIndex const esr,
                                localIndex const k )->void
  {
    localIndex const nodeListSize = elemsToNodes[er][esr].size(1);
    R1Tensor Xlocal[ElementRegionManager::maxNumNodesPerElem];

    R1Tensor & center = elemCenter[er][esr][k];
    center = 0.0;

    // TODO different center options
    for (localIndex a = 0; a < nodeListSize; ++a)
    {
      Xlocal[a] = X[elemsToNodes[er][esr][k][a]];
      center += Xlocal[a];
    }
    center /= nodeListSize;

    // TODO proper volumes for all shapes
    if( nodeListSize == 8 )
    {
        elemVolume[er][esr][k] = computationalGeometry::HexVolume(Xlocal);
    }
    else if( nodeListSize == 4)
    {
        elemVolume[er][esr][k] = computationalGeometry::TetVolume(Xlocal);
    }
    else if( nodeListSize == 6)
    {
        elemVolume[er][esr][k] = computationalGeometry::WedgeVolume(Xlocal);
    }
    else if ( nodeListSize == 5)
    {
        elemVolume[er][esr][k] = computationalGeometry::PyramidVolume(Xlocal);
    }
    else
    {
        GEOS_ERROR("GEOX does not support cells with " << nodeListSize << " nodes");
    }
  });

  r1_array & faceCenter = faceManager->getReference<r1_array>(FaceManager::viewKeyStruct::faceCenterString);
  array1d<array1d<localIndex>> const & faceToNodes = faceManager->nodeList();

  R1Tensor normal;
  for (localIndex kf = 0; kf < faceManager->size(); ++kf)
  {
    computationalGeometry::Centroid_3DPolygon(faceToNodes[kf], X, faceCenter[kf], normal);
  }
}

} // namespace geosx
