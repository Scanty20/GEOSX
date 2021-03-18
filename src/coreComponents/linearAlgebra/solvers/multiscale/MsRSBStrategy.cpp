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
 * @file InterpolationBuilderMsRsb.cpp
 */

#include "MsRSBStrategy.hpp"

#include "linearAlgebra/interfaces/InterfaceTypes.hpp"
#include "linearAlgebra/utilities/LAIHelperFunctions.hpp"
#include "linearAlgebra/utilities/TransposeOperator.hpp"

namespace geosx
{
namespace multiscale
{

template< typename LAI >
MsRSBStrategy< LAI >::MsRSBStrategy( string name,
                                     LinearSolverParameters::Multiscale params )
  : MultiscaleStrategy< LAI >( std::move( name ), std::move( params ) )
{}

template< typename LAI >
void MsRSBStrategy< LAI >::initializeFine( MeshLevel & mesh,
                                           DofManager const & dofManager,
                                           string const & fieldName )
{
  GEOSX_UNUSED_VAR( mesh, dofManager, fieldName )
}

template< typename LAI >
void MsRSBStrategy< LAI >::initializeCoarse( MultiscaleStrategy< LAI > const & fine_level )
{
  auto const & fine = dynamicCast< MsRSBStrategy< LAI > const & >( fine_level );
  GEOSX_UNUSED_VAR( fine )
}

template< typename VECTOR, typename MATRIX >
MATRIX makeMsrsbIterationMatrix( MATRIX const & fine_mat,
                                 localIndex const numDof,
                                 real64 const omega )
{
  // 1. Apply SC approximation
  MATRIX iterMat;
  LAIHelperFunctions::separateComponentFilter( fine_mat, iterMat, numDof );

  // 2. Compute -w * Dinv * Asc;
  VECTOR diag;
  iterMat.extractDiagonal( diag );
  diag.reciprocal();
  diag.scale( -omega );
  iterMat.leftScale( diag );

  // 3. Compute I - w * Dinv * Asc
  diag.set( 1.0 );
  iterMat.addDiagonal( diag );
  return iterMat;
}

template< typename LAI >
void MsRSBStrategy< LAI >::compute( Operator const & fine_operator )
{
  Matrix const & fine_mat = dynamicCast< Matrix const & >( fine_operator );
  localIndex const numDof = 3; // fine_mat.numLocalRows() / m_mesh.???
  Matrix const iterMat = makeMsrsbIterationMatrix< Vector >( fine_mat, numDof, m_params.msrsb.relaxation );

  for( integer iter = 0; iter < m_params.msrsb.maxSmoothingIter; ++iter )
  {
    /// TODO smoothing
  }

  /// Compute coarse operator
  fine_mat.multiplyPtAP( m_prolongation, m_matrix );

  /// Make restriction
  m_restriction = std::make_unique< TransposeOperator< LAI > >( m_prolongation );
}

// -----------------------
// Explicit Instantiations
// -----------------------
#ifdef GEOSX_USE_TRILINOS
template class MsRSBStrategy< TrilinosInterface >;
#endif

#ifdef GEOSX_USE_HYPRE
template class MsRSBStrategy< HypreInterface >;
#endif

#ifdef GEOSX_USE_PETSC
template class MsRSBStrategy< PetscInterface >;
#endif

} // namespace multiscale
} // namespace geosx
