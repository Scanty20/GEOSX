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
 * @file PreconditionerMultiscale.cpp
 */

#include "PreconditionerMultiscale.hpp"

#include "linearAlgebra/interfaces/InterfaceTypes.hpp"
#include "linearAlgebra/solvers/multiscale/MultiscaleStrategy.hpp"

namespace geosx
{

template< typename LAI >
PreconditionerMultiscale< LAI >::PreconditionerMultiscale( LinearSolverParameters params,
                                                           MeshLevel & mesh )
  : Base(),
  m_params( std::move( params.multiscale ) ),
  m_mesh( mesh ),
  m_initialized{ false }
{}

template< typename LAI >
void PreconditionerMultiscale< LAI >::setup( Matrix const & mat )
{
  Base::setup( mat );

  if( !m_initialized )
  {
    GEOSX_ERROR_IF( mat.dofManager() == nullptr, "PreconditionerMultiscale: DofManager is required" );
    createLevels( mat, *mat.dofManager() );
  }

  for( integer level_index = 1; level_index < m_params.maxLevels; ++level_index )
  {
    m_levels[level_index].strategy->compute( *m_levels[level_index-1].operator_ptr );
  }
}

template< typename LAI >
void PreconditionerMultiscale< LAI >::createLevels( Matrix const & mat,
                                                    DofManager const & dofManager )
{
  m_levels.clear();

  // create fine level
  m_levels.emplace_back();
  Level & fine = m_levels[0];
  fine.strategy = multiscale::MultiscaleStrategy< LAI >::createInstance( "level_0", m_params );
  fine.strategy->initializeFine( m_mesh, dofManager, m_params.fieldName );
  fine.operator_ptr = &mat;

  // create coarse levels
  for( integer level_index = 1; level_index < m_params.maxLevels; ++level_index )
  {
    m_levels.emplace_back();
    Level & coarse = m_levels[level_index];
    coarse.strategy = multiscale::MultiscaleStrategy< LAI >::createInstance( "level_" + std::to_string( level_index ), m_params );
    coarse.strategy->initializeCoarse( *m_levels[level_index - 1].strategy );
    coarse.operator_ptr = &coarse.strategy->getOperator();
    if( coarse.operator_ptr->numGlobalRows() <= m_params.minGlobalDof )
    {
      break;
    }
  }

  // create smoothers
  for( size_t level_index = 0; level_index < m_levels.size()-1; ++level_index )
  {
    Level & level = m_levels[level_index];
    // TODO: smoother options from input
    LinearSolverParameters smoother_params;
    smoother_params.preconditionerType = LinearSolverParameters::PrecondType::jacobi;
    if( m_params.preOrPostSmoothing == LinearSolverParameters::AMG::PreOrPost::pre
        || m_params.preOrPostSmoothing == LinearSolverParameters::AMG::PreOrPost::both )
    {
      level.presmoother = LAI::createPreconditioner( smoother_params );
    }
    if( m_params.preOrPostSmoothing == LinearSolverParameters::AMG::PreOrPost::post
        || m_params.preOrPostSmoothing == LinearSolverParameters::AMG::PreOrPost::both )
    {
      level.postsmoother = LAI::createPreconditioner( smoother_params );
    }
  }

  // create vectors
  for( size_t level_index = 0; level_index < m_levels.size(); ++level_index )
  {
    Level & level = m_levels[level_index];
    level.rhs.createWithLocalSize( level.operator_ptr->numLocalRows(), level.operator_ptr->getComm() );
    level.sol.createWithLocalSize( level.operator_ptr->numLocalRows(), level.operator_ptr->getComm() );
    level.tmp.createWithLocalSize( level.operator_ptr->numLocalRows(), level.operator_ptr->getComm() );
  }
}

template< typename LAI >
void PreconditionerMultiscale< LAI >::apply( Vector const & src,
                                             Vector & dst ) const
{
  // TODO: remove hardcoded V-cycle, abstract into a separate component

  m_levels[0].rhs.copy( src );

  // down phase
  for( size_t level_index = 0; level_index < m_levels.size() - 1; ++level_index )
  {
    Level const & fine = m_levels[level_index];
    Level const & coarse = m_levels[level_index + 1];
    if( fine.presmoother )
    {
      fine.presmoother->apply( fine.rhs, fine.sol );
      fine.operator_ptr->residual( fine.sol, fine.rhs, fine.rhs );
    }
    coarse.strategy->getRestriction().apply( fine.rhs, coarse.rhs );
  }

  // coarse level solve
  m_coarse_solver->apply( m_levels.back().rhs, m_levels.back().sol );

  // up phase
  for( size_t level_index = m_levels.size() - 1; level_index > 0; --level_index )
  {
    Level const & fine = m_levels[level_index - 1];
    Level const & coarse = m_levels[level_index];
    coarse.strategy->getProlongation().apply( coarse.sol, fine.tmp );
    fine.sol.axpy( 1.0, fine.tmp );
    fine.operator_ptr->residual( fine.tmp, fine.rhs, fine.rhs );
    if( fine.postsmoother )
    {
      fine.postsmoother->apply( fine.rhs, fine.tmp );
      fine.sol.axpy( 1.0, fine.tmp );
      fine.operator_ptr->residual( fine.tmp, fine.rhs, fine.rhs );
    }
  }

  dst.copy( m_levels[0].sol );
}

template< typename LAI >
PreconditionerMultiscale< LAI >::~PreconditionerMultiscale() = default;

// -----------------------
// Explicit Instantiations
// -----------------------
#ifdef GEOSX_USE_TRILINOS
template class PreconditionerMultiscale< TrilinosInterface >;
#endif

#ifdef GEOSX_USE_HYPRE
template class PreconditionerMultiscale< HypreInterface >;
#endif

#ifdef GEOSX_USE_PETSC
template class PreconditionerMultiscale< PetscInterface >;
#endif

} // namespace geosx
