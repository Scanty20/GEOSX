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
 * @file InterpolationBuilderMsRsb.hpp
 */

#ifndef GEOSX_LINEARALGEBRA_MSRSBSTRATEGY_HPP_
#define GEOSX_LINEARALGEBRA_MSRSBSTRATEGY_HPP_

#include "linearAlgebra/solvers/multiscale/MultiscaleStrategy.hpp"

namespace geosx
{
namespace multiscale
{

template< typename LAI >
class MsRSBStrategy : public MultiscaleStrategy< LAI >
{
public:

  /// Alias for base type
  using Base = MultiscaleStrategy< LAI >;

  /// Alias for vector type
  using Vector = typename Base::Vector;

  /// Alias for matrix type
  using Matrix = typename Base::Matrix;

  /// Alias for operator type
  using Operator = typename Base::Operator;

  explicit MsRSBStrategy( string name,
                          LinearSolverParameters::Multiscale params );

  virtual Operator const & getProlongation() const override
  {
    return m_prolongation;
  }

  virtual Operator const & getRestriction() const override
  {
    return *m_restriction;
  }

  virtual Operator const & getOperator() const override
  {
    return m_matrix;
  }

  virtual void initializeFine( MeshLevel & mesh,
                               DofManager const & dofManager,
                               string const & fieldName ) override;

  virtual void initializeCoarse( MultiscaleStrategy< LAI > const & fine_level ) override;

  virtual void compute( Operator const & fine_operator ) override;

private:

  using Base::m_params;

  struct Mesh
  {
    SparsityPattern< globalIndex > c2n;

    Matrix cellToNode;     ///< map from cells (rows) to nodes (cols)
    Matrix nodeToCell;     ///< map from nodes (rows) to cells (cols)
    Matrix nodeToBoundary; ///< map from nodes (rows) to boundaries (cols)
    Matrix cellToBoundary; ///< map from cells (rows) to boundaries (cols)

    array1d< localIndex > nodeToMesh; ///< map from
  };

  /// Mesh description at current level
  Mesh m_mesh;

  /// Prolongation matrix P
  Matrix m_prolongation;

  /// Restriction (kept as abstract operator to allow for memory efficiency, e.g. when R = P^T)
  std::unique_ptr< Operator > m_restriction;

  /// Level operator matrix
  Matrix m_matrix;
};

} // namespace multiscale
} // namespace geosx

#endif //GEOSX_LINEARALGEBRA_MSRSBSTRATEGY_HPP_
