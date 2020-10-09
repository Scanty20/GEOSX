/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file SolidMechanicsLagrangianFEMKernels.hpp
 */

#pragma once

#include "common/DataTypes.hpp"
#include "common/TimingMacros.hpp"
#include "constitutive/ConstitutiveBase.hpp"
#include "finiteElement/Kinematics.h"
#include "finiteElement/kernelInterface/ImplicitKernelBase.hpp"
#include "rajaInterface/GEOS_RAJA_Interface.hpp"

namespace geosx
{

namespace SolidMechanicsLagrangianFEMKernels
{

inline void velocityUpdate( arrayView2d< real64, nodes::ACCELERATION_USD > const & acceleration,
                            arrayView1d< real64 const > const & mass,
                            arrayView2d< real64, nodes::VELOCITY_USD > const & velocity,
                            R1Tensor const & gravityVector,
                            real64 const dt )
{
  GEOSX_MARK_FUNCTION;

  localIndex const N = acceleration.size( 0 );
  forAll< parallelDevicePolicy<> >( N, [=] GEOSX_DEVICE ( localIndex const i )
  {
    LvArray::tensorOps::scaledAdd< 3 >( velocity[ i ], acceleration[ i ], dt );
    LvArray::tensorOps::fill< 3 >( acceleration[ i ], 0 );
    LvArray::tensorOps::scaledAdd< 3 >( acceleration[ i ], gravityVector, mass[ i ] );
  } );
}

inline void velocityUpdate( arrayView2d< real64, nodes::ACCELERATION_USD > const & acceleration,
                            arrayView1d< real64 const > const & mass,
                            arrayView2d< real64, nodes::VELOCITY_USD > const & velocity,
                            real64 const dt,
                            real64 const massDamping,
                            SortedArrayView< localIndex const > const & indices )
{
  GEOSX_MARK_FUNCTION;

  forAll< parallelDevicePolicy<> >( indices.size(), [=] GEOSX_DEVICE ( localIndex const i )
  {
    localIndex const a = indices[ i ];
    LvArray::tensorOps::scaledAdd< 3 >( acceleration[ a ], velocity[ a ], - massDamping * mass[ a ] );
    LvArray::tensorOps::scale< 3 >( acceleration[ a ], 1.0 / ( mass[ a ] * ( 1 + 0.5 * dt * massDamping ) ) );
    LvArray::tensorOps::scaledAdd< 3 >( velocity[ a ], acceleration[ a ], dt );
  } );
}

inline void displacementUpdate( arrayView2d< real64 const, nodes::VELOCITY_USD > const & velocity,
                                arrayView2d< real64, nodes::INCR_DISPLACEMENT_USD > const & uhat,
                                arrayView2d< real64, nodes::TOTAL_DISPLACEMENT_USD > const & u,
                                real64 const dt )
{
  GEOSX_MARK_FUNCTION;

  localIndex const N = velocity.size( 0 );
  forAll< parallelDevicePolicy<> >( N, [=] GEOSX_DEVICE ( localIndex const i )
  {
    LvArray::tensorOps::scaledCopy< 3 >( uhat[ i ], velocity[ i ], dt );
    LvArray::tensorOps::add< 3 >( u[ i ], uhat[ i ] );
  } );
}


/**
 * @struct Structure to wrap templated function that implements the explicit time integration kernel.
 */
struct ExplicitKernel
{

  static inline real64
  CalculateSingleNodalForce( localIndex const k,
                             localIndex const targetNode,
                             localIndex const numQuadraturePoints,
                             arrayView4d< real64 const > const & dNdX,
                             arrayView2d< real64 const > const & detJ,
                             arrayView3d< real64 const, solid::STRESS_USD > const & stress,
                             R1Tensor & force )
  {
    GEOSX_MARK_FUNCTION;
    localIndex const & a = targetNode;

    //Compute Quadrature
    for( localIndex q = 0; q < numQuadraturePoints; ++q )
    {
      force[ 0 ] -= ( stress( k, q, 0 ) * dNdX( k, q, a, 0 ) +
                      stress( k, q, 5 ) * dNdX( k, q, a, 1 ) +
                      stress( k, q, 4 ) * dNdX( k, q, a, 2 ) ) * detJ( k, q );
      force[ 1 ] -= ( stress( k, q, 5 ) * dNdX( k, q, a, 0 ) +
                      stress( k, q, 1 ) * dNdX( k, q, a, 1 ) +
                      stress( k, q, 3 ) * dNdX( k, q, a, 2 ) ) * detJ( k, q );
      force[ 2 ] -= ( stress( k, q, 4 ) * dNdX( k, q, a, 0 ) +
                      stress( k, q, 3 ) * dNdX( k, q, a, 1 ) +
                      stress( k, q, 2 ) * dNdX( k, q, a, 2 ) ) * detJ( k, q );

    }//quadrature loop

    return 0;
  }

};


} // namespace SolidMechanicsLagrangianFEMKernels

} // namespace geosx
