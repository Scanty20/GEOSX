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
 * @file MultiphasePoroelasticKernel.hpp
 */

#ifndef GEOSX_PHYSICSSOLVERS_MULTIPHYSICS_MULTIPHASEPOROELASTICKERNEL_HPP_
#define GEOSX_PHYSICSSOLVERS_MULTIPHYSICS_MULTIPHASEPOROELASTICKERNEL_HPP_
#include "finiteElement/kernelInterface/ImplicitKernelBase.hpp"
#include "physicsSolvers/fluidFlow/CompositionalMultiphaseBase.hpp"

namespace geosx
{

namespace PoroelasticKernels
{

/**
 * @brief Implements kernels for solving quasi-static multiphase poromechanics.
 * @copydoc geosx::finiteElement::ImplicitKernelBase
 * @tparam NUM_NODES_PER_ELEM The number of nodes per element for the
 *                            @p SUBREGION_TYPE.
 * @tparam UNUSED An unused parameter since we are assuming that the test and
 *                trial space have the same number of support points.
 *
 * ### MultiphasePoroelastic Description
 * Implements the KernelBase interface functions required for solving the
 * quasi-static multiphase poromechanics problem using one of the
 * "finite element kernel application" functions such as
 * geosx::finiteElement::RegionBasedKernelApplication.
 *
 */
template< typename SUBREGION_TYPE,
          typename CONSTITUTIVE_TYPE,
          typename FE_TYPE >
class Multiphase :
  public finiteElement::ImplicitKernelBase< SUBREGION_TYPE,
                                            CONSTITUTIVE_TYPE,
                                            FE_TYPE,
                                            3,
                                            3 >
{
public:
  /// Alias for the base class;
  using Base = finiteElement::ImplicitKernelBase< SUBREGION_TYPE,
                                                  CONSTITUTIVE_TYPE,
                                                  FE_TYPE,
                                                  3,
                                                  3 >;

  /// Number of nodes per element...which is equal to the
  /// numTestSupportPointPerElem and numTrialSupportPointPerElem by definition.
  static constexpr int numNodesPerElem = Base::numTestSupportPointsPerElem;
  static constexpr int numMaxComponentsMultiphasePoroelastic = 3;
  using Base::numDofPerTestSupportPoint;
  using Base::numDofPerTrialSupportPoint;
  using Base::m_dofNumber;
  using Base::m_dofRankOffset;
  using Base::m_matrix;
  using Base::m_rhs;
  using Base::m_elemsToNodes;
  using Base::m_constitutiveUpdate;
  using Base::m_finiteElementSpace;

  /**
   * @brief Constructor
   * @copydoc geosx::finiteElement::ImplicitKernelBase::ImplicitKernelBase
   * @param inputGravityVector The gravity vector.
   */
  Multiphase( NodeManager const & nodeManager,
              EdgeManager const & edgeManager,
              FaceManager const & faceManager,
              localIndex const targetRegionIndex,
              SUBREGION_TYPE const & elementSubRegion,
              FE_TYPE const & finiteElementSpace,
              CONSTITUTIVE_TYPE & inputConstitutiveType,
              arrayView1d< globalIndex const > const & inputDispDofNumber,
              string const & inputFlowDofKey,
              globalIndex const rankOffset,
              real64 const (&inputGravityVector)[3],
              localIndex const numComponents,
              localIndex const numPhases,
              arrayView1d< string const > const & fluidModelNames,
              CRSMatrixView< real64, globalIndex const > const & inputMatrix,
              arrayView1d< real64 > const & inputRhs ):
    Base( nodeManager,
          edgeManager,
          faceManager,
          targetRegionIndex,
          elementSubRegion,
          finiteElementSpace,
          inputConstitutiveType,
          inputDispDofNumber,
          rankOffset,
          inputMatrix,
          inputRhs ),
    m_X( nodeManager.referencePosition()),
    m_disp( nodeManager.totalDisplacement()),
    m_uhat( nodeManager.incrementalDisplacement()),
    m_gravityVector{ inputGravityVector[0], inputGravityVector[1], inputGravityVector[2] },
    m_gravityAcceleration( LvArray::tensorOps::l2Norm< 3 >( inputGravityVector ) ),
    m_solidDensity( inputConstitutiveType.getDensity() ),
    m_fluidPhaseDensity( elementSubRegion.template getConstitutiveModel< constitutive::MultiFluidBase >( fluidModelNames[targetRegionIndex] ).phaseDensity() ),
    m_fluidPhaseDensityOld( elementSubRegion.template getReference< array2d< real64 > >( CompositionalMultiphaseBase::viewKeyStruct::phaseDensityOldString() ) ),
    m_dFluidPhaseDensity_dPressure( elementSubRegion.template getConstitutiveModel< constitutive::MultiFluidBase >( fluidModelNames[targetRegionIndex] ).dPhaseDensity_dPressure() ),
    m_dFluidPhaseDensity_dGlobalCompFraction( elementSubRegion.template getConstitutiveModel< constitutive::MultiFluidBase >( fluidModelNames[targetRegionIndex] ).dPhaseDensity_dGlobalCompFraction() ),
    m_fluidPhaseCompFrac( elementSubRegion.template getConstitutiveModel< constitutive::MultiFluidBase >( fluidModelNames[targetRegionIndex] ).phaseCompFraction() ),
    m_fluidPhaseCompFracOld( elementSubRegion.template getReference< array3d< real64 > >( CompositionalMultiphaseBase::viewKeyStruct::phaseComponentFractionOldString() ) ),
    m_dFluidPhaseCompFrac_dPressure( elementSubRegion.template getConstitutiveModel< constitutive::MultiFluidBase >( fluidModelNames[targetRegionIndex] ).dPhaseCompFraction_dPressure() ),
    m_fluidPhaseMassDensity( elementSubRegion.template getConstitutiveModel< constitutive::MultiFluidBase >( fluidModelNames[targetRegionIndex] ).phaseMassDensity() ),
    m_fluidPhaseSaturation( elementSubRegion.template getReference< array2d< real64 > >( CompositionalMultiphaseBase::viewKeyStruct::phaseVolumeFractionString() )),
    m_fluidPhaseSaturationOld( elementSubRegion.template getReference< array2d< real64 > >( CompositionalMultiphaseBase::viewKeyStruct::phaseVolumeFractionOldString() )),
    m_dFluidPhaseSaturation_dPressure( elementSubRegion.template getReference< array2d< real64 > >( CompositionalMultiphaseBase::viewKeyStruct::dPhaseVolumeFraction_dPressureString() )),
    m_dFluidPhaseSaturation_dGlobalCompDensity( elementSubRegion.template getReference< array3d< real64 > >(
                                                  CompositionalMultiphaseBase::viewKeyStruct:: dPhaseVolumeFraction_dGlobalCompDensityString() )),
    m_dGlobalCompFraction_dGlobalCompDensity( elementSubRegion.template getReference< array3d< real64 > >( CompositionalMultiphaseBase::viewKeyStruct::dGlobalCompFraction_dGlobalCompDensityString() )),
    m_dFluidPhaseCompFraction_dGlobalCompFraction( elementSubRegion.template getConstitutiveModel< constitutive::MultiFluidBase >(
                                                     fluidModelNames[targetRegionIndex] ).dPhaseCompFraction_dGlobalCompFraction() ),
    m_flowDofNumber( elementSubRegion.template getReference< array1d< globalIndex > >( inputFlowDofKey )),
    m_fluidPressure( elementSubRegion.template getReference< array1d< real64 > >( FlowSolverBase::viewKeyStruct::pressureString() ) ),
    m_deltaFluidPressure( elementSubRegion.template getReference< array1d< real64 > >( FlowSolverBase::viewKeyStruct::deltaPressureString() ) ),
    m_poroRef( elementSubRegion.template getReference< array1d< real64 > >( FlowSolverBase::viewKeyStruct::referencePorosityString() ) ),
    m_numComponents( numComponents ),
    m_numPhases( numPhases ),
    m_biotCoefficient( m_constitutiveUpdate.getBiotCoefficient() )
  {
    if( m_numComponents > numMaxComponentsMultiphasePoroelastic )
    {
      GEOSX_ERROR( "MultiphasePoroelastic solver allows at most three components at the moment" );
    }
  }

  //*****************************************************************************
  /**
   * @class StackVariables
   * @copydoc geosx::finiteElement::ImplicitKernelBase::StackVariables
   *
   * Adds a stack array for the displacement, incremental displacement, and the
   * constitutive stiffness.
   */
  struct StackVariables : public Base::StackVariables
  {
public:

    static constexpr int numDispDofPerElem =  Base::StackVariables::numRows;

    /// Constructor.
    GEOSX_HOST_DEVICE
    StackVariables():
      Base::StackVariables(),
            xLocal(),
            u_local(),
            uhat_local(),
            localFlowResidual{ 0.0 },
      localDispFlowJacobian{ {0.0} },
      localFlowDispJacobian{ {0.0} },
      localFlowFlowJacobian{ {0.0} },
      localFlowDofIndex{ 0 }
    {}

#if !defined(CALC_FEM_SHAPE_IN_KERNEL)
    /// Dummy
    int xLocal;
#else
    /// C-array stack storage for element local the nodal positions.
    real64 xLocal[numNodesPerElem][3];
#endif

    /// Stack storage for the element local nodal displacement
    real64 u_local[numNodesPerElem][numDofPerTrialSupportPoint];

    /// Stack storage for the element local nodal incremental displacement
    real64 uhat_local[numNodesPerElem][numDofPerTrialSupportPoint];

    real64 localFlowResidual[numMaxComponentsMultiphasePoroelastic];
    real64 localDispFlowJacobian[numDispDofPerElem][numMaxComponentsMultiphasePoroelastic + 1];
    real64 localFlowDispJacobian[numMaxComponentsMultiphasePoroelastic][numDispDofPerElem];
    real64 localFlowFlowJacobian[numMaxComponentsMultiphasePoroelastic][numMaxComponentsMultiphasePoroelastic + 1];

    /// C-array storage for the element local row degrees of freedom.
    globalIndex localFlowDofIndex[numMaxComponentsMultiphasePoroelastic+1];

  };
  //*****************************************************************************

  /**
   * @brief Copy global values from primary field to a local stack array.
   * @copydoc ::geosx::finiteElement::ImplicitKernelBase::setup
   *
   * For the Multiphase implementation, global values from the displacement,
   * incremental displacement, and degree of freedom numbers are placed into
   * element local stack storage.
   */
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  void setup( localIndex const k,
              StackVariables & stack ) const
  {
    for( localIndex a=0; a<numNodesPerElem; ++a )
    {
      localIndex const localNodeIndex = m_elemsToNodes( k, a );

      for( int i=0; i<3; ++i )
      {
#if defined(CALC_FEM_SHAPE_IN_KERNEL)
        stack.xLocal[a][i] = m_X[localNodeIndex][i];
#endif
        stack.u_local[a][i] = m_disp[localNodeIndex][i];
        stack.uhat_local[a][i] = m_uhat[localNodeIndex][i];
        stack.localRowDofIndex[a*3+i] = m_dofNumber[localNodeIndex]+i;
        stack.localColDofIndex[a*3+i] = m_dofNumber[localNodeIndex]+i;
      }
    }

    for( int flowDofIndex=0; flowDofIndex<numMaxComponentsMultiphasePoroelastic+1; ++flowDofIndex )
    {
      stack.localFlowDofIndex[flowDofIndex] = m_flowDofNumber[k] + flowDofIndex;
    }

  }

  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  void quadraturePointKernel( localIndex const k,
                              localIndex const q,
                              StackVariables & stack ) const
  {
    localIndex const NC = m_numComponents;
    localIndex const NP = m_numPhases;

    // Get displacement: (i) basis functions (N), (ii) basis function
    // derivatives (dNdX), and (iii) determinant of the Jacobian transformation
    // matrix times the quadrature weight (detJxW)
    real64 N[numNodesPerElem];
    real64 dNdX[numNodesPerElem][3];
    FE_TYPE::calcN( q, N );
    real64 const detJxW = m_finiteElementSpace.template getGradN< FE_TYPE >( k, q, stack.xLocal, dNdX );

    // Evaluate total stress tensor
    real64 strainIncrement[6] = {0};
    real64 totalStress[6];

    // --- Update effective stress tensor (stored in totalStress)
    typename CONSTITUTIVE_TYPE::KernelWrapper::DiscretizationOps stiffness;
    FE_TYPE::symmetricGradient( dNdX, stack.uhat_local, strainIncrement );
    m_constitutiveUpdate.smallStrainUpdate( k, q, strainIncrement, totalStress, stiffness );

    // --- Subtract pressure term
    real64 const biotTimesPressure = m_biotCoefficient * ( m_fluidPressure[k] + m_deltaFluidPressure[k] );
    totalStress[0] -= biotTimesPressure;
    totalStress[1] -= biotTimesPressure;
    totalStress[2] -= biotTimesPressure;

    // Evaluate fluid phase mass content
    // --- TODO: temporary solution -----------------------------------------------------------------------------------//
    //           see SinglePhasePoromechanicsKernel                                                                    //
    real64 const biotSkeletonModulusInverse = 0.0; //TODO: 1/N = 0 correct only for biotCoefficient = 1                //
    real64 const volumetricStrainNew = FE_TYPE::symmetricGradientTrace( dNdX, stack.u_local );                         //
    real64 const volumetricStrainOld = volumetricStrainNew - FE_TYPE::symmetricGradientTrace( dNdX, stack.uhat_local );//
    real64 const porosityOld = m_poroRef( k ) + m_biotCoefficient * volumetricStrainOld;// +  DeltaPoro                //
    real64 const dPorosity_dPressure = biotSkeletonModulusInverse;                                                     //
    real64 const dPorosity_dVolStrainIncrement =  m_biotCoefficient;                                                   //
    GEOSX_ERROR_IF_GT_MSG( fabs( m_biotCoefficient - 1.0 ),                                                            //
                           1e-10,                                                                                      //
                           "Correct only for Biot's coefficient equal to 1" );                                         //
    // --------------------------------------------------------------------------------------------------------------- //
    real64 const porosityNew = porosityOld
                               + m_biotCoefficient * (strainIncrement[0] + strainIncrement[1] + strainIncrement[2] )
                               + biotSkeletonModulusInverse * m_deltaFluidPressure[k];

    // Evaluate body force vector
    real64 bodyForce[3] = { m_gravityVector[0],
                            m_gravityVector[1],
                            m_gravityVector[2]};
    if( m_gravityAcceleration > 0.0 )
    {
      // Compute mixture density
      real64 mixtureDensity = m_fluidPhaseSaturation( k, 0 ) * m_fluidPhaseMassDensity( k, q, 0 );
      for( localIndex i = 1; i < NP; ++i )
      {
        mixtureDensity += m_fluidPhaseSaturation( k, i ) * m_fluidPhaseMassDensity( k, q, i );
      }
      mixtureDensity *= porosityNew;
      mixtureDensity += ( 1.0 - porosityNew ) * m_solidDensity( k, q );
      mixtureDensity *= detJxW;
      bodyForce[0] *= mixtureDensity;
      bodyForce[1] *= mixtureDensity;
      bodyForce[2] *= mixtureDensity;
    }

    // Assemble local jacobian and residual

    // --- Momentum balance
    for( localIndex i=0; i<6; ++i )
    {
      totalStress[i] *= -detJxW;
    }

    FE_TYPE::plusGradNajAijPlusNaFi( dNdX,
                                     totalStress,
                                     N,
                                     bodyForce,
                                     reinterpret_cast< real64 (&)[numNodesPerElem][3] >(stack.localResidual) );

    stiffness.template upperBTDB< numNodesPerElem >( dNdX, -detJxW, stack.localJacobian );

    for( integer a = 0; a < numNodesPerElem; ++a )
    {
      stack.localDispFlowJacobian[a*3+0][0] += dNdX[a][0] * m_biotCoefficient * detJxW;
      stack.localDispFlowJacobian[a*3+1][0] += dNdX[a][1] * m_biotCoefficient * detJxW;
      stack.localDispFlowJacobian[a*3+2][0] += dNdX[a][2] * m_biotCoefficient * detJxW;
    }

    if( m_gravityAcceleration > 0.0 )
    {
      // Assumptions: (  i) dMixtureDens_dVolStrain contribution is neglected
      //              ( ii) grains are assumed incompressible
      //              (iii) TODO add dMixtureDens_dPressure and dMixtureDens_dGlobalCompDensity
    }

    // --- Mass balance accumulation
    // --- --- sum contributions to component accumulation from each phase

    // --- --- temporary work arrays
    real64 dPhaseAmount_dC[numMaxComponentsMultiphasePoroelastic];
    real64 dPhaseCompFrac_dC[numMaxComponentsMultiphasePoroelastic];
    real64 componentAmount[numMaxComponentsMultiphasePoroelastic] = { 0.0 };

    for( localIndex ip = 0; ip < NP; ++ip )
    {
      real64 const phaseAmountNew = porosityNew * m_fluidPhaseSaturation( k, ip ) * m_fluidPhaseDensity( k, q, ip );
      real64 const phaseAmountOld = porosityOld * m_fluidPhaseSaturationOld( k, ip ) * m_fluidPhaseDensityOld( k, ip );

      real64 const dPhaseAmount_dP = dPorosity_dPressure * m_fluidPhaseSaturation( k, ip ) * m_fluidPhaseDensity( k, q, ip )
                                     + porosityNew * (m_dFluidPhaseSaturation_dPressure( k, ip ) * m_fluidPhaseDensity( k, q, ip )
                                                      + m_fluidPhaseSaturation( k, ip ) * m_dFluidPhaseDensity_dPressure( k, q, ip ) );

      // assemble density dependence
      applyChainRule( NC,
                      m_dGlobalCompFraction_dGlobalCompDensity[k],
                      m_dFluidPhaseDensity_dGlobalCompFraction[k][q][ip],
                      dPhaseAmount_dC );

      for( localIndex jc = 0; jc < NC; ++jc )
      {
        dPhaseAmount_dC[jc] = dPhaseAmount_dC[jc] * m_fluidPhaseSaturation( k, ip )
                              + m_fluidPhaseDensity( k, q, ip ) * m_dFluidPhaseSaturation_dGlobalCompDensity( k, ip, jc );
        dPhaseAmount_dC[jc] *= porosityNew;
      }

      // ic - index of component whose conservation equation is assembled
      // (i.e. row number in local matrix)
      real64 const fluidPhaseDensityTimesFluidPhaseSaturation = m_fluidPhaseDensity( k, q, ip ) * m_fluidPhaseSaturation( k, ip );
      for( localIndex ic = 0; ic < NC; ++ic )
      {
        real64 const phaseCompAmountNew = phaseAmountNew * m_fluidPhaseCompFrac( k, q, ip, ic );
        real64 const phaseCompAmountOld = phaseAmountOld * m_fluidPhaseCompFracOld( k, ip, ic );

        real64 const dPhaseCompAmount_dP = dPhaseAmount_dP * m_fluidPhaseCompFrac( k, q, ip, ic )
                                           + phaseAmountNew * m_dFluidPhaseCompFrac_dPressure( k, q, ip, ic );

        componentAmount[ic] = fluidPhaseDensityTimesFluidPhaseSaturation * m_fluidPhaseCompFrac( k, q, ip, ic );

        stack.localFlowResidual[ic] += ( phaseCompAmountNew - phaseCompAmountOld ) * detJxW;
        stack.localFlowFlowJacobian[ic][0] += dPhaseCompAmount_dP * detJxW;;

        // jc - index of component w.r.t. whose compositional var the derivative is being taken
        // (i.e. col number in local matrix)

        // assemble phase composition dependence
        applyChainRule( NC,
                        m_dGlobalCompFraction_dGlobalCompDensity[k],
                        m_dFluidPhaseCompFraction_dGlobalCompFraction[k][q][ip][ic],
                        dPhaseCompFrac_dC );

        for( localIndex jc = 0; jc < NC; ++jc )
        {
          real64 const dPhaseCompAmount_dC = dPhaseCompFrac_dC[jc] * phaseAmountNew
                                             + m_fluidPhaseCompFrac( k, q, ip, ic ) * dPhaseAmount_dC[jc];
          stack.localFlowFlowJacobian[ic][jc+1] += dPhaseCompAmount_dC  * detJxW;
        }
      }
    }
    for( localIndex ic = 0; ic < m_numComponents; ++ic )
    {
      for( integer a = 0; a < numNodesPerElem; ++a )
      {
        stack.localFlowDispJacobian[ic][a*3+0] += dPorosity_dVolStrainIncrement * componentAmount[ic] * dNdX[a][0] * detJxW;
        stack.localFlowDispJacobian[ic][a*3+1] += dPorosity_dVolStrainIncrement * componentAmount[ic] * dNdX[a][1] * detJxW;
        stack.localFlowDispJacobian[ic][a*3+2] += dPorosity_dVolStrainIncrement * componentAmount[ic] * dNdX[a][2] * detJxW;
      }
    }
  }

  /**
   * @copydoc geosx::finiteElement::ImplicitKernelBase::complete
   */
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  real64 complete( localIndex const k,
                   StackVariables & stack ) const
  {
    GEOSX_UNUSED_VAR( k );
    real64 maxForce = 0;

    CONSTITUTIVE_TYPE::KernelWrapper::DiscretizationOps::template fillLowerBTDB< numNodesPerElem >( stack.localJacobian );

    //int nFlowDof = m_numComponents + 1;

    constexpr int nUDof = numNodesPerElem * numDofPerTestSupportPoint;

    for( int localNode = 0; localNode < numNodesPerElem; ++localNode )
    {
      for( int dim = 0; dim < numDofPerTestSupportPoint; ++dim )
      {
        localIndex const dof = LvArray::integerConversion< localIndex >( stack.localRowDofIndex[numDofPerTestSupportPoint * localNode + dim] - m_dofRankOffset );
        if( dof < 0 || dof >= m_matrix.numRows() ) continue;
        m_matrix.template addToRowBinarySearchUnsorted< parallelDeviceAtomic >( dof,
                                                                                stack.localRowDofIndex,
                                                                                stack.localJacobian[numDofPerTestSupportPoint * localNode + dim],
                                                                                numNodesPerElem * numDofPerTrialSupportPoint );

        RAJA::atomicAdd< parallelDeviceAtomic >( &m_rhs[dof], stack.localResidual[numDofPerTestSupportPoint * localNode + dim] );
        maxForce = fmax( maxForce, fabs( stack.localResidual[numDofPerTestSupportPoint * localNode + dim] ) );

        m_matrix.template addToRowBinarySearchUnsorted< parallelDeviceAtomic >( dof,
                                                                                stack.localFlowDofIndex,
                                                                                stack.localDispFlowJacobian[numDofPerTestSupportPoint * localNode + dim],
                                                                                m_numComponents + 1 );

      }
    }

    localIndex const dof = LvArray::integerConversion< localIndex >( stack.localFlowDofIndex[0] - m_dofRankOffset );
    if( 0 <= dof && dof < m_matrix.numRows() )
    {
      for( localIndex i = 0; i < m_numComponents; ++i )
      {
        m_matrix.template addToRowBinarySearchUnsorted< serialAtomic >( dof + i,
                                                                        stack.localRowDofIndex,
                                                                        stack.localFlowDispJacobian[i],
                                                                        nUDof );
        m_matrix.template addToRow< serialAtomic >( dof + i,
                                                    stack.localFlowDofIndex,
                                                    stack.localFlowFlowJacobian[i],
                                                    m_numComponents + 1 );

        RAJA::atomicAdd< serialAtomic >( &m_rhs[dof+i], stack.localFlowResidual[i] );
      }
    }

    return maxForce;
  }



protected:
  /// The array containing the nodal position array.
  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const m_X;

  /// The rank-global displacement array.
  arrayView2d< real64 const, nodes::TOTAL_DISPLACEMENT_USD > const m_disp;

  /// The rank-global incremental displacement array.
  arrayView2d< real64 const, nodes::INCR_DISPLACEMENT_USD > const m_uhat;

  /// The gravity vector.
  real64 const m_gravityVector[3];
  real64 const m_gravityAcceleration;

  /// The rank global density
  arrayView2d< real64 const > const m_solidDensity;
  arrayView3d< real64 const > const m_fluidPhaseDensity;
  arrayView2d< real64 const > const m_fluidPhaseDensityOld;
  arrayView3d< real64 const > const m_dFluidPhaseDensity_dPressure;
  arrayView4d< real64 const > const m_dFluidPhaseDensity_dGlobalCompFraction;
  arrayView4d< real64 const > const m_fluidPhaseCompFrac;
  arrayView3d< real64 const > const m_fluidPhaseCompFracOld;
  arrayView4d< real64 const > const m_dFluidPhaseCompFrac_dPressure;

  arrayView3d< real64 const > const m_fluidPhaseMassDensity;
  arrayView3d< real64 const > const m_dFluidPhaseMassDensity_dPressure;
  arrayView4d< real64 const > const m_dFluidPhaseMassDensity_dGlobalCompFraction;

  arrayView2d< real64 const > const m_fluidPhaseSaturation;
  arrayView2d< real64 const > const m_fluidPhaseSaturationOld;
  arrayView2d< real64 const > const m_dFluidPhaseSaturation_dPressure;
  arrayView3d< real64 const > const m_dFluidPhaseSaturation_dGlobalCompFraction;
  arrayView3d< real64 const > const m_dFluidPhaseSaturation_dGlobalCompDensity;

  arrayView3d< real64 const > const m_dGlobalCompFraction_dGlobalCompDensity;

  arrayView5d< real64 const > const m_dFluidPhaseCompFraction_dGlobalCompFraction;

  /// The global degree of freedom number
  arrayView1d< globalIndex const > const m_flowDofNumber;

  /// The rank-global fluid pressure array.
  arrayView1d< real64 const > const m_fluidPressure;

  /// The rank-global delta-fluid pressure array.
  arrayView1d< real64 const > const m_deltaFluidPressure;

  /// The rank-global reference porosity array
  arrayView1d< real64 const > const m_poroRef;

  /// Number of components
  localIndex const m_numComponents;

  /// Number of phases
  localIndex const m_numPhases;

  /// Biot's coefficient
  real64 const m_biotCoefficient;
};


} // namespace PoroelasticKernels

} // namespace geosx

#include "finiteElement/kernelInterface/SparsityKernelBase.hpp"

#endif // GEOSX_PHYSICSSOLVERS_MULTIPHYSICS_MULTIPHASEPOROELASTICKERNEL_HPP_
