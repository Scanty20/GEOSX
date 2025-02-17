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
 * @file CO2Solubility.cpp
 */

#include "constitutive/fluid/PVTFunctions/CO2Solubility.hpp"

#include "constitutive/fluid/PVTFunctions/PVTFunctionHelpers.hpp"
#include "functions/FunctionManager.hpp"

namespace geosx
{

using namespace stringutilities;

namespace constitutive
{

namespace PVTProps
{

constexpr real64 T_K_f  = 273.15;
constexpr real64 P_Pa_f = 1e+5;
constexpr real64 P_c    = 73.773 * P_Pa_f;
constexpr real64 T_c    = 304.1282;
constexpr real64 Rgas   = 8.314467;
constexpr real64 V_c    = Rgas*T_c/P_c;

// these coefficients are in Table (A1) of Duan and Sun (2003)
constexpr real64 acoef[] =
{ 8.99288497e-2, -4.94783127e-1, 4.77922245e-2, 1.03808883e-2, -2.82516861e-2, 9.49887563e-2, 5.20600880e-4,
  -2.93540971e-4, -1.77265112e-3, -2.51101973e-5, 8.93353441e-5, 7.88998563e-5, -1.66727022e-2, 1.398, 2.96e-2 };

namespace detail
{

real64 ff( real64 const & T, real64 const & P, real64 const & V_r )
{
  // reduced pressure
  real64 const P_r = P*P_Pa_f/P_c;
  // reduced temperature
  real64 const T_r = (T_K_f+T)/T_c;

  // CO2 equation of state
  // see equation (A1) in Duan and Sun (2003)
  real64 const f_Z = 1.0
                     + ( acoef[0] + acoef[1]/(T_r * T_r) + acoef[2]/(T_r * T_r * T_r) )/V_r
                     + ( acoef[3] + acoef[4]/(T_r * T_r) + acoef[5]/(T_r * T_r * T_r) )/(V_r*V_r)
                     + ( acoef[6] + acoef[7]/(T_r * T_r) + acoef[8]/(T_r * T_r * T_r) )/(V_r*V_r*V_r*V_r)
                     + ( acoef[9] + acoef[10]/(T_r * T_r) + acoef[11]/(T_r * T_r * T_r) )/(V_r*V_r*V_r*V_r*V_r)
                     + acoef[12]/(T_r * T_r * T_r)/(V_r * V_r) * (acoef[13] + acoef[14]/(V_r * V_r)) * exp( -acoef[14]/(V_r * V_r)) - P_r * V_r / T_r;

  return f_Z;
}

real64 PWater( real64 const & T )
{
  // these coefficients are defined in Table (B1) of Duan and Sun (2003)
  constexpr real64 ccoef[] = { -38.640844, 5.8948420, 59.876516, 26.654627, 10.637097 };

  // H2O critical pressure (bars)
  real64 const P_c_w = 220.85;
  // H2O critical temperature (K)
  real64 const T_c_w = 647.29;
  real64 const tt = ( (T+T_K_f)-T_c_w )/T_c_w;
  // Empirical model for water pressure of equation (B1) of Duan and Sun (2003)
  real64 const x = (P_c_w*(T+T_K_f)/T_c_w)
                   * (1
                      + ccoef[0]*pow( -tt, 1.9 )
                      + ccoef[1]*tt
                      + ccoef[2]*tt*tt
                      + ccoef[3]*tt*tt*tt
                      + ccoef[4]*tt*tt*tt*tt);

  return x;
}

real64 logF( real64 const & T, real64 const & P, real64 const & V_r )
{
  // reduced pressure
  real64 const P_r = P*P_Pa_f/P_c;
  // reduced temperature
  real64 const T_r = (T_K_f+T)/T_c;
  real64 const Z   = P_r * V_r/T_r;

  // fugacity coefficient of CO2, equation (A6) of Duan and Sun (2003)
  real64 const log_f = Z - 1 - log( Z ) +
                       ( acoef[0] + acoef[1]/T_r/T_r + acoef[2]/T_r/T_r/T_r )/V_r
                       + ( acoef[3] + acoef[4]/T_r/T_r + acoef[5]/T_r/T_r/T_r )/2.0/V_r/V_r
                       + ( acoef[6] + acoef[7]/T_r/T_r + acoef[8]/T_r/T_r/T_r )/4.0/V_r/V_r/V_r/V_r
                       + ( acoef[9] + acoef[10]/T_r/T_r + acoef[11]/T_r/T_r/T_r )/5.0/V_r/V_r/V_r/V_r/V_r
                       + acoef[12]/2.0/T_r/T_r/T_r/acoef[14] * ( acoef[13] + 1.0 - (acoef[13] + 1.0 + acoef[14]/V_r/V_r) * exp( -acoef[14]/V_r/V_r ) );

  return log_f;
}

real64 Par( real64 const & T, real64 const & P, real64 const * cc )
{
  // "equation for the parameters", see equation (7) of Duan and Sun (2003)
  real64 x = cc[0]
             + cc[1]*T
             + cc[2]/T
             + cc[3]*T*T
             + cc[4]/(630.0-T)
             + cc[5]*P
             + cc[6]*P *log( T )
             + cc[7]*P/T
             + cc[8]*P/(630.0-T)
             + cc[9]*P*P/(630.0-T)/(630.0-T)
             + cc[10]*T *log( P );

  return x;
}

}

CO2Solubility::CO2Solubility( string_array const & inputPara,
                              string_array const & phaseNames,
                              string_array const & componentNames,
                              array1d< real64 > const & componentMolarWeight ):
  FlashModelBase( inputPara[1],
                  componentNames,
                  componentMolarWeight )
{
  GEOSX_THROW_IF( phaseNames.size() != 2,
                  "The CO2Solubility model is a two-phase model",
                  InputError );
  GEOSX_THROW_IF( componentNames.size() != 2,
                  "The CO2Solubility model is a two-component model",
                  InputError );

  char const * expectedCO2ComponentNames[] = { "CO2", "co2" };
  m_CO2Index = PVTFunctionHelpers::findName( componentNames, expectedCO2ComponentNames );
  GEOSX_THROW_IF( m_CO2Index < 0 || m_CO2Index >= componentNames.size(),
                  "Component CO2 is not found!",
                  InputError );

  char const * expectedWaterComponentNames[] = { "Water", "water" };
  m_waterIndex = PVTFunctionHelpers::findName( componentNames, expectedWaterComponentNames );
  GEOSX_THROW_IF( m_waterIndex < 0 || m_waterIndex >= componentNames.size(),
                  "Component Water/Brine is not found!",
                  InputError );

  char const * expectedGasPhaseNames[] = { "CO2", "co2", "gas", "Gas" };
  m_phaseGasIndex = PVTFunctionHelpers::findName( phaseNames,
                                                  expectedGasPhaseNames );
  GEOSX_THROW_IF( m_phaseGasIndex < 0 || m_phaseGasIndex >= phaseNames.size(),
                  "Phase co2/gas is not found!",
                  InputError );

  char const * expectedWaterPhaseNames[] = { "Water", "water", "Liquid", "liquid" };
  m_phaseLiquidIndex = PVTFunctionHelpers::findName( phaseNames,
                                                     expectedWaterPhaseNames );
  GEOSX_THROW_IF( m_phaseLiquidIndex < 0 || m_phaseLiquidIndex >= phaseNames.size(),
                  "Phase water/liquid is not found!",
                  InputError );

  makeTable( inputPara );
}

void CO2Solubility::makeTable( string_array const & inputPara )
{
  // initialize the (p,T) coordinates
  PTTableCoordinates tableCoords;
  PVTFunctionHelpers::initializePropertyTable( inputPara, tableCoords );

  // initialize salinity and tolerance
  GEOSX_THROW_IF( inputPara.size() < 9,
                  "Invalid property input!",
                  InputError );

  real64 tolerance = 1e-9;
  real64 salinity = 0.0;
  try
  {
    salinity = stod( inputPara[8] );
    if( inputPara.size() >= 10 )
    {
      tolerance = stod( inputPara[9] );
    }
  }
  catch( const std::invalid_argument & e )
  {
    GEOSX_THROW( "Invalid property argument:" + string( e.what()),
                 InputError );
  }

  array1d< real64 > values( tableCoords.nPressures() * tableCoords.nTemperatures() );
  calculateCO2Solubility( tolerance, tableCoords, salinity, values );

  FunctionManager & functionManager = FunctionManager::getInstance();
  m_CO2SolubilityTable = dynamicCast< TableFunction * >( functionManager.createChild( "TableFunction", "CO2SolubilityTable" ) );
  m_CO2SolubilityTable->setTableCoordinates( tableCoords.getCoords() );
  m_CO2SolubilityTable->setTableValues( values );
  m_CO2SolubilityTable->reInitializeFunction();
  m_CO2SolubilityTable->setInterpolationMethod( TableFunction::InterpolationType::Linear );
}

void CO2Solubility::calculateCO2Solubility( real64 const & tolerance,
                                            PTTableCoordinates const & tableCoords,
                                            real64 const & salinity,
                                            array1d< real64 > const & values )
{
  // Interaction parameters, see Table 2 of Duan and Sun (2003)
  constexpr real64 mu[] =
  { 28.9447706, -0.0354581768, -4770.67077, 1.02782768e-5, 33.8126098, 9.04037140e-3,
    -1.14934031e-3, -0.307405726, -0.0907301486, 9.32713393e-4, 0 };
  constexpr real64 lambda[] = { -0.411370585, 6.07632013e-4, 97.5347708, 0, 0, 0, 0, -0.0237622469, 0.0170656236, 0, 1.41335834e-5 };
  constexpr real64 zeta[] = { 3.36389723e-4, -1.98298980e-5, 0, 0, 0, 0, 0, 2.12220830e-3, -5.24873303e-3, 0, 0 };

  localIndex const nPressures = tableCoords.nPressures();
  localIndex const nTemperatures = tableCoords.nTemperatures();

  for( localIndex i = 0; i < nPressures; ++i )
  {
    real64 const P = tableCoords.getPressure( i ) / P_Pa_f;

    for( localIndex j = 0; j < nTemperatures; ++j )
    {
      real64 const T = tableCoords.getTemperature( j );

      // compute reduced volume by solving the CO2 equation of state
      real64 V_r = 0.0;
      CO2SolubilityFunction( tolerance, T, P, V_r, &detail::ff );

      // compute equation (6) of Duan and Sun (2003)
      real64 const logK = detail::Par( T+T_K_f, P, mu )
                          - detail::logF( T, P, V_r )
                          + 2*detail::Par( T+T_K_f, P, lambda ) * salinity
                          + detail::Par( T+T_K_f, P, zeta ) * salinity * salinity;

      // mole fraction of CO2 in vapor phase, equation (4) of Duan and Sun (2003)
      real64 const y_CO2 = (P - detail::PWater( T ))/P;
      values[j*nPressures+i] = y_CO2 * P / exp( logK );
    }
  }
}

void CO2Solubility::CO2SolubilityFunction( real64 const & tolerance,
                                           real64 const & T,
                                           real64 const & P,
                                           real64 & V_r,
                                           real64 (*f)( real64 const & x1, real64 const & x2, real64 const & x3 ) )
{
  constexpr real64 dx = 1e-10;
  int count = 0;
  real64 Vr_int = 0.05;

  V_r = 0.75*Rgas*(T_K_f+T)/(P*P_Pa_f)*(1/V_c);

  // iterate until the solution of the CO2 equation of state is found
  for(;; )
  {
    if( V_r < 0.0 )
    {
      V_r = Vr_int;
      Vr_int += 0.05;
    }

    real64 const v0 = (*f)( T, P, V_r );
    real64 const v1 = (*f)( T, P, V_r+dx );
    real64 const dre = -v0/((v1-v0)/dx);

    if( fabs( dre ) < tolerance )
    {
      break;
    }

    GEOSX_THROW_IF( count > 50,
                    "CO2Solubility NR convergence fails! " << "dre = " << dre << ", tolerance = " << tolerance,
                    InputError );

    count++;
    V_r += dre;
  }
}


CO2Solubility::KernelWrapper CO2Solubility::createKernelWrapper()
{
  return KernelWrapper( m_componentNames,
                        m_componentMolarWeight,
                        m_CO2SolubilityTable,
                        m_CO2Index,
                        m_waterIndex,
                        m_phaseGasIndex,
                        m_phaseLiquidIndex );
}

REGISTER_CATALOG_ENTRY( FlashModelBase, CO2Solubility, string_array const &, string_array const &, string_array const &, array1d< real64 > const & )

} // end namespace PVTProps

} // namespace constitutive

} // end namespace geosx
