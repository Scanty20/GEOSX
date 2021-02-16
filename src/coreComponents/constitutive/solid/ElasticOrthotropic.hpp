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
 *  @file ElasticOrthotropic.hpp
 */

#ifndef GEOSX_CONSTITUTIVE_SOLID_ELASTICORTHOTROPIC_HPP_
#define GEOSX_CONSTITUTIVE_SOLID_ELASTICORTHOTROPIC_HPP_

#include "SolidBase.hpp"
#include "constitutive/ExponentialRelation.hpp"
#include "LvArray/src/tensorOps.hpp"
#include "SolidModelDiscretizationOpsOrthotropic.hpp"

namespace geosx
{

namespace constitutive
{

/**
 * @class ElasticOrthotropicUpdates
 *
 * Class to provide elastic orthotropic material updates that
 * may be called from a kernel function.
 */
class ElasticOrthotropicUpdates : public SolidBaseUpdates
{
public:
  using DiscretizationOps = SolidModelDiscretizationOpsOrthotropic;

  /**
   * @brief Constructor
   * @param[in] c11 The 11 component of the Voigt stiffness tensor.
   * @param[in] c12 The 12 component of the Voigt stiffness tensor.
   * @param[in] c13 The 13 component of the Voigt stiffness tensor.
   * @param[in] c22 The 22 component of the Voigt stiffness tensor.
   * @param[in] c23 The 23 component of the Voigt stiffness tensor.
   * @param[in] c33 The 33 component of the Voigt stiffness tensor.
   * @param[in] c44 The 44 component of the Voigt stiffness tensor.
   * @param[in] c55 The 55 component of the Voigt stiffness tensor.
   * @param[in] c66 The 66 component of the Voigt stiffness tensor.
   * @param[in] stress The ArrayView holding the stress data for each quadrature
   *                   point.
   */
  ElasticOrthotropicUpdates( arrayView1d< real64 const > const & c11,
                             arrayView1d< real64 const > const & c12,
                             arrayView1d< real64 const > const & c13,
                             arrayView1d< real64 const > const & c22,
                             arrayView1d< real64 const > const & c23,
                             arrayView1d< real64 const > const & c33,
                             arrayView1d< real64 const > const & c44,
                             arrayView1d< real64 const > const & c55,
                             arrayView1d< real64 const > const & c66,
                             arrayView3d< real64, solid::STRESS_USD > const & stress ):
    SolidBaseUpdates( stress ),
    m_c11( c11 ),
    m_c12( c12 ),
    m_c13( c13 ),
    m_c22( c22 ),
    m_c23( c23 ),
    m_c33( c33 ),
    m_c44( c44 ),
    m_c55( c55 ),
    m_c66( c66 )
  {}

  /// Deleted default constructor
  ElasticOrthotropicUpdates() = delete;

  /// Default copy constructor
  ElasticOrthotropicUpdates( ElasticOrthotropicUpdates const & ) = default;

  /// Default move constructor
  ElasticOrthotropicUpdates( ElasticOrthotropicUpdates && ) = default;

  /// Deleted copy assignment operator
  ElasticOrthotropicUpdates & operator=( ElasticOrthotropicUpdates const & ) = delete;

  /// Deleted move assignment operator
  ElasticOrthotropicUpdates & operator=( ElasticOrthotropicUpdates && ) =  delete;


  GEOSX_HOST_DEVICE
  virtual void smallStrainNoState( localIndex const k,
                                   real64 const ( &voigtStrain )[ 6 ],
                                   real64 ( &stress )[ 6 ] ) const override final;

  GEOSX_HOST_DEVICE
  virtual void smallStrain( localIndex const k,
                            localIndex const q,
                            real64 const ( &voigtStrainInc )[ 6 ] ) const override final;

  GEOSX_HOST_DEVICE
  virtual void hypoElastic( localIndex const k,
                            localIndex const q,
                            real64 const ( &Ddt )[ 6 ],
                            real64 const ( &Rot )[ 3 ][ 3 ] ) const override final;

  GEOSX_HOST_DEVICE
  virtual void hyperElastic( localIndex const k,
                             real64 const (&FmI)[3][3],
                             real64 ( &stress )[ 6 ] ) const override final;

  GEOSX_HOST_DEVICE
  virtual void hyperElastic( localIndex const k,
                             localIndex const q,
                             real64 const (&FmI)[3][3] ) const override final;

  GEOSX_FORCE_INLINE
  GEOSX_HOST_DEVICE
  virtual void getStiffness( localIndex const k,
                             localIndex const q,
                             real64 (& c)[6][6] ) const override final
  {
    GEOSX_UNUSED_VAR( q );
    memset( c, 0, sizeof( c ) );
    c[0][0] = m_c11[k];
    c[0][1] = m_c12[k];
    c[0][2] = m_c13[k];
    c[1][0] = c[0][1];
    c[1][1] = m_c22[k];
    c[1][2] = m_c23[k];
    c[2][0] = c[0][2];
    c[2][1] = c[1][2];
    c[2][2] = m_c33[k];
    c[3][3] = m_c44[k];
    c[4][4] = m_c55[k];
    c[5][5] = m_c66[k];
  }

  GEOSX_FORCE_INLINE
  GEOSX_HOST_DEVICE
  void setDiscretizationOps( localIndex const k,
                             localIndex const q,
                             DiscretizationOps & discOps ) const
  {
    GEOSX_UNUSED_VAR( q )
    discOps.m_c11 = m_c11[k];
    discOps.m_c13 = m_c13[k];
    discOps.m_c33 = m_c33[k];
    discOps.m_c44 = m_c44[k];
    discOps.m_c66 = m_c66[k];
  }


  GEOSX_HOST_DEVICE
  virtual real64 calculateStrainEnergyDensity( localIndex const k,
                                               localIndex const q ) const override final
  {
    GEOSX_UNUSED_VAR( k, q );
    GEOSX_ERROR( "Not implemented" );
    return 0;
  }

private:
  /// A reference to the ArrayView holding c11 for each element.
  arrayView1d< real64 const > const m_c11;

  /// A reference to the ArrayView holding c12 for each element.
  arrayView1d< real64 const > const m_c12;

  /// A reference to the ArrayView holding c13 for each element.
  arrayView1d< real64 const > const m_c13;

  /// A reference to the ArrayView holding c22 for each element.
  arrayView1d< real64 const > const m_c22;

  /// A reference to the ArrayView holding c23 for each element.
  arrayView1d< real64 const > const m_c23;

  /// A reference to the ArrayView holding c33 for each element.
  arrayView1d< real64 const > const m_c33;

  /// A reference to the ArrayView holding c44 for each element.
  arrayView1d< real64 const > const m_c44;

  /// A reference to the ArrayView holding c55 for each element.
  arrayView1d< real64 const > const m_c55;

  /// A reference to the ArrayView holding c66 for each element.
  arrayView1d< real64 const > const m_c66;
};


GEOSX_FORCE_INLINE
GEOSX_HOST_DEVICE
void
ElasticOrthotropicUpdates::
  smallStrainNoState( localIndex const k,
                      real64 const ( &voigtStrain )[ 6 ],
                      real64 ( & stress )[ 6 ] ) const
{
  real64 const c12temp = ( m_c11[k] - 2.0 * m_c66[k] );
  stress[0] = m_c11[k] * voigtStrain[0] +  c12temp * voigtStrain[1] + m_c13[k]*voigtStrain[2];
  stress[1] =  c12temp * voigtStrain[0] + m_c11[k] * voigtStrain[1] + m_c13[k]*voigtStrain[2];
  stress[2] = m_c13[k] * voigtStrain[0] + m_c13[k] * voigtStrain[1] + m_c33[k]*voigtStrain[2];

  stress[3] = m_c44[k]*voigtStrain[3];
  stress[4] = m_c44[k]*voigtStrain[4];
  stress[5] = m_c66[k]*voigtStrain[5];
}


GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void
ElasticOrthotropicUpdates::
  smallStrain( localIndex const k,
               localIndex const q,
               real64 const ( &voigtStrainInc )[ 6 ] ) const
{
  real64 const temp = m_c11[ k ] * ( voigtStrainInc[ 0 ] + voigtStrainInc[ 1 ] ) + m_c13[ k ] * voigtStrainInc[ 2 ];
  m_stress( k, q, 0 ) += -2.0 * m_c66[ k ] * voigtStrainInc[ 1 ] + temp;
  m_stress( k, q, 1 ) += -2.0 * m_c66[ k ] * voigtStrainInc[ 0 ] + temp;
  m_stress( k, q, 2 ) = m_stress( k, q, 2 ) + m_c13[ k ] * ( voigtStrainInc[ 0 ] + voigtStrainInc[ 1 ] ) + m_c33[ k ] * voigtStrainInc[ 2 ];
  m_stress( k, q, 3 ) = m_stress( k, q, 3 ) + m_c44[ k ] * voigtStrainInc[ 3 ];
  m_stress( k, q, 4 ) = m_stress( k, q, 4 ) + m_c44[ k ] * voigtStrainInc[ 4 ];
  m_stress( k, q, 5 ) = m_stress( k, q, 5 ) + m_c66[ k ] * voigtStrainInc[ 5 ];
}

GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void
ElasticOrthotropicUpdates::
  hypoElastic( localIndex const k,
               localIndex const q,
               real64 const ( &Ddt )[ 6 ],
               real64 const ( &Rot )[ 3 ][ 3 ] ) const
{
  smallStrain( k, q, Ddt );
  real64 temp[ 6 ];
  LvArray::tensorOps::Rij_eq_AikSymBklAjl< 3 >( temp, Rot, m_stress[ k ][ q ] );
  LvArray::tensorOps::copy< 6 >( m_stress[ k ][ q ], temp );
}

GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void
ElasticOrthotropicUpdates::
  hyperElastic( localIndex const GEOSX_UNUSED_PARAM( k ),
                real64 const (&GEOSX_UNUSED_PARAM( FmI ))[3][3],
                real64 ( & )[ 6 ] ) const
{
  GEOSX_ERROR( "ElasticOrthotropicKernelWrapper::HyperElastic() is not implemented!" );
}

GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void
ElasticOrthotropicUpdates::
  hyperElastic( localIndex const GEOSX_UNUSED_PARAM( k ),
                localIndex const GEOSX_UNUSED_PARAM( q ),
                real64 const (&GEOSX_UNUSED_PARAM( FmI ))[3][3] ) const
{
  GEOSX_ERROR( "ElasticOrthotropicKernelWrapper::HyperElastic() is not implemented!" );
}



/**
 * @class ElasticOrthotropic
 *
 * Class to provide a linear elastic transverse isotropic material response.
 */
class ElasticOrthotropic : public SolidBase
{
public:

  /// @typedef Alias for ElasticOrthotropicUpdates
  using KernelWrapper = ElasticOrthotropicUpdates;

  /**
   * @brief constructor
   * @param[in]name name of the instance in the catalog
   * @param[in]parent the group which contains this instance
   */
  ElasticOrthotropic( string const & name, Group * const parent );

  /**
   * Destructor
   */
  virtual ~ElasticOrthotropic() override;

  /**
   * @name Static Factory Catalog members and functions
   */
  ///@{

  /// string name to use for this class in the catalog
  static constexpr auto m_catalogNameString = "ElasticOrthotropic";

  /**
   * @return A string that is used to register/lookup this class in the registry
   */
  static string catalogName() { return m_catalogNameString; }

  virtual string getCatalogName() const override { return catalogName(); }
  ///@}

  /**
   * @struct Set of "char const *" and keys for data specified in this class.
   */
  struct viewKeyStruct : public SolidBase::viewKeyStruct
  {
    /// string/key for transverse youngs modulus
    static constexpr auto defaultYoungsModulusTransverse     = "defaultYoungsModulusTransverse";

    /// string/key for axial Young's modulus
    static constexpr auto defaultYoungsModulusAxial          = "defaultYoungsModulusAxial";

    /// string/key for transverse Poisson's Ratio
    static constexpr auto defaultPoissonRatioTransverse      = "defaultPoissonRatioTransverse";

    /// string/key for axial Poisson's Ratio
    static constexpr auto defaultPoissonRatioAxialTransverse = "defaultPoissonRatioAxialTransverse";

    /// string/key for transverse shear modulus
    static constexpr auto defaultShearModulusAxialTransverse = "defaultShearModulusAxialTransverse";

    /// string/key for default c11 component of Voigt stiffness tensor
    static constexpr auto defaultC11 = "defaultC11";

    /// string/key for default c12 component of Voigt stiffness tensor
    static constexpr auto defaultC12 = "defaultC12";

    /// string/key for default c13 component of Voigt stiffness tensor
    static constexpr auto defaultC13 = "defaultC13";

    /// string/key for default c22 component of Voigt stiffness tensor
    static constexpr auto defaultC22 = "defaultC22";

    /// string/key for default c23 component of Voigt stiffness tensor
    static constexpr auto defaultC23 = "defaultC23";

    /// string/key for default c33 component of Voigt stiffness tensor
    static constexpr auto defaultC33 = "defaultC33";

    /// string/key for default c44 component of Voigt stiffness tensor
    static constexpr auto defaultC44 = "defaultC44";

    /// string/key for default c55 component of Voigt stiffness tensor
    static constexpr auto defaultC55 = "defaultC55";

    /// string/key for default c66 component of Voigt stiffness tensor
    static constexpr auto defaultC66 = "defaultC66";

    /// string/key for c11 component of Voigt stiffness tensor
    static constexpr auto c11 = "c11";

    /// string/key for c12 component of Voigt stiffness tensor
    static constexpr auto c12 = "c12";

    /// string/key for c13 component of Voigt stiffness tensor
    static constexpr auto c13 = "c13";

    /// string/key for c22 component of Voigt stiffness tensor
    static constexpr auto c22 = "c22";

    /// string/key for c23 component of Voigt stiffness tensor
    static constexpr auto c23 = "c23";

    /// string/key for c33 component of Voigt stiffness tensor
    static constexpr auto c33 = "c33";

    /// string/key for c44 component of Voigt stiffness tensor
    static constexpr auto c44 = "c44";

    /// string/key for c55 component of Voigt stiffness tensor
    static constexpr auto c55 = "c55";

    /// string/key for c66 component of Voigt stiffness tensor
    static constexpr auto c66 = "c66";
  };

  /**
   * @brief Getter for default transverse Young's modulus
   * @return The value of the default transverse Young's modulus.
   */
  real64 getDefaultYoungsModulusTransverse()
  {
    return m_defaultYoungsModulusTransverse;
  }

  /**
   * @brief Setter for the default transverse Young's modulus.
   * @param[in] input New value for the default transverse Young's modulus
   */
  void setDefaultYoungsModulusTransverse( real64 const input )
  {
    m_defaultYoungsModulusTransverse = input;
  }

  /**
   * @brief Getter for default axial Young's modulus
   * @return The value of the default axial Young's modulus.
   */
  real64 getDefaultYoungsModulusAxial()
  {
    return m_defaultYoungsModulusAxial;
  }

  /**
   * @brief Setter for the default axial Young's modulus.
   * @param[in] input New value for the default axial Young's modulus
   */
  void setDefaultYoungsModulusAxial( real64 const input )
  {
    m_defaultYoungsModulusAxial = input;
  }


  /**
   * @brief Getter for default transverse Poisson's ratio
   * @return The value of the default transverse Poisson's ratio.
   */
  real64 getDefaultPoissonsRatioTransverse()
  {
    return m_defaultPoissonTransverse;
  }

  /**
   * @brief Setter for the default transverse Poisson's ratio.
   * @param[in] input New value for the default transverse Poisson's ratio
   */
  void setDefaultPoissonsRatioTransverse( real64 const input )
  {
    m_defaultPoissonTransverse = input;
  }


  /**
   * @brief Getter for default axial Poisson's ratio
   * @return The value of the default axial/transverse Poisson's modulus.
   */
  real64 getDefaultPoissonsRatioAxialTransverse()
  {
    return m_defaultPoissonAxialTransverse;
  }

  /**
   * @brief Setter for the default axial Poisson's modulus.
   * @param[in] input New value for the default axial/transverse Poisson's
   *             modulus
   */
  void setDefaultPoissonsRatioAxialTransverse( real64 const input )
  {
    m_defaultPoissonAxialTransverse = input;
  }


  /**
   * @brief Getter for default axial/transverse Shear modulus
   * @return The value of the default axial/transverse Shear modulus.
   */
  real64 getDefaultShearModulusAxialTransverse()
  {
    return m_defaultShearModulusAxialTransverse;
  }

  /**
   * @brief Setter for the default axial/transverse Shear modulus.
   * @param[in] input New value for the default axial/transverse Shear modulus
   */
  void setDefaultShearModulusAxialTransverse( real64 const input )
  {
    m_defaultShearModulusAxialTransverse = input;
  }

  /**
   * @brief Const-Getter for 11 component of Voigt stiffness tensor.
   * @return reference to immutable 11 component of Voigt stiffness tensor.
   */
  arrayView1d< real64 const > getC11() const { return m_c11; }

  /**
   * @brief Getter for 11 component of Voigt stiffness tensor.
   * @return reference to mutable 11 component of Voigt stiffness tensor.
   */
  arrayView1d< real64 > getC11() { return m_c11; }


  /**
   * @brief Const-Getter for 13 component of Voigt stiffness tensor.
   * @return reference to immutable 13 component of Voigt stiffness tensor.
   */
  arrayView1d< real64 const > getC13() const { return m_c13; }

  /**
   * @brief Getter for 13 component of Voigt stiffness tensor.
   * @return reference to mutable 13 component of Voigt stiffness tensor.
   */
  arrayView1d< real64 > getC13() { return m_c13; }

  /**
   * @brief Const-Getter for 33 component of Voigt stiffness tensor.
   * @return reference to immutable 33 component of Voigt stiffness tensor.
   */
  arrayView1d< real64 const > getC33() const { return m_c33; }

  /**
   * @brief Getter for 33 component of Voigt stiffness tensor.
   * @return reference to mutable 33 component of Voigt stiffness tensor.
   */
  arrayView1d< real64 > getC33() { return m_c33; }

  /**
   * @brief Const-Getter for 44 component of Voigt stiffness tensor.
   * @return reference to immutable 44 component of Voigt stiffness tensor.
   */
  arrayView1d< real64 const > getC44() const { return m_c44; }

  /**
   * @brief Getter for 44 component of Voigt stiffness tensor.
   * @return reference to mutable 44 component of Voigt stiffness tensor.
   */
  arrayView1d< real64 > getC44() { return m_c44; }

  /**
   * @brief Const-Getter for 66 component of Voigt stiffness tensor.
   * @return reference to immutable 66 component of Voigt stiffness tensor.
   */
  arrayView1d< real64 const > getC66() const { return m_c66; }

  /**
   * @brief Getter for 66 component of Voigt stiffness tensor.
   * @return reference to mutable 66 component of Voigt stiffness tensor.
   */
  arrayView1d< real64 > getC66() { return m_c66; }

  /**
   * @brief Create a instantiation of the
   *        ElasticOrthotropicUpdates class that refers to the
   *        data in this.
   * @return An instantiation of ElasticOrthotropicUpdates.
   */
  ElasticOrthotropicUpdates createKernelUpdates()
  {
    return ElasticOrthotropicUpdates( m_c11,
                                      m_c12,
                                      m_c13,
                                      m_c22,
                                      m_c23,
                                      m_c33,
                                      m_c44,
                                      m_c55,
                                      m_c66,
                                      m_stress );
  }

  /**
   * @brief Construct an update kernel for a derived type.
   * @tparam UPDATE_KERNEL The type of update kernel from the derived type.
   * @tparam PARAMS The parameter pack to hold the constructor parameters for
   *   the derived update kernel.
   * @param constructorParams The constructor parameter for the derived type.
   * @return An @p UPDATE_KERNEL object.
   */
  template< typename UPDATE_KERNEL, typename ... PARAMS >
  UPDATE_KERNEL createDerivedKernelUpdates( PARAMS && ... constructorParams )
  {
    return UPDATE_KERNEL( std::forward< PARAMS >( constructorParams )...,
                          m_c11,
                          m_c12,
                          m_c13,
                          m_c22,
                          m_c23,
                          m_c33,
                          m_c44,
                          m_c55,
                          m_c66,
                          m_stress );
  }


protected:
  virtual void postProcessInput() override;

private:

  /// The default value of the transverse Young's modulus for any new
  /// allocations.
  real64 m_defaultYoungsModulusTransverse;

  /// The default value of the axial Young's modulus for any new
  /// allocations.
  real64 m_defaultYoungsModulusAxial;

  /// The default value of the transverse Poisson's ratio for any new
  /// allocations.
  real64 m_defaultPoissonTransverse;

  /// The default value of the axial/transverse Poisson's ratio for any new
  /// allocations.
  real64 m_defaultPoissonAxialTransverse;

  /// The default value of the axial/transverse Shear modulus for any new
  /// allocations.
  real64 m_defaultShearModulusAxialTransverse;

  /// The default value of the 11 component of the Voigt stiffness tensor
  /// for any new allocations.
  real64 m_defaultC11;

  /// The default value of the 12 component of the Voigt stiffness tensor
  /// for any new allocations.
  real64 m_defaultC12;

  /// The default value of the 13 component of the Voigt stiffness tensor
  /// for any new allocations.
  real64 m_defaultC13;

  /// The default value of the 22 component of the Voigt stiffness tensor
  /// for any new allocations.
  real64 m_defaultC22;

  /// The default value of the 23 component of the Voigt stiffness tensor
  /// for any new allocations.
  real64 m_defaultC23;

  /// The default value of the 33 component of the Voigt stiffness tensor
  /// for any new allocations.
  real64 m_defaultC33;

  /// The default value of the 44 component of the Voigt stiffness tensor
  /// for any new allocations.
  real64 m_defaultC44;

  /// The default value of the 55 component of the Voigt stiffness tensor
  /// for any new allocations.
  real64 m_defaultC55;

  /// The default value of the 66 component of the Voigt stiffness tensor
  /// for any new allocations.
  real64 m_defaultC66;

  /// The 11 component of the Voigt stiffness tensor.
  array1d< real64 > m_c11;

  /// The 12 component of the Voigt stiffness tensor.
  array1d< real64 > m_c12;

  /// The 13 component of the Voigt stiffness tensor.
  array1d< real64 > m_c13;

  /// The 22 component of the Voigt stiffness tensor.
  array1d< real64 > m_c22;

  /// The 23 component of the Voigt stiffness tensor.
  array1d< real64 > m_c23;

  /// The 33 component of the Voigt stiffness tensor.
  array1d< real64 > m_c33;

  /// The 44 component of the Voigt stiffness tensor.
  array1d< real64 > m_c44;

  /// The 55 component of the Voigt stiffness tensor.
  array1d< real64 > m_c55;

  /// The 66 component of the Voigt stiffness tensor.
  array1d< real64 > m_c66;

};

} /* namespace constitutive */

} /* namespace geosx */

#endif /* GEOSX_CONSTITUTIVE_SOLID_ELASTICORTHOTROPIC_HPP_ */

