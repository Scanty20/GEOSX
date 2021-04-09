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

#ifndef GEOSX_PHYSICSSOLVERS_SIMPLEPDE_LAPLACE_BASE_HPP
#define GEOSX_PHYSICSSOLVERS_SIMPLEPDE_LAPLACE_BASE_HPP

#include "common/EnumStrings.hpp"   // facilities for enum-string conversion (for reading enum values from XML input)
#include "physicsSolvers/SolverBase.hpp"  // an abstraction class shared by all physics solvers

namespace geosx
{

// Like most physics solvers, the Laplace solver derives from a generic SolverBase class.
// The base class is densely Doxygen-commented and worth a look if you have not done so already.
// Most important system assembly steps, linear and non-linear resolutions, and time-stepping mechanisms
// are implemented at the SolverBase class level and can thus be used in Laplace without needing reimplementation.

//START_SPHINX_INCLUDE_02
class LaplaceBase : public SolverBase
{
public:
  // The default nullary constructor is disabled to avoid compiler auto-generation:
  LaplaceBase() = delete;

  // The constructor needs a user-defined "name" and a parent Group (to place this instance in the tree structure of classes)
  LaplaceBase( const string & name,
               Group * const parent );

  // Destructor
  virtual ~LaplaceBase() override;

  // This method ties properties with their supporting mesh
  virtual void registerDataOnMesh( Group & meshBodies ) override final;

//END_SPHINX_INCLUDE_02

  // Choice of transient treatment options (steady, backward, forward Euler scheme):
  //START_SPHINX_INCLUDE_01
  enum class TimeIntegrationOption : integer
  {
    SteadyState,
    ImplicitTransient,
    ExplicitTransient
  };
  //END_SPHINX_INCLUDE_01


  // This structure stores ``dataRepository::ViewKey`` objects
  // used as binding between the input XML tags and source
  // code variables (here, timeIntegrationOption and fieldVarName)
  //START_SPHINX_INCLUDE_04
  struct viewKeyStruct : public SolverBase::viewKeyStruct
  {
    dataRepository::ViewKey timeIntegrationOption = { "timeIntegrationOption" };
    dataRepository::ViewKey fieldVarName = { "fieldName" };
  } laplaceBaseViewKeys;
  //END_SPHINX_INCLUDE_04

protected:     // TODO consider making these private

  // These two classes are specific to the Laplace solver:
  string m_fieldName;      // User-defined name of the physical quantity we wish to solve for (such as "Temperature", etc.)
  TimeIntegrationOption m_timeIntegrationOption;      // Choice of transient treatment (SteadyState, ImplicitTransient)

};


/* REGISTERING NEW ENUMS:
   --------------------------
   In order to register an enumeration type with the Data Repository and have its value read from input,
   we must define stream insertion/extraction operators. This is a common task, so GEOSX provides
   a facility for automating it. Upon including ``common/EnumStrings.hpp``, we can call the following macro
   at the namespace scope (in this case, right after the ``LaplaceBase`` class definition is complete):
 */
//START_SPHINX_INCLUDE_05
ENUM_STRINGS( LaplaceBase::TimeIntegrationOption, "SteadyState", "ImplicitTransient", "ExplicitTransient" )
//END_SPHINX_INCLUDE_05

} /* namespace geosx */

#endif /* GEOSX_PHYSICSSOLVERS_SIMPLEPDE_LAPLACE_BASE_HPP_ */
