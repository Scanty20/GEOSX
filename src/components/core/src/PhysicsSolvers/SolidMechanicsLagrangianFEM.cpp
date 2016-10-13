/*
 * NewtonianMechanics.cpp
 *
 *  Created on: Dec 4, 2014
 *      Author: rrsettgast
 */

#include "SolidMechanicsLagrangianFEM.hpp"

#include <vector>
#include <math.h>

// #include "RAJA/RAJA.hxx"
#include "dataRepository/ManagedGroup.hpp"
#include "common/DataTypes.hpp"


namespace geosx
{

namespace dataRepository
{
namespace keys
{
std::string const K = "K";
std::string const Ey = "Ey";
std::string const rho = "rho";
std::string const area = "area";
std::string const barLength = "barLength";
std::string const nElements = "nElements";
}
}

using namespace dataRepository;

SolidMechanics_LagrangianFEM::SolidMechanics_LagrangianFEM( const std::string& name,
                                                            ManagedGroup * const parent ) :
  SolverBase( name, parent )
{}



SolidMechanics_LagrangianFEM::~SolidMechanics_LagrangianFEM()
{
  // TODO Auto-generated destructor stub
}


void SolidMechanics_LagrangianFEM::FillDocumentationNode( dataRepository::ManagedGroup * const group )
{
  cxx_utilities::DocumentationNode * const docNode = this->getDocumentationNode();

  SolverBase::FillDocumentationNode( group );

  docNode->setName(this->CatalogName());    // If this method lived in Managed groups, this could be done automatically
  docNode->setSchemaType("Node");
  
  docNode->AllocateChildNode( "nElements",
                              keys::nElements,
                              -1,
                              "int32",
                              "int32",
                              "number of elements",
                              "number of elements",
                              "10",
                              "",
                              1,
                              0 );

  docNode->AllocateChildNode( "Ey",
                              keys::Ey,
                              -1,
                              "real64",
                              "real64",
                              "Elastic Young's Modulus",
                              "Elastic Young's Modulus",
                              "1.0e9",
                              "",
                              1,
                              0 );

  docNode->AllocateChildNode( "rho",
                              keys::rho,
                              -1,
                              "real64",
                              "real64",
                              "Initial Density",
                              "Initial Density",
                              "2600.0",
                              "",
                              1,
                              0 );

  docNode->AllocateChildNode( "area",
                              keys::area,
                              -1,
                              "real64",
                              "real64",
                              "cross section area",
                              "cross section area",
                              "1.0",
                              "",
                              1,
                              0 );

  docNode->AllocateChildNode( "barLength",
                              "barLength",
                              -1,
                              "real64",
                              "real64",
                              "reference length",
                              "reference length",
                              "1.0",
                              "",
                              1,
                              0 );
}

void SolidMechanics_LagrangianFEM::ReadXML( pugi::xml_node const & solverNode )
{
  SolverBase::ReadXML(solverNode);

  // cxx_utilities::DocumentationNode docVar;

  /*
  *(this->getData<int>(keys::nElements)) = solverNode.attribute("nElements").as_int(10);



  *(this->getData<real64>(keys::Ey)) = solverNode.attribute("Ey").as_double(10.0e9);

  *(this->getData<real64>(keys::rho)) = solverNode.attribute("rho").as_double(2650.0);


  *(this->getData<real64>(keys::area)) = solverNode.attribute("area").as_double(1.0);
  *(this->getData<real64>(keys::barLength)) = solverNode.attribute("barLength").as_double(1.0);
  */

  // Test auto-read:
  cxx_utilities::DocumentationNode * const docNode = this->getDocumentationNode();
  for( auto const & subDocNode : docNode.m_child )
  {
    std::string childType = subDocNode.getSchemaType();
    if (strcmp(childType, "real64") == 0)
    {
      real64 defVal = atof(subDocNode.getDefault());
      real64 xmlVal = solverNode.attribute(subDocNode.getStringKey()).as_double(defVal);
      *(this->getData<real64>(subDocNode.getStringKey())) = xmlVal;
    }
    else if (strcmp(childType, "int32") == 0)
    {
      int32 defVal = atoi(subDocNode.getDefault());
      int32 xmlVal = solverNode.attribute(subDocNode.getStringKey()).as_int(defVal);
      *(this->getData<int64>(subDocNode.getStringKey())) = xmlVal;
    }
    else if (strcmp(childType, "uint64") == 0)
    {
      uint32 defVal = atol(subDocNode.getDefault());
      uint32 xmlVal = solverNode.attribute(subDocNode.getStringKey()).as_uint(defVal);
      *(this->getData<uint64>(subDocNode.getStringKey())) = xmlVal;
    }
    else if (strcmp(childType, "string") == 0)
    {
      string defVal = subDocNode.getDefault();
      string xmlVal = solverNode.attribute(subDocNode.getStringKey()).as_double(defVal);
      *(this->getData<string>(subDocNode.getStringKey())) = xmlVal.empty() ? defVal : xmlVal;
    }
    else if (strcmp(childType, "real64_array") == 0)
    {
      string defVal = subDocNode.getDefault();
      std::vector<real64> xmlVal;
      solverNode.attribute(subDocNode.getStringKey()).load_double_array(xmlVal, defVal);
      // TODO: Figure out how to store these values (they may have different lengths)
    }
    else if (strcmp(childType, "int32_array") == 0)
    {
      string defVal = subDocNode.getDefault();
      std::vector<int32> xmlVal;
      solverNode.attribute(subDocNode.getStringKey()).load_int_array(xmlVal, defVal);
      // TODO: Figure out how to store these values (they may have different lengths)
    }
    else if (strcmp(childType, "uint32_array") == 0)
    {
      string defVal = subDocNode.getDefault();
      std::vector<uint32> xmlVal;
      solverNode.attribute(subDocNode.getStringKey()).load_uint_array(xmlVal, defVal);
      // TODO: Figure out how to store these values (they may have different lengths)
    }
    else if (strcmp(childType, "string_array") == 0)
    {
      string defVal = subDocNode.getDefault();
      std::vector<string> xmlVal;
      solverNode.attribute(subDocNode.getStringKey()).load_string_array(xmlVal, defVal);
      // TODO: Figure out how to store these values (they may have different lengths)
    }


    this->RegisterViewWrapper<subDocNode.getDataType()>(subDocNode.getStringKey());
  }


}


void SolidMechanics_LagrangianFEM::BuildDataStructure( ManagedGroup * const domain )
{
  SolverBase::BuildDataStructure( domain );

  ManagedGroup& nodes = domain->RegisterGroup<ManagedGroup>(keys::FEM_Nodes);
  ManagedGroup& elems = domain->RegisterGroup<ManagedGroup>(keys::FEM_Elements);

  nodes.RegisterViewWrapper<real64_array>(keys::TotalDisplacement);
  nodes.RegisterViewWrapper<real64_array>(keys::IncrementalDisplacement);
  nodes.RegisterViewWrapper<real64_array>(keys::Velocity);
  nodes.RegisterViewWrapper<real64_array>(keys::Acceleration);

  elems.RegisterViewWrapper<real64_array>(keys::Strain);
  elems.RegisterViewWrapper(keys::Force,  rtTypes::TypeIDs::real64_array_id );
  elems.RegisterViewWrapper<real64>(keys::Ey);
  elems.RegisterViewWrapper<real64_array>(keys::K);

  nodes.RegisterViewWrapper<real64_array>(keys::ReferencePosition);
  nodes.RegisterViewWrapper<real64_array>(keys::Mass);

  /*
  // Lagrange solver parameters
  this->RegisterViewWrapper<int>(keys::nElements);
  this->RegisterViewWrapper<real64>(keys::Ey);
  this->RegisterViewWrapper<real64>(keys::rho);
  this->RegisterViewWrapper<real64>(keys::area);
  this->RegisterViewWrapper<real64>(keys::barLength);
  */

  // Test auto-registration:
  RegisterDocumentationNodes();

}


void SolidMechanics_LagrangianFEM::Initialize( dataRepository::ManagedGroup& domain )
{
  ManagedGroup& nodes = domain.GetGroup<ManagedGroup >(keys::FEM_Nodes);
  ManagedGroup& elems = domain.GetGroup<ManagedGroup >(keys::FEM_Elements);

  int& nElements = *(this->getData<int>(keys::nElements));
  real64& Ey = *(this->getData<real64>(keys::Ey));
  real64& rho = *(this->getData<real64>(keys::rho));
  real64& area = *(this->getData<real64>(keys::area));
  real64& barLength = *(this->getData<real64>(keys::barLength));

  // HACK
  nodes.resize(nElements+1);
  elems.resize(nElements);

  ViewWrapper<real64_array>::rtype    X = nodes.getData<real64_array>(keys::ReferencePosition);
  ViewWrapper<real64_array>::rtype mass = nodes.getData<real64_array>(keys::Mass);
  ViewWrapper<real64_array>::rtype K = elems.getData<real64_array>(keys::K);


  std::cout<<"sound speed = "<<sqrt(Ey/rho)<<std::endl;
//  std::cout<<1.0/0.0<<std::endl;
  for( localIndex a=0 ; a<nodes.size() ; ++a )
  {
    X[a] = a * ( barLength/(nElements+1));
  }

  for( localIndex k=0 ; k<elems.size() ; ++k )
  {
    double dx = barLength / nElements;
    mass[k] += rho * area * dx / 2;
    mass[k+1] += rho * area * dx / 2;
    K[k] = Ey*area*dx;
  }

}

void SolidMechanics_LagrangianFEM::TimeStep( real64 const& time_n,
                                             real64 const& dt,
                                             const int cycleNumber,
                                             ManagedGroup& domain )
{
  TimeStepExplicit( time_n, dt, cycleNumber, domain );
}

void SolidMechanics_LagrangianFEM::TimeStepExplicit( real64 const& time_n,
                                                     real64 const& dt,
                                                     const int cycleNumber,
                                                     ManagedGroup& domain )
{
  ManagedGroup& nodes = domain.GetGroup<ManagedGroup>(keys::FEM_Nodes);
  ManagedGroup& elems = domain.GetGroup<ManagedGroup>(keys::FEM_Elements);

  localIndex const numNodes = nodes.size();
  localIndex const numElems = elems.size();

  ViewWrapper<real64_array>::rtype          X = nodes.getData<real64_array>(keys::ReferencePosition);
  ViewWrapper<real64_array>::rtype          u = nodes.getData<real64_array>(keys::TotalDisplacement);
  ViewWrapper<real64_array>::rtype       uhat = nodes.getData<real64_array>(keys::IncrementalDisplacement);
  ViewWrapper<real64_array>::rtype       vel  = nodes.getData<real64_array>(keys::Velocity);
  ViewWrapper<real64_array>::rtype       acc  = nodes.getData<real64_array>(keys::Acceleration);
  ViewWrapper<real64_array>::rtype_const mass = nodes.getWrapper<real64_array>(keys::Mass).data();

  ViewWrapper<real64_array>::rtype    Felem = elems.getData<real64_array>(keys::Force);
  ViewWrapper<real64_array>::rtype   Strain = elems.getData<real64_array>(keys::Strain);
  ViewWrapper<real64_array>::rtype_const  K = elems.getData<real64_array>(keys::K);


//  ViewWrapper<real64_array>::rtype          X2 = nodes.GetData(keys::ReferencePosition);


  Integration::OnePoint( acc, vel, dt/2, numNodes );
  vel[0] = 1.0;
  Integration::OnePoint( vel, uhat, u, dt, numNodes );

  for( localIndex a=0 ; a<numElems ; ++a )
  {
    acc[a] = 0.0;
  }

  for( localIndex k=0 ; k<numElems ; ++k )
  {
    Strain[k] = ( u[k+1] - u[k] ) / ( X[k+1] - X[k] );
    Felem[k] = K[k] * Strain[k];
    acc[k]   += Felem[k];
    acc[k+1] -= Felem[k];
  }

  for( localIndex a=0 ; a<numNodes ; ++a )
  {
    acc[a] /= mass[a];
  }
  Integration::OnePoint( acc, vel, dt/2, numNodes );
  vel[0] = 1.0;


  printf(" %6.5f : ", time_n + dt );
  for( localIndex a=0 ; a<numNodes ; ++a )
  {
    printf(" %4.2f ",vel[a] );
//    std::cout<<vel[a]<<std::endl;
  }
  printf("\n" );

  (void) time_n;
  (void) cycleNumber;
}


REGISTER_CATALOG_ENTRY( SolverBase, SolidMechanics_LagrangianFEM, std::string const &, ManagedGroup * const )
} /* namespace ANST */
