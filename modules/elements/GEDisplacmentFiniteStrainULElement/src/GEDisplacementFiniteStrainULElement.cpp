#include "Marmot/GEDisplacementFiniteStrainULElement.h"
#include "Marmot/Marmot.h"
#include "Marmot/MarmotFiniteElement.h"

namespace Marmot::Elements::Registration {

#define CONCAT( a, b ) a##b
  enum GEDisplacementFiniteStrainULElementCode {

    /* TAG explanation
     * XXXXXX
     * ||||||_    6: formulation
     * |||||__    5: formulation
     * ||||___    4: type of element
     * |||____    3: active fields
     * ||_____    2: number of nodes
     * |______    1: number of nodes
     *
     * formulation:     01 Update Lagrangian
     *
     * active fields:   1: displacement
     *
     *
     * type of element: 2: 2D full integration, plane stress
     *                  3: 3D full integration,
     *                  5: 2D red. integration, plane stress
     *                  6: 3D red. integration
     *                  7: 2D full integration, plane strain
     *                  8: 2D red. integration, plane strain
     *                  9: 2D red. integration, axisymmetric
     *                  1: 2D full integration, axisymmetric
     * */

    CPE8RUL  = CONCAT( 1193, 83801 ),
    C3D8UL   = CONCAT( 1193, 83301 ),
    C3D20UL  = CONCAT( 1193, 203301 ),
    C3D20RUL = CONCAT( 1193, 203601 ),
    CX8RUL   = CONCAT( 1193, 83901 ),
    CX8UL    = CONCAT( 1193, 83101 ),

  };
#undef CONCAT

  template < class T,
             Marmot::FiniteElement::Quadrature::IntegrationTypes integrationType,
             typename T::SectionType                             sectionType >
  MarmotLibrary::MarmotElementFactory::elementFactoryFunction makeFactoryFunction()
  {
    return []( int elementID ) -> MarmotElement* { return new T( elementID, integrationType, sectionType ); };
  }

  using namespace MarmotLibrary;
  using namespace Marmot::FiniteElement::Quadrature;

  const static bool CX8RUL_isRegistered = MarmotElementFactory::
    registerElement( "CX8RUL",
                     GEDisplacementFiniteStrainULElementCode::CX8RUL,
                     makeFactoryFunction< AxiSymmetricGEDisplacementFiniteStrainULElement< 8 >,
                                          ReducedIntegration,
                                          AxiSymmetricGEDisplacementFiniteStrainULElement< 8 >::PlaneStrain >() );

  const static bool CX8UL_isRegistered = MarmotElementFactory::
    registerElement( "CX8UL",
                     GEDisplacementFiniteStrainULElementCode::CX8UL,
                     makeFactoryFunction< AxiSymmetricGEDisplacementFiniteStrainULElement< 8 >,
                                          FullIntegration,
                                          AxiSymmetricGEDisplacementFiniteStrainULElement< 8 >::PlaneStrain >() );

  const static bool CPE8RGradientEnhancedMicropolar_isRegistered = MarmotElementFactory::
    registerElement( "CPE8RUL",
                     GEDisplacementFiniteStrainULElementCode::CPE8RUL,
                     makeFactoryFunction< GEDisplacementFiniteStrainULElement< 2, 8 >,
                                          ReducedIntegration,
                                          GEDisplacementFiniteStrainULElement< 2, 8 >::PlaneStrain >() );

  const static bool C3D8UL_isRegistered = MarmotLibrary::MarmotElementFactory::
    registerElement( "C3D8UL", GEDisplacementFiniteStrainULElementCode::C3D8UL, []( int elementID ) -> MarmotElement* {
      return new GEDisplacementFiniteStrainULElement<
        3,
        8 >( elementID,
             Marmot::FiniteElement::Quadrature::IntegrationTypes::FullIntegration,
             GEDisplacementFiniteStrainULElement< 3, 8 >::SectionType::Solid );
    } );

  const static bool C3D20RUL_isRegistered = MarmotLibrary::MarmotElementFactory::
    registerElement( "C3D20RUL",
                     GEDisplacementFiniteStrainULElementCode::C3D20RUL,
                     []( int elementID ) -> MarmotElement* {
                       return new GEDisplacementFiniteStrainULElement<
                         3,
                         20 >( elementID,
                               Marmot::FiniteElement::Quadrature::IntegrationTypes::ReducedIntegration,
                               GEDisplacementFiniteStrainULElement< 3, 20 >::SectionType::Solid );
                     } );

  const static bool C3D20UL_isRegistered = MarmotLibrary::MarmotElementFactory::
    registerElement( "C3D20UL",
                     GEDisplacementFiniteStrainULElementCode::C3D20UL,
                     []( int elementID ) -> MarmotElement* {
                       return new GEDisplacementFiniteStrainULElement<
                         3,
                         20 >( elementID,
                               Marmot::FiniteElement::Quadrature::IntegrationTypes::FullIntegration,
                               GEDisplacementFiniteStrainULElement< 3, 20 >::SectionType::Solid );
                     } );

} // namespace Marmot::Elements::Registration
