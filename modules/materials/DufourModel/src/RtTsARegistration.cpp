#include "Marmot/MarmotMaterialRegistrationHelper.h"
#include "Marmot/RtTsA.h"

namespace Marmot::Materials {

  namespace Registration {

    constexpr int RtTsACode = 11930000 + 20;

    using namespace MarmotLibrary;

    const static bool RtTsARegistered = MarmotMaterialFactory::
      registerMaterial( RtTsACode, "RtTsA", makeDefaultMarmotMaterialFactoryFunction< class RtTsA >() );
  } // namespace Registration
} // namespace Marmot::Materials
