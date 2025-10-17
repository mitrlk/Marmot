#include "Marmot/DufourModel.h"
#include "Marmot/MarmotMaterialRegistrationHelper.h"

namespace Marmot::Materials {

  namespace Registration {

    constexpr int DufourModelCode = 119300112;

    using namespace MarmotLibrary;

    const static bool DufourModelRegistered = MarmotMaterialFactory::
      registerMaterial( DufourModelCode,
                        "DUFOURMODEL",
                        makeDefaultMarmotMaterialFactoryFunction< class DufourModel >() );
  } // namespace Registration
} // namespace Marmot::Materials
