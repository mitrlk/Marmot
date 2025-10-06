#include "Marmot/DufourModel.h"
#include "Marmot/MarmotMaterialRegistrationHelper.h"

namespace Marmot::Materials {

  namespace Registration {

    constexpr int DufourModelCode = 11930000 + 20;

    using namespace MarmotLibrary;

    const static bool DufourModelRegistered = MarmotMaterialFactory::
      registerMaterial( DufourModelCode,
                        "DufourModel",
                        makeDefaultMarmotMaterialFactoryFunction< class DufourModel >() );
  } // namespace Registration
} // namespace Marmot::Materials
