
#include <satutils/constants.hpp>

namespace satutils {

std::size_t CodeLength(CodeId code)
{
  switch(code) {
    case GPSCA:
    case GPSL1C:
      return GPS_CA_CODE_LENGTH;
    case GPSL2CM:
      return GPS_L2CM_CODE_LENGTH;
    case GPSL2CL:
      return GPS_L2CL_CODE_LENGTH;
    case GPSL5I:
    case GPSL5Q:
    case GalileoE1OS:
    case GalileoE5A:
    case GalileoE5B:
    case GalileoE6CS:
    default:
      assert(false);
      return 0;
  }
}

} // namespace satutils
