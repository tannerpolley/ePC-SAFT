#include "ceres_backend.h"

namespace epcsaft::native::regression {

bool ceres_backend_compiled()
{
#ifdef EPCSAFT_HAS_CERES
    return true;
#else
    return false;
#endif
}

}  // namespace epcsaft::native::regression
