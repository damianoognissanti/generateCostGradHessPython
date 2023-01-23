#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Bruno_JExpBot2016 {

static constexpr std::array<sunindextype, 7> dxdotdw_colptrs_Bruno_JExpBot2016_ = {
    0, 3, 5, 8, 11, 13, 16
};

void dxdotdw_colptrs_Bruno_JExpBot2016(SUNMatrixWrapper &dxdotdw){
    dxdotdw.set_indexptrs(gsl::make_span(dxdotdw_colptrs_Bruno_JExpBot2016_));
}
} // namespace model_Bruno_JExpBot2016
} // namespace amici
