#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Bruno_JExpBot2016 {

static constexpr std::array<sunindextype, 16> dxdotdw_rowvals_Bruno_JExpBot2016_ = {
    0, 2, 3, 2, 3, 1, 2, 5, 1, 3, 4, 4, 5, 4, 5, 6
};

void dxdotdw_rowvals_Bruno_JExpBot2016(SUNMatrixWrapper &dxdotdw){
    dxdotdw.set_indexvals(gsl::make_span(dxdotdw_rowvals_Bruno_JExpBot2016_));
}
} // namespace model_Bruno_JExpBot2016
} // namespace amici
