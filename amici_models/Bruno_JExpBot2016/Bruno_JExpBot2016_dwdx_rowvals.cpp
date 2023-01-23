#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Bruno_JExpBot2016 {

static constexpr std::array<sunindextype, 6> dwdx_rowvals_Bruno_JExpBot2016_ = {
    0, 2, 3, 1, 4, 5
};

void dwdx_rowvals_Bruno_JExpBot2016(SUNMatrixWrapper &dwdx){
    dwdx.set_indexvals(gsl::make_span(dwdx_rowvals_Bruno_JExpBot2016_));
}
} // namespace model_Bruno_JExpBot2016
} // namespace amici
