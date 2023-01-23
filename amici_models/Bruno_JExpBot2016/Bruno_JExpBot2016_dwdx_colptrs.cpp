#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Bruno_JExpBot2016 {

static constexpr std::array<sunindextype, 8> dwdx_colptrs_Bruno_JExpBot2016_ = {
    0, 1, 3, 4, 4, 5, 5, 6
};

void dwdx_colptrs_Bruno_JExpBot2016(SUNMatrixWrapper &dwdx){
    dwdx.set_indexptrs(gsl::make_span(dwdx_colptrs_Bruno_JExpBot2016_));
}
} // namespace model_Bruno_JExpBot2016
} // namespace amici
