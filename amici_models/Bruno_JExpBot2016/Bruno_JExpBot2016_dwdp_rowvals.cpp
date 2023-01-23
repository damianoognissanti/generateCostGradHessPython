#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Bruno_JExpBot2016 {

static constexpr std::array<sunindextype, 12> dwdp_rowvals_Bruno_JExpBot2016_ = {
    5, 0, 1, 2, 3, 4, 5, 0, 1, 2, 3, 4
};

void dwdp_rowvals_Bruno_JExpBot2016(SUNMatrixWrapper &dwdp){
    dwdp.set_indexvals(gsl::make_span(dwdp_rowvals_Bruno_JExpBot2016_));
}
} // namespace model_Bruno_JExpBot2016
} // namespace amici
