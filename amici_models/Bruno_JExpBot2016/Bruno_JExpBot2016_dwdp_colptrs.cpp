#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Bruno_JExpBot2016 {

static constexpr std::array<sunindextype, 24> dwdp_colptrs_Bruno_JExpBot2016_ = {
    0, 0, 0, 0, 0, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 12, 12, 12, 12, 12, 12
};

void dwdp_colptrs_Bruno_JExpBot2016(SUNMatrixWrapper &dwdp){
    dwdp.set_indexptrs(gsl::make_span(dwdp_colptrs_Bruno_JExpBot2016_));
}
} // namespace model_Bruno_JExpBot2016
} // namespace amici
