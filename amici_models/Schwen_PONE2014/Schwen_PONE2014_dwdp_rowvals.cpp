#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Schwen_PONE2014 {

static constexpr std::array<sunindextype, 20> dwdp_rowvals_Schwen_PONE2014_ = {
    0, 1, 12, 13, 1, 13, 4, 5, 12, 13, 5, 13, 6, 7, 3, 2, 8, 9, 10, 11
};

void dwdp_rowvals_Schwen_PONE2014(SUNMatrixWrapper &dwdp){
    dwdp.set_indexvals(gsl::make_span(dwdp_rowvals_Schwen_PONE2014_));
}
} // namespace model_Schwen_PONE2014
} // namespace amici
