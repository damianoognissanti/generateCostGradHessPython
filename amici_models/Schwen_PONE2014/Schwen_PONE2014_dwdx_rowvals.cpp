#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Schwen_PONE2014 {

static constexpr std::array<sunindextype, 20> dwdx_rowvals_Schwen_PONE2014_ = {
    0, 1, 2, 12, 13, 0, 12, 1, 13, 4, 6, 12, 5, 7, 13, 8, 10, 9, 11, 3
};

void dwdx_rowvals_Schwen_PONE2014(SUNMatrixWrapper &dwdx){
    dwdx.set_indexvals(gsl::make_span(dwdx_rowvals_Schwen_PONE2014_));
}
} // namespace model_Schwen_PONE2014
} // namespace amici
