#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Schwen_PONE2014 {

static constexpr std::array<sunindextype, 32> dxdotdw_rowvals_Schwen_PONE2014_ = {
    0, 1, 3, 0, 2, 4, 0, 10, 0, 10, 0, 1, 3, 0, 2, 4, 3, 5, 4, 6, 3, 5, 4, 6, 1, 5, 9, 2, 6, 9, 7, 8
};

void dxdotdw_rowvals_Schwen_PONE2014(SUNMatrixWrapper &dxdotdw){
    dxdotdw.set_indexvals(gsl::make_span(dxdotdw_rowvals_Schwen_PONE2014_));
}
} // namespace model_Schwen_PONE2014
} // namespace amici
