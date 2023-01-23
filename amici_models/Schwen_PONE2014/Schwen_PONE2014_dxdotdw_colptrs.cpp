#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Schwen_PONE2014 {

static constexpr std::array<sunindextype, 15> dxdotdw_colptrs_Schwen_PONE2014_ = {
    0, 3, 6, 8, 10, 13, 16, 18, 20, 22, 24, 27, 30, 31, 32
};

void dxdotdw_colptrs_Schwen_PONE2014(SUNMatrixWrapper &dxdotdw){
    dxdotdw.set_indexptrs(gsl::make_span(dxdotdw_colptrs_Schwen_PONE2014_));
}
} // namespace model_Schwen_PONE2014
} // namespace amici
