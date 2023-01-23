#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Schwen_PONE2014 {

static constexpr std::array<sunindextype, 12> dwdx_colptrs_Schwen_PONE2014_ = {
    0, 5, 7, 9, 12, 15, 17, 19, 19, 19, 19, 20
};

void dwdx_colptrs_Schwen_PONE2014(SUNMatrixWrapper &dwdx){
    dwdx.set_indexptrs(gsl::make_span(dwdx_colptrs_Schwen_PONE2014_));
}
} // namespace model_Schwen_PONE2014
} // namespace amici
