#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Schwen_PONE2014 {

static constexpr std::array<sunindextype, 28> dwdp_colptrs_Schwen_PONE2014_ = {
    0, 0, 0, 4, 6, 10, 12, 13, 14, 15, 16, 17, 18, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20
};

void dwdp_colptrs_Schwen_PONE2014(SUNMatrixWrapper &dwdp){
    dwdp.set_indexptrs(gsl::make_span(dwdp_colptrs_Schwen_PONE2014_));
}
} // namespace model_Schwen_PONE2014
} // namespace amici
