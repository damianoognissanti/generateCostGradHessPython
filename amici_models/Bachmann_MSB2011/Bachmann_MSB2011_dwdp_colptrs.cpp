#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Bachmann_MSB2011 {

static constexpr std::array<sunindextype, 70> dwdp_colptrs_Bachmann_MSB2011_ = {
    0, 3, 4, 9, 11, 14, 15, 19, 21, 22, 23, 27, 28, 29, 29, 30, 39, 46, 51, 53, 56, 57, 58, 59, 60, 64, 68, 70, 70, 70, 70, 70, 70, 70, 70, 70, 70, 70, 70, 70, 70, 70, 70, 70, 70, 70, 70, 70, 70, 70, 70, 70, 70, 70, 70, 70, 70, 70, 70, 70, 70, 70, 70, 70, 70, 70, 70, 70, 70, 70
};

void dwdp_colptrs_Bachmann_MSB2011(SUNMatrixWrapper &dwdp){
    dwdp.set_indexptrs(gsl::make_span(dwdp_colptrs_Bachmann_MSB2011_));
}
} // namespace model_Bachmann_MSB2011
} // namespace amici
