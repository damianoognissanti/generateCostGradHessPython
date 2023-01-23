#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Bachmann_MSB2011 {

static constexpr std::array<sunindextype, 26> dwdx_colptrs_Bachmann_MSB2011_ = {
    0, 1, 6, 12, 16, 21, 24, 25, 30, 32, 33, 36, 37, 38, 39, 40, 41, 43, 45, 46, 47, 48, 49, 50, 52, 60
};

void dwdx_colptrs_Bachmann_MSB2011(SUNMatrixWrapper &dwdx){
    dwdx.set_indexptrs(gsl::make_span(dwdx_colptrs_Bachmann_MSB2011_));
}
} // namespace model_Bachmann_MSB2011
} // namespace amici
