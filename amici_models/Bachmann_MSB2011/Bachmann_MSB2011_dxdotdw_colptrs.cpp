#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Bachmann_MSB2011 {

static constexpr std::array<sunindextype, 37> dxdotdw_colptrs_Bachmann_MSB2011_ = {
    0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 19, 21, 23, 25, 27, 29, 31, 32, 34, 36, 38, 40, 42, 43, 44, 45, 46, 47, 49, 51, 53, 55, 57, 58, 59, 60, 61
};

void dxdotdw_colptrs_Bachmann_MSB2011(SUNMatrixWrapper &dxdotdw){
    dxdotdw.set_indexptrs(gsl::make_span(dxdotdw_colptrs_Bachmann_MSB2011_));
}
} // namespace model_Bachmann_MSB2011
} // namespace amici
