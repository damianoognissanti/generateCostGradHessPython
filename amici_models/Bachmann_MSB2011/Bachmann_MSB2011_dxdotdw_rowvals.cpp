#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Bachmann_MSB2011 {

static constexpr std::array<sunindextype, 61> dxdotdw_rowvals_Bachmann_MSB2011_ = {
    0, 1, 0, 1, 1, 2, 1, 3, 2, 4, 3, 4, 0, 2, 0, 3, 0, 4, 5, 6, 7, 6, 7, 8, 9, 8, 9, 9, 10, 8, 10, 11, 11, 12, 12, 13, 13, 14, 14, 15, 15, 16, 16, 17, 17, 17, 18, 18, 19, 19, 20, 20, 21, 21, 22, 22, 23, 23, 24, 24, 24
};

void dxdotdw_rowvals_Bachmann_MSB2011(SUNMatrixWrapper &dxdotdw){
    dxdotdw.set_indexvals(gsl::make_span(dxdotdw_rowvals_Bachmann_MSB2011_));
}
} // namespace model_Bachmann_MSB2011
} // namespace amici
