#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Bachmann_MSB2011 {

static constexpr std::array<sunindextype, 60> dwdx_rowvals_Bachmann_MSB2011_ = {
    0, 1, 2, 3, 10, 12, 4, 6, 9, 10, 12, 13, 5, 7, 10, 12, 8, 9, 10, 12, 13, 3, 4, 9, 10, 1, 6, 7, 8, 11, 12, 13, 14, 15, 16, 26, 17, 18, 19, 20, 21, 22, 23, 13, 24, 27, 28, 29, 30, 31, 32, 33, 0, 2, 3, 4, 5, 12, 13, 34
};

void dwdx_rowvals_Bachmann_MSB2011(SUNMatrixWrapper &dwdx){
    dwdx.set_indexvals(gsl::make_span(dwdx_rowvals_Bachmann_MSB2011_));
}
} // namespace model_Bachmann_MSB2011
} // namespace amici
