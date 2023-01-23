#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Bachmann_MSB2011 {

static constexpr std::array<sunindextype, 70> dwdp_rowvals_Bachmann_MSB2011_ = {
    13, 23, 25, 13, 17, 18, 19, 20, 21, 16, 22, 23, 24, 25, 25, 2, 3, 4, 5, 3, 4, 9, 0, 1, 6, 7, 8, 10, 11, 35, 0, 2, 3, 4, 5, 12, 13, 33, 35, 0, 2, 3, 4, 5, 12, 13, 27, 28, 29, 30, 31, 26, 32, 33, 34, 35, 13, 12, 15, 14, 9, 10, 12, 13, 1, 6, 7, 8, 16, 26
};

void dwdp_rowvals_Bachmann_MSB2011(SUNMatrixWrapper &dwdp){
    dwdp.set_indexvals(gsl::make_span(dwdp_rowvals_Bachmann_MSB2011_));
}
} // namespace model_Bachmann_MSB2011
} // namespace amici
