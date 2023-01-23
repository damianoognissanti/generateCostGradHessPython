#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Brannmark_JBC2010 {

static constexpr std::array<sunindextype, 14> dwdp_rowvals_Brannmark_JBC2010_ = {
    1, 1, 2, 3, 4, 5, 5, 6, 7, 8, 8, 10, 9, 11
};

void dwdp_rowvals_Brannmark_JBC2010(SUNMatrixWrapper &dwdp){
    dwdp.set_indexvals(gsl::make_span(dwdp_rowvals_Brannmark_JBC2010_));
}
} // namespace model_Brannmark_JBC2010
} // namespace amici
