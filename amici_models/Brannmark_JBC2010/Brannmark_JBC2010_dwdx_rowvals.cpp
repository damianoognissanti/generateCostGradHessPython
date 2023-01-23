#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Brannmark_JBC2010 {

static constexpr std::array<sunindextype, 15> dwdx_rowvals_Brannmark_JBC2010_ = {
    1, 2, 3, 4, 6, 8, 5, 8, 7, 8, 9, 10, 10, 5, 11
};

void dwdx_rowvals_Brannmark_JBC2010(SUNMatrixWrapper &dwdx){
    dwdx.set_indexvals(gsl::make_span(dwdx_rowvals_Brannmark_JBC2010_));
}
} // namespace model_Brannmark_JBC2010
} // namespace amici
