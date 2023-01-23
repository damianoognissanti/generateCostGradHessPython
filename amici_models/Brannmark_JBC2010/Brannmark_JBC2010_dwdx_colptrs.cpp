#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Brannmark_JBC2010 {

static constexpr std::array<sunindextype, 10> dwdx_colptrs_Brannmark_JBC2010_ = {
    0, 1, 3, 6, 8, 9, 10, 12, 13, 15
};

void dwdx_colptrs_Brannmark_JBC2010(SUNMatrixWrapper &dwdx){
    dwdx.set_indexptrs(gsl::make_span(dwdx_colptrs_Brannmark_JBC2010_));
}
} // namespace model_Brannmark_JBC2010
} // namespace amici
