#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Brannmark_JBC2010 {

static constexpr std::array<sunindextype, 13> dwdw_colptrs_Brannmark_JBC2010_ = {
    0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1
};

void dwdw_colptrs_Brannmark_JBC2010(SUNMatrixWrapper &dwdw){
    dwdw.set_indexptrs(gsl::make_span(dwdw_colptrs_Brannmark_JBC2010_));
}
} // namespace model_Brannmark_JBC2010
} // namespace amici
