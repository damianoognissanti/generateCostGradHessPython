#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Brannmark_JBC2010 {

static constexpr std::array<sunindextype, 1> dwdw_rowvals_Brannmark_JBC2010_ = {
    1
};

void dwdw_rowvals_Brannmark_JBC2010(SUNMatrixWrapper &dwdw){
    dwdw.set_indexvals(gsl::make_span(dwdw_rowvals_Brannmark_JBC2010_));
}
} // namespace model_Brannmark_JBC2010
} // namespace amici
