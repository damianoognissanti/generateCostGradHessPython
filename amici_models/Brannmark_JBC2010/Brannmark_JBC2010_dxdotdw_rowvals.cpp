#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Brannmark_JBC2010 {

static constexpr std::array<sunindextype, 22> dxdotdw_rowvals_Brannmark_JBC2010_ = {
    0, 1, 0, 1, 1, 2, 2, 3, 3, 4, 0, 2, 0, 4, 5, 6, 5, 6, 7, 8, 7, 8
};

void dxdotdw_rowvals_Brannmark_JBC2010(SUNMatrixWrapper &dxdotdw){
    dxdotdw.set_indexvals(gsl::make_span(dxdotdw_rowvals_Brannmark_JBC2010_));
}
} // namespace model_Brannmark_JBC2010
} // namespace amici
