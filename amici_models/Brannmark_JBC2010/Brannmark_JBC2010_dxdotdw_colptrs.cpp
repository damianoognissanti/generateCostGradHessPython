#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Brannmark_JBC2010 {

static constexpr std::array<sunindextype, 13> dxdotdw_colptrs_Brannmark_JBC2010_ = {
    0, 0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22
};

void dxdotdw_colptrs_Brannmark_JBC2010(SUNMatrixWrapper &dxdotdw){
    dxdotdw.set_indexptrs(gsl::make_span(dxdotdw_colptrs_Brannmark_JBC2010_));
}
} // namespace model_Brannmark_JBC2010
} // namespace amici
