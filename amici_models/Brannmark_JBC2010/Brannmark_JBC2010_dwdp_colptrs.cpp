#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Brannmark_JBC2010 {

static constexpr std::array<sunindextype, 22> dwdp_colptrs_Brannmark_JBC2010_ = {
    0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 12, 13, 14, 14, 14, 14, 14, 14, 14
};

void dwdp_colptrs_Brannmark_JBC2010(SUNMatrixWrapper &dwdp){
    dwdp.set_indexptrs(gsl::make_span(dwdp_colptrs_Brannmark_JBC2010_));
}
} // namespace model_Brannmark_JBC2010
} // namespace amici
