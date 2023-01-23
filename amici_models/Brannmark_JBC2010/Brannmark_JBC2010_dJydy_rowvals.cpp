#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Brannmark_JBC2010 {

static constexpr std::array<std::array<sunindextype, 1>, 3> dJydy_rowvals_Brannmark_JBC2010_ = {{
    {0}, 
    {0}, 
    {0}, 
}};

void dJydy_rowvals_Brannmark_JBC2010(SUNMatrixWrapper &dJydy, int index){
    dJydy.set_indexvals(gsl::make_span(dJydy_rowvals_Brannmark_JBC2010_[index]));
}
} // namespace model_Brannmark_JBC2010
} // namespace amici
