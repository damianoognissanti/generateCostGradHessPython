#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Schwen_PONE2014 {

static constexpr std::array<std::array<sunindextype, 1>, 4> dJydy_rowvals_Schwen_PONE2014_ = {{
    {0}, 
    {0}, 
    {0}, 
    {0}, 
}};

void dJydy_rowvals_Schwen_PONE2014(SUNMatrixWrapper &dJydy, int index){
    dJydy.set_indexvals(gsl::make_span(dJydy_rowvals_Schwen_PONE2014_[index]));
}
} // namespace model_Schwen_PONE2014
} // namespace amici
