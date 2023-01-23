#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Bruno_JExpBot2016 {

static constexpr std::array<std::array<sunindextype, 1>, 6> dJydy_rowvals_Bruno_JExpBot2016_ = {{
    {0}, 
    {0}, 
    {0}, 
    {0}, 
    {0}, 
    {0}, 
}};

void dJydy_rowvals_Bruno_JExpBot2016(SUNMatrixWrapper &dJydy, int index){
    dJydy.set_indexvals(gsl::make_span(dJydy_rowvals_Bruno_JExpBot2016_[index]));
}
} // namespace model_Bruno_JExpBot2016
} // namespace amici
