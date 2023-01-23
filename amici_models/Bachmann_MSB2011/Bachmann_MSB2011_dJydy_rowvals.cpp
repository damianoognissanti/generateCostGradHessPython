#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Bachmann_MSB2011 {

static constexpr std::array<std::array<sunindextype, 1>, 20> dJydy_rowvals_Bachmann_MSB2011_ = {{
    {0}, 
    {0}, 
    {0}, 
    {0}, 
    {0}, 
    {0}, 
    {0}, 
    {0}, 
    {0}, 
    {0}, 
    {0}, 
    {0}, 
    {0}, 
    {0}, 
    {0}, 
    {0}, 
    {0}, 
    {0}, 
    {0}, 
    {0}, 
}};

void dJydy_rowvals_Bachmann_MSB2011(SUNMatrixWrapper &dJydy, int index){
    dJydy.set_indexvals(gsl::make_span(dJydy_rowvals_Bachmann_MSB2011_[index]));
}
} // namespace model_Bachmann_MSB2011
} // namespace amici
