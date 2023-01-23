#include "amici/model.h"
#include "wrapfunctions.h"
#include "Bachmann_MSB2011.h"

namespace amici {
namespace generic_model {

std::unique_ptr<amici::Model> getModel() {
    return std::unique_ptr<amici::Model>(
        new amici::model_Bachmann_MSB2011::Model_Bachmann_MSB2011());
}


} // namespace generic_model

} // namespace amici
