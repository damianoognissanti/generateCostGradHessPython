#include "amici/model.h"
#include "wrapfunctions.h"
#include "Bruno_JExpBot2016.h"

namespace amici {
namespace generic_model {

std::unique_ptr<amici::Model> getModel() {
    return std::unique_ptr<amici::Model>(
        new amici::model_Bruno_JExpBot2016::Model_Bruno_JExpBot2016());
}


} // namespace generic_model

} // namespace amici
