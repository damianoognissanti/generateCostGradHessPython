#include "amici/model.h"
#include "wrapfunctions.h"
#include "Brannmark_JBC2010.h"

namespace amici {
namespace generic_model {

std::unique_ptr<amici::Model> getModel() {
    return std::unique_ptr<amici::Model>(
        new amici::model_Brannmark_JBC2010::Model_Brannmark_JBC2010());
}


} // namespace generic_model

} // namespace amici
