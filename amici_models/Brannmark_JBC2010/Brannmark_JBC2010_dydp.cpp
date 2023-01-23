#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <array>
#include <algorithm>

#include "Brannmark_JBC2010_x.h"
#include "Brannmark_JBC2010_p.h"
#include "Brannmark_JBC2010_k.h"
#include "Brannmark_JBC2010_h.h"
#include "Brannmark_JBC2010_w.h"

namespace amici {
namespace model_Brannmark_JBC2010 {

void dydp_Brannmark_JBC2010(realtype *dydp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const int ip, const realtype *w, const realtype *tcl, const realtype *dtcldp){
    switch(ip) {
        case 15:
            dydp[0] = IRiP + IRp;
            break;
        case 16:
            dydp[1] = IRSiP;
            break;
        case 17:
            dydp[2] = IRSiP;
            break;
    }
}

} // namespace model_Brannmark_JBC2010
} // namespace amici
