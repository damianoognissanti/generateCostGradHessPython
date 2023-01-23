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
#include "Brannmark_JBC2010_dwdw.h"

namespace amici {
namespace model_Brannmark_JBC2010 {

void dwdw_Brannmark_JBC2010(realtype *dwdw, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl){
    dflux_v1_v_0_dinsulin = 1.0*IR*k1a;  // dwdw[0]
}

} // namespace model_Brannmark_JBC2010
} // namespace amici
