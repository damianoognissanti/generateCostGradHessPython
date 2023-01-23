#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <array>
#include <algorithm>

#include "Schwen_PONE2014_x.h"
#include "Schwen_PONE2014_p.h"
#include "Schwen_PONE2014_k.h"
#include "Schwen_PONE2014_w.h"
#include "Schwen_PONE2014_dwdx.h"

namespace amici {
namespace model_Schwen_PONE2014 {

void dwdx_Schwen_PONE2014(realtype *dwdx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl){
    dflux_v1_ka1_dIns = 1.0*Rec1*ka1;  // dwdx[0]
    dflux_v2_ka2fold_dIns = 1.0*Rec2*ka1*ka2fold;  // dwdx[1]
    dflux_v3_v_2_dIns = 1.0*kon_unspec;  // dwdx[2]
    dflux_v13_v_12_dIns = 1.0*Rec1*ka1;  // dwdx[3]
    dflux_v14_v_13_dIns = 1.0*Rec2*ka1*ka2fold;  // dwdx[4]
    dflux_v1_ka1_dRec1 = 1.0*Ins*ka1;  // dwdx[5]
    dflux_v13_v_12_dRec1 = 1.0*Ins*ka1;  // dwdx[6]
    dflux_v2_ka2fold_dRec2 = 1.0*Ins*ka1*ka2fold;  // dwdx[7]
    dflux_v14_v_13_dRec2 = 1.0*Ins*ka1*ka2fold;  // dwdx[8]
    dflux_v5_kd1_dIR1 = 1.0*kd1;  // dwdx[9]
    dflux_v7_v_6_dIR1 = 1.0*kin;  // dwdx[10]
    dflux_v13_v_12_dIR1 = -1.0*kd1;  // dwdx[11]
    dflux_v6_kd2fold_dIR2 = 1.0*kd1*kd2fold;  // dwdx[12]
    dflux_v8_v_7_dIR2 = 1.0*kin2;  // dwdx[13]
    dflux_v14_v_13_dIR2 = -1.0*kd1*kd2fold;  // dwdx[14]
    dflux_v9_v_8_dIR1in = 1.0*kout;  // dwdx[15]
    dflux_v11_v_10_dIR1in = 1.0*kout_frag;  // dwdx[16]
    dflux_v10_v_9_dIR2in = 1.0*kout2;  // dwdx[17]
    dflux_v12_v_11_dIR2in = 1.0*kout_frag;  // dwdx[18]
    dflux_v4_v_3_dBoundUnspec = 1.0*koff_unspec;  // dwdx[19]
}

} // namespace model_Schwen_PONE2014
} // namespace amici
