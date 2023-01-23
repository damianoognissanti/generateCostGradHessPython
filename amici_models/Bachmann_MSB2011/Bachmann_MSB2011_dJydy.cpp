#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <algorithm>

#include "Bachmann_MSB2011_p.h"
#include "Bachmann_MSB2011_k.h"
#include "Bachmann_MSB2011_y.h"
#include "Bachmann_MSB2011_sigmay.h"
#include "Bachmann_MSB2011_my.h"
#include "Bachmann_MSB2011_dJydy.h"

namespace amici {
namespace model_Bachmann_MSB2011 {

void dJydy_Bachmann_MSB2011(realtype *dJydy, const int iy, const realtype *p, const realtype *k, const realtype *y, const realtype *sigmay, const realtype *my){
    switch(iy) {
        case 0:
            dJydy[0] = (-1.0*std::log(mobservable_CISRNA_foldA)/M_LN10 + 1.0*std::log(observable_CISRNA_foldA)/M_LN10)/(observable_CISRNA_foldA*std::pow(sigma_observable_CISRNA_foldA, 2)*M_LN10);
            break;
        case 1:
            dJydy[0] = (-1.0*std::log(mobservable_CISRNA_foldB)/M_LN10 + 1.0*std::log(observable_CISRNA_foldB)/M_LN10)/(observable_CISRNA_foldB*std::pow(sigma_observable_CISRNA_foldB, 2)*M_LN10);
            break;
        case 2:
            dJydy[0] = (-1.0*std::log(mobservable_CISRNA_foldC)/M_LN10 + 1.0*std::log(observable_CISRNA_foldC)/M_LN10)/(observable_CISRNA_foldC*std::pow(sigma_observable_CISRNA_foldC, 2)*M_LN10);
            break;
        case 3:
            dJydy[0] = (-1.0*std::log(mobservable_CIS_abs)/M_LN10 + 1.0*std::log(observable_CIS_abs)/M_LN10)/(observable_CIS_abs*std::pow(sigma_observable_CIS_abs, 2)*M_LN10);
            break;
        case 4:
            dJydy[0] = (-1.0*std::log(mobservable_CIS_au)/M_LN10 + 1.0*std::log(observable_CIS_au)/M_LN10)/(observable_CIS_au*std::pow(sigma_observable_CIS_au, 2)*M_LN10);
            break;
        case 5:
            dJydy[0] = (-1.0*std::log(mobservable_CIS_au1)/M_LN10 + 1.0*std::log(observable_CIS_au1)/M_LN10)/(observable_CIS_au1*std::pow(sigma_observable_CIS_au1, 2)*M_LN10);
            break;
        case 6:
            dJydy[0] = (-1.0*std::log(mobservable_CIS_au2)/M_LN10 + 1.0*std::log(observable_CIS_au2)/M_LN10)/(observable_CIS_au2*std::pow(sigma_observable_CIS_au2, 2)*M_LN10);
            break;
        case 7:
            dJydy[0] = (-1.0*std::log(mobservable_SHP1_abs)/M_LN10 + 1.0*std::log(observable_SHP1_abs)/M_LN10)/(observable_SHP1_abs*std::pow(sigma_observable_SHP1_abs, 2)*M_LN10);
            break;
        case 8:
            dJydy[0] = (-1.0*std::log(mobservable_SOCS3RNA_foldA)/M_LN10 + 1.0*std::log(observable_SOCS3RNA_foldA)/M_LN10)/(observable_SOCS3RNA_foldA*std::pow(sigma_observable_SOCS3RNA_foldA, 2)*M_LN10);
            break;
        case 9:
            dJydy[0] = (-1.0*std::log(mobservable_SOCS3RNA_foldB)/M_LN10 + 1.0*std::log(observable_SOCS3RNA_foldB)/M_LN10)/(observable_SOCS3RNA_foldB*std::pow(sigma_observable_SOCS3RNA_foldB, 2)*M_LN10);
            break;
        case 10:
            dJydy[0] = (-1.0*std::log(mobservable_SOCS3RNA_foldC)/M_LN10 + 1.0*std::log(observable_SOCS3RNA_foldC)/M_LN10)/(observable_SOCS3RNA_foldC*std::pow(sigma_observable_SOCS3RNA_foldC, 2)*M_LN10);
            break;
        case 11:
            dJydy[0] = (-1.0*std::log(mobservable_SOCS3_abs)/M_LN10 + 1.0*std::log(observable_SOCS3_abs)/M_LN10)/(observable_SOCS3_abs*std::pow(sigma_observable_SOCS3_abs, 2)*M_LN10);
            break;
        case 12:
            dJydy[0] = (-1.0*std::log(mobservable_SOCS3_au)/M_LN10 + 1.0*std::log(observable_SOCS3_au)/M_LN10)/(observable_SOCS3_au*std::pow(sigma_observable_SOCS3_au, 2)*M_LN10);
            break;
        case 13:
            dJydy[0] = (-1.0*std::log(mobservable_STAT5_abs)/M_LN10 + 1.0*std::log(observable_STAT5_abs)/M_LN10)/(observable_STAT5_abs*std::pow(sigma_observable_STAT5_abs, 2)*M_LN10);
            break;
        case 14:
            dJydy[0] = (-1.0*std::log(mobservable_pEpoR_au)/M_LN10 + 1.0*std::log(observable_pEpoR_au)/M_LN10)/(observable_pEpoR_au*std::pow(sigma_observable_pEpoR_au, 2)*M_LN10);
            break;
        case 15:
            dJydy[0] = (-1.0*std::log(mobservable_pJAK2_au)/M_LN10 + 1.0*std::log(observable_pJAK2_au)/M_LN10)/(observable_pJAK2_au*std::pow(sigma_observable_pJAK2_au, 2)*M_LN10);
            break;
        case 16:
            dJydy[0] = (-1.0*mobservable_pSTAT5B_rel + 1.0*observable_pSTAT5B_rel)/std::pow(sigma_observable_pSTAT5B_rel, 2);
            break;
        case 17:
            dJydy[0] = (-1.0*std::log(mobservable_pSTAT5_au)/M_LN10 + 1.0*std::log(observable_pSTAT5_au)/M_LN10)/(observable_pSTAT5_au*std::pow(sigma_observable_pSTAT5_au, 2)*M_LN10);
            break;
        case 18:
            dJydy[0] = (-1.0*std::log(mobservable_tSHP1_au)/M_LN10 + 1.0*std::log(observable_tSHP1_au)/M_LN10)/(observable_tSHP1_au*std::pow(sigma_observable_tSHP1_au, 2)*M_LN10);
            break;
        case 19:
            dJydy[0] = (-1.0*std::log(mobservable_tSTAT5_au)/M_LN10 + 1.0*std::log(observable_tSTAT5_au)/M_LN10)/(observable_tSTAT5_au*std::pow(sigma_observable_tSTAT5_au, 2)*M_LN10);
            break;
    }
}

} // namespace model_Bachmann_MSB2011
} // namespace amici
