#include "Brannmark_JBC2010.h"
#include <array>

namespace amici {

namespace model_Brannmark_JBC2010 {

std::array<const char*, 21> parameterNames = {
    "k1a", // p[0]
"k1aBasic", // p[1]
"k1b", // p[2]
"k1c", // p[3]
"k1d", // p[4]
"k1e", // p[5]
"k1f", // p[6]
"k1g", // p[7]
"k1r", // p[8]
"k21", // p[9]
"k22", // p[10]
"k3", // p[11]
"k_IRSiP_DosR", // p[12]
"km2", // p[13]
"km3", // p[14]
"observableParameter1_IR1_P", // p[15]
"observableParameter1_IRS1_P", // p[16]
"observableParameter1_IRS1_P_DosR", // p[17]
"noiseParameter1_IR1_P", // p[18]
"noiseParameter1_IRS1_P", // p[19]
"noiseParameter1_IRS1_P_DosR", // p[20]
};

std::array<const char*, 4> fixedParameterNames = {
    "insulin_time_1", // k[0]
"insulin_time_2", // k[1]
"insulin_dose_1", // k[2]
"insulin_dose_2", // k[3]
};

std::array<const char*, 9> stateNames = {
    "IR", // x_rdata[0]
"IRins", // x_rdata[1]
"IRp", // x_rdata[2]
"IRiP", // x_rdata[3]
"IRi", // x_rdata[4]
"IRS", // x_rdata[5]
"IRSiP", // x_rdata[6]
"X", // x_rdata[7]
"Xp", // x_rdata[8]
};

std::array<const char*, 3> observableNames = {
    "IR1_P", // y[0]
"IRS1_P", // y[1]
"IRS1_P_DosR", // y[2]
};

std::array<const ObservableScaling, 3> observableScalings = {
    ObservableScaling::lin, // y[0]
ObservableScaling::lin, // y[1]
ObservableScaling::lin, // y[2]
};

std::array<const char*, 12> expressionNames = {
    "insulin", // w[0]
"flux_v1_v_0", // w[1]
"flux_v2_v_1", // w[2]
"flux_v3_v_2", // w[3]
"flux_v4_v_3", // w[4]
"flux_v5_v_4", // w[5]
"flux_v6_v_5", // w[6]
"flux_v7_v_6", // w[7]
"flux_v8_v_7", // w[8]
"flux_v9_v_8", // w[9]
"flux_v10_v_9", // w[10]
"flux_v11_v_10", // w[11]
};

std::array<const char*, 21> parameterIds = {
    "k1a", // p[0]
"k1aBasic", // p[1]
"k1b", // p[2]
"k1c", // p[3]
"k1d", // p[4]
"k1e", // p[5]
"k1f", // p[6]
"k1g", // p[7]
"k1r", // p[8]
"k21", // p[9]
"k22", // p[10]
"k3", // p[11]
"k_IRSiP_DosR", // p[12]
"km2", // p[13]
"km3", // p[14]
"observableParameter1_IR1_P", // p[15]
"observableParameter1_IRS1_P", // p[16]
"observableParameter1_IRS1_P_DosR", // p[17]
"noiseParameter1_IR1_P", // p[18]
"noiseParameter1_IRS1_P", // p[19]
"noiseParameter1_IRS1_P_DosR", // p[20]
};

std::array<const char*, 4> fixedParameterIds = {
    "insulin_time_1", // k[0]
"insulin_time_2", // k[1]
"insulin_dose_1", // k[2]
"insulin_dose_2", // k[3]
};

std::array<const char*, 9> stateIds = {
    "IR", // x_rdata[0]
"IRins", // x_rdata[1]
"IRp", // x_rdata[2]
"IRiP", // x_rdata[3]
"IRi", // x_rdata[4]
"IRS", // x_rdata[5]
"IRSiP", // x_rdata[6]
"X", // x_rdata[7]
"Xp", // x_rdata[8]
};

std::array<const char*, 3> observableIds = {
    "IR1_P", // y[0]
"IRS1_P", // y[1]
"IRS1_P_DosR", // y[2]
};

std::array<const char*, 12> expressionIds = {
    "insulin", // w[0]
"flux_v1_v_0", // w[1]
"flux_v2_v_1", // w[2]
"flux_v3_v_2", // w[3]
"flux_v4_v_3", // w[4]
"flux_v5_v_4", // w[5]
"flux_v6_v_5", // w[6]
"flux_v7_v_6", // w[7]
"flux_v8_v_7", // w[8]
"flux_v9_v_8", // w[9]
"flux_v10_v_9", // w[10]
"flux_v11_v_10", // w[11]
};

std::array<int, 9> stateIdxsSolver = {
    0, 1, 2, 3, 4, 5, 6, 7, 8
};

std::array<bool, 4> rootInitialValues = {
    true, true, true, true
};

} // namespace model_Brannmark_JBC2010

} // namespace amici
