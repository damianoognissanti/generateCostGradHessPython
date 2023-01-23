#include "Schwen_PONE2014.h"
#include <array>

namespace amici {

namespace model_Schwen_PONE2014 {

std::array<const char*, 27> parameterNames = {
    "ini_R1", // p[0]
"ini_R2fold", // p[1]
"ka1", // p[2]
"ka2fold", // p[3]
"kd1", // p[4]
"kd2fold", // p[5]
"kin", // p[6]
"kin2", // p[7]
"koff_unspec", // p[8]
"kon_unspec", // p[9]
"kout", // p[10]
"kout2", // p[11]
"kout_frag", // p[12]
"observableParameter1_observable_IR1", // p[13]
"observableParameter2_observable_IR1", // p[14]
"observableParameter1_observable_IR2", // p[15]
"observableParameter2_observable_IR2", // p[16]
"observableParameter1_observable_IRsum", // p[17]
"observableParameter2_observable_IRsum", // p[18]
"observableParameter1_observable_Insulin", // p[19]
"observableParameter2_observable_Insulin", // p[20]
"observableParameter3_observable_Insulin", // p[21]
"observableParameter4_observable_Insulin", // p[22]
"noiseParameter1_observable_IR1", // p[23]
"noiseParameter1_observable_IR2", // p[24]
"noiseParameter1_observable_IRsum", // p[25]
"noiseParameter1_observable_Insulin", // p[26]
};

std::array<const char*, 1> fixedParameterNames = {
    "init_Ins", // k[0]
};

std::array<const char*, 11> stateNames = {
    "Insulin", // x_rdata[0]
"Receptors_low", // x_rdata[1]
"Receptors_high", // x_rdata[2]
"IR_low", // x_rdata[3]
"IR_high", // x_rdata[4]
"IR_low_internalized", // x_rdata[5]
"IR_high_internalized", // x_rdata[6]
"Uptake1", // x_rdata[7]
"Uptake2", // x_rdata[8]
"InsulinFragments", // x_rdata[9]
"BoundUnspec", // x_rdata[10]
};

std::array<const char*, 4> observableNames = {
    "", // y[0]
"", // y[1]
"", // y[2]
"", // y[3]
};

std::array<const ObservableScaling, 4> observableScalings = {
    ObservableScaling::log10, // y[0]
ObservableScaling::log10, // y[1]
ObservableScaling::log10, // y[2]
ObservableScaling::log10, // y[3]
};

std::array<const char*, 14> expressionNames = {
    "flux_v1_ka1", // w[0]
"flux_v2_ka2fold", // w[1]
"flux_v3_v_2", // w[2]
"flux_v4_v_3", // w[3]
"flux_v5_kd1", // w[4]
"flux_v6_kd2fold", // w[5]
"flux_v7_v_6", // w[6]
"flux_v8_v_7", // w[7]
"flux_v9_v_8", // w[8]
"flux_v10_v_9", // w[9]
"flux_v11_v_10", // w[10]
"flux_v12_v_11", // w[11]
"flux_v13_v_12", // w[12]
"flux_v14_v_13", // w[13]
};

std::array<const char*, 27> parameterIds = {
    "ini_R1", // p[0]
"ini_R2fold", // p[1]
"ka1", // p[2]
"ka2fold", // p[3]
"kd1", // p[4]
"kd2fold", // p[5]
"kin", // p[6]
"kin2", // p[7]
"koff_unspec", // p[8]
"kon_unspec", // p[9]
"kout", // p[10]
"kout2", // p[11]
"kout_frag", // p[12]
"observableParameter1_observable_IR1", // p[13]
"observableParameter2_observable_IR1", // p[14]
"observableParameter1_observable_IR2", // p[15]
"observableParameter2_observable_IR2", // p[16]
"observableParameter1_observable_IRsum", // p[17]
"observableParameter2_observable_IRsum", // p[18]
"observableParameter1_observable_Insulin", // p[19]
"observableParameter2_observable_Insulin", // p[20]
"observableParameter3_observable_Insulin", // p[21]
"observableParameter4_observable_Insulin", // p[22]
"noiseParameter1_observable_IR1", // p[23]
"noiseParameter1_observable_IR2", // p[24]
"noiseParameter1_observable_IRsum", // p[25]
"noiseParameter1_observable_Insulin", // p[26]
};

std::array<const char*, 1> fixedParameterIds = {
    "init_Ins", // k[0]
};

std::array<const char*, 11> stateIds = {
    "Ins", // x_rdata[0]
"Rec1", // x_rdata[1]
"Rec2", // x_rdata[2]
"IR1", // x_rdata[3]
"IR2", // x_rdata[4]
"IR1in", // x_rdata[5]
"IR2in", // x_rdata[6]
"Uptake1", // x_rdata[7]
"Uptake2", // x_rdata[8]
"InsulinFragments", // x_rdata[9]
"BoundUnspec", // x_rdata[10]
};

std::array<const char*, 4> observableIds = {
    "observable_IR1", // y[0]
"observable_IR2", // y[1]
"observable_IRsum", // y[2]
"observable_Insulin", // y[3]
};

std::array<const char*, 14> expressionIds = {
    "flux_v1_ka1", // w[0]
"flux_v2_ka2fold", // w[1]
"flux_v3_v_2", // w[2]
"flux_v4_v_3", // w[3]
"flux_v5_kd1", // w[4]
"flux_v6_kd2fold", // w[5]
"flux_v7_v_6", // w[6]
"flux_v8_v_7", // w[7]
"flux_v9_v_8", // w[8]
"flux_v10_v_9", // w[9]
"flux_v11_v_10", // w[10]
"flux_v12_v_11", // w[11]
"flux_v13_v_12", // w[12]
"flux_v14_v_13", // w[13]
};

std::array<int, 11> stateIdxsSolver = {
    0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10
};

std::array<bool, 0> rootInitialValues = {
    
};

} // namespace model_Schwen_PONE2014

} // namespace amici
