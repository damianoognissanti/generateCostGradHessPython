#include "Bruno_JExpBot2016.h"
#include <array>

namespace amici {

namespace model_Bruno_JExpBot2016 {

std::array<const char*, 23> parameterNames = {
    "init_b10", // p[0]
"init_bcar", // p[1]
"init_bcry", // p[2]
"init_ohb10", // p[3]
"init_zea", // p[4]
"k5", // p[5]
"kb1", // p[6]
"kb2", // p[7]
"kc1", // p[8]
"kc2", // p[9]
"kc4", // p[10]
"k5_multiplier", // p[11]
"kb1_multiplier", // p[12]
"kb2_multiplier", // p[13]
"kc1_multiplier", // p[14]
"kc2_multiplier", // p[15]
"kc4_multiplier", // p[16]
"noiseParameter1_ob10", // p[17]
"noiseParameter1_obcar", // p[18]
"noiseParameter1_obcry", // p[19]
"noiseParameter1_obio", // p[20]
"noiseParameter1_oohb10", // p[21]
"noiseParameter1_ozea", // p[22]
};

std::array<const char*, 0> fixedParameterNames = {
    
};

std::array<const char*, 7> stateNames = {
    "beta-carotin", // x_rdata[0]
"cry", // x_rdata[1]
"beta-10", // x_rdata[2]
"beta-io", // x_rdata[3]
"OH-beta-10", // x_rdata[4]
"OH-beta-io", // x_rdata[5]
"zea", // x_rdata[6]
};

std::array<const char*, 6> observableNames = {
    "", // y[0]
"", // y[1]
"", // y[2]
"", // y[3]
"", // y[4]
"", // y[5]
};

std::array<const ObservableScaling, 6> observableScalings = {
    ObservableScaling::lin, // y[0]
ObservableScaling::lin, // y[1]
ObservableScaling::lin, // y[2]
ObservableScaling::lin, // y[3]
ObservableScaling::lin, // y[4]
ObservableScaling::lin, // y[5]
};

std::array<const char*, 6> expressionNames = {
    "flux_v1_ReactionName", // w[0]
"flux_v2_ReactionName", // w[1]
"flux_v3_ReactionName", // w[2]
"flux_v4_ReactionName", // w[3]
"flux_v5_ReactionName", // w[4]
"flux_v6_ReactionName", // w[5]
};

std::array<const char*, 23> parameterIds = {
    "init_b10", // p[0]
"init_bcar", // p[1]
"init_bcry", // p[2]
"init_ohb10", // p[3]
"init_zea", // p[4]
"k5", // p[5]
"kb1", // p[6]
"kb2", // p[7]
"kc1", // p[8]
"kc2", // p[9]
"kc4", // p[10]
"k5_multiplier", // p[11]
"kb1_multiplier", // p[12]
"kb2_multiplier", // p[13]
"kc1_multiplier", // p[14]
"kc2_multiplier", // p[15]
"kc4_multiplier", // p[16]
"noiseParameter1_ob10", // p[17]
"noiseParameter1_obcar", // p[18]
"noiseParameter1_obcry", // p[19]
"noiseParameter1_obio", // p[20]
"noiseParameter1_oohb10", // p[21]
"noiseParameter1_ozea", // p[22]
};

std::array<const char*, 0> fixedParameterIds = {
    
};

std::array<const char*, 7> stateIds = {
    "bcar", // x_rdata[0]
"bcry", // x_rdata[1]
"b10", // x_rdata[2]
"bio", // x_rdata[3]
"ohb10", // x_rdata[4]
"ohbio", // x_rdata[5]
"zea", // x_rdata[6]
};

std::array<const char*, 6> observableIds = {
    "ob10", // y[0]
"obcar", // y[1]
"obcry", // y[2]
"obio", // y[3]
"oohb10", // y[4]
"ozea", // y[5]
};

std::array<const char*, 6> expressionIds = {
    "flux_v1_ReactionName", // w[0]
"flux_v2_ReactionName", // w[1]
"flux_v3_ReactionName", // w[2]
"flux_v4_ReactionName", // w[3]
"flux_v5_ReactionName", // w[4]
"flux_v6_ReactionName", // w[5]
};

std::array<int, 7> stateIdxsSolver = {
    0, 1, 2, 3, 4, 5, 6
};

std::array<bool, 0> rootInitialValues = {
    
};

} // namespace model_Bruno_JExpBot2016

} // namespace amici
