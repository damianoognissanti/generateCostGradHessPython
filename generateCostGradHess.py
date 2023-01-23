# Generates ODE solutions
import os
import numpy as np
import pypesto
import pypesto.petab
import pandas as pd
import amici

# Select a model to Test
#modelName = "Test_Cvijoviclab2022simple"
#modelName = "Test_Cvijoviclab2022"
modelName = "Boehm_JProteomeRes2014"
#modelName = "Bachmann_MSB2011"
longModelName = "model_" + modelName

inputFolderModel = "Data"
outputFolderBase = ""
inputFolder = os.path.join(outputFolderBase, longModelName)
outputFolder = os.path.join(outputFolderBase, longModelName)

outputFileNameCost = os.path.join(outputFolder,"Cost.csv")
outputFileNameGrad = os.path.join(outputFolder,"Grad.csv")
outputFileNameHess = os.path.join(outputFolder,"Hess.csv")

# Remove file if exists
if os.path.exists(outputFileNameCost):
    os.remove(outputFileNameCost)
if os.path.exists(outputFileNameGrad):
    os.remove(outputFileNameGrad)
if os.path.exists(outputFileNameHess):
    os.remove(outputFileNameHess)

# the yaml configuration file links to all needed files
yamlConfig = os.path.join(inputFolderModel, longModelName, modelName + ".yaml")

# import the model from yaml-file
importer = pypesto.petab.PetabImporter.from_yaml(yamlConfig)
petabProblem = importer.petab_problem

# open CSV-files with the Parameters
parameterPathScaled = os.path.join(inputFolder, "Params.csv")
parameterDFScaled = pd.read_csv(parameterPathScaled)

# open measurementData to extract timepoints
measurementDataPath = os.path.join(inputFolderModel, longModelName, "measurementData_" + modelName + ".tsv")
measurementDataDF = pd.read_csv(measurementDataPath, sep='\t', header=0)
timepoints = np.sort(np.unique(measurementDataDF.time.values))

# open experimental conditions
experimentalConditionPath = os.path.join(inputFolderModel, longModelName, "experimentalCondition_" + modelName + ".tsv")
experimentalConditionDF = pd.read_csv(experimentalConditionPath, sep='\t', header=0)

# create and set options for the objective function
obj = importer.create_objective()
obj.fim_for_hess = True
#obj.amici_solver.setMaxSteps(10000)
#obj.amici_solver.setRelativeTolerance(1e-15)
#obj.amici_solver.setAbsoluteTolerance(1e-15)
#obj.amici_solver.setLinearSolver(1)

# create model
model = obj.amici_model

# extract parameter names (different order and other names than in parameters file)
modelParameterNames = model.getParameterNames()

# extract expcond parameter names (different order and other names than in parameters file)
expCondModelParameterNames = model.getFixedParameterNames()

# Extract all problem IDs to make sure we use the correct order in Python
AllIds = petabProblem.x_ids
HessId = []
for row in AllIds:
    for col in AllIds:
        HessId.append(row+col)

# Create empty dataframes for Cost, Grad, Hess
costCols = ['Cost']
gradCols = AllIds
hessCols = HessId
CostDF = pd.DataFrame(columns=costCols)
GradDF = pd.DataFrame(columns=gradCols)
HessDF = pd.DataFrame(columns=hessCols)

# Csv Observables header strings
headerStringObs = ['time','observableId','calculatedObs']
# Csv ODE header strings
headerStringSol = np.array(model.getStateIds())
# Prepend the string 'time'
headerStringSol = np.insert(headerStringSol,0,'time')

for testCaseIndex in parameterDFScaled['Id']:
    print("Row " + str(testCaseIndex))
    ### Generate Cost Grad Hess
    parameterVector = []
    # To make sure order is the same
    for parid in AllIds:
        parameterValue = parameterDFScaled[parid][testCaseIndex]
        parameterVector.append(parameterValue)
    #    
    ret = obj(
            parameterVector,
            mode="mode_fun",
            sensi_orders=(0, 1, 2),
            return_dict=True
    )
    #
    PythonCost = ret['fval']
    PythonGrad = ret['grad']
    PythonHess = ret['hess']
    #
    outputCost = PythonCost
    outputGrad = PythonGrad
    #
    # Flatten Hessian
    outputHess = []
    hessRows, hessCols = PythonHess.shape
    for row in range(0,hessRows):
        for col in range(0,hessCols):
            outputHess.append(PythonHess[row][col])
    #
    CostDF.loc[len(CostDF.index)] = outputCost
    GradDF.loc[len(GradDF.index)] = outputGrad
    HessDF.loc[len(HessDF.index)] = outputHess


CostDF.to_csv(outputFileNameCost, header=CostDF.columns, index_label='Id', sep=',', mode='a')
GradDF.to_csv(outputFileNameGrad, header=GradDF.columns, index_label='Id', sep=',', mode='a')
HessDF.to_csv(outputFileNameHess, header=HessDF.columns, index_label='Id', sep=',', mode='a')
