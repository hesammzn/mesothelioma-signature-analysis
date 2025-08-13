from SigProfilerMatrixGenerator.scripts import SigProfilerMatrixGeneratorFunc as matGen

#Running SigMatrixGenerator
matrices = matGen.SigProfilerMatrixGeneratorFunc(
    "MEDUSA",                                                               #Project name
    "GRCh37",                                                               #Genome
    "C:/Users/hesam/Desktop/irp5/277/Sig/SigGenerator/Clonal/",             #Path to the input table
    exome=True                                                              #True if the data is WES
)
print("SigProfiler matrix generation complete.")