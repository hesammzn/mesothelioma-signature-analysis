from SigProfilerMatrixGenerator.scripts import CNVMatrixGenerator as scna

# Set parameters
file_type = "ASCAT"
input_file = "C:/Users/hesam/Desktop/irp5/277/Sig_CNV/Subclonal/Input_Subclonal_CNV.txt"
output_path = "C:/Users/hesam/Desktop/irp5/277/Sig_CNV/Subclonal/"
project = "277"

# Generate CNV matrix
scna.generateCNVMatrix(file_type, input_file, project, output_path)
