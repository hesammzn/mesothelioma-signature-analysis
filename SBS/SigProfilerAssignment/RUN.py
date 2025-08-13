#SigProfilerAssignment
from SigProfilerAssignment import Analyzer as Analyze

Analyze.cosmic_fit(
    samples="C:/Users/hesam/Desktop/irp5/277/Sig/SigGenerator/Subclonal/output/SBS/MEDUSA.SBS96.exome",                #Input file path
    output="C:/Users/hesam/Desktop/irp5/277/Sig/SigAssign_3.4_2/Subclonal/",                                           #Output path
    input_type="matrix",                    #Input file type
    context_type="96",                      #Number of rows (contents)
    collapse_to_SBS96=True,                 #Define if the context has 96 rows (contents)
    cosmic_version=3.3,                     #Define the reference COSMIC signature version
    exome=True,                             #True if the data is WES
    genome_build="GRCh37",                  #Genome version
)