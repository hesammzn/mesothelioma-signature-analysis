### **Hesam Moazzen**  

`MSc Bioinformatics, University of Leicester`


### Overview
In this repository, scripts for extracting and visualizing SBS, CNV and Chromothripsis signatures are published.

---

### Utilized tools:

**SBS:**
SigProfilerMatrixGenerator
SigProfilerAssignment
CloneSig
MutationalPatterns

**CNV:**
panConusig
Sigminer
SigProfilerMatrixGenerator
SigProfilerPlotting

**Chromothripsis:**
mclust and hdp packages in R

---

### Repository Structure

There are 3 folders in this repository:

- **SBS**  
  Contains all of the scripts for extracting SBS signatures with 3 methods: SigProfiler tools, CloneSig and MutationaPatterns.

- **CNV**  
  Contains scripts for extracting CNV signatures with panCounsig-Sigminer and SigProfiler tools.

- **Chromothripsis**  
  Includes the code for analyzing chromothripsis events based on Maclachlan et al., paper.

```
├───Chromothripsis
│   └───Figures
├───CNV
│   ├───panConusig
│   ├───Sigminer
│   │   └───Figures
│   ├───SigProfilerMatrixGenerator
│   └───SigProfilerPlotting
│       └───Figures
└───SBS
    ├───CloneSig
    │   └───Figures
    ├───MutationalPatterns
    │   └───Figures
    ├───SigProfilerAssignment
    │   └───Figures
    └───SigProfilerMatrixGenerator
```

---

### References

**SigProfiler tools:** https://cancer.sanger.ac.uk/signatures/tools/

**CloneSig:** https://github.com/judithabk6/clonesig

**MutationalPatterns:** https://github.com/UMCUGenetics/MutationalPatterns

**Sigminer:** https://github.com/ShixiangWang/sigminer

**panConusig:** https://github.com/UCL-Research-Department-of-Pathology/panConusig

**Chromothripsis:** Maclachlan, K.H., Rustad, E.H., Derkach, A. et al. Copy number signatures predict chromothripsis and clinical outcomes in newly diagnosed multiple myeloma. Nat Commun 12, 5172 (2021). https://doi.org/10.1038/s41467-021-25469-8
