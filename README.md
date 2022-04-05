# HIV time to infection estimation

Calculate next-generation sequence diversity scores including
- average pairwise distance diversity score  
- ambiguity score (fraction of polymorphic sites) 

based on the formulas in [Carlisle 2019](https://pubmed.ncbi.nlm.nih.gov/30835266/#:~:text=Conclusions%3A%20Viral%20diversity%20based%20on,absolute%20error%20of%20%3C1%20year.) and [Puller 2017](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1005775#:~:text=The%20error%20increased%20only%20slightly,contrast%20to%20most%20alternative%20biomarkers.).

# Install Software Dependencies
- prinseq
- bwa
- lofreq
- samtools
- pandas (python)

 The environment.yml file can be helpful. 
