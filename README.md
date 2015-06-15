# axtontalens
Simple TALEN assembly R script and associated files 

The script takes as arguments the target sequence string, a file of monomer sequences in FASTA format, and a set of recipes (validated with IDT software) for combining these units into fragments which can be synthesised by IDT and used to build functional TALENs. 

An example command to generate a TALEN capable of targeting the sequence GATCCGGATCAG might be:

"makeFragments.R GATCCGGATCAG monomers.fa fragment_defs.txt"
