Protease Prediction Plan 

Criteria as of 12/23/2025 

Utilizing the prosperous plus prediction tool (http://prosperousplus.unimelb-biotools.cloud.edu.au/) 

Derive a prediction of protease cleavage for all proteins that are from the human protein https://www.proteinatlas.org/humanproteome/brain/human+brain 

And defined as either “tissue enriched” or “tissue enhanced” 

Using the Uniprot ID, the FASTA sequence for the selected proteins (i.e., the amino acid sequence) can be pulled from here: https://www.uniprot.org/uniprotkb/P10636/entry - this is the one for MAPT (tau) 

CRITERIA: 

TWO proteases: A01.009 (Cathepsin D) AND C13.004 (Legumain) 

min_length = 25 

max_exceptions = 2 (total across BOTH proteases) 

score_threshold = 0.3 

 

To narrow it down for biomarker testing, take the amino acid sequence identified by your analysis, and then do a further analysis on the sequence including that and 5 amino acids on either side. So, as an example on tau, the algorithm identified 511-543 (SGYSSPGSPGTPGSRSRTPSLPTPPTREPKKVA), but we would want to evaluate 506-548 (KSGDRSGYSSPGSPGTPGSRSRTPSLPTPPTREPKKVAVVRTP).  

 

On this subsectioned sequence, we would then generate mutants: 

Every "S" and "T" residue would be converted to a "P" and then the analysis run for cathepsin L (C01.032) 

Every "S" and "T" residue would be converted to a "E" and then the analysis run for cathepsin B (C01.060). 

The math should be: (C01.032Normal-C01.03AlltoP) + (C01.06Normal – C01.060AlltoE) = PredictedChange. 

We want the highest predicted change – we decided the threshold should be the same as pTau217, which is 0.475. 

 

 

 