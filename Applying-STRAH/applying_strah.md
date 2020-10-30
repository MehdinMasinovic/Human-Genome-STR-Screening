# Applying STRAH on the human reference genome

In my thesis, I analysed short tandem repeats across the human genome. In total, the analysis includes a whole human genome screening of 310 short-tandem repeat types. 

STRAH was used to screen the human genome. As the dataset is ~15GB big, I decided to only provide the way STRAH was applied instead of the dataset. If someone is interested in the actual datasets, please do not hesitate to contact me via e-mail (mehdin.masinovic2@gmail.com).

The short-tandem repeat analysis can be grouped in five categories: 

1. Mononucleotide repeats (A, C, G, T)
2. Dinucleotide repeats (AC, AG, AT, CT, CG, CT)
3. Trinucleotide repeats (all DNA-triplet combinations, g.e. AAC, ACA, CAA, etc.)
4. Tetranucleotide repeats (all distinct DNA-quadruplet combinations, g.e. AACC, AAAC, CAAA, etc.)
5. DNA motifs (a set of DNA motifs that have been analysed in literature before, such as AAAAN, GGGGN, GGGNGGG, etc.)

In all categories, STRAH is used to analyze the human genome via the following function call:

```{R}
STR_analysis(nr.STRs, max.nr.STRs, nr.mismatch = 0, STR = letter, addToHs = 500, lens.grey = 0:5*1000,  reverse.comp = FALSE, pos_matrix, species=BSgenome.Hsapiens.UCSC.hg38::Hsapiens)
```

The main two functions of STRAH (STR_analysis() and STR_detection()) are included in this folder for reproducibility (in case STRAH is updated in the future). 

STR_analysis() is the main function in STRAH. All short-tandem repeat types are screened for using the abovementioned command, differing in the minimum repeat length. In all categories, we screen the human genome using the double-strand break map (DSB-map) of the work of Pratto et al. (https://science.sciencemag.org/content/346/6211/1256442). We define recombination maps as ± 500 base pairs up and downstream from the hotspot coordinates identified by Pratto et al.




