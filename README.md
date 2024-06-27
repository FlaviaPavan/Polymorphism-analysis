# Analysis_pi_data

Python script developed with the help of Sophie Galina to calculate the average nucleotide diversity (π) of different parts of miRNA hairpin precursors and miRNA targets. The script also performs permutation tests to assess the probability that the differences in π observed between the parts are different from chance. For example, each nucleotide associated with its π is permuted in the hairpin, and then the average π of each hairpin region is calculated. This operation is repeated a large number of times, making it possible to define a confidence interval.

> Authors : Sophie Gallina 1, Flavia Pavan 1

> Affiliations : 1. CNRS, Univ. Lille, UMR 8198 – Evo-Eco-Paleo, F-59000 Lille, France

**Corresponding author:** flavia.pavan@univ-lille.fr

*Command line to run the script :*

python3 Analysis_pi_data.py --in_pi PI_VALUES_FILE --in_bed HAIRPIN_POSITONS_FILE --out_dir OUTPUT_DIR_NAME  --simul NUMBER_OF_ITERATIONS  --plot --xmin 0 --xmax 0.05

> --in_bed : input bed file with 4 columns (CHROM START END CATEGORY) (category = eg stem, loop, mir, mir_star ...)

> --in_pi : input pi file with 3 colomns (CHROM POS PI_VALUE)

> --out_dir : output dir for results

> --simul : number of iterations, default=0

> --plot : default = no plot", action="store_const", const="Y", default="N")

> --xmin : min limit for x axis in density curves, default=0 (type=float)

> --xmax : max limit for x axis in density curves, default=0.5 (type=float)

> --skip_bed_dup : skip duplicated bed line (same seqid+start+end), default=Y", choices=["Y", "N"]

Output files : in out_dir folder

- observed_per_cat.csv : cat	n_sites	n_sites_polymorphes	mean_pi
    
- observed_per_item.csv : item	size	n_sites	n_sites_polymorphes	mean_pi
    
- observed_per_mir_pos.csv : pos	n_sites	n_sites_polymorphes	mean_pi
    
- simul_N.csv : cat	n_sites	n_sites_polymorphes	observed_mean_pi	n_iter	count0	count	simulated_mean_pi	percentile 2.5	percentile 97.5	percentile 5	percentile 95
    
- simul_N_category : N lines with the mean pi value for each simul iteration
    
- plot : one page for observed values as boxplot, one for simulated values as density curve, one page for each category with simulated values as density curve + observed value as vertical bar


    
