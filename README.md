# IBD_project
1. SNV Calling pipeline could be found in Call_SNP.
2. For the obtained SNV position, it is necessary to integrate the snps of each sample into a table.
3. Prepare the gff3 file for gene and protein annotation of related bacteria, and use Snp_to_gene.py to count and merge the number of SNVs in each gene.
4. Use Randomforest.pred-prob.simple--Shi.R for random forest analysis and drawing. randomforests_util.R and util.R should be placed in the same folder as Randomforest.pred-prob.simple--Shi.R
