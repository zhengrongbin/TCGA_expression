# data can be given by -i db/TCGA_XENA_exp_mat_with_pheno_survival.parquet or -dir db/cancer/, depending on whether only need one or several cancer data
# run basic function for genes, including tumor vs normal expression, distribution plot, survival analysis
for g in LINC00113 PBRM1 METTL3 METTL4 YTHDF1 
do
	python main.py basic -dir db/cancer/ -c all -p /Users/rongbinzheng/Documents/temp/liwei/out/Kidney -d -b -su -g $g -c TCGA-KIRC,TCGA-KIRP,TCGA-KICH
done

# do correlation for a pair gene
g1=LINC00113
g2=METTL3
python main.py corr -i  db/TCGA_XENA_exp_mat_with_pheno_survival.parquet -p ./test -g1 $g1 -g2 $g2 -c TCGA-KIRC,TCGA-KIRP,TCGA-KICH

# do correlation analysis for genes in pairwise
for g1 in LINC00113 PBRM1 METTL3 METTL4 YTHDF1
do
 for g2 in LINC00113 PBRM1 METTL3 METTL4 YTHDF1
	do
		if [[ g1 != g2 ]]; then
			python main.py corr -i db/TCGA_XENA_exp_mat_with_pheno_survival.parquet -p test -g1 $g1 -g2 $g2 -c TCGA-KIRC,TCGA-KIRP,TCGA-KICH
		fi
	done
done

# run advanced function, actually do correlation analysis in tumor and normal, and ranking score based on delta corr between tumor and normal
python main.py advanced -i db/TCGA_XENA_exp_mat_with_pheno_survival.parquet -p ./test -g PBRM1 -c all -j CTL
# partial correlation adjusted by CTL
