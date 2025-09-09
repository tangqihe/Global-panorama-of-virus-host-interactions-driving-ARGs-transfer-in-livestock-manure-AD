# 02Reads_module.sh
#!/bin/bash
#SBATCH --job-name=myjob_name
#SBATCH --chdir=/work/<USERNAME>
#SBATCH --output=/work/%u/%x-%j.log
#SBATCH --time=0-00:30:00
source deactivate
conda activate base
cd /media/kunyu/data6/05Zeev
conda activate tools
#-------------------------01SARG------------------------------------#
diamond makedb --in 04reads_based_analysis/01SARG/SARG.2.2.fasta -d 04reads_based_analysis/01SARG/SARG.2.2_nr
for i in 01Reads_fq/*.fastq
 do
     diamond blastx -d 04reads_based_analysis/01SARG/SARG.2.2_nr -q $i -o $i-SARG.txt --evalue 1e-5 --query-cover 75 --id 90 -k 1
     mv $i-SARG.txt 04reads_based_analysis/01SARG/
done
cd 04reads_based_analysis/01SARG
echo -e "SARG_family\t$(ls *-SARG.txt | sed 's/-SARG.txt//;s/^/count /' | tr '\n' '\t')" > Zeev_SARG.xls
for i in $(cut -f 1 SARG.2.2.txt | sort -u); do
  echo -ne "$i\t" >> Zeev_SARG.xls
  for file in *-SARG.txt; do
    prefix=${file%-SARG.txt}
    count=$(grep -c -F "$i" "$file")
    echo -ne "$count\t" >> Zeev_SARG.xls
  done
  echo >> Zeev_SARG.xls
done
awk -F'\t' 'BEGIN {OFS="\t"} NR==1 {print} NR>1 {sum=0; for (i=2; i<=NF; i++) {sum+=$i} if (sum!=0) {print}}' Zeev_SARG.xls > Zeev_SARG_final.xls
csvtk join -t --left-join -f 1 Zeev_SARG_final.xls SARG_mapping-20220906.txt > Zeev_SARG_mapping.xls
#-------------------------02CARD------------------------------------#
cd ..
cd ..
diamond makedb --in 04reads_based_analysis/02card/CARD3.1.4.fasta -d 04reads_based_analysis/02card/card3.1.4_nr
for i in 01Reads_fq/*.fastq
 do
     diamond blastx -d 04reads_based_analysis/02card/card3.1.4_nr -q $i -o $i-CARD.txt --evalue 1e-5 --query-cover 75 --id 90 -k 1
     mv $i-CARD.txt 04reads_based_analysis/02card/
done
cd 04reads_based_analysis/02card
echo -e "CARD_family\t$(ls *-CARD.txt | sed 's/-CARD.txt//;s/^/count /' | tr '\n' '\t')" > Zeev_CARD.xls
for i in $(cut -f 3 CARD3.1.4.txt | sort -u); do
  echo -ne "$i\t" >> Zeev_CARD.xls
  for file in *-CARD.txt; do
    prefix=${file%-CARD.txt}
    count=$(grep -c -F "$i" "$file")
    echo -ne "$count\t" >> Zeev_CARD.xls
  done
  echo >> Zeev_CARD.xls
done
awk -F'\t' 'BEGIN {OFS="\t"} NR==1 {print} NR>1 {sum=0; for (i=2; i<=NF; i++) {sum+=$i} if (sum!=0) {print}}' Zeev_CARD.xls > Zeev_CARD_final.xls
csvtk join -t --left-join -f 1 Zeev_CARD_final.xls CARD_mapping-20220905.txt > Zeev_CARD_mapping.xls
#-------------------------03bacmet------------------------------------#
cd ..
cd ..
diamond makedb --in 04reads_based_analysis/03bacmet/BacMet2.fasta -d 04reads_based_analysis/03bacmet/bacmet2_nr
for i in 01Reads_fq/*.fastq
 do
     diamond blastx -d 04reads_based_analysis/03bacmet/bacmet2_nr -q $i -o $i-BacMet2.txt --evalue 1e-5 --query-cover 75 --id 90 -k 1
     mv $i-BacMet2.txt 04reads_based_analysis/03bacmet/
done
cd 04reads_based_analysis/03bacmet
echo -e "BacMet2_family\t$(ls *-BacMet2.txt | sed 's/-BacMet2.txt//;s/^/count /' | tr '\n' '\t')" > Zeev_BacMet2.xls
for i in $(cut -f 1 BacMet2.txt | sort -u); do
  echo -ne "$i\t" >> Zeev_BacMet2.xls
  for file in *-BacMet2.txt; do
    prefix=${file%-BacMet2.txt}
    count=$(grep -c -F "$i" "$file")
    echo -ne "$count\t" >> Zeev_BacMet2.xls
  done
  echo >> Zeev_BacMet2.xls
done
awk -F'\t' 'BEGIN {OFS="\t"} NR==1 {print} NR>1 {sum=0; for (i=2; i<=NF; i++) {sum+=$i} if (sum!=0) {print}}' Zeev_BacMet2.xls > Zeev_BacMet2_final.xls
csvtk join -t --left-join -f 1 Zeev_BacMet2_final.xls BacMet21_EXP.753.mapping.txt > Zeev_BacMet2_mapping.xls
#-------------------------06victors------------------------------------#
cd ..
cd ..
diamond makedb --in 04reads_based_analysis/06victors/victors_pro.fasta -d 04reads_based_analysis/06victors/victors_nr
for i in 01Reads_fq/*.fastq
do
     diamond blastx -d 04reads_based_analysis/06victors/victors_nr -q $i -o $i-victors.txt --evalue 1e-5 --query-cover 75 --id 90 -k 1
     mv $i-victors.txt 04reads_based_analysis/06victors/
done
cd 04reads_based_analysis/06victors
echo -e "victors_family\t$(ls *-victors.txt | sed 's/-victors.txt//;s/^/count /' | tr '\n' '\t')" > Zeev_victors.xls
for i in $(cut -f 1 victors.txt | sort -u); do
  echo -ne "$i\t" >> Zeev_victors.xls
  for file in *-victors.txt; do
    prefix=${file%-victors.txt}
    count=$(grep -c -F "$i" "$file")
    echo -ne "$count\t" >> Zeev_victors.xls
  done
  echo >> Zeev_victors.xls
done
awk -F'\t' 'BEGIN {OFS="\t"} NR==1 {print} NR>1 {sum=0; for (i=2; i<=NF; i++) {sum+=$i} if (sum!=0) {print}}' Zeev_victors.xls > Zeev_victors_final.xls
csvtk join -t --left-join -f 1 Zeev_victors_final.xls victors_mapping-20220918.txt > Zeev_victors_mapping.xls
#-------------------------07kraken2_16S------------------------------------#
cd ..
cd ..
for i in 01Reads_fq/*.fastq
do
     kraken2 --db /media/kunyu/data2/db/metawrap/kraken2/16S_Greengenes_k2db $i --output $i.txt --report $i.report.txt --classified-out $i.16s.fasta --threads 20
     mv $i.report.txt 04reads_based_analysis/07kraken2/
     mv $i.16s.fasta 04reads_based_analysis/07kraken2/
     mv $i.txt 04reads_based_analysis/07kraken2/
done
cd 04reads_based_analysis/07kraken2
for file in *.fasta
do
  echo "$(basename $file .fasta): $(grep -c "^@" $file) reads" >> 16s_reads_number.txt
done
#-------------------------12Virus------------------------------------#
cd ..
cd ..
for i in 01Reads_fq/*.fastq
do
     kraken2 --db /media/kunyu/data2/db/metawrap/kraken2/k2_viral_20230314 $i --output $i.txt --report $i.report.txt --classified-out $i.viral.fasta --threads 30
     mv $i.report.txt 04reads_based_analysis/12virus/
     mv $i.viral.fasta 04reads_based_analysis/12virus/
     mv $i.txt 04reads_based_analysis/12virus/
done
cd 04reads_based_analysis/12virus
combine_kreports.py -r *.report.txt -o Zeev_kraken2_virus.txt
#----------------------------05MGEs-------------------------------#
cd ..
cd ..
for F in 01Reads_fq/*_1.fastq
do	
	R=${F%_*}_2.fastq
	BASE=${F##*/}
	SAMPLE=${BASE%_*}
	TMPDIR=/media/kunyu/data6/05Zeev/tmp coverm make -r 04reads_based_analysis/05MGEs/MGEs.fasta -1 $F -2 $R -o 04reads_based_analysis/05MGEs/ -t 30
	coverm filter -b 04reads_based_analysis/05MGEs/MGEs.fasta.$BASE.bam -o 04reads_based_analysis/05MGEs/$SAMPLE.bam --min-read-percent-identity 95 --min-read-aligned-percent 75 --threads 30
	coverm contig --bam-files 04reads_based_analysis/05MGEs/$SAMPLE.bam -o 04reads_based_analysis/05MGEs/$SAMPLE.rpkm.txt --trim-min 10 --trim-max 90 --min-read-percent-identity 95 --min-read-aligned-percent 75 -m rpkm -t 30
	coverm contig --bam-files 04reads_based_analysis/05MGEs/$SAMPLE.bam -o 04reads_based_analysis/05MGEs/$SAMPLE.count.txt --min-covered-fraction 0 -m count -t 30
	coverm contig --bam-files 04reads_based_analysis/05MGEs/$SAMPLE.bam -o 04reads_based_analysis/05MGEs/$SAMPLE.bases.txt --min-covered-fraction 0 -m covered_bases -t 30
done
cd 04reads_based_analysis/05MGEs
paste *.rpkm.txt >Zeev.rpkm.csv
paste *.count.txt >Zeev.count.csv
paste *.bases.txt >Zeev.bases.csv
#----------------------------14crassphage-------------------------------#
cd ..
cd ..
for F in 01Reads_fq/*_1.fastq
do	
	R=${F%_*}_2.fastq
	BASE=${F##*/}
	SAMPLE=${BASE%_*}
	TMPDIR=/media/kunyu/data6/05Zeev/tmp coverm make -r 04reads_based_analysis/14crassphage/crassphage.fasta -1 $F -2 $R -o 04reads_based_analysis/14crassphage/ -t 30
	coverm filter -b 04reads_based_analysis/14crassphage/crassphage.fasta.$BASE.bam -o 04reads_based_analysis/14crassphage/$SAMPLE.bam --min-read-percent-identity 95 --min-read-aligned-percent 75 --threads 30
	coverm contig --bam-files 04reads_based_analysis/14crassphage/$SAMPLE.bam -o 04reads_based_analysis/14crassphage/$SAMPLE.rpkm.txt --trim-min 10 --trim-max 90 --min-read-percent-identity 95 --min-read-aligned-percent 75 -m rpkm -t 30
	coverm contig --bam-files 04reads_based_analysis/14crassphage/$SAMPLE.bam -o 04reads_based_analysis/14crassphage/$SAMPLE.count.txt --min-covered-fraction 0 -m count -t 30
	coverm contig --bam-files 04reads_based_analysis/14crassphage/$SAMPLE.bam -o 04reads_based_analysis/14crassphage/$SAMPLE.bases.txt --min-covered-fraction 0 -m covered_bases -t 30
done
cd 04reads_based_analysis/14crassphage
paste *.rpkm.txt >Zeev.rpkm.csv
paste *.count.txt >Zeev.count.csv
paste *.bases.txt >Zeev.bases.csv
#--------------------------------metaphlan4-------------------------#
cd ..
cd ..
conda deactivate
conda activate metaphlan4
for F in 01Reads_fq/*_1.fastq
do
	R=${F%_*}_2.fastq
	BASE=${F##*/}
	SAMPLE=${BASE%_*}
	metaphlan $F,$R --input_type fastq --nproc 50 --bowtie2out 10metaphlan4/$SAMPLE.bowtie2out --output_file 10metaphlan4/$SAMPLE.profile.txt --bowtie2db /media/kunyu/data2/db/metaphlan4 --tmp_dir /media/kunyu/data6/05Zeev/tmp
done
for F in 10metaphlan4/*.bowtie2out
do
	BASE=${F##*/}
	SAMPLE=${BASE%.bowtie2out*}
	metaphlan $F --nproc 50 --input_type bowtie2out --output_file 10metaphlan4/$SAMPLE.profile.txt --bowtie2db /media/kunyu/data2/db/metaphlan4
done
merge_metaphlan_tables.py 10metaphlan4/*.txt > 10metaphlan4/merged_abundance_table.tsv
cd 10metaphlan4
Rscript /home/kunyu/miniconda3/envs/metaphlan4/lib/python3.7/site-packages/metaphlan/utils/calculate_diversity.R -f merged_abundance_table.tsv -d beta -m bray-curtis
Rscript /home/kunyu/miniconda3/envs/metaphlan4/lib/python3.7/site-packages/metaphlan/utils/calculate_diversity.R -f merged_abundance_table.tsv -d alpha -m richness
Rscript /home/kunyu/miniconda3/envs/metaphlan4/lib/python3.7/site-packages/metaphlan/utils/calculate_diversity.R -f merged_abundance_table.tsv -d alpha -m shannon
Rscript /home/kunyu/miniconda3/envs/metaphlan4/lib/python3.7/site-packages/metaphlan/utils/calculate_diversity.R -f merged_abundance_table.tsv -d alpha -m simpson
Rscript /home/kunyu/miniconda3/envs/metaphlan4/lib/python3.7/site-packages/metaphlan/utils/calculate_diversity.R -f merged_abundance_table.tsv -d alpha -m gini
conda deactivate
conda activate graphlan
export2graphlan.py --skip_rows 1,2 -i 10metaphlan4/merged_abundance_table.txt --tree 10metaphlan4/merged_abundance.tree.txt --annotation 10metaphlan4/merged_abundance.annot.txt --most_abundant 100 --abundance_threshold 1 --least_biomarkers 10 --annotations 5,6 --external_annotations 7 --min_clade_size 1
graphlan_annotate.py --annot 10metaphlan4/merged_abundance.annot.txt 10metaphlan4/merged_abundance.tree.txt 10metaphlan4/merged_abundance.xml
graphlan.py --dpi 300 10metaphlan4/merged_abundance.xml 10metaphlan4/merged_abundance.pdf --external_legends
