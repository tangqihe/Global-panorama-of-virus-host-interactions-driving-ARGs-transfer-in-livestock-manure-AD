#!/bin/bash
#SBATCH --job-name=myjob_name
#SBATCH --chdir=/work/<USERNAME>
#SBATCH --output=/work/%u/%x-%j.log
#SBATCH --time=0-00:30:00
#nohup bash 04Contigs_module.sh > 04Contigs_module.log 2>&1 &
source deactivate
conda activate base
conda activate tools
cd /media/kunyu/data6/05Zeev
#-------------------------Contigs_prodigal------------------------------------#
for F in 01Reads_fq/*_1.fastq
do
	R=${F%_*}_2.fastq
	BASE=${F##*/}
	SAMPLE=${BASE%_*}
	prodigal -i 05ASSEMBLY/$SAMPLE/final_assembly.fasta -o 05ASSEMBLY/$SAMPLE/$SAMPLE.genes -a 05ASSEMBLY/$SAMPLE/$SAMPLE.pro.fa -d 05ASSEMBLY/$SAMPLE/$SAMPLE.fa -p meta
	mv 05ASSEMBLY/"$SAMPLE"/"$SAMPLE".genes 05ASSEMBLY/ORFs_genes/
	mv 05ASSEMBLY/"$SAMPLE"/"$SAMPLE".pro.fa 05ASSEMBLY/ORFs_protein/
	mv 05ASSEMBLY/"$SAMPLE"/"$SAMPLE".fa 05ASSEMBLY/ORFs_nucl/	
done
#--------------------------------Zeev-SARG-------------------------#
cd 05ASSEMBLY
diamond makedb --in ORFs_nucl/01SARG/SARG.2.2.fasta -d ORFs_nucl/01SARG/SARG.2.2_nr
for F in ORFs_nucl/*.fa;
do
	BASE=${F##*/}	
	SAMPLE=${BASE%.*}
	diamond blastx -d ORFs_nucl/01SARG/SARG.2.2_nr -q $F -o $F-SARG.txt --evalue 1e-10 --query-cover 70 --id 80 -k 1
	mv $F-SARG.txt ORFs_nucl/01SARG/"$SAMPLE".fa.txt
	cut -f 1,2 ORFs_nucl/01SARG/"$SAMPLE".fa.txt | awk 'BEGIN{FS=OFS="\t"}{gsub("_","\t",$1)}1' | cut -f 1,2,3,4,5,6,8 | awk -F " " '{print $1"_"$2"_"$3"_"$4"_"$5"_"$6" "$7}' | sort -k 1n >ORFs_nucl/01SARG/"$SAMPLE"_ARGs.txt
	sort ORFs_nucl/01SARG/"$SAMPLE"_ARGs.txt >ORFs_nucl/01SARG/"$SAMPLE"_ARGs_contigs.txt
	cut -d' ' -f1 ORFs_nucl/01SARG/"$SAMPLE"_ARGs_contigs.txt | sort -u >ORFs_nucl/01SARG/"$SAMPLE"_ARGs_list.txt
	seqtk subseq /media/kunyu/data2/03Zeev/05ASSEMBLY/$SAMPLE/final_assembly.fasta ORFs_nucl/01SARG/"$SAMPLE"_ARGs_list.txt | sed "s/k141/$SAMPLE/g" >ORFs_nucl/01SARG/"$SAMPLE"_ARGs.fa
done
cat ORFs_nucl/01SARG/*_ARGs.fa >ORFs_nucl/01SARG/Zeev_SARG_ARCs.fa
#--------------------------------Zeev-CARD-------------------------#
diamond makedb --in ORFs_nucl/02card/CARD3.1.4.fasta -d ORFs_nucl/02card/card3.1.4_nr
for F in ORFs_nucl/*.fa;
do
	BASE=${F##*/}	
	SAMPLE=${BASE%.*}
	diamond blastx -d ORFs_nucl/02card/card3.1.4_nr -q $F -o $F-card.txt --evalue 1e-5 --query-cover 70 --id 80 -k 1
	mv $F-card.txt ORFs_nucl/02card/"$SAMPLE".fa.txt
	cut -f 1,2 ORFs_nucl/02card/"$SAMPLE".fa.txt | awk 'BEGIN{FS=OFS="\t"}{gsub("_","\t",$1)}1' | cut -f 1,2,3,4,5,6,8 | awk -F " " '{print $1"_"$2"_"$3"_"$4"_"$5"_"$6" "$7}' | sort -k 1n >ORFs_nucl/02card/"$SAMPLE"_ARGs.txt
	sort ORFs_nucl/02card/"$SAMPLE"_ARGs.txt >ORFs_nucl/02card/"$SAMPLE"_ARGs_contigs.txt
	cut -d' ' -f1 ORFs_nucl/02card/"$SAMPLE"_ARGs_contigs.txt | sort -u >ORFs_nucl/02card/"$SAMPLE"_ARGs_list.txt
	seqtk subseq /media/kunyu/data2/03Zeev/05ASSEMBLY/$SAMPLE/final_assembly.fasta ORFs_nucl/02card/"$SAMPLE"_ARGs_list.txt | sed "s/k141/$SAMPLE/g" >ORFs_nucl/02card/"$SAMPLE"_ARGs.fa
done
cat ORFs_nucl/02card/*_ARGs.fa >ORFs_nucl/02card/Zeev_card_ARCs.fa
cat ORFs_nucl/01SARG/Zeev_SARG_ARCs.fa ORFs_nucl/02card/Zeev_card_ARCs.fa >/media/kunyu/data2/03Zeev/07ARCs_analysis/Zeev_ARCs.fasta
cd /media/kunyu/data2/03Zeev/07ARCs_analysis
cd-hit-est -i Zeev_ARCs.fasta -o Zeev_ARCs_nr.fasta -c 1.0 -n 10 -M 0 -T 8
prodigal -i Zeev_ARCs_nr.fasta -o Zeev_ARCs_nr.genes -a Zeev_ARCs_nr.pro.fa -d Zeev_ARCs_nr.nul.fa -p meta
#--------------------------------Zeev-01SARG-------------------------#
diamond makedb --in 01SARG/SARG.2.2.fasta -d 01SARG/SARG.2.2_nr
diamond blastx -d 01SARG/SARG.2.2_nr -q Zeev_ARCs_nr.nul.fa -o 01SARG/Zeev_ARCs-SARG.txt --evalue 1e-10 --query-cover 70 --id 80 -k 1
cut -f 1,2 01SARG/Zeev_ARCs-SARG.txt | awk 'BEGIN{FS=OFS="\t"}{gsub("_","\t",$1)}1' | cut -f 1,2,3,4,5,6,8 | awk -F " " '{print $1"_"$2"_"$3"_"$4"_"$5"_"$6" "$7}' | sort -k 1n >01SARG/Zeev_ARCs-SARG_list.txt
cat 01SARG/Zeev_ARCs-SARG_list.txt | datamash -sW -g1 collapse 2 >01SARG/Zeev_ARCs-SARD_list2.txt
#--------------------------------Zeev-02card-------------------------#
diamond makedb --in 02card/CARD3.1.4.fasta -d 02card/card3.1.4_nr
diamond blastx -d 02card/card3.1.4_nr -q Zeev_ARCs_nr.nul.fa -o 02card/Zeev_ARCs-card.txt --evalue 1e-10 --query-cover 70 --id 80 -k 1
cut -f 1,2 02card/Zeev_ARCs-card.txt | awk 'BEGIN{FS=OFS="\t"}{gsub("_","\t",$1)}1' | cut -f 1,2,3,4,5,6,8 | awk -F " " '{print $1"_"$2"_"$3"_"$4"_"$5"_"$6" "$7}' | sort -k 1n >02card/Zeev_ARCs-card_list.txt
awk '{split ($2, T, "|"); $2 = T[4]}1' OFS="\t" 02card/Zeev_ARCs-card_list.txt | datamash -sW -g1 collapse 2 >02card/Zeev_ARCs-card_list2.txt
#--------------------------------Zeev-03bacmet-------------------------#
diamond makedb --in 03bacmet/BacMet2.fasta -d 03bacmet/bacmet2_nr
diamond blastx -d 03bacmet/bacmet2_nr -q Zeev_ARCs_nr.nul.fa -o 03bacmet/Zeev_ARCs-bacmet2.txt --evalue 1e-10 --query-cover 70 --id 80 -k 1
cut -f 1,2 03bacmet/Zeev_ARCs-bacmet2.txt | awk 'BEGIN{FS=OFS="\t"}{gsub("_","\t",$1)}1' | cut -f 1,2,3,4,5,6,8 | awk -F " " '{print $1"_"$2"_"$3"_"$4"_"$5"_"$6" "$7}' | sort -k 1n >03bacmet/Zeev_ARCs-bacmet2_list.txt
cat 03bacmet/Zeev_ARCs-bacmet2_list.txt | datamash -sW -g1 collapse 2 >03bacmet/Zeev_ARCs-bacmet2_list2.txt
#--------------------------------Zeev-04ICEs-------------------------#
diamond makedb --in 04ICEs/ICE3.0.fasta -d 04ICEs/ICE3.0_nr
diamond blastx -d 04ICEs/ICE3.0_nr -q Zeev_ARCs_nr.nul.fa -o 04ICEs/Zeev_ARCs-ICE.txt --evalue 1e-10 --query-cover 70 --id 80 -k 1
cut -f 1,2 04ICEs/Zeev_ARCs-ICE.txt | awk 'BEGIN{FS=OFS="\t"}{gsub("_","\t",$1)}1' | cut -f 1,2,3,4,5,6,8 | awk -F " " '{print $1"_"$2"_"$3"_"$4"_"$5"_"$6" "$7}' | sort -k 1n >04ICEs/Zeev_ARCs-ICE_list.txt
cat 04ICEs/Zeev_ARCs-ICE_list.txt | datamash -sW -g1 collapse 2 >04ICEs/Zeev_ARCs-ICE_list2.txt
#--------------------------------Zeev-05MGEs-------------------------#
makeblastdb -dbtype nucl -in 05MGEs/MGEs.fasta -input_type fasta -title MGEs -out 05MGEs/MGEs 
blastn -query Zeev_ARCs_nr.nul.fa -db 05MGEs/MGEs -out 05MGEs/Zeev_ARCs-MGEs.txt -evalue 1e-10 -perc_identity 80 -num_threads 10 -outfmt 6 -max_target_seqs 1
cut -f 1,2 05MGEs/Zeev_ARCs-MGEs.txt | awk 'BEGIN{FS=OFS="\t"}{gsub("_","\t",$1)}1' | cut -f 1,2,3,4,5,6,8 | awk -F " " '{print $1"_"$2"_"$3"_"$4"_"$5"_"$6" "$7}' | sort -k 1n >05MGEs/Zeev_ARCs-MGEs_list.txt
cat 05MGEs/Zeev_ARCs-MGEs_list.txt | datamash -sW -g1 collapse 2 >05MGEs/Zeev_ARCs-MGEs_list2.txt
#--------------------------------Zeev-06VFs-------------------------#
diamond makedb --in 06victors/victors_pro.fasta -d 06victors/victors_nr
diamond blastx -d 06victors/victors_nr -q Zeev_ARCs_nr.nul.fa -o 06victors/Zeev_ARCs-victors.txt --evalue 1e-10 --query-cover 70 --id 80 -k 1
cut -f 1,2 06victors/Zeev_ARCs-victors.txt | awk 'BEGIN{FS=OFS="\t"}{gsub("_","\t",$1)}1' | cut -f 1,2,3,4,5,6,8 | awk -F " " '{print $1"_"$2"_"$3"_"$4"_"$5"_"$6" "$7}' | sort -k 1n >06victors/Zeev_ARCs-victors-list.txt
cat 06victors/Zeev_ARCs-victors-list.txt | datamash -sW -g1 collapse 2 >06victors/Zeev_ARCs-victors-list2.txt
#--------------------------------Zeev-09CAT-------------------------#
CAT contigs -c Zeev_ARCs_nr.fasta -d /media/kunyu/data2/db/CAT_prepare_20210107/2021-01-07_CAT_database -t /media/kunyu/data2/db/CAT_prepare_20210107/2021-01-07_taxonomy -p Zeev_ARCs_nr.pro.fa
CAT add_names -i out.CAT.contig2classification.txt -o CAT_contigs_classification.txt -t /media/kunyu/data2/db/CAT_prepare_20210107/2021-01-07_taxonomy
CAT add_names -i out.CAT.ORF2LCA.txt -o CAT_ORFs_classification.txt -t /media/kunyu/data2/db/CAT_prepare_20210107/2021-01-07_taxonomy
mv out.CAT.contig2classification.txt 07CAT/
mv CAT_contigs_classification.txt 07CAT/
mv out.CAT.ORF2LCA.txt 07CAT/
mv CAT_ORFs_classification.txt 07CAT/
mv out.CAT.predicted_proteins.faa 07CAT/
mv out.CAT.predicted_proteins.gff 07CAT/
mv out.CAT.log 07CAT/
rm -rf out.CAT.alignment.diamond
#--------------------------------Zeev-08plasflow-------------------------#
conda deactivate
conda activate plasflow
PlasFlow.py --input Zeev_ARCs_nr.fasta --output 08plasflow/Zeev_ARCs_nr.plasflow.txt --threshold 0.7
#--------------------------------Zeev-08metawrap-------------------------#
conda deactivate
conda activate metawrap-env
cp Zeev_ARCs_nr.fasta 10contigs/
metawrap classify_bins -b 10contigs -o 11metawrap-classify -t 10
#--------------------------------Zeev-ARCs-coverm-------------------------#
conda deactivate
conda activate tools
cd ..
for F in 01Reads_fq/*_1.fastq
do	
	R=${F%_*}_2.fastq
	BASE=${F##*/}
	SAMPLE=${BASE%_*}
	TMPDIR=/media/kunyu/data2/03Zeev/tmp coverm make -r 07ARCs_analysis/Zeev_ARCs_nr.fasta -1 $F -2 $R -o 07ARCs_analysis/12coverm -t 30
	coverm filter -b 07ARCs_analysis/12coverm/Zeev_ARCs_nr.fasta.$BASE.bam -o 07ARCs_analysis/12coverm/$SAMPLE.bam --min-read-percent-identity 95 --min-read-aligned-percent 75 --threads 30
	coverm contig --bam-files 07ARCs_analysis/12coverm/$SAMPLE.bam -o 07ARCs_analysis/12coverm/$SAMPLE.rpkm.txt --trim-min 10 --trim-max 90 --min-read-percent-identity 95 --min-read-aligned-percent 75 -m rpkm -t 30
	coverm contig --bam-files 07ARCs_analysis/12coverm/$SAMPLE.bam -o 07ARCs_analysis/12coverm/$SAMPLE.count.txt --min-covered-fraction 0 -m count -t 30
	coverm contig --bam-files 07ARCs_analysis/12coverm/$SAMPLE.bam -o 07ARCs_analysis/12coverm/$SAMPLE.bases.txt --min-covered-fraction 0 -m covered_bases -t 30
done
cd 07ARCs_analysis/12coverm
paste *.rpkm.txt >Zeev_SSU.rpkm.txt
paste *.count.txt >Zeev_SSU.count.txt
paste *.bases.txt >Zeev_SSU.bases.txt
