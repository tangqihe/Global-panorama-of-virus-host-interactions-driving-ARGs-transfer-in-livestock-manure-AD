#!/bin/bash
#SBATCH --job-name=myjob_name
#SBATCH --chdir=/work/<USERNAME>
#SBATCH --output=/work/%u/%x-%j.log
#SBATCH --time=0-00:30:00
#nohup bash 05MAGs_module.sh > 05MAGs_module.log 2>&1 &
source deactivate
conda activate base
cd /media/kunyu/data5/07Compost_meta
#--------------------------------MAGs_dRep-------------------------#
conda activate dRep
for F in 01Reads_fq/*_1.fastq
do
	R=${F%_*}_2.fastq
	BASE=${F##*/}
	SAMPLE=${BASE%_*}
	cd /media/kunyu/data5/06AD_meta/06MAGs/02refine_bins/$SAMPLE/metawrap_50_20_bins
	rename 's/^/'$SAMPLE'-/g' *.fa
	cd /media/kunyu/data5/06AD_meta
	cp 06MAGs/02refine_bins/$SAMPLE/metawrap_50_20_bins/*.fa 06MAGs/03beforedRep/
done
dRep dereplicate 06MAGs/04dRep/ -g 06MAGs/03beforedRep/*.fa -comp 70 -con 10 -p 50
conda deactivate
conda activate gtdbtk
gtdbtk classify_wf --skip_ani_screen --genome_dir 06MAGs/04dRep/dereplicated_genomes/ -x fa --out_dir 06MAGs/06gtdbtk_classify/ --cpus 50
conda deactivate
conda activate metawrap-env
metawrap classify_bins -b 06MAGs/04dRep/dereplicated_genomes -o 06MAGs/05metawrap-classify -t 50
metawrap quant_bins -b 06MAGs/04dRep/dereplicated_genomes -o 06MAGs/05MAG_QUANTIFY -t 50 01Reads_fq/*fastq
#--------------------------------MAGs_ARGs-------------------------#
conda deactivate
conda activate tools
for F in 06MAGs/04dRep/dereplicated_genomes_clean/*.fa
do
	BASE=${F##*/}	
	SAMPLE=${BASE%.*}
	prodigal -i $F -o 06MAGs/ORFs_genes/$SAMPLE.genes -a 06MAGs/ORFs_pro/$SAMPLE.pro.fa -d 06MAGs/ORFs_nucl/$SAMPLE.fa -p meta
done
diamond makedb --in 06MAGs/ORFs_nucl/01SARG/SARG.2.2.fasta -d  06MAGs/ORFs_nucl/01SARG/SARG.2.2_nr
diamond makedb --in 06MAGs/ORFs_nucl/02card/CARD3.1.4.fasta -d 06MAGs/ORFs_nucl/02card/card3.1.4_nr
diamond makedb --in 06MAGs/ORFs_nucl/03bacmet/BacMet2.fasta -d 06MAGs/ORFs_nucl/03bacmet/bacmet2_nr
diamond makedb --in 06MAGs/ORFs_nucl/04ICEs/ICE3.0.fasta -d 06MAGs/ORFs_nucl/04ICEs/ICE3.0_nr
diamond makedb --in 06MAGs/ORFs_nucl/06victors/victors_pro.fasta -d 06MAGs/ORFs_nucl/06victors/victors_nr
diamond makedb --in 06MAGs/ORFs_nucl/08Ncyc/Ncyc.fasta -d  06MAGs/ORFs_nucl/08Ncyc/Ncyc_nr
diamond makedb --in 06MAGs/ORFs_nucl/09Scyc/SCycDB.fasta -d 06MAGs/ORFs_nucl/09Scyc/SCycDB_nr
diamond makedb --in 06MAGs/ORFs_nucl/10Sulfur/sulfur_database.fasta -d 06MAGs/ORFs_nucl/10Sulfur/sulfur_database_nr
diamond makedb --in 06MAGs/ORFs_nucl/11Propionate/Propionate_database.fasta -d 06MAGs/ORFs_nucl/11Propionate/Propionate_database_nr
for F in 06MAGs/ORFs_nucl/*.fa;
do
	BASE=${F##*/}	
	SAMPLE=${BASE%.*}
	diamond blastx -d 06MAGs/ORFs_nucl/01SARG/SARG.2.2_nr -q $F -o $F-SARG.txt --evalue 1e-10 --query-cover 70 --id 80 -k 1
	diamond blastx -d 06MAGs/ORFs_nucl/02card/card3.1.4_nr -q $F -o $F-card.txt --evalue 1e-5 --query-cover 70 --id 80 -k 1
	diamond blastx -d 06MAGs/ORFs_nucl/03bacmet/bacmet2_nr -q $F -o $F-bacmet2.txt --evalue 1e-10 --query-cover 70 --id 80 -k 1
	diamond blastx -d 06MAGs/ORFs_nucl/04ICEs/ICE3.0_nr -q $F -o $F-ICE.txt --evalue 1e-10 --query-cover 70 --id 80 -k 1
	diamond blastx -d 06MAGs/ORFs_nucl/06victors/victors_nr -q $F -o $F-victors.txt --evalue 1e-10 --query-cover 70 --id 80 -k 1
diamond blastx -d 06MAGs/ORFs_nucl/04methane/methane_nr -q $F -o $F-methane.txt --evalue 1e-10 --query-cover 70 --id 80 -k 1
diamond blastx -d 06MAGs/ORFs_nucl/08Ncyc/Ncyc_nr -q $F -o $F-Ncyc.txt --evalue 1e-10 --query-cover 70 --id 80 -k 1
diamond blastx -d 06MAGs/ORFs_nucl/09Scyc/SCycDB_nr -q $F -o $F-SCyc.txt --evalue 1e-10 --query-cover 70 --id 80 -k 1
diamond blastx -d 06MAGs/ORFs_nucl/10Sulfur/sulfur_database_nr -q $F -o $F-sulfur.txt --evalue 1e-10 --query-cover 70 --id 80 -k 1
diamond blastx -d 06MAGs/ORFs_nucl/11Propionate/Propionate_nr -q $F -o $F-Prop.txt --evalue 1e-10 --query-cover 70 --id 80 -k 1
diamond blastx -d /media/kunyu/data2/db/dbcan/CAZy -q $F -o $F-CAZy.txt --evalue 1e-10 --query-cover 70 --id 80 -k 1
	mv $F-methane.txt 06MAGs/ORFs_nucl/04methane/
	mv $F-SARG.txt 06MAGs/ORFs_nucl/01SARG/
	mv $F-card.txt 06MAGs/ORFs_nucl/02card/
	mv $F-bacmet2.txt 06MAGs/ORFs_nucl/03bacmet/
	mv $F-ICE.txt 06MAGs/ORFs_nucl/04ICEs/
	mv $F-victors.txt 06MAGs/ORFs_nucl/06victors/
	mv $F-Ncyc.txt 06MAGs/ORFs_nucl/08Ncyc/
	mv $F-SCyc.txt 06MAGs/ORFs_nucl/09Scyc/
	mv $F-sulfur.txt 06MAGs/ORFs_nucl/10Sulfur/
	mv $F-Prop.txt 06MAGs/ORFs_nucl/11Propionate/
	mv $F-CAZy.txt 06MAGs/ORFs_nucl/13CAZy/
done
makeblastdb -dbtype nucl -in 06MAGs/ORFs_nucl/05MGEs/MGEs.fasta -input_type fasta -title MGEs -out 06MAGs/ORFs_nucl/05MGEs/MGEs 
for F in 06MAGs/ORFs_nucl/*.fa;
do
	BASE=${F##*/}	
	SAMPLE=${BASE%.*}
	blastn -query $F -db 06MAGs/ORFs_nucl/05MGEs/MGEs -out $F-MGEs.txt -evalue 1e-10 -perc_identity 80 -num_threads 10 -outfmt 6 -max_target_seqs 1
	mv $F-MGEs.txt 06MAGs/ORFs_nucl/05MGEs/
done
#SARG数据处理
cd /media/kunyu/data5/07Compost_meta/06MAGs/ORFs_nucl/01SARG
mkdir 01
for file in *.fa-SARG.txt; do
     filename="${file%.fa-SARG.txt}"
    if [ -f SARG_risk_mapping-20220906.txt ]; then
        awk 'NR==FNR{a[$1]=$0;next} {print $0, a[$2]}' SARG_risk_mapping-20220906.txt "$file" > 01/"${filename}.fa-SARG.txt"; fi
 done
find 01 -type f -name "*.fa-SARG.txt" -exec bash -c 'for file; do 
    filename=$(basename "$file"); 
    column2=$(awk "{printf \"%s,\", \$2}" "$file" | sed "s/,$//"); 
    column14=$(awk "{printf \"%s,\", \$14}" "$file" | sed "s/,$//"); 
     column15=$(awk "{printf \"%s,\", \$15}" "$file" | sed "s/,$//"); 
    if [ -n "$column2" ]; then
echo "${filename}, ${column2}, ${column14}, ${column15}"; fi
done' _ {} + > 01SARG.txt
#CARD数据处理
cd ..
find 02card -type f -name "*.fa-card.txt" -exec bash -c ' for file; do
    filename=$(basename "$file")
    aro_codes=""; gene_names=""; while IFS= read -r line; do
        while IFS= read -r match; do
            if [[ $match =~ gb\|([^|]+)\|ARO:([^|]+)\|([^[:space:]]+) ]]; then
                gene_name="${BASH_REMATCH[3]}"; aro_code="${BASH_REMATCH[2]}"
                [ -n "$aro_codes" ] && aro_codes+=",";  [ -n "$gene_names" ] && gene_names+=","
                aro_codes+="ARO:${aro_code}"; gene_names+="$gene_name"; fi
        done <<< "$(echo "$line" | grep -oE "gb\|[^|]+\|ARO:[^|]+\|[^[:space:]]+")"
    done < "$file"
    [ -n "$aro_codes" ] && echo "${filename}" "${aro_codes}" "${gene_names}"
done' bash {} + > 02card/02CARD.txt
#重金属抗性基因表格处理
find 03bacmet -type f -name "*.fa-bacmet2.txt" -exec bash -c ' for file; do
    filename=$(basename "$file")
    column2=$(awk "{print \$2}" "$file" | paste -sd "," -)
    if [ -n "$column2" ]; then
        echo "${filename}, ${column2}"; fi
done' bash {} + > 03bacmet/03bacmet.txt
#ICEs表格处理
find 04ICEs -type f -name "*.fa-ICE.txt" -exec bash -c ' for file; do
    filename=$(basename "$file")
    column2=$(awk "{print \$2}" "$file" | paste -sd "," -)
    if [ -n "$column2" ]; then
        echo "${filename}, ${column2}"; fi
done' bash {} + > 04ICEs/04ICEs.txt
#MGEs表格处理
find 05MGEs -type f -name "*.fa-MGEs.txt" -exec bash -c ' for file; do
    filename=$(basename "$file")
    column2=$(awk "{print \$2}" "$file" | paste -sd "," -)
    if [ -n "$column2" ]; then
        echo "${filename}, ${column2}"; fi
done' bash {} + > 05MGEs/05MGEs.txt
#victors表格处理
find 06victors -type f -name "*.fa-victors.txt" -exec bash -c ' for file; do
    filename=$(basename "$file")
    column2=$(awk "{print \$2}" "$file" | paste -sd "," -)
    if [ -n "$column2" ]; then
        echo "${filename}, ${column2}"; fi
done' bash {} + > 06victors/06victors.txt
#SCyc数据处理
find 09Scyc -type f -name "*.fa-SCyc.txt" -exec bash -c ' for file; do
    filename=$(basename "$file")
    column2=$(awk "{print \$2}" "$file" | paste -sd "," -)
    if [ -n "$column2" ]; then
        echo "${filename}, ${column2}"; fi
done' bash {} + > 09Scyc/09Scyc.txt
#10Sulfur数据处理
find 10Sulfur -type f -name "*.fa-sulfur.txt" -exec bash -c ' for file; do
    filename=$(basename "$file")
    column2=$(awk "{print \$2}" "$file" | paste -sd "," -)
    if [ -n "$column2" ]; then
        echo "${filename}, ${column2}"; fi
done' bash {} + > 10Sulfur/10Sulfur.txt
#11Propionate数据处理
find 11Propionate -type f -name "*.fa-Prop.txt" -exec bash -c ' for file; do
    filename=$(basename "$file")
    column2=$(awk "{print \$2}" "$file" | paste -sd "," -)
    if [ -n "$column2" ]; then
        echo "${filename}, ${column2}"; fi
done' bash {} + > 11Propionate/11Propionate.txt
#methane数据处理
find 04methane -type f -name "*.fa-methane.txt" -exec bash -c ' for file; do
    filename=$(basename "$file")
    column2=$(awk "{print \$2}" "$file" | paste -sd "," -)
    if [ -n "$column2" ]; then
        echo "${filename}, ${column2}"; fi
done' bash {} + > 04methane/04methane.txt
#13CAZy数据处理
find 13CAZy -type f -name "*.fa-CAZy.txt" -exec bash -c ' for file; do
    filename=$(basename "$file")
    column2=$(awk "{print \$2}" "$file" | paste -sd "," -)
    if [ -n "$column2" ]; then
        echo "${filename}, ${column2}"; fi
done' bash {} + > 13CAZy/13CAZy.txt
#--------------------------------MAGs_RNAs-------------------------#
cd ..
cd ..
for F in 06MAGs/04dRep/dereplicated_genomes_clean/*.fa
do
	BASE=${F##*/}	
	SAMPLE=${BASE%.*}
	barrnap --outseq 06MAGs/06RNA/$SAMPLE.RNA.fa $F
done
#--------------------------------MAGs_coverM-------------------------#
cd 06MAGs/04dRep
for F in dereplicated_genomes_clean/*.fa
do
	BASE=${F##*/}
	SAMPLE=${BASE%.fa*}
 	sed "s/k141\|k119/$SAMPLE/g" $F >02MAGs_database/$SAMPLE.fa
done
cat 02MAGs_database/*.fa >/media/kunyu/data5/07Compost_meta/06MAGs/13coverm/Compost_MAGs_database.fa
cd /media/kunyu/data5/07Compost_meta
for F in 01Reads_fq/*_1.fastq
do
	R=${F%_*}_2.fastq
	BASE=${F##*/}
	SAMPLE=${BASE%_*}
	TMDIR=/media/kunyu/data5/07Compost_meta/tmp coverm make -r 06MAGs/13coverm/Compost_MAGs_database.fa -1 $F -2 $R -o 06MAGs/13coverm -t 50
	coverm filter -b 06MAGs/13coverm/Compost_MAGs_database.fa.$BASE.bam -o 06MAGs/13coverm/$SAMPLE.bam --min-read-percent-identity 95 --min-read-aligned-percent 75 --threads 50
	coverm genome --bam-files 06MAGs/13coverm/$SAMPLE.bam --genome-fasta-directory 06MAGs/04dRep/02MAGs_database/ -x .fa -o 06MAGs/13coverm/$SAMPLE.rpkm.txt --trim-min 10 --trim-max 90 --min-read-percent-identity 95 --min-read-aligned-percent 75 -m rpkm -t 50
	coverm genome --bam-files 06MAGs/13coverm/$SAMPLE.bam --genome-fasta-directory 06MAGs/04dRep/02MAGs_database/ -x .fa -o 06MAGs/13coverm/$SAMPLE.count.txt --min-covered-fraction 0 -m count -t 50
	coverm genome --bam-files 06MAGs/13coverm/$SAMPLE.bam --genome-fasta-directory 06MAGs/04dRep/02MAGs_database/ -x .fa -o 06MAGs/13coverm/$SAMPLE.bases.txt --min-covered-fraction 0 -m covered_bases -t 50
done
cd 06MAGs/13coverm
paste *.rpkm.txt >Compost_meta_MAGs.rpkm.txt
paste *.count.txt >Compost_meta_MAGs.count.txt
paste *.bases.txt >Compost_meta_MAGs.bases.txt
awk 'BEGIN { FS="\t"; OFS="\t" } { printf "%s", $1; for (i=2; i<=NF; i+=2) printf "%s%s", OFS, $i; printf "\n" }'  Compost_meta_MAGs.rpkm.txt >  00Compost_meta_MAGs.rpkm.txt
awk 'BEGIN { FS="\t"; OFS="\t" } { printf "%s", $1; for (i=2; i<=NF; i+=2) printf "%s%s", OFS, $i; printf "\n" }'  Compost_meta_MAGs.bases.txt >  00Compost_Meta_MAGs.bases.txt
awk 'BEGIN { FS="\t"; OFS="\t" } { printf "%s", $1; for (i=2; i<=NF; i+=2) printf "%s%s", OFS, $i; printf "\n" }'  Compost_meta_MAGs.count.txt >  00Compost_Meta_MAGs.count.txt
