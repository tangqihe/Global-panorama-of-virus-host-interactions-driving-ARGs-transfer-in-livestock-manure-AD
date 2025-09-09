#!/bin/bash
#SBATCH --job-name=myjob_name
#SBATCH --chdir=/work/<USERNAME>
#SBATCH --output=/work/%u/%x-%j.log
#SBATCH --time=0-00:30:00
#nohup bash 06Virus_module.sh > 06Virus_module.log 2>&1 &
source deactivate
conda activate base
cd /media/kunyu/data5/06AD_meta
#..................................contigs>5000bp....................................#
conda deactivate
conda activate tools
for F in 01Reads_fq/*_1.fastq
do
	R=${F%_*}_2.fastq
	BASE=${F##*/}
	SAMPLE=${BASE%_*}
	sed -E "s/k141|k119/$SAMPLE/g" 05ASSEMBLY/$SAMPLE/final_assembly.fasta | seqtk seq -L 5000 >06VCs/00Contigs/$SAMPLE.fa
done
cat 06VCs/00Contigs/*.fa >06VCs/00Contigs/AD_Meta_contigs5000bp.fasta
#--------------------------------Genomad-------------------------#
conda deactivate
conda activate genomad
genomad end-to-end --cleanup 06VCs/00Contigs/AD_Meta_contigs5000bp.fasta 06VCs/01genomad /media/kunyu/data2/db/genomad_db/genomad_db
cd-hit-est -i 06VCs/01genomad/AD_Meta_contigs5000bp_summary/AD_Meta_contigs5000bp_virus.fna -o 06VCs/AD_Meta_virus_nr.fasta -c 1.0 -n 10 -M 16000 -T 8
cd 06VCs
#....................................Clustergenome.......................................................#
conda deactivate
conda activate tools
makeblastdb -in AD_Meta_virus_nr.fasta -dbtype nucl -out 01Clustergenome/AD_Meta_virus_nr
blastn -query AD_Meta_virus_nr.fasta -db 01Clustergenome/AD_Meta_virus_nr -outfmt '6 std qlen slen' -max_target_seqs 10000 -out 01Clustergenome/AD_Meta_virus_nr.tsv -num_threads 32 
python 01Clustergenome/anicalc.py -i 01Clustergenome/AD_Meta_virus_nr.tsv -o 01Clustergenome/AD_Meta_virus_ani.tsv
python 01Clustergenome/aniclust.py --fna AD_Meta_virus_nr.fasta --ani 01Clustergenome/AD_Meta_virus_ani.tsv --out 01Clustergenome/AD_Meta_virus_clusters.tsv --min_ani 95 --min_tcov 85 --min_qcov 0
cut -f 1 01Clustergenome/AD_Meta_virus_clusters.tsv >01Clustergenome/AD_Meta_virus_clusters.txt
seqtk subseq AD_Meta_virus_nr.fasta 01Clustergenome/AD_Meta_virus_clusters.txt >01Clustergenome/AD_Meta_VCs.fasta
#....................................checkV.......................................................#
conda deactivate
conda activate vs2
checkv end_to_end 01Clustergenome/AD_Meta_VCs.fasta 03checkv -t 50 -d /media/kunyu/data2/db/checkv/checkv-db-v1.4
#--------------------------------Virus-coverm-------------------------#
cd ..
conda deactivate
conda activate tools
for F in 01Reads_fq/*_1.fastq
do	
	R=${F%_*}_2.fastq
	BASE=${F##*/}
	SAMPLE=${BASE%_*}
	TMDIR=/media/kunyu/data5/06AD_meta/tmp coverm make -r 06VCs/01Clustergenome/AD_Meta_VCs.fasta -1 $F -2 $R -o 06VCs/11coverM -t 30
	coverm filter -b 06VCs/11coverM/AD_Meta_VCs.fasta.$BASE.bam -o 06VCs/11coverM/$SAMPLE.bam --min-read-percent-identity 95 --min-read-aligned-percent 75 --threads 30
	coverm contig --bam-files 06VCs/11coverM/$SAMPLE.bam -o 06VCs/11coverM/$SAMPLE.rpkm.txt --trim-min 10 --trim-max 90 --min-read-percent-identity 95 --min-read-aligned-percent 75 -m rpkm -t 30
	coverm contig --bam-files 06VCs/11coverM/$SAMPLE.bam -o 06VCs/11coverM/$SAMPLE.count.txt --min-covered-fraction 0 -m count -t 30
	coverm contig --bam-files 06VCs/11coverM/$SAMPLE.bam -o 06VCs/11coverM/$SAMPLE.bases.txt --min-covered-fraction 0 -m covered_bases -t 30
done
cd 06VCs/11coverM
paste *.rpkm.txt >AD_Meta_virus.rpkm.txt
paste *.count.txt >AD_Meta_virus.count.txt
paste *.bases.txt >AD_Meta_virus.bases.txt
awk 'BEGIN { FS="\t"; OFS="\t" } { printf "%s", $1; for (i=2; i<=NF; i+=2) printf "%s%s", OFS, $i; printf "\n" }'  AD_Meta_virus.rpkm.txt >  00AD_Meta_virus.rpkm.txt
awk 'BEGIN { FS="\t"; OFS="\t" } { printf "%s", $1; for (i=2; i<=NF; i+=2) printf "%s%s", OFS, $i; printf "\n" }'  AD_Meta_virus.bases.txt >  00AD_Meta_virus.bases.txt
awk 'BEGIN { FS="\t"; OFS="\t" } { printf "%s", $1; for (i=2; i<=NF; i+=2) printf "%s%s", OFS, $i; printf "\n" }'  AD_Meta_virus.count.txt >  00AD_Meta_virus.count.txt
#....................................vibrant.......................................................#
cd ..
conda deactivate
conda activate vibrant
VIBRANT_run.py -i 01Clustergenome/AD_Meta_VCs.fasta -f nucl -t 50 -virome
#....................................ARGs.......................................................#
conda deactivate
conda activate tools
prodigal -i 01Clustergenome/AD_Meta_VCs.fasta -o 01Clustergenome/AD_Meta_VCs.genes -a 01Clustergenome/AD_Meta_VCs.faa -d 01Clustergenome/AD_Meta_VCs.fna -p meta
diamond makedb --in 06ARGs/SARG.2.2.fasta -d 06ARGs/SARG.2.2_nr
diamond blastx -d 06ARGs/SARG.2.2_nr -q 01Clustergenome/AD_Meta_VCs.fna -o 06ARGs/AD_Meta_VCs-SARG.txt --evalue 1e-10 --query-cover 70 --id 80 -k 1
diamond makedb --in 06ARGs/CARD3.1.4.fasta -d 06ARGs/card3.1.4_nr
diamond blastx -d 06ARGs/card3.1.4_nr -q 01Clustergenome/AD_Meta_VCs.fna -o 06ARGs/AD_Meta_VCs-card.txt --evalue 1e-10 --query-cover 70 --id 80 -k 1
#....................................Other.......................................................#
diamond makedb --in 06other/BacMet2.fasta -d 06other/BacMet2_nr
diamond blastx -d 06other/BacMet2_nr -q 01Clustergenome/AD_Meta_VCs.fna -o 06other/001AD_Meta_VCs-BacMet.txt --evalue 1e-10 --query-cover 70 --id 80 -k 1
diamond makedb --in 06other/Ncyc.fasta -d 06other/Ncyc_nr
diamond blastx -d 06other/Ncyc_nr -q 01Clustergenome/AD_Meta_VCs.fna -o 06other/001AD_Meta_VCs-Ncyc.txt --evalue 1e-10 --query-cover 70 --id 80 -k 1
diamond makedb --in 06other/SCycDB.fasta -d 06other/SCycDB_nr
diamond blastx -d 06other/SCycDB_nr -q 01Clustergenome/AD_Meta_VCs.fna -o 06other/001AD_Meta_VCs-SCycDB.txt --evalue 1e-10 --query-cover 70 --id 80 -k 1
diamond makedb --in 06other/victors_pro.fasta -d 06other/victors_pro_nr
diamond blastx -d 06other/victors_pro_nr -q 01Clustergenome/AD_Meta_VCs.fna -o 06other/001AD_Meta_VCs-victors.txt --evalue 1e-10 --query-cover 70 --id 80 -k 1
diamond blastx -d /media/kunyu/data2/db/dbcan/CAZy -q 01Clustergenome/AD_Meta_VCs.fna -o 001AD_Meta_VCs-CAZy.txt --evalue 1e-10 --query-cover 70 --id 80 -k 1
makeblastdb -dbtype nucl -in 06other/MGEs.fasta -input_type fasta -title MGEs -out 06other/MGEs 
blastn -query 01Clustergenome/AD_Meta_VCs.fna -db 06other/MGEs -out 001AD_Meta_VCs-MGEs.txt -evalue 1e-10 -perc_identity 80 -num_threads 10 -outfmt 6 -max_target_seqs 1
	mv 001AD_Meta_VCs-CAZy.txt  06other/
                   mv 001AD_Meta_VCs-MGEs.txt 06other/
diamond makedb --in 06other/Methane_mechanism_database.fasta -d 06other/methane_nr
diamond blastx -d 06other/methane_nr -q 01Clustergenome/AD_Meta_VCs.fna -o 06other/001AD_Meta_VCs-methane.txt --evalue 1e-10 --query-cover 70 --id 80 -k 1
#--------------------------------Genomad_classify-------------------------#
conda deactivate
conda activate genomad
genomad end-to-end --cleanup 01Clustergenome/AD_Meta_VCs.fasta 02genomad /media/kunyu/data2/db/genomad_db/genomad_db
#....................................CAT.......................................................#
conda deactivate
conda activate tools
CAT contigs -c 01Clustergenome/AD_Meta_VCs.fasta -d /media/kunyu/data2/db/CAT_prepare_20210107/2021-01-07_CAT_database -t /media/kunyu/data2/db/CAT_prepare_20210107/2021-01-07_taxonomy
CAT add_names -i out.CAT.contig2classification.txt -o CAT_contigs_classification.txt -t /media/kunyu/data2/db/CAT_prepare_20210107/2021-01-07_taxonomy
CAT add_names -i out.CAT.ORF2LCA.txt -o CAT_ORFs_classification.txt -t /media/kunyu/data2/db/CAT_prepare_20210107/2021-01-07_taxonomy
mv out.CAT.contig2classification.txt 02CAT/
mv CAT_contigs_classification.txt 02CAT/
mv out.CAT.ORF2LCA.txt 02CAT/
mv CAT_ORFs_classification.txt 02CAT/
mv out.CAT.predicted_proteins.faa 02CAT/
mv out.CAT.predicted_proteins.gff 02CAT/
mv out.CAT.log 02CAT/
rm -rf out.CAT.alignment.diamond
#....................................vcontact2.......................................................#
conda deactivate
conda activate vcontact2
seqtk seq -L 10000 01Clustergenome/AD_Meta_VCs.fasta >04vcontact2/AD_Meta_VCs_10kb.fasta
prodigal -i 04vcontact2/AD_Meta_VCs_10kb.fasta -o 04vcontact2/AD_Meta_VCs_10kb.genes -a 04vcontact2/AD_Meta_VCs_10kb.faa -p meta
vcontact2_gene2genome -p 04vcontact2/AD_Meta_VCs_10kb.faa -o 04vcontact2/AD_Meta_VCs_10kb_g2g.csv -s 'Prodigal-FAA'
vcontact2 -r 04vcontact2/AD_Meta_VCs_10kb.faa --rel-mode 'Diamond' --proteins-fp 04vcontact2/AD_Meta_VCs_10kb_g2g.csv --db 'ProkaryoticViralRefSeq201-Merged' --pcs-mode MCL --vcs-mode ClusterONE --c1-bin /home/kunyu/miniconda3/envs/vcontact2/bin/cluster_one-1.0.jar --output-dir 04vcontact2/vContactOut/
#....................................majority_rules_Refseq.......................................................#
conda deactivate
conda activate tools
diamond makedb --in /media/kunyu/data2/db/viral_refseq/viral.1.protein.faa -d 08majority_rules/viral_refseq
diamond blastp -q 01Clustergenome/AD_Meta_VCs.faa -d 08majority_rules/viral_refseq -o 08majority_rules/AD_Meta_VCs_blastp.txt --query-cover 50 --subject-cover 50 --evalue 1e-5 -k 1
csvtk filter -t -f "12>=50" 08majority_rules/AD_Meta_VCs_blastp.txt > 08majority_rules/AD_Meta_VCs_blastp_50score.txt
cut -f 2 08majority_rules/AD_Meta_VCs_blastp_50score.txt >08majority_rules/AD_Meta_VCs_blastp_accession.txt
rg -f 08majority_rules/AD_Meta_VCs_blastp_accession.txt /media/kunyu/data2/db/protein_taxID/prot.accession2taxid.FULL -n --no-line-number >08majority_rules/AD_Meta_VCs_refseq_accession2taxid.txt
cut -f 2 08majority_rules/AD_Meta_VCs_refseq_accession2taxid.txt >08majority_rules/AD_Meta_VCs_blastp_taxid.txt
taxonkit lineage 08majority_rules/AD_Meta_VCs_blastp_taxid.txt --data-dir /media/kunyu/data2/db/taxdump | taxonkit reformat -F --data-dir /media/kunyu/data2/db/taxdump | cut -f 1,3 >08majority_rules/AD_Meta_VCs_taxid_taxonomy.txt
paste 08majority_rules/AD_Meta_VCs_refseq_accession2taxid.txt 08majority_rules/AD_Meta_VCs_taxid_taxonomy.txt | cut -f 1,4 >08majority_rules/AD_Meta_VCs_accession_taxonomy.txt
awk -F'\t' 'FNR==NR{a[$1]=$2; next}; {if($2 in a) {print $0, "\t"a[$2];} else {print $0, "\tNA"}}' 08majority_rules/AD_Meta_VCs_accession_taxonomy.txt 08majority_rules/AD_Meta_VCs_blastp_50score.txt >08majority_rules/AD_Meta_VCs_blastp_viref_50_tax.txt
awk -F'\t' -v OFS='\t' '{sub(/_[0-9]+$/, "", $1); split($13, a, ";"); $13=a[1] ";" a[2] ";" a[3] ";" a[4] ";" a[5]; print $1, $13, 1}' 08majority_rules/AD_Meta_VCs_blastp_viref_50_tax.txt>08majority_rules/AD_Meta_VCs_family.txt
awk -F'\t' '{sums[$1"\t"$2] += $3} END {for (key in sums) print key"\t"sums[key]}' 08majority_rules/AD_Meta_VCs_family.txt > 08majority_rules/AD_Meta_VCs_family_stax.txt
cut -f 1,3 08majority_rules/AD_Meta_VCs_family_stax.txt | datamash -sW -g1 sum 2  >08majority_rules/AD_Meta_VCs_contigs_stax.txt
awk 'FNR==NR{a[$1]=$2; next}; {if($1 in a) {print $0, "\t"a[$1];} else {print $0, "\tNA"}}' 08majority_rules/AD_Meta_VCs_contigs_stax.txt 08majority_rules/AD_Meta_VCs_family_stax.txt >08majority_rules/AD_Meta_VCs_family_contigs.txt
#  MJ        ݱ AD_Meta_VCs_family.txt  ͳ  ÿ г  ִ       4  Ϊͳ ƴ      õ    ݱ 01_AD_Meta_VCs_family.txt  Ȼ 󽫸ñ vOTU г  ִ      ķ      
awk -F '\t' '
    {
        key = $0;         #         Ϊ ؼ   
        counts[key]++;    # ͳ  ÿ г  ֵĴ   
    }
    END {
        #           ʱ ļ 
        for (key in counts) {
            print key "\t" counts[key];   #      ʽΪ  ԭʼ       +  Ʊ   + ͳ ƽ  
        }
    }
' AD_Meta_VCs_family.txt > 01_AD_Meta_VCs_family.txt && \
awk -F '\t' '
    {
        key = $1;        #   ȡ  һ   ֶ   Ϊ ؼ   
        value = $4;      #   ȡ   ĸ  ֶ   Ϊ Ƚ ֵ
        category = $2;   #   ȡ ڶ    ֶ   Ϊ   
        
        #   ʼ     ֵΪ       
        if (!(key in max_values)) {
            max_values[key] = -999999;  #  趨һ    С ĳ ʼֵ
        }
        
        #  жϵ ǰ   Ƿ        
        if (value > max_values[key] || (value == max_values[key] && category == "Viruses")) {
            max_values[key] = value;    #        ֵ
            lines[key] = $0;            #    ±           
        }
    }
    END {
        #         Ľ         ļ 
        for (key in lines) {
            print lines[key];
        }
    }
' 01_AD_Meta_VCs_family_01 > 00_AD_Meta_VCs_family.txt
#....................................majority_rules_IMG/VR..............................................#
diamond blastp -q 01Clustergenome/AD_Meta_VCs.faa -d /media/kunyu/data2/db/IMG_VR/IMGVR_proteins -o 08majority_rules/01IMG_VR/AD_Meta_VCs_blastp_IMGVR.txt --query-cover 50 --subject-cover 50 --evalue 1e-5 -k 1
csvtk filter -t -f "12>=50" 08majority_rules/01IMG_VR/AD_Meta_VCs_blastp_IMGVR.txt> 08majority_rules/01IMG_VR/AD_Meta_VCs_blastp_IMGVR_50.txt
cut -f 1,2 08majority_rules/01IMG_VR/AD_Meta_VCs_blastp_IMGVR_50.txt | sed 's/|.*//g' > 08majority_rules/01IMG_VR/AD_Meta_VCs_IMGVR.txt
cut -f 2 08majority_rules/01IMG_VR/AD_Meta_VCs_IMGVR.txt > 08majority_rules/01IMG_VR/AD_Meta_VCs_IMGVR_uvig.txt
while read line; do
    if grep -w "$line" /media/kunyu/data2/db/IMG_VR/IMGVR_all_Sequence_information-high_confidence.tsv; then
        echo $line >> 08majority_rules/01IMG_VR/new_file2.txt
    fi
done < 08majority_rules/01IMG_VR/AD_Meta_VCs_IMGVR_uvig.txt >08majority_rules/01IMG_VR/AD_Meta_VCs_IMGVR_uvig_infor.txt
paste 08majority_rules/01IMG_VR/AD_Meta_VCs_IMGVR.txt 08majority_rules/01IMG_VR/AD_Meta_VCs_IMGVR_uvig_infor.txt >08majority_rules/01IMG_VR/AD_Meta_VCs_IMGVR_results.txt
cut -f 1,2,10,17,19 08majority_rules/01IMG_VR/AD_Meta_VCs_IMGVR_results.txt > 08majority_rules/01IMG_VR/01AD_VCs_IMGVR_results.txt
#....................................Refseq_genome.......................................................#
makeblastdb -dbtype nucl -in /media/kunyu/data2/db/viral_refseq/viral.1.1.genomic.fna -input_type fasta -title virus_genome -out 09blast2IMGVR/virus_genome -blastdb_version 5
blastn -query 01Clustergenome/AD_Meta_VCs.fasta -out 09blast2IMGVR/AD_Meta_VCs_blastnVG.txt -outfmt "6 qseqid sseqid length qlen slen qstart qend sstart send evalue pident staxids sscinames scomnames sblastnames bitscore salltitles qcovs qcovhsp stitle" -db 09blast2IMGVR/virus_genome -dust no -max_target_seqs 1 -perc_identity 95 -evalue 0.00001 -num_threads 20
#....................................IMG/VR.......................................................#
blastn -query 01Clustergenome/AD_Meta_VCs.fasta -out 09blast2IMGVR/AD_Meta_VCs_blastnIMGVR.txt -outfmt "6 qseqid sseqid length qlen slen qstart qend sstart send evalue pident staxids sscinames scomnames sblastnames bitscore salltitles qcovs qcovhsp stitle" -db /media/kunyu/data2/db/IMG_VR/IMGVR -dust no -max_target_seqs 1 -perc_identity 95 -evalue 0.00001 -num_threads 20
cut -f 2 09blast2IMGVR/AD_Meta_VCs_blastnIMGVR.txt >09blast2IMGVR/AD_Meta_IMGVR_VCblastn.txt
cat /media/kunyu/data2/db/IMG_VR/IMGVR_all_Host_information-high_confidence.tsv | csvtk grep -t -P 09blast2IMGVR/AD_Meta_IMGVR_VCblastn.txt >09blast2IMGVR/IMGVR_host_results.txt
#....................................Novelty.......................................................#
cd /media/kunyu/data5/06AD_meta/06VCs
conda deactivate
conda activate tools
seqkit seq -n 01Clustergenome/AD_Meta_VCs.faa | cut -c 1-5 > 10novelty/AD_Meta_VCs_names.txt
diamond blastp -q 01Clustergenome/AD_Meta_VCs.faa -d /media/kunyu/data2/db/IMG_VR/IMGVR_proteins -o 10novelty/01AD_Meta_VCs_IMGVR.txt --id 30 --query-cover 50 --evalue 1e-5 -k 1
diamond blastp -q 01Clustergenome/AD_Meta_VCs.faa -d 08majority_rules/viral_refseq -o 10novelty/02AD_Meta_VCs_Refseq.txt --id 30 --query-cover 50 --evalue 1e-5 -k 1
cut -f1 10novelty/01AD_Meta_VCs_IMGVR.txt 10novelty/02AD_Meta_VCs_Refseq.txt | sort -u > 10novelty/known_VCs.txt
