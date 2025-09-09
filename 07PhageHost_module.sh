#!/bin/bash
#SBATCH --job-name=myjob_name
#SBATCH --chdir=/work/<USERNAME>
#SBATCH --output=/work/%u/%x-%j.log
#SBATCH --time=0-00:30:00
#nohup bash 07PhageHost_module.sh > 07PhageHost_module.log 2>&1 &
source deactivate
conda activate base
cd /media/kunyu/data5/07Compost_meta
#....................................MAGs_crispr.......................................................#
conda activate tools
cd 06MAGs/04dRep
for F in dereplicated_genomes_clean/*.fa
do
	BASE=${F##*/}
	SAMPLE=${BASE%.fa*}
	java -cp CRT1.2-CLI.jar crt $F $F.out
	sed -i "s/k1\(19\|41\)/$SAMPLE/g" "$F.out"
	mv $F.out 01MAGs_crispr/
done
Rscript MAGs_CRT_spacers.R
cd 01MAGs_spacers
cat *.txt >MAGs_spacers.fasta
seqtk seq -L 1 MAGs_spacers.fasta > /media/kunyu/data5/07Compost_meta/06VCs/10host-link/01crispr/MAGs_spacers_nr.fasta
#....................................virus_tRNA.......................................................#
cd /media/kunyu/data5/07Compost_meta/06VCs
aragorn -t -fasta -wa -o 10host-link/02tRNA/Compost_Meta_VCs_tRNA.fasta 01Clustergenome/Compost_VCs.fasta
sed -e '/gene/d' 10host-link/02tRNA/Compost_Meta_VCs_tRNA.fasta | sed -e '/      /d'| seqtk seq -N -L 1 >10host-link/02tRNA/Compost_Meta_VCs_tRNA_nr.fasta
sed -i 's/ c\[/c\[/g; s/)\ \[/)\[/g' 10host-link/02tRNA/Compost_Meta_VCs_tRNA_nr.fasta
#....................................MAGs_database_used_for_tRNA.......................................................#
cd /media/kunyu/data5/07Compost_meta/06MAGs/04dRep
for F in dereplicated_genomes_clean/*.fa
do
	BASE=${F##*/}
	SAMPLE=${BASE%.fa*}
 	sed -E "s/k141|k119/$SAMPLE/g" $F >02MAGs_database/$SAMPLE.fa
done
cat 02MAGs_database/*.fa >/media/kunyu/data5/07Compost_meta/06VCs/10host-link/02tRNA/Compost_Meta_MAGs_database.fa
cd /media/kunyu/data5/07Compost_meta/06VCs
makeblastdb -dbtype nucl -in 10host-link/02tRNA/Compost_Meta_MAGs_database.fa -input_type fasta -title Compost_Meta_MAGs -out 10host-link/02tRNA/Compost_Meta_MAGs -blastdb_version 5
blastn -query 10host-link/02tRNA/Compost_Meta_VCs_tRNA_nr.fasta -out 10host-link/02tRNA/Compost_Meta_MAGs_hosts_tRNA.txt -outfmt "6 qseqid sseqid length qlen slen qstart qend sstart send evalue pident staxids bitscore salltitles qcovs qcovhsp stitle" -db 10host-link/02tRNA/Compost_Meta_MAGs -dust no -perc_identity 100 -evalue 0.0001 -num_threads 20 
makeblastdb -dbtype nucl -in 01Clustergenome/Compost_VCs.fasta -input_type fasta -title Compost_VCs -out 10host-link/01crispr/Compost_VCs -blastdb_version 5
blastn -task blastn-short -query 10host-link/01crispr/MAGs_spacers_nr.fasta -db 10host-link/01crispr/Compost_VCs -perc_identity 97 -out 10host-link/01crispr/Compost_VCs_MAGs_criprcas_host.txt -outfmt "6 qseqid sseqid length qlen slen qstart qend sstart send evalue pident bitscore mismatch qcovs qcovhsp" -num_threads 20
#....................................Homology.......................................................#
cp 10host-link/02tRNA/Compost_Meta_MAGs_database.fa 10host-link/03homology/
makeblastdb -dbtype nucl -in 10host-link/03homology/Compost_Meta_MAGs_database.fa -input_type fasta -title Compost_Meta_MAGs -out 10host-link/03homology/Compost_Meta_MAGs -blastdb_version 5
blastn -query 01Clustergenome/Compost_VCs.fasta -out 10host-link/03homology/Compost_Meta_Virus_homology.txt -outfmt "6 qseqid sseqid length qlen slen qstart qend sstart send evalue pident staxids sscinames scomnames sblastnames bitscore salltitles qcovs qcovhsp stitle" -db 10host-link/03homology/Compost_Meta_MAGs -dust no -perc_identity 80 -evalue 0.00001 -num_threads 40
#...................................................................................#
cd /media/kunyu/data5/07Compost_meta/06VCs/10host-link
awk -F'\t' '{
    if ($11 >= 97 && $14 >= 90 && $15 >= 90) {
        print
    }
}' 01crispr/Compost_VCs_MAGs_criprcas_host.txt | \
awk '!seen[$1,$2]++' | \
awk 'BEGIN { FS = OFS = "\t" } { $17 = $1; $16 = $1; print }'  | \
awk 'BEGIN { FS = OFS = "\t" } { gsub(/-[^_]*_/, "_", $16); sub(/_.*/, "", $1); print }' | \
awk 'BEGIN { FS = OFS = "\t" } $2 != $16' > 00Compost_VCs_MAGs_criprcas_host.txt

awk -F'\t' '{ if ($11 >= 80 && $3 >= 1000) print }' 03homology/Compost_Meta_Virus_homology.txt | \
awk -F'\t' '{ if ($4 != 0 && ($3 / $4) >= 0.5) print }' | \
awk '!seen[$1,$2]++' | \
awk 'BEGIN {FS=OFS="\t"} { sub(/_.*/, "", $2); gsub(/-[^_]*_/, "_", $20); print }' | \
awk -F'\t' '{ if ($1 != $20 && ($2 ~ /k119/ || !seen[$1,$2]++)) { print; } }' > 00Compost_Meta_Virus_homology.txt

awk 'BEGIN { OFS = "\t"; b = "" } { if ($0 !~ /^>t/) { b = $0 } else { if (b != "") { print_b(b) }; print } } END { if (b != "") { print_b(b) } } function print_b(content) { print $0, content }' | \
grep $'\t' | \
sed 's/)\s/)/g' | \
sed -E 's/>//g' | \
awk 'BEGIN { FS = OFS = "\t" } NR == FNR { mapping[$1] = $2; next } { if ($1 in mapping) { $1 = mapping[$1] } print }' - 02tRNA/Compost_Meta_MAGs_hosts_tRNA.txt | \
awk -F'\t' '{ if ($11 >= 100 && $15 >= 100 && $16 >= 100) { print } }' | \
awk '!seen[$1,$2]++' | \
awk 'BEGIN { FS = OFS = "\t" } { gsub(/-[^_]*_/, "_", $17); sub(/_.*/, "", $2); print }' | \
awk -F'\t' '{ if ($1 != $17 && ($2 ~ /k119/ || !seen[$1,$2]++)) { print; } }' > 00Compost_Meta_MAGs_hosts_tRNA.txt
