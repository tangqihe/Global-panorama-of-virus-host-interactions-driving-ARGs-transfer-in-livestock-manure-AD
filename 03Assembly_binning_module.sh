#!/bin/bash
#SBATCH --job-name=myjob_name
#SBATCH --chdir=/work/<USERNAME>
#SBATCH --output=/work/%u/%x-%j.log
#SBATCH --time=0-00:30:00
#nohup bash 03Assembly.sh > 03Assembly.log 2>&1 &
source deactivate
conda activate base
cd /media/kunyu/data5/07Compost_meta/00ass
conda activate metawrap-env
#-------------------------Assembly&binning------------------------------------#
for F in 01Reads_fq/*_1.fastq
do
	R=${F%_*}_2.fastq
	BASE=${F##*/}
	SAMPLE=${BASE%_*}
	metawrap assembly -1 $F -2 $R -m 900 -t 50 -o 05ASSEMBLY/$SAMPLE
	metawrap binning -o 06MAGs/01bins/$SAMPLE -t 50 -m 900 -a 05ASSEMBLY/$SAMPLE/final_assembly.fasta --metabat2 --concoct --maxbin2 01Reads_fq/$SAMPLE*.fastq --universal
	metawrap bin_refinement -o 06MAGs/02refine_bins/$SAMPLE -t 50 -m 900 -A 06MAGs/01bins/$SAMPLE/maxbin2_bins/ -B 06MAGs/01bins/$SAMPLE/metabat2_bins/ -C 06MAGs/01bins/$SAMPLE/concoct_bins/ -c 50 -x 20
        rm -rf 05ASSEMBLY/$SAMPLE/megahit
        rm -rf 06MAGs/01bins/$SAMPLE/work_files
        rm -rf 06MAGs/02refine_bins/$SAMPLE/work_files
done
