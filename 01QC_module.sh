# 01QC_module.sh
#!/bin/bash
#SBATCH --job-name=myjob_name
#SBATCH --chdir=/work/<USERNAME>
#SBATCH --output=/work/%u/%x-%j.log
#SBATCH --time=0-00:30:00
source deactivate
conda activate base
conda activate metawrap-env
cd /media/kunyu/data6/05Zeev
for F in 01RAW_READS/*_1.fastq; do 
	R=${F%_*}_2.fastq
	BASE=${F##*/}
	SAMPLE=${BASE%_*}
	metawrap read_qc -1 $F -2 $R -t 20 --skip-bmtagger -o 01READ_QC/$SAMPLE
	pigz $F $R
done
for i in 01READ_QC/*; do 
	b=${i#*/}
	mv ${i}/final_pure_reads_1.fastq 01Reads_fq/${b}_1.fastq
	mv ${i}/final_pure_reads_2.fastq 01Reads_fq/${b}_2.fastq
done
conda deactivate
for file in 01Reads_fq/*.fastq
do
  echo "$(basename $file .fastq): $(grep -c "^@" $file) reads" >> 01Reads_fq/Zeev_Meta_reads_number.txt
done
