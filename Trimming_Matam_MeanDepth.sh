#! /bin/bash

folder_out=Haslea_O/16S
folder_in=Haslea_O/data

#TRIMMING PART, cutadapt 4.0
mkdir -p $[folder_out}/trimming
for R1 in ${folder_in}/*R1*.fastq.gz; do
	R2=$(ls ${R1//R1/R2})
	name_R1=$(basename ${R1})
	name_R2=$(basename ${R2})
	cutadapt -u 9 -U 9 -m 250 -q 20 --pair-filter=both --trim-n -o $[folder_out}/trimming/Trim_${name_R1} -p $[folder_out}/trimming/Trim_${name_R2} ${R1} ${R2} 
done

#ASSEMBLE PART, matam Apr2020
conda activate matam 
mkdir -p ${folder_out}/Matam/DB
mkdir -p ${folder_out}/Matam/data/matam_assembly
for R1 in {folder_out}/trimming/Trim*R1.fastq.gz; do
	R2=$(ls "${R1//R1/R2}")
	tmp=$(basename "${R1}")
	name=$(echo "${tmp%%_R1*}")
	if [ ! -f ${folder_out}/Matam/data/matam_assembly/${name}_scaffolds.matam.fa ]; then
		zcat "${R1}" > ${folder_out}/${name}.fastq
		zcat "${R2}" >> ${folder_out}/${name}.fastq
		echo "${name}" >> ${folder_out}/Matam/sample_in.txt
		cd ${folder_out}/Matam/data/matam_assembly
		matam/bin/index_default_ssu_rrna_db.py -d $[folder_out}/Matam/${name} --max_memory 10000
		matam/bin/matam_assembly.py -d $[folder_out}/Matam/DB/SILVA_128_SSURef_NR95 -i $[folder_out}/${name}.fastq --cpu 8 --max_memory 10000 -v
	fi
done

#MEAN DEPTH COUNT
#bwa-mem2 2.2.1
#samtools 1.10
#python 3.8.10
#pandas 1.4.1
mkdir -p ${folder_out}/Matam/Depth
for file in ${folder_out}/Matam/data/matam_assembly/*_scaffolds.matam.fa; do
	tmp=$(basename "${file}")
	name=$(echo "${tmp%%_scaffolds*}")

	if [ ! -f ${folder_out}/Matam/Depth/${name}.meandepth ]; then
		R1=$(ls {folder_out}/trimming/Trim*R1*)
		R2=$(ls {folder_out}/trimming/Trim*R2*)
		bwa-mem2 index $file
		bwa-mem2 mem -t 8 -o ${folder_out}/${name}.sam $file "$R1" "$R2"
		samtools view -u ${folder_out}/${name}.sam | samtools sort -o ${folder_out}/${name}.bam
		rm ${folder_out}/${name}.sam
		samtools depth -a ${folder_out}/${name}.bam -o ${folder_out}/${name}.depth
		samtools flagstat ${folder_out}/${name}.bam > ${folder_out}/${name}.bam.stat
		sed -i "s/DepthA=pandas.*/DepthA\=pandas.read_csv(\"\Haslea_O\/16S\/Matam\/Depth\/${name}.depth\",sep\=\'\\\t\')/" ${folder_out}/Matam/mean.py
		sed -i "s/MeanA.to.*/MeanA.to_csv(\"\Haslea_O\/16S\/Matam\/Depth\/${name}.meandepth\",sep\=\'\\\t\',header\=False)/" ${folder_out}/Matam/mean.py 	
		python3 ${folder_out}/Matam/mean.py 
	fi
done
