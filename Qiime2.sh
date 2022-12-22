#! /bin/bash

folder_out=Haslea_O/16S
folder_in=Haslea_O/data
out=Haslea_O/16S/Qiime2
q_taxa=silva_132_99_16S_taxonomy_7.qza
q_fasta=silva_132_99_16S.qza
train_silva=silva-138-99-nb-classifier.qza

#IDENTIFICATION PART, qiime2 2021.2.0
mkdir -p ${out}/data
for file in $(ls ${folder_out}/Matam/data/matam_assembly/*_scaffolds.matam.fa); do
	tmp=$(basename ${file})
	name=$(echo "${tmp%%.*}")
	if [ ! -f ${out}/${name}.taxa.tsv ]; then
		qiime tools import --type 'FeatureData[Sequence]' --input-path $file --output-path ${out}/data/${name}.qza
		qiime feature-classifier classify-consensus-blast --i-query ${out}/data/${name}.qza --i-reference-reads $q_fasta --i-reference-taxonomy $q_taxa --o-classification ${name}.taxa --output-dir ${out}/tmp
		qiime tools export --input-path ${out}/${name}.taxa.qza --output-path ${out}
		mv ${out}/taxonomy.tsv ${out}/${name}.taxa.tsv
		rm -r ${out}/tmp
	fi
done

#MEAN DEPTH + TAXA in a tsv file
for tsv in $(ls ${out}/*tsv); do
	tmp=$(basename ${tsv})
	name=$(echo "${tmp%%.taxa*}")
	nameO=$(echo "${tmp%%_scaffolds*}")
	sed -i 's/ /_/g' "$tsv"
	for contig in $(egrep -v "Feature" "$tsv"| cut -f1); do
		depth=$(grep -w "^$contig" ${folder_out}/Matam/Depth/${nameO}.meandepth | cut -f2)
		echo -e $depth"\t"$(grep -w "^$contig" "$tsv" | cut -f2 | sed 's/;\s*/\t/g')>> ${folder_out}/${name}.tmp 
	done
	sed -i 's/ /;/g' ${folder_out}/${name}.tmp
	mv ${folder_out}/${name}.tmp ${folder_out}/Matam/Depth/${name}.meandepth.taxa.tsv
done
