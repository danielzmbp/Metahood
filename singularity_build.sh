#!/bin/bash

container_list=("pythonenv" "bandage" "bedtools" "blast" "bwasamtools" "cat" "concoct" "desman" "diamond" "drep" "fastp" "fasttree" "gtdbtk" "ip" "kofamscan" "krakenuniq" "mafft" "megahit" "metabat2" "multiqc" "prodigal" "trimal")

mkdir -p builds

for item in ${container_list[@]}; do
  sudo singularity build builds/${item}.sif container_recipes/singularity/Singularity.${item}
done
