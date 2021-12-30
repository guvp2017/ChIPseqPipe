#!/bin/bash  -l

usage(){
echo "
	usage $0 : <fastq_R1> <output_dir> <sample_name> <genome> <spikein>
	<output_dir> : root of the results
	<sample_name> : name of the results 
	<genome> : hg38 
	<spikein>: ecoli 
"
}

## handle input parameters
if [ $# -lt 4 ];then
	usage; exit 1;
fi
fq1="$1"; fq2=${fq1/_R1/_R2}
output_dir=`realpath $2`; output_dir=${output_dir%\/}
sample=$3;
genome=$4
spike=${5:-"No"}; spike=${spike/eco*/eco};
output=$output_dir/$sample; 
echo "$0  parameters
	output : $output
	sample : $sample
	genome : $genome
	spike : $spike
"
#if [ -d $output ];then echo "$output exists">&2; exit -1;	fi
echo "generate $output .. "
mkdir -p $output


## reference files
workdir=/mnt/rstor/SOM_GENE_BEG33/ChIP_seq/hg38/lscratch
bi=/mnt/rstor/SOM_GENE_BEG33/ref_files/hg38/bwa_indices/GRChg38.d1.vd1.fa
gn=/mnt/rstor/SOM_GENE_BEG33/ref_files/hg38/bwa_indices/GRChg38.genome
pc=/mnt/rstor/SOM_GENE_BEG33/software/Picard
igvtools=/mnt/rstor/SOM_GENE_BEG33/software/IGVTools_2.3.98
rose=/mnt/rstor/SOM_GENE_BEG33/software/ROSE 
bi_ecoli=/mnt/rstor/SOM_GENE_BEG33/ref_files/Escherichia_coli_K_12_MG1655/Sequence/BWAIndex/genome.fa
bi_dm3=/mnt/rstor/SOM_GENE_BEG33/ref_files/dm3/dm3.fa
bi_mm10=/mnt/rstor/SOM_GENE_BEG33/ref_files/mm10/bwa_indices/mm10.fa
chipseq_spikein=/mnt/rstor/SOM_GENE_BEG33/ChIP_seq/hg38/scripts/chipseq-spikein.sh

bi_spike="No"
if [ $spike == "eco" ];then
	bi_spike=$bi_ecoli
elif [ $spike == "mm10" ];then
	bi_spike=$bi_mm10
fi

## slurm paramegers
mem=64g
nproc=12 ; ## multi threads

## STEP1 : fastq QC and mapping 
#sbatch<<-EOF 
jb1=`sbatch<<-EOF | egrep -o -e "\b[0-9]+$" 
	#!/bin/bash -l
	#SBATCH -t 10:00:00
	#SBATCH -p smp
	#SBATCH -n 1
	#SBATCH -o $output/slurm-%j.out
	#SBATCH -e $output/slurm-%j.err
	#SBATCH -c $nproc
	#SBATCH --mem=$mem
	#SBATCH -J fastchip-mapping

	module load picard/2.11
	module load bwa/0.7.17
	module load samtools/1.8
	module load fastqc/0.11.9
	module load base ##perl
	#module load gcc/6.3.0  ## inactivate samtools
	#module load R

	## skip Fastqc for now
	# mkdir -p $output/Fastqc
	#fastqc -t $nproc $output/sub50m_R1.fastq -o $output/Fastqc
	if [ -f $output/sub50m.fastq ];then 
		echo "$output/sub50m.fastq exists. skip this step" >&2
	else
		gunzip -dc $fq1 | head -n 50000000 > $output/sub50m.fastq
		if [ -f $fq2 ];then 
			## skip fastqc and add fastq entries
			#fastqc -t $nproc $fq2 -o $output/Fastqc
			gunzip -dc $fq2 | head -n 50000000 >> $output/sub50m.fastq
		fi
	fi
	#unzip $output/Fastqc/${sample}*.zip -d $output/Fastqc

	## skip paring 
	if [ -f $output/${sample}.bam ];then 
		echo "$output/${sample}.bam exists. skip this step." >&2	
	else
		bwa mem -5SPM -t $nproc $bi $output/sub50m.fastq  | samtools view -h -q 30  > $output/${sample}.sam 
		java -Xmx60g  -Djava.io.tmpdir=$workdir -jar $pc/picard.jar SortSam -INPUT $output/${sample}.sam -OUTPUT $output/${sample}.bam -SORT_ORDER coordinate
		rm $output/${sample}.sam
		samtools index $output/${sample}.bam
		samtools flagstat $output/${sample}.bam > $output/${sample}.flagstat.txt
		java -Xmx$mem -jar $igvtools/igvtools.jar count $output/${sample}.bam $output/${sample}.tdf $bi
	fi

	if [[ -f $bi_spike  ]];then

		if [[ $spike == "No" || -f $output/${sample}_${spike}.bam ]];then
			echo "$output/${sample}_${spike}.bam exitst. skip this step" >&2
		else 
			bwa mem -5SPM -t $nproc $bi_spike $output/sub50m.fastq  | samtools view -h -q 30 > $output/${sample}_${spike}.sam
			java -Xmx60g  -Djava.io.tmpdir=$workdir -jar $pc/picard.jar SortSam -INPUT $output/${sample}_${spike}.sam -OUTPUT  $output/${sample}_${spike}.bam -SORT_ORDER coordinate
			rm $output/${sample}_${spike}.sam
			samtools index $output/${sample}_${spike}.bam
			samtools flagstat $output/${sample}_${spike}.bam > $output/${sample}_${spike}.flagstat.txt
			## put spikein function
		fi
		if [ -f  $output/SpikeIn/spike_map_summary ];then
			echo "$output/SpikeIn/spike_map_summary exists. skip this step" >&2
		else
			$chipseq_spikein $output $sample $genome $spike
		fi
	fi
EOF`

BLACKLIST=/mnt/rstor/SOM_GENE_BEG33/ChIP_seq/hg38/ref/Blacklist_ENCFF356LFX.bed 
## treat as single-reads 
mac3f="BAM"; #if [ -f $fq2 ];then mac3f="BAMPE"; fi
jb2="";
for p in 4 5 7 14; do
	pstr=`seq $p | xargs -I {} echo -n 0`;
	out="$output/MACS_p-$p"
	mkdir -p $out/ROSE/gff/ $out/ROSE/mappedGFF/
	if [ -f $out/${sample}.clean.bed ]; then
		echo " $out/${sample}.clean.bed exists! skip this step" >&2
		continue;
	fi
	sbatch<<-EOF 
		#!/bin/bash -l
		#SBATCH -t 24:00:00
		#SBATCH -n 1 
		#SBATCH -p smp
		#SBATCH -c $nproc
		#SBATCH -o $output/slurm-%j.out
		#SBATCH -e $output/slurm-%j.err
		#SBATCH --mem=$mem
		#SBATCH -J peakcall
		#SBATCH --dependency=afterok:$jb1

		module load miniconda3/4.9.2
		source activate "/home/sxg1131/.conda/envs/py368"
		module load gcc
		module load bedtools
		module load R/4.0.2
		module load base #perl

		macs3 callpeak -t $output/${sample}.bam -f $mac3f -g hs -n ${sample} -B -p 0.${pstr}1  --outdir $output/MACS_p-$p 
		grep -vE "(chrUn|chrM|HIV|random|CMV)" $out/${sample}_peaks.narrowPeak > $out/${sample}.clean.bed 
		bedtools intersect -a $out/${sample}.clean.bed -b $BLACKLIST -v > $out/${sample}.nobl.bed 
		findMotifsGenome.pl $out/${sample}_peaks.narrowPeak hg38 $out/${sample}_motifs
		rm $out/${sample}*.bdg 2> /dev/null
		conda deactivate

		module load miniconda3/4.9.2
		export PATH=/home/sxg1131/.conda/envs/py27/bin:$PATH # activate python 2.7.5 for ROSE2
		cd $rose
  		python ROSE_main.py -g $genome -i $out/${sample}.nobl.bed -r $output/${sample}.bam -o $out/ROSE -s 12500 -t 2500
		
	EOF
	jb2="$jb2:$jid"
done
jb2=${jb2%:}



