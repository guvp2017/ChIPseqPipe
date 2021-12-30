#!/usr/bin/env bash

module load bedtools
module load gcc
module load R/4.1.1
GENOME="hg38"
HOME="/mnt/rstor/SOM_GENE_BEG33"
DATA_HOME="${HOME}/ChIP_seq/${GENOME}/DATA"
SCRIPT_HOME="${HOME}/ChIP_seq/${GENOME}/scripts/multisample_analysis"

#runBedCovCompare Usage:
# batch submission recommended: sbatch --ntasks=4 --mem=32g --cpus-per-task=32 --gres=lscratch:200 -J bCComp
#
#options required:
# -sample1 		Sample 1 name, same as sample directory name under ChIP_seq/DATA/*/
# -name1		Short name for the 1st sample
# -sample2 		Sample 2 name, same as sample directory name under ChIP_seq/DATA/*/
# -name2		Short name for the 2nd sample
# -pValue  		MACS pValue, set to be the same across both samples being compared (example, p-7)
# -overlapname	name of the genomic feature being overlapped
# -overlapbed	bed filename (file location optional if in current directory) for genomic feature to be used for overlapping
#
#optional options
# -bedsub		bed that replaces merged bed
# -annotate		bed file with a 4th column for annotation of a feature when overlapping (ie, a gene name)
#
#future options for more flexibility:
# -mergeparameters (-d 1000) #currently hardcoded
# -spiked	#currently spike in mode exlusive

	output_dir=${name1}_${name2}_${pValue}_BCC_spike; echo output_dir: $output_dir
	mkdir $output_dir

### 1. MERGE beds
	if [ -f "$bedsub" ]
		then
			echo "Substituting in an external bedfile, not using MACS out from samples."
			bedcat=$bedsub; echo "bedcat: $bedcat"
			bedsort=${output_dir}/${bedsub%.*}.sort.bed; echo "bedsort: $bedsort"
			bedmerge=${output_dir}/${bedsub%.*}.merge.bed; echo "bedmerge: $bedmerge"
			
		else
		#build file paths to BEDS
		bed1path=${DATA_HOME}/${sample1}/MACS_${pValue}/${sample1}.clean.bed; echo "bed1path: $bed1path"
		bed2path=${DATA_HOME}/${sample2}/MACS_${pValue}/${sample2}.clean.bed; echo "bed2path: $bed2path"
		bedcat=${output_dir}/${name1}_${name2}.cat.bed
		bedsort=${output_dir}/${name1}_${name2}.sort.bed
		bedmergetoomanycolumns=${output_dir}/${name1}_${name2}.mergeallcol.bed
		bedmerge=${output_dir}/${name1}_${name2}.merge.bed
		
		#cat (concatenate beds, then sort them)
		echo "making bedcat: $bedcat"
		cat $bed1path $bed2path > $bedcat
	fi
	#sort
	echo "sorting bedcat to bedsort: $bedsort"
	sort -k1,1 -k2,2n $bedcat > $bedsort
	#merge using bedtools
	mergeparameters='-d 1000'  #this is kept simple for now so that the bed only has 3 columns going into bedCov steps
	echo "running bedtools merge with parameters: $mergeparameters"
	bedtools merge $mergeparameters -i $bedsort > $bedmergetoomanycolumns
	awk -F "\t" 'BEGIN{OFS="\t"}{print $1,$2,$3}' $bedmergetoomanycolumns > $bedmerge
	
### 2. Calculate BAM coverage in BEDS
#runBedCov
	#setup variables
	bam_1=${DATA_HOME}/${sample1}/${sample1}.bam
	bam_2=${DATA_HOME}/${sample2}/${sample2}.bam
	bed_file=$bedmerge
	out_file_count=${output_dir}/${name1}_${name2}.raw.bed
	out_file_spikescaled=${output_dir}/${name1}_${name2}.RRPM.bed

	if [ -f "$out_file_spikescaled" ]
	then
		echo "Coverage completed previously, skipping ahead to plots"
	else
		#run multicov!
		echo "bedtools multicov -bams $bam_1 $bam_2 -bed $bed_file > $out_file_count"
		bedtools multicov -bams $bam_1 $bam_2 -bed $bed_file > $out_file_count
		echo "raw read coverage mapped, now calculating $out_file_spikescaled"
		
		#fetch spike in dm3 or ecoli read counts for normlization (ChIP-Rx method)
		bam_path1=`dirname $bam_1` 
		bam_path2=`dirname $bam_2`

		spike_reads1=`cat $bam_path1/SpikeIn/spike_map_summary | sed -n '2p' | cut -f 3`
		spike_reads2=`cat $bam_path2/SpikeIn/spike_map_summary | sed -n '2p' | cut -f 3`

		#echos for debugging:
		echo $bam_1
		echo `cat $bam_path1/SpikeIn/spike_map_summary | sed -n '2p' | cut -f 3`
		echo spike_reads1
		echo $bam_2
		echo `cat $bam_path2/SpikeIn/spike_map_summary | sed -n '2p' | cut -f 3`
		echo spike_reads2

		awk -F "\t" -v spike_reads1="$spike_reads1" -v spike_reads2="$spike_reads2" 'BEGIN{OFS="\t"}{print $1,$2,$3,$4*1000000/spike_reads1,$5*1000000/spike_reads2}' $out_file_count > $out_file_spikescaled
		echo "...bedCov step complete!"
	fi
		
### 3. Get an overlapped feature
#make a bedtools count to split the cov file into 2 groups (SE or not, enhancer or promoter, etc.)

	echo "overlapping merge.bed with genomic locations of $overlapname"
	echo "overlapbed: $overlapbed"
	bedoverlap=${output_dir}/${name1}_${name2}.${overlapname}.bed
	echo "bedtools intersect -wa -c -a $bed_file -b $overlapbed > $bedoverlap"
	bedtools intersect -wa -c -a $bed_file -b $overlapbed > $bedoverlap
	
### Split and stitch RRPM bed!
	echo "split by overlap then stich the overlapping"
	bedoutside=${output_dir}/${name1}_${name2}.outside_${overlapname}.RRPM.bed
	bedinside=${output_dir}/${name1}_${name2}.inside_${overlapname}.RRPM.bed
	
	bedtools intersect -wa -v -a $out_file_spikescaled -b $overlapbed > $bedoutside
	bedtools intersect -wa -a $out_file_spikescaled -b $overlapbed > $bedinside
	
	bedoutsidestitch=${output_dir}/${name1}_${name2}.${overlapname}.out.stitch.bed
	stitchparameters='-d 2000 -c 4,5 -o sum,sum'  
	echo "Outside peaks, running bedtools merge with parameters: $stitchparameters"
	bedtools merge $stitchparameters -i $bedoutside > $bedoutsidestitch
	
	bedoverlapstitch=${output_dir}/${name1}_${name2}.${overlapname}.in.stitch.bed
	stitchparameters='-d 2000 -c 4,5 -o sum,sum'  
	echo "Inside peaks, running bedtools merge with parameters: $stitchparameters"
	bedtools merge $stitchparameters -i $bedinside > $bedoverlapstitch
	
	bedout_and_institch=${output_dir}/${name1}_${name2}.${overlapname}.outandinstitch.RRPM.bed
	bedsortinstitch=${output_dir}/${name1}_${name2}.${overlapname}.sortoutandinstitch.RRPM.bed
	
	echo "combining outside and inside stitched beds"
	cat $bedoutsidestitch $bedoverlapstitch > $bedout_and_institch
	#sort
	echo "sorting bedcat to bedsort: $bedsortinstitch"
	sort -k1,1 -k2,2n $bedout_and_institch > $bedsortinstitch
	
### 4. Run R script for making plots
	
	if [ -f "$annotate" ]
		then
			echo "annotating stitched bed"
			bedsortinstitch_anot=${output_dir}/${name1}_${name2}.${overlapname}.stitch.anot.RRPM.bed
			bedtools intersect -wa -loj -a $bedsortinstitch -b $annotate > $bedsortinstitch_anot
			rm $bedsortinstitch
			echo "Rscript $SCRIPT_HOME/plotBedCovComp_AnotFuse.R $out_file_spikescaled $bedoverlap $overlapname $output_dir $name1 $name2 $bedsortinstitch_anot"
			Rscript $SCRIPT_HOME/plotBedCovComp_AnotFuse.R $out_file_spikescaled $bedoverlap $overlapname $output_dir $name1 $name2 $bedsortinstitch_anot
		else
			echo "Rscript $SCRIPT_HOME/plotBedCovComp.R $out_file_spikescaled $bedoverlap $overlapname $output_dir $name1 $name2 "
			Rscript $SCRIPT_HOME/plotBedCovComp.R $out_file_spikescaled $bedoverlap $overlapname $output_dir $name1 $name2 
	fi

### 5. Cleanup, update permissions
	
	#clean up
	rm $bedout_and_institch $bedoutside $bedinside $bedoverlapstitch $bedoutsidestitch
	echo "all done, opening permissions for files in $output_dir"
	chgrp beg33 -R $output_dir