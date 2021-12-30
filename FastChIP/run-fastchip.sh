#!/bin/bash  -l

output=/mnt/rstor/SOM_GENE_BEG33/ChIP_seq/hg38/DATA ## default
home=/mnt/rstor/SOM_GENE_BEG33/ChIP_seq/hg38/scripts
runchipseq=$home/run-chipseq.sh
runchipseq_hichip=$home/run-chipseq-for-hic.sh
readxlsx=$home/read-xlsx.sh


usage(){
echo "usage: $0 <input.txt> [<output>]
function: 
	<input.txt>: (sample, groupdir, spikein, pairedinput)
	run bwa mapping, MACS peak finding, Homer motif finding, ROSE tools 
	<output>: optional output directory default $output
"
}
if [ $# -lt 1 ];then usage; exit 1; fi

output=${2:-$output}
#$readxlsx | awk -v OFS="\t" -v S=$1 'S==$2' > sample.txt
cat $1 | grep -vE "(^$|^#)"  | while read -r line;do
	a=( `echo $line ` );
	if [ ${#a[@]} -ne 4 ];then 
		echo "$line error" >&2
		exit 1
	fi
	sample=${a[0]};
	group=${a[1]};
	spike=${a[2]};
	pairedinput=${a[3]};
	genome=hg38


	#fq1=/mnt/rstor/SOM_GENE_BEG33/data/$group/${sample}_R1.fastq.gz 
	fq1=`ls /mnt/rstor/SOM_GENE_BEG33/data/$group/${sample}*_R1*.fastq.gz`

	if [ ! -f $fq1 ]; then
		echo "$fq1 not exits! now exit" >&2; exit 1;
	fi

	if [ -d $output/$sample ];then 
		echo "$output/$sample exits but overwrite!  ">&2;
	fi

	pairedinput_fq1=`ls /mnt/rstor/SOM_GENE_BEG33/data/$group/${pairedinput}*_R1*.fastq.gz 2> /dev/null`
	if [[ ${#pairedinput} -lt 3 || ! -f $pairedinput_fq1 ]];then pairedinput="No"; fi

	if [ -z ${sample##*HiChIP*} ];then
		## paired input not implemented in HiChIP
		bash $runchipseq_hichip $fq1 $output $sample $genome $spike 
	else
		bash $runchipseq $fq1 $output $sample $genome $spike $pairedinput
	fi
done
