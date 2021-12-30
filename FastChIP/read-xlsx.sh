#!/bin/bash

usage(){
echo "
usage $0 <filter> [<xlsx.file>]
	<filter> : group name (default no)
"
}
if [[ $# > 0 && $1 == "-h" ]];then usage; exit 1; fi
filter=${1:-no}
input=${2:-/mnt/rstor/SOM_GENE_BEG33/ChIP_seq/hg38/manage_samples/ChIP_seq_samples.xlsx}

module load gcc
module load R/4.1.1

## save submit files in a temp file
cat<<-'EOF' |  R --no-save --vanilla --args $input $filter |& tail -n+2 | grep "^>@" | cut -c 3- 
library("readxl")
library("dplyr")
args = commandArgs(TRUE)
input = args[1];
filter = args[2];
output = args[3];

print(paste0("reading ",input," | filter ",filter))
ChIP_samples = read_excel(input)
ChIP_samples$SpikeIn[ChIP_samples$SpikeIn==""]<-"No"
if( filter == "no" ){
	res=ChIP_samples %>%
        mutate( str=paste(SampleFiles,SequencingRun_GEO,SpikeIn,PairedInput,sep="\t"))
}else{
	res=ChIP_samples %>%
        filter( SequencingRun_GEO == !!filter ) %>%
        mutate( str=paste(SampleFiles,SequencingRun_GEO,SpikeIn,PairedInput,sep="\t"))
}
cat(paste0(">@",res$str),sep="\n")
EOF

