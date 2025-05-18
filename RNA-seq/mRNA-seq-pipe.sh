#!/bin/bash

R1=$1
R2=$2

for i in ./log ./mapping ./mapping/STAR_log ./ERCC-spike
do
    if [ ! -d $i ];then
        mkdir $i
    fi
done

cutadapt -j 15 -m 26 --max-n=0 -e 0.15 -q 20 --nextseq-trim=20 -O 6 -n 3 --pair-filter=any \
    -a AGATCGGAAGAGCACACGTCTG -A AGATCGGAAGAGCGTCGTGT \
    -o ${R1/.fastq.gz/.trimmed_1.fastq.gz} -p ${R2/.fastq.gz/.trimmed_1.fastq.gz} \
    $R1 $R2 > ./log/${R1/.R1.fastq.gz/}.log
#mv $R1 $R2 ~/data/HCT116/mRNA-seq
clumpify.sh in=${R1/.fastq.gz/.trimmed_1.fastq.gz} in2=${R2/.fastq.gz/.trimmed_1.fastq.gz} \
    out=${R1/.fastq.gz/.dedupe.fastq.gz} out2=${R2/.fastq.gz/.dedupe.fastq.gz} dedupe 2>> ./log/${R1/.R1.fastq.gz/}.log
rm ${R1/.fastq.gz/.trimmed_1.fastq.gz} ${R2/.fastq.gz/.trimmed_1.fastq.gz}

echo "ERCC-spike_in mapping:" >> ./log/${R1/.R1.fastq.gz/}.log
bowtie2 -p 15 --no-unal --end-to-end -L 16 -N 1 --mp 5 \
    --un-conc-gz ${R1/R1.fastq.gz/after_spike_in.fastq.gz} \
    -x ~/Genome/ERCC_spike_in/ERCC92 -1 ${R1/.fastq.gz/.dedupe.fastq.gz} \
    -2 ${R2/.fastq.gz/.dedupe.fastq.gz} > ./ERCC-spike/${R1/.R1.fastq.gz/}.sam 2>> ./log/${R1/.R1.fastq.gz/}.log
samtools sort -O BAM -@ 6 -o ./ERCC-spike/${R1/.R1.fastq.gz/}.sorted.bam ./ERCC-spike/${R1/.R1.fastq.gz/}.sam
rm ./ERCC-spike/${R1/.R1.fastq.gz/}.sam
rm ${R1/.fastq.gz/.dedupe.fastq.gz} ${R2/.fastq.gz/.dedupe.fastq.gz}

echo "hg38 mapping:" >> ./log/${R1/.R1.fastq.gz/}.log
ulimit -n 1000000
STAR --runThreadN 20 --genomeDir ~/Genome/hg38_UCSC \
    --readFilesIn ${R1/R1.fastq.gz/after_spike_in.fastq.1.gz} ${R1/R1.fastq.gz/after_spike_in.fastq.2.gz} \
    --readFilesCommand gunzip -c --alignEndsType Local --outFilterMatchNminOverLread 0.66 \
    --outFilterMatchNmin 15 --outFilterMismatchNmax 5 --outFilterMismatchNoverLmax 0.2 \
    --outFilterMultimapNmax 10 --outSAMmultNmax -1 --outReadsUnmapped None \
    --limitBAMsortRAM 8000000000 --outSAMtype BAM Unsorted --outFileNamePrefix ./mapping/${R1/.R1.fastq.gz/}_
cat ./mapping/${R1/.R1.fastq.gz/}_Log.final.out >> ./log/${R1/.R1.fastq.gz/}.log
mv ./mapping/${R1/.R1.fastq.gz/}_Log.final.out ./mapping/${R1/.R1.fastq.gz/}_Log.out \
    ./mapping/${R1/.R1.fastq.gz/}_Log.progress.out ./mapping/${R1/.R1.fastq.gz/}_SJ.out.tab ./mapping/STAR_log
rm ${R1/R1.fastq.gz/after_spike_in.fastq.1.gz} ${R1/R1.fastq.gz/after_spike_in.fastq.2.gz}

samtools sort -@ 10 --input-fmt-option 'filter=(flag == 99 || flag == 147 || flag == 83 || flag == 163 )' \
    -o ./mapping/${R1/.R1.fastq.gz/}.sorted.bam ./mapping/${R1/.R1.fastq.gz/}_Aligned.out.bam

rm ./mapping/${R1/.R1.fastq.gz/}_Aligned.out.bam


