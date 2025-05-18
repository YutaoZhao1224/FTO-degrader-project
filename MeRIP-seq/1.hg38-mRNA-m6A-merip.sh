#!/bin/bash
# how to use: bash 1.hg38-mRNA-m6A-merip.sh name.R1.fastq.gz name.R2.fastq.gz


R1=$1
R2=$2
for i in ./log ./mapping ./mapping/STAR_log ./probe_mapping ./mapping/bw
do
	if [ ! -d $i ];then
		mkdir $i
	fi
done
cutadapt -j 15 -m 26 --max-n=0 -e 0.15 -q 20 -m 30 --nextseq-trim=20 -O 6 --pair-filter=any \
	 -a AGATCGGAAGAGCACACGTCTG -A AGATCGGAAGAGCGTCGTGT \
	 -o ${R1/.fastq.gz/.trimmed_1.fastq.gz} -p ${R2/.fastq.gz/.trimmed_1.fastq.gz} \
	 $R1 $R2 > ./log/${R1/.R1.fastq.gz/}.log
clumpify.sh in=${R1/.fastq.gz/.trimmed_1.fastq.gz} in2=${R2/.fastq.gz/.trimmed_1.fastq.gz} \
	 out=${R1/.fastq.gz/.dedupe.fastq.gz} out2=${R2/.fastq.gz/.dedupe.fastq.gz} dedupe 2>> ./log/${R1/.R1.fastq.gz/}.log
cutadapt -j 15 -u -3 -U 3 --rename='{id}_{r2.cut_prefix} {comment}' \
	-o ${R1/.fastq.gz/.barcoded.fastq.gz} -p ${R2/.fastq.gz/.barcoded.fastq.gz} \
	 ${R1/.fastq.gz/.dedupe.fastq.gz} ${R2/.fastq.gz/.dedupe.fastq.gz} >> ./log/${R1/.R1.fastq.gz/}.log
rm ${R1/.fastq.gz/.trimmed_1.fastq.gz} ${R2/.fastq.gz/.trimmed_1.fastq.gz} ${R1/.fastq.gz/.dedupe.fastq.gz} ${R2/.fastq.gz/.dedupe.fastq.gz}
# save data to another folder
mv $1 $2 ~/data/Wen-long

# NEB m6A probe mapping
echo "m6A probe mapping:" >> ./log/${R1/.R1.fastq.gz/}.log
bowtie2 -p 15 --no-unal --end-to-end -L 16 -N 0 --un-conc-gz ${R1/R1.fastq.gz/after_probe.fastq.gz} \
	-x ~/Genome/tools_custom/Merip/NEB_merip_C-G-Luc/NEB-merip -1 ${R1/.fastq.gz/.barcoded.fastq.gz} -2 ${R2/.fastq.gz/.barcoded.fastq.gz} \
	2>> ./log/${R1/.R1.fastq.gz/}.log | samtools sort -@ 15 --input-fmt-option "filter=[NM]<=10" \
	-O BAM -o ./probe_mapping/${R1/.R1.fastq.gz/}.NEB_Probe.bam
rm ${R1/.fastq.gz/.barcoded.fastq.gz} ${R2/.fastq.gz/.barcoded.fastq.gz}
# hg38 maping
echo "hg38 mapping:" >> ./log/${R1/.R1.fastq.gz/}.log
ulimit -n 1000000
STAR --runThreadN 20 --genomeDir ~/Genome/hg38_UCSC \
       	--readFilesIn ${R1/R1.fastq.gz/after_probe.fastq.1.gz} ${R1/R1.fastq.gz/after_probe.fastq.2.gz} \
        --readFilesCommand gunzip -c --alignEndsType Local --outFilterMatchNminOverLread 0.66 \
	--outFilterMatchNmin 15 --outFilterMismatchNmax 5 --outFilterMismatchNoverLmax 0.2 \
	--outFilterMultimapNmax 10 --outSAMmultNmax -1 --outReadsUnmapped None --limitBAMsortRAM 8000000000 \
	--outSAMtype BAM Unsorted --outFileNamePrefix ./mapping/${R1/.R1.fastq.gz/}_
cat ./mapping/${R1/.R1.fastq.gz/}_Log.final.out >> ./log/${R1/.R1.fastq.gz/}.log
mv ./mapping/${R1/.R1.fastq.gz/}_Log.final.out ./mapping/${R1/.R1.fastq.gz/}_Log.out \
	    ./mapping/${R1/.R1.fastq.gz/}_Log.progress.out ./mapping/${R1/.R1.fastq.gz/}_SJ.out.tab ./mapping/STAR_log
rm ${R1/R1.fastq.gz/after_probe.fastq.1.gz} ${R1/R1.fastq.gz/after_probe.fastq.2.gz}

samtools sort -@ 10 --input-fmt-option 'filter=(flag == 99 || flag == 147 || flag == 83 || flag == 163 )' \
	    -o ./mapping/${R1/.R1.fastq.gz/}.sorted.bam ./mapping/${R1/.R1.fastq.gz/}_Aligned.out.bam
rm ./mapping/${R1/.R1.fastq.gz/}_Aligned.out.bam

samtools index ./mapping/${R1/.R1.fastq.gz/}.sorted.bam -@ 10
# make bigwig file
bamCoverage -p 20 -b ./mapping/${R1/.R1.fastq.gz/}.sorted.bam --normalizeUsing RPKM --binSize 10 -o ./mapping/bw/${R1/.R1.fastq.gz/}.bw 


### calculate efficiency of IP
echo "non-m6A C_Luc reads:" >> ./log/${R1/.R1.fastq.gz/}.log
echo $(samtools view ./probe_mapping/${R1/.R1.fastq.gz/}.NEB_Probe.bam | awk '$3=="C_Luc"' | wc -l) >> ./log/${R1/.R1.fastq.gz/}.log
echo "m6A G_Luc reads:" >> ./log/${R1/.R1.fastq.gz/}.log
echo $(samtools view ./probe_mapping/${R1/.R1.fastq.gz/}.NEB_Probe.bam | awk '$3=="G-Luc"' | wc -l) >> ./log/${R1/.R1.fastq.gz/}.log
