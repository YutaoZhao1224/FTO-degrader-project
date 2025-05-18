#!/bin/bash
IP=$1
Input=$2

for i in ./Peaks ./Peaks/${IP/.sorted.bam/}
do
	if [ ! -d $i ];then
		mkdir $i
	fi
done
samtools view -@ 10 --input-fmt-option 'filter=(flag == 99 || flag == 147 )' -o ${IP/.sorted.bam/.minus.sorted.bam} $IP &
samtools view -@ 10 --input-fmt-option 'filter=(flag == 83 || flag == 163 )' -o ${IP/.sorted.bam/.plus.sorted.bam} $IP &
samtools view -@ 10 --input-fmt-option 'filter=(flag == 99 || flag == 147 )' -o ${Input/.sorted.bam/.minus.sorted.bam} $Input &
samtools view -@ 10 --input-fmt-option 'filter=(flag == 83 || flag == 163 )' -o ${Input/.sorted.bam/.plus.sorted.bam} $Input &
# wait until all & jobs finish
#wait
# forward strand bedgraph
#bamCoverage -b $IP --binSize 10 -p 10 --normalizeUsing None --outFileFormat bedgraph --filterRNAstrand forward \
#	--normalizeUsing None -o ${IP/.sorted.bam/.forward.bedgraph} &
#bamCoverage -b $IP --binSize 10 -p 10 --normalizeUsing None --outFileFormat bedgraph --filterRNAstrand reverse \
#	--normalizeUsing None -o ${IP/.sorted.bam/.reverse.bedgraph} &
#bamCoverage -b $Input --binSize 10 -p 10 --normalizeUsing None --outFileFormat bedgraph --filterRNAstrand forward \
#	--normalizeUsing None -o ${Input/.sorted.bam/.forward.bedgraph} &
#bamCoverage -b $Input --binSize 10 -p 10 --normalizeUsing None --outFileFormat bedgraph --filterRNAstrand reverse \
#	--normalizeUsing None -o ${Input/.sorted.bam/.reverse.bedgraph} &
wait
# call peaks
macs2 callpeak -t ${IP/.sorted.bam/.plus.sorted.bam} -c ${Input/.sorted.bam/.plus.sorted.bam} --keep-dup all \
	-f BAM --extsize 100 -g 100e6 -m 4 50 --outdir ./Peaks/${IP/.sorted.bam/} --nomodel -q 0.01 -n ${IP/.sorted.bam/.plus} &
macs2 callpeak -t ${IP/.sorted.bam/.minus.sorted.bam} -c ${Input/.sorted.bam/.minus.sorted.bam} --keep-dup all \
	-f BAM --extsize 100 -g 100e6 -m 4 50 --outdir ./Peaks/${IP/.sorted.bam/} --nomodel -q 0.01 -n ${IP/.sorted.bam/.minus} &
wait

#cp ./Peaks/${IP/.sorted.bam/}/${IP/.sorted.bam/.plus}_peaks.narrowPeak ./Peaks
#cp ./Peaks/${IP/.sorted.bam/}/${IP/.sorted.bam/.minus}_peaks.narrowPeak ./Peaks

# make peak.bed
cut -f 1-5 ./Peaks/${IP/.sorted.bam/}/${IP/.sorted.bam/.plus}_peaks.narrowPeak | awk -v OFS="\t" '{print $0,"+"}' \
	> ./Peaks/${IP/.sorted.bam/}/${IP/.sorted.bam/.plus}_peaks.bed
cut -f 1-5 ./Peaks/${IP/.sorted.bam/}/${IP/.sorted.bam/.minus}_peaks.narrowPeak | awk -v OFS="\t" '{print $0,"-"}' \
	> ./Peaks/${IP/.sorted.bam/}/${IP/.sorted.bam/.minus}_peaks.bed
cat ./Peaks/${IP/.sorted.bam/}/${IP/.sorted.bam/.plus}_peaks.bed ./Peaks/${IP/.sorted.bam/}/${IP/.sorted.bam/.minus}_peaks.bed \
	> ./Peaks/${IP/.sorted.bam/}/${IP/.sorted.bam/}_peaks.bed
# intersect with transcriptome to remove background
bedtools intersect -a ./Peaks/${IP/.sorted.bam/}/${IP/.sorted.bam/}_peaks.bed -b ~/Genome/preparsed/hg38-genebody.bed -wa -s -u \
	| sort -k 1,1 -k2,2n > ./Peaks/${IP/.sorted.bam/}/${IP/.sorted.bam/}.peaks.bed
rm ./Peaks/${IP/.sorted.bam/}/${IP/.sorted.bam/.plus}_peaks.bed ./Peaks/${IP/.sorted.bam/}/${IP/.sorted.bam/.minus}_peaks.bed \
	./Peaks/${IP/.sorted.bam/}/${IP/.sorted.bam/}_peaks.bed


#rm ${IP/.sorted.bam/.plus.sorted.bam} ${IP/.sorted.bam/.minus.sorted.bam} ${Input/.sorted.bam/.plus.sorted.bam} ${Input/.sorted.bam/.minus.sorted.bam}
