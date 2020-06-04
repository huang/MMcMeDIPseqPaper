#!/bin/bash



#--1 trimming and quality control --
mkdir trimmed
for sample in MMC_2pos MMC_4pos MMC_5pos MMC_2neg MMC_4neg MMC_5neg  MMC_2pos_input MMC_4pos_input MMC_5pos_input MMC_2neg_input MMC_4neg_input MMC_5neg_input; do    
    java -jar /home/jhuang/Tools/Trimmomatic-0.36/trimmomatic-0.36.jar SE -threads 16 raw_data/${sample}_R1.fastq.gz trimmed/${sample}_R1.fastq.gz ILLUMINACLIP:/home/jhuang/Tools/Trimmomatic-0.36/adapters/TruSeq3-SE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:35 AVGQUAL:20 HEADCROP:3 CROP:80;
    fastqc trimmed/${sample}_R1.fastq.gz;
done


#--2 filtered to include only de-duplicated uniquely mapped reads resulting from sam files (under host env) --
mkdir alg
#fetchChromSizes mm10 > mm10.chrom.sizes
for sample in MMC_2pos MMC_4pos MMC_5pos MMC_2neg MMC_4neg MMC_5neg  MMC_2pos_input MMC_4pos_input MMC_5pos_input MMC_2neg_input MMC_4neg_input MMC_5neg_input; do
#--2.1 alignment--
bwa mem -M -t 12 genome_mm10/mouse_genome.fa trimmed/${sample}_R1.fastq.gz | samtools view -h -bS - > alg/${sample}.bam        

#--2.2 sort--
samtools sort alg/${sample}.bam alg/${sample}.sorted
rm alg/${sample}.bam

#--2.3 index sorted--
samtools index alg/${sample}.sorted.bam 
samtools view -F 260 alg/${sample}.sorted.bam > alg/${sample}_filter.sam  #mapped reads

#--2.4 convert to BED--
#Usage: sam2bedgff.pl path2sam sample species read_length chrom_list single_end? bed? path2output
#to be run on filtered sam files. convert to bed/gff with only UNIQUE READS INCLUDED
./sam2bedgff.pl alg/${sample}_filter.sam ${sample} Mouse 84 mm10.chrom.sizes 1 0 alg/${sample}.sam   #add headers to the sam and convert to bam, using them as the input of step4.
./sam2bedgff.pl alg/${sample}_filter.sam ${sample} Mouse 84 mm10.chrom.sizes 1 1 alg/${sample}.bed   
done


#--3 transfer bed to wig and bw (under host env) --
for sample in MMC_2pos MMC_4pos MMC_5pos MMC_2neg MMC_4neg MMC_5neg  MMC_2pos_input MMC_4pos_input MMC_5pos_input MMC_2neg_input MMC_4neg_input MMC_5neg_input; do
#medips_qc.R $NAME $SPECIES $REFNAME $WINDOW_SIZE 0 $PATH2OUTPUT
./medips_qc.R --args ${sample} Mouse mm10 800 0 /mnt/Seagate_Corona/Data_Nicole_MeDIP/medip_source_code/alg/     #-->from ${sample}.bed --> ${sample}_rpkm.wig
#grep '  +' alg/${sample}.bed > alg/${sample}_for.bed
#./medips_qc.R --args ${sample}_for Mouse mm10 800 0 /mnt/Seagate_Corona/Data_Nicole_MeDIP/medip_source_code/alg/
#grep '  -' alg/${sample}.bed > alg/${sample}_rev.bed
#./medips_qc.R --args ${sample}_rev Mouse mm10 800 0 /mnt/Seagate_Corona/Data_Nicole_MeDIP/medip_source_code/alg/
done
conda activate medip  #(under medip env)
for sample in MMC_2pos MMC_4pos MMC_5pos MMC_2neg MMC_4neg MMC_5neg  MMC_2pos_input MMC_4pos_input MMC_5pos_input MMC_2neg_input MMC_4neg_input MMC_5neg_input; do
wigToBigWig -clip alg/${sample}_rpkm.wig mm10.chrom.sizes alg/${sample}_rpkm.bw
rm alg/${sample}_filter.sam
#rm alg/${sample}_for.bed
#rm alg/${sample}_rev.bed
done



#--4 call diff reps --
diffReps.pl --treatment alg/MMC_2neg.bed alg/MMC_4neg.bed alg/MMC_5neg.bed --btr alg/MMC_2neg_input.bed alg/MMC_4neg_input.bed alg/MMC_5neg_input.bed  --control alg/MMC_2pos.bed alg/MMC_4pos.bed alg/MMC_5pos.bed --bco alg/MMC_2pos_input.bed alg/MMC_4pos_input.bed alg/MMC_5pos_input.bed    --nsd broad --gname mm10 --report MMC_neg_vs_pos.mm10_800_0.0001.diff.gt.txt --meth gt --window 800 --frag 250 --nproc 14 --pval 0.0001;
for sample in 800; do
grep "Up" MMC_neg_vs_pos.mm10_${sample}_0.0001.diff.gt.txt > ${sample}.up.gt.txt;
grep "Down" MMC_neg_vs_pos.mm10_${sample}_0.0001.diff.gt.txt > ${sample}.down.gt.txt;
awk 'BEGIN {OFS="\t"} {print $1,$2,$3,"diffreps_"NR,$12}' ${sample}.up.gt.txt > ${sample}.up.gt.bed;
awk 'BEGIN {OFS="\t"} {print $1,$2,$3,"diffreps_"NR,$12}' ${sample}.down.gt.txt > ${sample}.down.gt.bed;
cat ${sample}.up.gt.bed ${sample}.down.gt.bed > ${sample}.diff.gt.bed;
cut f1-3 -d$'\t' ${sample}.diff.gt.bed > f1-3.txt;
cut f5-5 -d$'\t' ${sample}.diff.gt.bed > f5.txt;
paste -d$'\t' f1-3 f5.txt > f1-3_5.bed;
rm ${sample}.up.gt.txt;
rm ${sample}.down.gt.txt;
rm ${sample}.up.gt.bed;
rm ${sample}.down.gt.bed;
done


#--5 prepare bed files using 100k bin size for circos plot drawing --
#-- ${sample}.sam --> ${sample}.bam --> ${sample}_100k.bedgraph --> ${sample}_100k.bed --
for sample in MMC_2pos MMC_4pos MMC_5pos MMC_2neg MMC_4neg MMC_5neg  MMC_2pos_input MMC_4pos_input MMC_5pos_input MMC_2neg_input MMC_4neg_input MMC_5neg_input; do
samtools view -H alg/${sample}.sorted.bam > alg/${sample}.header
cat alg/${sample}.sam >> alg/${sample}.header
mv alg/${sample}.header alg/${sample}.sam
samtools view -h -bS alg/${sample}.sam > alg/${sample}.bam
bamCoverage --bam alg/${sample}.bam -o alg/${sample}_100k.bedgraph --binSize 100000 --normalizeUsing RPKM --effectiveGenomeSize 2725521391 --outFileFormat bedgraph
grep -v -P "\\t0\\n" alg/${sample}_100k.bedgraph | cut -f4-4 - | sort -g > ${sample}_100k.bed
rm alg/${sample}.sam
done


#--6 drawing circos plot using cicos.conf --
circos -conf circos.conf


#--7 merge replicates together --
bigWigMerge MMC_2neg_input_rpkm.bw MMC_4neg_input_rpkm.bw MMC_5neg_input_rpkm.bw MMC_neg_input.bedGraph
sort -k1,1 -k2,2n MMC_neg_input.bedGraph > MMC_neg_input.sorted.bedGraph
bedGraphToBigWig MMC_neg_input.sorted.bedGraph mm10.chrom.sizes MMC_neg_input.bw
bigWigMerge MMC_2pos_input_rpkm.bw MMC_4pos_input_rpkm.bw MMC_5pos_input_rpkm.bw MMC_pos_input.bedGraph
sort -k1,1 -k2,2n MMC_pos_input.bedGraph > MMC_pos_input.sorted.bedGraph
bedGraphToBigWig MMC_pos_input.sorted.bedGraph mm10.chrom.sizes MMC_pos_input.bw
bigWigMerge MMC_2neg_rpkm.bw MMC_4neg_rpkm.bw MMC_5neg_rpkm.bw MMC_neg.bedGraph
sort -k1,1 -k2,2n MMC_neg.bedGraph > MMC_neg.sorted.bedGraph
bedGraphToBigWig MMC_neg.sorted.bedGraph mm10.chrom.sizes MMC_neg.bw
bigWigMerge MMC_2pos_rpkm.bw MMC_4pos_rpkm.bw MMC_5pos_rpkm.bw MMC_pos.bedGraph
sort -k1,1 -k2,2n MMC_pos.bedGraph > MMC_pos.sorted.bedGraph
bedGraphToBigWig MMC_pos.sorted.bedGraph mm10.chrom.sizes MMC_pos.bw



