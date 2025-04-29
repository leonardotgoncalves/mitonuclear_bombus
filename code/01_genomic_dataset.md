## 1. Illumina dataset processing and mapping

### 1.1. Read processing
```bash

for i in ./*_1.fastq.gz; do

R1=$(basename ${i})
R2=${R1%_1.fastq.gz}_2.fastq.gz
OUT=${R1%_1.fastq.gz}.fastq.gz

java -jar trimmomatic-0.39.jar PE -threads ${SLURM_NTASKS} \
  $R1 $R2 -baseout $OUT \
  ILLUMINACLIP:TruSeq3-PE-2.fa:2:30:10:8:true SLIDINGWINDOW:4:20 MINLEN:50 LEADING:3 TRAILING:3

done

```
### 1.2. Mapping Illumina reads to a high-quality assembly
```bash
FASTA=hortorum.fasta

# Index reference genome
bwa index -a bwtsw ${FASTA}
samtools faidx ${FASTA}

# Align trimmed reads to reference and sort by query name
for i in ./reads/${FASTA%.fasta}/*_1P.fastq.gz; do

R1=${i}
R2=${R1%_1P.fastq.gz}_2P.fastq.gz

TMPDIR=./

# Align and sort
bwa mem -M -t ${SLURM_NTASKS} ./ref/${FASTA} ${R1} ${R2} | \
  samtools sort @ ${SLURM_NTASKS} -O BAM -T ${TMPDIR} -o ${R1%_1P.fastq.gz}_sorted.bam

done

```
## 2. Pacbio dataset processing and mapping

### 2.1. Read processing
```bash
# Conservative approach: discards the read if an adapter is found

for i in hifi_reads/*.fastq.gz; do

SAMPLE=$(basename ${i})
echo ${SAMPLE}

cutadapt -b ATCTCTCTCAACAACAACAACGGAGGAGGAGGAAAAGAGAGAGAT \
  -b ATCTCTCTCTTTTCCTCCTCCTCCGTTGTTGTTGTTGAGAGAGAT --discard-trimmed \
  -o hifi_reads/${SAMPLE%.fastq.gz}_clean.fastq.gz ${i}

done
```
### 2.2. Mapping Pacbio reads to a high-quality assembly
```bash
# Independent runs, depending on the reference genome
# Outputs minimap2 results into a sorted bam

for i in hifi_reads/*.fq.gz; do

SAMPLE=$(basename ${i})
REFERENCE=ref/pascuorum.fasta
echo ${SAMPLE}

minimap2 -ax map-pb ${REFERENCE} ${SAMPLE%.fastq.gz}_clean.fastq.gz | \
  samtools sort -@ ${SLURM_NTASKS} -O BAM -T ./ -o ${SAMPLE%.fastq.gz}_sorted.bam

done
```
## 3. Obtaining mapping statistics
```bash
for i in ./bam/*.bam; do

# Mean read depth
mean_depth=$(samtools depth -a "${i}" | awk '{c++;s+=$3}END{if (c>0) print s/c; else print 0}')

# Breadth of coverage
breadth=$(samtools depth -a "${i}" | awk '{c++; if($3>0) total+=1}END{if (c>0) print (total/c)*100; else print 0}')

echo "${i},${mean_depth},{$breadth}"

done

```
## 4. Generating a consensus sequence
```bash
for i in ./sorted_bam/*.bam; do

samtools consensus -m simple -f fasta -o ${i%.bam}.fasta ${i}

done

```
