## 1. BLASTn searches
```bash
#!/bin/bash

# Usage: bash autoblast.sh bellicosus

# To use this script, make a folder named DB and index the blastn databases
# inside the DB folder using the command "makeblastdb -in <ref>.fasta -dbtype nucl -out <dbname>"

# Check if the database argument was provided
if [ -z "$1" ]; then
  echo "Error: Please provide a database."
  exit 1
fi

# Assign the database argument to a variable
db_name="$1"

echo "Searching in blast database: $db_name"

# Loop through each query file
for query_file in queries/*; do
    
  # Get the basename of the query file
  query_name=$(basename "$query_file")

  # Create the output filename
  output="${query_name%.*}_${db_name}.txt"
    
  # Blast    
  echo "Searching for ${query_name%.*}..."
  blastn -query "${query_file}" -db "DB/${db_name}" -out "output/${output}"

done
```

## 2. Sequence alignment
```bash
# Run MACSE
for i in 1_unaligned/*.fasta; do

file=$(basename "$i")

# Translation table (-gc_def) is 1 for nuclear-encoded sequences, 5 for mt-encoded sequences
java -jar macse.jar -prog alignSequences \
  -seq  1_unaligned/${file} \
  -out_NT 2_aligned_NT/${file} \
  -out_AA 3_aligned_AA/${file} \
  -gc_def 1

done

# Run Gblocks for the nucleotide dataset
for i in 2_aligned_NT/*.fasta; do
Gblocks ${i} -b5=h -t=c
done

# Run Gblocks for the protein dataset
for i in 3_aligned_AA/*.fasta; do
Gblocks ${i} -b5=h -t=p
done

```

## 3. TargetP
```bash
# Applied for identifying mt-targeting cleavage sites of the N-mt dataset
for i in 3_aligned_AA/*; do

targetp -fasta ${i} -org non-pl -format short -prefix targetp

done
```
