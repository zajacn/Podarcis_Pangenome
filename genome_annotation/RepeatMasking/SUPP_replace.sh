#!/bin/bash

cat family-*.fa.refiner_cons > output.fasta
grep 'RepeatScout:' families.stk | sed 's"#=GF DE    RepeatModeler Generated - rnd-1_""g' | sed 's/, RepeatScout/\t(RepeatScout Family = /g' | sed  's/\[//g' | sed  's/\]/)/g' | sed 's/://g' | sed 's/RS//g' > replacement.table.txt

# Define the table and fasta file
table="replacement.table.txt"        # Your table with family info
fasta="output.fasta"      # Your FASTA file with the headers to be replaced

# Process the table into a format for easy lookup
awk -F"\t" '{print $1"\t"$2}' $table > table_lookup.txt

# Loop through the FASTA file and modify the headers
awk -F ' ' '
BEGIN {
  # Load the lookup table into an array
  while ((getline < "table_lookup.txt") > 0) {
    split($0, arr, "\t");
    table[arr[1]] = arr[2];
  }
}
{
  # Check if the line is a header (starting with ">")
  if ($0 ~ /^>/) {
    # Extract family name
    split($0, parts, " ");
    family = substr(parts[1], 2);  # Remove ">"
    
    # Replace header with the corresponding value from the table
    if (family in table) {
      print ">" family " " table[family];  # Print new header
    } else {
      print $0;  # Print original header if no match found
    }
  } else {
    print $0;  # For sequence lines, print as is
  }
}
' $fasta > consensi-refined.fa


rm output.fasta table_lookup.txt replacement.table.txt 
