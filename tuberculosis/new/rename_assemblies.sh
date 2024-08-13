#!/bin/bash

# Loop through each directory in the current directory
for dir in */; do
  # Remove the trailing slash from the directory name
  sampleid=$(basename "$dir" .fa)

  # Define the path to the contigs.fa file
  contigs_file="${dir}contigs.fa"

  # Define the new file name with the SAMPLEID
  new_file_name="${dir}${sampleid}.fa"

  # Check if contigs.fa exists
  if [ -f "$contigs_file" ]; then
    # Rename contigs.fa to SAMPLEID.fa
    mv "$contigs_file" "$new_file_name"
    echo "Renamed $contigs_file to $new_file_name"
  else
    echo "contigs.fa not found in $dir"
  fi
done


#!/bin/bash

# Loop through each directory ending with .fa in the current directory
for dir in *.fa/; do
  # Extract the base name of the directory without the .fa suffix
  new_dir=$(basename "$dir" .fa)

  # Rename the directory by removing the .fa suffix
  mv "$dir" "${new_dir}/"
  echo "Renamed directory $dir to ${new_dir}/"
done


