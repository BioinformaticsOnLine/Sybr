#!/usr/bin/env bash

# Check if resolution arguments are provided
if [ $# -eq 0 ]; then
    echo "Usage: $0 <reference_species> [resolution1 resolution2 resolution3 ...]"
    echo "Example: $0 Adineta_vaga 100 300 500"
    echo "If no resolutions provided, default: 100 300 500"
    exit 1
fi

# First argument is the reference species name
REF_SPECIES="$1"
shift

# Get resolutions from arguments or use defaults
if [ $# -gt 0 ]; then
    RESOLUTIONS=("$@")
else
    RESOLUTIONS=(100 300 500)
fi

# Output file
OUT="all_sequence_lengths.txt"
> "$OUT"

# Process each FASTA file
for fasta in *.fa *.fasta *.fna; do
    [ -e "$fasta" ] || continue
    
    species=$(basename "$fasta")
    species=${species%%.*}
    
    # Extract sequence lengths and store in temporary file
    temp_file=$(mktemp)
    seqkit fx2tab -n -l "$fasta" | awk -v sp="$species" '{
        print $0 "\t" sp
    }' > "$temp_file"
    
    # Check if this is the reference species
    if [ "$species" = "$REF_SPECIES" ]; then
        # For reference species, create entries for each resolution
        for res in "${RESOLUTIONS[@]}"; do
            awk -v sp="$species" -v res="$res" '{
                chr++
                print chr "\t" $2 "\t" sp ":" res "k"
            }' "$temp_file" >> "$OUT"
        done
    else
        # For other species, create single entries
        awk -v sp="$species" '{
            chr++
            print chr "\t" $2 "\t" sp
        }' "$temp_file" >> "$OUT"
    fi
    
    # Clean up temporary file
    rm "$temp_file"
done

echo "Output written to: $OUT"
echo "Reference species '$REF_SPECIES' processed with resolutions: ${RESOLUTIONS[*]}k"
