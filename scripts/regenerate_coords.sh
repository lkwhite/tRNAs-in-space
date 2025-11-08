#!/bin/bash
# Regenerate all tRNA coordinate files from R2DT outputs
# This script runs R2DT on FASTA files and then processes the JSON outputs

set -e  # Exit on error

echo "========================================="
echo "tRNAs-in-space: Regenerate Coordinates"
echo "========================================="
echo ""

# Check if Docker is available
if ! command -v docker &> /dev/null; then
    echo "Error: Docker is required but not installed."
    echo "Please install Docker from https://www.docker.com/get-started"
    exit 1
fi

# Check if R2DT image is available
if ! docker images | grep -q "rnacentral/r2dt"; then
    echo "Pulling R2DT Docker image..."
    docker pull rnacentral/r2dt
fi

# Process each organism
for fasta in fastas/*.fa; do
    basename=$(basename "$fasta" .fa)
    json_dir="outputs/${basename}_jsons"
    output_tsv="outputs/${basename}_global_coords.tsv"

    echo "Processing $basename..."

    # Run R2DT if JSON directory doesn't exist
    if [ ! -d "$json_dir" ]; then
        echo "  Running R2DT annotation..."
        docker run --rm -v "$(pwd):/data" rnacentral/r2dt \
            r2dt.py gtrnadb draw "/data/$fasta" "/data/$json_dir"
    else
        echo "  Using existing R2DT outputs in $json_dir"
    fi

    # Generate coordinates
    echo "  Generating global coordinates..."
    python scripts/trnas_in_space.py "$json_dir" "$output_tsv"

    echo "  ✓ Generated $output_tsv"
    echo ""
done

echo "========================================="
echo "✓ All datasets regenerated successfully!"
echo "========================================="
