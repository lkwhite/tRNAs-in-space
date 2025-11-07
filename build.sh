#!/bin/bash
# Build and validation script for tRNAs-in-space project

set -e  # Exit on error

echo "========================================="
echo "tRNAs-in-space Build & Validation"
echo "========================================="
echo ""

# Check Python version
echo "1. Checking Python version..."
python --version
echo ""

# Check dependencies
echo "2. Checking dependencies..."
if ! python -c "import pandas, numpy" 2>/dev/null; then
    echo "   Installing dependencies from requirements.txt..."
    pip install -q -r requirements.txt
    echo "   ✓ Dependencies installed"
else
    echo "   ✓ Dependencies already satisfied"
fi
echo ""

# Validate main script syntax
echo "3. Validating Python syntax..."
python -m py_compile scripts/trnas_in_space.py
echo "   ✓ scripts/trnas_in_space.py syntax OK"
python -m py_compile test_trnas_in_space.py
echo "   ✓ test_trnas_in_space.py syntax OK"
echo ""

# Run tests
echo "4. Running tests..."
python test_trnas_in_space.py
echo ""

# Check output files
echo "5. Validating output files..."
for file in outputs/*.tsv; do
    if [ -f "$file" ]; then
        lines=$(wc -l < "$file")
        echo "   ✓ $(basename "$file"): $lines lines"
    fi
done
echo ""

echo "========================================="
echo "✓ Build completed successfully!"
echo "========================================="
