#!/bin/bash
set -euo pipefail

# Only run in remote Claude Code environment (web)
if [ "${CLAUDE_CODE_REMOTE:-}" != "true" ]; then
  exit 0
fi

echo "========================================="
echo "Setting up tRNAs-in-space environment..."
echo "========================================="
echo ""

# Install Python dependencies
if [ -f "$CLAUDE_PROJECT_DIR/requirements.txt" ]; then
  echo "📦 Installing Python dependencies from requirements.txt..."
  pip install -q -r "$CLAUDE_PROJECT_DIR/requirements.txt"
  echo "✓ Dependencies installed successfully"
  echo ""
fi

# Set PYTHONPATH to include scripts directory
if [ -n "${CLAUDE_ENV_FILE:-}" ]; then
  echo 'export PYTHONPATH="$CLAUDE_PROJECT_DIR:$CLAUDE_PROJECT_DIR/scripts"' >> "$CLAUDE_ENV_FILE"
fi

# Verify Python version
echo "🐍 Python version:"
python --version
echo ""

# Quick syntax check
echo "🔍 Running quick syntax validation..."
python -m py_compile scripts/trnas_in_space.py 2>/dev/null && echo "  ✓ scripts/trnas_in_space.py syntax OK" || echo "  ⚠ scripts/trnas_in_space.py has syntax errors"
python -m py_compile test_trnas_in_space.py 2>/dev/null && echo "  ✓ test_trnas_in_space.py syntax OK" || echo "  ⚠ test_trnas_in_space.py has syntax errors"
echo ""

# Display available commands
echo "📋 Available commands:"
echo "  make test      - Run all tests"
echo "  make lint      - Check code quality"
echo "  make format    - Auto-format code"
echo "  make build     - Run full build"
echo "  make all       - Format, lint, test, and build"
echo "  ./build.sh     - Run build script directly"
echo ""

echo "========================================="
echo "✓ Environment setup complete!"
echo "========================================="
