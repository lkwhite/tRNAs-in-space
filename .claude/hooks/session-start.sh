#!/bin/bash
set -euo pipefail

# Only run in remote Claude Code environment (web)
if [ "${CLAUDE_CODE_REMOTE:-}" != "true" ]; then
  exit 0
fi

echo "Setting up tRNAs-in-space environment..."

# Install Python dependencies
if [ -f "$CLAUDE_PROJECT_DIR/requirements.txt" ]; then
  echo "Installing Python dependencies from requirements.txt..."
  pip install -q -r "$CLAUDE_PROJECT_DIR/requirements.txt"
  echo "✓ Dependencies installed successfully"
fi

# Set PYTHONPATH to include scripts directory
if [ -n "${CLAUDE_ENV_FILE:-}" ]; then
  echo 'export PYTHONPATH="$CLAUDE_PROJECT_DIR:$CLAUDE_PROJECT_DIR/scripts"' >> "$CLAUDE_ENV_FILE"
fi

echo "✓ Environment setup complete"
