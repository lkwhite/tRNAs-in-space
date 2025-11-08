.PHONY: help install test lint format clean build check all

# Default target
help:
	@echo "tRNAs-in-space Development Tasks"
	@echo "================================="
	@echo ""
	@echo "Available targets:"
	@echo "  make install    - Install all dependencies"
	@echo "  make test       - Run all tests"
	@echo "  make lint       - Run all linters"
	@echo "  make format     - Format code with black and isort"
	@echo "  make check      - Run linters without making changes"
	@echo "  make build      - Run the build script"
	@echo "  make clean      - Remove generated files"
	@echo "  make all        - Run format, lint, test, and build"
	@echo ""

# Install dependencies
install:
	@echo "Installing dependencies..."
	pip install -r requirements.txt
	@echo "✓ Dependencies installed"

# Run tests
test:
	@echo "Running tests..."
	@python test_trnas_in_space.py
	@echo ""
	@echo "Running tests with pytest..."
	@pytest test_trnas_in_space.py -v
	@echo "✓ All tests passed"

# Run tests with coverage
test-cov:
	@echo "Running tests with coverage..."
	@pytest test_trnas_in_space.py --cov=scripts --cov-report=term --cov-report=html
	@echo "✓ Coverage report generated in htmlcov/"

# Format code
format:
	@echo "Formatting code with black..."
	@black scripts/ test_trnas_in_space.py
	@echo "Sorting imports with isort..."
	@isort scripts/ test_trnas_in_space.py
	@echo "✓ Code formatted"

# Check code style without changes
check:
	@echo "Checking code formatting..."
	@black --check --diff scripts/ test_trnas_in_space.py || true
	@echo ""
	@echo "Checking import order..."
	@isort --check-only --diff scripts/ test_trnas_in_space.py || true
	@echo ""
	@echo "Running flake8..."
	@flake8 scripts/ test_trnas_in_space.py --max-line-length=100 --statistics || true
	@echo ""
	@echo "Running mypy..."
	@mypy scripts/ test_trnas_in_space.py --ignore-missing-imports || true

# Lint code
lint:
	@echo "Running flake8..."
	@flake8 scripts/ test_trnas_in_space.py --count --select=E9,F63,F7,F82 --show-source --statistics
	@flake8 scripts/ test_trnas_in_space.py --count --exit-zero --max-complexity=10 --max-line-length=100 --statistics
	@echo ""
	@echo "Type checking with mypy..."
	@mypy scripts/ test_trnas_in_space.py --ignore-missing-imports || true
	@echo "✓ Linting complete"

# Build and validate
build:
	@echo "Running build script..."
	@chmod +x build.sh
	@./build.sh
	@echo "✓ Build complete"

# Clean generated files
clean:
	@echo "Cleaning generated files..."
	@rm -rf __pycache__ scripts/__pycache__ .pytest_cache .mypy_cache .coverage htmlcov/ *.pyc
	@find . -type d -name "__pycache__" -exec rm -rf {} + 2>/dev/null || true
	@find . -type f -name "*.pyc" -delete
	@echo "✓ Clean complete"

# Run everything
all: format lint test build
	@echo ""
	@echo "========================================"
	@echo "✓ All tasks completed successfully!"
	@echo "========================================"
