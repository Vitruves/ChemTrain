# Makefile for ChemTrain project

# Compiler settings
CXX := g++
CC := gcc
CMAKE := cmake
NINJA := ninja
UV := uv

# Default build type
BUILD_TYPE ?= Release

# Build directory
BUILD_DIR := build

# Number of parallel jobs for compilation
JOBS ?= $(shell nproc)

# Executable name
EXECUTABLE := chemtrain

.PHONY: all clean configure build release debug relwithdebinfo deps test format lint install uninstall compiledb help

# Default target
all: release

variance_analyzer: scripts/variance_analyzer.c
	@echo "Building variance analyzer..."
	@cd scripts && gcc -o variance_analyzer variance_analyzer.c -O3 -fopenmp -lm
	
# Clean build directory
clean:
	@echo "Cleaning build directory..."
	@rm -rf $(BUILD_DIR)

# Create build directory (ensured by targets needing it)
ensure_build_dir:
	@mkdir -p $(BUILD_DIR)

# Add a new target for generating the C registry header if it doesn't exist already
ensure_cregistry:
	@if [ ! -f src/cregistry.h ]; then \
		echo "Creating C registry header..."; \
		touch src/cregistry.h; \
	fi

# Configure CMake (internal helper)
# This target is not meant to be called directly, but used by build types
configure: ensure_build_dir ensure_cregistry
	@echo "Configuring CMake ($(BUILD_TYPE))..."
	@cd $(BUILD_DIR) && $(CMAKE) -GNinja \
		-DCMAKE_BUILD_TYPE=$(BUILD_TYPE) \
		-DCMAKE_CXX_COMPILER=$(CXX) \
		-DCMAKE_C_COMPILER=$(CC) \
		..

# Build target (internal helper)
# This target is not meant to be called directly, but used by build types
build: configure
	@echo "Building ChemTrain ($(BUILD_TYPE))..."
	@cd $(BUILD_DIR) && $(NINJA) -j$(JOBS)

# --- User-facing Build Targets ---
release:
	@$(MAKE) build BUILD_TYPE=Release JOBS=$(JOBS)

debug:
	@$(MAKE) build BUILD_TYPE=Debug JOBS=$(JOBS)

relwithdebinfo:
	@$(MAKE) build BUILD_TYPE=RelWithDebInfo JOBS=$(JOBS)

# --- Other Targets ---

# Install dependencies using uv
deps:
	@echo "Installing Python dependencies..."
	@$(UV) pip install -r requirements.txt

# Run tests
# Note: Tests should probably use the currently configured BUILD_TYPE
test: ensure_build_dir
	@echo "Running tests (using BUILD_TYPE=$(BUILD_TYPE))..."
	@cd $(BUILD_DIR) && $(CMAKE) .. -DCMAKE_BUILD_TYPE=$(BUILD_TYPE) # Reconfigure if needed
	@cd $(BUILD_DIR) && $(NINJA) -j$(JOBS) test # Build the 'test' target if defined in CMake, or run ctest
# If ctest is preferred directly:
# test: build
# 	@echo "Running tests (using BUILD_TYPE=$(BUILD_TYPE))..."
# 	@cd $(BUILD_DIR) && ctest --output-on-failure -C $(BUILD_TYPE)

# Format code
format:
	@echo "Formatting code..."
	@find src tests examples \( -iname *.hpp -o -iname *.cpp -o -iname *.c -o -iname *.h \) -exec clang-format -i {} \;

# Static analysis
lint:
	@echo "Running static analysis (requires compile_commands.json)..."
	@if [ ! -f compile_commands.json ]; then $(MAKE) compiledb; fi
	@find src tests examples \( -iname *.hpp -o -iname *.cpp -o -iname *.c -o -iname *.h \) -print0 | xargs -0 clang-tidy -p $(BUILD_DIR)

# Install the built binaries
# Note: Installs the currently configured BUILD_TYPE
install: build
	@echo "Installing ChemTrain ($(BUILD_TYPE))..."
	@cd $(BUILD_DIR) && $(NINJA) install

# Uninstall
uninstall: ensure_build_dir
	@echo "Uninstalling ChemTrain..."
	@if [ -f $(BUILD_DIR)/install_manifest.txt ]; then \
		xargs rm -f < $(BUILD_DIR)/install_manifest.txt; \
	else \
		echo "No install manifest found, trying 'ninja uninstall'... (may fail)"; \
		cd $(BUILD_DIR) && $(NINJA) uninstall; \
	fi

# Generate compile_commands.json for IDE support
# Uses the current BUILD_TYPE setting
compiledb: ensure_build_dir
	@echo "Generating compile_commands.json (using BUILD_TYPE=$(BUILD_TYPE))..."
	@cd $(BUILD_DIR) && $(CMAKE) -GNinja \
		-DCMAKE_BUILD_TYPE=$(BUILD_TYPE) \
		-DCMAKE_EXPORT_COMPILE_COMMANDS=ON \
		..
	@ln -sf $(BUILD_DIR)/compile_commands.json .

# Help target
help:
	@echo "ChemTrain Build System"
	@echo
	@echo "Usage: make [TARGET] [OPTIONS]"
	@echo
	@echo "Available targets:"
	@echo "  all            - Build ChemTrain (Release mode)"
	@echo "  release        - Build with optimizations"
	@echo "  debug          - Build with debug symbols"
	@echo "  relwithdebinfo - Build with both optimizations and debug symbols"
	@echo "  clean          - Remove build directory and compile_commands.json link"
	@echo "  deps           - Install Python dependencies using uv"
	@echo "  test           - Build and run tests for the current BUILD_TYPE"
	@echo "  format         - Format code using clang-format"
	@echo "  lint           - Run static analysis using clang-tidy"
	@echo "  install        - Install built binaries for the current BUILD_TYPE"
	@echo "  uninstall      - Remove installed files (best effort)"
	@echo "  compiledb      - Generate compile_commands.json for IDE support"
	@echo "  help           - Display this help message"
	@echo
	@echo "Build options (can be set on the command line):"
	@echo "  BUILD_TYPE=<?> - Set build type (Release, Debug, RelWithDebInfo). Default: Release"
	@echo "  JOBS=<?>       - Number of parallel jobs (Default: system cores)"
	@echo
	@echo "Example usage:"
	@echo "  make release JOBS=8     # Build release version with 8 parallel jobs"
	@echo "  make debug              # Build debug version"
	@echo "  make clean && make      # Clean and rebuild release version"
	@echo "  make test BUILD_TYPE=Debug # Build debug and run tests"