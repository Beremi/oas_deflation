SHELL := /bin/bash
.SHELLFLAGS := -eu -o pipefail -c

BUILD_DIR ?= ../oas_deflation-build/release
FAST_BUILD_DIR ?= ../oas_deflation-build/fast
BUILD_TYPE ?= Release
JOBS ?= $(shell nproc)
THREADS ?= 4
MKLROOT ?= /opt/intel/oneapi/mkl/latest
VTK_DIR ?= /usr/lib/cmake/vtk
USE_VTK ?= ON
OAS_CXX_FLAGS ?= -DEIGEN_MKL_NO_DIRECT_CALL
SOLVER ?= EigenLDLT

OAS_BIN := $(BUILD_DIR)/bin/OAS
DOGBONE_ARCHIVE := data/archives/Dogbone.zip
DOGBONE_CASE := data/cases/Dogbone
DOGBONE_MASTER := $(DOGBONE_CASE)/master.inp

.PHONY: configure build configure-fast clean-build dogbone dogbone-solver extract-dogbone sync-upstream status agent-context

configure:
	MKLROOT="$(MKLROOT)" cmake -S . -B "$(BUILD_DIR)" -G "Unix Makefiles" \
		-DCMAKE_BUILD_TYPE="$(BUILD_TYPE)" \
		-DCMAKE_C_COMPILER_LAUNCHER=ccache \
		-DCMAKE_CXX_COMPILER_LAUNCHER=ccache \
		-DCMAKE_CXX_FLAGS="$(OAS_CXX_FLAGS)" \
		-DUSE_PARDISO=ON \
		-DUSE_CHOLMOD=ON \
		-DUSE_VTK="$(USE_VTK)" \
		-DUSE_DOCS=OFF \
		-DVTK_DIR:PATH="$(VTK_DIR)"

build: configure
	cmake --build "$(BUILD_DIR)" --parallel "$(JOBS)"

configure-fast:
	MKLROOT="$(MKLROOT)" cmake -S . -B "$(FAST_BUILD_DIR)" -G "Unix Makefiles" \
		-DCMAKE_BUILD_TYPE="$(BUILD_TYPE)" \
		-DCMAKE_C_COMPILER_LAUNCHER=ccache \
		-DCMAKE_CXX_COMPILER_LAUNCHER=ccache \
		-DCMAKE_CXX_FLAGS="$(OAS_CXX_FLAGS)" \
		-DUSE_PARDISO=OFF \
		-DUSE_CHOLMOD=OFF \
		-DUSE_VTK=OFF \
		-DUSE_DOCS=OFF

clean-build:
	cmake -E rm -rf "$(BUILD_DIR)" "$(FAST_BUILD_DIR)"

extract-dogbone: $(DOGBONE_MASTER)

$(DOGBONE_MASTER): $(DOGBONE_ARCHIVE)
	mkdir -p data/cases
	unzip -q -o "$(DOGBONE_ARCHIVE)" -d data/cases

dogbone: build extract-dogbone
	MKLROOT="$(MKLROOT)" MKL_NUM_THREADS="$(THREADS)" OMP_NUM_THREADS="$(THREADS)" \
		"$(OAS_BIN)" -j "$(THREADS)" "$(DOGBONE_MASTER)"

dogbone-solver: extract-dogbone
	perl -i -pe 's/^solver_type\s+\S+/solver_type\t$(SOLVER)/' "$(DOGBONE_CASE)/solver.inp"
	$(MAKE) dogbone

sync-upstream:
	git fetch upstream
	git rebase upstream/master

status:
	git status --short --branch
	git remote -v
	@test -x "$(OAS_BIN)" && ls -lh "$(OAS_BIN)" || true

agent-context:
	scripts/update-agent-context.sh
