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
USE_AMGCL ?= ON
USE_HYPRE ?= OFF
OAS_CXX_FLAGS ?= -DEIGEN_MKL_NO_DIRECT_CALL
SOLVER ?= EigenLDLT
AMGCL_TOL ?= 1e-6
AMGCL_MAXIT ?= 500
AMGCL_EPS_STRONG ?= 0
AMGCL_RELAX ?= 1
AMGCL_BLOCK_SIZE ?= 1
AMGCL_COARSE_ENOUGH ?= 0
AMGCL_NPRE ?= 1
AMGCL_NPOST ?= 1
AMGCL_NCYCLE ?= 1
AMGCL_NULLSPACE ?= 1
AMGCL_REUSE_INITIAL_GUESS ?= 0
AMGCL_BLOCK_RELAX ?= ilu0
DFGMRES_TOL ?= 1e-6
DFGMRES_TRUE_TOL ?= $(DFGMRES_TOL)
DFGMRES_MAXIT ?= 500
DFGMRES_RESTART ?= 80
DFGMRES_N ?= 20
DFGMRES_EPS ?= 1e-3
DFGMRES_PRECONDITIONER ?= amgcl
HYPRE_TOL ?= 1e-6
HYPRE_MAXIT ?= 500
HYPRE_COARSEN_TYPE ?= 8
HYPRE_INTERP_TYPE ?= 6
HYPRE_STRONG_THRESHOLD ?= 0.5
HYPRE_NODAL ?= 4
HYPRE_RELAX_TYPE ?= 6
HYPRE_P_MAX ?= 4
PROFILE_DIR ?= $(DOGBONE_CASE)/results
PROFILE_SOLVER ?= $(shell printf '%s' '$(SOLVER)' | tr '[:upper:]' '[:lower:]')
PROFILE_STAMP ?= $(shell date +%Y%m%d-%H%M%S)
OUT_DIR ?= results/dogbone-$(PROFILE_SOLVER)-$(PROFILE_STAMP)

OAS_BIN := $(BUILD_DIR)/bin/OAS
DOGBONE_ARCHIVE := data/archives/Dogbone.zip
DOGBONE_CASE := data/cases/Dogbone
DOGBONE_MASTER := $(DOGBONE_CASE)/master.inp
TSN65_CASE := data/cases/TS-N_65
TSN65_MASTER := $(TSN65_CASE)/master.inp

.PHONY: configure build configure-fast clean-build dogbone dogbone-solver dogbone-profile dogbone-amgcl-prefix dogbone-hypre-prefix dogbone-dfgmres-prefix dogbone-dfgmres-profile dogbone-dfgmres-hypre-sweep dogbone-amgcl-profile dogbone-amgcl-sweep dogbone-amgcl-tolerance-sweep tsn65-amgcl-prefix tsn65-hypre-prefix linear-profile-report extract-dogbone sync-upstream status agent-context

configure:
	MKLROOT="$(MKLROOT)" cmake -S . -B "$(BUILD_DIR)" -G "Unix Makefiles" \
		-DCMAKE_BUILD_TYPE="$(BUILD_TYPE)" \
		-DCMAKE_C_COMPILER_LAUNCHER=ccache \
		-DCMAKE_CXX_COMPILER_LAUNCHER=ccache \
		-DCMAKE_CXX_FLAGS="$(OAS_CXX_FLAGS)" \
		-DUSE_PARDISO=ON \
		-DUSE_CHOLMOD=ON \
		-DUSE_AMGCL="$(USE_AMGCL)" \
		-DUSE_HYPRE="$(USE_HYPRE)" \
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
		-DUSE_AMGCL="$(USE_AMGCL)" \
		-DUSE_HYPRE="$(USE_HYPRE)" \
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

dogbone-profile: build extract-dogbone
	scripts/run-dogbone-profile.sh "$(DOGBONE_CASE)" "$(DOGBONE_MASTER)" "$(OAS_BIN)" "$(MKLROOT)" "$(THREADS)" "$(SOLVER)"
	$(MAKE) linear-profile-report PROFILE_DIR="$(PROFILE_DIR)" OUT_DIR="$(OUT_DIR)" SOLVER="$(SOLVER)"

dogbone-amgcl-prefix: build extract-dogbone
	scripts/run-oas-profile.sh "$(DOGBONE_CASE)" "$(DOGBONE_MASTER)" "$(OAS_BIN)" "$(THREADS)" "AmgclCGElastic" \
		"Dogbone AMGCL Elastic Prefix Linear-Solve Profile" Dogbone \
		"results/dogbone-amgcl-elastic-prefix-$(PROFILE_STAMP)" \
		total_time=0.01 max_iterations=5 min_iterations=5 limit_tolerance=10 \
		linear_solver_replay_dump=1 linear_solver_replay_limit=5 \
		amgcl_tolerance="$(AMGCL_TOL)" amgcl_max_iterations="$(AMGCL_MAXIT)" \
		amgcl_eps_strong="$(AMGCL_EPS_STRONG)" amgcl_relax="$(AMGCL_RELAX)" \
		amgcl_block_size="$(AMGCL_BLOCK_SIZE)" amgcl_coarse_enough="$(AMGCL_COARSE_ENOUGH)" \
		amgcl_npre="$(AMGCL_NPRE)" amgcl_npost="$(AMGCL_NPOST)" amgcl_ncycle="$(AMGCL_NCYCLE)" \
		amgcl_near_nullspace="$(AMGCL_NULLSPACE)" amgcl_elastic_full_lift=1 amgcl_use_block_backend=1 \
		amgcl_block_relaxation="$(AMGCL_BLOCK_RELAX)" \
		amgcl_reuse_initial_guess="$(AMGCL_REUSE_INITIAL_GUESS)"

dogbone-hypre-prefix: build extract-dogbone
	scripts/run-oas-profile.sh "$(DOGBONE_CASE)" "$(DOGBONE_MASTER)" "$(OAS_BIN)" "$(THREADS)" "HypreBoomerAMGCG" \
		"Dogbone hypre BoomerAMG-CG Prefix Linear-Solve Profile" Dogbone \
		"results/dogbone-hypre-boomeramg-prefix-$(PROFILE_STAMP)" \
		total_time=0.01 max_iterations=5 min_iterations=5 limit_tolerance=10 \
		linear_solver_replay_dump=1 linear_solver_replay_limit=5 \
		hypre_tolerance="$(HYPRE_TOL)" hypre_max_iterations="$(HYPRE_MAXIT)" \
		hypre_coarsen_type="$(HYPRE_COARSEN_TYPE)" hypre_interp_type="$(HYPRE_INTERP_TYPE)" \
		hypre_strong_threshold="$(HYPRE_STRONG_THRESHOLD)" hypre_nodal="$(HYPRE_NODAL)" \
		hypre_relax_type="$(HYPRE_RELAX_TYPE)" hypre_p_max="$(HYPRE_P_MAX)" \
		hypre_use_dof_functions=0 hypre_use_interp_vectors=0 hypre_interp_vec_variant=2

dogbone-dfgmres-prefix: build extract-dogbone
	scripts/run-oas-profile.sh "$(DOGBONE_CASE)" "$(DOGBONE_MASTER)" "$(OAS_BIN)" "$(THREADS)" "DeflatedFGMRES" \
		"Dogbone DeflatedFGMRES Prefix Linear-Solve Profile" Dogbone \
		"results/dogbone-dfgmres-prefix-$(PROFILE_STAMP)" \
		total_time=0.01 max_iterations=5 min_iterations=5 limit_tolerance=10 \
		linear_solver_replay_dump=1 linear_solver_replay_limit=5 \
		dfgmres_tolerance="$(DFGMRES_TOL)" dfgmres_true_tolerance="$(DFGMRES_TRUE_TOL)" \
		dfgmres_max_iterations="$(DFGMRES_MAXIT)" dfgmres_restart="$(DFGMRES_RESTART)" \
		dfgmres_deflation_vectors="$(DFGMRES_N)" dfgmres_deflation_eps="$(DFGMRES_EPS)" \
		dfgmres_collect_newton_steps=1 dfgmres_preconditioner="$(DFGMRES_PRECONDITIONER)" \
		amgcl_tolerance="$(AMGCL_TOL)" amgcl_max_iterations="$(AMGCL_MAXIT)" \
		amgcl_eps_strong="$(AMGCL_EPS_STRONG)" amgcl_relax="$(AMGCL_RELAX)" \
		amgcl_block_size="$(AMGCL_BLOCK_SIZE)" amgcl_coarse_enough="$(AMGCL_COARSE_ENOUGH)" \
		amgcl_npre="$(AMGCL_NPRE)" amgcl_npost="$(AMGCL_NPOST)" amgcl_ncycle="$(AMGCL_NCYCLE)" \
		amgcl_near_nullspace="$(AMGCL_NULLSPACE)" amgcl_elastic_full_lift=1 amgcl_use_block_backend=1 \
		amgcl_backend=hybrid amgcl_block_relaxation="$(AMGCL_BLOCK_RELAX)"

dogbone-dfgmres-profile: build extract-dogbone
	scripts/run-oas-profile.sh "$(DOGBONE_CASE)" "$(DOGBONE_MASTER)" "$(OAS_BIN)" "$(THREADS)" "DeflatedFGMRES" \
		"Dogbone DeflatedFGMRES Linear-Solve Profile" Dogbone \
		"results/dogbone-dfgmres-N$(DFGMRES_N)-tol$(DFGMRES_TOL)-$(PROFILE_STAMP)" \
		dfgmres_tolerance="$(DFGMRES_TOL)" dfgmres_true_tolerance="$(DFGMRES_TRUE_TOL)" \
		dfgmres_max_iterations="$(DFGMRES_MAXIT)" dfgmres_restart="$(DFGMRES_RESTART)" \
		dfgmres_deflation_vectors="$(DFGMRES_N)" dfgmres_deflation_eps="$(DFGMRES_EPS)" \
		dfgmres_collect_newton_steps=1 dfgmres_preconditioner="$(DFGMRES_PRECONDITIONER)" \
		amgcl_tolerance="$(AMGCL_TOL)" amgcl_max_iterations="$(AMGCL_MAXIT)" \
		amgcl_eps_strong="$(AMGCL_EPS_STRONG)" amgcl_relax="$(AMGCL_RELAX)" \
		amgcl_block_size="$(AMGCL_BLOCK_SIZE)" amgcl_coarse_enough="$(AMGCL_COARSE_ENOUGH)" \
		amgcl_npre="$(AMGCL_NPRE)" amgcl_npost="$(AMGCL_NPOST)" amgcl_ncycle="$(AMGCL_NCYCLE)" \
		amgcl_near_nullspace="$(AMGCL_NULLSPACE)" amgcl_elastic_full_lift=1 amgcl_use_block_backend=1 \
		amgcl_backend=hybrid amgcl_block_relaxation="$(AMGCL_BLOCK_RELAX)"

dogbone-dfgmres-hypre-sweep:
	$(MAKE) build USE_VTK=OFF USE_AMGCL=ON USE_HYPRE=ON THREADS="$(THREADS)"
	$(MAKE) extract-dogbone
	scripts/run-dogbone-dfgmres-hypre-sweep.py --threads "$(THREADS)" --oas-bin "$(OAS_BIN)" --mklroot "$(MKLROOT)"

dogbone-amgcl-profile: build extract-dogbone
	scripts/run-oas-profile.sh "$(DOGBONE_CASE)" "$(DOGBONE_MASTER)" "$(OAS_BIN)" "$(THREADS)" "AmgclCG" \
		"Dogbone AmgclCG Linear-Solve Profile" Dogbone \
		"results/dogbone-amgclcg-$(PROFILE_STAMP)" \
		amgcl_tolerance="$(AMGCL_TOL)" amgcl_max_iterations="$(AMGCL_MAXIT)" \
		amgcl_eps_strong="$(AMGCL_EPS_STRONG)" amgcl_relax="$(AMGCL_RELAX)" \
		amgcl_block_size="$(AMGCL_BLOCK_SIZE)" amgcl_coarse_enough="$(AMGCL_COARSE_ENOUGH)" \
		amgcl_npre="$(AMGCL_NPRE)" amgcl_npost="$(AMGCL_NPOST)" amgcl_ncycle="$(AMGCL_NCYCLE)" \
		amgcl_near_nullspace="$(AMGCL_NULLSPACE)" amgcl_block_relaxation="$(AMGCL_BLOCK_RELAX)" \
		amgcl_reuse_initial_guess="$(AMGCL_REUSE_INITIAL_GUESS)"

dogbone-amgcl-sweep: build extract-dogbone
	scripts/sweep-amgcl.py --mode dogbone-prefix --threads "$(THREADS)" --oas-bin "$(OAS_BIN)" \
		--tolerances "1e-6" --max-iterations "500,1000" --eps-strong "0,0.02,0.08" \
		--relax "0.8,1" --npre "1,2" --npost "1,2" --near-nullspace "1" \
		--reuse-initial-guess "0" \
		--limit "$${LIMIT:-0}"

dogbone-amgcl-tolerance-sweep: build extract-dogbone
	scripts/sweep-amgcl.py --mode dogbone-full --threads "$(THREADS)" --oas-bin "$(OAS_BIN)" \
		--tolerances "$${TOLS:-1e-4,1e-5,1e-6}" --max-iterations "$${AMGCL_MAXITS:-1000}" \
		--eps-strong "$${AMGCL_EPS_STRONGS:-0}" --relax "$${AMGCL_RELAXES:-1}" \
		--npre "$${AMGCL_NPRES:-1}" --npost "$${AMGCL_NPOSTS:-1}" --near-nullspace "$${AMGCL_NULLSPACES:-1}" \
		--reuse-initial-guess "$${AMGCL_REUSE_INITIAL_GUESSES:-0}"

tsn65-amgcl-prefix: build
	OAS_TIMEOUT="$${OAS_TIMEOUT:-45m}" scripts/run-oas-profile.sh "$(TSN65_CASE)" "$(TSN65_MASTER)" "$(OAS_BIN)" "$(THREADS)" "AmgclCGElastic" \
		"TS-N_65 AMGCL Elastic $(THREADS)-thread First-Step Prefix Linear-Solve Profile" TS-N_65 \
		"results/tsn65-amgcl-elastic-$(THREADS)t-firststep-5it-$(PROFILE_STAMP)" \
		total_time=5.000000e-03 max_iterations=5 min_iterations=5 limit_tolerance=10 \
		linear_solver_profile_matrix_delta=0 linear_solver_replay_dump=0 \
		amgcl_tolerance="$(AMGCL_TOL)" amgcl_max_iterations="$(AMGCL_MAXIT)" \
		amgcl_eps_strong="$(AMGCL_EPS_STRONG)" amgcl_relax="$(AMGCL_RELAX)" \
		amgcl_block_size="$(AMGCL_BLOCK_SIZE)" amgcl_coarse_enough="$(AMGCL_COARSE_ENOUGH)" \
		amgcl_npre="$(AMGCL_NPRE)" amgcl_npost="$(AMGCL_NPOST)" amgcl_ncycle="$(AMGCL_NCYCLE)" \
		amgcl_near_nullspace="$(AMGCL_NULLSPACE)" amgcl_elastic_full_lift=1 amgcl_use_block_backend=1 \
		amgcl_block_relaxation="$(AMGCL_BLOCK_RELAX)" \
		amgcl_reuse_initial_guess="$(AMGCL_REUSE_INITIAL_GUESS)"

tsn65-hypre-prefix: build
	OAS_TIMEOUT="$${OAS_TIMEOUT:-45m}" scripts/run-oas-profile.sh "$(TSN65_CASE)" "$(TSN65_MASTER)" "$(OAS_BIN)" "$(THREADS)" "HypreBoomerAMGCG" \
		"TS-N_65 hypre BoomerAMG-CG $(THREADS)-thread First-Step Prefix Linear-Solve Profile" TS-N_65 \
		"results/tsn65-hypre-boomeramg-$(THREADS)t-firststep-5it-$(PROFILE_STAMP)" \
		total_time=5.000000e-03 max_iterations=5 min_iterations=5 limit_tolerance=10 \
		linear_solver_profile_matrix_delta=0 linear_solver_replay_dump=0 \
		hypre_tolerance="$(HYPRE_TOL)" hypre_max_iterations="$(HYPRE_MAXIT)" \
		hypre_coarsen_type="$(HYPRE_COARSEN_TYPE)" hypre_interp_type="$(HYPRE_INTERP_TYPE)" \
		hypre_strong_threshold="$(HYPRE_STRONG_THRESHOLD)" hypre_nodal="$(HYPRE_NODAL)" \
		hypre_relax_type="$(HYPRE_RELAX_TYPE)" hypre_p_max="$(HYPRE_P_MAX)" \
		hypre_use_dof_functions=0 hypre_use_interp_vectors=0 hypre_interp_vec_variant=2

linear-profile-report:
	python3 scripts/analyze-linear-profile.py \
		--profile-dir "$(PROFILE_DIR)" \
		--out-dir "$(OUT_DIR)" \
		--title "Dogbone $(SOLVER) Linear-Solve Profile" \
		--case Dogbone \
		--threads "$(THREADS)"

sync-upstream:
	git fetch upstream
	git rebase upstream/master

status:
	git status --short --branch
	git remote -v
	@test -x "$(OAS_BIN)" && ls -lh "$(OAS_BIN)" || true

agent-context:
	scripts/update-agent-context.sh
