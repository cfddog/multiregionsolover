# -----------------------------------------------------------------------------
# 兼容入口：真正的构建规则在仓库根目录 makefile
#   src/   只放源码（*.f90）
#   build/ 编译中间文件（*.o / *.mod / *.d）
#   bin/   可执行文件（MultiPhysMixSolver + 软链 opencfd-ec1.16a.out）
# 旧的 `cd src && make` / `cd src && make clean` 仍然可用（转发到 ../makefile）。
# -----------------------------------------------------------------------------
.DEFAULT_GOAL := default

default all:
	@$(MAKE) --no-print-directory -C .. $@

clean:
	@$(MAKE) --no-print-directory -C .. clean

info:
	@$(MAKE) --no-print-directory -C .. info

.PHONY: default all clean info

