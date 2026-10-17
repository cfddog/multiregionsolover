# =============================================================================
# MultiPhysMixSolver -- 顶层 makefile
#
#   目录布局（源码与中间产物彻底分离）：
#     src/      仅 Fortran 源码 (*.f90)
#     build/    编译中间文件 (*.o / *.mod / *.d)，由本 makefile 自动创建
#     bin/      可执行文件 (MultiPhysMixSolver + 兼容软链 opencfd-ec1.16a.out)
#     cases/    算例；进入算例目录运行 ../../bin/MultiPhysMixSolver
#
#   用法：
#     make                 # 增量编译 -> bin/MultiPhysMixSolver
#     make -j4             # 并行编译（模块依赖由 -MMD 自动收集，可安全并行）
#     make clean           # 删除 build/ 与 bin/
#     make info            # 打印编译器与选项
#
#   兼容：cd src && make 仍然可用（src/makefile 转发到本文件）。
# =============================================================================

f77 = mpif90
opt = -g -std=legacy -ffree-line-length-none -freal-4-real-8 -O3 
#opt=-O3 -openmp

SRCDIR   = src
BUILDDIR = build
BINDIR   = bin
EXE      = MultiPhysMixSolver

srs=  sub_modules.f90 \
      multiphysmixsolver.f90 \
      sub_interfaces.f90 sub_init.f90 sub_Residual.f90 \
      sub_flux_split.f90 sub_boundary.f90 sub_turbulence_BL.f90 \
      sub_turbulence_SA.f90 sub_turbulence_NewSA.f90 sub_turbulence_SST.f90 \
      sub_time_advance.f90  sub_update_buffer_mpi.f90 \
      sub_Post.f90  sub_time_acceleraction.f90 \
      sub_LU_SGS.f90  \
      sub_IO.f90  sub_geometry.f90 \
      sub_partation_mpi.f90 sub_convert_inp.f90 \
      sub_read_parameter.f90 sub_scheme.f90 sub_comput_dw.f90 \
	  sub_debug.f90 sub_filtering.f90 sub_boundary_user.f90  \
	  sub_udate_Meshlink_mpi.f90  \
	  sub_Finite_Difference1.f90  sub_Finite_Difference2.f90 \
	  sub_restart.f90 \
	  sub_Post_timeAverage.f90 sub_limitflow.f90 sub_multi_region.f90 \
	  sub_lowspeed.f90 sub_porous.f90 sub_lowspeed_ac.f90 sub_porous_ac.f90

OBJS = $(addprefix $(BUILDDIR)/,$(srs:.f90=.o))

# -J: 模块文件写入 build/；-I: 从 build/ 查找模块
FLAGS = $(opt) -J$(BUILDDIR) -I$(BUILDDIR)

TARGET = $(BINDIR)/$(EXE)

default all: $(TARGET)

$(TARGET): $(OBJS) | $(BINDIR)
	$(f77) $(opt) -o $@ $(OBJS)
	@ln -sf $(EXE) $(BINDIR)/opencfd-ec1.16a.out   # 兼容旧名（旧算例脚本/软链）
	@echo "==> $(TARGET)  (+ $(BINDIR)/opencfd-ec1.16a.out)"

$(BUILDDIR)/%.o: $(SRCDIR)/%.f90 | $(BUILDDIR)
	$(f77) $(FLAGS) -c $< -o $@

$(BUILDDIR) $(BINDIR):
	@mkdir -p $@

clean:
	@rm -rf $(BUILDDIR) $(BINDIR)
	@echo "cleaned: $(BUILDDIR)/  $(BINDIR)/"

info:
	@echo "f77     = $(f77)"
	@echo "FLAGS   = $(FLAGS)"
	@echo "TARGET  = $(TARGET)"
	@echo "OBJS    = $(words $(OBJS)) object(s) in $(BUILDDIR)/"

# 文件级模块依赖（由 util/gen_moddep.sh 生成）：改过 `use` 关系后请重新生成，
# 这样 `make -j`（包括首次构建）才会按正确顺序编译。
include $(SRCDIR)/modules.dep

.PHONY: default all clean info
