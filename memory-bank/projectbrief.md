# Project Brief — OpenCFD-EC 1.16a 多块/多区域求解框架

## 1. 项目目标
在可压缩密度基基线 **OpenCFD-EC 1.16a**（Li Xinliang，LHD/中科院力学所）之上，扩展
**多块、多区域混合求解**能力，用于流固/流多孔/高低速耦合等工程问题：

| 块类型码 | 名称 | 求解器 | 状态 |
|---|---|---|---|
| 0 | `BLOCK_FLUID` | 可压缩密度基 N-S（RANS 框架、LU-SGS） | 继承基线，完整 |
| 1 | `BLOCK_SOLID` | 固体导热 FVM（Tw/Qw 热边界） | 新增 |
| 2 | `BLOCK_LOWSPEED` | 不可压低速（SIMPLE / SIMPLEC / AC-FV + Rhie–Chow） | 新增 |
| 3 | `BLOCK_POROUS` | 多孔介质（Darcy–Brinkman–Forchheimer + LTNE 双温度） | 新增 |

显示式**接口码**（写在 `bc3d.inp`，标识“块的类型对”，不含方向）：
11 流体↔固体、12 流体↔低速、13 低速↔固体、14 固固、15 流流、
16 低速↔多孔、17 固↔多孔（未实现）、18 多多（未实现）、19 可压↔多孔。

## 2. 代码布局
```
src/                     单一可执行的主程序 + 各子模块（含 makefile）
  opencfd_ec3d_v1.16a.f90  主程序（含交错耦合分派 Iflag_Couple_Scheme==1）
  sub_modules.f90          精度/常量/块类型/接口码/全局控制量/全局数组
  sub_read_parameter.f90   namelist / material.in / porous.inp / solid_bc.inp 读取 + MPI 广播
  sub_IO.f90              网格读写、bc3d 转换、VTK/单文件 SI 输出、界面码运行时提升
  sub_convert_inp.f90      bc3d.inp / bc3d_interface.inp → .inc；连接描述符 L1..L3
  sub_geometry.f90         几何量、网格质量、多块连接检查
  sub_multi_region.f90     多区域耦合例程 + 交错分段驱动（本仓库扩展核心）
  sub_lowspeed.f90         低速 SIMPLE 求解器
  sub_lowspeed_ac.f90      低速 AC（人工压缩）求解器 + 温度方程
  sub_porous.f90           多孔 LTNE 求解器
  sub_porous_ac.f90        多孔 AC 求解器
cases/                   算例（每个含 Mesh3d.x/bc3d.inp/bc3d_interface.inp/material.in/control.ec…）
docs/                    工作日志、程序能力总结、进度摘要、AC-FV 后续工作
memory-bank/             本记忆库（projectbrief/activeContext/progress/worklog/constraints）
util/  run_lid/  sample_code/
```

## 3. 构建与运行
```bash
cd src && make            # mpif90，opt=-O3 -freal-4-real-8 -std=legacy …
                          # 产物：src/opencfd-ec1.16a.out
cd ../cases/<case> && mpirun -np 1 ./opencfd-ec1.16a.out
```
（算例目录里通常把 `opencfd-ec1.16a.out` 软链到 `src/` 的产物。）

## 4. 算例关键输入文件
| 文件 | 作用 |
|---|---|
| `control.ec` | **7 组 namelist**：`$freestream_ec`（来流/参考态/气体物性）、`$flow_ec`（时间步/格式/模型/网格并行/输出）、`$lowspeed_ec`（`LS_*`）、`$ac_ec`（`AC_*`，低速与多孔共用）、`$solid_ec`（固体 GS 控制）、`$porous_ec`（`Porous_*`）、`$couple_ec`（交错耦合调度）。**每组可选、变量可选**（全部有代码默认值）；旧式单组 `control_ec`（128 变量写一起）仍完全兼容，新旧并存时新组优先；组存在但拼错变量名 ⇒ 报错 `stop 1`。详见 `docs/control.ec-说明.md`，模板 `control.ec.template` |
| `Mesh3d.x` / `Mesh3d.dat` | 网格（`Mesh_File_Format=2` / `0`）|
| `bc3d.inp` | 物理边界 + 接口码（每面：i1 i2 j1 j2 k1 k2 bc）|
| `bc3d_interface.inp` | 块间对接关系（界面行可用 `-1` + 续行给 nb1；配对也可由几何自动识别）|
| `material.in` | 块类型表 + 各块 (rho, Cp, k)（固体/多孔骨架用）|
| `porous.inp` | 多孔块：`Block_no eps dp(m) hv(W/m³/K)` |
| `solid_bc.inp` | 固体热边界：`Block_no / Num_faces / face Tw Qw`（Tw>0 等温，Tw<0 给定热流）|

## 5. 面编号约定（**全仓库通用，务必牢记**）
`1=i-, 2=j-, 3=k-, 4=i+, 5=j+, 6=k+`。
对接描述符 `L1,L2,L3`：本块第 k 维与邻块第 `|Lk|` 维相连，符号表正/反向连接
（`sub_convert_inp.f90::Convert_bc`）。

## 6. 参考文档
- `docs/工作日志.md`（830+ 行，逐次开发的权威记录，含历史 commit 说明）
- `docs/程序能力与算例考核总结.md`（能力→算例→验证状态，可作为功能清单）
- `docs/AC-FV-后续工作.md`（低速 AC 专题）
- `docs/进度_fluid_solid三算例.md`（2026-09-18 中断点记录）

## 7. 参考算例：`cases/fluid_solid` 拓扑（2026-09-18 扩展，重要）
- Block1 = 67×101×2 薄板，Block2 = 101×201×2 通道；界面 = **Block1 j=1 ↔ Block2 i=101**
  （共形 67 节点）；网格单位 mm（`Lscale=0.001`）。
- 两块界面法向在**物理上都是 y**（Block1 的 j 与 Block2 的 i 都沿 y）
  ⇒ `U(2)=切向 x、U(3)=法向 y`，物理公式与旧 `(j-,j+)` 拓扑一致，只需换下标。
- 运行时连接描述符：Block1 侧 `L=(-2,1,3)`；Block2（可压）侧 `L=(2,-1,3)`
  ⇒ 本拓扑 = **gas i+ (face 4) ↔ 对侧 j- (face 2)，切向反向**。
- 该拓扑已支持：码 11（`run_staggered_fluid_solid`）、12（`run_staggered_highlow`）、
  19（`run_staggered_fluid_porous` 的 `pair42` 分支，含发汗壁 `wall_bound_blowing`）。
- 三个算例目录：`run_case1_solid`(11) / `run_case2_lowspeed`(12) / `run_case3_porous`(19)。
- 界面码写显式 11/12/19 于 `bc3d.inp`；`bc3d_interface.inp` 沿用 Pointwise 原结构（`-1`+续行给 nb1）。
