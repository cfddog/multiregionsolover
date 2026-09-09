# fluid_solid / run_conj_11 — code-11 可压流体↔固体 真共轭验收

在 fluid_solid 网格上启用**真正的 code-11 共轭界面**(不手工改 bc3d.inp):
`src/sub_IO.f90::align_interface_to_bc_msg` 会把“物理墙码 2 + bc3d_interface 真实配对”
的耦合面在运行时提升为 `BC_Interface_FluidSolid(11)`;同时把 LS/AC 边界几处严格的
`bc<0` 判断改为 `is_interface_bc`,使 couple/边界/solid 统一把它当接口处理。

运行:`Iflag_Couple_Scheme=1`(code-11 CHT 交错驱动),Kstep_Couple_Comp=200、
warm=1、outer=3、Tol_Couple_Tw=1e-2。

结果(run.log):外轮 1 固体侧 T_w≈800 K → 外轮 2 `max|ΔT_w|=5.0e-3 < 1e-2`
**收敛**;界面 q_w≈29–40 W/m²(近热平衡)。couple 共轭分支真实触发、界面温度往返
正常,无 NaN。
