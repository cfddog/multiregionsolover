# fluid_solid / run_staggered_11 — code-11 可压流体↔固体 CHT 分段交错试跑

将 staggered 分段交错用于 BLOCK_FLUID↔BLOCK_SOLID(code-11)。外层:可压块按调度推进
(`Kstep_Couple_Comp=200`、`Niter_Couple_Warm=2`、减半至 `Kstep_Couple_Min=50`),
固体段用 `solid_solver_one_block` 全收敛,返回固体侧界面 T_w。

试跑(5 轮外层)调度/固体求解/界面输出均正常、无 NaN;T_w 单调向平衡推进
(每轮 |dT_w| 见 run.log),最终 T_w≈558 K,未达 `Tol_Couple_Tw` 前已用完外轮数。
> 注意:本网格(code-11 自检几何)的跨类界面在 `bc_msg2` 中按物理墙式(bc=2+
> 邻接字段)表示,`couple_fluid_solid_interfaces` 的 is_interface 过滤不触发共轭反馈;
> 因此本目录主要验证交错调度与代码路径。真正 code-11 显式标记界面的共轭算例
> 需另建网格/接口文件后验收。
