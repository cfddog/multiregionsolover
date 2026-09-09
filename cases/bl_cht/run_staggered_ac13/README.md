# bl_cht / run_staggered_ac13 — code-13 低速(AC)流体↔固体 CHT 分段交错试跑

BLOCK_LOWSPEED(LS_Algorithm=3, AC+温度方程)↔ BLOCK_SOLID(code-13)。LS-AC 段每轮
一次“至收敛”(AC_Max_Iter=10000;流场收敛后 `ac_lowspeed_energy` 解 U(5) 温度),
固体段全 GS 收敛,返回固体侧界面 T_w。

试跑:调度与 AC/固体求解路径正常(第 1 轮固体 GS 收敛、T_w=800 K 来自该自检网格的
固体物理等温壁 BC;第 2 轮 LS-AC 开始推进)。
> 注意:同 fluid_solid/run_staggered_11,该自检网格的 code-13 界面按物理墙式表示,
> 未走 couple_fluid_solid_interfaces 的共轭热平衡分支;本目录主要验证交错调度与
> AC 温度耦合的数据路径。真共轭验收需 code-13 显式标记界面的网格。
