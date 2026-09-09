# bl_cht / run_conj_ac13 — code-13 低速(AC)流体↔固体 真共轭验收

在 bl_cht 网格上启用真正的 code-13 共轭界面(同 code-11 的运行时提升机制:
`align_interface_to_bc_msg` 把“墙码 2+配对”面提升为 `BC_Interface_LowSolid(13)`)。
流体侧为 BLOCK_LOWSPEED + **LS_Algorithm=3(AC)**(流场+`ac_lowspeed_energy` 温度)。

运行:`Iflag_Couple_Scheme=1`(code-13 CHT 交错驱动),LS AC 每外轮至收敛(AC 20 步
即收敛),固体段全 GS;outer=3、Tol_Couple_Tw=1e-2。

结果(run.log):T_w 逐轮 337→387→424 K 单调上升、q_w≈1e5 W/m²,共轭热交换真实
进行;3 轮用完时 max|ΔT_w|≈43 K(未达 Tol,需更多外轮/更长 AC 收敛)。无 NaN。
