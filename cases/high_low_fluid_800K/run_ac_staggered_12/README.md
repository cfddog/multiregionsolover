# high_low_fluid_800K / run_ac_staggered_12 — code-12 高低速流体交错(AC)试跑

可压高速区(Block2/3/4, Ma=3, T_inf=800 K)与 BLOCK_LOWSPEED(Block1, code-12 接口)
的双流体**匹配式块 Gauss–Seidel 分段交错**(`Iflag_Couple_Scheme=1`):
每轮先跑高速区(暖机 `Kstep_Couple_Comp=300`、`Niter_Couple_Warm=1`、之后减半至
`Kstep_Couple_Min=50`),再跑低速区一次至收敛(LS_Algorithm=3, AC);每步用既有
`couple_highlow_fluid_face` 做双向单位换算 ghost。低速 AC 求解器在流场收敛后还会
求解温度(U(5), `ac_lowspeed_energy`, LS_k=4.55>0)。外层判敛用界面 T/p/|u|
(`Tol_Couple_Tw/Tol_Couple_p/Tol_Couple_u`)。

结果(见 run.log / iface_highlow.dat):3 轮外层(高速 300+150+50 步)正常调度、无 NaN;
低速段 AC 3000 伪步收敛到 res_q≈4.3e5;界面 T_w≈440–447 K、p_w≈2.77e4 Pa、
|u|_w≈O(1e-3) m/s;外层 dT≈3.5 K / dp≈1.1e2 Pa(未到 `Tol_Couple_p`,需更多轮次)。

复现:`mpirun -np 1 ./opencfd-ec1.16a.out`(编译好的 binary 在 src/)。
