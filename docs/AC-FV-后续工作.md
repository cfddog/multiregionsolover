# AC-FV(虚拟压缩法)后续工作备忘

分支:`feat/lowspeed-ac-fv`(从 main `2857573` 切出)
今日完成(2026-09-08,提交 `4996011` / `fe6d707` / `5612cdd`):
- AC 求解器 `src/sub_lowspeed_ac.f90`(LS_Algorithm=3):cell-center FV、
  van-Leer MUSCL + AC-Rusanov、显式边界通量、粘性面梯度 Laplacian、
  伪时间隐式 LU-SGS(定常)+ 双时间源项(已接线,DTS 未做专项考核);
  控制量 AC_beta/AC_CFL/AC_CFLv/AC_Max_Iter/AC_Tol/AC_w/AC_Print。
- 考核:channel_ac(Poiseuille Re_H=50,充分发展 u_max 1.495 vs 1.5,
  剖面 RMS≈0.1%);lidcavity_ac Re=100(40² 均匀,u-RMS≈0.021)与
  Re=1000(80² 均匀/轻度拉伸/双曲正切,均匀 u-RMS≈0.044,
  首层 0.0034 加密后 ≈0.034)。

## 明确列为后续工作(本期不实施)
1. 无粘通量改进:把 Rusanov 对角耗散升级为 **AC 特征分裂 / FVS(AUSM、Steger–Warming 型)**
   以降低 Re 较大时强剪切区的数值耗散;
2. **更高阶重构**:当前为 2 阶 van-Leer MUSCL;可上 WENO5/3 阶 upwind(AC 变量逐分量重构),
   边界处降阶处理;
3. Re=1000 方腔底部低涡区残差成因判定:建议 160² 加密做网格收敛,并与 1、2 的耗散
   改进对照(底部 y≈0.06–0.10 对首层间距不敏感 → 疑似角涡耗散/总点数而非首层)。
4. LU-SGS-DTS 非定常专项考核(槽道起动/方腔阶跃;动量 BDF2 源项与 Un/Un1 已接线);
5. T 被动标量(CHT/bl_cht)求解;压力入口复测(porous_plug/darcy);
6. 多区域接口(12/13/16)与高/低速混合构型的 LS_Algorithm=3 兼容性复跑;
7. AC 的 MPI/多块缓冲语义复核(目前按单进程/单块对验证)。
