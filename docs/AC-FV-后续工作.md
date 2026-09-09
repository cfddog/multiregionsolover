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

## AC-porous 续接(2026-09-09 新增,与 `docs/工作日志.md` 2026-09-09 条目联动)
8. **AC 版 LTNE 双温度标量对**:实现 `src/sub_porous_ac.f90::ac_porous_ltne`(AC 面通量的
   Tf 对流-扩散 + hv 源 + `(1−ε)k_s` 骨架导热 + 骨架面热边界 ghost),验证
   `porous_ltne_1d` 两温度闭式解(目标 ≤0.5%)。✅ **已完成(2026-09-10)**:
   `cases/porous_ltne_1d/run_ac[_g101|_g201]`,`nx=401` 时 Tf RMS=1.22 K / rel≈0.19%、
   Ts RMS=1.17 K、能量缺口 0.28% —— 与 SIMPLE 参考(0.18%)一致;图
   `cases/porous_ltne_1d/validation_ltne_ac{,_101,_201}.png`。
9. **plug/darcy <1% 定量收口**:统一压差拟合窗口与出口压力口径,溯源 AC 与 SIMPLE
   (u_bar 8.06e-3 vs 7.85e-3 ≈ +2.7%)差异;随后做 LTNE 热态 plug 定性/定量。
   ✅ **完成(2026-09-11)**:根因=AC 的 β 配置——`β=AC_beta·Uref²`,Uref 误取
   `LS_U_in=0.01`(压差驱动下仅作标度)使存档运行实为 β=0.2(日志 `beta=0.2`),而非
   文档假设的 β≈500;小 β 使 AC-Rusanov 连续方程数值耗散过大、粗网格稳态解偏移。
   plug AC `LS_U_in=0.01→0.5`(β≈500)后解≈理想(40×20:u_bar=7.958e-3、
   dp/dx=−3.1259e-2≈−Δp/L,5 172 步收敛);SIMPLE@40×20 自身距理想 −1.44%,故
   同网格两者差 ~1.4%;**80×40 同网格验收 AC vs SIMPLE:u_bar +0.20%、dp/dx −0.21%
   → <1% ✅**。记录:`cases/porous_plug_ac/hstudy/`(h20/h40/h80 收敛表+图)、
   `docs/工作日志.md` 2026-09-11 条目。
10. **beavers_joseph 接口(16)AC 双块复测**(阶段 2)✅/部分:粗网格
    (`cases/beavers_joseph/run_ac`,fluid 81×21 / porous 81×25)在 `LS_Algorithm=3` 下
    可收敛(res_m≈2.3e-4@3e5 外步),dp/dx 双块一致、Darcy 平台 −1.5%,剖面 L2≈3%;
    但界面滑移 u_i 偏高 ~15%(线性外推口径,单元层面偏高 ~6%),与 plug/darcy 的 AC
    速度偏高 +2.7%(待办 9)属同一系统性差,细网格收口待其解决后复核;另把“阻力主导
    下 AC 收敛慢于 SIMPLE”量化(本条目收敛曲线 `run_ac/residual_history.txt`)
    记入考核总结 §方法对比。
11. 补充实验教训:阻力主导工况需大 β(`AC_beta=2000,β≈U²·2000` 才在 ~1e5 步内收敛);
    建议后续加入阻力相关对角预条件 / β 自适应。
