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
12. **接口 19(可压↔多孔)分段交错 + AC 试跑**(2026-09-11):`cases/high_low_fluid_800K/
    run_ac_staggered` —— 在 4 块 code-19 网格(Block1 按 interface19_phaseA 改 POROUS)
    上用 `LS_Algorithm=3`(AC-porous)+ `Iflag_Couple_Scheme=1` 跑“1000 步暖机 3 轮后
    减半精修”交错调度。代码改动(向后兼容):`run_staggered_fluid_porous` 的 code-19
    配对搜索推广到全部 FLUID 块(原只查第一个 FLUID 块,此网格界面在 Block3);
    多孔段改经 `solver_one_block` 分派(LS_Algorithm=3 → porous AC)。12 轮外层调度/
    数据流正常、无 NaN;T_w 300→306–308 K、q_w 2.09e6→0.76e6 W/m² 单调收敛趋势,
    但未达 `Tol_Couple_Tw`(无骨架热沉+hv=0,属 staggered_tw README 已知限制)。
    AC-porous 段 20 000 伪步未达 AC_Tol(退化 LS 网格),每段跑满 AC_Max_Iter。
13. **分段交错推广到 11/12/13(AC)**:`Iflag_Couple_Scheme=1` 改为 `run_staggered_multiregion`
    自动按块类型分派:19(原 FLUID↔POROUS)、11(FLUID↔SOLID CHT)、13(LOWSPEED(AC)↔SOLID
    CHT)、12(FLUID↔LOWSPEED 匹配式块 GS)。新增:LS-AC 温度方程 `ac_lowspeed_energy`
    (LS_k>0,重建面质量流+`lowspeed_energy` GS)、判敛 `Tol_Couple_p/Tol_Couple_u`。
    - ✅ 12 试跑:`cases/high_low_fluid_800K/run_ac_staggered_12`(高速 300→减半与 LS-AC
      至收敛,界面 T/p/|u| 度量正常,无 NaN)。
    - ~ 11/13 调度路径试跑:`cases/fluid_solid/run_staggered_11`、`cases/bl_cht/
      run_staggered_ac13`(跑通但这两个自检网格的跨类界面以物理墙式 bc=2 表示,未触发
      couple_fluid_solid_interfaces 共轭分支;真共轭验收需 code-11/13 显式界面网格)。
    (下一步修改建议按优先级详见 `docs/工作日志.md` 2026-09-11 追加②末尾列表。)
14. **P0 完成(2026-09-11 追加③)**:code-11/13 真共轭验收采用“运行时界面提升”——
    `align_interface_to_bc_msg` 把“墙码 2+真实配对”的耦合面按块类型提升为接口码
    11/13(等),LS/AC 侧 `<0` 判断统一为 `is_interface_bc`;验收:code-11
    `fluid_solid/run_conj_11` 外轮 2 收敛(T_w≈800 K);code-13 `bl_cht/run_conj_ac13`
    (LS-AC)T_w 逐轮上升、q_w≈1e5 W/m² 共轭演化正常。README/日志已同步。

## 多块同类界面耦合收口（2026-09-15 新增，与 `docs/工作日志.md` 2026-09-15 条目联动）

15. ✅ **同类（LOWSPEED↔LOWSPEED）界面通量单值共享**：AC 界面原由两侧各自算通量，因
    Rusanov 耗散项 `−0.5·lam·(qR−qL)` 的 `lam` 对法向为**偶函数** ⇒ `F(−n) ≠ −F(n)`，
    差值恰为 `lam·Δq`；两侧传入的 (L,R) 状态顺序又相同（qA→qB）⇒ **通量不抵消**，每个
    界面有 O(1) 级虚假质量/动量源。已实现 owner 计算一次 + 对侧只读：
    `ac_build_ifc_registry` / `ac_ifc_find` / `ac_ifc_src_cell`（复用 `Umessage_send_mpi`
    的 L1/L2/L3+P 映射 ⇒ 旋转界面一致）/ `ac_ifc_flux_own` / `ac_ifc_flux_nb` /
    `ac_ifc_face` / `ac_interface_fluxes` + `FSH` 注册表；配对不上（id=0）自动退回原路径。
    **限制（2026-09-15 记录）**：注册表按**单个 mesh** 建立、切换 mesh 时重建（并清零
    `FSH`）；同类界面总在同一 mesh 内，故对现有算例（2/4 块槽道、BJ、bl_cht、fluid_solid）
    严格成立。若将来出现"同类型界面跨 mesh 交替求解"，应改为一次性遍历全部 mesh 建表、
    按 `ac_ifc_msh(id)` 取 owner（源码 `lows_ac_work` 中已有对应注释）。
    **验收**：`cases/channel_ac_2blk` 两块均 `res_q≈1e-7`（22570/18947 步）；
    `cases/channel_ac_4blk`（2×2 角点构型）blk1/2/3 `≈1e-7`、blk4 `8.3e-7`（撞 60000 内迭代上限）。
16. ✅ **耦合压力锚定**：`ac_have_pressure_bc` 改为**网格级**——串联同类低速块中任一有压力
    边界即视为全体有；否则只有速度入口的入口块每伪步 `ac_zero_mean_p` 会向界面**注入压力
    跳变**，典型症状“先收敛后发散”（blk1 res_q 1e-7 → 6.4e-2）。修复后入口块正常收敛。
    另注：内迭代数必须足以让 `AC_Tol` 判据（`iter>=20`）能触发，否则残差停在平台无意义。
17. ✅ **多块"残差 ~1e-7 但解非物理"定性（2026-09-15 收口）**：现象=邻块首列（紧贴界面
    的 `i=1`）出现反号速度（2 块 blk2 `u_min=−2.72`；4 块 blk4 `u_min=−14.4`、`u_max=2.11`）
    而残差 ~1e-7。**根因不是 VTK 输出**（`POINTS` 反推校验：只写内点、i 最快、块坐标正确；
    单块自校验 `u_max=1.49546` 正常），**也不是共享逻辑**（诊断打印证实邻块
    `map=(80,1,1)` 读到的 `FSH` 与 owner 写入**逐位相同**），而是**块间调度**：见 §20。
    已按 §20 重跑，非物理极值全部消失。
18. ~ **2D 三界面交汇验收（`cases/2D_T_tube`）**：**文件职责更正**——`bc3d_interface.inp`
    只服务**跨类型** block 对接；**同类块**的对接与物理边界全在 `bc3d.inp`（转换器
    `sub_convert_inp.f90::Convert_bc` 自动生成 `ib1..ke1/nb1/face1/L1..L3` 到 `bc3d.inc`），
    单一块类型时该文件不参与（缺省打印 `not found, skip interface`；实测删掉后
    `channel_ac_2blk` 仍解析出同类界面并正常收敛）。实测该算例**能加载能跑**
    （`run_probe.log`：4 块 61×61/61×101/61×101/101×61 节点、`Num_Cell=21600`、
    `shared same-class interfaces = 3`、100 外步无报错、能写 VTK），但 **`res_q≈3e2`、
    `res_m≈7.7e1`、`u_min≈−7.5`、数百个负速度单元 ⇒ 当前设置/几何未收口**（与接口文件无关）。
    待办：(i) 逐面核对 `bc3d.inp` 物理边界码与出口压力锚点；(ii) 检查 T 形交汇处 3 界面
    的角点处理；(iii) 调参（`AC_CFL=2/AC_CFLv=0.1` 远小于槽道 20；`AC_Flux=3`）；
    (iv) 文献数据需用户指定。**3D 方管算例仓库内不存在**。共形多块验收已完成：
    `channel_ac_2blk`/`channel_ac_4blk`（残差 1e-11/1e-14、`u_max`=1.49546/1.49503、
    界面两侧 Σu 差 0.12%、镜像严格对称 ⇒ **无角点病态**）。
19. [ ] **收窄 `Ts_send/recv_mpi` 过滤条件**：仍用 `is_interface_bc`（含 11–19），固/多孔
    温度缓冲跑在跨类共轭面上、与 `couple_*` 重复写 ghost；应按同类面（-1）过滤，并对
    `fluid_solid`(11) / `bl_cht`(13) / `beavers_joseph`(16) 回归后再定稿。

20. ✅ **多块 AC 必须"小内迭代 + 多外步"（同步调度）—— 已定案并写入算例**（2026-09-15）：
    `AC_Max_Iter` 大 + `t_end=1`（把某块一次迭代到收敛）会让该块求解期间**界面通量被冻结**：
    邻块（`i=1` 首列，紧贴界面）会被锁进**假解**——残差仍到 1e-7（该列"自洽假平衡"满足 max
    残差判据）但场量非物理（2 块 blk2 `u_min=−2.72`；4 块 blk4 `u_min=−14.4`,`u_max=2.11`）。
    正确用法：`AC_Max_Iter=20` + `t_end=3000`（各块每外步只推进 20 个伪步、界面每外步刷新）
    ⇒ 2 块 `res=4.8e-11/1.3e-11`、4 块 `res=8.5e-14…1.8e-14`，`u_max`=1.49439/1.49546（单块基准
    1.49546）、界面两侧 Σu 差 0.12%、出口 Σu=40.00000（解析）、**0 个负速度单元**。
    已更新 `cases/channel_ac_2blk|4blk` 的 `control.ec` 与 README，并在 `sub_lowspeed_ac.f90`
    保留默认关闭的界面诊断开关 `ac_dbg`（`ac_dbg=1` 打印 owner/neighbour 的 `FSH` 读写槽、
    两侧状态与 `(mBlock,ksub)->(id,is_owner)` 解析，用于此类问题快速定位）。
    **验收提示**：多块不能只看 max 残差，必须同时检查"界面两侧剖面/通量一致性 + 场量极值"。
    同类多块问题的排查顺序建议：①VTK 口径（`POINTS` 反推单元中心自校验）；②界面槽映射
    （`ac_dbg`）；③调度（先只改 `AC_Max_Iter/t_end` 做对照）。

## 圆柱绕流（半域+对称）标模：AC 阻力 vs Rogers 与无粘 Cp 解析验证（2026-09-16 新增）

21. ✅ **单块外流标模 `cases/cylinder_re40_half`**（上半圆柱 + 对称面 y=0，O 型贴体网格，
    `LS_Algorithm=3`）。粘性 Re=40：加密后 `Cd=1.5536`（壁面 2 阶积分）/`1.5582`（远场动量
    平衡），与 Rogers & Kwak（NASA TM-101051 Table 4）的 **1.549** 差 +0.3%、与
    Dennis & Chang 的 1.522 差 +2.1%；`Cd_p=1.0086`（文献 1.011/0.998）、`Cp_front=1.159`
    （文献 1.147/1.144）、分离角 125.8°/后驻点 54.2°（文献 53.0~54.8°）。
    **新增两种独立阻力口径**（壁面积分 + 远场圆弧动量平衡，且控制面半径自洽）写进
    `postprocess.py`，可直接复用于其它外流算例。
22. ✅ **无粘（`If_viscous=0`）逐点 Cp 解析验证**（`verify_cp_inviscid.py`）：解析
    `Cp=1-4sin²θ`（θ 自后驻点）。前/上半面 θ≳130° 一致到 0.5~3%、吸力峰 +0.10(3.3%)、
    前驻点 +2.6% ⇒ 壁面压力/通量/远场/对称 BC 正确；**但后驻点区（θ≲20°）Cp 误差达 −0.61**
    （数值"死水区"），并伴随 **首层单元总压偏低（H=0.185 vs 0.502，首层速度 −8%）**，
    表现为虚假阻力 `Cd≈0.02~0.04`（d'Alembert 被破坏）。⇒ 现有 MUSCL2+AC-Rusanov 的主要
    弱点在"**强减速/回流区 + 近壁层**"，与 T 形管交汇处同类，是后续改进（WENO5/限幅器/
    近壁重构/特征远场）的验收靶子。
23. ⚠ **暴露的三个鲁棒性问题（待修，均与多块/接口无关）**：
    (a) `AC_Flux=3`(AUSM+) 在无粘圆柱上 2000 步内 NaN；改用 `AC_Flux=1` 正常；
    (b) `LS_Algorithm=1`（低速 SIMPLE）同网格 **2 步内发散**（弧形进出口 + 壁/对称相交角点）；
    (c) 远场距离敏感：R_out 25D→10D 时 Cd 由 1.57 涨到 1.67（+6%），提示远场 BC 反射
    （Rogers 用特征型无反射 BC）——与文献同域对比前需先实现特征远场。
24. 🔧 代码（向后兼容）：`ac_set_ghost_cell` 中 `BC_Wall + If_viscous=0` 改为**滑移壁**
    （只反转法向分量），与可压求解器同一约定；`If_viscous=1` 分支逐位不变。


## 近壁/强减速区精度改进（2026-09-17 新增，与 `docs/工作日志.md` 同日条目联动）

25. ✅ **壁面压力二阶重构 `AC_WallP`（默认 1）**：壁/对称面压力由"取相邻单元值"（法向零梯度，
    隐含 `∂q/∂n=0`）改为通过面做内部场线性重构 `q_wall=1.5q₁−0.5q₂`（新例程
    `ac_wall_pressure` + `ac_wallface` 面映射）。平壁/边界层 `∂q/∂n=0` 时与旧式等价；曲壁上
    保留离心项 `∂q/∂n=u_t²/R`。**效果（无粘圆柱）**：max|ΔCp| 0.6106→**0.5413**（−11%）、
    RMS 0.2129→0.2072、前驻点 ΔCp +0.026→**−0.003**、表面速度比 0.7850→0.8071、Cd_p 0.0393→0.0364；
    **代价**：吸力峰（θ≈80–90°）ΔCp +0.097→+0.168（该处主因是首层速度亏损 −8%，与壁压无关）。
    粘性圆柱 Re=40：Cd 1.5706→**1.5679**、Cd_p 1.0152→1.0128、Cd_f/Cp_stag 不变或略优；
    方腔 u-RMS 0.02105→**0.02056**。`AC_WallP=0` 可退回旧式。
26. ✅ **壁/对称面切向状态单侧二阶重构 `AC_WallRecon`（默认 1）**：镜像 ghost 使切向 MUSCL 斜率
    恒为 0 ⇒ 壁面侧状态退化为一阶（原缺陷）；改为单侧内点线性外推 `a_eff=2c−d`，法向分量仍用
    镜像 ghost（无穿透性/质量通量不变），新增 `ac_tang_merge`。**效果**：无粘圆柱中性；
    槽道 u_max 1.495464→**1.495718**、2 块槽道 1.49439/1.49546→**1.49465/1.49571**（更接近解析 1.5），
    残差与界面一致性不变。对切向梯度大的壁面/角点（T 形管交汇、方管、圆柱尾驻点）意义最大。
27. ⚠ **Rusanov 耗散逐场分裂 `AC_MomDiss`（默认 0，可选 1/2；`AC_MomFrac` 默认 0.2）**：AC 特征值
    {u_n,u_n,u_n±c}，质量方程保留 `|u_n|+c`（压力耦合/防棋盘阻尼），`=1` 时动量方程改用
    `max(|u_n|)+AC_MomFrac·max(c)`。**动机（已量化）**：原实现把 `λ=max(|u_n|+c)` 统一施加于四方程，
    低马赫下给动量额外加了 ~√β·Δx 的人工粘性；无粘圆柱 0<θ<18° 的壁面 Cp 平台（+0.39 vs 解析 +1.0）、
    对称线后方 r=1.7R/2.1R 处 u=0.17/0.46（解析 0.67/0.78）即"数值死水区"表现。
    **实测结论（重要，决定默认值）**：`=1` 使**壁面主导/粘性内流明显改善**——槽道 u_max 误差
    0.285%→0.153%、剖面 RMS 0.184%→0.130%（收敛步数不变 39377）；粘性圆柱 Re=40 Cd
    1.5679→**1.5589**（Rogers 1.549）、Cd_p 1.0128→1.0066——但**无粘欧拉解变差**
    （max|ΔCp| 0.5413→0.8712 @Frac=0.2、0.6209 @Frac=0.5），因为无粘时数值耗散是唯一的稳定机制。
    因此**默认保持 0**（既有验证结果不变），按算例开启 1 使用。纯 `|u_n|`（`=2`）**不稳定**
    （驻点处 u_n→0 使耗散消失，无粘圆柱 ~15000 伪步发散），已作为负结果记录在代码注释中。
    槽道残差轨迹与旧式逐点重合 ⇒ **`=1` 无收敛代价**。
28. ✅ **粘性面距离按面法向投影**（`ac_viscous_res`/`dist_ac`）：斜/畸变网格上按指标线距离会高估
    法向导数 1/cosθ，现投影到面单位法向（`B%ni/nj/nk`）。正交网格恒等 ⇒ 槽道/方腔/圆柱逐位不变。
    另注：壁面 `d=0.5|x₂−x₀|` 在"节点镜像"约定下恰为 2×(壁到首层中心距)，故壁面剪切原本已是二阶。
29. ✅ **移除 `AC_Flux=3`（AUSM+）+ `ac_check_controls()` 校验**：原因 (a) AC 无物理声速，AUSM+ 的
    马赫分裂与 LU-SGS 的 `c=√(u_n²+β)` 口径不一致；(b) 实现缺 AUSM+-up 的低马赫压力耗散项 p_u；
    (c) 实测无粘圆柱 2000 步 NaN。现仅支持 1(Rusanov)/2(Steger–Warming)，越界值告警回落为 1；
    `AC_Recon/AC_Limiter/AC_WenoBlend/AC_WallRecon/AC_WallP/AC_MomDiss/AC_MomFrac` 同步校验；
    启动打印完整 AC 控制量。`sub_porous_ac.f90` 的 AC 路径同步接线（新增 `ac_wall_pressure` 调用与新
    数组分配；`ac_wall_pressure` 内部对未分配数组做提前返回保护）。
30. [ ] **遗留（下一步）**：(i) 吸力峰首层速度亏损 −8%（及由此导致的 Cl 偏差）需靠**近壁切向
    重构/曲线壁面处理**继续改（当前壁压手段已用尽）；(ii) `AC_MomFrac` 的**最优值扫描**（0.1/0.2/0.3/0.5）
    与"死水区"指标（尾部 Cp 平台、对称线速度剖面、Cd_inviscid）联动；(iii) WENO5/WENO3 在 AC 共置格式
    下仍不稳（WENO5+blend0.1 → NaN；WENO3 收敛慢），需要真正的低马赫限幅/混合策略，不建议直接使用；
    (iv) 特征型远场 BC（消除 R_out=10D 的 +6% 反射）。


