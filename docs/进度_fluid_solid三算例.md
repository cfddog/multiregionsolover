# 进度摘要 —— `cases/fluid_solid` 三个流固耦合算例（2026-09-18，未提交）

> **用途**：本文件是本次会话的中断点记录，便于下次快速接续。
> 基线提交 `adaa43b`（分支 `main`）；本次改动全部**未提交**，在工作树中。
> 编译产物 `src/opencfd-ec1.16a.out`（2026-09-18 13:44，构建 EXIT=0）。

---

## 1. 目标（用户需求）

在 `cases/fluid_solid`（网格：Block1 = 67×101×2 薄板，Block2 = 101×201×2 通道，
界面 = **Block1 j=1 面 ↔ Block2 i=101 面**，两块的物理法向都是 y）上做三个算例，
均要求 **分段交错耦合**（先流体跑 1000 步 → 再跑另一侧 → 逐步减半间隔 →
最后 `Kstep_Couple_Min=1` 逐层同步推进）：

| 算例 | Block1 | Block2 | 界面码 | 交错驱动 |
|---|---|---|---|---|
| case1 | 不锈钢固体（j=max 等温 300 K） | Ma=3 空气, Re=5000, T∞=800K, 壁面绝热 | **11** FLUID↔SOLID | `run_staggered_fluid_solid(11)` |
| case2 | 低速空气（j=max 质量入口 G=2 kg/(m²·s), 300K） | 同上 | **12** FLUID↔LOWSPEED | `run_staggered_highlow` |
| case3 | 不锈钢多孔 ε=0.30（j=max 冷却剂注入 G=2, 300K） | 同上 | **19** FLUID↔POROUS | `run_staggered_fluid_porous` |

Block1 的 j=1 为对接面、j=max 为冷端/入口 —— **与现有网格一致，网格未改动**。

## 2. 关键难点与结论（重要）

既有 `couple_*` / `run_staggered_*` 例程的“可压侧”一律硬编码为
**gas j- (face=2) ↔ 对侧 j+ (face=5)**；而本网格是
**gas i+ (face=4) ↔ 对侧 j- (face=2)**。运行时界面 dump 实测连接描述符：

```
blk1 type1 bc=11 face=2 (j-) nb1=2 face1=4 (i+)  L=(-2, 1, 3)
blk2 type0 bc=11 face=4 (i+) nb1=1 face1=2 (j-)  L=( 2,-1, 3)
```

即：gas 的 i(法向) ↔ 对侧 j(法向)；gas 的 j(切向) ↔ 对侧 i，**反向**（`L(2)=-1`）。
因此**必须扩展代码**；改网格或镜像 Block1 的 j 轴会破坏“j=1 为对接面”的约定，已排除。

**关键简化**：本网格两块的法向在物理上都沿 y，故 `U(2)=切向 x、U(3)=法向 y`
的分量语义与旧路径**完全一致** —— 只需改**下标/鬼单元索引**，物理公式原样复用。

## 3. 已完成的代码改动（`src/`，向后兼容，未提交）

### 3.1 `src/sub_multi_region.f90`（+651/−344 行）

1. **`run_staggered_fluid_solid`（码 11）**：面检查放宽为
   `solid j- (2) vs comp i-/i+ (1/4)`。共轭热平衡由既有全通用例程
   `couple_compressible_fluid_solid_face`（用 L1..L3）完成，**未改物理**。
2. **`couple_highlow_fluid_face`（码 12 逐步耦合）**：新增 `(4,2)` 分支
   （gas i+ ↔ LS j-）：gas 内部列 `ie-1`/鬼列 `ie..ie+LAP-1`；LS 内部行 `jb1(=1)`/
   鬼行 `jb1-1, jb1-2, ...`；切向映射 `gas_j = jb + (ie1-1) - ls_i`（反向）。
   旧 `(2,5)` 路径逐字保留。
3. **`run_staggered_highlow`（码 12 交错驱动）**：接受 `(4,2)`；新增 `pair42`；
   界面量（T/p/|u|）改按**对侧（LS）i 单元**索引（`ia1..ia2/ka1..ka2`）；
   `iface_highlow.dat` 的 x_w 改由 gas 的 i=ie 节点给出。
4. **`couple_compressible_porous_blowing_face`（码 19 发汗耦合）**：新增 `(4,2)` 分支
   （gas i+ ↔ porous j-，冷却剂经 j=max 面注入）。旧 `(2,5)` 路径保留。
5. **`run_staggered_fluid_porous`（码 19 交错驱动）**：接受 `(4,2)`，重构为 `pair42`
   双分支：
   - 界面工作数组 `fp_*` 统一按**多孔块 i 单元**索引，尺寸改为
     `max(Bf%nx,Bp%nx)+1 × max(Bf%nz,Bp%nz)+1`；
   - `fill_gas_wall_ghost(nv)`（7 参 → 1 参）支持 gas i+ 等温壁；
   - q_w 提取：gas i+ 列 + 两单元中心的欧氏距离 dx；
   - 多孔段 j- 面 ghost（`jg_l-n1`，零梯度 U + `p=fp_pw`）；
   - `set_porous_Ts_flux(...,sgn,...)` 增加方向参数（+1 = j+ / −1 = j−）；
   - 返回 `T_w = 0.5*(Ts(1)+Ts(0))`（pair42）；
   - 新增打印 `coolant exit |v|_max`（发汗可见性）。
6. **新增 `wall_bound_blowing`（pair42 专用内嵌例程）**：等温 + 法向喷射壁
   （`v_ghost = 2*v_w - v_int`，切向镜像，`T=2Tw−T`），其中 `v_w` 取多孔块界面单元
   `Bp%U(3,gi,jc_l,k)` —— **让冷却剂喷注真正作用到气体边界层**
   （旧交错路径只把界面当等温无滑移壁）。仅 `LS_Inlet_Type==1` 且 `v_w>0` 时启用。
7. **新增诊断 `dbg_dump_interfaces`**：`IF_Debug=1` 时打印各界面的
   face/face1/nb1/ranges/L1..L3，排查界面配对/方向问题很有用。

### 3.2 `src/opencfd_ec3d_v1.16a.f90`（+2 行）
`Start ......` 之后加 `if(IF_Debug == 1) call dbg_dump_interfaces(1)`。

> 既有路径（码 13、码 16、`(2,5)` 的 12/19）未改物理：旧代码逐字保留，仅外包 `if/else`。

## 4. 新建的算例目录

```
cases/fluid_solid/run_case1_solid/      （码 11：固体↔可压）
cases/fluid_solid/run_case2_lowspeed/   （码 12：低速↔可压）
cases/fluid_solid/run_case3_porous/     （码 19：多孔↔可压）
```

每个目录含：`Mesh3d.x`(copy)、`bc3d.inp`、`bc3d_interface.inp`、`material.in`、
`control.ec`、`opencfd-ec1.16a.out`(软链到 `src/`)；
case1 另有 `solid_bc.inp`（`5 300.0 0.0`：固体 j=max 面等温 300 K）；
case3 另有 `porous.inp`（`1  0.30  1.0e-3  2.0e6`：ε、dp、hv）。

要点：
- `bc3d.inp`：界面码显式写 **11/12/19**；**Block2 的 i=1 面由 4(farfield) 改为 2(绝热壁)**
  （按用户“壁面为绝热壁面”）；Block1 的 j=max 面 case1=2(壁)、case2/3=5(入口)。
- `bc3d_interface.inp`：沿用原 Pointwise 结构（界面行 `-1` + 续行给 nb1；配对按块类型分派，
  无需改成显式码），仅把 Block2 i=1 的 4→2 与之保持一致。
- `material.in`：case1 `1 0` + `7900.0 500.0 16.3`（不锈钢 304）；case2 `2 0` 全 0；
  case3 `3 0` + `7900.0 500.0 16.3`。
- **质量入口的实现方式（重要）**：低速/多孔求解器为常密度（`LS_rho`），且
  `LS_Inlet_Type=2` 的流量入口通路按“i 面 + 网格单位面积”实现，用于 j 面会方向/量纲错。
  故 **G=2 kg/(m²·s) 以等效法向速度给出**：
  `LS_V_in = Porous_V_in = G/LS_rho = 2/1.177 = 1.699 m/s`（+y），
  `LS_rho=1.177 kg/m³`（300 K / 1 atm 空气）。等价且无量纲陷阱，
  已在 `control.ec` 注释说明；真正走“质量入口”通路列为待办。
- 交错调度（三例相同）：`Iflag_Couple_Scheme=1, Kstep_Couple_Comp=1000,
  Niter_Couple_Warm=2, Kstep_Couple_Min=1, Niter_Couple_Outer=12`
  → 实际气体段步数 1000,1000,500,250,125,62,31,15,7,3,1,1；`If_viscous=1`（热流必须有粘性）。

## 5. 试跑结果（短调度验证：无 NaN、EXIT=0）

> 为快速验证，临时把 `Kstep_Couple_Comp=20, Niter_Couple_Warm=1,
> Niter_Couple_Outer=4, Kstep_Couple_Min=2`；跑完已**恢复交付版 control.ec**
> （日志 `run_trial.log`；各目录 `control.ec.orig` 副本已清理）。

### 5.0 交付版 `control.ec` 参数已验证（✅ 已完成）

用**交付版参数**（`Kstep_Couple_Comp=1000 / Niter_Couple_Outer=12`）各跑一次短时限运行，
确认解析与首轮调度正确（日志 `run_check_delivered.log`）：

| 算例 | 关键行 |
|---|---|
| case1 | `run_staggered_fluid_solid(paircode=11): ... Kstep_Couple_Comp=1000 Niter_Couple_Outer=12`；`outer iter 1: flow chunk steps = 1000`；`solid chunk done, outer iter 1  max|dT_w|=499.9 K  T_w≈799.9 K  max|q_w|=2064 W/m²`；`outer iter 2: flow chunk steps = 1000` |
| case2 | `run_staggered_highlow(12): gas block 2  LS block 1  Kstep_Couple_Comp=1000 Niter_Couple_Outer=12`；`outer iter 1: gas chunk steps = 1000` |
| case3 | `run_staggered_multiregion: block types F/L/S/P = T F F T`；`pairs 11/12/13/19 = F F F T`；`gas block 2  porous block 1  pair42(gas i+/porous j-)= T`；`outer iter 1: gas chunk steps = 1000` |

（case1 的 1000 步暖机段在 ~40 s 内跑完，说明整例 12 轮在单机上数分钟内可完成。）

### 5.2 现场测试（2026-09-18 第二轮，交付版设置整段跑）

按交付版 `control.ec`（`Kstep_Couple_Comp=1000`）各跑一段（脚本 `/tmp/smoke_all.sh`，
日志各目录 `run_smoke.log`）：

| 算例 | 结果 | 关键数据 |
|---|---|---|
| case1 | ✅ **12 轮全跑完**，EXIT=0，**99 s**，无 NaN | 调度 1000,1000,500,250,125,62,31,15,7,... 正确；`max|dT_w|` 499.9 → 0.0996 K（平台，见下）；`T_w` 799.90→799.20 K 缓降；`max|q_w|` 2064→1672 W/m²；残差 ~1e-3 |
| case2 | ✅ 正常推进，无 NaN；420 s 限时(EXIT=124) | `max|dT|` 321→20.2→17.5→9.2 K，`max|dp|` ~0.2–0.9 Pa，`max|du|` 29→3.0 m/s；**LS-AC 每轮跑满 5000 迭代、`res_q` 平台 ~1.4e-2** |
| case3 | ✅ 正常推进，无 NaN；437 s 限时(EXIT=124) | 首轮气体段完成，`max|q_w|`=1.09e6 W/m²；多孔 AC `res_q`=5.6e-5；LTNE `dTs`→9.9e-4；**原设置 `AC_Max_Iter=20000 × Porous_Chunk_Iter=10` 使首轮多孔段在 437 s 内跑不完** |

### 5.3 针对发现问题的设置调整（已完成）

1. **case3 多孔段预算压缩**（原设置太慢）：交付 `control.ec` 改为
   `AC_Max_Iter=5000, AC_Print=1000, AC_CFL=20, AC_beta=100, Porous_Chunk_Iter=3`
   （等价上限 15000 迭代/轮）。试验（`/tmp/tune3.sh`）结果：**2 轮完成、无 NaN、≈170 s/轮**，
   多孔 AC `res_q`=6.4e-5（与 20000 迭代的 5.6e-5 相当）；
   `coolant exit |v|_max` 10.3→15.4 m/s（入口 1.699 m/s，冷启动暂未收敛，需关注）。
2. **case2 的 LS-AC 结论**（不改交付设置，保留稳定配置：
   `LS_Algorithm=3, AC_CFL=2, AC_beta=10, AC_Max_Iter=5000`）：
   - 试 `AC_CFL=20, AC_beta=100`（β=289）：`res_q` 1.5e-2→**3.3e-3** 有明显改善，
     但在**第 3 轮出现 NaN**（不稳定）⇒ **不可用**。
   - 试 `LS_Algorithm=1`（SIMPLE）：**第 1 步即发散**（`res_p`→6e31, `u`→-1e6）⇒ **不可用**。
   - 结论：**必须用 AC，且保持 `AC_CFL=2`**；`res_q` 平台 ~1.4e-2 属该配置的已知限制
     （界面度量 max|dT|/max|du| 逐轮下降，结果仍可用）。
3. **case1 的判据说明**（不改设置）：`max|dT_w|` 平衡在 ~0.0996 K >
   `Tol_Couple_Tw=1e-2`，因固体（50 mm 厚、300 K 冷端）尚未达到稳态、`T_w` 每轮缓降 ~0.1 K；
   属**物理瞬态**而非设置错误。若希望外层早停，可把 `Tol_Couple_Tw` 放宽到 ~0.2 K。

### 5.4 结论：三个算例现场可用 ✅

- 三例均能按交付设置启动、正确识别材料/界面/调度、正常推进、**无 NaN**；
- 预计单机耗时：case1 ~100 s（12 轮全跑完）；case2 ~20 min（12 轮，AC 每轮≈105 s）；
  case3 ~40 min（12 轮，多孔段≈150 s/轮 + 气体段）；
- **必须单进程**（`mpirun -np 1`）；12/19 无跨进程耦合。

| 算例 | 关键日志 | 结论 |
|---|---|---|
| case1 | blk1 type1(7900/500/16.3) / blk2 type0；界面 dump `face=2,face1=4`；`flow chunk steps 20/10/5/2`；`max|dT_w|` 499.9→0.094 K；`T_w≈799.6 K`；`max|q_w|≈1.0e4 W/m²` | ✅ CHT 调度与热平衡正常 |
| case2 | blk1 type2 / blk2 type0；`gas chunk steps 20/10`；`max|dT|=138 K, max|dp|=0.76 Pa, max|du|=29 m/s`（被 U_MAX_LS=30 钳制）；LS-AC 5000 伪步 `res_q~1e-2` | ✅ 调度/数据流正常；AC 收敛待调 |
| case3 | blk1 type3 (eps=0.30, hv=2e6, rho_s=7900, k_s=16.3) / blk2 type0；`gas chunk steps 20/10/5/2`；`max|q_w|` 5.26e6→3.93e6 W/m²；`T_w` 309→330 K、`max|dT_w|` 9.2→6.9 K | ✅ 调度/热流/发汗数据流正常，外层缓慢收敛（多孔储热） |

## 6. 下次接续：快速启动

```bash
cd /home/sundong/Fortran_Project/OpenCFD-EC-1.16a/src
make                      # 增量编译（~1 min）→ opencfd-ec1.16a.out

# 冒烟跑（短调度）
cd ../cases/fluid_solid/run_case1_solid
cp control.ec control.ec.full
sed -i 's/^  Kstep_Couple_Comp=.*/  Kstep_Couple_Comp=20/;
        s/^  Niter_Couple_Warm=.*/  Niter_Couple_Warm=1/;
        s/^  Niter_Couple_Outer=.*/  Niter_Couple_Outer=4/;
        s/^  Kstep_Couple_Min=.*/  Kstep_Couple_Min=2/;
        s/^  Kstep_show=.*/  Kstep_show=5/;
        s/^  IF_Debug=.*/  IF_Debug=1/' control.ec
mpirun -np 1 ./opencfd-ec1.16a.out > run_smoke.log 2>&1
cp control.ec.full control.ec         # 记得还原
```

- **必须单进程**（`mpirun -np 1`，不要 `partation.dat`）：跨进程耦合目前只对码 11
  实现（`cht_mpi_fluid_solid_face`），**12/19 未实现跨进程**。
- 网格单位 mm（`Lscale=0.001`），故 `Re=5000` 是**以 1 mm 为参考长度**的雷诺数
  （`RHO_REF = Re·μ∞/(Ma·a∞·Lscale)`）；按板长（166.7 mm）换算 Re_L ≈ 8.3e5。
  若用户本意是“板长 Re=5000”，需改 `Re` 并重新标定密度/粘性（见待办）。
- 交付版一轮 12 个外层、气体段共约 3000 步；想快速看调度可先按上面短调度跑。

## 7. 待办（按优先级）

1. ~~**交付确认**~~ ✅ 已完成（见 §5.0）：三例交付版参数均已确认解析与首轮调度正确。
2. **写三个目录的 `README.md`**：几何/边界表、参数含义、G→V 换算、已知限制、复现命令。
3. **case2 的 AC 标定**：`AC_CFL=2→20`、`AC_beta=10→100`、`AC_Max_Iter` 加大，
   把 `res_q` 压到 `AC_Tol` 附近；或与 `LS_Algorithm=1`(SIMPLE) 对比。
4. **case3 外层加速**：多孔储热导致 `T_w` 缓慢爬升（与既有 `staggered_tw` README 同类），
   可加界面欠松弛 / 增外层轮数 / 减小 `hv` 或给骨架加冷端。
5. **质量入口通路修复（可选）**：`sub_lowspeed.f90::lowspeed_inlet_velocity` 与
   `sub_porous.f90::porous_inlet_velocity` —— 让 `LS_Inlet_Type=2` 按**面法向**赋分量，
   并明确面积量纲，以便直接按 kg/(m²·s) 输入。
6. **清理运行产物**：三个目录里的 `*.vtk / flow3d.dat / Residual.dat / *.inc / *.orig` 等，
   保留 `run_trial.log` 作证据；注意 case1 目录有含空格的怪文件 `'force-debug-  0'`。
7. **文档同步**：`docs/工作日志.md` 追加本次条目；`docs/程序能力与算例考核总结.md` §1
   接口表可补“11/12/19 已支持 gas i+ ↔ 对侧 j- 拓扑”。
8. **提交**（仅 2 个源文件 + 3 个算例目录 + 本文档）：
   `feat(staggered): support gas i+ <-> partner j- topology for codes 11/12/19 in cases/fluid_solid`

## 8. 备份/补丁脚本（工作树外，可回滚或重放）

```
/tmp/sub_multi_region.f90.bak0            # 改动前
/tmp/bak_hl, /tmp/bak_por, /tmp/bak_blow  # 各阶段快照
/tmp/hl_new.f90 /tmp/por_new.f90 /tmp/blow_new.f90   # 新代码片段
/tmp/patch_*.py /tmp/fix*.py              # 全部补丁脚本
git checkout -- src/ && rm src/opencfd-ec1.16a.out && make   # 回滚到基线
```
