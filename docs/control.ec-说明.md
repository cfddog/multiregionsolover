# control.ec 分组说明（7 组 namelist）

> 变更时间：2026-09-23　代码位置：`src/sub_read_parameter.f90`、`src/sub_modules.f90`
> 模板文件：`control.ec.template`（仓库根）

## 1. 为什么要改

`$control_ec` 一个组里塞了 **128 个变量**，而**全部 128 个都有代码默认值**。
实测几个代表算例：文件里 70~75% 的赋值行与默认值完全相同（可整行删除）。

| 算例 | 行数 | 赋值条数 | == 默认值 | 必须保留 |
|---|---|---|---|---|
| fluid_solid/run_case1_solid | 121 | 112 | 75 | 33 |
| fluid_solid/run_case2_lowspeed | 128 | 114 | 75 | 35 |
| fluid_solid/run_case3_porous | 129 | 113 | 69 | 40 |
| channel_ac | 93 | 91 | 65 | 22 |
| lidcavity | 89 | 87 | 66 | 17 |

因此：**按物理/功能域分成 7 组**，且**每组可选、每个变量可选**，算例只需写与默认值不同的项。

## 2. 七个组（变量归属，一个不漏）

| 组 | 变量数 | 用途 | 变量 |
|---|---|---|---|
| `$freestream_ec` | 24 | 来流条件 / 参考态 / 气体物性 | `Ma, Re, AoA, AoS, p_outlet, gamma, PrL, PrT, T_inf, Twall, Lscale, Ref_S, Ref_L, Centroid, Cood_Y_UP, Kt_inf, Wt_inf, IF_TurboMachinary, Ref_medium_usrdef, Turbo_P0, Turbo_T0, Turbo_L0, Turbo_w, Turbo_Periodic_seta` |
| `$flow_ec` | 49 (+3) | 流场控制：时间步/格式/模型/网格并行/限制器/输出/重启 | `t_end, Kstep_save, Kstep_show, Kstep_average, Kstep_smooth, Kstep_init_smooth, CFL, dt_global, dtmax, dtmin, Iflag_local_dt, Time_Method, Step_Inner_Limit, Res_Inner_Limit, w_LU, If_Residual_smoothing, If_dtime_mesh, Iflag_Scheme, Iflag_Flux, IFlag_Reconstruction, Bound_Scheme, IF_Scheme_Positivity, IFLAG_LIMIT_FLOW, Iflag_turbulence_model, If_viscous, Iflag_init, MUT_MAX, CP1_NSA, CP2_NSA, Ldmin, Ldmax, Lpmin, Lpmax, Lumax, LSAmax, Mesh_File_Format, Num_Mesh, Pre_Step_Mesh, NUM_THREADS, IF_Debug, Pdebug, Periodic_dX, Periodic_dY, Periodic_dZ, IF_Innerflow, Iflag_savefile, Iflag_vtk_onefile, Iflag_vtk_SI, Iflag_bc_check, **Iflag_restart, Kstep_restart, Iflag_flow_node**` |
| `$lowspeed_ec` | 21 | 低速（不可压）块 | `LS_rho, LS_mu, LS_k, LS_Cp, LS_T_ref, LS_Inlet_Type, LS_U_in, LS_V_in, LS_W_in, LS_Mdot_in, LS_P_in, LS_P_out, LS_T_wall, LS_U_lid, LS_alpha_p, LS_alpha_u, LS_alpha_T, LS_Max_Iter, LS_Tol, LS_Scheme, LS_Algorithm` |
| `$ac_ec` | 15 | 人工压缩(AC)求解器（低速与多孔共用） | `AC_Max_Iter, AC_Print, AC_beta, AC_CFL, AC_CFLv, AC_Tol, AC_w, AC_Flux, AC_Recon, AC_Limiter, AC_WenoBlend, AC_WallRecon, AC_WallP, AC_MomDiss, AC_MomFrac` |
| `$solid_ec` | 4 | 固体导热块（**新增**：原为代码内硬编码） | `Solid_GS_Omega=1.7, Solid_Max_Iter=20000, Solid_Min_Iter=5, Solid_Tol=1e-9` |
| `$porous_ec` | 8 | 多孔介质块 | `Porous_T_ref, Porous_alpha_Ts, Porous_Max_Iter, Porous_Tol, Porous_U_in, Porous_V_in, Porous_W_in, Porous_T_in` |
| `$couple_ec` | 12 | 跨区域（分段交错）耦合调度 | `Iflag_Couple_Scheme, Kstep_Couple_Comp, Niter_Couple_Outer, Niter_Couple_Warm, Kstep_Couple_Min, Twall_Couple_Init, Tol_Couple_Tw, Tol_Couple_p, Tol_Couple_u, Iflag_Couple_WallFlux, Porous_Chunk_Iter, Iflag_Couple_Restart` |
| 旧组 `$control_ec` | 128(+4) | 历史单组，**完全保留** | 全部变量，`Solid_*` 追加在末尾 |

合计：现有 128 个变量全部覆盖，另加 4 个固体控制量。

**不在 control.ec 里**（仍是独立文件，因为是逐块/逐面数据）：
- `material.in`：块类型表 + 每块 `(rho, Cp, k)`
- `solid_bc.inp`：固体/多孔逐面热边界 `Tw / Qw / htc`
- `porous.inp`：多孔块 `eps, dp, hv`
- `bc3d.inp` / `bc3d_interface.inp`：边界与块间对接

## 3. 读取规则（兼容 + 严格）

1. 先**扫描** `control.ec`，判定实际出现了哪些组（大小写不敏感；`$` 与 `&` 均可；
   `!` 之后视为注释；**只有“首个非空字符为 `$`/`&`”的行才算组头**，
   所以正文/注释里提到组名不会误判）。
2. 若有 `$control_ec` → **按原样读**（119 个既有算例零改动，行为逐位不变）。
3. 再按 `freestream → flow → lowspeed → ac → solid → porous → couple` 顺序读实际存在的组；
   **分组值覆盖旧组值**（支持“一次迁移一组”的渐进迁移）。
4. 某组**存在但解析失败**（变量名未知/拼写错、缺少 `$end` 或 `/`）⇒ **立即报错并 `stop 1`**，
   打印组名与 `iostat`，**绝不静默回落到默认值**。
5. 启动时终端与 `output_para.out` 末尾会回显各组来源（`T`=来自文件，`F`=用默认值），
   并打印固体 GS 四个控制量，便于回归比对。

## 4. 新算例怎么写（示例）

```ec
$freestream_ec
  Ma=3.0d0, Re=5000.d0, T_inf=800.d0, Lscale=0.001d0
$end
$flow_ec
  t_end=1.0d6, Kstep_save=500, Kstep_show=50,
  Mesh_File_Format=2, dt_global=1.0d0, dtmax=1000.d0,
  IFLAG_LIMIT_FLOW=0, Iflag_vtk_onefile=1, Iflag_vtk_SI=1
$end
$couple_ec
  Iflag_Couple_Scheme=1, Niter_Couple_Warm=2, Niter_Couple_Outer=12, Tol_Couple_Tw=1.0d-2
$end
```
上面 15 行即可替代原来 121 行的 `run_case1_solid/control.ec`（值完全等价）。
写新算例时直接从 `control.ec.template` 拷贝，删掉不需要的组即可。

## 5. 迁移建议（本次不强制）

本次只做**分组**，`cases/**/control.ec` 一律不改（靠双读兼容）。若要精简某个算例：

1. 只删除“与默认值相同”的行 —— 注意像 `LS_Max_Iter=2500`（默认 5000）、
   `AC_Max_Iter=5000`（默认 40000）、`Niter_Couple_Outer=12`（默认 60）、
   `Tol_Couple_Tw=1e-2`（默认 2e-2）、`Kstep_save=500`（默认 1000）**看着像废话其实必须留**。
2. 用 `output_para.out` 做**逐项 diff** 验证：新旧两次运行的 `output_para.out`
   除末尾 4 行“组来源回显”外必须完全一致。

## 6. 同步改动

| 文件 | 改动 |
|---|---|
| `src/sub_modules.f90` | 新增 `Solid_GS_Omega/Solid_Max_Iter/Solid_Min_Iter/Solid_Tol` |
| `src/sub_read_parameter.f90` | 1 组 → 7 组 namelist；新增扫描/严格报错/组来源回显；新增 4 个固体默认值 |
| `src/sub_multi_region.f90` | `solid_solver_one_block` 的硬编码常量改为读 `$solid_ec`（默认值与原值相同） |
| `util/readflow3d-ver2.5.f90`、`util/readflow3d-ver2.4a.f90` | 同步双读；并补齐 `Lscale` 与 `LS_*/AC_*/Porous_*/couple` 声明（原先这两个后处理工具**读不了**含这些变量的算例文件，属历史遗留问题，现已修复） |

## 7. 验证记录（2026-09-23）

- `make` EXIT=0，无告警；两个 util 工具 `gfortran` 编译 EXIT=0。
- **旧算例 A/B**（`cases/solid_1d`、`cases/channel_ac`，旧二进制 vs 新二进制）：
  `output_para.out` 除新增 4 行外**逐字节一致**。
- **新格式**：3 组写法取值正确；`$freestream_ec` 覆盖旧组生效（Ma 0.1→0.3）。
- **严格模式**：组内变量名写错 ⇒ 报错指出组名并以非零状态退出。
- **模板整体**：7 组同时存在 ⇒ 全部解析成功；`Solid_Max_Iter=30` 生效
  （日志 `GS iterations: 31` = 30+1）；算例跑完并正常输出。
- **后处理工具**：新旧两种 `control.ec` 都能正确取到 `Ma/Re/T_inf` 并跑完（EXIT=0）。

## 8. 重启后的耦合状态处理（2026-09-25 新增，`Iflag_Couple_Restart`）

**背景**：交错耦合（`Iflag_Couple_Scheme=1`）跑完 `Niter_Couple_Outer` 轮就正常退出（`t_end` 在该模式下无效）。
续算时若把整段调度（暖机 + 每轮步数减半）从头再播，并且界面量从 `control.ec` 初值重播，
就是"白跑几轮 + 界面跳变"。现在把**耦合状态**写进重启文件，由 `Iflag_Couple_Restart` 决定怎么续。

- `field_restart.dat` 版本 `iver` 1 → **2**：所有块之后追加一条 trailer
  （`mode/conv/pair/iter/nf,nk` + `max|dT_w|/tol` + 界面量 `T_w, q_w, u, p_w`）。
  **旧文件（iver=1）仍可读**，耦合状态视为"未知" ⇒ 保持旧行为。
- `Iflag_Couple_Restart`（`$couple_ec`，默认 **0**）：

  | 取值 | 含义 |
  |---|---|
  | `0` | **自动**：上次结束时已满足 `Tol_Couple_Tw` ⇒ 不再播交错，直接切到**逐步强耦合**（`Iflag_Couple_Scheme` 内部置 0、`t_end` 重新生效）一直跑到 `t_end`；否则继续交错，但**界面量从重启恢复**（跨重启零跳变） |
  | `1` | 重启后强制逐步强耦合 |
  | `2` | 重启后强制走交错驱动（旧行为） |
  | `-1` | 完全忽略重启里的耦合状态（不判定、不恢复）＝ 2026-09-25 之前的行为（回归测试用） |

- 切到逐步强耦合时会**自动限制固体每步 GS 预算**：`Solid_Max_Iter=min(现值,20)`、
  `Solid_Tol=max(现值,1e-6)`（实测默认 20000/1e-9 时**每步**固体求解约 15–20 s，无法长跑）；
  松弛因子与物理不变，实际取值会打印并回显到 `output_para.out`。
- 该路径退出前会**补写一次重启文件**（否则 `t_end` 落在两次 `Kstep_save` 之间就没有文件可续）。
- case1 这类"`max|dT_w|` 平台 ~0.1 K > `Tol_Couple_Tw=1e-2`"的算例永远不会自动判定为已收敛：
  要么放宽 `Tol_Couple_Tw`，要么直接 `Iflag_Couple_Restart=1`。

