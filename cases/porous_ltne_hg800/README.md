# porous_ltne_hg800 — Step-2（解耦两步：高温高速流对流供热的多孔壁 LTNE）

沿用 `cases/porous_ltne_1d` 的全部多孔/冷却剂设定（L=10 mm、ε=0.3、k_s=13.4、
dp=1 mm、hv=2e6 W/m³K、冷端 G=2 kg/(m²·s) 300 K 空气注入 + 骨架 Robin hc=31.4、
骨架材质不变），仅把**热端定热流 q=2e6 换成由高速流体提供的对流热边界**：
`q = hg·(T_rec − Ts(L))`，其中 `hg=1.48e3 W/(m²·K)`、`T_rec=2240 K`
（由 `cases/high_low_fluid_800K` Step-1 标定：Ma=3、T∞=800 K）。
`solid_bc.inp` 面 4 用对流(Robin)型：`4  0.0 0.0 1.48e3 2240.0`。

## 结果（nx=401，150000 SIMPLE 步收敛，res_p~7e-17）
| 站位 | Tf 数值 | Tf 解析 | Ts 数值 | Ts 解析 | Ts−Tf |
|---|---|---|---|---|---|
| y=0    | 302.0 | 302.1 | 465.7 | 466.7 | 163.7 |
| y=L/2  | 605.1 | 605.9 | 662.0 | 662.8 | 56.9  |
| y=L    | 1059.3| 1059.4| 1198.4| 1198.7| 139.1 |

- Tf RMS 0.84 K、Ts RMS 0.76 K（相对量程 ≈0.1 %），与（对流热端）解析解吻合。
- 热端壁面 Ts(L)≈1198 K、冷却剂出口 Tf(L)≈1059 K。
- 热量自洽：`q_hot=hg(T_rec−Ts(L))≈1.54e6 W/m²` vs `G·cp·ΔTf+hc·(Ts(0)−300)≈1.53e6`
  （相对缺口 0.46%）。
- 与 Step-1 的 2.14e6 W/m² 差异属**解耦误差**：Step-1 的壁面替身温度≈797 K，而真正
  多孔壁热面为 1198 K，驱动温差变小。全收敛需"可压↔多孔（码 19）"共轭迭代。

复现：`python3 gen_case.py 401 2 2 && mpirun -np 1 ./opencfd-ec1.16a.out && conda run -n num_python python3 validate_conv.py`。
图：`validation_ltne_hg800.png`。
