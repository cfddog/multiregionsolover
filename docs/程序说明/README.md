# OpenCFD-EC 1.16a 程序说明（LaTeX 源）

本目录是**程序说明文档**的 LaTeX 源码，最终 PDF 为
`docs/OpenCFD-EC-1.16a-程序说明.pdf`（根目录另有一份副本）。

## 编译

```bash
cd docs/程序说明
./build.sh          # 生成 main.pdf，并拷贝到 docs/OpenCFD-EC-1.16a-程序说明.pdf
./build.sh clean    # 清理中间文件
```

依赖：`texlive-xetex`、`texlive-lang-chinese`（`ctexart`）、Fandol 字体
（`texlive-fonts-extra`，本机路径 `/usr/share/texlive/texmf-dist/fonts/opentype/public/fandol`）。
**必须用 `xelatex`**（或 `lualatex`），`pdflatex` 无法处理中文。

## 目录结构

| 文件 | 内容 |
|---|---|
| `main.tex` | 主文件（导言区、标题、目录、`\input` 各章） |
| `chap01_overview.tex` | 概述：目标、块类型、接口码、面编号、单位与无量纲化、目录结构 |
| `chap02_inputs.tex` | 输入文件与算例结构（control.ec / 网格 / bc3d / material / porous / solid_bc） |
| `chap03_compressible.tex` | 可压高速求解器与参数 |
| `chap04_lowspeed.tex` | 低速不可压求解器（SIMPLE/SIMPLEC/AC）与参数 |
| `chap05_porous.tex` | 多孔介质（DBF + LTNE）与参数 |
| `chap06_solid.tex` | 固体导热求解器与参数 |
| `chap07_coupling.tex` | 跨区域耦合：接口码、同步耦合、分段交错、重启耦合状态 |
| `chap08_output_restart.tex` | 输出文件与重启（含 `iver=2` trailer） |
| `chap09_case_fluid_solid.tex` | 三个 fluid_solid 耦合算例的完整设置 |
| `chap10_build_run.tex` | 构建、运行、并行与集群注意事项 |
| `chap11_validation.tex` | 验证/回归方法与已知问题 |
| `91_appendix_params.tex` | 附录 A：`control.ec` 全部参数总表（默认值） |
| `92_appendix_codes.tex` | 附录 B：面编号 / 接口码 / 物理边界码 / 文件格式速查 |
| `93_changelog.tex` | 附录 C：更新记录（**新功能必须追加在此**） |

## 维护约定（重要）

程序每新增/修改一个功能，必须：

1. 更新对应章节（物理模型、控制参数表、示例）；
2. 若新增参数：同步附录 A 参数表；
3. 在 `93_changelog.tex` 追加一条记录（日期 + 功能 + 涉及源码 + 验证算例）；
4. 重新运行 `./build.sh` 并确认编译无错误（见 `build.log`）。
