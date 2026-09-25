#!/bin/bash
# 编译 OpenCFD-EC 程序说明（LaTeX/中文）。
# 需要 texlive-xetex + texlive-lang-chinese（ctex）+ Fandol 字体（texlive-fonts-extra）。
# 用法：  ./build.sh            # 生成 main.pdf
#         ./build.sh clean      # 清理中间文件
set -e
cd "$(dirname "$0")"
if [ "$1" = "clean" ]; then
  rm -f main.aux main.log main.out main.toc main.pdf
  echo "cleaned."
  exit 0
fi
if ! command -v xelatex >/dev/null 2>&1; then
  echo "ERROR: xelatex not found (install texlive-xetex)." >&2; exit 1
fi
# 两次编译以生成目录/书签；第三次可省略但更保险
xelatex -interaction=nonstopmode -halt-on-error main.tex > build.log 2>&1 || {
  echo "xelatex FAILED - see build.log (tail):"; tail -30 build.log; exit 1; }
xelatex -interaction=nonstopmode -halt-on-error main.tex >> build.log 2>&1 || {
  echo "xelatex (2nd pass) FAILED - see build.log (tail):"; tail -30 build.log; exit 1; }
cp -f main.pdf "../OpenCFD-EC-1.16a-程序说明.pdf"
echo "OK: docs/程序说明/main.pdf"
echo "    docs/OpenCFD-EC-1.16a-程序说明.pdf"
