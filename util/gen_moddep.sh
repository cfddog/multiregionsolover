#!/bin/bash
# ---------------------------------------------------------------------------
# 生成 src/modules.dep：Fortran `module` / `use` 关系 -> 文件级 make 依赖。
#
# 仓库根 makefile 通过 `include src/modules.dep` 使用它，因此：
#   * `make -j` 可以并行编译（首次构建也安全）；
#   * 不产生 .d 中间文件，build/ 里只有 *.o / *.mod。
#
# 触发时机：新增/删除源文件，或改动了 `use`/`module` 关系之后，运行
#   bash util/gen_moddep.sh
# 再把 src/modules.dep 一起提交。
# ---------------------------------------------------------------------------
set -e

R="$(cd "$(dirname "$0")/.." && pwd)"
cd "$R"

LC_ALL=C awk '
FNR==1 { file=FILENAME; sub(/^src\//,"",file); sub(/\.f90$/,"",file) }
{
  line=tolower($0); sub(/!.*/,"",line)
  if (line ~ /^[ \t]*module[ \t]+[a-z0-9_]+/ &&
      line !~ /module[ \t]+(procedure|subroutine|function)/) {
    split(line,a,/module[ \t]+/); m=a[2]; sub(/[^a-z0-9_].*/,"",m)
    if (m!="") def[m]=file
  }
  if (line ~ /^[ \t]*use[ \t]+[a-z0-9_]+/) {
    split(line,a,/use[ \t]+/); u=a[2]; sub(/[^a-z0-9_].*/,"",u)
    if (u!="" && u!="iso_fortran_env" && u!="iso_c_binding" && u!="omp_lib") used[file" "u]=1
  }
}
END {
  for (k in used) {
    split(k,kk," "); f=kk[1]; u=kk[2]
    if (u in def) { d=def[u]; if (d!=f) dep[f]=dep[f]" "d }
  }
  for (f in dep) {
    printf("$(BUILDDIR)/%s.o:", f)
    split(dep[f],arr," ")
    for (i=1;i<=length(arr);i++) if (arr[i]!="" && !(arr[i] in seen)) { seen[arr[i]]=1; printf(" $(BUILDDIR)/%s.o", arr[i]) }
    printf("\n")
    delete seen
  }
}' src/*.f90 | sort > src/modules.dep

echo "wrote src/modules.dep ($(wc -l < src/modules.dep) rules)"
