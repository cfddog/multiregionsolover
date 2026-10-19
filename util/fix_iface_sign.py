#!/usr/bin/env python3
"""Check / repair the (i,j,k) sign convention of bc3d*.inp interface lines.

Convention (documented in src/sub_convert_inp.f90::Convert_bc)
-------------------------------------------------------------
Each face line holds three node ranges

        (ib,ie) (jb,je) (kb,ke)  type

For a connection face (type < 0) -- both on the local line and on the
neighbour (continuation) line -- the three pairs must be

  * exactly ONE collapsed pair       st == ed       -> the face NORMAL
  * exactly ONE fully negative pair  st<0 and ed<0  -> the coupling edge
  * exactly ONE fully positive pair  st>0 and ed>0  -> other tangential edge

The i/j/k correspondence is then unique: negative<->negative,
positive<->positive, normal<->normal, and the order inside a direction
follows |st| -> |ed| (start to start).

Two positive tangential ranges make the correspondence ambiguous (the old
reader then produced e.g. a "k-k" link on a 2-D case).

bc3d.inp never carries the connection: for same-class pairs it either uses
type 11..19 (e.g. 14 = solid-solid) with NO pairing line, or the generic -1
with the pairing line.  A type-only face must therefore be written with plain
ascending positive ranges -- the sign is reserved for the pairing descriptor
of bc3d_interface.inp (single source of truth for block connections).

Use
---
    python3 util/fix_iface_sign.py --check cases/*/bc3d*.inp
    python3 util/fix_iface_sign.py --fix   cases/*/bc3d*.inp

--fix marks, on both sides of every offending pair, the last tangential
direction (k, else j, else i) negative -- or keeps an already marked one --
so that the pairing obeys the convention.  Physical faces (type > 0) are
never touched.
"""
import glob
import sys

INTERFACE_CODES = range(11, 20)      # cross-class codes written by gen_bc3d


def read_records(path):
    """-> [[block, face_line_index, face_fields, pair_fields_or_None]]

    A continuation (pairing) line follows a face when
      * bc3d.inp           : bc < 0
      * bc3d_interface.inp : bc < 0 or 11 <= bc <= 19   (is_interface_bc)
    """
    txt = open(path, errors='ignore').read().splitlines()
    iface_file = path.endswith('bc3d_interface.inp')
    recs = []
    blk = 0
    i = 0
    while i < len(txt):
        v = txt[i].split()
        if len(v) == 3:
            try:
                [int(x) for x in v]
                blk += 1
            except ValueError:
                pass
        if len(v) == 7:
            try:
                f = [int(x) for x in v]
            except ValueError:
                i += 1
                continue
            pair = None
            idx = i
            cont = f[6] < 0 or (iface_file and f[6] in INTERFACE_CODES)
            if cont:
                p = txt[i + 1].split()
                if len(p) == 7:
                    pair = [int(x) for x in p]
                    i += 1
            recs.append([blk, idx, f, pair])
        i += 1
    return recs


def cls(f):
    """-> (collapsed_dim, [negative dims], [positive dims], [mixed dims])"""
    coll = [k for k in range(3) if f[2 * k] == f[2 * k + 1]]
    neg = [k for k in range(3) if f[2 * k] < 0 and f[2 * k + 1] < 0]
    pos = [k for k in range(3) if f[2 * k] > 0 and f[2 * k + 1] > 0]
    mixed = [k for k in range(3) if (f[2 * k] < 0) != (f[2 * k + 1] < 0)]
    if len(coll) == 1:
        for lst in (neg, pos):
            if coll[0] in lst:
                lst.remove(coll[0])
    return (coll[0] if len(coll) == 1 else None, neg, pos, mixed)


def extent(f, d):
    return abs(f[2 * d + 1] - f[2 * d])


def compliant(f):
    c, neg, pos, mixed = cls(f)
    return c is not None and len(neg) == 1 and len(pos) == 1 and not mixed


def partner_map(f, p):
    """local tangential dim -> neighbour tangential dim (or None)

    Used when at least one side is not compliant: match by node extent.
    """
    c = cls(f)[0]
    c1 = cls(p)[0]
    if c is None or c1 is None:
        return None
    T = [d for d in range(3) if d != c]
    T1 = [d for d in range(3) if d != c1]
    out = {}
    for d in T:
        m = [e for e in T1 if extent(p, e) == extent(f, d)]
        if len(m) != 1:
            return None
        out[d] = m[0]
    if len(set(out.values())) != 2:
        return None
    return out


def sign_map(f, p):
    """class-based correspondence: normal<->normal, neg<->neg, pos<->pos

    Returns a dict (dim -> dim) when both sides are compliant, else None.
    This is the rule Convert_bc applies, so two compliant faces are already
    unambiguous even when both tangential extents are equal.
    """
    if not (compliant(f) and compliant(p)):
        return None
    cf, nf, pf, _ = cls(f)
    cp, npp, pp, _ = cls(p)
    m = {cf: cp, nf[0]: npp[0], pf[0]: pp[0]}
    return m if len(set(m.values())) == 3 else None


def extents_agree(f, p, m):
    return all(extent(f, d) == extent(p, m[d]) for d in m)

def repaint(f, d):
    """make tangential dim d negative, every other tangential positive"""
    c = cls(f)[0]
    for e in range(3):
        if e == d:
            f[2 * e], f[2 * e + 1] = -abs(f[2 * e]), -abs(f[2 * e + 1])
        elif e != c:
            f[2 * e], f[2 * e + 1] = abs(f[2 * e]), abs(f[2 * e + 1])


def fix_pair(f, p):
    """fix one connection pair in place -> (changed, message)"""
    c, negf, posf, mixf = cls(f)
    c1, negp, posp, mixp = cls(p)
    if c is None or c1 is None or mixf or mixp:
        return False, 'unusable face (no unique collapsed pair / mixed signs)'
    mp = partner_map(f, p)
    if mp is None:
        return False, 'cannot match the tangential directions by extent'
    if compliant(f) and compliant(p):
        if mp[negf[0]] == negp[0]:
            return False, ''                     # already consistent
        dm, dm1 = negf[0], negp[0]               # marks disagree -> re-mark
    else:
        if len(negf) == 1:
            dm = negf[0]
        elif len(negp) == 1:
            cand = [d for d, e in mp.items() if e == negp[0]]
            dm = cand[0] if len(cand) == 1 else max(mp.keys())
        else:
            dm = max(mp.keys())                  # k, else j, else i
        dm1 = mp[dm]
    repaint(f, dm)
    repaint(p, dm1)
    return True, ''


def process(path, apply_fix):
    recs = read_records(path)
    changed = bad = 0
    msgs = []
    for rec in recs:
        blk, idx, f, p = rec
        if p is None:
            # bc3d.inp type-only interface face (11..19, e.g. 14 = solid-solid):
            # it carries NO pairing line, so it must not carry directional
            # signs either -- write plain ascending positive ranges.
            if f[6] in INTERFACE_CODES:
                want = []
                for d in range(3):
                    a, b = f[2 * d], f[2 * d + 1]
                    want += [min(abs(a), abs(b)), max(abs(a), abs(b))]
                if want != f[:6]:
                    if not apply_fix:
                        bad += 1
                        msgs.append('  blk%d: type-only face %s coded %d carries signs'
                                    % (blk, f[:6], f[6]))
                    else:
                        f[:6] = want
                        changed += 1
            continue
        if compliant(f) and compliant(p):
            m = sign_map(f, p)
            if m is not None and extents_agree(f, p, m):
                continue                          # consistent, nothing to do
            if not apply_fix:
                bad += 1
                msgs.append('  blk%d: compliant but inconsistent pair %s / %s'
                            % (blk, f[:6], p[:6]))
                continue
        if not apply_fix:
            bad += 1
            msgs.append('  blk%d: local %s / nb %s' % (blk, f[:6], p[:6]))
            continue
        ok, msg = fix_pair(f, p)
        if ok:
            changed += 1
        else:
            bad += 1
            msgs.append('  blk%d: %s' % (blk, msg))
    if apply_fix and changed:
        rewrite(path, recs)
    return changed, bad, msgs


def format_like(orig, fields):
    """re-emit `fields` keeping the right edge of every original column"""
    toks = orig.split()
    if len(toks) != len(fields):
        return ''.join('%8d' % x for x in fields)
    ends = []
    j = 0
    for t in toks:
        k = orig.index(t, j)
        j = k + len(t)
        ends.append(j)
    widths = [ends[0]]
    for k in range(1, len(ends)):
        widths.append(ends[k] - ends[k - 1])
    widths = [max(w, len(str(x))) for w, x in zip(widths, fields)]
    return ''.join(str(x).rjust(w) for x, w in zip(fields, widths))


def rewrite(path, recs):
    """replace only the lines whose content changed (format preserved)"""
    lines = open(path, errors='ignore').read().splitlines()
    new_lines = {}
    for blk, idx, f, p in recs:
        orig = lines[idx].split()
        if [int(x) for x in orig] != f:
            new_lines[idx] = format_like(lines[idx], f)
        if p is not None:
            origp = lines[idx + 1].split()
            if [int(x) for x in origp] != p:
                new_lines[idx + 1] = format_like(lines[idx + 1], p)
    if not new_lines:
        return 0
    out = [new_lines.get(i, ln) for i, ln in enumerate(lines)]
    open(path, 'w').write('\n'.join(out) + '\n')
    return len(new_lines)


def main(argv):
    apply_fix = '--fix' in argv
    paths = []
    for a in argv:
        if a.startswith('-'):
            continue
        paths += sorted(glob.glob(a))
    total_changed = total_bad = 0
    for path in paths:
        changed, bad, msgs = process(path, apply_fix)
        total_changed += changed
        total_bad += bad
        flag = 'FIXED' if (apply_fix and changed) else ('BAD  ' if bad else 'OK   ')
        print('%s %-52s changed=%d bad=%d' % (flag, path, changed, bad))
        for m in msgs[:6]:
            print(m)
    print('total changed=%d, still bad=%d' % (total_changed, total_bad))
    return 0 if total_bad == 0 else 1


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))

