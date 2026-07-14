#!/usr/bin/env python3
"""
ad_fix_varargs.py -- wrap printf-family arguments in spval() at sites
flagged by -Werror=conditionally-supported (sfloat passed through
varargs prints garbage otherwise).

Usage: ad_fix_varargs.py <build.log> <srcdir>

Reads 'file:line:col: error: passing objects ... through ...' from the
log, finds the printf/fprintf/sprintf/snprintf call at each site, and
wraps every argument after the format with spval(...). spval is an
identity template for non-sfloat types, so blanket wrapping is safe.
Idempotent: args already wrapped are skipped.
"""

import re, sys, os

FUNCS = {"printf": 1, "fprintf": 2, "sprintf": 2, "snprintf": 3}

def split_args(s):
    """split argument list at top-level commas, respecting (), "" and ''"""
    args, depth, cur, i = [], 0, "", 0
    instr = None
    while i < len(s):
        c = s[i]
        if instr:
            cur += c
            if c == "\\": cur += s[i+1]; i += 1
            elif c == instr: instr = None
        elif c in "\"'": instr = c; cur += c
        elif c in "([{": depth += 1; cur += c
        elif c in ")]}": depth -= 1; cur += c
        elif c == "," and depth == 0: args.append(cur); cur = ""
        else: cur += c
        i += 1
    args.append(cur)
    return args

def find_call(text, lo, hi):
    """innermost printf-family call whose span overlaps [lo,hi];
       returns (name, open_idx, close_idx) or None"""
    best = None
    for m in re.finditer(r"\b(v?f?s?n?printf)\s*\(", text):
        name = m.group(1).lstrip("v")
        if name not in FUNCS: continue
        # find matching close paren
        depth, j, instr = 0, m.end()-1, None
        while j < len(text):
            c = text[j]
            if instr:
                if c == "\\": j += 1
                elif c == instr: instr = None
            elif c in "\"'": instr = c
            elif c == "(": depth += 1
            elif c == ")":
                depth -= 1
                if depth == 0: break
            j += 1
        if m.start() <= hi and j >= lo:       # call span overlaps line
            if best is None or m.start() > best[1]:
                best = (name, m.end()-1, j)   # innermost enclosing call
    return best

def fix_file(path, lines_flagged):
    src = open(path).read()
    offsets = [0]
    for line in src.split("\n"): offsets.append(offsets[-1]+len(line)+1)
    edits = []
    done_spans = set()
    for ln in sorted(lines_flagged):
        lo, hi = offsets[ln-1], offsets[ln]-1
        call = find_call(src, lo, hi)
        if call is None:
            print(f"  !! no call found {path}:{ln}"); continue
        name, op, cl = call
        if (op,cl) in done_spans: continue
        done_spans.add((op,cl))
        inner = src[op+1:cl]
        args = split_args(inner)
        nfmt = FUNCS[name]
        changed = False
        for k in range(nfmt, len(args)):
            a = args[k].strip()
            if not a or a.startswith("spval(") or a.startswith('"'): continue
            lead = args[k][:len(args[k])-len(args[k].lstrip())]
            args[k] = f"{lead}spval({a})"
            changed = True
        if changed:
            edits.append((op+1, cl, ",".join(args)))
    # apply edits back-to-front
    for op, cl, rep in sorted(edits, reverse=True):
        src = src[:op] + rep + src[cl:]
    if edits:
        open(path, "w").write(src)
    return len(edits)

def main():
    log, srcdir = sys.argv[1], sys.argv[2]
    flagged = {}
    pat = re.compile(r"\.\./([\w./]+):(\d+):\d+: error: passing objects")
    for line in open(log):
        m = pat.search(line)
        if m:
            flagged.setdefault(m.group(1), set()).add(int(m.group(2)))
    total = 0
    for f, lines in sorted(flagged.items()):
        n = fix_file(os.path.join(srcdir, f), lines)
        print(f"{f}: {n} call(s) wrapped")
        total += n
    print(f"total {total} calls fixed in {len(flagged)} files")

if __name__ == "__main__":
    main()
