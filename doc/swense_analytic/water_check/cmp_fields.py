import os, re, sys
def read_field(path):
    txt = open(path, errors='replace').read()
    m = re.search(r'internalField\s+nonuniform\s+List<(\w+)>\s*\n?(\d+)\s*\n?\(', txt)
    if not m:
        m2 = re.search(r'internalField\s+uniform\s+([^;]+);', txt)
        return ('uniform', m2.group(1).strip() if m2 else None)
    n = int(m.group(2)); start = m.end()
    body = txt[start:]
    # tokens between parentheses per entry (vectors) or scalars
    body = body[:body.index(')') if ')' in body else len(body)]
    if m.group(1) == 'vector':
        vals = [float(x) for x in re.findall(r'[-+0-9.eE]+', body.replace('(', ' ').replace(')', ' '))]
    else:
        vals = [float(x) for x in re.findall(r'[-+0-9.eE]+', body)]
    return ('nonuniform', vals)
a, b = sys.argv[1], sys.argv[2]
files = sorted(f for f in os.listdir(a) if os.path.isfile(os.path.join(a, f)) and f not in ('uniform',))
worst = 0.0
for f in files:
    pb = os.path.join(b, f)
    if not os.path.exists(pb):
        print('missing in B:', f); continue
    ta, va = read_field(os.path.join(a, f))
    tb, vb = read_field(pb)
    if ta != 'nonuniform' or tb != 'nonuniform':
        print(f, 'skipped (non-field):', ta, tb); continue
    if len(va) != len(vb):
        print(f, 'length differs', len(va), len(vb)); worst = float('inf'); continue
    d = max((abs(x - y) for x, y in zip(va, vb)), default=0.0)
    worst = max(worst, d)
    print(f'{f:24s} n={len(va):7d} max abs diff = {d:.3e}')
print('WORST', worst)
