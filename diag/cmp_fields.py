# max |on - off| of the internal fields of two cases at one time
# python3 cmp_fields.py caseOn caseOff time field...
import sys, numpy as np
def internal(f):
    s = open(f).read(); i = s.index('internalField'); line = s[i:s.index('\n', i)]
    if 'nonuniform' not in line:
        return np.fromstring(line.split('uniform')[1].replace('(', ' ').replace(')', ' ').replace(';', ' '), sep=' ')
    a = s.index('(', s.index('List', i)); b = s.index('\n)', a)
    return np.fromstring(s[a+1:b].replace('(', ' ').replace(')', ' '), sep=' ')
on, off, t = sys.argv[1:4]
for f in sys.argv[4:]:
    A = internal('%s/%s/%s' % (on, t, f)); B = internal('%s/%s/%s' % (off, t, f))
    d = np.abs(A - B).max(); print('%-12s t %s  max|on-off| %.3e  (max|off| %.3e)' % (f, t, d, np.abs(B).max()))
