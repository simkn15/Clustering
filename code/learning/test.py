def r(tp, fn):
    return tp / (tp + fn)

def p(tp, fp):
    return tp / (tp + fp)

def f(r, p):
    return ((2 * (r * p)) / (p + r))

r1 = r(7., 1)
p1 = p(7., 2)
r2 = r(15., 2)
p2 = p(15, 3)

print f(r1, p1)
print f(r2, p2)
