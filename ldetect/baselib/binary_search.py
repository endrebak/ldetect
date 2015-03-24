import bisect

# From docs.python.org

def index(a, x):
    'Locate the leftmost value exactly equal to x'
    i = bisect.bisect_left(a, x)
    if i != len(a) and a[i] == x:
        return i
    raise ValueError('index: Not found!')

def find_lt_ind(a, x):
    'Find rightmost value less than x'
    i = bisect.bisect_left(a, x)
    if i:
        return i-1
    raise ValueError('find_lt: Not found!')

def find_lt(a, x):
    i = find_lt_ind(a, x)
    return a[i]

def find_le_ind(a, x):
    'Find rightmost value less than or equal to x'
    i = bisect.bisect_right(a, x)
    if i:
        return i-1
    raise ValueError('find_le: Not found!')

def find_le(a, x):
    i = find_le_ind(a, x)
    return a[i]

def find_gt_ind(a, x):
    'Find leftmost value greater than x'
    i = bisect.bisect_right(a, x)
    if i != len(a):
        return i
    raise ValueError('find_gt: Not found!')

def find_gt(a, x):
    i = find_gt_ind(a, x)
    return a[i]

def find_ge_ind(a, x):
    'Find leftmost item greater than or equal to x'
    i = bisect.bisect_left(a, x)
    if i != len(a):
        return i
    print('a: '+repr(a))
    print('x: '+repr(x))
    raise ValueError('find_ge: Not found!')

def find_ge(a, x):
    i = find_ge_ind(a, x)
    return a[i]
