"""
utils for math/stat operations
"""

def sum_arrays(a, b):
    '''
    calcs per-element sum of two arrays (lists)
    e.g., list a[1,2] + list b[3,4] --> [3,6]
    '''
    return [x + y for x, y in zip(a, b)]


def sum_array_list(arrays):
    '''
    calcs per-element sum for a 
    list of arrays (lists)
    e.g., [[1,2],[3,4]] --> [3,6]
    '''
    return [sum(x) for x in zip(*arrays)]


def average_array_list(arrays):
    '''
    calcs per-element average for
    a list of arrays (lists)
    e.g., [[1,2],[3,4]] --> [1.5,3]
    '''

    n = len(arrays)
    return [float(x) / float(n) for x in sum_array_list(arrays)]
