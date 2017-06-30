import os
from collections import Counter, OrderedDict
from math import log, erfc, sqrt


gzopen = lambda f: gzip.open(f, mode="rt") if f.endswith(".gz") else open(f)


def touch(fname, mode=0o666, dir_fd=None, **kwargs):
    """https://stackoverflow.com/questions/1158076/implement-touch-using-python"""
    flags = os.O_CREAT | os.O_APPEND
    with os.fdopen(os.open(fname, flags=flags, mode=mode, dir_fd=dir_fd)) as f:
        os.utime(f.fileno() if os.utime in os.supports_fd else fname,
            dir_fd=None if os.supports_fd else dir_fd, **kwargs)


def index_of_list_items(lists):
    """Find the highest index for items among list of lists and return ordered dictionary sorted by
    decreasing index value.

    Args:
        lists (list): list of lists

    Returns:
        OrderedDict of list items, sorted by depth with deepest list item first

    Example:
        >>> lineages = [['1', '131567', '99999999', '2'],
                        ['1', '131567', '99999999', '2', '1224', '1236', '72274', '135621', '286'],
                        ['1', '131567', '99999999', '2', '1224'],
                        ['1', '131567', '99999999', '2', '1224', '1236']]
        >>> index_of_list_items(lineages)
        OrderedDict([('286', 8),
                     ('135621', 7),
                     ('72274', 6),
                     ('1236', 5),
                     ('1224', 4),
                     ('2', 3),
                     ('99999999', 2),
                     ('131567', 1),
                     ('1', 0)])

    """
    indexes = {}
    for l in lists:
        for i, item in enumerate(l):
            if item in indexes:
                if indexes[item] < i:
                    indexes[item] = i
            else:
                indexes[item] = i
    return OrderedDict(sorted(indexes.items(), key=lambda t: t[1], reverse=True))


def nettleton_pvalue(items, key):
    """Calculate pvalue based on Nettleton result.

    Args:
        items (list): list of items
        key (string): key within items you're testing

    Returns:
        float

    Raises:
        AssertionError when key is not present in items
    """
    assert key in items
    if len(items) == 1:
        return 1

    item_counts = Counter(items)
    if item_counts.most_common()[0][0] == key:
        try:
            # take second entries count value
            max_count = item_counts.most_common()[1][1]
        except IndexError:
            max_count = 0
    else:
        # take first entries count value
        max_count = item_counts.most_common()[0][1]
    if item_counts[key] <= max_count:
        return 1
    else:
        try:
            t = 2 * (max_count * log((2 * max_count) / (max_count + item_counts[key])) \
                  + (item_counts[key] * log((2 * item_counts[key] / \
                                                         (max_count + item_counts[key])))))
        except ValueError:
            # max_count = 0
            t = 2 * (item_counts[key] * log((2 * item_counts[key] / \
                                                         (max_count + item_counts[key]))))
        return erfc(sqrt(t / 2))
