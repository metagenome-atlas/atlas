def gen_names_for_range(N, prefix="", start=1):
    """generates a range of IDS with leading zeros so sorting will be ok"""
    n_leading_zeros = len(str(N))
    format_int = prefix + "{:0" + str(n_leading_zeros) + "d}"
    return [format_int.format(i) for i in range(start, N + start)]
