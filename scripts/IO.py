

def cat_reads(f1, f2, outF):
    with open(outF, 'w') as catted:
        for fpath in [f1, f2]:
            with open(fpath, 'r') as infile:
                for line in infile:
                    catted.write(line)


def interleave_reads(f1, f2, outF):
    with open(outF, 'w') as interleaved:
        with open(f1, 'r') as infile1:
            with open(f2, 'r') as infile2:
                for l1, l2 in zip(infile1, infile2):
                    interleaved.write(l1)
                    interleaved.write(l2)
