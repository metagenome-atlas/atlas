import os



def gunzip(path, dst):
    subprocess.call(['gunzip', '<', path, '>', dst])


if __name__ == '__main__':
    gunzip('./test.fastq.gz', '.')
