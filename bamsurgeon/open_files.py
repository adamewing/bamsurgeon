import pysam

compression_threads = 1

def open_aligment_file(*args):
    return pysam.AlignmentFile(*args, threads=compression_threads)
