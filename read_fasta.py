"""
Read a FASTA/FASTQ file.
See https://github.com/lh3/readfq
"""

def read_fasta(fin):
    "Read a FASTA/FASTQ file."
    last = None # this is a buffer keeping the last unprocessed line
    while True: # mimic closure; is it a bad idea?
        if not last: # the first record or a record following a fastq
            for line in fin: # search for the start of the next record
                if line[0] in '>@': # fasta/q header line
                    last = line[:-1] # save this line
                    break
        if not last:
            break
        xs = last[1:].split(None, 1)
        if len(xs) == 1:
            name = xs[0]
            bx = None
        else:
            name, bx = xs
            bx = bx[5:] if bx.startswith("BX:Z:") else None
        seqs = []
        last = None
        for line in fin: # read the sequence
            if line[0] in '@+>':
                last = line[:-1]
                break
            seqs.append(line[:-1])
        if not last or last[0] != '+': # this is a fasta record
            yield name, ''.join(seqs), bx, None # yield a fasta record
            if not last:
                break
        else: # this is a fastq record
            seq, leng, seqs = ''.join(seqs), 0, []
            for line in fin: # read the quality
                seqs.append(line[:-1])
                leng += len(line) - 1
                if leng >= len(seq): # have read enough quality
                    last = None
                    yield name, seq, bx, ''.join(seqs) # yield a fastq record
                    break
            if last: # reach EOF before reading enough quality
                yield name, seq, bx, None # yield a fasta record instead
                break
