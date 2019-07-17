#!/usr/bin/env pypy3
"""
Hash k-mers and integers
Written by Shaun Jackman @sjackman
Edits by Lauren Coombe @lcoombe
"""

def hash_int(key, mask=0xffffffffffffffff):
    """
    Hash a 64-bit integer (invertible).
    See https://gist.github.com/lh3/59882d6b96166dfc3d8d
    """
    assert 0 <= key < 0x10000000000000000
    key = (~key + (key << 21)) & mask # key = (key << 21) - key - 1
    key = key ^ key >> 24
    key = ((key + (key << 3)) + (key << 8)) & mask # key * 265
    key = key ^ key >> 14
    key = ((key + (key << 2)) + (key << 4)) & mask # key * 21
    key = key ^ key >> 28
    key = (key + (key << 31)) & mask
    assert 0 <= key < 0x10000000000000000
    return key

def unhash_int(key, mask=0xffffffffffffffff):
    """
    Invert hash_int.
    https://gist.githubusercontent.com/lh3/974ced188be2f90422cc/raw/55fbbb63e489328fd9d1641897954ca997b65951/inthash.c
    """
    assert 0 <= key < 0x10000000000000000

    # Invert key = key + (key << 31)
    tmp = (key - (key << 31))
    key = (key - (tmp << 31)) & mask

    # Invert key = key ^ (key >> 28)
    tmp = key ^ key >> 28
    key = key ^ tmp >> 28

    # Invert key *= 21
    key = (key * 14933078535860113213) & mask

    # Invert key = key ^ (key >> 14)
    tmp = key ^ key >> 14
    tmp = key ^ tmp >> 14
    tmp = key ^ tmp >> 14
    key = key ^ tmp >> 14

    # Invert key *= 265
    key = (key * 15244667743933553977) & mask

    # Invert key = key ^ (key >> 24)
    tmp = key ^ key >> 24
    key = key ^ tmp >> 24

    # Invert key = (~key) + (key << 21)
    tmp = ~key
    tmp = ~(key - (tmp << 21))
    tmp = ~(key - (tmp << 21))
    key = ~(key - (tmp << 21)) & mask

    assert 0 <= key < 0x10000000000000000
    return key

# Translate the ACGT to ASCII 0-3.
TRANSLATE_ACGT_0123 = str.maketrans("ACGT", "\0\1\2\3")

# Translate the ASCII 0-3 to ACGT.
TRANSLATE_0123_ACGT = str.maketrans("\0\1\2\3", "ACGT")

# Translate the ACGT to TGCA.
TRANSLATE_ACGT_TGCA = str.maketrans("ACGT", "TGCA")

# Complement nucleotides.
TRANSLATE_COMPLEMENT = str.maketrans(
    "ACGTUNMRWSYKVHDBacgtunmrwsykvhdb",
    "TGCAANKYWSRMBDHVtgcaankywsrmbdhv")

def reverse_complement(seq):
    "Return the reverse complement of this sequence."
    return seq[::-1].translate(TRANSLATE_COMPLEMENT)

def kmer_to_int(kmer):
    "Convert a k-mer to an integer."
    x = 0
    for c in kmer.translate(TRANSLATE_ACGT_0123):
        assert 0 <= ord(c) < 4
        x <<= 2
        x += ord(c)
    return x

def int_to_kmer(x, k):
    "Convert an integer to a k-mer."
    assert x >= 0
    xs = k * [None]
    for i in reversed(range(k)):
        xs[i] = chr(x & 3)
        x >>= 2
    assert x == 0
    return str.join("", xs).translate(TRANSLATE_0123_ACGT)

def hash_kmer(kmer):
    "Hash a k-mer to an integer."
    return min(hash_int(kmer_to_int(kmer[0])), hash_int(kmer_to_int(reverse_complement(kmer[0])))), kmer[1]


def unhash_kmer(x, k):
    "Unhash an integer to a k-mer."
    return int_to_kmer(unhash_int(x), k)
