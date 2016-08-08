# Like the naive algorithm but we break out of the inner loop as soon as our
# mismatch budget exceeds the maximum allowed Hamming distance.
def naive_approx_hamming(p, t, maxHammingDistance=1):
    occurrences = []
    for i in range(0, len(t) - len(p) + 1): # for all alignments
        nmm = 0
        for j in range(0, len(p)):          # for all characters
            if t[i+j] != p[j]:               # does it match?
                nmm += 1                     # mismatch
                if nmm > maxHammingDistance:
                    break                    # exceeded maximum distance
        if nmm <= maxHammingDistance:
            # approximate match; return pair where first element is the
            # offset of the match and second is the Hamming distance
            occurrences.append((i, nmm))
    return occurrences
