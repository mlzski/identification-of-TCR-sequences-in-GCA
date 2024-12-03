from anarci import anarci


def Align_sequences(sequences):

    try:
        # Hand the list of sequences to the anarci function. Number them with the IMGT scheme
        results = anarci(sequences, scheme="imgt", output=False)
    except Exception as e:
        print(e)
        return ("", 0)

    # Unpack the results
    numbering, alignment_details, hit_tables = results

    if numbering[0] is None:
        aligned_CDR3 = ''
        domain_size = 0
        return (aligned_CDR3, domain_size)

    assert len(numbering) == len(alignment_details) == len(hit_tables
                                                           ) == len(sequences)

    # Iterate over the sequences
    for i in range(len(sequences)):
        domain_size = len(numbering[i])
        # Iterate over the domains
        for j in range(len(numbering[i])):
            domain_numbering, start_index, end_index = numbering[i][j]
            aligned_seq = [tup[1] for tup in domain_numbering]
            # Extract CDR3 region (105->117 positions together with anchor position (104))
            aligned_CDR3 = aligned_seq[103:117]
            aligned_CDR3 = ''.join(aligned_CDR3)

    return (aligned_CDR3, domain_size)
