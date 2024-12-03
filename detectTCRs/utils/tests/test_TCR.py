from detectTCRs.utils.TCR import TCR

V = 'TRBV7-8*01'
J = 'TRBJ2-3*01'
CDR3 = 'CASSLNRLGLDTQYF'

test_tcr = TCR(V, J, CDR3)
assert test_tcr.v == V
assert test_tcr.j == J
assert test_tcr.cdr3 == CDR3

assert test_tcr.full(
) == 'GAGVSQSPRYKVAKRGQDVALRCDPISGHVSLFWYQQALGQGPEFLTYFQNEAQLDKSGLPSDRFFAERPEGSVSTLKIQRTQQEDSAVYLCASSLNRLGLDTQYFGPGTRLTVL'
assert test_tcr.aligned() == 'CASSLNRLGLDTQY'
