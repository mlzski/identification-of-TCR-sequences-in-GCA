# TCR class
from detectTCRs.utils.cdr3_to_full_seq import to_full_seq
from detectTCRs.utils.align import Align_sequences


class TCR:

    def __init__(self, V, J, CDR3):
        self.v = V
        self.j = J
        self.cdr3 = CDR3

    def full(self):
        return str(to_full_seq(self.v, self.j, self.cdr3)[0])

    def aligned(self):
        sequence = [("human0:TRB", self.full())]
        return Align_sequences(sequence)[0]
