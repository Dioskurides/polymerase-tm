import pytest
from polymerase_tm import tm, primer_hairpin, primer_dimer

def test_primer3_tm():
    seq = "ATGTCCCTGCTCTTCTCTCGATGCAA"
    # Internal owczarzy
    tm_int = tm(seq, method="owczarzy")
    # Primer3
    tm_p3 = tm(seq, method="primer3")
    
    print(f"Internal Tm: {tm_int}")
    print(f"Primer3 Tm: {tm_p3}")
    
    # They should be roughly similar (within a few degrees)
    assert abs(tm_int - tm_p3) < 5

def test_primer3_hairpin():
    seq = "AAAAAAAATTTTTTTT" # perfect hairpin
    hps = primer_hairpin(seq)
    assert len(hps) > 0
    assert "dg" in hps[0]
    assert hps[0]["dg"] < 0 # strong hairpin

def test_primer3_dimer():
    fwd = "ATGCATGCATGC"
    rev = "GCATGCATGCAT"
    res = primer_dimer(fwd, rev)
    assert "hetero_dimer_dg" in res
    assert "risk_level" in res
