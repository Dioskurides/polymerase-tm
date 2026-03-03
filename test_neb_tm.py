import math
from polymerase_tm.constants import NN_PARAMS

def calc_neb_tm(seq: str, use_sch=True, use_sl=False) -> float:
    dH = 200.0  # cal/mol
    dS = -5.7   # cal/K*mol
    
    # AT penalty
    if seq[0] in "AT":
        dH += 2200.0
        dS += 6.9
    if seq[-1] in "AT":
        dH += 2200.0
        dS += 6.9
        
    for i in range(len(seq) - 1):
        pair = seq[i:i+2]
        h, s = NN_PARAMS[pair]
        dH += h * 1000.0  # kcal to cal
        dS += s
        
    R = 1.987
    Ct = 500e-9
    salt_M = 0.150
    N = len(seq)
    
    # SantaLucia Entropy salt correction
    if use_sl:
        dS += 0.368 * (N - 1) * math.log(salt_M)
        
    tm_raw = (dH) / (dS + R * math.log(Ct / 4)) - 273.15
    
    # Schildkraut Tm salt correction
    if use_sch:
        # 16.6 * log10(saltM)
        tm_raw += 16.6 * math.log10(salt_M)
        
    return tm_raw

rev = "CACCATATGCGGTGTG"
fwd = "ATGGCGATGATTACGAATTC"  # Right side 20nt
fwd17 = "GCGATGATTACGAATTC"     # Right side 17nt

for name, seq in [("REV", rev), ("FWD 20nt", fwd), ("FWD 17nt", fwd17)]:
    print(f"=== {name} ({len(seq)}nt) ===")
    
    tm1 = calc_neb_tm(seq, use_sch=True, use_sl=False)
    print(f"Schildkraut (Tm+16.6logNa) : {tm1:.4f} -> {math.round(tm1) if hasattr(math, 'round') else math.floor(tm1 + 0.5)}")
    
    tm2 = calc_neb_tm(seq, use_sch=False, use_sl=True)
    print(f"SantaLucia (dS+0.368lnNa)  : {tm2:.4f} -> {math.floor(tm2 + 0.5)}")
