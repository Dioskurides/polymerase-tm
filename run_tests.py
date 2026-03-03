# -*- coding: utf-8 -*-
import math
from polymerase_tm.constants import NN_PARAMS
from polymerase_tm.core import calc_nn_raw, owczarzy_correction

# Check our NN initiation values
print("Our initiation (cal):")
print(f"  dH_init = {NN_PARAMS['init']['dH']}")
print(f"  dS_init = {NN_PARAMS['init']['dS']}")
print(f"  dH_AT = {NN_PARAMS['AT_penalty']['dH']}")
print(f"  dS_AT = {NN_PARAMS['AT_penalty']['dS']}")

print()
print("NEB JS values (from extraction):")
print("  dH_init = 200 cal/mol")
print("  dS_init = -5.7 cal/molK")
print("  AT_penalty dH = 2200 cal/mol per end")
print("  AT_penalty dS = 6.9 cal/molK per end")

# Raw NN calc step by step for REV = CACCATATGCGGTGTG
seq = "CACCATATGCGGTGTG"
print(f"\n=== Step-by-step for REV ({seq}) ===")

# Count AT terminals
has_5prime_AT = seq[0] in "AT"
has_3prime_AT = seq[-1] in "AT"
print(f"5' base: {seq[0]} -> AT terminal: {has_5prime_AT}")
print(f"3' base: {seq[-1]} -> AT terminal: {has_3prime_AT}")

# Compute NN sum
dH = NN_PARAMS["init"]["dH"]
dS = NN_PARAMS["init"]["dS"]
if has_5prime_AT:
    dH += NN_PARAMS["AT_penalty"]["dH"]
    dS += NN_PARAMS["AT_penalty"]["dS"]
if has_3prime_AT:
    dH += NN_PARAMS["AT_penalty"]["dH"]
    dS += NN_PARAMS["AT_penalty"]["dS"]

print(f"After init+penalties: dH={dH}, dS={dS}")

for i in range(len(seq) - 1):
    pair = seq[i:i+2]
    if pair in NN_PARAMS:
        dH += NN_PARAMS[pair]["dH"]
        dS += NN_PARAMS[pair]["dS"]
        print(f"  {pair}: dH += {NN_PARAMS[pair]['dH']}, dS += {NN_PARAMS[pair]['dS']} -> dH={dH}, dS={dS}")
    else:
        print(f"  {pair}: NOT FOUND!")

R = 1.987
Ct = 500e-9
raw_tm = (dH / (dS + R * math.log(Ct / 4))) - 273.15
print(f"\nraw Tm = {dH} / ({dS} + {R} * ln({Ct/4})) - 273.15")
print(f"raw Tm = {raw_tm:.4f}")

# Owczarzy correction
corrected = owczarzy_correction(raw_tm, seq, 150)
print(f"Owczarzy 2004 (150mM): {corrected:.4f}")
print(f"JS Math.round: {math.floor(corrected + 0.5)}")

# Try with NEB init values
print("\n=== With NEB init values ===")
dH_neb = 200  # 0.2 kcal = 200 cal
dS_neb = -5.7
# NEB: AT penalty +2200 cal +6.9 cal/K per end
if has_5prime_AT:
    dH_neb += 2200
    dS_neb += 6.9
if has_3prime_AT:
    dH_neb += 2200
    dS_neb += 6.9

for i in range(len(seq) - 1):
    pair = seq[i:i+2]
    if pair in NN_PARAMS:
        dH_neb += NN_PARAMS[pair]["dH"]
        dS_neb += NN_PARAMS[pair]["dS"]

raw_neb = (dH_neb / (dS_neb + R * math.log(Ct / 4))) - 273.15
corrected_neb = owczarzy_correction(raw_neb, seq, 150)
print(f"NEB init: dH={dH_neb}, dS={dS_neb}")
print(f"NEB raw Tm = {raw_neb:.4f}")
print(f"NEB corrected (150mM) = {corrected_neb:.4f}")
print(f"JS Math.round: {math.floor(corrected_neb + 0.5)}")
