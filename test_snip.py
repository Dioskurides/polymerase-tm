from polymerase_tm import (
    batch_tm,
    optimal_binding_length,
    check_pair,
    pcr_protocol,
    reverse_complement,
)

print("--- BATCH TM ---")
results = batch_tm(['ATCGATCGATCG', 'GCGCGCGCGCGC', 'AATTCCGGAATT'])
for r in results:
    print(f"{r['sequence']}: Tm={r['tm']} degC, GC={r['gc_pct']}%")

print("\n--- OPTIMAL BINDING LENGTH ---")
result = optimal_binding_length('ATGTCCCTGCTCTTCTCTCGATGCAA', target_tm=65)
print(f"{result['binding_seq']} ({result['length']} nt, Tm={result['tm']})")

print("\n--- CHECK PAIR ---")
pair = check_pair('ATGTCCCTGCTCTTCTCTCGATGCAA', 'GTGCCTCCGAGCCAGCACC')
print(f"Ta={pair['ta']}, compatible={pair['compatible']}")
if pair["additive"]["recommended"]:
    print(f"Use {pair['additive']['additive']} ({pair['additive']['concentration']})")

print("\n--- PCR PROTOCOL ---")
protocol = pcr_protocol(
    'ATGTCCCTGCTCTTCTCTCGATGCAA',
    'GTGCCTCCGAGCCAGCACC',
    template='ATGTCCCTGCTCTTCTCTCGATGCAAGTGCCTCCGAGCCAGCACC',
)
for step in protocol['cycling']:
    print(f"{step['step']:25s} {step['temp']} degC  {step['time']}")
