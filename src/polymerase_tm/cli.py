"""Command-line interface for polymerase-tm."""

from __future__ import annotations

import argparse
import sys
import textwrap

EXAMPLES = """\
examples:
  Single primer Tm (default: Q5 polymerase):
    polymerase-tm ATCGATCGATCG

  Primer pair Ta:
    polymerase-tm ATGTCCCTGCTCTTCTCTCGATGCAA GTGCCTCCGAGCCAGCACC

  Different polymerase:
    polymerase-tm -p taq ATGTCCCTGCTCTTCTCTCGATGCAA GTGCCTCCGAGCCAGCACC

  Override buffer (e.g. when not using master mix):
    polymerase-tm --buffer thermopol ATGTCCCTGCTCTTCTCTCGATGCAA

  Direct salt concentration (mM):
    polymerase-tm --salt 50 ATGTCCCTGCTCTTCTCTCGATGCAA

  Ta with 3%% DMSO correction:
    polymerase-tm --dmso 3 ATGTCCCTGCTCTTCTCTCGATGCAA GTGCCTCCGAGCCAGCACC

  DMSO analysis with GenBank template:
    polymerase-tm --dmso-check --template plasmid.gbk FWD_SEQ REV_SEQ

  List all supported polymerases / buffers:
    polymerase-tm --list
    polymerase-tm --list-buffers

python API (new in v0.5):
  from polymerase_tm import primer_dimer, restriction_scan, primer_quality
  from polymerase_tm import gibson_overlaps, list_buffers, RESTRICTION_ENZYMES

  # Buffer / salt override
  tm(seq, buffer="thermopol")    # NEB buffer name
  tm(seq, salt_mM=50)            # direct mM value

  # Primer dimer check
  primer_dimer(fwd, rev)  -> {risk_level, max_score, ...}

  # Restriction scan (~120 NEB enzymes, or pass names/custom dict)
  restriction_scan(seq)                            # all enzymes
  restriction_scan(seq, enzymes=["EcoRI","BamHI"]) # by name
  restriction_scan(seq, enzymes={"My": "AACGTT"})  # custom

  # Primer quality score (0-100)
  primer_quality(seq)  -> {score, grade, issues, ...}

  # Gibson/HiFi Assembly overlap primer design
  gibson_overlaps(fwd_bind, rev_bind, left_seq, right_seq)

algorithm:
  Tm:  SantaLucia (1998) nearest-neighbor + Owczarzy (2004) salt correction
  Ta:  Polymerase-specific rules (e.g. Q5: min(Tm1,Tm2)+1, cap 72 degC)
  DMSO: -0.6 degC per 1%% (v/v)

  Verified against the official NEB Tm Calculator with 0 degC deviation.
  https://tmcalculator.neb.com/
"""


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(
        prog="polymerase-tm",
        description=(
            "Compute primer Tm/Ta for 22 NEB polymerases, plus primer dimer check, "
            "restriction site scan (~120 NEB enzymes), quality scoring, and "
            "Gibson Assembly overlap design. Reproduces the NEB Tm Calculator "
            "(SantaLucia 1998 + Owczarzy 2004) with polymerase-specific buffer "
            "conditions and Ta rules."
        ),
        epilog=EXAMPLES,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "primers",
        nargs="*",
        metavar="SEQUENCE",
        help=(
            "One or two DNA primer sequences. Use the binding region only "
            "(no Gibson/overlap overhangs). With one sequence, only Tm is "
            "calculated. With two, both Tm values and the recommended Ta "
            "are shown."
        ),
    )
    parser.add_argument(
        "-p", "--polymerase",
        default="q5",
        metavar="KEY",
        help=(
            "Polymerase product key (default: q5). Each polymerase uses "
            "its own buffer salt concentration and Ta calculation rule. "
            "Use --list to see all 22 supported products."
        ),
    )
    parser.add_argument(
        "--dmso",
        type=float,
        default=0,
        metavar="PCT",
        help=(
            "DMSO percentage (v/v) for Ta correction. Applies -0.6 degC "
            "per 1%% DMSO to the final Ta. Common values: 2-5%%."
        ),
    )
    parser.add_argument(
        "--dmso-check",
        action="store_true",
        help=(
            "Run a comprehensive DMSO recommendation analysis. Checks "
            "primer hairpins, template GC content, and GC-rich hotspots. "
            "Optionally provide --template for amplicon-level analysis."
        ),
    )
    parser.add_argument(
        "--template",
        default=None,
        metavar="FILE",
        help=(
            "GenBank (.gbk/.gb) template file for DMSO amplicon analysis. "
            "The amplicon is extracted between the primer binding sites and "
            "analysed for GC hotspots and secondary structures."
        ),
    )
    parser.add_argument(
        "--buffer",
        default=None,
        metavar="NAME",
        help=(
            "NEB buffer name to override the polymerase default. "
            "Use --list-buffers to see all available buffer names."
        ),
    )
    parser.add_argument(
        "--salt",
        type=int,
        default=None,
        metavar="mM",
        help=(
            "Direct monovalent salt concentration in mM. Overrides "
            "both --polymerase and --buffer. Useful when using a "
            "non-NEB buffer or custom reaction conditions."
        ),
    )
    parser.add_argument(
        "--list",
        action="store_true",
        dest="list_poly",
        help="List all 22 supported NEB polymerases with their buffer and Ta parameters.",
    )
    parser.add_argument(
        "--list-buffers",
        action="store_true",
        dest="list_buffers",
        help="List all 15 NEB buffers with their effective salt concentrations.",
    )
    from polymerase_tm import __version__  # type: ignore
    parser.add_argument(
        "--version",
        action="version",
        version=f"%(prog)s {__version__}",
    )

    args = parser.parse_args(argv)

    from polymerase_tm import (  # type: ignore
        tm,
        ta,
        list_polymerases,
        list_buffers,
        dmso_recommendation,
        print_dmso_report,
        pcr_protocol,
    )

    if args.list_buffers:
        print(f"\n  {'Buffer':25s} {'Salt (mM)':>10s}")
        print("  " + "-" * 38)
        for b in list_buffers():
            print(f"  {b['name']:25s} {b['salt_mM']:>7d} mM")
        print()
        return

    if args.list_poly:
        print(f"\n  {'Key':<25s} {'Salt':>6s} {'Conc':>6s}  {'Ta rule':<22s}  Description")
        print("  " + "-" * 95)
        for p in list_polymerases():
            ta_rule = p.get("ta_rule", "")
            print(
                f"  {p['key']:<25s} {p['buffer_salt_mM']:>4d} mM"
                f" {p['primer_conc_nM']:>4d} nM  {ta_rule:<22s}  {p['description']}"
            )
        print()
        return

    if not args.primers:
        parser.print_help()
        sys.exit(0)

    seq1 = args.primers[0].strip().upper()
    poly = args.polymerase

    if len(args.primers) == 1:
        # Single primer -- just Tm
        t = tm(seq1, polymerase=poly, buffer=args.buffer, salt_mM=args.salt)
        gc = (seq1.count("G") + seq1.count("C")) / len(seq1) * 100
        print(f"\n  Primer:     {seq1}")
        print(f"  Length:     {len(seq1)} nt")
        print(f"  GC:         {gc:.1f} %")
        print(f"  Polymerase: {poly}")
        print(f"  Tm:         {t} degC\n")

    elif len(args.primers) >= 2:
        seq2 = args.primers[1].strip().upper()
        result_ta, t1, t2 = ta(seq1, seq2, polymerase=poly, dmso_pct=args.dmso,
                               buffer=args.buffer, salt_mM=args.salt)

        gc1 = (seq1.count("G") + seq1.count("C")) / len(seq1) * 100
        gc2 = (seq2.count("G") + seq2.count("C")) / len(seq2) * 100

        print(f"\n  Primer 1:   {seq1}")
        print(f"              {len(seq1)} nt, GC {gc1:.1f} %, Tm = {t1} degC")
        print(f"  Primer 2:   {seq2}")
        print(f"              {len(seq2)} nt, GC {gc2:.1f} %, Tm = {t2} degC")
        print(f"  Polymerase: {poly}")
        if args.dmso > 0:
            print(f"  DMSO:       {args.dmso} %")
        print(f"  Ta:         {result_ta} degC")

        # Auto additive recommendation
        from polymerase_tm import _additive_recommendation  # type: ignore
        additive = _additive_recommendation(seq1, seq2, polymerase=poly)
        if additive["recommended"]:
            print(f"\n  [!] RECOMMENDATION: {additive['additive']} ({additive['concentration']})")
            for reason in additive["reasons"]:
                print(f"      - {reason}")
            if "note" in additive:
                print(f"      Note: {additive['note']}")
        else:
            print(f"\n  [OK] No additive required")
            
        print()

        # Generate PCR protocol if template is provided
        if args.template:
            try:
                from Bio import SeqIO
                record = SeqIO.read(args.template, "genbank")  # Assume genbank or fasta can be parsed? SeqIO.parse might be better, we'll try genbank first as it was standard
                template_seq = str(record.seq)
                
                protocol = pcr_protocol(seq1, seq2, polymerase=poly, dmso_pct=args.dmso, template=template_seq)
                amp_len = protocol.get("amplicon_length")
                
                print(f"  [PCR CYCLING PROTOCOL]")
                if amp_len:
                    print(f"  Expected Amplicon: {amp_len} bp")
                else:
                    print(f"  Expected Amplicon: Not found (check primer orientation/binding)")
                
                for step in protocol["cycling"]:
                    print(f"    - {step['step']:20s} {step['temp']:>2d} degC, {step['time']}")
                print(f"  Total Time: ~{protocol['total_time_min']} min\n")
                
            except ImportError:
                print("  [WARNING] Biopython is required to read template files. Install with: pip install biopython")
            except Exception as e:
                print(f"  [WARNING] Could not read template file: {e}")

        if args.dmso_check:
            report = dmso_recommendation(
                fwd_bind=seq1,
                rev_bind=seq2,
                template_file=args.template,
            )
            print_dmso_report(report)


if __name__ == "__main__":
    main()
