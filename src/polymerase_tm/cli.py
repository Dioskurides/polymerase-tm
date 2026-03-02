"""Command-line interface for polymerase-tm."""

from __future__ import annotations

import argparse
import sys


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(
        prog="polymerase-tm",
        description=(
            "Polymerase Tm Calculator -- compute primer Tm and Ta for NEB polymerases."
        ),
    )
    parser.add_argument(
        "primers",
        nargs="*",
        metavar="SEQUENCE",
        help="One or two primer sequences (binding region only, no overhangs).",
    )
    parser.add_argument(
        "-p", "--polymerase",
        default="q5",
        help="Polymerase key (default: q5). Use --list to see all options.",
    )
    parser.add_argument(
        "--dmso",
        type=float,
        default=0,
        metavar="PCT",
        help="DMSO percentage (v/v) for Ta correction (-0.6 degC per 1 %%).",
    )
    parser.add_argument(
        "--dmso-check",
        action="store_true",
        help="Run DMSO recommendation analysis.",
    )
    parser.add_argument(
        "--template",
        default=None,
        metavar="FILE",
        help="GenBank template file for DMSO amplicon analysis (requires biopython).",
    )
    parser.add_argument(
        "--list",
        action="store_true",
        dest="list_poly",
        help="List all available polymerases and exit.",
    )

    args = parser.parse_args(argv)

    from polymerase_tm import (
        tm,
        ta,
        list_polymerases,
        dmso_recommendation,
        print_dmso_report,
    )

    if args.list_poly:
        print(f"\n  {'Key':<25s} {'Salt':>6s} {'Conc':>6s}  Description")
        print("  " + "-" * 80)
        for p in list_polymerases():
            print(
                f"  {p['key']:<25s} {p['buffer_salt_mM']:>4d} mM"
                f" {p['primer_conc_nM']:>4d} nM  {p['description']}"
            )
        print()
        return

    if not args.primers:
        parser.error("At least one primer sequence is required (or use --list).")

    seq1 = args.primers[0].strip().upper()
    poly = args.polymerase

    if len(args.primers) == 1:
        # Single primer -- just Tm
        t = tm(seq1, polymerase=poly)
        print(f"\n  Primer:     {seq1}")
        print(f"  Length:     {len(seq1)} nt")
        gc = (seq1.count("G") + seq1.count("C")) / len(seq1) * 100
        print(f"  GC:         {gc:.0f} %")
        print(f"  Polymerase: {poly}")
        print(f"  Tm:         {t} degC\n")

    elif len(args.primers) >= 2:
        seq2 = args.primers[1].strip().upper()
        result_ta, t1, t2 = ta(seq1, seq2, polymerase=poly, dmso_pct=args.dmso)

        print(f"\n  Primer 1:   {seq1}")
        print(f"              {len(seq1)} nt, Tm = {t1} degC")
        print(f"  Primer 2:   {seq2}")
        print(f"              {len(seq2)} nt, Tm = {t2} degC")
        print(f"  Polymerase: {poly}")
        if args.dmso > 0:
            print(f"  DMSO:       {args.dmso} %")
        print(f"  Ta:         {result_ta} degC\n")

        if args.dmso_check:
            report = dmso_recommendation(
                fwd_bind=seq1,
                rev_bind=seq2,
                template_file=args.template,
            )
            print_dmso_report(report)


if __name__ == "__main__":
    main()
