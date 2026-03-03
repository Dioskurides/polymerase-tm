"""Command-line interface for polymerase-tm."""

from __future__ import annotations

import argparse
import sys
import textwrap

EXAMPLES = """\
examples:
  Single primer Tm (default: Q5 polymerase):
    polymerase-tm ATCGATCGATCG

  Primer pair Ta + PCR protocol:
    polymerase-tm ATGTCCCTGCTCTTCTCTCGATGCAA GTGCCTCCGAGCCAGCACC

  Different polymerase:
    polymerase-tm -p taq ATGTCCCTGCTCTTCTCTCGATGCAA GTGCCTCCGAGCCAGCACC

  Override buffer (e.g. when not using master mix):
    polymerase-tm --buffer thermopol ATGTCCCTGCTCTTCTCTCGATGCAA

  Direct salt concentration (mM):
    polymerase-tm --salt 50 ATGTCCCTGCTCTTCTCTCGATGCAA

  Ta with 3%% DMSO correction:
    polymerase-tm --dmso 3 ATGTCCCTGCTCTTCTCTCGATGCAA GTGCCTCCGAGCCAGCACC

  PCR protocol with template (auto amplicon length & extension time):
    polymerase-tm FWD_SEQ REV_SEQ --template plasmid.gbk

  DMSO analysis with GenBank template:
    polymerase-tm --dmso-check --template plasmid.gbk FWD_SEQ REV_SEQ

  Virtual agarose gel:
    polymerase-tm FWD_SEQ REV_SEQ --template plasmid.gbk --plot-gel gel.png
    polymerase-tm --plot-gel gel.png --ladder 100bp --plot-gel-sizes 150 400
    polymerase-tm --plot-gel gel.png --plot-gel-sizes 3000 3000 --topology coiled nicked

  Gel physics (custom agarose/voltage/time):
    polymerase-tm --plot-gel gel.png --plot-gel-sizes 1500 --agarose 1.5 --voltage 110 --time 90

  List all supported polymerases / buffers:
    polymerase-tm --list
    polymerase-tm --list-buffers

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
        "--plot-gel",
        nargs="?",
        const="virtual_gel.png",
        default=None,
        metavar="PATH",
        help="Generate a virtual agarose gel image if an amplicon size is known (default: virtual_gel.png).",
    )
    parser.add_argument(
        "--plot-gel-sizes",
        nargs="+",
        type=int,
        default=None,
        metavar="BP",
        help="Provide additional amplicon sizes (in bp) to render them in separate lanes next to the ladder.",
    )
    parser.add_argument(
        "--ladder",
        default="1kb_plus",
        choices=["1kb_plus", "1kb", "100bp", "50bp", "low_mw", "pcr_marker"],
        help="Choose the DNA ladder for the virtual gel. Default: 1kb_plus",
    )
    parser.add_argument(
        "--topology",
        nargs="+",
        choices=["linear", "coiled", "nicked"],
        default=["linear"],
        help="Topological conformation of the amplicons/plasmids. Use a list of topologies matching --plot-gel-sizes. Default: linear",
    )
    parser.add_argument(
        "--agarose",
        type=float,
        default=1.0,
        metavar="PCT",
        help="Agarose percentage (w/v) for virtual gel physics. Default: 1.0",
    )
    parser.add_argument(
        "--voltage",
        type=float,
        default=110.0,
        metavar="V",
        help="Voltage for virtual gel physics. Default: 110.0",
    )
    parser.add_argument(
        "--time",
        type=float,
        default=60.0,
        metavar="MIN",
        help="Runtime in minutes for virtual gel physics. Default: 60.0",
    )
    parser.add_argument(
        "--gel-length",
        type=float,
        default=12.0,
        metavar="CM",
        help="Total length of the virtual gel in cm. Default: 12.0",
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
        help="List all 17 NEB buffers with their effective salt concentrations.",
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
        plot_virtual_gel,
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
        # Check if the user just wants to render a virtual gel without any primers
        if args.plot_gel and getattr(args, "plot_gel_sizes", None):
            try:
                from polymerase_tm.gel import plot_virtual_gel  # type: ignore
                plot_virtual_gel(
                    args.plot_gel_sizes, 
                    ladder_name=args.ladder, 
                    output_path=args.plot_gel,
                    agarose_pct=args.agarose,
                    voltage=args.voltage,
                    time_min=args.time,
                    gel_length_cm=args.gel_length,
                    amplicon_topologies=args.topology
                )
                print(f"\n  [GEL] Saved virtual gel to: {args.plot_gel}\n")
            except ImportError as e:
                print(f"\n  [GEL] {e}\n")
            sys.exit(0)
            
        parser.print_help()
        sys.exit(0)

    seq1 = args.primers[0].strip().upper()
    poly = args.polymerase

    # Validate DNA sequences
    for i, seq in enumerate(args.primers, 1):
        invalid = set(seq.strip().upper()) - set("ATGC")
        if invalid:
            print(
                f"\n  [ERROR] Primer {i} contains invalid character(s): "
                f"{', '.join(sorted(invalid))}.\n"
                f"  Only A, T, G, C are allowed (binding region only, no overhangs).\n"
            )
            sys.exit(1)

    if len(args.primers) == 1:
        if args.template:
            print("\n  [WARNING] --template is ignored for single-primer Tm calculation.\n")
        if args.dmso_check:
            print("\n  [WARNING] --dmso-check requires a primer pair. Ignored.\n")
        # Single primer -- just Tm
        t = tm(seq1, polymerase=poly, buffer=args.buffer, salt_mM=args.salt)
        gc = (seq1.count("G") + seq1.count("C")) / len(seq1) * 100
        print(f"\n  Primer:     {seq1}")
        print(f"  Length:     {len(seq1)} nt")
        print(f"  GC:         {gc:.1f} %")
        print(f"  Polymerase: {poly}")
        if args.buffer:
            print(f"  Buffer:     {args.buffer}")
        if args.salt:
            print(f"  Salt:       {args.salt} mM")
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
        if args.buffer:
            print(f"  Buffer:     {args.buffer}")
        if args.salt:
            print(f"  Salt:       {args.salt} mM")
        if args.dmso > 0:
            print(f"  DMSO:       {args.dmso} %")
        print(f"  Ta:         {result_ta} degC")

        # Auto additive recommendation
        from polymerase_tm import additive_recommendation  # type: ignore
        additive = additive_recommendation(seq1, seq2, polymerase=poly)
        if additive["recommended"]:
            print(f"\n  [!] RECOMMENDATION: {additive['additive']} ({additive['concentration']})")
            for reason in additive["reasons"]:
                print(f"      - {reason}")
            if "note" in additive:
                print(f"      Note: {additive['note']}")
        else:
            print(f"\n  [OK] No additive required")
            
        print()

        # Generate PCR protocol (always for primer pairs; template is optional)
        template_seq = None
        template_topology = "linear"
        
        if args.template:
            try:
                from Bio import SeqIO  # type: ignore
                record = SeqIO.read(args.template, "genbank")
                template_seq = str(record.seq)
                template_topology = record.annotations.get("topology", "linear")
            except ImportError:
                print("  [WARNING] Biopython is required to read template files. Install with: pip install polymerase-tm[bio]\n")
            except Exception as e:
                print(f"  [WARNING] Could not read template file: {e}\n")

        protocol = pcr_protocol(seq1, seq2, polymerase=poly, dmso_pct=args.dmso,
                                template=template_seq, buffer=args.buffer,
                                salt_mM=args.salt)
        amp_len = protocol.get("amplicon_length")
        
        print(f"  [PCR CYCLING PROTOCOL]")
        if amp_len:
            print(f"  Expected Amplicon: {amp_len} bp (linear)")
        elif template_seq:
            template_len = len(template_seq)
            print(f"  Expected Amplicon: Not found (check primer orientation/binding)")
            print(f"  Template Size: {template_len} bp ({template_topology})")
        
        for step in protocol["cycling"]:
            temp = step['temp']
            temp_str = f"{temp:>2d}" if isinstance(temp, (int, float)) else f"{temp}"
            print(f"    - {step['step']:20s} {temp_str} degC, {step['time']}")
        print(f"  Total Time: ~{protocol['total_time_min']} min\n")
        
        # Virtual gel (requires template or --plot-gel-sizes)
        sizes_to_plot = []
        tops_to_plot = []
        
        if amp_len:
            sizes_to_plot.append(amp_len)
            tops_to_plot.append("linear")  # PCR products are always linear
        elif template_seq:
            sizes_to_plot.append(len(template_seq))
            # Biopython uses "circular", map to our "coiled" for plasmids
            tops_to_plot.append("coiled" if template_topology == "circular" else template_topology)
            
        if getattr(args, "plot_gel_sizes", None):
            sizes_to_plot.extend(args.plot_gel_sizes)
            
            # Extend topologies from args, pad with linear if needed
            custom_tops = args.topology if hasattr(args, "topology") else []
            for i in range(len(args.plot_gel_sizes)):
                if i < len(custom_tops):
                    tops_to_plot.append(custom_tops[i])
                else:
                    tops_to_plot.append("linear")
            
        if args.plot_gel and sizes_to_plot:
            try:
                plot_virtual_gel(
                    sizes_to_plot, 
                    ladder_name=args.ladder, 
                    output_path=args.plot_gel,
                    agarose_pct=args.agarose,
                    voltage=args.voltage,
                    time_min=args.time,
                    gel_length_cm=args.gel_length,
                    amplicon_topologies=tops_to_plot
                )
                print(f"  [GEL] Saved virtual gel to: {args.plot_gel} (Ladder: {args.ladder})\n")
            except ImportError as e:
                print(f"  [GEL] {e}\n")

        if args.dmso_check:
            report = dmso_recommendation(
                fwd_bind=seq1,
                rev_bind=seq2,
                template_file=args.template,
            )
            print_dmso_report(report)


if __name__ == "__main__":
    main()
