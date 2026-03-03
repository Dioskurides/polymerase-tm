import sys
from pathlib import Path

# Add src to path to import gel module
src_path = str(Path(__file__).resolve().parent.parent / "src")
if src_path not in sys.path:
    sys.path.insert(0, src_path)

from polymerase_tm.gel import _get_migration_distance_cm, LADDERS

def test_migration():
    for name, ladder in LADDERS.items():
        print(f"\n=== LADDER: {name} ===")
        print(f"{'Band (bp)':>10s} | {'Distance (cm)':>15s} | {'Gap to prev':>15s}")
        print("-" * 48)
        prev_dist = None
        for bp, intensity in ladder:
            dist = _get_migration_distance_cm(bp, 1.0, 110.0, 60.0)
            if prev_dist is None:
                gap = 0.0
            else:
                gap = dist - prev_dist
                
            print(f"{bp:10d} | {dist:15.2f} cm | {gap:15.2f} cm")
            
            prev_dist = dist

if __name__ == "__main__":
    test_migration()
