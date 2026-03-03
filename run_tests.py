import sys
from pathlib import Path
import pytest

sys.path.insert(0, str(Path(__file__).parent / "src"))

if __name__ == "__main__":
    sys.exit(pytest.main(sys.argv[1:] if len(sys.argv) > 1 else ["tests"]))
