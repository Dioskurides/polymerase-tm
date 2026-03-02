# Contributing to polymerase-tm

Thank you for your interest in contributing! This guide will help you get started.

## Development Setup

```bash
# Clone the repository
git clone https://github.com/Dioskurides/polymerase-tm.git
cd polymerase-tm

# Install in editable mode
pip install -e .

# Verify installation
polymerase-tm --version
```

## How to Contribute

### Reporting Bugs

- Open an [issue](https://github.com/Dioskurides/polymerase-tm/issues/new?template=bug_report.md)
- Include the primer sequence(s), polymerase used, and expected vs. actual result
- If possible, compare with the [NEB Tm Calculator](https://tmcalculator.neb.com/)

### Suggesting Features

- Open an [issue](https://github.com/Dioskurides/polymerase-tm/issues/new?template=feature_request.md) describing your use case
- Explain how the feature would help your workflow

### Submitting Changes

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/my-feature`)
3. Make your changes
4. Test against the NEB Tm Calculator for any Tm/Ta changes
5. Commit with a descriptive message
6. Push and open a Pull Request

## Code Style

- Follow PEP 8
- Add docstrings with NumPy-style parameters
- Include type hints for function signatures

## Testing

All Tm/Ta values should be verified against the [NEB Tm Calculator](https://tmcalculator.neb.com/) before submitting.

```bash
# Quick verification
python -c "from polymerase_tm import tm; print(tm('ATGTCCCTGCTCTTCTCTCGATGCAA'))"
# Expected: 72
```

## Project Structure

```
polymerase-tm/
├── src/polymerase_tm/
│   ├── __init__.py    # Core module (Tm, Ta, automation functions)
│   └── cli.py         # Command-line interface
├── pyproject.toml     # Package metadata
└── README.md          # Documentation
```

## Questions?

Open an issue or start a discussion -- we're happy to help!
