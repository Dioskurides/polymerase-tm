# Security Policy

## Supported Versions

| Version | Supported          |
| ------- | ------------------ |
| 0.4.x   | :white_check_mark: |
| < 0.4   | :x:                |

## Reporting a Vulnerability

If you discover a security vulnerability, please report it responsibly:

1. **Do not** open a public issue
2. Email the maintainer directly or use GitHub's [private vulnerability reporting](https://github.com/Dioskurides/polymerase-tm/security/advisories/new)
3. Include a description of the vulnerability and steps to reproduce

We will respond within 48 hours and work on a fix.

## Scope

This package performs local DNA sequence calculations. It does not:
- Make network requests
- Handle authentication or credentials
- Process untrusted binary data

The primary security concern is malicious input via CSV files (`from_csv`). We sanitize inputs but recommend validating CSV sources.
