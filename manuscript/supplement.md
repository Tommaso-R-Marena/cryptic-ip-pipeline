# Supplement

## S1. Configuration

The full default merged configuration is available via
`crypticip --help` and `config/default.yaml`. Run-specific hashes are
recorded in every output's `metadata.config_hash` so audit trails are
trivial.

## S2. Tier definitions

| Tier  | Filter pass count | Score   | APBS         |
|-------|-------------------|---------|--------------|
| 1     | 6 / 6             | ≥ 0.65  | ok           |
| 2     | ≥ 5 / 6           | ≥ 0.55  | any          |
| 3     | any               | ≥ 0.40  | any          |
| Reject| any               | < 0.40  | any          |

## S3. Validation set details

See `config/validation.yaml`.

## S4. Hardware

The "fast" screening mode (no APBS) processes a yeast-sized proteome
(~6 000 entries) in ~2 hours on 8 cores. The "full" mode is closer to
8–24 hours per organism on the same hardware.
