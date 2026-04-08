# SPECIFICITY README

This document describes the specificity checking feature of the primer evaluation tool.

## Overview

The tool performs DNA primer pair evaluation including specificity checking against a template sequence.

## Features

### Basic Primer Evaluation
- Tm calculation
- GC content analysis
- Hairpin structure assessment
- Homodimer/heterodimer analysis
- 3'-end stability check

### Specificity Checking
- Primer-template alignment via sliding window
- Mismatch tolerance configuration
- 3'-end specificity (critical for DNA polymerase extension)
- Potential non-specific amplification product prediction
- Multi-binding site detection

## Usage

### Python API

```python
from primer_eval import Primer3Validator

validator = Primer3Validator()

# Basic analysis
analysis = validator.analyze_primer_pair("ATGCCCTGAGCTAAAGCTG", "GTGAGCTTTGTCTCGGTGA")

# With specificity check
analysis = validator.analyze_primer_pair_with_template(
    "ATGCCCTGAGCTAAAGCTG",
    "GTGAGCTTTGTCTCGGTGA",
    "ATCGATCGATCG..."
)
```

### CLI

```bash
python -m primer_eval.cli <forward> <reverse> --template <sequence_or_file>
```
