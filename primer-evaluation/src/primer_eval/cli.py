#!/usr/bin/env python3
"""CLI entry point for primer evaluation."""

import logging
import sys
from argparse import ArgumentParser

from primer_eval.validator import Primer3Validator, read_template_file

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
logger = logging.getLogger(__name__)


def main():
    parser = ArgumentParser(description="Evaluate DNA primer pairs for PCR experiments")
    parser.add_argument("forward_primer", help="Forward primer sequence (5'-3')")
    parser.add_argument("reverse_primer", help="Reverse primer sequence (5'-3')")
    parser.add_argument(
        "--template",
        help="Template sequence or FASTA/GenBank file path for specificity check",
    )
    parser.add_argument(
        "--max-mismatches", type=int, default=3,
        help="Maximum allowed mismatches (default: 3)",
    )
    parser.add_argument(
        "--allow-3prime-mismatches", type=int, default=1,
        help="Allowed 3'-end mismatches (default: 1)",
    )
    args = parser.parse_args()

    forward = str(args.forward_primer).upper()
    reverse = str(args.reverse_primer).upper()

    template_seq = None
    if args.template:
        import os
        if os.path.exists(args.template):
            try:
                template_seq = read_template_file(args.template)
            except Exception as e:
                logger.error("Failed to read template file: %s", e)
                sys.exit(1)
        else:
            template_seq = args.template.upper()

    validator = Primer3Validator()

    try:
        if template_seq:
            result = validator.analyze_primer_pair_with_template(
                forward, reverse, template_seq,
                args.max_mismatches, args.allow_3prime_mismatches,
            )
        else:
            result = validator.analyze_primer_pair(forward, reverse)

        print(validator.generate_report(result))
        warnings = len(result.warnings)
        print(f"\nAnalysis completed. Warnings: {warnings}")

    except ValueError as e:
        logger.error("Validation error: %s", e)
        sys.exit(1)
    except RuntimeError as e:
        logger.error("Analysis error: %s", e)
        sys.exit(1)
    except Exception as e:
        logger.error("Unexpected error: %s", e)
        sys.exit(1)


if __name__ == "__main__":
    main()
