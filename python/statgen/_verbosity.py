import sys

_LEVELS = ("quiet", "info")
_VERBOSITY = "info"


def set_verbosity(level: str) -> None:
    global _VERBOSITY
    if level not in _LEVELS:
        raise ValueError("verbosity must be one of: quiet, info")
    _VERBOSITY = level


def get_verbosity() -> str:
    return _VERBOSITY


def info(message: str) -> None:
    if _VERBOSITY == "info":
        print(message, file=sys.stderr)


def report_variant_swaps(source_path, shard_label, num_swaps, num_ambiguous, num_disabled):
    context = f"{source_path}: shard {shard_label}"
    if num_swaps:
        info(f"{context}: {num_swaps} successful allele swap(s)")
    if num_ambiguous:
        info(f"{context}: {num_ambiguous} strand-ambiguous A/T or C/G swap(s)")
    if num_disabled:
        info(f"{context}: {num_disabled} indel row(s) unmatched because indel swaps are disabled")
