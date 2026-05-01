CANONICAL_CHR_ORDER = [str(i) for i in range(1, 23)] + ["X"]
CHR_RANK = {c: i for i, c in enumerate(CANONICAL_CHR_ORDER)}


def validate_requested_shards(shards, available_labels, where: str) -> list[str]:
    labels = list(available_labels)
    if shards is None:
        return labels
    if isinstance(shards, str):
        raise ValueError(f"{where}: shards must be a non-empty list of unique canonical contig labels")

    requested = list(shards)
    if not requested:
        raise ValueError(f"{where}: shards must be a non-empty list of unique canonical contig labels")

    seen = set()
    prev_rank = None
    for label in requested:
        if label not in CHR_RANK:
            raise ValueError(f"{where}: unsupported shard label {label!r}; expected canonical labels 1-22 or X")
        if label in seen:
            raise ValueError(f"{where}: duplicate shard label {label!r} in shards")
        rank = CHR_RANK[label]
        if prev_rank is not None and rank <= prev_rank:
            raise ValueError(f"{where}: shards must be in canonical subsequence order")
        seen.add(label)
        prev_rank = rank
        if label not in labels:
            raise ValueError(f"{where}: requested shard {label!r} is not present")
    return requested
