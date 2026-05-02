# AGENTS.md

## Environment

Work from the repository root. Use the `Makefile` for common tasks (`make install`, `make fixtures`, `make test`).

## Source of truth

Follow `spec/SPEC.md` and all documents in `spec/` as the source of truth. Do not invent undocumented behavior. When spec and code disagree, align implementation to the spec unless explicitly instructed otherwise. Spec changes should prioritize minimalism and no duplication — say a thing once, in the most appropriate document.

## Change discipline

Prefer minimal, surgical edits. Avoid unrelated refactors, formatting churn, and public-interface changes unless required. Keep I/O separate from core logic. Preserve existing module organization.

## Testing discipline

Run the narrowest relevant tests first. Add or update tests for behavior changes. Never claim tests passed unless they were actually run. If tests could not be run, say so clearly.
Use the `statgen` conda environment: `conda activate statgen && pytest tests/`
Octave tests run automatically when Octave is on the path and are skipped otherwise.

## Schema and format stability

Preserve portable file formats (`ld_shard` binary layout, BIM column order, sumstats TSV columns, `metadata.json` fields) unless explicitly instructed otherwise. If formats change, regenerate fixtures with `make fixtures`.

## Domain guardrails

Do not normalize contig labels, swap alleles, or silently drop unmatched variants. `statgen` consumes post-harmonization inputs as-is. Fail early on ambiguous chromosome naming, allele convention, or checksum mismatch.

## Reproducibility

Do not hardcode local paths or machine-specific assumptions. Fixture generation must remain deterministic (`make fixtures` produces bit-identical output on repeated runs).

## Commits

Do not commit and do not stage changes unless explicitly asked. If instructions say "implement phases X–Y with one commit per phase" that covers all those commits. Prefer informative commit messages with a short description.
