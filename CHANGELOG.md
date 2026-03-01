# Changelog

## 1.0.0 (Unreleased)

### Fixed
- **SA-mode loop counting now resolves duplicate QNAME records correctly** for MAT_PAT merged BAMs (or any pipeline that produces multiple alignments per molecule). Previously, SA-mode used a "first seen QNAME wins" shortcut, which could silently bias allele counts and label true allele-biased loops as balanced.
- **Code quality**: Removed unused variables and reduced cyclomatic complexity in core counting and I/O modules.

### Details
- In SA loop mode, LOPHOS now counts **at most one piece of evidence per QNAME per loop** and selects the best alignment per QNAME using a deterministic score: `(MAPQ, AS, -NM)`.
- If the best score is tied across conflicting alleles, the read is counted as **ambiguous**.
- **Refactored** peak and loop counting to improve maintainability:
  - Extracted per-QNAME resolution logic (`_update_per_qname_mates`, `_update_per_qname_peaks`) to eliminate nested functions and reduce cognitive load.
  - Reduced cyclomatic complexity in `_counts_for_single_loop_mates`, `_counts_for_single_loop_sa`, and `count_peaks`.
  - Decomposed `write_loops()` into focused normalization helpers (`_normalize_informative_counts`, `_normalize_evidence_columns`, etc.).
- All changes maintain **100% backward compatibility** — outputs and CLI behavior unchanged.

### Added
- Unit tests covering SA-mode per-QNAME resolution logic.
- Comprehensive **Code Quality & Development** section in README documenting linting, formatting, and testing tools.
- Detailed project structure documentation with references to refactored modules.
