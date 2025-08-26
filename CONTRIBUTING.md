# Contributing to LOPHOS

Thanks for your interest in improving LOPHOS! Contributions of all kinds are welcome: bug reports, feature ideas, documentation, tests, and code.

> By participating, you agree to abide by our [Code of Conduct](CODE_OF_CONDUCT.md).

---

## Ways to Contribute

* **Report bugs** — open an issue with steps to reproduce, expected vs. actual behavior, and environment info.
* **Request features** — describe the biological use case, desired CLI/outputs, and any relevant references.
* **Improve docs** — clarify wording, add examples, or expand explanations for biologists.
* **Contribute code** — small, focused pull requests with tests are easiest to review.

---

## Development Setup

LOPHOS targets **Python ≥ 3.10**.

```bash
# 1) Clone
git clone https://github.com/<org-or-user>/lophos.git
cd lophos

# 2) Create and activate a virtual env (choose one)
python -m venv .venv && source .venv/bin/activate
# or: conda create -n lophos python=3.10 -y && conda activate lophos
# or: pyenv virtualenv 3.10.14 lophos && pyenv local lophos

# 3) Install with dev extras and hooks
pip install -e ".[dev]"
pre-commit install

# 4) Verify
python -m lophos --help
pytest -q
```

**Makefile shortcuts** (optional):

```bash
make fmt   # ruff --fix + black
make lint  # ruff + black --check + mypy
make test  # pytest
```

Equivalent raw commands:

```bash
ruff check . --fix; black .
ruff check .; black --check .; mypy src/lophos
pytest -q
```

---

## Pull Request (PR) Workflow

1. **Create a feature branch** from `main`:

   ```bash
   git checkout -b feat/<short-description>
   ```
2. **Keep changes focused** and add tests for new behavior.
3. Run the full suite locally before pushing:

   ```bash
   make fmt && make lint && make test
   ```
4. **Open a PR** with:

   * What/why summary (biologist-friendly language encouraged)
   * Any CLI changes (flags, defaults)
   * Notes on performance or memory if relevant
5. Collaborate through review; please be responsive to suggestions.

**Commit messages**: Use clear, descriptive messages. Conventional style is welcome (e.g., `feat:`, `fix:`, `docs:`) but not required.

---

## Coding Guidelines

* **Typing**: add type annotations; keep `mypy` green. If third‑party stubs are missing, prefer adding `types-...` packages or isolate `# type: ignore[reason]` narrowly.
* **Formatting**: `black` (line length configured in `pyproject.toml`).
* **Linting**: `ruff` (imports ordering, unused vars, etc.).
* **Imports**: standard lib → third‑party → local; no relative imports that jump directories.
* **Logging**: use `lophos.utils.logging.get_logger()`; avoid `print` in library code.
* **Errors**: raise informative exceptions; avoid silent failures.
* **Tests**: prefer small, deterministic unit tests; add integration tests for CLI flows.
* **Docs**: write concise docstrings (one‑line summary + args/returns). Update `README.md` if user-facing behavior changes.

---

## Project Structure (reference)

```
src/lophos/
  cli.py                   # Typer CLI
  core/                    # counts, stats, calls, validation (placeholder), APA, motif
  io/                      # BAM / BED / BEDPE / YAML loaders
  report/                  # writers and summary
tests/                     # unit + integration tests
docs/                      # user docs
examples/                  # example configs/notebooks
```

### Common contribution areas

* **RG mapping tokens**: If your pipeline uses different maternal/paternal identifiers, update `src/lophos/constants.py` and add a unit test.
* **Loop counting**: If mates in your data do not share `RG`, propose/implement explicit mate lookup in `core/counts_loops.py` (with tests).
* **Validation**: Extend `core/validate_local.py` to compute real local backgrounds and Z‑scores. Include a switch to keep `local` behavior backward compatible.
* **Outputs**: Add new format writers in `report/writers.py` (e.g., TSV/CSV/BigBed) and document them.

---

## Testing

* Run `pytest -q` for quick feedback; use markers for slow/IO tests.
* Keep tests self‑contained; when possible, use small synthetic BAM/BED/BEDPE fixtures.
* Add regression tests for any bugfix.

---

## Performance

* Profile before optimizing. If you introduce significant speedups, include a brief note and (if possible) a micro‑benchmark script in `examples/`.
* Avoid unnecessary data copies; prefer iterators/generators for large scans; consider chunking when reading huge files.

---

## Security & Privacy

Do **not** post sensitive or patient‑identifying data in issues/PRs. Use synthetic or minimized examples, or share through approved private channels.

---

## License

By contributing, you agree that your contributions are licensed under the repository’s **MIT License**.

---

## Acknowledgments

LOPHOS development benefits from scientific guidance by **Dr. Joanna Borkowska** and **Prof. Dariusz Plewczyński** (Structural and Functional Genomics Laboratory). We appreciate contributions from the broader community.
