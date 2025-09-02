# from your lophos repo root
git checkout -b fix/bedpe-extra-cols

# (you’ve already edited src/lophos/io/bedpe.py to the new version)

# format + lint + type-check + tests
ruff check . --fix
black .
mypy src/lophos
pytest -q

# (optional but recommended) bump patch version
# edit pyproject.toml:  version = "0.1.2"
# and set src/lophos/__init__.py:
#   __all__ = []
#   __version__ = "0.1.2"

# run hooks once more
pre-commit run --all-files



git add src/lophos/io/bedpe.py pyproject.toml src/lophos/__init__.py
git commit -m "io: accept >10 BEDPE columns; keep extras as extra1..N (fix length mismatch)"
git push -u origin fix/bedpe-extra-cols


git checkout main
git merge --no-ff fix/bedpe-extra-cols -m "Merge fix: tolerant BEDPE reader"
git push
git tag v0.1.2 -m "BEDPE reader: keep extra columns; fix length mismatch"
git push --tags
