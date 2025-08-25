from pathlib import Path
from typing import Any

import yaml

def load_yaml_if_exists(path: Path) -> dict[str, Any]:
    if path is None or not path.exists():
        return {}
    with path.open() as fh:
        return yaml.safe_load(fh) or {}
