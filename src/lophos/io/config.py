from pathlib import Path
from typing import Dict, Any
import yaml

def load_yaml_if_exists(path: Path) -> Dict[str, Any]:
    if path is None or not path.exists():
        return {}
    with path.open() as fh:
        return yaml.safe_load(fh) or {}
