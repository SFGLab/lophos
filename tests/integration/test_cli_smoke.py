import subprocess, sys

def test_cli_help():
    cmd = [sys.executable, "-m", "lophos", "--help"]
    assert subprocess.run(cmd, check=False).returncode == 0
