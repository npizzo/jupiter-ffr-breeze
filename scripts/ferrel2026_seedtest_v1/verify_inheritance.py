#!/usr/bin/env python3
"""Source-only audit; no model data or scientific calculation."""
import hashlib
import json
from pathlib import Path

here = Path(__file__).resolve().parent
repo = here.parent.parent
manifest = json.loads((here / "driver_inheritance.json").read_text())
base = repo / manifest["base_directory"]
driver = (base / "run.jl").read_text()
assert hashlib.sha256(driver.encode()).hexdigest() == manifest["base_driver_sha256"]
expected = driver
for edit in manifest["edits"]:
    assert expected.count(edit["old"]) == 1, "Ambiguous inherited-driver edit"
    expected = expected.replace(edit["old"], edit["new"])
assert expected == (here / "run.jl").read_text(), "Unreviewed driver change"
for name, digest in manifest["copied_components"].items():
    assert hashlib.sha256((base / name).read_bytes()).hexdigest() == digest
    assert (here / name).read_bytes() == (base / name).read_bytes(), name
print(json.dumps({"status": "PASS_SOURCE_INHERITANCE_ONLY", "scope": manifest["scope"],
                  "base_driver_sha256": manifest["base_driver_sha256"],
                  "new_driver_sha256": hashlib.sha256(expected.encode()).hexdigest(),
                  "copied_components": manifest["copied_components"]}, indent=2))
