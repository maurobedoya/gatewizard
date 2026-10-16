"""Write packmol-memgen --available_lipids_all into the lipid catalog."""

from __future__ import annotations

import json
import os
import shutil
import subprocess
from pathlib import Path

DEST = Path(__file__).resolve().parents[1] / "gatewizard" / "data" / "packmol_memgen_lipids.json"


def main() -> None:
    exe = shutil.which("packmol-memgen")
    if not exe:
        raise SystemExit("packmol-memgen was not found on PATH")
    env = dict(os.environ)
    env["COLUMNS"] = "240"
    proc = subprocess.run(
        [exe, "--available_lipids_all"],
        capture_output=True,
        text=True,
        env=env,
    )
    text = "\n".join(part for part in (proc.stdout, proc.stderr) if part)
    if "Charge" not in text:
        raise SystemExit(proc.stderr or "packmol-memgen did not print a lipid table")
    lines = text.splitlines()
    version = ""
    for line in lines:
        if "VERSION:" in line:
            version = line.split("VERSION:", 1)[1].strip()
            break
    rows = []
    for line in lines:
        parts = line.split()
        if len(parts) < 3:
            continue
        try:
            charge = int(parts[1])
        except ValueError:
            continue
        name = parts[0]
        if name in {"Lipid", "-----"}:
            continue
        note = ""
        body = line
        if "***" in line:
            body, note = line.split("***", 1)
            note = note.strip()
        full_name = " ".join(body.split()[2:])
        rows.append(
            {"name": name, "charge": charge, "full_name": full_name, "note": note}
        )
    rows.sort(key=lambda row: row["name"])
    DEST.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "source": "packmol-memgen --available_lipids_all",
        "packmol_memgen_version": version,
        "count": len(rows),
        "lipids": rows,
    }
    DEST.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    notes: dict[str, int] = {}
    for row in rows:
        key = row["note"] or "core"
        notes[key] = notes.get(key, 0) + 1
    print(f"{version} {len(rows)} {DEST}")
    for key, count in sorted(notes.items()):
        print(f"{count} {key}")


if __name__ == "__main__":
    main()
