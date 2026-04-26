#!/usr/bin/env python3
"""Resolve and verify CODEX channel→marker mapping per slide (config-driven).

Panel-source priority is set by the dataset config:

  panel:
    source: qupath_v2 | csv | inline | <constant-key from scripts/panels.PANELS>
    qupath_project_path: <relative-to-preproc-dataset-root>
    csv_path: <relative-or-absolute>
    inline: [DAPI, ...]
    fallback_constant: <constant-key>

  panel_verify:
    expected_dapi_at_index_0: bool
    expected_total_channels: int | null
    cycle_structure_constant: <key from scripts/panels.PANEL_CYCLES> | null

Verifies length, DAPI position (when configured), and per-channel reporter
consistency (when a cycle constant is named). Writes:

  index.json[slides][sid].codex.markers / .panel_source
  preproc/<dataset>/panels/<sid>.csv
  preproc/<dataset>/panels/_audit.csv

Updates patches/<sid>.h5 attrs `channels` if H5 already exists.
"""

from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

import h5py
import numpy as np

from _common import ROOT, get_logger, load_config, load_index, save_index
from panels import PANELS, PANEL_CYCLES, canonical_fluor

log = get_logger("panels")


# ---------------------------------------------------------------- panel sources

def _qupath_project_map(project_path: Path) -> dict[str, str]:
    """Return {sample_id_normalized: entry_id_str} for QuPath v2 projects."""
    proj = project_path / "project.qpproj"
    if not proj.exists():
        log.warning("QuPath project file not found: %s", proj)
        return {}
    d = json.loads(proj.read_text())
    out: dict[str, str] = {}
    for img in d.get("images", []):
        eid = img.get("entryID") or img.get("imageID") or img.get("id")
        name = img.get("imageName") or img.get("name") or ""
        if not eid or not name or name.endswith("_HE"):
            continue
        norm = name
        for suf in ("_A", "_B"):
            if norm.endswith(suf):
                norm = norm[: -len(suf)]
        if norm.endswith("_v6"):
            norm = norm[: -len("_v6")]
        out[norm] = str(eid)
    return out


def _qupath_panel(project_path: Path, entry_id: str) -> list[str] | None:
    sj = project_path / "data" / entry_id / "server.json"
    if not sj.exists():
        return None
    d = json.loads(sj.read_text())

    def find_chans(o):
        if isinstance(o, dict):
            if "channels" in o and isinstance(o["channels"], list):
                return o["channels"]
            for v in o.values():
                r = find_chans(v)
                if r:
                    return r
        elif isinstance(o, list):
            for x in o:
                r = find_chans(x)
                if r:
                    return r
        return None

    chans = find_chans(d) or []
    return [c.get("name") for c in chans if c.get("name") is not None]


def _csv_panel(csv_path: Path) -> list[str]:
    """Read marker names from a CSV. Tries common column names."""
    with csv_path.open() as f:
        reader = csv.DictReader(f)
        rows = list(reader)
    for col in ("marker", "Marker", "name", "Name", "channel", "Channel", "Target"):
        if rows and col in rows[0]:
            return [r[col].strip() for r in rows if r.get(col)]
    raise SystemExit(f"{csv_path}: cannot find a marker column (tried marker/name/channel/Target)")


# ---------------------------------------------------------------- verification

def _verify(sid: str, markers: list[str], ome_channels: list[str],
            *, expected_dapi: bool, expected_total: int | None,
            cycle_key: str | None) -> list[str]:
    warnings: list[str] = []
    if len(markers) != len(ome_channels):
        raise SystemExit(
            f"{sid}: panel size mismatch — markers={len(markers)} OME channels={len(ome_channels)}"
        )
    if expected_total is not None and len(markers) != expected_total:
        raise SystemExit(f"{sid}: panel size {len(markers)} != expected_total_channels {expected_total}")
    if expected_dapi and markers and markers[0] != "DAPI":
        raise SystemExit(f"{sid}: ch0 must be DAPI, got {markers[0]!r}")

    if cycle_key:
        cycles = PANEL_CYCLES.get(cycle_key)
        if cycles is None:
            raise SystemExit(f"{sid}: cycle_structure_constant {cycle_key!r} not in panels.PANEL_CYCLES")
        if len(cycles) != len(markers):
            raise SystemExit(f"{sid}: cycle structure {cycle_key!r} length {len(cycles)} != panel {len(markers)}")
        for i, (cyc, expected) in enumerate(cycles):
            ome_canon = canonical_fluor(ome_channels[i])
            if ome_canon != expected:
                raise SystemExit(
                    f"{sid}: ch{i} fluorophore mismatch — panel expects {expected} "
                    f"(cycle {cyc}, marker {markers[i]!r}), OME has {ome_channels[i]!r}"
                )

    return warnings


def _write_panel_csv(sid: str, markers: list[str], ome_channels: list[str], source: str,
                     cycle_key: str | None, out_dir: Path) -> Path:
    out_dir.mkdir(parents=True, exist_ok=True)
    fp = out_dir / f"{sid}.csv"
    cycles = PANEL_CYCLES.get(cycle_key) if cycle_key else None
    with fp.open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["ch_idx", "marker", "cycle", "fluorophore_expected", "fluorophore_ome", "source"])
        for i, m in enumerate(markers):
            cyc, expected = (cycles[i] if cycles else (None, None))
            w.writerow([i, m, cyc if cyc is not None else "", expected or "", ome_channels[i], source])
    return fp


def _update_h5_channels(sid: str, markers: list[str], patches_dir: Path) -> bool:
    h5 = patches_dir / f"{sid}.h5"
    if not h5.exists():
        return False
    with h5py.File(h5, "a") as f:
        f.attrs["channels"] = np.array(markers, dtype=h5py.string_dtype())
    return True


# ---------------------------------------------------------------- main

def _resolve_panel(sid: str, ome_channels: list[str], dc) -> tuple[list[str] | None, str]:
    spec = dc.section("panel")
    source = spec.get("source", "auto")

    # 1) QuPath v2 project
    if source in ("qupath_v2", "auto") and spec.get("qupath_project_path"):
        proj = (dc.preproc_dir / spec["qupath_project_path"]).resolve()
        qmap = _qupath_project_map(proj)
        if sid in qmap:
            ml = _qupath_panel(proj, qmap[sid])
            if ml and len(ml) == len(ome_channels):
                return ml, "qupath_v2"

    # 2) CSV
    if source in ("csv", "auto") and spec.get("csv_path"):
        csv_path = Path(spec["csv_path"])
        if not csv_path.is_absolute():
            csv_path = (dc.preproc_dir / csv_path).resolve()
        if csv_path.exists():
            ml = _csv_panel(csv_path)
            if ml and len(ml) == len(ome_channels):
                return ml, f"csv:{csv_path.name}"

    # 3) Inline
    if source == "inline" and spec.get("inline"):
        ml = list(spec["inline"])
        if len(ml) == len(ome_channels):
            return ml, "inline"

    # 4) Named constant (fallback)
    fallback = spec.get("fallback_constant") or (source if source in PANELS else None)
    if fallback and fallback in PANELS and len(PANELS[fallback]) == len(ome_channels):
        return list(PANELS[fallback]), f"constant:{fallback}"

    return None, "unresolved"


def main() -> int:
    dc = load_config()
    index = load_index(dc)
    slides = index.get("slides", {})
    if not slides:
        log.error("index.json empty — run 01_index first")
        return 1

    pv = dc.section("panel_verify")
    expected_dapi = bool(pv.get("expected_dapi_at_index_0"))
    expected_total = pv.get("expected_total_channels")
    cycle_key = pv.get("cycle_structure_constant")

    audit_rows: list[dict] = []
    for sid in sorted(slides):
        codex = slides[sid].get("codex")
        if not codex:
            continue
        ome_channels = list(codex.get("channels") or [])
        markers, source = _resolve_panel(sid, ome_channels, dc)
        if markers is None:
            log.error("%s: cannot resolve panel (n_ome=%d)", sid, len(ome_channels))
            audit_rows.append({"sid": sid, "source": "none", "n_ome": len(ome_channels), "warnings": "unresolved"})
            continue

        warnings = _verify(sid, markers, ome_channels,
                           expected_dapi=expected_dapi, expected_total=expected_total, cycle_key=cycle_key)
        for w in warnings:
            log.warning("%s: %s", sid, w)

        codex["markers"] = markers
        codex["panel_source"] = source
        csv_path = _write_panel_csv(sid, markers, ome_channels, source, cycle_key, dc.panels_dir)
        h5_updated = _update_h5_channels(sid, markers, dc.patches_dir)
        log.info(
            "%s: source=%s ch0..3=%s ch%d..%d=%s csv=%s h5_updated=%s",
            sid, source, markers[:4], len(markers)-2, len(markers)-1, markers[-2:],
            csv_path.relative_to(ROOT), h5_updated,
        )
        audit_rows.append({"sid": sid, "source": source, "n_ome": len(ome_channels),
                           "warnings": "; ".join(warnings) or "ok"})

    save_index(index, dc)
    audit_fp = dc.panels_dir / "_audit.csv"
    audit_fp.parent.mkdir(parents=True, exist_ok=True)
    with audit_fp.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["sid", "source", "n_ome", "warnings"])
        w.writeheader()
        w.writerows(audit_rows)
    log.info("audit → %s", audit_fp)
    return 0


if __name__ == "__main__":
    sys.exit(main())
