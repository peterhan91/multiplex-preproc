#!/usr/bin/env python3
"""Unpack the downloaded archives for the active dataset.

Driven by `configs/datasets/<DATASET>.yaml`:

  unpack:
    archives:
      - {src: HE_scans.zip, dst: he_raw, type: zip}
      - {src: foo.zip,      dst: qupath_v2, type: zip-metadata,
         include_suffixes: ["server.json","project.qpproj","summary.json"]}
    visium_glob: "GBM_*.tar.gz"   # nullable

Resolves: src is relative to <root>/raw/<dataset>/; dst relative to
<root>/preproc/<dataset>/.
"""

from __future__ import annotations

import sys
import tarfile
import zipfile
from pathlib import Path

from _common import ROOT, get_logger, load_config

log = get_logger("unpack")


def unzip(src: Path, dst: Path) -> None:
    dst.mkdir(parents=True, exist_ok=True)
    if any(dst.iterdir()):
        log.info("skip unzip %s — dst %s already populated", src.name, dst)
        return
    try:
        log.info("unzip %s → %s", src.name, dst)
        with zipfile.ZipFile(src) as z:
            z.extractall(dst)
    except (zipfile.BadZipFile, EOFError) as e:
        log.warning("  archive %s not yet complete (%s) — skipping", src.name, e)


def unzip_metadata(src: Path, dst: Path, include_suffixes: list[str]) -> None:
    dst.mkdir(parents=True, exist_ok=True)
    if any(dst.iterdir()):
        log.info("skip metadata-unzip %s — dst %s already populated", src.name, dst)
        return
    try:
        log.info("unzip-metadata %s → %s (suffixes %s)", src.name, dst, include_suffixes)
        with zipfile.ZipFile(src) as z:
            for info in z.infolist():
                if info.is_dir():
                    continue
                if any(info.filename.endswith(s) for s in include_suffixes):
                    z.extract(info, dst)
    except (zipfile.BadZipFile, EOFError) as e:
        log.warning("  archive %s not yet complete (%s) — skipping", src.name, e)


def untar(src: Path, dst: Path) -> None:
    dst.mkdir(parents=True, exist_ok=True)
    if any(dst.iterdir()):
        log.info("skip untar %s — dst %s already populated", src.name, dst)
        return
    try:
        log.info("untar %s → %s", src.name, dst)
        with tarfile.open(src, "r:gz") as t:
            t.extractall(dst)
    except (tarfile.ReadError, EOFError) as e:
        log.warning("  archive %s not yet complete (%s) — skipping", src.name, e)


def main() -> int:
    dc = load_config()
    raw = dc.raw_dir
    pre = dc.preproc_dir
    if not raw.exists():
        log.error("raw dir does not exist: %s", raw)
        return 1
    spec = dc.section("unpack")
    archives = spec.get("archives") or []
    visium_glob = spec.get("visium_glob")

    for entry in archives:
        src = (raw / entry["src"]).resolve()
        dst = (pre / entry["dst"]).resolve()
        if not src.exists():
            log.warning("archive missing: %s — skipping", src.name)
            continue
        kind = entry.get("type", "zip")
        if kind == "zip":
            unzip(src, dst)
        elif kind == "zip-metadata":
            unzip_metadata(src, dst, entry.get("include_suffixes") or [])
        elif kind == "tgz":
            untar(src, dst)
        else:
            log.warning("unknown archive type %r for %s", kind, src.name)

    if visium_glob:
        visium_dir = pre / "visium"
        for tgz in sorted(raw.glob(visium_glob)):
            stem = tgz.name[:-len(".tar.gz")] if tgz.name.endswith(".tar.gz") else tgz.stem
            untar(tgz, visium_dir / stem)

    log.info("unpack done — dataset %s", dc.dataset_id)
    return 0


if __name__ == "__main__":
    sys.exit(main())
