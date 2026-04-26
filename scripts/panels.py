"""Panel registry — named cohort-level panels + cycle structures.

Indexed by `<key>` so dataset configs reference them as e.g.
``panel.fallback_constant: greenwald43`` and
``panel_verify.cycle_structure_constant: greenwald43_cycles``.

To add a new published panel: append a constant + register it in PANELS /
PANEL_CYCLES below. Per-slide panel resolution (the QuPath / CSV / inline
sources) is handled by 01b_panels.py — this module is just the cohort-level
fallback & verification reference.
"""

from __future__ import annotations

# ============================================================ Greenwald 2024 *Cell*
# Verified against:
# - paper/mmc1.pdf Table S3 (40 antibodies + DAPI)
# - preproc/greenwald-gbm-codex/qupath_v2/CODEX_qupath_project/data/<id>/server.json
# OME-XML reporter order within each cycle: Atto550 → AF488 (FITC) → AF647 (Cy5).
GREENWALD_PANEL = [
    "DAPI",                              # ch00 — nuclear stain
    "CD69",   "CD45",   "CD3",           # ch01-03  cycle  1
    "Ki67",   "CD90",   "CD11c",         # ch04-06  cycle  2
    "HLA-DR", "CD8",    "CD4",           # ch07-09  cycle  3
    "CD31",   "FN1",    "SOX10",         # ch10-12  cycle  4
    "PDPN",   "VIM",    "PLP1",          # ch13-15  cycle  5
    "CD279",  "GFAP",   "EGFR",          # ch16-18  cycle  6  (CD279 = PD-1)
    "DCX",    "SOX4",   "CD14",          # ch19-21  cycle  7
    "NeuN",   "AQP4",   "OLIG2",         # ch22-24  cycle  8
    "GLUT1",  "PDGFRA", "CA9",           # ch25-27  cycle  9
    "p53",    "BCAN",   "CD44",          # ch28-30  cycle 10
    "NDRG1",  "GAP43",  "CD163",         # ch31-33  cycle 11
    "APOD",   "SOX2",                    # ch34-35  cycle 12 (no Atto550)
    "APOE",   "S100B",                   # ch36-37  cycle 13 (no AF488)
    "CHI3L1", "FOXP3",  "CD1c",          # ch38-40  cycle 14 (FOXP3, CD1c not in Table S3)
    "MAP2",   "CD19",                    # ch41-42  cycle 15 (no AF488)
]
assert len(GREENWALD_PANEL) == 43

GREENWALD_PANEL_CYCLES: list[tuple[int, str]] = [
    (0,  "DAPI"),
    (1,  "Atto550"), (1,  "AF488"), (1,  "AF647"),
    (2,  "Atto550"), (2,  "AF488"), (2,  "AF647"),
    (3,  "Atto550"), (3,  "AF488"), (3,  "AF647"),
    (4,  "Atto550"), (4,  "AF488"), (4,  "AF647"),
    (5,  "Atto550"), (5,  "AF488"), (5,  "AF647"),
    (6,  "Atto550"), (6,  "AF488"), (6,  "AF647"),
    (7,  "Atto550"), (7,  "AF488"), (7,  "AF647"),
    (8,  "Atto550"), (8,  "AF488"), (8,  "AF647"),
    (9,  "Atto550"), (9,  "AF488"), (9,  "AF647"),
    (10, "Atto550"), (10, "AF488"), (10, "AF647"),
    (11, "Atto550"), (11, "AF488"), (11, "AF647"),
    (12, "AF488"),   (12, "AF647"),
    (13, "Atto550"), (13, "AF647"),
    (14, "Atto550"), (14, "AF488"), (14, "AF647"),
    (15, "Atto550"), (15, "AF647"),
]
assert len(GREENWALD_PANEL_CYCLES) == 43

# ============================================================ Registries
PANELS: dict[str, list[str]] = {
    "greenwald43": GREENWALD_PANEL,
}

PANEL_CYCLES: dict[str, list[tuple[int, str]]] = {
    "greenwald43_cycles": GREENWALD_PANEL_CYCLES,
    "greenwald43":        GREENWALD_PANEL_CYCLES,  # alias
}


# ============================================================ Fluorophore aliases

FLUOR_ALIASES: dict[str, set[str]] = {
    "DAPI":    {"DAPI", "Hoechst"},
    "Atto550": {"ATTO 550", "Atto550", "ATTO550"},
    "AF488":   {"FITC", "AF488", "Alexa 488", "Alexa Fluor 488"},
    "AF647":   {"Cy5", "AF647", "Alexa 647", "Alexa Fluor 647"},
}


def canonical_fluor(name: str | None) -> str | None:
    """Map an OME-XML fluorophore label → canonical reporter name."""
    if not name:
        return None
    n = name.strip()
    n_lower = n.lower()
    for canon, aliases in FLUOR_ALIASES.items():
        if n in aliases or n_lower in {a.lower() for a in aliases}:
            return canon
    return None
