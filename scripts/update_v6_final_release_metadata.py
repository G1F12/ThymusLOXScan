"""Refresh v6.0-final release inventory, state, and checksums.

This utility updates documentation metadata only. It excludes the checksum
index, manifest, and state record from the aggregate snapshot digest to avoid
self-reference, while listing the manifest and state hashes in SHA256SUMS.
"""

from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
RELEASE = ROOT / "release" / "v6.0-final"
MANIFEST = RELEASE / "RELEASE_MANIFEST.tsv"
CHECKSUMS = RELEASE / "checksums" / "SHA256SUMS.txt"
STATE = RELEASE / "governance" / "repository_version_state_v6_final.json"
POLISH_START = "32ab4f3d3b7981a0460e40185b5345cdcb23b7fa"
SPECIAL = {
    "checksums/SHA256SUMS.txt",
    "RELEASE_MANIFEST.tsv",
    "governance/repository_version_state_v6_final.json",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def release_files() -> dict[str, Path]:
    return {
        path.relative_to(RELEASE).as_posix(): path
        for path in RELEASE.rglob("*")
        if path.is_file()
    }


def snapshot_digest(files: dict[str, Path]) -> str:
    lines = [
        f"{sha256(path)}  {rel}\n"
        for rel, path in sorted(files.items())
        if rel not in SPECIAL
    ]
    return hashlib.sha256("".join(lines).encode("utf-8")).hexdigest()


def write_state(files: dict[str, Path], manifest_hash: str | None) -> None:
    state = {
        "schema": "thymusloxscan.release_state.v6",
        "record_role": "presentation_polish_snapshot_state",
        "source_branch": "human-relevance-and-wetlab-plan",
        "source_head_at_polish_start": POLISH_START,
        "dirty_tree_at_build": True,
        "release_snapshot_checksum_sha256": snapshot_digest(files),
        "snapshot_checksum_scope": (
            "all release assets excluding this state record, RELEASE_MANIFEST.tsv, "
            "and checksums/SHA256SUMS.txt to avoid self-reference"
        ),
        "final_manuscript_checksum_sha256": sha256(
            RELEASE / "manuscript" / "LOX_thymus_aging_public_preprint_v6_final.md"
        ),
        "release_manifest_checksum_sha256": manifest_hash,
        "intended_tag": "v6.0-final",
        "release_commit": None,
        "release_tag_status": "PROPOSED_NOT_CREATED",
        "presentation_polish_commit": None,
        "state_binding_status": "AWAITING_LATER_IMMUTABLE_RELEASE_FREEZE",
        "scientific_computations_executed_during_polish": 0,
        "numerical_outputs_changed_during_polish": 0,
        "figures_in_manuscript": 3,
        "final_snapshot_status": "PRESENTATION_POLISHED_WITH_DOCUMENTED_WARNINGS",
    }
    STATE.write_text(
        json.dumps(state, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
        newline="\n",
    )


def figure_row(rel: str, path: Path) -> dict[str, str]:
    return {
        "release_path": rel,
        "classification": "FIGURE",
        "public_role": "publication figure",
        "scientific_scope": "frozen transcript-level visual summary",
        "source_path": f"results/figures/v6_final/{path.name}",
        "source_authority": "frozen authoritative tables",
        "generation_method": "presentation-only assembly without scientific recomputation",
        "size_bytes": str(path.stat().st_size),
        "sha256": sha256(path),
        "inclusion_reason": "Reader-facing visualization of frozen accepted outputs",
    }


def refresh_manifest(files: dict[str, Path]) -> None:
    with MANIFEST.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = list(reader.fieldnames or [])
        rows = {row["release_path"]: row for row in reader}

    for rel, path in files.items():
        if rel.startswith("figures/"):
            rows[rel] = figure_row(rel, path)
        elif rel in rows and rel not in SPECIAL:
            rows[rel]["size_bytes"] = str(path.stat().st_size)
            rows[rel]["sha256"] = sha256(path)

    rows["checksums/SHA256SUMS.txt"].update(
        {"size_bytes": "SELF_REFERENTIAL", "sha256": "SELF_OMITTED"}
    )
    rows["RELEASE_MANIFEST.tsv"].update(
        {"size_bytes": "SELF_REFERENTIAL", "sha256": "SEE_CHECKSUM_INDEX"}
    )
    rows["governance/repository_version_state_v6_final.json"].update(
        {"size_bytes": "SEE_CHECKSUM_INDEX", "sha256": "SEE_CHECKSUM_INDEX"}
    )

    missing = set(files) - set(rows)
    if missing:
        raise RuntimeError(f"Release files missing from manifest: {sorted(missing)}")
    orphan = set(rows) - set(files)
    if orphan:
        raise RuntimeError(f"Manifest paths missing from release: {sorted(orphan)}")

    ordered = sorted(rows.values(), key=lambda row: row["release_path"].lower())
    with MANIFEST.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(ordered)


def write_checksums(files: dict[str, Path]) -> None:
    lines = [
        f"{sha256(path)}  {rel}\n"
        for rel, path in sorted(files.items())
        if rel != "checksums/SHA256SUMS.txt"
    ]
    CHECKSUMS.write_text("".join(lines), encoding="utf-8", newline="\n")


def main() -> None:
    files = release_files()
    write_state(files, manifest_hash=None)
    files = release_files()
    refresh_manifest(files)
    manifest_hash = sha256(MANIFEST)
    write_state(release_files(), manifest_hash=manifest_hash)
    write_checksums(release_files())
    print(f"Release assets: {len(release_files())}")
    print(f"Manifest SHA-256: {manifest_hash}")
    print(f"Snapshot SHA-256: {json.loads(STATE.read_text())['release_snapshot_checksum_sha256']}")


if __name__ == "__main__":
    main()
