from __future__ import annotations

import json
from dataclasses import asdict, dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Optional

from mtase_motif import __version__


SCHEMA_VERSION = 1


@dataclass
class DbManifest:
    schema_version: int = SCHEMA_VERSION
    created_at: str = field(
        default_factory=lambda: datetime.now(timezone.utc).isoformat(timespec="seconds")
    )
    mtase_motif_version: str = __version__
    providers: dict[str, Any] = field(default_factory=dict)


@dataclass
class ProviderArtifact:
    path: str
    sha256: str
    url: Optional[str] = None


@dataclass
class ProviderManifest:
    schema_version: int = SCHEMA_VERSION
    provider: str = ""
    version: str = ""
    fetched_at: str = field(
        default_factory=lambda: datetime.now(timezone.utc).isoformat(timespec="seconds")
    )
    license: str = "unknown"
    urls: list[str] = field(default_factory=list)
    artifacts: list[ProviderArtifact] = field(default_factory=list)
    notes: Optional[str] = None


def load_manifest(path: Path) -> DbManifest:
    data = json.loads(path.read_text())
    return DbManifest(**data)


def write_manifest(path: Path, manifest: DbManifest) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(asdict(manifest), indent=2, sort_keys=True) + "\n")


def load_provider_manifest(path: Path) -> ProviderManifest:
    data = json.loads(path.read_text())
    artifacts = [ProviderArtifact(**a) for a in data.get("artifacts", [])]
    return ProviderManifest(
        schema_version=int(data.get("schema_version", SCHEMA_VERSION)),
        provider=str(data.get("provider", "")),
        version=str(data.get("version", "")),
        fetched_at=str(data.get("fetched_at", "")),
        license=str(data.get("license", "unknown")),
        urls=list(data.get("urls", [])),
        artifacts=artifacts,
        notes=data.get("notes"),
    )


def write_provider_manifest(path: Path, manifest: ProviderManifest) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(asdict(manifest), indent=2, sort_keys=True) + "\n")
