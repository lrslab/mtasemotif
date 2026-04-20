from __future__ import annotations

from abc import ABC, abstractmethod
from pathlib import Path
from typing import Optional

from mtase_motif.config import DbLayout
from mtase_motif.db.manifest import ProviderManifest


class Provider(ABC):
    name: str

    def __init__(self, layout: DbLayout) -> None:
        self.layout = layout

    @abstractmethod
    def fetch(self, *, source: Optional[Path] = None) -> ProviderManifest:
        """Fetch provider data into db_dir. Returns the updated provider manifest."""

    @abstractmethod
    def index(self) -> ProviderManifest:
        """Build indexes for already-fetched data. Returns the updated provider manifest."""

    @abstractmethod
    def status(self) -> dict[str, object]:
        """Return a lightweight status summary for db status output."""

    @abstractmethod
    def validate(self) -> None:
        """Raise MtaseMotifError if the provider payload is incomplete."""
