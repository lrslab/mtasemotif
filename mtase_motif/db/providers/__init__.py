from __future__ import annotations

from typing import Type

from mtase_motif.db.base import Provider
from mtase_motif.db.providers.pfam import PfamProvider
from mtase_motif.db.providers.rebase import RebaseProvider
from mtase_motif.db.providers.tigrfams import TigrfamsProvider


PROVIDERS: dict[str, Type[Provider]] = {
    "pfam": PfamProvider,
    "tigrfams": TigrfamsProvider,
    "rebase": RebaseProvider,
}


def list_providers() -> list[str]:
    return sorted(PROVIDERS.keys())


def get_provider_class(name: str) -> Type[Provider]:
    try:
        return PROVIDERS[name]
    except KeyError:
        raise ValueError(f"Unknown provider: {name}. Expected one of: {', '.join(list_providers())}")

