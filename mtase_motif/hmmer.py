from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Literal


@dataclass(frozen=True)
class DomainHit:
    source: str
    model_acc: str
    model_name: str
    i_evalue: float


def parse_domtblout(
    domtbl: Path,
    evalue_thresh: float,
    *,
    source: str,
    program: Literal["hmmscan", "hmmsearch"] = "hmmscan",
) -> dict[str, list[DomainHit]]:
    hits_by_query: dict[str, list[DomainHit]] = {}
    if not domtbl.exists():
        return hits_by_query
    for raw in domtbl.read_text().splitlines():
        if not raw or raw.startswith("#"):
            continue
        fields = raw.split(maxsplit=22)
        if len(fields) < 22:
            continue

        if program == "hmmscan":
            model_name = fields[0]
            model_acc_raw = fields[1]
            query_id = fields[3]
        elif program == "hmmsearch":
            model_name = fields[3]
            model_acc_raw = fields[4]
            query_id = fields[0]
        else:
            raise ValueError(f"Unsupported HMMER domtbl program: {program}")
        try:
            i_eval = float(fields[12])
        except ValueError:
            continue
        if i_eval > evalue_thresh:
            continue

        model_acc = (
            model_acc_raw.split(".", 1)[0]
            if model_acc_raw not in ("-", "")
            else model_acc_raw
        )
        hits_by_query.setdefault(query_id, []).append(
            DomainHit(source=source, model_acc=model_acc, model_name=model_name, i_evalue=i_eval)
        )
    return hits_by_query
