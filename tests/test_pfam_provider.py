from __future__ import annotations

import gzip
from pathlib import Path

from mtase_motif.config import DbLayout
from mtase_motif.db.providers.pfam import PfamProvider


class _FakeResponse:
    def __init__(self, payload: bytes) -> None:
        self._payload = payload

    def __enter__(self) -> "_FakeResponse":
        return self

    def __exit__(self, exc_type, exc, tb) -> None:
        return None

    def read(self) -> bytes:
        return self._payload


def _gzip_text(text: str) -> bytes:
    return gzip.compress(text.encode("utf-8"))


def test_pfam_provider_fetch_remote_subset(monkeypatch, tmp_path: Path) -> None:
    payloads = {
        "PF01555": _gzip_text(
            "HMMER3/f [3.3.2 | Nov 2020]\nNAME  N6_N4_Mtase\nACC   PF01555.25\n//\n"
        ),
        "PF00145": _gzip_text(
            "HMMER3/f [3.3.2 | Nov 2020]\nNAME  DNA_methylase\nACC   PF00145.18\n//\n"
        ),
    }

    def fake_urlopen(req):
        url = req.full_url
        for accession, payload in payloads.items():
            if accession in url:
                return _FakeResponse(payload)
        raise AssertionError(f"unexpected URL: {url}")

    monkeypatch.setattr("urllib.request.urlopen", fake_urlopen)

    layout = DbLayout(db_dir=tmp_path / "db")
    layout.ensure_dirs()
    provider = PfamProvider(layout)
    pm = provider.fetch()

    subset = layout.pfam_dir / "subset.hmm"
    text = subset.read_text()
    assert subset.exists()
    assert "PF01555.25" in text
    assert "PF00145.18" in text
    assert pm.version == "PF01555.25,PF00145.18"
    assert len(pm.urls) == 2
    assert pm.notes is not None and "InterPro API" in pm.notes


def test_pfam_provider_index_uses_existing_subset(monkeypatch, tmp_path: Path) -> None:
    layout = DbLayout(db_dir=tmp_path / "db")
    layout.ensure_dirs()
    subset = layout.pfam_dir / "subset.hmm"
    subset.write_text(
        "HMMER3/f [3.3.2 | Nov 2020]\nNAME  N6_N4_Mtase\nACC   PF01555.25\n//\n"
        "HMMER3/f [3.3.2 | Nov 2020]\nNAME  DNA_methylase\nACC   PF00145.18\n//\n"
    )
    calls = []

    def fake_run_cmd(argv, **kwargs):
        calls.append(list(argv))
        for ext in (".h3f", ".h3i", ".h3m", ".h3p"):
            Path(str(subset) + ext).write_text("")
        return None

    monkeypatch.setattr("mtase_motif.db.providers.pfam.which", lambda name: "/tmp/hmmpress")
    monkeypatch.setattr("mtase_motif.db.providers.pfam.run_cmd", fake_run_cmd)

    provider = PfamProvider(layout)
    pm = provider.index()

    assert calls == [["/tmp/hmmpress", "-f", str(subset)]]
    assert pm.version == "PF01555.25,PF00145.18"
    assert pm.notes == "Ran hmmpress on existing subset.hmm (2 models)."
    assert any(a.path.endswith("subset.hmm.h3f") for a in pm.artifacts)
