from __future__ import annotations

import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path

from mtase_motif.util import MtaseMotifError


@dataclass(frozen=True)
class SearchDb:
    kind: str
    db_prefix: Path


def search_rebase(
    *,
    query_faa: Path,
    out_tsv: Path,
    mmseqs_db: Path,
    blast_db: Path,
    threads: int = 1,
    max_target_seqs: int = 50,
) -> None:
    out_tsv.parent.mkdir(parents=True, exist_ok=True)
    header = "query\ttarget\tevalue\tbits\tpident\tqcov\ttcov\n"

    if not query_faa.exists() or query_faa.stat().st_size == 0:
        out_tsv.write_text(header)
        return

    db = pick_db(mmseqs_db=mmseqs_db, blast_db=blast_db)
    if db.kind == "mmseqs":
        run_mmseqs(query_faa, db.db_prefix, out_tsv, header, threads)
        return
    if db.kind == "blast":
        run_blast(
            query_faa,
            db.db_prefix,
            out_tsv,
            header,
            threads,
            max_target_seqs=max_target_seqs,
        )
        return
    raise MtaseMotifError(f"Unknown DB kind: {db.kind}")


def pick_db(*, mmseqs_db: Path, blast_db: Path) -> SearchDb:
    mmseqs_exe = shutil.which("mmseqs")
    blastp_exe = shutil.which("blastp")

    mmseqs_present = mmseqs_exe is not None and _mmseqs_db_exists(mmseqs_db)
    blast_present = blastp_exe is not None and _blast_db_exists(blast_db)

    if mmseqs_present:
        return SearchDb(kind="mmseqs", db_prefix=mmseqs_db)
    if blast_present:
        return SearchDb(kind="blast", db_prefix=blast_db)

    raise MtaseMotifError(
        "No usable REBASE protein search DB found.\n"
        "Expected either:\n"
        f"- mmseqs + MMseqs DB at {mmseqs_db}\n"
        f"- blastp + BLAST DB at {blast_db}"
    )


def run_mmseqs(query_faa: Path, target_db: Path, out_tsv: Path, header: str, threads: int) -> None:
    mmseqs = shutil.which("mmseqs")
    if mmseqs is None:
        raise MtaseMotifError("mmseqs not found in PATH (required for MMseqs REBASE search)")

    work_dir = out_tsv.parent / "mmseqs_work"
    work_dir.mkdir(parents=True, exist_ok=True)
    query_db = work_dir / "query_db"
    result_db = work_dir / "result_db"
    tmp_dir = work_dir / "tmp"
    tmp_dir.mkdir(parents=True, exist_ok=True)
    raw_out = work_dir / "convertalis.tsv"

    _run([mmseqs, "createdb", str(query_faa), str(query_db)])
    _run(
        [
            mmseqs,
            "search",
            str(query_db),
            str(target_db),
            str(result_db),
            str(tmp_dir),
            "-a",
            "--threads",
            str(threads),
        ]
    )
    _run(
        [
            mmseqs,
            "convertalis",
            str(query_db),
            str(target_db),
            str(result_db),
            str(raw_out),
            "--format-output",
            "query,target,evalue,bits,pident,qcov,tcov",
        ]
    )

    with out_tsv.open("w") as out:
        out.write(header)
        if not raw_out.exists():
            raise MtaseMotifError(f"MMseqs convertalis output missing: {raw_out}")
        out.write(raw_out.read_text())


def run_blast(
    query_faa: Path,
    blast_db: Path,
    out_tsv: Path,
    header: str,
    threads: int,
    *,
    max_target_seqs: int,
) -> None:
    blastp = shutil.which("blastp")
    if blastp is None:
        raise MtaseMotifError("blastp not found in PATH (required for BLAST REBASE search)")

    proc = subprocess.run(
        [
            blastp,
            "-query",
            str(query_faa),
            "-db",
            str(blast_db),
            "-max_target_seqs",
            str(max_target_seqs),
            "-num_threads",
            str(threads),
            "-outfmt",
            "6 qseqid sseqid evalue bitscore pident qcovs slen length",
        ],
        capture_output=True,
        text=True,
    )
    if proc.returncode != 0:
        raise MtaseMotifError(proc.stderr.strip() or "blastp failed")

    with out_tsv.open("w") as out:
        out.write(header)
        for line in proc.stdout.splitlines():
            if not line.strip():
                continue
            parts = line.split("\t")
            if len(parts) != 8:
                continue
            query, target, evalue, bits, pident, qcov_raw, slen_raw, length_raw = parts
            try:
                qcov = float(qcov_raw)
                slen = float(slen_raw)
                alen = float(length_raw)
            except ValueError:
                continue
            tcov = 0.0 if slen <= 0 else (100.0 * alen / slen)
            out.write(f"{query}\t{target}\t{evalue}\t{bits}\t{pident}\t{qcov:g}\t{tcov:g}\n")


def _mmseqs_db_exists(prefix: Path) -> bool:
    return prefix.exists() or Path(str(prefix) + ".dbtype").exists()


def _blast_db_exists(prefix: Path) -> bool:
    return any(Path(str(prefix) + ext).exists() for ext in (".pin", ".psq", ".phr"))


def _run(argv: list[str]) -> None:
    proc = subprocess.run(argv, capture_output=True, text=True)
    if proc.returncode != 0:
        raise MtaseMotifError(f"Command failed: {' '.join(argv)}\n{proc.stderr}")
