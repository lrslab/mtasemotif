from __future__ import annotations

import gzip
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Iterator, Optional, TextIO


def extract_dna_entity_sequences(path: Path) -> list[str]:
    """
    Best-effort extraction of DNA polymer entity sequences from an mmCIF file.

    This intentionally avoids heavy dependencies; it parses the `_entity_poly` loop and
    returns canonicalized sequences (A/C/G/T/N) for entities whose type indicates DNA.
    """
    dna: list[str] = []
    with _open_text(path) as f:
        stream = _TokenStream(_iter_tokens(f))
        while stream.has_next():
            tok = stream.next()
            if tok != "loop_":
                continue
            cols: list[str] = []
            while stream.has_next() and stream.peek().startswith("_"):
                cols.append(stream.next())
            if not cols:
                continue
            if not any(c.startswith("_entity_poly.") for c in cols):
                _skip_loop_rows(stream, width=len(cols))
                continue

            idx_type = _col_index(cols, "_entity_poly.type")
            if idx_type is None:
                _skip_loop_rows(stream, width=len(cols))
                continue

            idx_seq = _col_index(cols, "_entity_poly.pdbx_seq_one_letter_code_can")
            if idx_seq is None:
                idx_seq = _col_index(cols, "_entity_poly.pdbx_seq_one_letter_code")
            if idx_seq is None:
                _skip_loop_rows(stream, width=len(cols))
                continue

            for row in _iter_loop_rows(stream, width=len(cols)):
                poly_type = row[idx_type]
                if not _is_dna_polymer_type(poly_type):
                    continue
                raw_seq = row[idx_seq]
                if raw_seq in {".", "?", ""}:
                    continue
                seq = _clean_polymer_seq(raw_seq)
                if seq:
                    dna.append(seq)
    return dna


def _open_text(path: Path) -> TextIO:
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", errors="ignore")
    return path.open("r", encoding="utf-8", errors="ignore")


def _is_dna_polymer_type(value: str) -> bool:
    t = (value or "").strip().lower()
    return "deoxyribonucleotide" in t or t in {"dna", "polydeoxyribonucleotide"}


def _clean_polymer_seq(raw: str) -> str:
    out: list[str] = []
    in_parens = 0
    for ch in raw:
        if ch == "(":
            in_parens += 1
            continue
        if ch == ")" and in_parens > 0:
            in_parens -= 1
            continue
        if in_parens:
            continue
        if ch.isspace():
            continue
        base = ch.upper()
        if base == "U":
            base = "T"
        if base in {"A", "C", "G", "T", "N"}:
            out.append(base)
        elif "A" <= base <= "Z":
            out.append("N")
    return "".join(out)


def _col_index(cols: list[str], name: str) -> Optional[int]:
    try:
        return cols.index(name)
    except ValueError:
        return None


def _skip_loop_rows(stream: "_TokenStream", *, width: int) -> None:
    for _row in _iter_loop_rows(stream, width=width):
        pass


def _iter_loop_rows(stream: "_TokenStream", *, width: int) -> Iterator[list[str]]:
    if width <= 0:
        return
    while stream.has_next():
        nxt = stream.peek()
        if nxt in {"loop_", "stop_"} or nxt.startswith("_") or nxt.startswith("data_"):
            break
        row: list[str] = []
        for _i in range(width):
            if not stream.has_next():
                return
            row.append(stream.next())
        if len(row) == width:
            yield row


def _iter_tokens(handle: TextIO) -> Iterator[str]:
    multiline: Optional[list[str]] = None
    for raw in handle:
        line = raw.rstrip("\n")
        if multiline is not None:
            if line.startswith(";"):
                yield "\n".join(multiline)
                multiline = None
                continue
            multiline.append(line)
            continue

        if not line:
            continue
        if line.lstrip().startswith("#"):
            continue
        if line.startswith(";"):
            multiline = [line[1:]]
            continue
        for tok in _split_line_tokens(line):
            yield tok
    if multiline is not None:
        yield "\n".join(multiline)


def _split_line_tokens(line: str) -> list[str]:
    tokens: list[str] = []
    i = 0
    n = len(line)
    while i < n:
        while i < n and line[i].isspace():
            i += 1
        if i >= n:
            break
        ch = line[i]
        if ch in {"'", '"'}:
            quote = ch
            i += 1
            start = i
            while i < n and line[i] != quote:
                i += 1
            tokens.append(line[start:i])
            i += 1
            continue
        start = i
        while i < n and not line[i].isspace():
            i += 1
        tok = line[start:i]
        if tok.startswith("#"):
            break
        tokens.append(tok)
    return tokens


@dataclass
class _TokenStream:
    source: Iterable[str]

    def __post_init__(self) -> None:
        self._it = iter(self.source)
        self._buffer: list[str] = []

    def has_next(self) -> bool:
        try:
            self.peek()
            return True
        except StopIteration:
            return False

    def peek(self) -> str:
        if not self._buffer:
            self._buffer.append(next(self._it))
        return self._buffer[0]

    def next(self) -> str:
        if self._buffer:
            return self._buffer.pop(0)
        return next(self._it)

