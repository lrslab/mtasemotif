from __future__ import annotations

import re


_LABEL_ORDER = ("m5C", "m6A", "m4C")
_PAREN_LABELS = {
    "4": "m4C",
    "5": "m5C",
    "6": "m6A",
}
_TEXT_PATTERNS: tuple[tuple[str, tuple[str, ...]], ...] = (
    ("m5C", ("M5C", "5MC", "5-METHYLCYTOSINE", "5 METHYLCYTOSINE")),
    ("m4C", ("M4C", "4MC", "N4-METHYLCYTOSINE", "N4 METHYLCYTOSINE")),
    ("m6A", ("M6A", "6MA", "N6-METHYLADENINE", "N6 METHYLADENINE")),
)


def methylation_labels(value: str) -> set[str]:
    text = value.strip()
    if not text:
        return set()

    upper = text.upper()
    labels: set[str] = set()
    for label, tokens in _TEXT_PATTERNS:
        if any(token in upper for token in tokens):
            labels.add(label)

    for digit in re.findall(r"\((\d)\)", upper):
        label = _PAREN_LABELS.get(digit)
        if label is not None:
            labels.add(label)

    return labels


def normalize_methylation(value: str) -> str:
    text = value.strip()
    if not text:
        return ""

    labels = methylation_labels(text)
    if not labels:
        return text
    return "/".join(label for label in _LABEL_ORDER if label in labels)
