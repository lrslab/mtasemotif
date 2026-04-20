from __future__ import annotations

import json
import csv
import tempfile
import unittest
from pathlib import Path

from mtase_motif.structure_motif.foldseek import FoldseekHit
from mtase_motif.structure_motif.template_panel import infer_template_panel_motif
from mtase_motif.structure_motif.training import build_structure_label_db
from mtase_motif.structure_motif.validation import validate_structure_predictions


MMCIF_WITH_DNA = """data_test
loop_
_entity_poly.entity_id
_entity_poly.type
_entity_poly.pdbx_seq_one_letter_code_can
1 polydeoxyribonucleotide GATC
#
"""

MMCIF_NO_PANEL_DNA = """data_test
loop_
_entity_poly.entity_id
_entity_poly.type
_entity_poly.pdbx_seq_one_letter_code_can
1 polydeoxyribonucleotide NNNN
#
"""


class StructureMotifTests(unittest.TestCase):
    def test_template_panel_infers_motif_from_mmcif(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            (root / "1abc.cif").write_text(MMCIF_WITH_DNA)
            panel = infer_template_panel_motif(
                [FoldseekHit(target="1abc", evalue=1e-30, bits=120.0, alntmscore=0.74)],
                mmcif_mirror_dir=root,
                mmcif_cache_dir=root / "cache",
                allow_download=False,
                top_hits=1,
                min_alntmscore=0.5,
                kmin=4,
                kmax=4,
                panel_size=1,
            )
            self.assertIsNotNone(panel)
            assert panel is not None
            self.assertEqual(panel.motif_iupac, "GATC")
            self.assertEqual(panel.templates_used, 1)
            self.assertEqual(panel.confidence, "low")

    def test_template_panel_uses_contributing_hits_for_provenance(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            (root / "1abc.cif").write_text("# no DNA here\n")
            (root / "2def.cif").write_text(MMCIF_WITH_DNA)
            (root / "3ghi.cif").write_text(MMCIF_WITH_DNA)
            panel = infer_template_panel_motif(
                [
                    FoldseekHit(target="1abc", evalue=1e-40, bits=300.0, alntmscore=0.99),
                    FoldseekHit(target="2def", evalue=1e-20, bits=120.0, alntmscore=0.59),
                    FoldseekHit(target="3ghi", evalue=1e-18, bits=118.0, alntmscore=0.58),
                ],
                mmcif_mirror_dir=root,
                mmcif_cache_dir=root / "cache",
                allow_download=False,
                top_hits=3,
                min_alntmscore=0.5,
                kmin=4,
                kmax=4,
                panel_size=1,
            )
            self.assertIsNotNone(panel)
            assert panel is not None
            self.assertEqual(panel.best_target, "2def")
            self.assertEqual(panel.contributing_targets, ("2def", "3ghi"))
            self.assertEqual(panel.templates_used, 2)
            self.assertEqual(panel.confidence, "low")

    def test_build_structure_label_db_writes_labels_and_summary(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            mmcif = root / "1abc.cif"
            mmcif.write_text(MMCIF_WITH_DNA)
            manifest = root / "manifest.tsv"
            manifest.write_text(
                "\t".join(
                    ["target", "motif_iupac", "methylation", "family_id", "split", "structure_path"]
                )
                + "\n"
                + "\t".join(["1abc_model", "GATC", "m6A", "fam1", "train", str(mmcif)])
                + "\n"
                + "\t".join(["2def_model", "CCWGG", "m5C", "fam2", "test", ""])
                + "\n"
            )
            out_dir = root / "training"
            summary = build_structure_label_db(
                manifest,
                out_dir,
                allow_download=False,
                kmin=4,
                kmax=4,
                panel_size=1,
            )
            self.assertEqual(summary.total_rows, 2)
            self.assertEqual(summary.label_transfer_ready, 2)
            self.assertEqual(summary.template_panel_ready, 1)
            self.assertEqual(summary.missing_structure, 1)

            labels = (out_dir / "labels.tsv").read_text().splitlines()
            self.assertEqual(labels[0], "target\tmotif_iupac\tmethylation")
            self.assertIn("1abc_model\tGATC\tm6A", labels[1:])

            examples = (out_dir / "training_examples.tsv").read_text()
            self.assertIn("self_panel_iupac", examples)
            self.assertIn("GATC", examples)

            training_summary = json.loads((out_dir / "training_summary.json").read_text())
            self.assertEqual(training_summary["unique_families"], 2)
            self.assertEqual(training_summary["split_counts"]["train"], 1)
            self.assertEqual(training_summary["split_counts"]["test"], 1)

    def test_build_structure_label_db_requires_successful_panel_construction(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            mmcif = root / "1abc.cif"
            mmcif.write_text(MMCIF_NO_PANEL_DNA)
            manifest = root / "manifest.tsv"
            manifest.write_text(
                "\t".join(
                    ["target", "motif_iupac", "methylation", "family_id", "split", "structure_path"]
                )
                + "\n"
                + "\t".join(["1abc_model", "GATC", "m6A", "fam1", "train", str(mmcif)])
                + "\n"
            )
            out_dir = root / "training"
            summary = build_structure_label_db(
                manifest,
                out_dir,
                allow_download=False,
                kmin=4,
                kmax=4,
                panel_size=1,
            )

            self.assertEqual(summary.total_rows, 1)
            self.assertEqual(summary.label_transfer_ready, 1)
            self.assertEqual(summary.template_panel_ready, 0)
            with (out_dir / "training_examples.tsv").open(newline="", encoding="utf-8") as handle:
                row = next(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(row["usable_for_template_panel"], "False")

    def test_validate_structure_predictions_reports_match_classes(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            gold = root / "gold.tsv"
            pred = root / "pred.tsv"
            gold.write_text(
                "\t".join(["candidate_id", "motif_iupac", "methylation", "split"])
                + "\n"
                + "\t".join(["cand_exact", "GATC", "m6A", "validation"])
                + "\n"
                + "\t".join(["cand_rc", "AGC", "m5C", "validation"])
                + "\n"
                + "\t".join(["cand_compatible", "AAN", "m6A", "validation"])
                + "\n"
                + "\t".join(["cand_missing", "CCWGG", "m5C", "validation"])
                + "\n"
            )
            pred.write_text(
                "\t".join(["candidate_id", "motif_iupac", "methylation", "method", "confidence"])
                + "\n"
                + "\t".join(["cand_exact", "GATC", "m6A", "foldseek_structure", "high"])
                + "\n"
                + "\t".join(["cand_rc", "GCT", "m5C", "foldseek_structure", "medium"])
                + "\n"
                + "\t".join(["cand_compatible", "AAA", "m6A", "foldseek_template_panel", "low"])
                + "\n"
            )

            out_dir = root / "validation"
            summary = validate_structure_predictions(
                gold,
                pred,
                out_dir,
                split_filter={"validation"},
            )

            self.assertEqual(summary.total_gold, 4)
            self.assertEqual(summary.resolved_predictions, 3)
            self.assertEqual(summary.exact_matches, 1)
            self.assertEqual(summary.reverse_complement_matches, 1)
            self.assertEqual(summary.compatible_matches, 3)
            self.assertEqual(summary.methylation_matches, 3)
            self.assertEqual(summary.unresolved_predictions, 1)
            self.assertGreater(summary.mean_positional_overlap, 0.0)
            self.assertTrue((out_dir / "validation_results.tsv").exists())
            self.assertTrue((out_dir / "validation_summary.json").exists())


if __name__ == "__main__":
    unittest.main()
