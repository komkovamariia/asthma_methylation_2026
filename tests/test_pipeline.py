"""Regression checks for reference sequences, sample mapping and launch boundaries."""

import csv
import hashlib
import json
import os
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def run_script(path, *args, env=None):
    return subprocess.run(
        [sys.executable, str(path), *map(str, args)],
        env=env,
        check=False,
        capture_output=True,
        text=True,
    )


class PipelineTests(unittest.TestCase):
    def setUp(self):
        self.workspace = tempfile.TemporaryDirectory()
        self.addCleanup(self.workspace.cleanup)
        self.temp = Path(self.workspace.name)
        self.project = self.temp / "project"
        shutil.copytree(ROOT / "second_try", self.project)
        self.raw = self.temp / "raw, with spaces"
        self.raw.mkdir()

    def add_fastq_pair(self, sample, lane=1, nested=False):
        directory = self.raw / "nested" if nested else self.raw
        directory.mkdir(exist_ok=True)
        paths = []
        for read in (1, 2):
            path = directory / f"{sample}_S1_L{lane:03d}_R{read}_001.fastq.gz"
            path.write_bytes(b"fixture")
            paths.append(path)
        return paths

    def make_sheet(self):
        return run_script(
            self.project / "make_methylseq_samplesheet.py",
            "--fastq_dir",
            self.raw,
            "--out",
            self.project / "samplesheet_methylseq.csv",
        )

    def launch_env(self):
        reference = self.temp / "reference"
        (reference / "Bisulfite_Genome").mkdir(parents=True)
        fasta = reference / "targets.fasta"
        fasta.write_text(">target\nACGT\n")
        binaries = self.temp / "bin"
        binaries.mkdir()
        nextflow = binaries / "nextflow"
        nextflow.write_text(
            "#!/usr/bin/env python3\nimport json, os, sys\n"
            "with open(os.environ['CAPTURE'], 'w') as f: json.dump(sys.argv[1:], f)\n"
            "sys.exit(int(os.environ.get('MOCK_STATUS', '0')))\n"
        )
        nextflow.chmod(0o755)
        return {
            **os.environ,
            "PROJECT_DIR": str(self.project),
            "RESULTS_DIR": str(self.temp / "results"),
            "WORK_DIR": str(self.temp / "work"),
            "REFERENCE_FASTA": str(fasta),
            "BISMARK_INDEX_DIR": str(reference),
            "CAPTURE": str(self.temp / "argv.json"),
            "PATH": str(binaries) + os.pathsep + os.environ["PATH"],
        }

    def test_reference_sequences(self):
        script = self.project / "create_tz_reference.py"
        result = run_script(script)
        self.assertEqual(result.returncode, 0, result.stderr)
        fasta = self.project / "reference/tz_amplicons/tz_targets.fasta"
        self.assertEqual(
            hashlib.sha256(fasta.read_bytes()).hexdigest(),
            "4b5c85316c3a0c2645fe95e7f48714289b7e01358740dfc3d7638e58acd34d01",
        )
        self.assertEqual(fasta.read_text().count(">"), 19)

    def test_metadata_matches_supplied_design(self):
        sheet = self.temp / "all.csv"
        with sheet.open("w", newline="") as handle:
            writer = csv.writer(handle)
            writer.writerow(("sample", "fastq_1", "fastq_2"))
            for sample in range(66, 0, -1):
                writer.writerow((f"hla{sample}", "r1", "r2"))
                writer.writerow((f"hla{sample}", "lane2r1", "lane2r2"))
        output = self.temp / "new/subdir/metadata.tsv"
        result = run_script(
            self.project / "make_tz_metadata.py",
            "--samplesheet_csv",
            sheet,
            "--out_tsv",
            output,
        )
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(
            output.read_bytes(), (ROOT / "postanalysis/sample_sheet.tsv").read_bytes()
        )

    def test_recursive_discovery_keeps_lanes(self):
        self.add_fastq_pair("hla1", lane=1)
        self.add_fastq_pair("hla1", lane=2, nested=True)
        self.add_fastq_pair("hla66")
        self.add_fastq_pair("hla67")
        self.add_fastq_pair("Undetermined")
        result = self.make_sheet()
        self.assertEqual(result.returncode, 0, result.stderr)
        with (self.project / "samplesheet_methylseq.csv").open() as handle:
            rows = list(csv.DictReader(handle))
        self.assertEqual([row["sample"] for row in rows].count("hla1"), 2)
        self.assertEqual({row["sample"] for row in rows}, {"hla1", "hla66"})
        self.assertTrue(all(Path(row["fastq_2"]).is_file() for row in rows))

    def test_missing_mate_fails(self):
        _, mate = self.add_fastq_pair("hla1")
        mate.unlink()
        result = self.make_sheet()
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("Missing or empty FASTQ", result.stderr)
        self.assertFalse((self.project / "samplesheet_methylseq.csv").exists())

    def test_empty_fastq_fails(self):
        read, _ = self.add_fastq_pair("hla1")
        read.write_bytes(b"")
        self.assertNotEqual(self.make_sheet().returncode, 0)

    def test_launch_filters_samples_and_preserves_options(self):
        for sample in ("hla1", "hla4", "hla5", "hla10"):
            self.add_fastq_pair(sample)
        self.assertEqual(self.make_sheet().returncode, 0)
        env = self.launch_env()
        result = subprocess.run(
            ["bash", str(self.project / "run_methylseq_nfcore.sh")],
            env=env,
            text=True,
            check=False,
            capture_output=True,
        )
        self.assertEqual(result.returncode, 0, result.stderr)
        with (self.project / "samplesheet_methylseq.filtered.csv").open() as handle:
            rows = list(csv.DictReader(handle))
        self.assertEqual({r["sample"] for r in rows}, {"hla1", "hla10"})
        argv = json.loads((self.temp / "argv.json").read_text())
        self.assertEqual(argv[argv.index("-r") + 1], "4.2.0")
        for flag in (
            "--skip_deduplication",
            "--non_directional",
            "--cytosine_report",
            "-resume",
        ):
            self.assertIn(flag, argv)
        self.assertIn(str(self.temp / "results"), argv)

    def test_nextflow_failure_reaches_caller(self):
        self.add_fastq_pair("hla1")
        self.assertEqual(self.make_sheet().returncode, 0)
        env = self.launch_env()
        env["MOCK_STATUS"] = "17"
        result = subprocess.run(
            ["bash", str(self.project / "run_methylseq_nfcore.sh")],
            env=env,
            check=False,
            capture_output=True,
        )
        self.assertEqual(result.returncode, 17)

    def test_empty_filtered_sheet_does_not_launch(self):
        self.add_fastq_pair("hla4")
        self.assertEqual(self.make_sheet().returncode, 0)
        env = self.launch_env()
        result = subprocess.run(
            ["bash", str(self.project / "run_methylseq_nfcore.sh")],
            env=env,
            check=False,
            capture_output=True,
        )
        self.assertNotEqual(result.returncode, 0)
        self.assertFalse((self.temp / "argv.json").exists())

    def test_optional_settings_reach_nextflow(self):
        self.add_fastq_pair("hla1")
        self.assertEqual(self.make_sheet().returncode, 0)
        env = self.launch_env()
        bed = self.temp / "targets.bed"
        bed.write_text("target\t0\t4\n")
        env.update(NON_DIRECTIONAL="0", RESUME="0", TARGET_BED=str(bed))
        result = subprocess.run(
            ["bash", str(self.project / "run_methylseq_nfcore.sh")],
            env=env,
            check=False,
            capture_output=True,
        )
        self.assertEqual(result.returncode, 0, result.stderr)
        argv = json.loads((self.temp / "argv.json").read_text())
        self.assertNotIn("--non_directional", argv)
        self.assertNotIn("-resume", argv)
        self.assertEqual(argv[argv.index("--target_regions_file") + 1], str(bed))

    def test_hg38_discovery_keeps_nested_lanes_and_partition(self):
        self.add_fastq_pair("hla1", lane=1)
        self.add_fastq_pair("hla1", lane=2, nested=True)
        self.add_fastq_pair("hla4")
        env = self.launch_env()
        reference = Path(env["BISMARK_INDEX_DIR"])
        (reference / "hg38.fa").write_text(">chr1\nACGT\n")
        diagnostic = self.temp / "diagnostic"
        shutil.copytree(ROOT / "second_try_after", diagnostic)
        sbatch = self.temp / "bin/sbatch"
        sbatch.write_text(
            "#!/usr/bin/env python3\nimport json, os, sys\n"
            "with open(os.environ['CAPTURE'], 'w') as f: json.dump(sys.argv[1:], f)\n"
        )
        sbatch.chmod(0o755)
        env.update(
            PROJECT_DIR=str(diagnostic),
            INPUT_FASTQ_DIR=str(self.raw),
            REF_DIR=str(reference),
            PARTITION="infinite",
        )
        result = subprocess.run(
            ["bash", str(diagnostic / "run_realign_hg38.sh")],
            env=env,
            check=False,
            capture_output=True,
        )
        self.assertEqual(result.returncode, 0, result.stderr)
        with (diagnostic / "samplesheet_realign_hg38.tsv").open() as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual(len(rows), 2)
        self.assertEqual({row["sample"] for row in rows}, {"hla1"})
        argv = json.loads((self.temp / "argv.json").read_text())
        self.assertIn("--partition=infinite", argv)
        self.assertIn("--array=1-2%10", argv)

    def test_hg38_array_isolates_lane_outputs(self):
        env = self.launch_env()
        reference = Path(env["BISMARK_INDEX_DIR"])
        diagnostic = self.temp / "diagnostic"
        shutil.copytree(ROOT / "second_try_after", diagnostic)
        sheet = self.temp / "realign.tsv"
        with sheet.open("w", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t")
            writer.writerow(("sample", "fastq_1", "fastq_2"))
            for lane in (1, 2):
                reads = self.add_fastq_pair("hla1", lane=lane)
                writer.writerow(("hla1", *map(str, reads)))
        bismark = self.temp / "bin/bismark"
        bismark.write_text(
            "#!/usr/bin/env python3\nimport sys\nfrom pathlib import Path\n"
            "a=sys.argv[1:]\no=Path(a[a.index('--output_dir')+1])\n"
            "b=a[a.index('--basename')+1]\n"
            "(o/(b+'_pe.bam')).write_bytes(b'bam fixture')\n"
            "(o/(b+'_PE_report.txt')).write_text('Mapping efficiency: 10%\\n')\n"
        )
        bismark.chmod(0o755)
        samtools = self.temp / "bin/samtools"
        samtools.write_text(
            "#!/usr/bin/env python3\nimport sys\nfrom pathlib import Path\n"
            "a=sys.argv[1:]\n"
            "if a[0]=='sort': Path(a[a.index('-o')+1]).write_bytes(b'sorted fixture')\n"
            "elif a[0]=='index': Path(a[-1]+'.bai').write_bytes(b'index fixture')\n"
            "else: print('stats fixture')\n"
        )
        samtools.chmod(0o755)
        env.update(
            SAMPLESHEET_TSV=str(sheet),
            ALIGN_DIR=str(self.temp / "align"),
            TMP_DIR=str(self.temp / "tmp"),
            REF_DIR=str(reference),
            SLURM_JOB_ID="test",
            SLURM_CPUS_PER_TASK="16",
        )
        for lane in (1, 2):
            env["SLURM_ARRAY_TASK_ID"] = str(lane)
            result = subprocess.run(
                ["bash", str(diagnostic / "slurm_realign_hg38.sh")],
                env=env,
                check=False,
                capture_output=True,
            )
            self.assertEqual(result.returncode, 0, result.stderr)
            output = self.temp / f"align/hla1/task_{lane}/hla1.hg38.sorted.bam"
            self.assertTrue(output.is_file())
            self.assertTrue(Path(str(output) + ".bai").is_file())
        self.assertEqual(list((self.temp / "tmp").iterdir()), [])


if __name__ == "__main__":
    unittest.main()
