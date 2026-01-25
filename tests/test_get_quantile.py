# Apache License Version 2.0
# Copyright 2022 Xin Huang
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#    http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import builtins
import subprocess
import textwrap

import numpy as np
import pandas as pd
import pytest

from sstar.get_quantile import (
    get_quantile,
    _cal_quantile,
    _summary,
    _ms2vcf,
    _run_ms_simulation,
    _run_ms_simulation_worker,
)


# ------------------------------------------------------
# ORIGINAL TEST (UNCHANGED)
# ------------------------------------------------------


@pytest.fixture
def data():
    pytest.model = "./examples/models/BonoboGhost_4K19_no_introgression.yaml"
    pytest.exp_quantile = "./tests/results/test.quantile.exp.summary"


def test_get_quantile(data):
    get_quantile(
        model=pytest.model,
        ms_dir="./ext/msdir",
        seeds=[1, 2, 3],
        N0=1000,
        nsamp=22,
        nreps=20000,
        ref_index=4,
        ref_size=20,
        tgt_index=3,
        tgt_size=2,
        mut_rate=1.2e-8,
        rec_rate=0.7e-8,
        seq_len=40000,
        snp_num_range=[25, 30, 5],
        output_dir="./tests/results/simulation",
        thread=2,
    )
    res = open("./tests/results/simulation/quantile.summary.txt").read()
    exp_res = open(pytest.exp_quantile).read()
    assert res == exp_res


# ------------------------------------------------------
# EXTRA TESTS
# ------------------------------------------------------


def test_get_quantile_raises_if_ref_equals_tgt(tmp_path, data):
    outdir = tmp_path / "error_case"

    with pytest.raises(
        Exception,
        match="The reference population should be different from the target population",
    ):
        get_quantile(
            model=pytest.model,
            ms_dir="./ext/msdir",
            seeds=None,
            N0=1000,
            nsamp=22,
            nreps=10,
            ref_index=3,
            ref_size=20,
            tgt_index=3,
            tgt_size=2,
            mut_rate=1.2e-8,
            rec_rate=0.7e-8,
            seq_len=1000,
            snp_num_range=[10, 10, 1],
            output_dir=str(outdir),
            thread=1,
        )


def test_cal_quantile_basic(tmp_path):
    in_file = tmp_path / "scores.txt"
    out_file = tmp_path / "quantile.txt"

    df = pd.DataFrame(
        {
            "chrom": ["1", "1", "1", "1"],
            "start": [0, 0, 100, 100],
            "end": [100, 100, 200, 200],
            "S*_score": [1.0, 2.0, 3.0, 4.0],
        }
    )
    df.to_csv(in_file, sep="\t", index=False)

    _cal_quantile(str(in_file), str(out_file), snp_num=42)
    out = pd.read_csv(out_file, sep="\t")

    assert set(out.columns) == {"S*_score", "SNP_num", "quantile"}
    assert (out["SNP_num"].unique() == [42]).all()
    assert np.isclose(out["quantile"].min(), 0.5)
    assert np.isclose(out["quantile"].max(), 0.995)


def test_summary_aggregates_quantiles(tmp_path):
    sub1 = tmp_path / "10"
    sub2 = tmp_path / "20"
    sub1.mkdir()
    sub2.mkdir()

    q1 = pd.DataFrame(
        {"S*_score": [0.1, 0.2], "SNP_num": [10, 10], "quantile": [0.5, 0.6]}
    )
    q2 = pd.DataFrame(
        {"S*_score": [0.3, 0.4], "SNP_num": [20, 20], "quantile": [0.7, 0.8]}
    )
    q1.to_csv(sub1 / "sim.quantile", sep="\t", index=False)
    q2.to_csv(sub2 / "sim.quantile", sep="\t", index=False)

    _summary(str(tmp_path), rec_rate=1e-8)
    out = pd.read_csv(tmp_path / "quantile.summary.txt", sep="\t")

    assert "log(local_recomb_rate)" in out.columns
    assert set(out["SNP_num"]) == {10, 20}
    assert np.isclose(out["log(local_recomb_rate)"].iloc[0], -8.0)


def test_ms2vcf_roundtrip(tmp_path):
    ms_file = tmp_path / "sim.ms"
    vcf_file = tmp_path / "sim.vcf"

    ms_content = textwrap.dedent(
        """\
        ms 4 1 -t 5 -r 5 1000

        //
        positions: 0.1 0.5
        0100
        0011
        1100
        0011
        """
    )
    ms_file.write_text(ms_content)

    _ms2vcf(str(ms_file), str(vcf_file), nsamp=4, seq_len=1000, ploidy=2)
    lines = [l for l in vcf_file.read_text().splitlines() if not l.startswith("#")]

    assert len(lines) == 2
    fields = lines[0].split("\t")
    assert fields[0] == "1"
    assert fields[3] == "A"
    assert fields[4] == "T"
    assert all("|" in gt for gt in fields[9:])


def test_run_ms_simulation_ref_lt_tgt_importerror(tmp_path, data, monkeypatch):
    output_dir = tmp_path / "ms_sim"
    output_dir.mkdir()

    real_import = builtins.__import__

    def fake_import(name, *args, **kwargs):
        if name == "pytest_cov.embed":
            raise ImportError
        return real_import(name, *args, **kwargs)

    monkeypatch.setattr(builtins, "__import__", fake_import)

    _run_ms_simulation(
        model=pytest.model,
        ms_dir="./ext/msdir",
        N0=1000,
        nsamp=4,
        nreps=1,
        ref_index=1,
        ref_size=4,
        tgt_index=2,
        tgt_size=2,
        seq_len=1000,
        snp_num_range=[10, 5, 1],
        output_dir=str(output_dir),
        thread=1,
        seeds=[1, 2, 3],
    )

    assert (output_dir / "sim.ref.list").exists()
    assert (output_dir / "sim.tgt.list").exists()


# ------------------------------------------------------
# UPDATED TEST — FULL PHASED / UNPHASED COVERAGE
# ------------------------------------------------------


def test_run_ms_simulation_worker_simple(tmp_path, monkeypatch):
    """
    Fully covers phased and unphased logic in _run_ms_simulation_worker.
    """

    class FakeQueue:
        def __init__(self, values):
            self.values = list(values)

        def get(self):
            if self.values:
                return self.values.pop(0)
            raise StopIteration

    class FakeOutQueue:
        def __init__(self):
            self.items = []

        def put(self, value):
            self.items.append(value)

    def run_case(is_phased):
        output_dir = tmp_path / ("phased" if is_phased else "unphased")
        output_dir.mkdir()
        (output_dir / "10").mkdir()

        (output_dir / "rates.combination").write_text("x\n")
        (output_dir / "sim.ref.list").write_text("")
        (output_dir / "sim.tgt.list").write_text("")

        captured = []

        # Intercept the "bash run_ms.sh" call
        def fake_check_call(cmd, *a, **k):
            captured.append(cmd)
            return 0

        # Intercept the "sstar score ..." call (now uses subprocess.run)
        def fake_run(cmd, *a, **k):
            captured.append(cmd)
            return subprocess.CompletedProcess(cmd, 0)

        monkeypatch.setattr("sstar.get_quantile.subprocess.check_call", fake_check_call)
        monkeypatch.setattr("sstar.get_quantile.subprocess.run", fake_run)
        monkeypatch.setattr("sstar.get_quantile._ms2vcf", lambda *a, **k: None)
        monkeypatch.setattr("sstar.get_quantile._cal_quantile", lambda *a, **k: None)

        in_q = FakeQueue([10])
        out_q = FakeOutQueue()

        with pytest.raises(StopIteration):
            _run_ms_simulation_worker(
                in_q,
                out_q,
                str(output_dir),
                str(output_dir / "rates.combination"),
                "/usr/bin/ms",
                4,
                1,
                1000,
                "-I 2 2 2",
                str(output_dir / "sim.ref.list"),
                str(output_dir / "sim.tgt.list"),
                seeds=[1, 2, 3],
                is_phased=is_phased,
            )

        score_cmd = [c for c in captured if isinstance(c, list) and c[0] == "sstar"][0]
        assert ("--phased" in score_cmd) is is_phased
        assert out_q.items == ["Finished"]

    run_case(False)
    run_case(True)

