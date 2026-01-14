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
    f1 = open("./tests/results/simulation/quantile.summary.txt", "r")
    res = f1.read()
    f1.close()

    f2 = open(pytest.exp_quantile, "r")
    exp_res = f2.read()
    f2.close()

    assert res == exp_res


# ------------------------------------------------------
# EXTRA TESTS 
# ------------------------------------------------------

def test_get_quantile_raises_if_ref_equals_tgt(tmp_path, data):
    """
    seeds=None branch + ref_index == tgt_index error.
    """
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
            tgt_index=3,  # same index -> triggers exception
            tgt_size=2,
            mut_rate=1.2e-8,
            rec_rate=0.7e-8,
            seq_len=1000,
            snp_num_range=[10, 10, 1],
            output_dir=str(outdir),
            thread=1,
        )


def test_cal_quantile_basic(tmp_path):
    """
    Direct unit test of _cal_quantile using synthetic S* data.
    """
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
    """
    Direct unit test of _summary with two fake quantile folders.
    """
    sub1 = tmp_path / "10"
    sub2 = tmp_path / "20"
    sub1.mkdir()
    sub2.mkdir()

    q1 = pd.DataFrame(
        {
            "S*_score": [0.1, 0.2],
            "SNP_num": [10, 10],
            "quantile": [0.5, 0.6],
        }
    )
    q2 = pd.DataFrame(
        {
            "S*_score": [0.3, 0.4],
            "SNP_num": [20, 20],
            "quantile": [0.7, 0.8],
        }
    )
    q1.to_csv(sub1 / "sim.quantile", sep="\t", index=False)
    q2.to_csv(sub2 / "sim.quantile", sep="\t", index=False)

    _summary(str(tmp_path), rec_rate=1e-8)

    out_file = tmp_path / "quantile.summary.txt"
    assert out_file.exists()

    out = pd.read_csv(out_file, sep="\t")

    assert "log(local_recomb_rate)" in out.columns
    assert set(out["SNP_num"]) == {10, 20}
    assert np.isclose(out["log(local_recomb_rate)"].iloc[0], -8.0)


def test_ms2vcf_roundtrip(tmp_path):
    """
    Tests ms → VCF conversion using a tiny synthetic ms file.
    """
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

    with open(ms_file, "w") as f:
        f.write(ms_content)

    _ms2vcf(str(ms_file), str(vcf_file), nsamp=4, seq_len=1000, ploidy=2)

    assert vcf_file.exists()

    with open(vcf_file, "r") as f:
        lines = [l.rstrip("\n") for l in f]

    data_lines = [l for l in lines if not l.startswith("#")]
    assert len(data_lines) == 2  # two variants

    chrom, pos, _id, ref, alt, qual, flt, info, fmt, *samples = data_lines[0].split("\t")

    assert chrom == "1"
    assert ref == "A"
    assert alt == "T"
    assert fmt == "GT"
    assert len(samples) == 2  # two diploid individuals
    assert all("|" in gt for gt in samples)
    assert set("".join(samples)) <= {"0", "1", "|"}


def test_run_ms_simulation_ref_lt_tgt_importerror(tmp_path, data, monkeypatch):
    """
    Covers:
      - ref_index < tgt_index branch (writing sim.ref.list / sim.tgt.list)
      - ImportError path for pytest_cov.embed
    """
    output_dir = tmp_path / "ms_sim"
    output_dir.mkdir()

    real_import = builtins.__import__

    def fake_import(name, *args, **kwargs):
        if name == "pytest_cov.embed":
            raise ImportError("forced for coverage")
        return real_import(name, *args, **kwargs)

    monkeypatch.setattr(builtins, "__import__", fake_import)

    # snp_num_range chosen so np.arange(...) is empty → no worker jobs
    snp_num_range = [10, 5, 1]  # start > stop

    _run_ms_simulation(
        model=pytest.model,
        ms_dir="./ext/msdir",
        N0=1000,
        nsamp=4,
        nreps=1,
        ref_index=1,     # < tgt_index → hit elif branch
        ref_size=4,
        tgt_index=2,
        tgt_size=2,
        seq_len=1000,
        snp_num_range=snp_num_range,
        output_dir=str(output_dir),
        thread=1,
        seeds=[1, 2, 3],
    )

    ref_list = output_dir / "sim.ref.list"
    tgt_list = output_dir / "sim.tgt.list"
    assert ref_list.exists()
    assert tgt_list.exists()

    ref_lines = ref_list.read_text().strip().splitlines()
    tgt_lines = tgt_list.read_text().strip().splitlines()

    # ref_size and tgt_size are halved internally (4→2, 2→1)
    assert len(ref_lines) == 2
    assert len(tgt_lines) == 1


def test_run_ms_simulation_worker_simple(tmp_path, monkeypatch):
    """
    Compact test for _run_ms_simulation_worker:
      - writes run_ms
      - calls subprocess.call, _ms2vcf, _cal_quantile
      - puts 'Finished' into out_queue
    """

    class FakeQueue:
        def __init__(self, values):
            self.values = list(values)

        def get(self):
            if self.values:
                return self.values.pop(0)
            raise StopIteration  # break worker loop

    class FakeOutQueue:
        def __init__(self):
            self.items = []

        def put(self, value):
            self.items.append(value)

    output_dir = tmp_path / "worker"
    output_dir.mkdir()
    snp_num = 10

    # Ensure subdir exists so open(.../10/run_ms, 'w') works
    (output_dir / str(snp_num)).mkdir(parents=True, exist_ok=True)

    rates = output_dir / "rates.combination"
    rates.write_text("dummy\n")

    ref_list = output_dir / "sim.ref.list"
    tgt_list = output_dir / "sim.tgt.list"
    ref_list.write_text("")
    tgt_list.write_text("")

    calls = {"subprocess": 0, "ms2vcf": 0, "cal": 0}

    def fake_subprocess_call(cmd, *args, **kwargs):
        calls["subprocess"] += 1
        return 0

    def fake_ms2vcf(ms_file, vcf_file, nsamp_arg, seq_len_arg):
        calls["ms2vcf"] += 1

    def fake_cal_quantile(score_file, quantile_file, snp_num_arg):
        calls["cal"] += 1

    monkeypatch.setattr("sstar.get_quantile.subprocess.call", fake_subprocess_call)
    monkeypatch.setattr("sstar.get_quantile._ms2vcf", fake_ms2vcf)
    monkeypatch.setattr("sstar.get_quantile._cal_quantile", fake_cal_quantile)

    ms_exec = "/usr/bin/ms"  # dummy string
    nsamp = 4
    nreps = 1
    seq_len = 1000
    ms_params = "-I 2 2 2"

    in_queue = FakeQueue([snp_num])
    out_queue = FakeOutQueue()

    with pytest.raises(StopIteration):
        _run_ms_simulation_worker(
            in_queue,
            out_queue,
            str(output_dir),
            str(rates),
            ms_exec,
            nsamp,
            nreps,
            seq_len,
            ms_params,
            str(ref_list),
            str(tgt_list),
            seeds=[1, 2, 3],  # seeds not None branch
        )

    assert out_queue.items == ["Finished"]
    assert calls["ms2vcf"] >= 1
    assert calls["cal"] >= 1

