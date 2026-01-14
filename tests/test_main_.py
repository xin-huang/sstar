import pytest
from unittest.mock import patch

from sstar.__main__ import main  # import from the package, not __main__ of pytest


def test_score_subcommand_calls_runner():
    # Patch the function where it is defined: in sstar.__main__
    with patch("sstar.__main__._run_score") as mock_run:
        main([
            "score",
            "--vcf", "test.vcf",
            "--ref", "ref.txt",
            "--tgt", "tgt.txt",
            "--output", "out.txt",
        ])
        mock_run.assert_called_once()
        call_args = mock_run.call_args[0][0]
        assert call_args.vcf == "test.vcf"
        assert call_args.ref_ind == "ref.txt"
        assert call_args.tgt_ind == "tgt.txt"


def test_score_subcommand_phased_flag():
    """
    Ensure that the --phased flag on the CLI is parsed and propagated
    to _run_score (and therefore to cal_s_star).
    """
    with patch("sstar.__main__._run_score") as mock_run:
        main(
            [
                "score",
                "--vcf",
                "test.vcf",
                "--ref",
                "ref.txt",
                "--tgt",
                "tgt.txt",
                "--output",
                "out.txt",
                "--phased",
            ]
        )

        # Verify that _run_score was called exactly once
        mock_run.assert_called_once()

        # Extract the argparse Namespace passed into _run_score
        args = mock_run.call_args[0][0]

        # Confirm the CLI flag is correctly set
        assert args.phased is True



def test_quantile_subcommand_calls_runner():
    with patch("sstar.__main__._run_quantile") as mock_run:
        main([
            "quantile",
            "--model", "model.yaml",
            "--ms-dir", "/tmp/ms",
            "--N0", "10000",
            "--nsamp", "20",
            "--nreps", "10",
            "--ref-index", "1",
            "--ref-size", "10",
            "--tgt-index", "2",
            "--tgt-size", "10",
            "--mut-rate", "1e-8",
            "--rec-rate", "1e-8",
            "--seq-len", "100000",
            "--snp-num-range", "50", "200", "10",
            "--output-dir", "quant_out",
        ])
        mock_run.assert_called_once()


def test_threshold_subcommand_calls_runner():
    with patch("sstar.__main__._run_threshold") as mock_run:
        main([
            "threshold",
            "--score", "scores.txt",
            "--sim-data", "sim.txt",
            "--quantile", "0.99",
            "--output", "thr.txt",
        ])
        mock_run.assert_called_once()
        args = mock_run.call_args[0][0]
        assert args.score == "scores.txt"
        assert args.sim_data == "sim.txt"


def test_matchrate_subcommand_calls_runner():
    with patch("sstar.__main__._run_match_pct") as mock_run:
        main([
            "matchrate",
            "--vcf", "x.vcf",
            "--ref", "r.txt",
            "--tgt", "t.txt",
            "--src", "s.txt",
            "--output", "out.txt",
            "--score", "sc.txt",
        ])
        mock_run.assert_called_once()


def test_tract_subcommand_calls_runner():
    with patch("sstar.__main__._run_tract") as mock_run:
        main([
            "tract",
            "--threshold", "thr.txt",
            "--match-rate", "m1.txt",
            "--output-prefix", "op",
        ])
        mock_run.assert_called_once()
        args = mock_run.call_args[0][0]
        assert args.threshold == "thr.txt"
        assert args.match_pct == ["m1.txt"]
        assert args.output == "op"

