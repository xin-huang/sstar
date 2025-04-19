import pytest
import argparse
from unittest.mock import patch
from sstar.parsers.score_parser import add_score_parser, run_score


def test_add_score_parser():
    """Test that add_score_parser correctly adds the score subcommand."""
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest="subcommand")
    add_score_parser(subparsers)
    
    # Test with required arguments
    args = parser.parse_args([
        'score',
        '--vcf', 'data.vcf',
        '--ref', 'reference.txt',
        '--tgt', 'target.txt',
        '--output', 'results.txt'
    ])
    
    assert args.subcommand == 'score'
    assert args.vcf == 'data.vcf'
    assert args.ref_ind == 'reference.txt'
    assert args.tgt_ind == 'target.txt'
    assert args.output == 'results.txt'
    assert args.anc_allele is None  # Default value
    assert args.win_len == 50000  # Default value
    assert args.win_step == 10000  # Default value
    assert args.thread == 1  # Default value
    assert args.match_bonus == 5000  # Default value
    assert args.max_mismatch == 5  # Default value
    assert args.mismatch_penalty == -10000  # Default value
    
    # Test with all arguments
    args = parser.parse_args([
        'score',
        '--vcf', 'data.vcf',
        '--ref', 'reference.txt',
        '--tgt', 'target.txt',
        '--output', 'results.txt',
        '--anc-allele', 'ancestors.bed',
        '--win-len', '100000',
        '--win-step', '20000',
        '--thread', '8',
        '--match-bonus', '10000',
        '--max-mismatch', '10',
        '--mismatch-penalty', '-20000'
    ])
    
    assert args.anc_allele == 'ancestors.bed'
    assert args.win_len == 100000
    assert args.win_step == 20000
    assert args.thread == 8
    assert args.match_bonus == 10000
    assert args.max_mismatch == 10
    assert args.mismatch_penalty == -20000


@patch('sstar.calc_s_star.calc_s_star')
def test_run_score(mock_calc_s_star):
    """Test that run_score calls calc_s_star with correct arguments."""
    # Setup mock args
    args = argparse.Namespace(
        vcf='data.vcf',
        ref_ind='reference.txt',
        tgt_ind='target.txt',
        anc_allele='ancestors.bed',
        win_len=100000,
        win_step=20000,
        output='results.txt',
        thread=8,
        match_bonus=10000,
        max_mismatch=10,
        mismatch_penalty=-20000
    )
    
    # Run the function
    run_score(args)
    
    # Check that calc_s_star was called with the correct arguments
    mock_calc_s_star.assert_called_once_with(
        vcf='data.vcf',
        ref_ind_file='reference.txt',
        tgt_ind_file='target.txt',
        anc_allele_file='ancestors.bed',
        win_len=100000,
        win_step=20000,
        output='results.txt',
        thread=8,
        match_bonus=10000,
        max_mismatch=10,
        mismatch_penalty=-20000
    )
