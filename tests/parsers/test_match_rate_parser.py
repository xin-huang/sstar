import pytest
import argparse
from unittest.mock import patch
from sstar.parsers.match_rate_parser import add_match_rate_parser, run_match_rate


def test_add_match_rate_parser():
    """Test that add_match_rate_parser correctly adds the matchrate subcommand."""
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest="subcommand")
    add_match_rate_parser(subparsers)
    
    # Test with required arguments
    args = parser.parse_args([
        'matchrate',
        '--vcf', 'data.vcf',
        '--ref', 'reference.txt',
        '--tgt', 'target.txt',
        '--src', 'source.txt',
        '--output', 'results.txt',
        '--score', 'scores.txt'
    ])
    
    assert args.subcommand == 'matchrate'
    assert args.vcf == 'data.vcf'
    assert args.ref_ind == 'reference.txt'
    assert args.tgt_ind == 'target.txt'
    assert args.src_ind == 'source.txt'
    assert args.output == 'results.txt'
    assert args.score == 'scores.txt'
    assert args.thread == 1  # Default value
    assert args.anc_allele is None  # Default value
    assert args.mapped_region_file is None  # Default value
    
    # Test with all arguments
    args = parser.parse_args([
        'matchrate',
        '--vcf', 'data.vcf',
        '--ref', 'reference.txt',
        '--tgt', 'target.txt',
        '--src', 'source.txt',
        '--output', 'results.txt',
        '--thread', '4',
        '--anc-allele', 'ancestors.bed',
        '--mapped-region', 'mapped.bed',
        '--score', 'scores.txt'
    ])
    
    assert args.thread == 4
    assert args.anc_allele == 'ancestors.bed'
    assert args.mapped_region_file == 'mapped.bed'


@patch('sstar.calc_match_rate.calc_match_pct')
def test_run_match_rate(mock_calc_match_pct):
    """Test that run_match_rate calls calc_match_pct with correct arguments."""
    # Setup mock args
    args = argparse.Namespace(
        vcf='data.vcf',
        ref_ind='reference.txt',
        tgt_ind='target.txt',
        src_ind='source.txt',
        anc_allele='ancestors.bed',
        output='results.txt',
        thread=4,
        score='scores.txt',
        mapped_region_file='mapped.bed'
    )
    
    # Run the function
    run_match_rate(args)
    
    # Check that calc_match_pct was called with the correct arguments
    mock_calc_match_pct.assert_called_once_with(
        vcf='data.vcf',
        ref_ind_file='reference.txt',
        tgt_ind_file='target.txt',
        src_ind_file='source.txt',
        anc_allele_file='ancestors.bed',
        output='results.txt',
        thread=4,
        score_file='scores.txt',
        mapped_region_file='mapped.bed'
    )
