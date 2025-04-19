import pytest
import argparse
from sstar.parsers.common_args import (
    add_common_args, add_mapped_args, add_window_args,
    add_ref_ind_args, add_tgt_ind_args, add_src_ind_args, add_score_args
)


def test_add_common_args():
    """Test that add_common_args adds expected arguments with correct defaults."""
    parser = argparse.ArgumentParser()
    add_common_args(parser)
    
    # Parse with minimal required arguments
    args = parser.parse_args(['--output', 'output.txt'])
    assert args.output == 'output.txt'
    assert args.thread == 1  # Default value
    assert args.anc_allele is None  # Default value
    
    # Parse with all arguments
    args = parser.parse_args([
        '--anc-allele', 'ancestors.bed',
        '--output', 'results.txt',
        '--thread', '8'
    ])
    assert args.anc_allele == 'ancestors.bed'
    assert args.output == 'results.txt'
    assert args.thread == 8


def test_add_mapped_args():
    """Test that add_mapped_args adds the mapped-region argument."""
    parser = argparse.ArgumentParser()
    add_mapped_args(parser)
    
    args = parser.parse_args(['--mapped-region', 'mapped_regions.bed'])
    assert args.mapped_region_file == 'mapped_regions.bed'


def test_add_window_args():
    """Test that add_window_args adds window arguments with correct defaults."""
    parser = argparse.ArgumentParser()
    add_window_args(parser)
    
    # Test defaults
    args = parser.parse_args([])
    assert args.win_len == 50000
    assert args.win_step == 10000
    
    # Test custom values
    args = parser.parse_args(['--win-len', '100000', '--win-step', '20000'])
    assert args.win_len == 100000
    assert args.win_step == 20000


def test_add_ref_ind_args():
    """Test that add_ref_ind_args adds the required reference argument."""
    parser = argparse.ArgumentParser()
    add_ref_ind_args(parser)
    
    args = parser.parse_args(['--ref', 'reference.txt'])
    assert args.ref_ind == 'reference.txt'
    
    # Should fail without required argument
    with pytest.raises(SystemExit):
        parser.parse_args([])


def test_add_tgt_ind_args():
    """Test that add_tgt_ind_args adds the required target argument."""
    parser = argparse.ArgumentParser()
    add_tgt_ind_args(parser)
    
    args = parser.parse_args(['--tgt', 'target.txt'])
    assert args.tgt_ind == 'target.txt'
    
    # Should fail without required argument
    with pytest.raises(SystemExit):
        parser.parse_args([])


def test_add_src_ind_args():
    """Test that add_src_ind_args adds the required source argument."""
    parser = argparse.ArgumentParser()
    add_src_ind_args(parser)
    
    args = parser.parse_args(['--src', 'source.txt'])
    assert args.src_ind == 'source.txt'
    
    # Should fail without required argument
    with pytest.raises(SystemExit):
        parser.parse_args([])


def test_add_score_args():
    """Test that add_score_args adds the required score argument."""
    parser = argparse.ArgumentParser()
    add_score_args(parser)
    
    args = parser.parse_args(['--score', 'scores.txt'])
    assert args.score == 'scores.txt'
    
    # Should fail without required argument
    with pytest.raises(SystemExit):
        parser.parse_args([])
