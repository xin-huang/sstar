import pytest
import argparse
from unittest.mock import patch
from sstar.parsers.quantile_parser import add_quantile_parser, run_quantile


def test_add_quantile_parser():
    """Test that add_quantile_parser correctly adds the quantile subcommand."""
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest="subcommand")
    add_quantile_parser(subparsers)
    
    # Test with required arguments
    args = parser.parse_args([
        'quantile',
        '--model', 'demographic.yaml',
        '--ms-dir', '/path/to/ms',
        '--N0', '10000',
        '--nsamp', '100',
        '--nreps', '1000',
        '--ref-index', '1',
        '--ref-size', '50',
        '--tgt-index', '2',
        '--tgt-size', '50',
        '--mut-rate', '1e-8',
        '--rec-rate', '1e-8',
        '--seq-len', '50000',
        '--snp-num-range', '100', '500', '100',
        '--output-dir', 'output_directory'
    ])
    
    assert args.subcommand == 'quantile'
    assert args.model == 'demographic.yaml'
    assert args.ms_dir == '/path/to/ms'
    assert args.N0 == 10000
    assert args.nsamp == 100
    assert args.nreps == 1000
    assert args.ref_index == 1
    assert args.ref_size == 50
    assert args.tgt_index == 2
    assert args.tgt_size == 50
    assert args.mut_rate == 1e-8
    assert args.rec_rate == 1e-8
    assert args.seq_len == 50000
    assert args.snp_num_range == [100, 500, 100]
    assert args.output_dir == 'output_directory'
    assert args.thread == 1  # Default value
    assert args.seeds is None  # Default value
    
    # Test with optional seeds and thread
    args = parser.parse_args([
        'quantile',
        '--model', 'demographic.yaml',
        '--ms-dir', '/path/to/ms',
        '--N0', '10000',
        '--nsamp', '100',
        '--nreps', '1000',
        '--seeds', '123', '456', '789',
        '--ref-index', '1',
        '--ref-size', '50',
        '--tgt-index', '2',
        '--tgt-size', '50',
        '--mut-rate', '1e-8',
        '--rec-rate', '1e-8',
        '--seq-len', '50000',
        '--snp-num-range', '100', '500', '100',
        '--output-dir', 'output_directory',
        '--thread', '8'
    ])
    
    assert args.seeds == [123, 456, 789]
    assert args.thread == 8


def test_seeds_must_be_three():
    """Test that --seeds requires exactly 3 integer values."""
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest="subcommand")
    add_quantile_parser(subparsers)
    
    # Basic required args
    basic_args = [
        'quantile',
        '--model', 'demographic.yaml',
        '--ms-dir', '/path/to/ms',
        '--N0', '10000',
        '--nsamp', '100',
        '--nreps', '1000',
        '--ref-index', '1',
        '--ref-size', '50',
        '--tgt-index', '2',
        '--tgt-size', '50',
        '--mut-rate', '1e-8',
        '--rec-rate', '1e-8',
        '--seq-len', '50000',
        '--snp-num-range', '100', '500', '100',
        '--output-dir', 'output_directory',
    ]
    
    # Test with 1 seed (should fail)
    with pytest.raises(SystemExit):
        parser.parse_args(basic_args + ['--seeds', '123'])
    
    # Test with 2 seeds (should fail)
    with pytest.raises(SystemExit):
        parser.parse_args(basic_args + ['--seeds', '123', '456'])
    
    # Test with 4 seeds (should fail)
    with pytest.raises(SystemExit):
        parser.parse_args(basic_args + ['--seeds', '123', '456', '789', '012'])
    
    # Test with exactly 3 seeds (should pass)
    args = parser.parse_args(basic_args + ['--seeds', '123', '456', '789'])
    assert args.seeds == [123, 456, 789]


@patch('sstar.get_quantile.get_quantile')
def test_run_quantile(mock_get_quantile):
    """Test that run_quantile calls get_quantile with correct arguments."""
    # Setup mock args
    args = argparse.Namespace(
        model='demographic.yaml',
        ms_dir='/path/to/ms',
        N0=10000,
        nsamp=100,
        nreps=1000,
        ref_index=1,
        ref_size=50,
        tgt_index=2,
        tgt_size=50,
        mut_rate=1e-8,
        rec_rate=1e-8,
        seq_len=50000,
        snp_num_range=[100, 500, 100],
        output_dir='output_directory',
        thread=8,
        seeds=[123, 456, 789]
    )
    
    # Run the function
    run_quantile(args)
    
    # Check that get_quantile was called with the correct arguments
    mock_get_quantile.assert_called_once_with(
        model='demographic.yaml',
        ms_dir='/path/to/ms',
        N0=10000,
        nsamp=100,
        nreps=1000,
        ref_index=1,
        ref_size=50,
        tgt_index=2,
        tgt_size=50,
        mut_rate=1e-8,
        rec_rate=1e-8,
        seq_len=50000,
        snp_num_range=[100, 500, 100],
        output_dir='output_directory',
        thread=8,
        seeds=[123, 456, 789]
    )
