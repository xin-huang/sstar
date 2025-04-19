import pytest
import argparse
from unittest.mock import patch
from sstar.parsers.tract_parser import add_tract_parser, run_tract


def test_add_tract_parser():
    """Test that add_tract_parser correctly adds the tract subcommand."""
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest="subcommand")
    add_tract_parser(subparsers)
    
    # Test with required arguments and one match-rate file
    args = parser.parse_args([
        'tract',
        '--threshold', 'threshold.txt',
        '--match-rate', 'match_src1.txt',
        '--output-prefix', 'output'
    ])
    
    assert args.subcommand == 'tract'
    assert args.threshold == 'threshold.txt'
    assert args.match_pct == ['match_src1.txt']
    assert args.output == 'output'
    assert args.diff == 0  # Default value
    
    # Test with two match-rate files and custom diff
    args = parser.parse_args([
        'tract',
        '--threshold', 'threshold.txt',
        '--match-rate', 'match_src1.txt', 'match_src2.txt',
        '--output-prefix', 'output',
        '--diff', '0.1'
    ])
    
    assert args.threshold == 'threshold.txt'
    assert args.match_pct == ['match_src1.txt', 'match_src2.txt']
    assert args.output == 'output'
    assert args.diff == 0.1


def test_match_rate_limits():
    """Test that match-rate requires between 1 and 2 arguments."""
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest="subcommand")
    add_tract_parser(subparsers)
    
    base_args = [
        'tract',
        '--threshold', 'threshold.txt',
        '--output-prefix', 'output'
    ]
    
    # No match-rate file should be allowed (match_pct=None)
    args = parser.parse_args(base_args)
    assert args.match_pct is None
    
    # Test with empty match-rate (should fail)
    with pytest.raises(SystemExit):
        parser.parse_args(base_args + ['--match-rate'])
    
    # Test with 3 match-rate files (should fail)
    with pytest.raises(argparse.ArgumentTypeError):
        parser.parse_args(base_args + ['--match-rate', 'file1.txt', 'file2.txt', 'file3.txt'])
    
    # Test with 1 file (should pass)
    args = parser.parse_args(base_args + ['--match-rate', 'file1.txt'])
    assert args.match_pct == ['file1.txt']
    
    # Test with 2 files (should pass)
    args = parser.parse_args(base_args + ['--match-rate', 'file1.txt', 'file2.txt'])
    assert args.match_pct == ['file1.txt', 'file2.txt']


@patch('sstar.get_tract.get_tract')
def test_run_tract(mock_get_tract):
    """Test that run_tract calls get_tract with correct arguments."""
    args = argparse.Namespace(
        threshold='threshold.txt',
        match_pct=['match_src1.txt', 'match_src2.txt'],
        output='output',
        diff=0.1
    )
    
    run_tract(args)
    
    # Check that get_tract was called with the correct arguments
    mock_get_tract.assert_called_once_with(
        threshold_file='threshold.txt',
        match_pct_files=['match_src1.txt', 'match_src2.txt'],
        output_prefix='output',
        diff=0.1
    )


@patch('sstar.get_tract.get_tract')
def test_run_tract_no_match_files(mock_get_tract):
    """Test that run_tract calls get_tract correctly with no match files."""
    args = argparse.Namespace(
        threshold='threshold.txt',
        match_pct=None,
        output='output',
        diff=0
    )
    
    run_tract(args)
    
    mock_get_tract.assert_called_once_with(
        threshold_file='threshold.txt',
        match_pct_files=None,
        output_prefix='output',
        diff=0
    )
