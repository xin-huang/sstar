import pytest
import argparse
from unittest.mock import patch
from sstar.parsers.threshold_parser import add_threshold_parser, run_threshold


def test_add_threshold_parser():
    """Test that add_threshold_parser correctly adds the threshold subcommand."""
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest="subcommand")
    add_threshold_parser(subparsers)
    
    # Test with required arguments
    args = parser.parse_args([
        'threshold',
        '--sim-data', 'simulated.txt',
        '--score', 'scores.txt',
        '--quantile', '0.99',
        '--output', 'results.txt'
    ])
    
    assert args.subcommand == 'threshold'
    assert args.sim_data == 'simulated.txt'
    assert args.score == 'scores.txt'
    assert args.quantile == 0.99
    assert args.output == 'results.txt'
    assert args.recomb_rate == 1e-8  # Default value
    assert args.recomb_map is None  # Default value
    assert args.k == 8  # Default value
    
    # Test with all arguments
    args = parser.parse_args([
        'threshold',
        '--sim-data', 'simulated.txt',
        '--score', 'scores.txt',
        '--recomb-rate', '2e-8',
        '--recomb-map', 'recombination.map',
        '--quantile', '0.95',
        '--output', 'results.txt',
        '--k', '12'
    ])
    
    assert args.recomb_rate == 2e-8
    assert args.recomb_map == 'recombination.map'
    assert args.quantile == 0.95
    assert args.k == 12


@patch('sstar.calc_threshold.calc_threshold')
def test_run_threshold(mock_calc_threshold):
    """Test that run_threshold calls calc_threshold with correct arguments."""
    args = argparse.Namespace(
        sim_data='simulated.txt',
        score='scores.txt',
        recomb_rate=2e-8,
        recomb_map='recombination.map',
        quantile=0.95,
        output='results.txt',
        k=12
    )
    
    run_threshold(args)
    
    # Check that calc_threshold was called with the correct arguments
    mock_calc_threshold.assert_called_once_with(
        simulated_data='simulated.txt',
        score_file='scores.txt',
        recomb_rate=2e-8,
        recomb_map='recombination.map',
        quantile=0.95,
        output='results.txt',
        k=12
    )
