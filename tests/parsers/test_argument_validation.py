import argparse
import pytest
from sstar.parsers.argument_validation import required_length


def test_required_length_valid_inputs():
    """Test that required_length accepts valid numbers of arguments."""
    # Create a parser with an argument using required_length(2, 4)
    parser = argparse.ArgumentParser()
    parser.add_argument('--test', nargs='+', action=required_length(2, 4))
    
    # Test minimum number of arguments
    args = parser.parse_args(['--test', 'arg1', 'arg2'])
    assert args.test == ['arg1', 'arg2']
    
    # Test maximum number of arguments
    args = parser.parse_args(['--test', 'arg1', 'arg2', 'arg3', 'arg4'])
    assert args.test == ['arg1', 'arg2', 'arg3', 'arg4']
    
    # Test in-between number of arguments
    args = parser.parse_args(['--test', 'arg1', 'arg2', 'arg3'])
    assert args.test == ['arg1', 'arg2', 'arg3']


def test_required_length_invalid_inputs():
    """Test that required_length rejects invalid numbers of arguments."""
    parser = argparse.ArgumentParser()
    parser.add_argument('--test', nargs='+', action=required_length(2, 4))
    
    # Test too few arguments
    with pytest.raises(argparse.ArgumentTypeError, match="requires between 2 and 4 arguments"):
        parser.parse_args(['--test', 'arg1'])
    
    # Test too many arguments
    with pytest.raises(argparse.ArgumentTypeError, match="requires between 2 and 4 arguments"):
        parser.parse_args(['--test', 'arg1', 'arg2', 'arg3', 'arg4', 'arg5'])


def test_required_length_equal_bounds():
    """Test that required_length works correctly when min equals max."""
    parser = argparse.ArgumentParser()
    parser.add_argument('--test', nargs='+', action=required_length(3, 3))
    
    # Test exact number of arguments
    args = parser.parse_args(['--test', 'arg1', 'arg2', 'arg3'])
    assert args.test == ['arg1', 'arg2', 'arg3']
    
    # Test too few arguments
    with pytest.raises(argparse.ArgumentTypeError, match="requires between 3 and 3 arguments"):
        parser.parse_args(['--test', 'arg1', 'arg2'])
    
    # Test too many arguments
    with pytest.raises(argparse.ArgumentTypeError, match="requires between 3 and 3 arguments"):
        parser.parse_args(['--test', 'arg1', 'arg2', 'arg3', 'arg4'])  
