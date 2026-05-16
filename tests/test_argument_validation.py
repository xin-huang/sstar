import argparse
import pytest

from sstar.parsers.argument_validation import (
    non_negative_int,
    non_positive_int,
    positive_int,
)


def test_argument_types():
    assert positive_int("1") == 1
    assert non_negative_int("0") == 0
    assert non_positive_int("-1") == -1


@pytest.mark.parametrize(
    ("value", "func"),
    [
        ("0", positive_int),
        ("-1", positive_int),
        ("-1", non_negative_int),
        ("1", non_positive_int),
    ],
)
def test_argument_types_raises(value, func):
    with pytest.raises(argparse.ArgumentTypeError):
        func(value)
