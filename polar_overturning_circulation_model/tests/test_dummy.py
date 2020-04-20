import pytest
from polar_overturning_circulation_model.dummy import dummy_foo


def test_dummy():
    assert dummy_foo(4) == (4 + 4)
