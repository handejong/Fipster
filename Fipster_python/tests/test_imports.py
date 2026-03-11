from fipster import FIP_signal, Sweepset, __version__
from Fipster import FIP_signal as LegacyFIPSignal


def test_public_package_exports():
    assert FIP_signal is LegacyFIPSignal
    assert Sweepset.__name__ == "Sweepset"
    assert isinstance(__version__, str)
