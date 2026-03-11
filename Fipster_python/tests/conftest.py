import numpy as np
import pandas as pd
import pytest

from fipster import FIP_signal


@pytest.fixture
def signal_object():
    signal = FIP_signal()
    signal.hasdata = True
    signal.framerate = 10
    signal.signal = np.zeros((2, 100, 2), dtype=float)
    signal.raw_ref = np.zeros((2, 100, 2), dtype=float)
    timeline = np.arange(100, dtype=float) / signal.framerate
    trace = np.linspace(1.0, 2.0, 100)
    signal.signal[0, :, 0] = trace * 1.5 + 0.1
    signal.signal[0, :, 1] = trace * 1.2 + 0.2
    signal.raw_ref[0, :, 0] = trace
    signal.raw_ref[0, :, 1] = trace * 0.8
    signal.signal[1, :, :] = timeline[:, None]
    signal.raw_ref[1, :, :] = timeline[:, None]
    signal.nr_signals = 2
    signal.labels = ["green", "red"]
    signal.filename = "synthetic.mat"
    signal.external_signal = [False, False]
    signal.smooth = [False, False]
    signal.logAI = pd.DataFrame(
        {
            "Time": np.arange(0.005, 1.005, 0.005),
            "TTL 1": np.where((np.arange(200) >= 20) & (np.arange(200) < 60), 5.0, 0.0),
            "TTL 2": np.where((np.arange(200) >= 100) & (np.arange(200) < 180), 5.0, 0.0),
        }
    )
    return signal
