import numpy as np


def test_get_data_short_recording_detrend_does_not_raise(signal_object):
    signal_object.signal = signal_object.signal[:, :5, :]
    signal_object.raw_ref = signal_object.raw_ref[:, :5, :]
    signal_object.signal[1, :, :] = np.arange(5, dtype=float)[:, None] / signal_object.framerate
    signal_object.raw_ref[1, :, :] = np.arange(5, dtype=float)[:, None] / signal_object.framerate
    signal_object.detrend = True

    data = signal_object.get_data()

    assert data.shape == (2, 5, 2)


def test_get_ref_short_recording_detrend_does_not_raise(signal_object):
    signal_object.signal = signal_object.signal[:, :5, :]
    signal_object.raw_ref = signal_object.raw_ref[:, :5, :]
    signal_object.signal[1, :, :] = np.arange(5, dtype=float)[:, None] / signal_object.framerate
    signal_object.raw_ref[1, :, :] = np.arange(5, dtype=float)[:, None] / signal_object.framerate

    ref = signal_object.get_ref(detrend=True)

    assert ref.shape == (2, 5, 2)


def test_check_settings_removes_invalid_key_cleanly(signal_object):
    signal_object.settings["invalid"] = 123

    try:
        signal_object.check_settings()
    except Exception as exc:
        assert "invalid" in str(exc)

    assert "invalid" not in signal_object.settings


def test_find_signal_switch_supports_more_than_two_channels(signal_object):
    signal_object.signal = np.zeros((2, 10, 3), dtype=float)
    signal_object.raw_ref = np.zeros((2, 10, 3), dtype=float)
    signal_object.signal[0, 5, 0] = 10
    signal_object.raw_ref[0, 6, 0] = 8
    signal_object.nr_signals = 3

    index, size = signal_object.find_signal_switch()

    assert len(index) == 2
    assert len(size) == 2
