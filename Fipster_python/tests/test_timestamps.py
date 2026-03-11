def test_derive_timestamps_stores_filtered_stamps(signal_object):
    stamps = signal_object.derive_timestamps(TTL=1, min_length=0.15, max_length=0.5)

    assert len(stamps) == 1
    assert len(signal_object.timestamps["TTL 1_onset"]) == 1
    assert len(signal_object.timestamps["TTL 1_offset"]) == 1


def test_quick_peri_derives_each_ttl_independently(signal_object, monkeypatch):
    captured = []

    def fake_peri_event(stamps, window=10, from_ref=False):
        captured.append(stamps.copy())

        class DummySweepset:
            def __init__(self):
                self.settings = {}

            def make_figure(self):
                return None

        return DummySweepset()

    monkeypatch.setattr(signal_object, "peri_event", fake_peri_event)

    output = signal_object.quick_peri(TTL=[1, 2], stamps_min=0.05, stamps_max=1.0)

    assert len(output) == 2
    assert len(captured) == 2
    assert captured[0][0] != captured[1][0]
