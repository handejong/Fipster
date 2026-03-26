from .core import FIP_signal


def _load_signal(argv=None):
    import sys

    args = sys.argv[1:] if argv is None else argv
    if not args or not (args[-1].endswith(".mat") or args[-1].endswith(".csv")):
        return None

    signal = FIP_signal(filename=args[-1])
    signal.facecolor = "k"
    signal.axcolor = "w"
    return signal


def main(argv=None):
    import matplotlib.pyplot as plt

    signal = _load_signal(argv)
    if signal is not None:
        signal.plot()
        plt.show()


def shell_main(argv=None):
    import code
    import matplotlib.pyplot as plt

    signal = _load_signal(argv)
    if signal is None:
        return

    plt.ion()
    signal.plot()

    banner = (
        "Fipster interactive shell\n"
        "Loaded recording as `signal`.\n"
        "Matplotlib interactive mode is enabled.\n"
        "Use `signal.plot()` or normal `plt` commands to open additional figures."
    )

    try:
        from IPython import start_ipython
    except ImportError:
        print(banner)
        code.interact(banner=banner, local={"signal": signal, "plt": plt})
        return

    print(banner)
    start_ipython(argv=["--matplotlib"], user_ns={"signal": signal, "plt": plt})
