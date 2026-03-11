from .core import FIP_signal


def main(argv=None):
    import matplotlib.pyplot as plt
    import sys

    args = sys.argv[1:] if argv is None else argv

    if args and args[-1].endswith(".mat"):
        signal = FIP_signal(filename=args[-1])
        signal.facecolor = "k"
        signal.axcolor = "w"

        signal.plot()
        plt.show()
