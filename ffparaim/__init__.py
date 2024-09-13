# __init__.py

try:
    from ._version import __version__, __version_tuple__
except ImportError:
    __version__ = "0.0.0.post0"
    __version_tuple__ = (0, 0, 0, "a-dev")

from .ffparaim import FFparAIM
