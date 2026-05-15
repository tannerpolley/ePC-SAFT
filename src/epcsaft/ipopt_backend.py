"""Compatibility module for the optional IPOPT backend bridge."""

from __future__ import annotations

import sys

from ._optional_backends import ipopt as _ipopt

sys.modules[__name__] = _ipopt
