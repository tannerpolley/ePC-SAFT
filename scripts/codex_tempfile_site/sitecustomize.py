"""Sandbox-safe tempfile behavior for local Codex build subprocesses.

Windows Codex sandbox tokens can fail to create child directories inside
``tempfile.mkdtemp()`` directories made with the stdlib's restrictive 0o700
mode. Build backends such as scikit-build-core rely on nested temporary
directories, so opt in to a less restrictive mode for local build commands.
"""

from __future__ import annotations

import errno
import os
import sys
import tempfile


def _sandbox_safe_mkdtemp(suffix=None, prefix=None, dir=None):
    prefix, suffix, dir, output_type = tempfile._sanitize_params(prefix, suffix, dir)

    names = tempfile._get_candidate_names()
    if output_type is bytes:
        names = map(os.fsencode, names)

    for _ in range(tempfile.TMP_MAX):
        name = next(names)
        path = os.path.join(dir, prefix + name + suffix)
        sys.audit("tempfile.mkdtemp", path)
        try:
            os.mkdir(path, 0o777)
        except FileExistsError:
            continue
        except PermissionError:
            if os.name == "nt" and os.path.isdir(dir) and os.access(dir, os.W_OK):
                continue
            raise
        return os.path.abspath(path)

    raise FileExistsError(errno.EEXIST, "No usable temporary directory name found")


if os.environ.get("EPCSAFT_SANDBOX_SAFE_TEMPFILE") == "1":
    tempfile.mkdtemp = _sandbox_safe_mkdtemp
