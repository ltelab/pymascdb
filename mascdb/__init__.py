# -----------------------------------------------------------------------------.
# Copyright (c) 2021-2025 MASCDB developers
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the MIT License.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#
# You should have received a copy of the MIT License
# along with this program.  If not, see <https://opensource.org/license/mit/>.
# -----------------------------------------------------------------------------.
"""MASCDB."""
import contextlib
import os
from importlib.metadata import PackageNotFoundError, version

from mascdb.pd_sns_accessor import SeabornAccessor  # noqa: F401

package_dir = os.path.dirname(os.path.realpath(__file__))

# Get version
with contextlib.suppress(PackageNotFoundError):
    __version__ = version("mascdb")
