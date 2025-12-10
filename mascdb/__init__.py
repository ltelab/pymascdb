import contextlib
import os
from importlib.metadata import PackageNotFoundError, version

from mascdb.pd_sns_accessor import SeabornAccessor  # noqa: F401

package_dir = os.path.dirname(os.path.realpath(__file__))

# Get version
with contextlib.suppress(PackageNotFoundError):
    __version__ = version("disdrodb")
