import sys
from .version_string import package_versions

__version__ = '0.0.1'

assert sys.version_info >= (3, 7), "Python 3.7+ required."