from __future__ import absolute_import

from . import _load_data
from . import _gdca

if not hasattr(_load_data, '__pythran__') or not hasattr(_gdca, '__pythran__'):
    import warnings
    warnings.warn(RuntimeWarning('Modules were not compiled, pyGaussDCA will be slow.'))

from .gaussdca import run
