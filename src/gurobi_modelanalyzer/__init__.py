__version__ = "v3.0.0"

from .common import _config

from .results_analyzer import (
    kappa_explain,
    angle_explain,
    matrix_bitmap,
    converttofractions,
)

from .solcheck import SolCheck

from .scaling import scale_model

set_env = _config.set_env
