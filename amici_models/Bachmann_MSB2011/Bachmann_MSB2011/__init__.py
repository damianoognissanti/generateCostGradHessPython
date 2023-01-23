"""AMICI-generated module for model Bachmann_MSB2011"""

import amici
from pathlib import Path

# Ensure we are binary-compatible, see #556
if '0.15.0' != amici.__version__:
    raise amici.AmiciVersionError(
        f'Cannot use model `Bachmann_MSB2011` in {Path(__file__).parent}, '
        'generated with amici==0.15.0, '
        f'together with amici=={amici.__version__} '
        'which is currently installed. To use this model, install '
        'amici==0.15.0 or re-import the model with the amici '
        'version currently installed.'
    )

from Bachmann_MSB2011._Bachmann_MSB2011 import *

__version__ = '0.1.0'
