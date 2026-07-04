"""pdl_protocol: reusable core of the THEMIS PDL measurement protocol.

Quick API (see docs/API.md):
    from pdl_protocol import load_encounter, member_mask, shell_contrast, compute_s
    from pdl_protocol.spectral import spectral_metrics, classify_spectrum

Install:  pip install -e .        (from the repository root)
"""
__version__ = "1.2.0"

from .core import load_encounter, member_mask, shell_contrast, load_enc, shell_dn, KB
from .coords import shue_r0, shue_alpha, shue_r, jelinek_r, compute_s
from .psub import pmap, list_files
from .finder import find_candidates, analyse_encounter
from .contrast import shell_profile, population_contrast, paired_coordinate_check
from . import spectral, config

__all__ = ["load_encounter", "member_mask", "shell_contrast", "compute_s",
           "shue_r0", "shue_alpha", "shue_r", "jelinek_r",
           "pmap", "list_files", "spectral", "config", "KB", "__version__",
           "find_candidates", "analyse_encounter",
           "shell_profile", "population_contrast", "paired_coordinate_check"]
