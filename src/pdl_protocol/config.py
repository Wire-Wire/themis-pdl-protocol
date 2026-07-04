"""Path configuration. Set the PDL_ROOT environment variable to relocate everything;
default = the repository root. The scripts keep their original analysis-machine path
strings as provenance; P() maps them onto this configurable layout."""
import os

ROOT = os.environ.get("PDL_ROOT", os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")))  # repo root (src/pdl_protocol -> src -> repo)
DATA = os.path.join(ROOT, "data")
SUBSTRATE = os.path.join(DATA, "substrate")
OMNI_CACHE = os.path.join(DATA, "omni_cache")
OUTPUTS = os.path.join(ROOT, "outputs")
COMMITTED = os.path.join(ROOT, "committed_outputs")
PIPE = os.path.join(ROOT, "scripts", "reproduce")
EVENTS = os.path.join(DATA, "events")


def _j(base, rest):
    return os.path.join(base, *rest.split("/")) if rest else base


_MAP = [
    ("h:/0mssl/review/01_current__rebuild/substrate", SUBSTRATE),
    ("h:/0mssl/review/01_current__rebuild/runs", OUTPUTS),
    ("h:/0mssl/review/repair/option3/omni_cache", OMNI_CACHE),
    ("h:/0mssl/review/repair/option3", PIPE),  # legacy event npz now live in data/events

    ("h:/0mssl/0mssl519", os.path.join(DATA, "external_0mssl519")),
    ("h:/0mssl/review", ROOT),
]


def P(orig):
    """Map an original analysis-machine path onto the configurable repo layout."""
    s = str(orig).replace('\\', "/")
    low = s.lower()
    for pre, base in _MAP:
        if low.startswith(pre):
            return _j(base, s[len(pre):].lstrip("/"))
    return orig
