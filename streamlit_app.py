"""Entry point for Streamlit Community Cloud (AC-7).

The package lives in ``src/``, which is not on the import path when Streamlit
Cloud runs this file from the repository root, so it is added explicitly.
"""

from __future__ import annotations

import sys
from pathlib import Path

_SRC = Path(__file__).resolve().parent / "src"
if str(_SRC) not in sys.path:
    sys.path.insert(0, str(_SRC))


def main() -> None:
    """Configure the page and delegate rendering to :mod:`app.pages`."""
    import streamlit as st

    from app.pages import render_app
    from tp53.constants import APP_TITLE

    st.set_page_config(page_title=APP_TITLE, layout="wide")
    render_app()


if __name__ == "__main__":
    main()
else:
    # Streamlit imports this module directly when serving.
    main()
