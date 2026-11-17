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
    """Run the app. Full UI is implemented in phase 3; this is the stub."""
    import streamlit as st

    from tp53.constants import APP_TITLE

    st.set_page_config(page_title=APP_TITLE, layout="wide")
    st.title(APP_TITLE)
    st.info("The app UI is under construction — see PLAN.md phase 3.")


if __name__ == "__main__":
    main()
else:
    # Streamlit imports this module directly when serving.
    main()
