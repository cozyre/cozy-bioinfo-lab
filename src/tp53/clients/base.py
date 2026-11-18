"""Shared HTTP plumbing for all remote clients (§5, §8).

Every network client takes an injected requests-like session so tests can
replace it with a mock and never call a live service. All calls use a timeout
and exponential-backoff retries; failures become exceptions from this module
so the UI can show one consistent error message per stage.
"""

from __future__ import annotations

import time
from collections.abc import Callable
from typing import Any

import requests

from tp53 import config


class RemoteServiceError(RuntimeError):
    """A remote service call failed after exhausting its retry budget."""


class RemoteTimeout(RemoteServiceError):
    """A remote service call timed out (per request or overall job deadline)."""


class RemoteHTTPError(RemoteServiceError):
    """A remote service returned an HTTP error status.

    Attributes:
        status_code: HTTP status of the last failing response, if known.
    """

    def __init__(self, message: str, status_code: int | None = None) -> None:
        super().__init__(message)
        self.status_code = status_code


def request_with_retry(
    session: Any,
    method: str,
    url: str,
    *,
    timeout: int = config.REQUEST_TIMEOUT_S,
    max_retries: int = config.HTTP_MAX_RETRIES,
    backoff_factor: float = config.HTTP_BACKOFF_FACTOR_S,
    sleep: Callable[[float], None] = time.sleep,
    **kwargs: Any,
) -> requests.Response:
    """Send one HTTP request, retrying transient failures with backoff.

    Retries timeouts, connection errors, HTTP 429 and 5xx responses. Other
    HTTP errors (e.g. a bad accession) fail immediately because retrying them
    cannot change the outcome. Tests inject ``sleep`` to avoid real delays.

    Args:
        session: injected requests-like session exposing ``request()``.
        method: HTTP method, e.g. ``"GET"`` or ``"POST"``.
        url: absolute URL.
        timeout: per-request timeout in seconds.
        max_retries: total number of attempts.
        backoff_factor: base for exponential backoff (factor ** attempt).
        sleep: sleep function, injectable for offline tests.
        **kwargs: forwarded to ``session.request`` (params, data, ...).

    Returns:
        The successful response.

    Raises:
        RemoteTimeout: every attempt timed out.
        RemoteHTTPError: non-retryable HTTP error, or retries exhausted.
        RemoteServiceError: connection failure, or retries exhausted.
    """
    last_error: RemoteServiceError | None = None
    for attempt in range(max_retries):
        try:
            response = session.request(method, url, timeout=timeout, **kwargs)
            response.raise_for_status()
            return response
        except requests.Timeout as exc:
            last_error = RemoteTimeout(f"{method} {url} timed out: {exc}")
        except requests.ConnectionError as exc:
            last_error = RemoteServiceError(f"{method} {url} connection failed: {exc}")
        except requests.HTTPError as exc:
            response = exc.response
            status = getattr(response, "status_code", None)
            message = f"{method} {url} returned HTTP {status}: {exc}"
            if status is not None and status not in config.HTTP_RETRYABLE_STATUSES:
                raise RemoteHTTPError(message, status) from exc
            last_error = RemoteHTTPError(message, status)
        except requests.RequestException as exc:
            last_error = RemoteServiceError(f"{method} {url} failed: {exc}")

        if attempt + 1 < max_retries:
            sleep(backoff_factor**attempt)
    assert last_error is not None  # max_retries >= 1, every failure sets it
    raise last_error
