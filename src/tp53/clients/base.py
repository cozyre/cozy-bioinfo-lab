"""Shared HTTP plumbing for all remote clients (§5, §8).

Every network client takes an injected requests-like session so tests can
replace it with a mock and never call a live service. All calls use a timeout
and exponential-backoff retries; failures become exceptions from this module
so the UI can show one consistent error message per stage.
"""

from __future__ import annotations

import time
from collections.abc import Callable
from dataclasses import dataclass
from typing import Any, Literal

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
    retry: bool = True,
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
        retry: whether transient failures may be retried. Set to ``False``
            for non-idempotent submissions whose retry could create a second
            remote job when the first response was lost after acceptance.
        **kwargs: forwarded to ``session.request`` (params, data, ...).

    Returns:
        The successful response.

    Raises:
        RemoteTimeout: every attempt timed out.
        RemoteHTTPError: non-retryable HTTP error, or retries exhausted.
        RemoteServiceError: connection failure, or retries exhausted.
    """
    last_error: RemoteServiceError | None = None
    for attempt in range(max_retries if retry else 1):
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

        if retry and attempt + 1 < max_retries:
            sleep(backoff_factor**attempt)
    assert last_error is not None  # max_retries >= 1, every failure sets it
    raise last_error


@dataclass(frozen=True)
class PollOutcome:
    """Result of one check of an asynchronous remote job (§5, §9.1).

    Services classify their own statuses; the shared polling loop only needs
    to know whether the job is still pending, finished, or failed.
    """

    state: Literal["pending", "done", "failed"]
    value: str = ""
    error: str = ""

    @classmethod
    def pending(cls) -> PollOutcome:
        """The job has not finished yet."""
        return cls("pending")

    @classmethod
    def complete(cls, value: str) -> PollOutcome:
        """The job finished; ``value`` is the fetched result text."""
        return cls("done", value=value)

    @classmethod
    def failed(cls, message: str) -> PollOutcome:
        """The job failed; ``message`` is the complete error message."""
        return cls("failed", error=message)


def poll_until(
    check: Callable[[], PollOutcome],
    *,
    label: str,
    timeout_s: int,
    interval_s: int,
    sleep: Callable[[float], None] = time.sleep,
    clock: Callable[[], float] = time.monotonic,
) -> str:
    """Call ``check`` until the job finishes, fails, or the deadline passes.

    The deadline is checked before every request and the wait is clamped to
    the remaining time, so a job that finishes just in time is not preceded by
    a full poll interval of idle waiting. ``sleep`` and ``clock`` are injected
    so offline tests never wait in real time.

    Args:
        check: zero-argument function performing one poll and interpreting it.
        label: job description used in the timeout message.
        timeout_s: overall wall-clock deadline in seconds.
        interval_s: minimum seconds between two checks.
        sleep: sleep function, injectable for offline tests.
        clock: monotonic clock, injectable for offline tests.

    Returns:
        Result text from :meth:`PollOutcome.complete`.

    Raises:
        RemoteServiceError: ``check`` reported a failed job.
        RemoteTimeout: the deadline passed before the job completed.
    """
    start = clock()
    while clock() - start < timeout_s:
        outcome = check()
        if outcome.state == "done":
            return outcome.value
        if outcome.state == "failed":
            raise RemoteServiceError(outcome.error)
        remaining = timeout_s - (clock() - start)
        if remaining <= 0:
            break
        sleep(min(interval_s, remaining))
    raise RemoteTimeout(f"{label} did not finish within {timeout_s} s.")
