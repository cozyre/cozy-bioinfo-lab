"""Offline tests for the HTTP clients (M2: success, timeout and error paths).

No test in this module contacts a live service: all HTTP goes through a fake
session injected into the client, and sleeps/clocks are injected where a job
is polled.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pytest
import requests

from tp53 import config
from tp53.clients import (
    RemoteHTTPError,
    RemoteServiceError,
    RemoteTimeout,
    request_with_retry,
)
from tp53.clients.ebi import run_clustalo
from tp53.clients.ncbi import (
    fetch_protein,
    poll_blast,
    submit_blastp,
    submit_blastp_staggered,
)
from tp53.clients.uniprot import fetch_query_info

DATA_DIR = Path(__file__).parent / "data"


class FakeResponse:
    """Minimal requests.Response stand-in."""

    def __init__(
        self,
        text: str = "",
        *,
        json_data: Any = None,
        status_code: int = 200,
    ) -> None:
        self.text = text
        self._json_data = json_data
        self.status_code = status_code

    def raise_for_status(self) -> None:
        if self.status_code >= 400:
            raise requests.HTTPError(
                f"HTTP {self.status_code}", response=self
            )

    def json(self) -> Any:
        if self._json_data is None:
            raise ValueError("No JSON object could be decoded")
        return self._json_data


class FakeSession:
    """Queue-based session recording every call made by a client."""

    def __init__(self, responses: list[Any]) -> None:
        self.responses = list(responses)
        self.calls: list[tuple[str, str, dict[str, Any]]] = []

    def request(self, method: str, url: str, **kwargs: Any) -> FakeResponse:
        self.calls.append((method, url, kwargs))
        item = self.responses.pop(0)
        if isinstance(item, Exception):
            raise item
        if isinstance(item, FakeResponse):
            return item
        return FakeResponse(item)


class FakeClock:
    """Monotonic clock advanced by the injected sleep function."""

    def __init__(self) -> None:
        self.now = 0.0

    def __call__(self) -> float:
        return self.now

    def sleep(self, seconds: float) -> None:
        self.now += seconds


@pytest.fixture(autouse=True)
def no_real_sleep(monkeypatch: pytest.MonkeyPatch) -> None:
    """Retry backoff must never actually wait in the offline suite."""
    monkeypatch.setattr("tp53.clients.base.time.sleep", lambda _seconds: None)


class TestRetryHelper:
    def test_success_first_attempt(self) -> None:
        session = FakeSession(["ok"])
        response = request_with_retry(session, "GET", "https://example.org")
        assert response.text == "ok"
        assert len(session.calls) == 1
        assert session.calls[0][2]["timeout"] == config.REQUEST_TIMEOUT_S

    def test_retry_on_timeout_then_success(self) -> None:
        session = FakeSession(
            [
                requests.Timeout("slow"),
                requests.Timeout("slow"),
                "ok",
            ]
        )
        response = request_with_retry(session, "GET", "https://example.org")
        assert response.text == "ok"
        assert len(session.calls) == config.HTTP_MAX_RETRIES

    def test_timeout_exhausted(self) -> None:
        session = FakeSession([requests.Timeout("slow")] * config.HTTP_MAX_RETRIES)
        with pytest.raises(RemoteTimeout, match="timed out"):
            request_with_retry(session, "GET", "https://example.org")
        assert len(session.calls) == config.HTTP_MAX_RETRIES

    def test_server_error_retried(self) -> None:
        session = FakeSession(
            [FakeResponse(status_code=503)] * config.HTTP_MAX_RETRIES
        )
        with pytest.raises(RemoteHTTPError) as excinfo:
            request_with_retry(session, "GET", "https://example.org")
        assert excinfo.value.status_code == 503
        assert len(session.calls) == config.HTTP_MAX_RETRIES

    def test_client_error_not_retried(self) -> None:
        session = FakeSession([FakeResponse(status_code=404)])
        with pytest.raises(RemoteHTTPError) as excinfo:
            request_with_retry(session, "GET", "https://example.org")
        assert excinfo.value.status_code == 404
        assert len(session.calls) == 1

    def test_connection_error(self) -> None:
        session = FakeSession(
            [requests.ConnectionError("offline")] * config.HTTP_MAX_RETRIES
        )
        with pytest.raises(RemoteServiceError, match="connection failed"):
            request_with_retry(session, "GET", "https://example.org")
        assert len(session.calls) == config.HTTP_MAX_RETRIES


class TestUniProtClient:
    def test_real_fixture_json(self) -> None:
        payload = json.loads((DATA_DIR / "uniprot_P04637.json").read_text())
        session = FakeSession([FakeResponse(json_data=payload)])
        info = fetch_query_info("P04637", session)
        assert info.record.accession == "P04637"
        assert info.record.organism == "Homo sapiens"
        assert info.record.sequence == payload["sequence"]["value"]
        assert "TP53" in info.record.description
        assert any(
            feature.type == "DNA binding"
            and (feature.start, feature.end) == (102, 292)
            for feature in info.features
        )
        assert all(
            feature.type in config.UNIPROT_FEATURE_TYPES
            for feature in info.features
        )
        assert session.calls[0][1] == f"{config.UNIPROT_BASE_URL}/P04637.json"

    def test_fuzzy_feature_boundaries_are_skipped(self) -> None:
        payload = {
            "primaryAccession": "Q12345",
            "organism": {"scientificName": "Test species"},
            "sequence": {"value": "M" * 30},
            "features": [
                {
                    "type": "Domain",
                    "description": "Fuzzy",
                    "location": {
                        "start": {"value": 1, "modifier": "UNKNOWN"},
                        "end": {"value": 10},
                    },
                },
                {
                    "type": "Domain",
                    "description": "Exact",
                    "location": {"start": {"value": 12}, "end": {"value": 20}},
                },
                {"type": "Helix", "location": {"start": {"value": 2}, "end": {"value": 5}}},
            ],
        }
        session = FakeSession([FakeResponse(json_data=payload)])
        info = fetch_query_info("Q12345", session)
        assert [(f.name, f.start, f.end) for f in info.features] == [
            ("Exact", 12, 20)
        ]

    def test_invalid_json(self) -> None:
        session = FakeSession([FakeResponse("<html>not json</html>")])
        with pytest.raises(RemoteServiceError, match="invalid JSON"):
            fetch_query_info("P04637", session)

    def test_missing_sequence(self) -> None:
        session = FakeSession([FakeResponse(json_data={"primaryAccession": "X"})])
        with pytest.raises(RemoteServiceError, match="no sequence"):
            fetch_query_info("X", session)

    def test_http_error(self) -> None:
        session = FakeSession([FakeResponse(status_code=404)])
        with pytest.raises(RemoteHTTPError) as excinfo:
            fetch_query_info("NOPE", session)
        assert excinfo.value.status_code == 404
        assert len(session.calls) == 1


class TestNCBIClient:
    def test_fetch_protein(self) -> None:
        fasta = ">NP_035770.2 cellular tumor antigen p53 [Mus musculus]\nMEEPQSDPS\n"
        session = FakeSession([FakeResponse(fasta)])
        record = fetch_protein("NP_035770.2", session, email="user@example.org")
        assert record.id == "NP_035770.2"
        assert record.organism == "Mus musculus"
        assert record.description == "cellular tumor antigen p53"
        assert record.sequence == "MEEPQSDPS"
        params = session.calls[0][2]["params"]
        assert params["id"] == "NP_035770.2"
        assert params["email"] == "user@example.org"
        assert params["tool"] == config.NCBI_TOOL_NAME
        assert params["retmode"] == "text"

    def test_fetch_protein_empty_response(self) -> None:
        session = FakeSession([FakeResponse("not a FASTA record")])
        with pytest.raises(RemoteServiceError, match="no FASTA record"):
            fetch_protein("NP_000000", session, email="user@example.org")

    def test_submit_blastp(self) -> None:
        session = FakeSession(["RID = ABC123\nRTOE = 42\n"])
        rid = submit_blastp(
            "MEEPQSDPS",
            "Mus musculus",
            config.BlastParams(),
            "user@example.org",
            session,
        )
        assert rid == "ABC123"
        method, url, kwargs = session.calls[0]
        assert (method, url) == ("POST", config.NCBI_BLAST_URL)
        data = kwargs["data"]
        assert data["CMD"] == "Put"
        assert data["PROGRAM"] == "blastp"
        assert data["DATABASE"] == "refseq_protein"
        assert data["ENTREZ_QUERY"] == '"Mus musculus"[Organism]'
        assert data["HITLIST_SIZE"] == str(config.BLAST_HITLIST_SIZE)
        assert data["email"] == "user@example.org"

    def test_submit_blastp_staggered(self) -> None:
        sleeps: list[float] = []
        session = FakeSession(["RID = RID1\n", "RID = RID2\n"])
        rids = submit_blastp_staggered(
            [("Homo sapiens", "MEEPQSDPS"), ("Mus musculus", "MEEPQSDPS")],
            config.BlastParams(),
            "user@example.org",
            session,
            sleep=sleeps.append,
        )
        assert rids == {"Homo sapiens": "RID1", "Mus musculus": "RID2"}
        assert sleeps == [config.BLAST_STAGGER_DELAY_S]
        assert [
            call[2]["data"]["ENTREZ_QUERY"] for call in session.calls
        ] == ['"Homo sapiens"[Organism]', '"Mus musculus"[Organism]']

    def test_submit_blastp_ncbi_error(self) -> None:
        html = '<div class="error msInf">Message ID#32 Error: Query contains no data</div>'
        session = FakeSession([FakeResponse(html)])
        with pytest.raises(RemoteServiceError, match="Message ID#32"):
            submit_blastp("X", "Mus musculus", config.BlastParams(), "a@b.c", session)

    def test_poll_blast_returns_xml(self) -> None:
        xml = (DATA_DIR / "blast_tp53_five_species.xml").read_text()
        clock = FakeClock()
        session = FakeSession(["Status=WAITING\n", FakeResponse(xml)])
        result = poll_blast(
            "ABC123",
            session,
            timeout_s=120,
            email="user@example.org",
            clock=clock,
            sleep=clock.sleep,
        )
        assert result == xml
        assert len(session.calls) == 2
        assert clock.now == pytest.approx(config.BLAST_POLL_INTERVAL_S)
        assert all(
            call[2]["params"]["RID"] == "ABC123"
            and call[2]["params"]["email"] == "user@example.org"
            for call in session.calls
        )

    def test_poll_blast_timeout(self) -> None:
        clock = FakeClock()
        session = FakeSession(["Status=WAITING\n"] * 5)
        with pytest.raises(RemoteTimeout, match="did not finish"):
            poll_blast(
                "ABC123",
                session,
                timeout_s=config.BLAST_POLL_INTERVAL_S + 1,
                clock=clock,
                sleep=clock.sleep,
            )
        assert len(session.calls) == 2

    def test_poll_blast_failed_job(self) -> None:
        session = FakeSession(["Status=FAILED\nMessage ID#51 Error: bad RID"])
        with pytest.raises(RemoteServiceError, match="failed"):
            poll_blast("ABC123", session, timeout_s=120, clock=FakeClock())

    def test_poll_blast_http_error(self) -> None:
        session = FakeSession([FakeResponse(status_code=500)] * config.HTTP_MAX_RETRIES)
        with pytest.raises(RemoteHTTPError) as excinfo:
            poll_blast("ABC123", session, timeout_s=120, clock=FakeClock())
        assert excinfo.value.status_code == 500


class TestEBIClient:
    def test_run_clustalo_success(self) -> None:
        fasta = ">Homo\nMEEPQSDPS\n>Mus\nMEEPQSDPT\n"
        clock = FakeClock()
        session = FakeSession(
            [
                "clustalo-R2026-0001",
                "RUNNING",
                "FINISHED",
                ">Homo\nMEEPQSDPS\n>Mus\nMEEPQSDPT\n",
            ]
        )
        aligned = run_clustalo(
            fasta,
            "user@example.org",
            session,
            timeout_s=30,
            clock=clock,
            sleep=clock.sleep,
        )
        assert aligned.startswith(">Homo")
        assert [call[1] for call in session.calls] == [
            f"{config.EBI_CLUSTALO_URL}/run",
            f"{config.EBI_CLUSTALO_URL}/status/clustalo-R2026-0001",
            f"{config.EBI_CLUSTALO_URL}/status/clustalo-R2026-0001",
            f"{config.EBI_CLUSTALO_URL}/result/clustalo-R2026-0001/out",
        ]
        assert session.calls[0][2]["data"]["email"] == "user@example.org"
        assert session.calls[0][2]["data"]["sequence"] == fasta
        assert clock.now == pytest.approx(config.EBI_POLL_INTERVAL_S)

    def test_run_clustalo_job_failed(self) -> None:
        session = FakeSession(["job-1", "FAILED"])
        with pytest.raises(RemoteServiceError, match="failed"):
            run_clustalo(">a\nAAAA\n", "a@b.c", session, 30, clock=FakeClock())

    def test_run_clustalo_timeout(self) -> None:
        clock = FakeClock()
        session = FakeSession(["job-1"] + ["RUNNING"] * 5)
        with pytest.raises(RemoteTimeout, match="did not finish"):
            run_clustalo(
                ">a\nAAAA\n",
                "a@b.c",
                session,
                2 * config.EBI_POLL_INTERVAL_S,
                clock=clock,
                sleep=clock.sleep,
            )

    def test_run_clustalo_no_job_id(self) -> None:
        session = FakeSession(["   "])
        with pytest.raises(RemoteServiceError, match="no job ID"):
            run_clustalo(">a\nAAAA\n", "a@b.c", session, 30, clock=FakeClock())

    def test_run_clustalo_submission_http_error(self) -> None:
        session = FakeSession([FakeResponse(status_code=503)] * config.HTTP_MAX_RETRIES)
        with pytest.raises(RemoteHTTPError) as excinfo:
            run_clustalo(">a\nAAAA\n", "a@b.c", session, 30, clock=FakeClock())
        assert excinfo.value.status_code == 503
