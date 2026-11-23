from tp53 import config
from tp53.clients.uniprot import fetch_query_info
from tp53.clients.ncbi import submit_blastp, poll_blast
from tp53.clients.ebi import run_clustalo
import requests

session = requests.Session()
info = fetch_query_info("P04637", session)
rid = submit_blastp(info.record.sequence, "Mus musculus", config.BlastParams(), "austin.412024020@civitas.ukrida.ac.id", session)
xml = poll_blast(rid, session, timeout_s=config.BLAST_MAX_WAIT_S, email="you@example.org")
fasta = run_clustalo(">Homo\n" + info.record.sequence + "\n", "austin.412024020@civitas.ukrida.ac.id", session, config.EBI_MAX_WAIT_S)