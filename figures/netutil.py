"""Shared network helper: the UniProt / PDBe / InterPro / Ensembl endpoints all time out
intermittently, so every fetch in this package goes through one retrying opener."""
import time, urllib.request

# some endpoints (geneontology.org among them) reject the default python-urllib agent
USER_AGENT = "Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) LehrPDI-figures/1.0"

def urlopen_retry(url, tries=5, timeout=40, backoff=2, headers=None):
    """Fetch `url`, retrying with exponential backoff. Returns raw bytes."""
    req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT, **(headers or {})})
    for attempt in range(tries):
        try:
            return urllib.request.urlopen(req, timeout=timeout).read()
        except Exception:
            if attempt == tries - 1:
                raise
            time.sleep(backoff ** attempt)

def urlpost_retry(url, fields, tries=5, timeout=90, backoff=2, headers=None):
    """POST form-encoded `fields`, retrying with exponential backoff. Returns raw bytes.

    Needed wherever a query carries more identifiers than fit in a URL -- STRING's
    enrichment endpoint takes the whole detected proteome at once."""
    import urllib.parse
    data = urllib.parse.urlencode(fields).encode()
    req = urllib.request.Request(url, data=data,
                                 headers={"User-Agent": USER_AGENT, **(headers or {})})
    for attempt in range(tries):
        try:
            return urllib.request.urlopen(req, timeout=timeout).read()
        except Exception:
            if attempt == tries - 1:
                raise
            time.sleep(backoff ** attempt)


_urlopen_retry = urlopen_retry          # backwards-compatible alias


# ---------------------------------------------------------------- package-relative paths
import os as _os

PKG_DIR = _os.path.dirname(_os.path.abspath(__file__))

def pkg_path(*parts) -> str:
    """Resolve a path relative to this package directory, so the modules' data files are found
    no matter what the notebook's working directory is."""
    return _os.path.join(PKG_DIR, *parts)
