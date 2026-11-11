"""Download global bilateral MONTHLY trade flows from UN Comtrade,
HS chapter 27 at 4-digit level (2701-2716), 2016-2024.

Monthly goods data (typeCode C, freq M, classification HS), all reporters x all
partners, cmdCode=2701..2716, flows M and X.

Sharding: one CSV per reporter-year (12 monthly periods per API call) in
raw/, e.g. hs27m_156_2023.csv. If a year would hit the 250k-record cap it is
re-fetched quarter by quarter, and a capped quarter month by month.
Existing non-empty files are skipped -> fully resumable.

Descriptions are NOT repeated in every row (includeDesc=false); code lookup
tables are saved once to reference/.

API key: plain-text single line in .comtrade_key (never printed/logged).

Usage:
    python download_hs27.py                 # everything
    python download_hs27.py --years 2016,2017
    python download_hs27.py --only 156,842  # reporter codes
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import pandas as pd
import requests

BASE_DIR = Path(__file__).resolve().parent
RAW_DIR = BASE_DIR / "raw"
REF_DIR = BASE_DIR / "reference"
KEY_FILE = BASE_DIR / ".comtrade_key"

YEARS = list(range(2016, 2025))
HS4_CODES = ",".join(f"27{i:02d}" for i in range(1, 17))  # 2701..2716
QUARTERS = {1: (1, 2, 3), 2: (4, 5, 6), 3: (7, 8, 9), 4: (10, 11, 12)}
API_URL = "https://comtradeapi.un.org/data/v1/get/C/M/HS"
REF_BASE = "https://comtradeapi.un.org/files/v1/app/reference"
MAX_RECORDS = 250_000
CALL_GAP_SECONDS = 1.2
RETRIES = 4


def load_key() -> str:
    if not KEY_FILE.exists():
        sys.exit(f"API key file not found: {KEY_FILE}")
    key = KEY_FILE.read_text().strip()
    if not key:
        sys.exit(f"API key file is empty: {KEY_FILE}")
    return key


def save_reference_tables() -> pd.DataFrame:
    """Download code lookup tables once; return reporter list."""
    REF_DIR.mkdir(parents=True, exist_ok=True)
    tables = {
        "reporters": f"{REF_BASE}/Reporters.json",
        "partners": f"{REF_BASE}/partnerAreas.json",
        "hs_codes": f"{REF_BASE}/HS.json",
    }
    for name, url in tables.items():
        out = REF_DIR / f"{name}.json"
        if not out.exists():
            r = requests.get(url, timeout=60)
            r.raise_for_status()
            out.write_text(json.dumps(r.json()["results"]))
            print(f"reference table saved: {out.name}", flush=True)
    df = pd.DataFrame(json.loads((REF_DIR / "reporters.json").read_text()))
    df = df[df["isGroup"] == False]  # noqa: E712  (drop aggregates like EU)
    return df[["id", "text"]].astype({"id": int})


def api_call(key: str, reporter: int, periods: str) -> pd.DataFrame | None:
    """One API call; returns None after exhausting retries.

    Raises RuntimeError("cap") if the record cap was hit (caller must split)."""
    params = {
        "reporterCode": reporter,
        "period": periods,
        "cmdCode": HS4_CODES,
        "flowCode": "M,X",
        "maxRecords": MAX_RECORDS,
        "format": "JSON",
        "breakdownMode": "classic",
        "includeDesc": "false",
    }
    headers = {"Ocp-Apim-Subscription-Key": key}
    for attempt in range(1, RETRIES + 1):
        try:
            r = requests.get(API_URL, params=params, headers=headers, timeout=180)
            if r.status_code == 429:
                wait = 30 * attempt
                print(f"    429 rate-limited, sleeping {wait}s", flush=True)
                time.sleep(wait)
                continue
            r.raise_for_status()
            payload = r.json()
            if payload.get("count", 0) >= MAX_RECORDS:
                raise RuntimeError("cap")
            return pd.DataFrame(payload.get("data", []))
        except RuntimeError:
            raise
        except (requests.RequestException, ValueError) as exc:
            wait = 15 * attempt
            print(f"    attempt {attempt}/{RETRIES} failed ({exc}); "
                  f"retry in {wait}s", flush=True)
            time.sleep(wait)
    return None


def fetch_months(key: str, reporter: int, year: int,
                 months: tuple[int, ...]) -> pd.DataFrame | None:
    """Fetch a batch of months in one call; on cap, split month by month."""
    periods = ",".join(f"{year}{m:02d}" for m in months)
    try:
        return api_call(key, reporter, periods)
    except RuntimeError:
        if len(months) == 1:
            print(f"    WARNING: single month {year}-{months[0]:02d} still "
                  f"capped (reporter {reporter}) — data truncated", flush=True)
            return None
        print(f"    cap hit on {periods}, splitting", flush=True)
        parts = []
        for m in months:
            part = fetch_months(key, reporter, year, (m,))
            if part is None:
                return None
            parts.append(part)
            time.sleep(CALL_GAP_SECONDS)
        return pd.concat(parts, ignore_index=True)


def fetch_year(key: str, reporter: int, year: int) -> pd.DataFrame | None:
    """Whole year in one call; on cap, fall back to quarters (then months)."""
    all_months = tuple(range(1, 13))
    try:
        return api_call(key, reporter,
                        ",".join(f"{year}{m:02d}" for m in all_months))
    except RuntimeError:
        print(f"    cap hit, splitting {year} into quarters", flush=True)
        parts = []
        for q in (1, 2, 3, 4):
            part = fetch_months(key, reporter, year, QUARTERS[q])
            if part is None:
                return None
            parts.append(part)
            time.sleep(CALL_GAP_SECONDS)
        return pd.concat(parts, ignore_index=True)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--years", help="comma list, default 2016-2024")
    ap.add_argument("--only", help="comma list of reporter codes (shard/debug)")
    args = ap.parse_args()

    years = ([int(y) for y in args.years.split(",")] if args.years else YEARS)
    key = load_key()
    RAW_DIR.mkdir(parents=True, exist_ok=True)

    reporters = save_reference_tables()
    if args.only:
        keep = {int(c) for c in args.only.split(",")}
        reporters = reporters[reporters["id"].isin(keep)]
    total = len(reporters) * len(years)
    print(f"{len(reporters)} reporters x {len(years)} years = {total} tasks",
          flush=True)

    done = skipped = empty = failed = 0
    for year in years:
        for _, rep in reporters.iterrows():
            code, name = rep["id"], rep["text"]
            out = RAW_DIR / f"hs27m_{code}_{year}.csv"
            if out.exists() and out.stat().st_size > 0:
                skipped += 1
                continue
            df = fetch_year(key, code, year)
            if df is None:
                failed += 1
                print(f"  FAILED  {name} ({code}) {year}", flush=True)
                continue
            df.to_csv(out, index=False)
            done += 1
            if df.empty:
                empty += 1
            print(f"  [{done + skipped}/{total}] {name} ({code}) {year}: "
                  f"{'empty' if df.empty else f'{len(df)} rows'}", flush=True)
            time.sleep(CALL_GAP_SECONDS)

    print(f"\nfinished: {done} downloaded ({empty} empty), "
          f"{skipped} skipped (already present), {failed} failed", flush=True)
    if failed:
        sys.exit(1)


if __name__ == "__main__":
    main()
