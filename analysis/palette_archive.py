"""Fetch the DOI-pinned Steen localization archive and inventory every extracted file.

Run on a POSIX data volume. Resume downloads in the same directory.
"""
import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import fcntl
import hashlib
import json
import os
from pathlib import Path
import time
import zlib
from urllib.error import HTTPError, URLError
from urllib.request import Request, urlopen

URL = "https://zenodo.org/api/records/10960858/files/dye_analysis.rar/content"
SIZE = 38_156_495_377
MD5 = "4d46f75d5cd23313925f044298c89ddf"
BLOCK = 16 * 1024**2


def digest(path, algorithm):
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, algorithm).hexdigest()


def fetch(directory, workers):
    directory.mkdir(parents=True, exist_ok=True)
    dest = directory / "dye_analysis.rar"
    part = dest.with_suffix(".rar.part")
    journal = dest.with_suffix(".ranges")
    with (directory / ".download.lock").open("a") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        if dest.exists():
            if dest.stat().st_size != SIZE or digest(dest, "md5") != MD5:
                raise ValueError("Existing archive differs from the DOI-pinned source")
            return
        total = (SIZE + BLOCK - 1) // BLOCK
        done = set()
        if part.exists() and part.stat().st_size == SIZE and journal.exists():
            done = {int(i) for i in journal.read_text().splitlines() if i.isdigit() and int(i) < total}
        else:
            journal.write_text("")
        with part.open("r+b" if part.exists() else "w+b") as data, journal.open("a", buffering=1) as ledger:
            data.truncate(SIZE)
            def transfer(i):
                start, end = i * BLOCK, min((i + 1) * BLOCK, SIZE) - 1
                req = Request(URL, headers={"Range": f"bytes={start}-{end}",
                                           "User-Agent": "Blink-Atlas-reproducer/2.0"})
                for attempt in range(10):
                    try:
                        with urlopen(req, timeout=120) as src:
                            if src.status != 206 or src.headers.get("Content-Range") != f"bytes {start}-{end}/{SIZE}":
                                raise OSError("The server returned a different byte range")
                            pos = start
                            while pos <= end:
                                buf = src.read(min(512 * 1024, end + 1 - pos))
                                if not buf:
                                    raise OSError("Truncated range")
                                view = memoryview(buf)
                                while view:
                                    written = os.pwrite(data.fileno(), view, pos)
                                    if written == 0:
                                        raise OSError("Unable to write the downloaded bytes")
                                    pos += written
                                    view = view[written:]
                        return i
                    except (OSError, HTTPError, URLError) as error:
                        if attempt == 9:
                            raise
                        pause = min(60, 2 ** attempt)
                        if isinstance(error, HTTPError) and error.headers.get("Retry-After","").isdigit():
                            pause = max(pause, int(error.headers["Retry-After"]))
                        print(f"Range {i} retry {attempt+1} after {error}", flush=True)
                        time.sleep(pause)
            with ThreadPoolExecutor(workers) as pool:
                tasks = [pool.submit(transfer, i) for i in range(total) if i not in done]
                last = time.monotonic()
                for future in as_completed(tasks):
                    i = future.result()
                    ledger.write(f"{i}\n")
                    done.add(i)
                    if time.monotonic() - last > 20 or len(done) == total:
                        print(f"Received {len(done)}/{total} ranges", flush=True)
                        last = time.monotonic()
            os.fsync(data.fileno())
        if digest(part, "md5") != MD5:
            raise ValueError("Source MD5 failed. Inspect the range ledger before retrying.")
        part.replace(dest)
        print("Verified archive MD5", MD5, flush=True)


def inventory(directory, output):
    root = directory / "extracted"
    import h5py
    from collections import Counter
    sizes = json.loads((directory / ".rar-member-sizes.json").read_text(encoding="utf8"))
    present = {p.relative_to(root).as_posix() for p in root.rglob("*") if p.is_file()}
    if present != set(sizes):
        raise ValueError("Extracted paths differ from the verified archive's member list")
    records = []
    for name, size in sorted(sizes.items()):
        rel = Path(name)
        path = root / rel
        parts = rel.parts[1:]
        if path.stat().st_size != size:
            raise ValueError(f"Extracted member has incorrect length: {name}")
        sha, crc = hashlib.sha256(), 0
        with path.open("rb") as stream:
            for block in iter(lambda: stream.read(4 * 1024**2), b""):
                sha.update(block)
                crc = zlib.crc32(block, crc)
        rec = {"path": name, "bytes": size, "sha256": sha.hexdigest(), "crc32": f"{crc:08x}",
               "type": path.suffix.lower().lstrip(".") or "other"}
        if path.suffix.lower() == ".hdf5":
            with h5py.File(path, "r") as f:
                rec["datasets"] = {name: {"shape": list(obj.shape), "fields": list(obj.dtype.names or [])}
                                   for name,obj in f.items() if isinstance(obj, h5py.Dataset)}
        if parts:
            rec["section"] = parts[0]
        records.append(rec)
    if not any(r["type"] == "hdf5" for r in records):
        raise ValueError("No HDF5 data found in extracted archive")
    result = {"source": URL, "record": "https://zenodo.org/records/10960858",
              "license": "CC BY 4.0", "archive_bytes": SIZE, "archive_md5": MD5,
              "extracted_bytes": sum(sizes.values()), "counts": dict(Counter(r["type"] for r in records)), "files": records}
    output.parent.mkdir(parents=True, exist_ok=True)
    tmp = output.with_suffix(".tmp")
    tmp.write_text(json.dumps(result,separators=(",",":")) + "\n", encoding="utf8")
    tmp.replace(output)


def extract(directory):
    import libarchive
    import rarfile
    archive = directory / "dye_analysis.rar"
    if digest(archive, "md5") != MD5:
        raise ValueError("Archive checksum differs from Zenodo")
    root = (directory / "extracted").resolve()
    root.mkdir(parents=True, exist_ok=True)
    sizes, checksums, pending = {}, {}, set()
    # Read CRCs from this verified archive. A prior extraction or an independent
    # range reader may already have delivered some of the members.
    with rarfile.RarFile(archive) as rar:
        for item in rar.infolist():
            dest = (root / item.filename).resolve()
            if not dest.is_relative_to(root) or dest == root:
                raise ValueError(f"Invalid archive path: {item.filename}")
            if item.isdir():
                dest.mkdir(parents=True, exist_ok=True)
                continue
            if not item.is_file():
                raise ValueError(f"Unexpected archive member: {item.filename}")
            rel = dest.relative_to(root).as_posix()
            if rel in sizes or item.CRC is None:
                raise ValueError(f"Duplicate file or missing checksum: {rel}")
            sizes[rel], checksums[rel] = item.file_size, item.CRC
            if dest.is_file() and dest.stat().st_size == item.file_size:
                crc = 0
                with dest.open("rb") as stream:
                    for block in iter(lambda: stream.read(4 * 1024**2), b""):
                        crc = zlib.crc32(block, crc)
                if crc == item.CRC:
                    continue
            pending.add(rel)
    print(f"Member CRCs verified for {len(sizes)-len(pending)} existing files; extracting {len(pending)}.", flush=True)
    with libarchive.file_reader(str(archive)) as entries:
        for entry in entries:
            if not pending:
                break
            dest = (root / entry.pathname).resolve()
            if not dest.is_relative_to(root) or dest == root:
                raise ValueError(f"Invalid archive path: {entry.pathname}")
            if entry.isdir:
                dest.mkdir(parents=True, exist_ok=True)
                continue
            if not entry.isfile:
                raise ValueError(f"Unexpected entry type: {entry.pathname}")
            rel = dest.relative_to(root).as_posix()
            if rel not in pending:
                continue
            dest.parent.mkdir(parents=True, exist_ok=True)
            tmp = dest.with_name(dest.name + ".extracting")
            crc = 0
            with tmp.open("wb") as out:
                for chunk in entry.get_blocks():
                    out.write(chunk)
                    crc = zlib.crc32(chunk, crc)
            if tmp.stat().st_size != sizes[rel] or crc != checksums[rel]:
                raise ValueError(f"Archive member failed size/CRC32: {rel}")
            tmp.replace(dest)
            pending.remove(rel)
    if pending:
        raise ValueError("Some archive members could not be extracted: " + ", ".join(sorted(pending)))
    listing = directory / ".rar-member-sizes.json"
    tmp = listing.with_suffix(".tmp")
    tmp.write_text(json.dumps(sizes, sort_keys=True) + "\n", encoding="utf8")
    tmp.replace(listing)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=["fetch", "extract", "inventory"])
    parser.add_argument("directory", type=Path)
    parser.add_argument("--workers", type=int, default=8)
    parser.add_argument("--output", type=Path, default=Path("docs/data/palette-manifest.json"))
    args = parser.parse_args()
    if args.action == "fetch":
        fetch(args.directory, args.workers)
    elif args.action == "extract":
        extract(args.directory)
    else:
        inventory(args.directory, args.output)
