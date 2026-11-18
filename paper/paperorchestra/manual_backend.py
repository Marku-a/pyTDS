"""File-queue Gemini stand-in for PaperOrchestra.

When GEMINI_BACKEND=manual, PaperOrchestra's genai client is replaced by this
one. Every generate_content() call is written to $PO_QUEUE/req_NNNN.md (text
parts; images/PDFs saved next to it; large repeated texts stored once under
blobs/ and referenced by hash). The call then blocks until a response file
$PO_QUEUE/resp_NNNN.txt appears, whose text is returned as response.text.

This lets an external model (here: Claude in a Claude Code session) act as the
LLM behind PaperOrchestra's unchanged agents, prompts and control flow.
"""
import hashlib
import json
import os
import threading
import time
from pathlib import Path

_lock = threading.Lock()
_counter = [0]
BLOB_MIN = 1500  # chars; texts longer than this are deduplicated into blobs/


class _Resp:
    def __init__(self, text):
        self.text = text


def _queue() -> Path:
    q = Path(os.environ.get("PO_QUEUE", "po_queue")).resolve()
    (q / "blobs").mkdir(parents=True, exist_ok=True)
    (q / "files").mkdir(parents=True, exist_ok=True)
    return q


def _text_block(q: Path, text: str) -> str:
    if len(text) < BLOB_MIN:
        return text
    h = hashlib.sha1(text.encode()).hexdigest()[:12]
    p = q / "blobs" / f"{h}.txt"
    if not p.exists():
        p.write_text(text)
    first = " ".join(text[:160].split())
    return f"[[BLOB {h} | {len(text)} chars | starts: {first} ...]]"


def _file_block(q: Path, n: int, k: int, data: bytes, mime: str) -> str:
    ext = {"application/pdf": "pdf", "image/png": "png", "image/jpeg": "jpg"}.get(mime, "bin")
    h = hashlib.sha1(data).hexdigest()[:12]
    p = q / "files" / f"{h}.{ext}"
    if not p.exists():
        p.write_bytes(data)
    return f"[[FILE files/{p.name} | {mime} | {len(data)} bytes]]"


def _serialize(q: Path, n: int, contents) -> list:
    out, k = [], [0]

    def walk(x):
        if x is None:
            return
        if isinstance(x, str):
            out.append(_text_block(q, x))
        elif isinstance(x, (list, tuple)):
            for y in x:
                walk(y)
        elif hasattr(x, "parts") and x.parts is not None:  # types.Content
            walk(list(x.parts))
        elif hasattr(x, "inline_data") or hasattr(x, "text"):  # types.Part
            if getattr(x, "text", None):
                out.append(_text_block(q, x.text))
            elif getattr(x, "inline_data", None) is not None:
                k[0] += 1
                out.append(_file_block(q, n, k[0], x.inline_data.data, x.inline_data.mime_type))
            elif getattr(x, "file_data", None) is not None:
                out.append(f"[[FILE-URI {x.file_data.file_uri}]]")
        elif hasattr(x, "save"):  # PIL image
            import io
            b = io.BytesIO()
            x.save(b, format="PNG")
            k[0] += 1
            out.append(_file_block(q, n, k[0], b.getvalue(), "image/png"))
        else:
            out.append(_text_block(q, str(x)))

    walk(contents)
    return out


class _Models:
    def generate_content(self, model=None, contents=None, config=None):
        q = _queue()
        with _lock:
            _counter[0] += 1
            n = _counter[0]
            while (q / f"req_{n:04d}.md").exists():
                _counter[0] += 1
                n = _counter[0]
            req = q / f"req_{n:04d}.md"
            req.write_text("")  # reserve number
        meta = {"model": model}
        sys_inst = None
        if config is not None:
            sys_inst = getattr(config, "system_instruction", None)
            tools = getattr(config, "tools", None)
            if tools:
                meta["tools"] = [
                    "google_search" if getattr(t, "google_search", None) is not None else str(t)
                    for t in tools
                ]
            for a in ("temperature", "response_mime_type"):
                v = getattr(config, a, None)
                if v is not None:
                    meta[a] = v
        parts = []
        if sys_inst is not None:
            parts.append("## SYSTEM INSTRUCTION\n" + "\n".join(_serialize(q, n, sys_inst)))
        parts.append("## USER CONTENTS\n" + "\n\n---\n\n".join(_serialize(q, n, contents)))
        req.write_text(f"# REQUEST {n:04d}\n```json\n{json.dumps(meta)}\n```\n\n" + "\n\n".join(parts))
        resp = q / f"resp_{n:04d}.txt"
        t0 = time.time()
        while not resp.exists():
            time.sleep(2)
            if time.time() - t0 > float(os.environ.get("PO_TIMEOUT", "21600")):
                raise TimeoutError(f"no response for request {n}")
        time.sleep(0.5)  # let the writer finish
        return _Resp(resp.read_text())


class ManualClient:
    def __init__(self):
        self.models = _Models()
