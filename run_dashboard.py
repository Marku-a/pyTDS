"""
pyTDS Test Dashboard
====================
Runs the full test suite and writes results to an Excel spreadsheet.

Usage:
    python run_dashboard.py

Output:
    tests reports/test_dashboard_<timestamp>.xlsx
"""

import sys
import os
import datetime

import pytest
import openpyxl
from openpyxl.styles import (
    PatternFill, Font, Alignment, Border, Side
)
from openpyxl.utils import get_column_letter


# ── Pytest plugin — collect results during the run ─────────────────────────────

class ResultCollector:
    """Minimal pytest plugin that captures each test outcome."""

    def __init__(self):
        self.results = []   # list of dicts

    def pytest_runtest_logreport(self, report):
        if report.when != "call":
            # Also capture setup failures (e.g. import errors)
            if report.when == "setup" and report.failed:
                self.results.append({
                    "file":   report.nodeid.split("::")[0],
                    "name":   report.nodeid.split("::")[-1],
                    "status": "ERROR",
                    "reason": str(report.longrepr),
                })
            return

        status = "PASSED" if report.passed else (
            "FAILED" if report.failed else "SKIPPED"
        )
        reason = ""
        if report.failed:
            reason = str(report.longrepr)
            # Keep only the short assertion message (last non-empty line)
            lines = [l.strip() for l in reason.splitlines() if l.strip()]
            reason = lines[-1] if lines else reason

        self.results.append({
            "file":   report.nodeid.split("::")[0],
            "name":   report.nodeid.split("::")[-1],
            "status": status,
            "reason": reason,
        })


# ── Excel writer ───────────────────────────────────────────────────────────────

# Colour palette
GREEN_FILL   = PatternFill("solid", fgColor="C6EFCE")
RED_FILL     = PatternFill("solid", fgColor="FFC7CE")
YELLOW_FILL  = PatternFill("solid", fgColor="FFEB9C")
HEADER_FILL  = PatternFill("solid", fgColor="1F4E79")
SECTION_FILL = PatternFill("solid", fgColor="D6E4F0")

GREEN_FONT   = Font(color="276221", bold=False)
RED_FONT     = Font(color="9C0006", bold=False)
YELLOW_FONT  = Font(color="9C6500", bold=False)
HEADER_FONT  = Font(color="FFFFFF", bold=True, size=11)
SECTION_FONT = Font(color="1F4E79", bold=True)

THIN = Side(style="thin", color="BBBBBB")
THIN_BORDER = Border(left=THIN, right=THIN, top=THIN, bottom=THIN)

STATUS_ICON = {"PASSED": "✅", "FAILED": "❌", "SKIPPED": "⏭", "ERROR": "💥"}


def _cell(ws, row, col, value="", fill=None, font=None, bold=False,
          align="left", wrap=False, border=True):
    c = ws.cell(row=row, column=col, value=value)
    if fill:
        c.fill = fill
    if font:
        c.font = font
    elif bold:
        c.font = Font(bold=True)
    c.alignment = Alignment(horizontal=align, vertical="center",
                             wrap_text=wrap)
    if border:
        c.border = THIN_BORDER
    return c


def write_excel(results, out_path):
    wb = openpyxl.Workbook()
    ws = wb.active
    ws.title = "Test Dashboard"

    now = datetime.datetime.now().strftime("%Y-%m-%d  %H:%M")
    n_passed  = sum(1 for r in results if r["status"] == "PASSED")
    n_failed  = sum(1 for r in results if r["status"] == "FAILED")
    n_error   = sum(1 for r in results if r["status"] == "ERROR")
    n_skipped = sum(1 for r in results if r["status"] == "SKIPPED")
    n_total   = len(results)

    # ── Title row ─────────────────────────────────────────────────────────────
    ws.merge_cells("A1:E1")
    c = ws["A1"]
    c.value = f"pyTDS — Test Dashboard   |   {now}"
    c.fill   = HEADER_FILL
    c.font   = Font(color="FFFFFF", bold=True, size=13)
    c.alignment = Alignment(horizontal="center", vertical="center")
    ws.row_dimensions[1].height = 28

    # ── Summary row ───────────────────────────────────────────────────────────
    ws.merge_cells("A2:E2")
    summary_color = GREEN_FILL if n_failed == 0 and n_error == 0 else RED_FILL
    summary_font  = GREEN_FONT if n_failed == 0 and n_error == 0 else Font(color="9C0006", bold=True)
    overall = "ALL PASSED" if n_failed == 0 and n_error == 0 else f"{n_failed + n_error} FAILED"
    c = ws["A2"]
    c.value = (f"  {overall}   |   "
               f"✅ {n_passed} passed   "
               f"❌ {n_failed + n_error} failed   "
               f"⏭ {n_skipped} skipped   "
               f"Total: {n_total}")
    c.fill = summary_color
    c.font = summary_font
    c.alignment = Alignment(horizontal="center", vertical="center")
    ws.row_dimensions[2].height = 22

    # ── Column headers ────────────────────────────────────────────────────────
    headers = ["#", "Test File", "Test Name", "Status", "Failure Reason"]
    col_widths = [5, 22, 48, 10, 60]

    for col, (h, w) in enumerate(zip(headers, col_widths), start=1):
        _cell(ws, 3, col, h, fill=HEADER_FILL, font=HEADER_FONT, align="center")
        ws.column_dimensions[get_column_letter(col)].width = w

    ws.row_dimensions[3].height = 18

    # ── Test rows, grouped by file ────────────────────────────────────────────
    row = 4
    current_file = None
    counter = 0

    for r in results:
        # Section header when file changes
        if r["file"] != current_file:
            current_file = r["file"]
            ws.merge_cells(f"A{row}:E{row}")
            _cell(ws, row, 1, f"  {current_file}",
                  fill=SECTION_FILL, font=SECTION_FONT)
            # Overwrite merged cells' borders
            for col in range(2, 6):
                ws.cell(row=row, column=col).fill = SECTION_FILL
            ws.row_dimensions[row].height = 16
            row += 1

        counter += 1
        status = r["status"]
        icon   = STATUS_ICON.get(status, "?")

        fill = (GREEN_FILL  if status == "PASSED"  else
                RED_FILL    if status in ("FAILED", "ERROR") else
                YELLOW_FILL)
        font = (GREEN_FONT  if status == "PASSED"  else
                RED_FONT    if status in ("FAILED", "ERROR") else
                YELLOW_FONT)

        _cell(ws, row, 1, counter,         align="center")
        _cell(ws, row, 2, r["file"])
        _cell(ws, row, 3, r["name"])
        _cell(ws, row, 4, f"{icon} {status}", fill=fill, font=font, align="center")
        _cell(ws, row, 5, r["reason"],     fill=fill if r["reason"] else None,
              font=font if r["reason"] else None, wrap=True)

        # Auto-height for long failure messages
        if r["reason"]:
            lines = max(1, len(r["reason"]) // 80 + 1)
            ws.row_dimensions[row].height = max(15, lines * 14)
        else:
            ws.row_dimensions[row].height = 15

        row += 1

    # Freeze panes below header
    ws.freeze_panes = "A4"

    wb.save(out_path)
    return out_path


# ── Main ───────────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    project_root = os.path.dirname(os.path.abspath(__file__))
    tests_dir    = os.path.join(project_root, "tests")
    reports_dir  = os.path.join(project_root, "tests reports")
    os.makedirs(reports_dir, exist_ok=True)
    timestamp    = datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
    out_path     = os.path.join(reports_dir, f"test_dashboard_{timestamp}.xlsx")

    collector = ResultCollector()

    print("Running tests...")
    exit_code = pytest.main(
        [tests_dir, "-v", "--tb=short", "--no-header", "-q"],
        plugins=[collector],
    )

    print(f"\nWriting dashboard to: {out_path}")
    write_excel(collector.results, out_path)

    n_failed = sum(1 for r in collector.results if r["status"] in ("FAILED", "ERROR"))
    n_total  = len(collector.results)
    print(f"Dashboard saved — {n_total - n_failed}/{n_total} passed")

    sys.exit(exit_code)
