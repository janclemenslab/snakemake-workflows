"""FastHTML dashboard for Snakemake controller logs."""

from __future__ import annotations

from html import escape
import os
from pathlib import Path
import sys

if __package__ in {None, ""}:
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from fasthtml.common import *  # noqa: F403
from monsterui.all import *  # noqa: F403

from fab_commands.monitor import (
    apply_squeue_state,
    controller_job_log_context,
    controller_logs,
    fetch_squeue_output,
    job_rows,
    parse_controller_log,
    read_log_tail,
    resolve_log_paths,
)


PROJECT_DIR = Path(os.environ.get("FAB_MONITOR_PROJECT_DIR", ".")).resolve()
DEFAULT_USER = os.environ.get("FAB_MONITOR_USER", "").strip()
DEFAULT_REMOTE_HOST = os.environ.get("FAB_MONITOR_REMOTE_HOST", "").strip()
SORTABLE_COLUMNS = {
    "state": "state",
    "rule": "rule",
    "wildcards": "wildcards",
    "elapsed": "elapsed",
    "node_or_reason": "node_or_reason",
    "log_file": "log_file",
}
LOG_BROWSER_CSS = """
@keyframes fab-monitor-spin {
    to { transform: rotate(360deg); }
}

.fab-monitor-loading-overlay {
    position: fixed;
    inset: 0;
    z-index: 1200;
    display: flex;
    align-items: center;
    justify-content: center;
    padding: 1rem;
    background: rgba(15, 23, 42, 0.38);
}

.fab-monitor-loading-window {
    width: min(28rem, 96vw);
    border: 1px solid #cbd5e1;
    border-radius: 1rem;
    background: white;
    box-shadow: 0 24px 80px rgba(15, 23, 42, 0.22);
    padding: 1rem 1.25rem;
}

.fab-monitor-spinner {
    width: 1.5rem;
    height: 1.5rem;
    border-radius: 999px;
    border: 3px solid #cbd5e1;
    border-top-color: #0f766e;
    animation: fab-monitor-spin 0.8s linear infinite;
    flex: none;
}

.fab-monitor-log-tabs {
    display: flex;
    gap: 0.5rem;
    overflow-x: auto;
    padding-bottom: 0.25rem;
}

.fab-monitor-log-tab {
    min-width: 12rem;
    border: 1px solid #cbd5e1;
    border-radius: 0.85rem;
    background: #f8fafc;
    padding: 0.75rem 0.9rem;
    text-align: left;
    transition: border-color 120ms ease, background-color 120ms ease, box-shadow 120ms ease;
}

.fab-monitor-log-tab.is-active {
    border-color: #0f766e;
    background: #ecfeff;
    box-shadow: inset 0 0 0 1px rgba(15, 118, 110, 0.12);
}

.fab-monitor-log-panel {
    display: none;
}

.fab-monitor-log-panel.is-active {
    display: block;
}

.fab-monitor-log-pre {
    width: 100%;
    overflow: auto;
    border-radius: 0.75rem;
    border: 1px solid #cbd5e1;
    padding: 0.875rem;
    white-space: pre-wrap;
    font-family: ui-monospace, SFMono-Regular, Menlo, monospace;
    font-size: 12px;
    line-height: 1.5;
    background: white;
}
"""
LOG_BROWSER_SCRIPT = """
window.fabMonitorDashboard = window.fabMonitorDashboard || {
  pollPausedUntil: 0,
  stateOrder: ['running', 'pending', 'submitted', 'finished', 'failed'],
  collectDashboardValues() {
    const values = {};
    ['controls-form', 'jobs-filter-form'].forEach((formId) => {
      const form = document.getElementById(formId);
      if (!form) return;
      new FormData(form).forEach((value, key) => {
        values[key] = value;
      });
    });
    return values;
  },
  refreshDashboardShell() {
    if (typeof htmx === 'undefined') return;
    const target = document.getElementById('dashboard-shell');
    if (!target) return;
    htmx.ajax('GET', '/dashboard-shell', {
      target: '#dashboard-shell',
      swap: 'outerHTML',
      values: this.collectDashboardValues()
    });
  },
  showLoading() {
    const container = document.getElementById('modal-container');
    const template = document.getElementById('job-modal-loading-template');
    if (!container || !template) return;
    container.innerHTML = template.innerHTML;
  },
  pausePolling(ms = 8000) {
    this.pollPausedUntil = Math.max(this.pollPausedUntil, Date.now() + ms);
  },
  shouldPoll() {
    return Date.now() >= this.pollPausedUntil;
  },
  jobRows() {
    return Array.from(document.querySelectorAll('#jobs-card tbody tr[data-job-row="true"]'));
  },
  sortStateValues(values) {
    const unique = Array.from(new Set(values.filter(Boolean)));
    unique.sort((left, right) => {
      const leftIndex = this.stateOrder.indexOf(left);
      const rightIndex = this.stateOrder.indexOf(right);
      if (leftIndex !== -1 || rightIndex !== -1) {
        const normalizedLeft = leftIndex === -1 ? Number.MAX_SAFE_INTEGER : leftIndex;
        const normalizedRight = rightIndex === -1 ? Number.MAX_SAFE_INTEGER : rightIndex;
        if (normalizedLeft !== normalizedRight) {
          return normalizedLeft - normalizedRight;
        }
      }
      return left.localeCompare(right);
    });
    return unique;
  },
  sortTextValues(values) {
    return Array.from(new Set(values.filter(Boolean))).sort((left, right) =>
      left.localeCompare(right, undefined, { sensitivity: 'base' })
    );
  },
  setSelectOptions(select, allLabel, values, selectedValue) {
    if (!select) return '';
    const nextSelected = values.includes(selectedValue) ? selectedValue : '';
    const optionSignature = JSON.stringify([[ '', allLabel ], ...values.map((value) => [value, value])]);
    if (select.dataset.optionSignature !== optionSignature) {
      select.innerHTML = '';
      select.appendChild(new Option(allLabel, ''));
      values.forEach((value) => {
        select.appendChild(new Option(value, value));
      });
      select.dataset.optionSignature = optionSignature;
    }
    select.value = nextSelected;
    return nextSelected;
  },
  syncJobsFilterOptions() {
    const rows = this.jobRows();
    const stateSelect = document.getElementById('state_filter');
    const ruleSelect = document.getElementById('rule_filter');
    if (!stateSelect || !ruleSelect) {
      return { state: '', rule: '' };
    }

    let stateValue = stateSelect.value || '';
    let ruleValue = ruleSelect.value || '';

    const statesForRule = ruleValue
      ? rows.filter((row) => row.dataset.rule === ruleValue).map((row) => row.dataset.baseState)
      : rows.map((row) => row.dataset.baseState);
    stateValue = this.setSelectOptions(
      stateSelect,
      'All states',
      this.sortStateValues(statesForRule),
      stateValue
    );

    const rulesForState = stateValue
      ? rows.filter((row) => row.dataset.baseState === stateValue).map((row) => row.dataset.rule)
      : rows.map((row) => row.dataset.rule);
    ruleValue = this.setSelectOptions(
      ruleSelect,
      'All rules',
      this.sortTextValues(rulesForState),
      ruleValue
    );

    const adjustedStatesForRule = ruleValue
      ? rows.filter((row) => row.dataset.rule === ruleValue).map((row) => row.dataset.baseState)
      : rows.map((row) => row.dataset.baseState);
    stateValue = this.setSelectOptions(
      stateSelect,
      'All states',
      this.sortStateValues(adjustedStatesForRule),
      stateValue
    );

    const adjustedRulesForState = stateValue
      ? rows.filter((row) => row.dataset.baseState === stateValue).map((row) => row.dataset.rule)
      : rows.map((row) => row.dataset.rule);
    ruleValue = this.setSelectOptions(
      ruleSelect,
      'All rules',
      this.sortTextValues(adjustedRulesForState),
      ruleValue
    );

    const clearButton = document.getElementById('jobs-clear-filters');
    if (clearButton) {
      clearButton.hidden = !(stateValue || ruleValue);
    }

    return { state: stateValue, rule: ruleValue };
  },
  renderJobsTable({ resetPage = false } = {}) {
    const jobsCard = document.getElementById('jobs-card');
    if (!jobsCard) return;

    const pageInput = document.getElementById('page');
    const limitInput = document.getElementById('limit');
    if (resetPage && pageInput) {
      pageInput.value = '1';
    }

    const { state, rule } = this.syncJobsFilterOptions();
    const rows = this.jobRows();
    const matchedRows = rows.filter((row) => {
      if (state && row.dataset.baseState !== state) return false;
      if (rule && row.dataset.rule !== rule) return false;
      return true;
    });

    let limit = parseInt(limitInput?.value || '100', 10);
    if (!Number.isFinite(limit) || limit < 1) {
      limit = 100;
    }

    let page = parseInt(pageInput?.value || '1', 10);
    if (!Number.isFinite(page) || page < 1) {
      page = 1;
    }

    const pageCount = Math.max(1, Math.ceil(matchedRows.length / limit));
    if (page > pageCount) {
      page = pageCount;
    }
    if (pageInput) {
      pageInput.value = String(page);
    }

    const startIndex = matchedRows.length ? (page - 1) * limit : 0;
    const endIndex = Math.min(startIndex + limit, matchedRows.length);
    const visibleRows = new Set(matchedRows.slice(startIndex, endIndex));

    rows.forEach((row) => {
      row.hidden = !visibleRows.has(row);
    });

    const summary = document.getElementById('jobs-summary-text');
    if (summary) {
      const text = matchedRows.length
        ? `Showing ${startIndex + 1}-${endIndex} of ${matchedRows.length}`
        : 'Showing 0 jobs';
      summary.textContent =
        matchedRows.length === rows.length ? text : `${text} | ${rows.length} total`;
    }

    const pageStatus = document.getElementById('jobs-page-status');
    if (pageStatus) {
      pageStatus.textContent = `Page ${page} of ${pageCount} | ${matchedRows.length} jobs`;
    }

    const prevButton = document.getElementById('jobs-prev-page');
    if (prevButton) {
      prevButton.disabled = page <= 1 || matchedRows.length === 0;
    }
    const nextButton = document.getElementById('jobs-next-page');
    if (nextButton) {
      nextButton.disabled = page >= pageCount || matchedRows.length === 0;
    }

    document.querySelectorAll('[data-page-size-button]').forEach((button) => {
      const isActive = button.dataset.pageSizeButton === String(limit);
      button.classList.toggle('uk-btn-primary', isActive);
      button.classList.toggle('uk-btn-ghost', !isActive);
    });

    const emptyState = document.getElementById('jobs-empty-state');
    if (emptyState) {
      emptyState.hidden = matchedRows.length > 0;
    }

    const tableWrapper = document.getElementById('jobs-table-wrapper');
    if (tableWrapper) {
      tableWrapper.hidden = matchedRows.length === 0;
    }
  },
  applyJobsFilters() {
    this.pausePolling(2500);
    this.closeTransientModal();
    this.renderJobsTable({ resetPage: true });
  },
  clearJobsFilters() {
    this.pausePolling(2500);
    this.closeTransientModal();
    const stateSelect = document.getElementById('state_filter');
    const ruleSelect = document.getElementById('rule_filter');
    if (stateSelect) stateSelect.value = '';
    if (ruleSelect) ruleSelect.value = '';
    this.renderJobsTable({ resetPage: true });
  },
  changeJobsPage(delta) {
    this.pausePolling(2500);
    this.closeTransientModal();
    const pageInput = document.getElementById('page');
    if (pageInput) {
      const currentPage = parseInt(pageInput.value || '1', 10);
      pageInput.value = String(Math.max(1, currentPage + delta));
    }
    this.renderJobsTable();
  },
  setJobsPageSize(size) {
    this.pausePolling(2500);
    this.closeTransientModal();
    const limitInput = document.getElementById('limit');
    const pageInput = document.getElementById('page');
    if (limitInput) limitInput.value = String(size);
    if (pageInput) pageInput.value = '1';
    this.renderJobsTable();
  },
  closeTransientModal() {
    const container = document.getElementById('modal-container');
    if (!container) return;
    container.innerHTML = '';
  },
  selectLogTab(rootId, tabId) {
    const root = document.getElementById(rootId);
    if (!root) return;
    root.querySelectorAll('[data-log-tab]').forEach((button) => {
      const active = button.dataset.tabId === tabId;
      button.classList.toggle('is-active', active);
      button.setAttribute('aria-selected', active ? 'true' : 'false');
    });
    root.querySelectorAll('[data-log-panel]').forEach((panel) => {
      const active = panel.dataset.tabId === tabId;
      panel.classList.toggle('is-active', active);
      panel.hidden = !active;
    });
  }
};

if (!window.fabMonitorDashboardBound) {
  window.fabMonitorDashboardBound = true;
  document.addEventListener('DOMContentLoaded', () => {
    window.fabMonitorDashboard.renderJobsTable();
  });
  document.body.addEventListener('htmx:afterSwap', (event) => {
    if (event.target && event.target.id === 'dashboard-shell') {
      window.fabMonitorDashboard.renderJobsTable();
    }
  });
}
"""


def _env_int(name: str, default: int, *, minimum: int = 0) -> int:
    try:
        value = int(os.environ.get(name, default))
    except (TypeError, ValueError):
        value = default
    return max(minimum, value)


DEFAULT_LIMIT = _env_int("FAB_MONITOR_LIMIT", 100, minimum=10)
DEFAULT_REFRESH = _env_int("FAB_MONITOR_REFRESH", 5, minimum=0)
DEFAULT_HOST = os.environ.get("FAB_MONITOR_BIND_HOST", "127.0.0.1").strip() or "127.0.0.1"
DEFAULT_PORT = _env_int("FAB_MONITOR_PORT", 2718, minimum=1)
PAGE_SIZE_OPTIONS = [50, 100, 200]
DASHBOARD_FORM_INCLUDE = "#controls-form, #jobs-filter-form"
DASHBOARD_REQUEST_SYNC = "#dashboard-shell:replace"


app, rt = fast_app(  # noqa: F405
    pico=False,
    live=False,
    title="Snakemake Monitor",
    hdrs=Theme.slate.headers(),  # noqa: F405
)


def _available_logs() -> list[Path]:
    return controller_logs(PROJECT_DIR)


def _normalize_limit(limit: int | str | None) -> int:
    try:
        value = int(limit if limit is not None else DEFAULT_LIMIT)
    except (TypeError, ValueError):
        value = DEFAULT_LIMIT
    return max(10, min(value, 200))


def _normalize_page(page: int | str | None) -> int:
    try:
        value = int(page if page is not None else 1)
    except (TypeError, ValueError):
        value = 1
    return max(1, value)


def _normalize_refresh(refresh: int | str | None) -> int:
    try:
        value = int(refresh if refresh is not None else DEFAULT_REFRESH)
    except (TypeError, ValueError):
        value = DEFAULT_REFRESH
    return max(0, min(value, 60))


def _normalize_state_filter(state_filter: str | None) -> str:
    return (state_filter or "").strip().lower()


def _normalize_rule_filter(rule_filter: str | None) -> str:
    return (rule_filter or "").strip()


def _selected_log_path(controller_log: str | None) -> str:
    logs = _available_logs()
    if controller_log:
        candidate = Path(controller_log).expanduser()
        if candidate.exists():
            return str(candidate.resolve())
    if logs:
        return str(logs[-1])
    return ""


def _load_run(
    controller_log: str,
    user: str,
    remote_host: str,
    *,
    smk_jobid: str = "",
):
    selected_log = _selected_log_path(controller_log)
    if not selected_log:
        return None, "", ""

    run = parse_controller_log(Path(selected_log))
    queue_warning = ""
    queue_user = (user or "").strip()
    queue_host = (remote_host or "").strip()
    needs_live_slurm = (
        run.job_needs_live_slurm_state(smk_jobid)
        if smk_jobid
        else run.needs_live_slurm_state()
    )

    if needs_live_slurm and queue_user and queue_host:
        squeue_output = fetch_squeue_output(user=queue_user, host=queue_host)
        if squeue_output:
            apply_squeue_state(run, squeue_output)
        else:
            queue_warning = "Could not fetch live SLURM state."
    elif needs_live_slurm and queue_user and not queue_host:
        queue_warning = "Set an SSH host to enrich jobs with live SLURM state."

    return run, queue_warning, selected_log


def _elapsed_sort_value(value: str) -> int:
    if not value or value == "-":
        return -1
    days = 0
    clock = value
    if "-" in value:
        day_text, clock = value.split("-", maxsplit=1)
        if day_text.isdigit():
            days = int(day_text)
    parts = clock.split(":")
    if len(parts) != 3 or not all(part.isdigit() for part in parts):
        return -1
    hours, minutes, seconds = (int(part) for part in parts)
    return (((days * 24) + hours) * 60 + minutes) * 60 + seconds


def _sort_value(row: dict[str, str], sort_by: str):
    if sort_by == "elapsed":
        return (_elapsed_sort_value(row["elapsed"]), row["smk_jobid"])
    return (row.get(sort_by, "").lower(), row["smk_jobid"])


def _sorted_rows(run, sort_by: str, sort_dir: str):
    rows = job_rows(run, limit=0)
    if sort_by in SORTABLE_COLUMNS:
        rows = sorted(rows, key=lambda row: _sort_value(row, sort_by), reverse=sort_dir == "desc")
    return rows


def _row_base_state(row: dict[str, str]) -> str:
    return row["state"].split("/", maxsplit=1)[0].strip().lower()


def _available_states(rows: list[dict[str, str]]) -> list[str]:
    preferred_order = ["running", "pending", "submitted", "finished", "failed"]
    states = {_row_base_state(row) for row in rows if row["state"] and row["state"] != "-"}
    ordered = [state for state in preferred_order if state in states]
    ordered.extend(sorted(state for state in states if state not in preferred_order))
    return ordered


def _available_rules(rows: list[dict[str, str]]) -> list[str]:
    rules = {row["rule"] for row in rows if row["rule"] and row["rule"] != "-"}
    return sorted(rules, key=str.casefold)


def _filter_rows(
    rows: list[dict[str, str]],
    *,
    state_filter: str,
    rule_filter: str,
) -> list[dict[str, str]]:
    filtered = rows
    if state_filter:
        filtered = [row for row in filtered if _row_base_state(row) == state_filter]
    if rule_filter:
        filtered = [row for row in filtered if row["rule"] == rule_filter]
    return filtered


def _warning_card(title: str, messages: list[str], tone: str):
    tone_cls = {
        "warn": "border-orange-300 bg-orange-50",
        "danger": "border-red-300 bg-red-50",
    }.get(tone, "border-slate-300 bg-slate-50")
    return Card(  # noqa: F405
        Ul(*[Li(message) for message in messages], cls="list-disc pl-5 space-y-1"),  # noqa: F405
        header=H4(title, cls="text-sm font-semibold"),  # noqa: F405
        cls=f"h-full border {tone_cls}",
        body_cls="space-y-3",
    )


def _state_badge(state: str):
    base_state = state.split("/", maxsplit=1)[0].strip().lower()
    tone_cls = {
        "failed": "border-red-300 bg-red-50 text-red-700",
        "error": "border-red-300 bg-red-50 text-red-700",
        "finished": "border-emerald-300 bg-emerald-50 text-emerald-700",
    }.get(base_state, "border-slate-300 bg-slate-50 text-slate-700")
    return CodeSpan(  # noqa: F405
        state,
        cls=f"rounded border px-2 py-0.5 text-xs font-medium {tone_cls}",
    )


def _controls(
    controller_log: str,
    user: str,
    remote_host: str,
    limit: int,
    page: int,
    refresh: int,
    sort_by: str,
    sort_dir: str,
    controls_open: bool = False,
):
    logs = list(reversed(_available_logs()))

    refresh_options = [
        ("Off", 0),
        ("3s", 3),
        ("5s", 5),
        ("10s", 10),
        ("30s", 30),
    ]

    def native_select(
        *,
        select_id: str,
        name: str,
        options: list[tuple[str, str, bool]],
        title: str,
        class_name: str,
        style: str = "",
        disabled: bool = False,
    ):
        options_html = "".join(
            f'<option value="{escape(value, quote=True)}"'
            f'{" selected" if selected else ""}>{escape(label)}</option>'
            for value, label, selected in options
        )
        attrs = {
            "id": select_id,
            "name": name,
            "class": class_name,
            "title": title,
            "onfocus": "window.fabMonitorDashboard.pausePolling(12000);",
            "onpointerdown": "window.fabMonitorDashboard.pausePolling(12000);",
            "onchange": (
                "window.fabMonitorDashboard.pausePolling(2500);"
                "window.fabMonitorDashboard.closeTransientModal();"
                "document.getElementById('page').value='1';"
            ),
            "hx-get": "/dashboard-shell",
            "hx-trigger": "change",
            "hx-target": "#dashboard-shell",
            "hx-swap": "outerHTML",
            "hx-include": DASHBOARD_FORM_INCLUDE,
            "hx-sync": DASHBOARD_REQUEST_SYNC,
        }
        if style:
            attrs["style"] = style
        if disabled:
            attrs["disabled"] = "disabled"
        attr_html = " ".join(
            f'{key}="{escape(value, quote=True)}"' for key, value in attrs.items()
        )
        return NotStr(f"<select {attr_html}>{options_html}</select>")  # noqa: F405

    def native_input(
        *,
        input_id: str,
        name: str,
        value: str,
        placeholder: str,
        title: str,
        class_name: str = "uk-input uk-form-small w-full",
    ):
        attrs = {
            "id": input_id,
            "name": name,
            "type": "text",
            "value": value,
            "placeholder": placeholder,
            "title": title,
            "autocomplete": "off",
            "class": class_name,
            "onfocus": "window.fabMonitorDashboard.pausePolling(12000);",
            "onpointerdown": "window.fabMonitorDashboard.pausePolling(12000);",
            "onchange": (
                "window.fabMonitorDashboard.pausePolling(2500);"
                "window.fabMonitorDashboard.closeTransientModal();"
                "document.getElementById('page').value='1';"
            ),
            "hx-get": "/dashboard-shell",
            "hx-trigger": "change",
            "hx-target": "#dashboard-shell",
            "hx-swap": "outerHTML",
            "hx-include": DASHBOARD_FORM_INCLUDE,
            "hx-sync": DASHBOARD_REQUEST_SYNC,
        }
        attr_html = " ".join(
            f'{key}="{escape(attr_value, quote=True)}"' for key, attr_value in attrs.items()
        )
        return NotStr(f"<input {attr_html}>")  # noqa: F405

    run_select = Div(
        DivFullySpaced(  # noqa: F405
            P("Run", cls="text-sm font-medium"),
            Button(  # noqa: F405
                UkIcon("refresh-cw", width=14, height=14),  # noqa: F405
                submit=False,
                cls=(ButtonT.ghost, ButtonT.sm),  # noqa: F405
                title="Refresh controller log list",
                hx_get="/controls?controls_open=1",
                hx_target="#controls-card",
                hx_swap="outerHTML",
                hx_include="#controls-form",
            ),
        ),
        native_select(
            select_id="controller_log",
            name="controller_log",
            options=[
                (str(path), path.name, str(path) == controller_log) for path in logs
            ],
            title="Select controller log",
            class_name="uk-select uk-form-small w-full",
            style="width: 100%;",
            disabled=not logs,
        ),
        cls="space-y-2 min-w-0",
        style="flex: 2.2 1 24rem; min-width: 18rem;",
    )

    user_input = Div(
        P("User", cls="text-sm font-medium"),
        native_input(
            input_id="user",
            name="user",
            value=user,
            placeholder="abcd1234",
            title="SSH username for live SLURM state",
            class_name="uk-input uk-form-small w-full",
        ),
        cls="space-y-2 min-w-0",
        style="flex: 0.75 1 8.5rem; min-width: 8rem;",
    )

    host_input = Div(
        P("Host", cls="text-sm font-medium"),
        native_input(
            input_id="remote_host",
            name="remote_host",
            value=remote_host,
            placeholder="rosa.hpc.uni-oldenburg.de",
            title="SSH host for live SLURM state",
        ),
        cls="space-y-2 min-w-0",
        style="flex: 1.2 1 15rem; min-width: 12rem;",
    )

    refresh_select = Div(
        P("Refresh", cls="text-sm font-medium"),
        native_select(
            select_id="refresh",
            name="refresh",
            options=[
                (str(seconds), label, seconds == refresh)
                for label, seconds in refresh_options
            ],
            title="Auto-refresh interval",
            class_name="uk-select uk-form-small",
            style="width: 5.5rem;",
        ),
        cls="space-y-2",
        style="flex: 0 0 5.5rem;",
    )

    form = Form(  # noqa: F405
        Input(type="hidden", id="page", name="page", value=str(page)),  # noqa: F405
        Input(type="hidden", id="sort_by", name="sort_by", value=sort_by),  # noqa: F405
        Input(type="hidden", id="sort_dir", name="sort_dir", value=sort_dir),  # noqa: F405
        Input(type="hidden", id="limit", name="limit", value=str(limit)),  # noqa: F405
        Div(
            run_select,
            user_input,
            host_input,
            refresh_select,
            cls="flex flex-wrap xl:flex-nowrap gap-2 items-end pt-2",
        ),
        id="controls-form",
        cls="space-y-0",
    )

    summary_line = DivFullySpaced(  # noqa: F405
        Div(
            H2(PROJECT_DIR.name, cls="text-lg font-semibold tracking-tight"),  # noqa: F405
            P(Path(controller_log).name if controller_log else "No controller log selected", cls="text-xs opacity-70"),
            cls="space-y-0.5",
        ),
        P("Options", cls="text-sm opacity-60"),
    )

    return Card(  # noqa: F405
        Details(  # noqa: F405
            Summary(summary_line, cls="cursor-pointer list-none"),  # noqa: F405
            form,
            open=controls_open,
            cls="group",
        ),
        id="controls-card",
        body_cls="space-y-0 py-3",
    )


def _summary(run, queue_warning: str, selected_log: str):
    if run is None:
        return Card(  # noqa: F405
            P("No controller logs found under `log/slurm`.", cls="text-sm"),
            cls="border border-amber-300 bg-amber-50",
            body_cls="space-y-0 py-3",
        )

    counts = run.counts()
    progress = "-"
    if run.progress_total:
        progress = f"{run.progress_done}/{run.progress_total} ({run.progress_percent}%)"

    compact_row = DivFullySpaced(  # noqa: F405
        Div(
            H3(Path(selected_log).name, cls="text-base font-semibold"),  # noqa: F405
            P(
                f"{run.started_at} | {run.state} | {progress}",
                cls="text-sm opacity-75",
            ),
            cls="space-y-0.5",
        ),
        Div(
            Span("running ", Strong(str(counts["running"]))),  # noqa: F405
            Span("pending ", Strong(str(counts["submitted"] - counts["running"]))),  # noqa: F405
            Span("failed ", Strong(str(counts["failed"]))),  # noqa: F405
            cls="flex flex-wrap gap-4 text-sm whitespace-nowrap",
        ),
    )

    blocks = [compact_row]
    meta = []
    if run.slurm_run_id:
        meta.append(f"SLURM run ID: {run.slurm_run_id}")
    if run.status_line:
        meta.append(run.status_line)
    if meta:
        blocks.append(P(" | ".join(meta), cls="text-sm opacity-70 break-all"))
    warnings = []
    if queue_warning:
        warnings.append(_warning_card("Live Queue", [queue_warning], "warn"))
    error_messages = [message for message in run.errors if message and message != "None"]
    if error_messages:
        warnings.append(_warning_card("Recent Errors", error_messages[-5:], "danger"))
    if warnings:
        blocks.append(Grid(*warnings, cols_min=1, cols_lg=2, cls="gap-2"))  # noqa: F405
    return Card(*blocks, body_cls="space-y-2 py-3")  # noqa: F405


def _job_modal_loading_template():
    return Div(  # noqa: F405
        Div(
            Div(
                Button(  # noqa: F405
                    "Close",
                    submit=False,
                    cls=(ButtonT.ghost, ButtonT.sm),  # noqa: F405
                    onclick="window.fabMonitorDashboard.closeTransientModal();",
                ),
                cls="flex justify-end",
            ),
            Div(
                Span(cls="fab-monitor-spinner", aria_hidden="true"),  # noqa: F405
                Div(
                    P("Loading logs", cls="text-base font-semibold"),
                    P(
                        "Fetching controller, rule, and SLURM log output.",
                        cls="text-sm opacity-70",
                    ),
                    cls="space-y-1",
                ),
                cls="flex items-center gap-3",
                role="status",
                aria_live="polite",
            ),
            cls="fab-monitor-loading-window space-y-3",
            onclick="event.stopPropagation();",
        ),
        cls="fab-monitor-loading-overlay",
        onclick="window.fabMonitorDashboard.closeTransientModal();",
    )


def _log_link(row: dict[str, str]):
    if row["log_file"] == "-":
        return CodeSpan("-", cls="text-xs")  # noqa: F405
    return A(  # noqa: F405
        row["log_file"],
        href="#",
        cls="font-mono text-xs hover:underline",
        title=row["log_file"],
        onclick="window.fabMonitorDashboard.showLoading();",
        hx_get=f"/job-modal?smk_jobid={row['smk_jobid']}",
        hx_target="#modal-container",
        hx_swap="innerHTML",
        hx_include="#controls-form",
    )


def _sort_header(label: str, key: str, sort_by: str, sort_dir: str):
    marker = ""
    next_dir = "asc"
    if sort_by == key:
        marker = " ^" if sort_dir == "asc" else " v"
        next_dir = "desc" if sort_dir == "asc" else "asc"
    return Th(  # noqa: F405
        A(  # noqa: F405
            f"{label}{marker}",
            href="#",
            cls="hover:underline whitespace-nowrap",
            onclick=(
                "window.fabMonitorDashboard.closeTransientModal();"
                f"document.getElementById('sort_by').value='{key}';"
                f"document.getElementById('sort_dir').value='{next_dir}';"
            ),
            hx_get="/dashboard-shell",
            hx_target="#dashboard-shell",
            hx_swap="outerHTML",
            hx_include=DASHBOARD_FORM_INCLUDE,
            hx_sync=DASHBOARD_REQUEST_SYNC,
        )
    )


def _pager_button(label: str, *, button_id: str, delta: int, disabled: bool = False):
    classes = (ButtonT.secondary, ButtonT.sm)
    if disabled:
        return Button(label, submit=False, cls=classes, disabled=True, id=button_id)  # noqa: F405
    return Button(  # noqa: F405
        label,
        submit=False,
        id=button_id,
        cls=classes,
        onclick=f"window.fabMonitorDashboard.changeJobsPage({delta});",
    )


def _page_size_button(size: int, current: int):
    classes = (ButtonT.primary if size == current else ButtonT.ghost, ButtonT.sm)
    return Button(  # noqa: F405
        str(size),
        submit=False,
        cls=classes,
        data_page_size_button=str(size),
        onclick=f"window.fabMonitorDashboard.setJobsPageSize({size});",
    )


def _jobs_filter_form(
    rows: list[dict[str, str]],
    *,
    state_filter: str,
    rule_filter: str,
):
    state_rows = _filter_rows(rows, state_filter="", rule_filter=rule_filter)
    rule_rows = _filter_rows(rows, state_filter=state_filter, rule_filter="")

    clear_button = Button(  # noqa: F405
        "Clear",
        submit=False,
        id="jobs-clear-filters",
        hidden=not (state_filter or rule_filter),
        cls=(ButtonT.ghost, ButtonT.sm),  # noqa: F405
        onclick="window.fabMonitorDashboard.clearJobsFilters();",
    )

    def native_select(
        *,
        select_id: str,
        name: str,
        width_cls: str,
        title: str,
        options: list[tuple[str, str, bool]],
    ):
        options_html = "".join(
            f'<option value="{escape(value, quote=True)}"'
            f'{" selected" if selected else ""}>{escape(label)}</option>'
            for value, label, selected in options
        )
        attrs = {
            "id": select_id,
            "name": name,
            "class": f"uk-select uk-form-small {width_cls}",
            "title": title,
            "onfocus": "window.fabMonitorDashboard.pausePolling(12000);",
            "onpointerdown": "window.fabMonitorDashboard.pausePolling(12000);",
            "onchange": (
                "window.fabMonitorDashboard.pausePolling(2500);"
                "window.fabMonitorDashboard.closeTransientModal();"
                "document.getElementById('page').value='1';"
                "window.fabMonitorDashboard.applyJobsFilters();"
            ),
        }
        attr_html = " ".join(
            f'{key}="{escape(value, quote=True)}"'
            for key, value in attrs.items()
        )
        return NotStr(f"<select {attr_html}>{options_html}</select>")  # noqa: F405

    state_select = native_select(
        select_id="state_filter",
        name="state_filter",
        width_cls="w-36",
        title="Filter by state",
        options=[("", "All states", not state_filter)]
        + [(state, state, state == state_filter) for state in _available_states(state_rows)],
    )
    rule_select = native_select(
        select_id="rule_filter",
        name="rule_filter",
        width_cls="min-w-44 max-w-72",
        title="Filter by rule",
        options=[("", "All rules", not rule_filter)]
        + [(rule, rule, rule == rule_filter) for rule in _available_rules(rule_rows)],
    )

    return Form(  # noqa: F405
        state_select,
        rule_select,
        clear_button,
        id="jobs-filter-form",
        cls="flex items-center gap-2 flex-nowrap overflow-x-auto",
    )


def _jobs_table(
    run,
    limit: int,
    page: int,
    sort_by: str,
    sort_dir: str,
    state_filter: str,
    rule_filter: str,
):
    if run is None:
        return Div()  # noqa: F405

    all_rows = _sorted_rows(run, sort_by, sort_dir)
    filtered_rows = _filter_rows(
        all_rows,
        state_filter=state_filter,
        rule_filter=rule_filter,
    )
    filter_form = _jobs_filter_form(
        all_rows,
        state_filter=state_filter,
        rule_filter=rule_filter,
    )
    if not all_rows:
        return Card(  # noqa: F405
            P("No jobs were parsed from the selected controller log.", cls="text-sm"),
            header=DivFullySpaced(  # noqa: F405
                H3("Jobs", cls="text-base font-semibold"),  # noqa: F405
                filter_form,
            ),
            body_cls="space-y-0 py-2",
        )

    total_rows = len(filtered_rows)
    page_count = max(1, (total_rows + limit - 1) // limit)
    current_page = min(max(1, page), page_count)
    start = (current_page - 1) * limit
    stop = start + limit
    visible_jobids = {
        row["smk_jobid"]
        for row in filtered_rows[start:stop]
    }

    body = []
    for row in all_rows:
        body.append(
            Tr(  # noqa: F405
                Td(_state_badge(row["state"]), cls="whitespace-nowrap"),  # noqa: F405
                Td(row["rule"], cls="font-medium"),
                Td(row["wildcards"], cls="max-w-xl text-sm"),
                Td(row["elapsed"], cls="whitespace-nowrap text-sm"),
                Td(row["node_or_reason"], cls="max-w-sm text-sm"),
                Td(_log_link(row), cls="max-w-sm"),
                data_job_row="true",
                data_base_state=_row_base_state(row),
                data_rule=row["rule"],
                data_smk_jobid=row["smk_jobid"],
                hidden=row["smk_jobid"] not in visible_jobids,
            )
        )

    table_body = (
        Div(
            Table(  # noqa: F405
                Thead(  # noqa: F405
                    Tr(  # noqa: F405
                        _sort_header("State", "state", sort_by, sort_dir),
                        _sort_header("Rule", "rule", sort_by, sort_dir),
                        _sort_header("Wildcards", "wildcards", sort_by, sort_dir),
                        _sort_header("Elapsed", "elapsed", sort_by, sort_dir),
                        _sort_header("Node / Reason", "node_or_reason", sort_by, sort_dir),
                        _sort_header("Log", "log_file", sort_by, sort_dir),
                    )
                ),
                Tbody(*body),  # noqa: F405
            ),
            id="jobs-table-wrapper",
            hidden=not filtered_rows,
            cls="overflow-x-auto",
        )
    )

    return Card(  # noqa: F405
        table_body,
        P("No jobs match the current filters.", id="jobs-empty-state", cls="text-sm", hidden=bool(filtered_rows)),
        footer=DivFullySpaced(  # noqa: F405
            Div(
                _pager_button("Previous", button_id="jobs-prev-page", delta=-1, disabled=current_page <= 1),
                P(
                    f"Page {current_page} of {page_count} | {total_rows} jobs",
                    cls="text-xs opacity-70",
                    id="jobs-page-status",
                ),
                _pager_button("Next", button_id="jobs-next-page", delta=1, disabled=current_page >= page_count),
                cls="flex items-center gap-2 flex-wrap",
            ),
            Div(
                P("Rows per page", cls="text-xs opacity-70"),
                *[_page_size_button(size, limit) for size in PAGE_SIZE_OPTIONS],
                cls="flex items-center gap-2 flex-wrap",
            ),
        ),
        header=DivFullySpaced(  # noqa: F405
            Div(
                H3("Jobs", cls="text-base font-semibold"),  # noqa: F405
                P(
                    (
                        f"Showing {start + 1}-{min(stop, total_rows)} of {total_rows}"
                        if filtered_rows
                        else "Showing 0 jobs"
                    )
                    + (f" | {len(all_rows)} total" if total_rows != len(all_rows) else ""),
                    cls="text-xs opacity-70",
                    id="jobs-summary-text",
                ),
                cls="space-y-0.5",
            ),
            filter_form,
        ),
        id="jobs-card",
        body_cls="space-y-0 py-2",
    )


def _log_viewer_height_style(content: str) -> str:
    line_count = max(1, len(content.splitlines()))
    estimated_rem = max(12.0, min(42.0, 8.0 + (line_count * 0.22)))
    return f"height: min(72vh, {estimated_rem:.1f}rem); min-height: 12rem; max-height: 72vh;"


def _log_source_kind(path: Path) -> str:
    return "SLURM log" if "/log/slurm/" in str(path).replace("\\", "/") else "Rule log"


def _build_log_sources(run, row: dict[str, str], job) -> list[dict[str, str]]:
    sources = []
    for title, content in controller_job_log_context(
        run.controller_log,
        job,
        max_lines=200,
        max_chars=40000,
    ):
        sources.append(
            {
                "title": title,
                "subtitle": str(run.controller_log),
                "content": content,
                "kind": "Controller",
                "tone_cls": "border border-sky-200 bg-sky-50",
            }
        )

    for path in resolve_log_paths(PROJECT_DIR, row["log"]):
        sources.append(
            {
                "title": path.name,
                "subtitle": str(path),
                "content": read_log_tail(path, max_lines=250, max_chars=40000),
                "kind": _log_source_kind(path),
                "tone_cls": "",
            }
        )

    return sources


def _tab_label(title: str, seen_labels: dict[str, int]) -> str:
    count = seen_labels.get(title, 0) + 1
    seen_labels[title] = count
    if count == 1:
        return title
    return f"{title} {count}"


def _log_browser(log_sources: list[dict[str, str]], *, smk_jobid: str):
    if not log_sources:
        return Card(  # noqa: F405
            P("No readable log files for this job.", cls="text-sm"),
            cls="border border-amber-300 bg-amber-50",
        )

    root_id = f"job-log-tabs-{smk_jobid}"
    label_counts: dict[str, int] = {}
    buttons = []
    panels = []

    for index, source in enumerate(log_sources, start=1):
        tab_id = f"{root_id}-tab-{index}"
        active = index == 1
        label = _tab_label(source["title"], label_counts)
        line_count = len(source["content"].splitlines())
        buttons.append(
            Button(  # noqa: F405
                Div(
                    P(label, cls="text-sm font-medium leading-tight"),
                    P(
                        f"{source['kind']} | {line_count} lines",
                        cls="text-xs opacity-70",
                    ),
                    cls="space-y-1",
                ),
                submit=False,
                cls=f"fab-monitor-log-tab{' is-active' if active else ''}",
                role="tab",
                aria_selected="true" if active else "false",
                title=source["subtitle"],
                data_log_tab="true",
                data_tab_id=tab_id,
                onclick=f"window.fabMonitorDashboard.selectLogTab('{root_id}', '{tab_id}');",
            )
        )
        panels.append(
            Div(  # noqa: F405
                Card(
                    Pre(  # noqa: F405
                        source["content"],
                        cls="fab-monitor-log-pre",
                        style=_log_viewer_height_style(source["content"]),
                    ),
                    header=DivFullySpaced(  # noqa: F405
                        Div(
                            P(source["title"], cls="text-sm font-medium"),
                            P(source["subtitle"], cls="text-xs opacity-70 break-all"),
                            cls="space-y-1",
                        ),
                        P(f"{line_count} lines", cls="text-xs opacity-70 whitespace-nowrap"),
                    ),
                    body_cls="space-y-3",
                    cls=source["tone_cls"],
                ),
                cls=f"fab-monitor-log-panel{' is-active' if active else ''}",
                data_log_panel="true",
                data_tab_id=tab_id,
                hidden=not active,
            )
        )

    return Div(  # noqa: F405
        Div(
            Div(
                P("Log Sources", cls="text-xs font-medium uppercase tracking-wide opacity-60"),
                P(
                    f"{len(log_sources)} sources available. Select a tab to switch logs.",
                    cls="text-sm opacity-75",
                ),
                cls="space-y-1",
            ),
            Div(*buttons, cls="fab-monitor-log-tabs"),
            cls="sticky top-0 z-10 space-y-3 border-b border-slate-200 bg-white/95 pb-3 backdrop-blur",
        ),
        Div(*panels, cls="space-y-3"),
        id=root_id,
        cls="space-y-3",
    )


def _dashboard_shell(
    controller_log: str,
    user: str,
    remote_host: str,
    limit: int,
    page: int,
    refresh: int,
    sort_by: str,
    sort_dir: str,
    state_filter: str,
    rule_filter: str,
):
    run, queue_warning, selected_log = _load_run(controller_log, user, remote_host)
    poll_attrs = {}
    if refresh > 0 and (run is None or not run.is_terminal()):
        poll_attrs = {
            "hx_get": "/dashboard-shell",
            "hx_trigger": f"every {refresh}s",
            "hx_target": "this",
            "hx_swap": "outerHTML",
            "hx_include": DASHBOARD_FORM_INCLUDE,
            "hx_sync": DASHBOARD_REQUEST_SYNC,
            "hx-on::before-request": "if (!window.fabMonitorDashboard.shouldPoll()) event.preventDefault();",
        }

    return Div(  # noqa: F405
        _summary(run, queue_warning, selected_log),
        _jobs_table(run, limit, page, sort_by, sort_dir, state_filter, rule_filter),
        id="dashboard-shell",
        cls="space-y-2",
        **poll_attrs,
    )


@rt("/")  # noqa: F405
def index(
    controller_log: str = "",
    user: str = DEFAULT_USER,
    remote_host: str = DEFAULT_REMOTE_HOST,
    limit: int = DEFAULT_LIMIT,
    page: int = 1,
    refresh: int = DEFAULT_REFRESH,
    sort_by: str = "",
    sort_dir: str = "asc",
    state_filter: str = "",
    rule_filter: str = "",
):
    limit = _normalize_limit(limit)
    page = _normalize_page(page)
    refresh = _normalize_refresh(refresh)
    controller_log = _selected_log_path(controller_log)
    state_filter = _normalize_state_filter(state_filter)
    rule_filter = _normalize_rule_filter(rule_filter)

    return Title(f"Snakemake Monitor | {PROJECT_DIR.name}"), Container(  # noqa: F405
        Style(LOG_BROWSER_CSS),  # noqa: F405
        Script(LOG_BROWSER_SCRIPT),  # noqa: F405
        _controls(controller_log, user, remote_host, limit, page, refresh, sort_by, sort_dir),
        _dashboard_shell(
            controller_log,
            user,
            remote_host,
            limit,
            page,
            refresh,
            sort_by,
            sort_dir,
            state_filter,
            rule_filter,
        ),
        Div(id="modal-container"),  # noqa: F405
        Div(_job_modal_loading_template(), id="job-modal-loading-template", hidden=True),  # noqa: F405
        cls="space-y-3 py-3",
    )


@rt("/controls")  # noqa: F405
def controls(
    controller_log: str = "",
    user: str = DEFAULT_USER,
    remote_host: str = DEFAULT_REMOTE_HOST,
    limit: int = DEFAULT_LIMIT,
    page: int = 1,
    refresh: int = DEFAULT_REFRESH,
    sort_by: str = "",
    sort_dir: str = "asc",
    state_filter: str = "",
    rule_filter: str = "",
    controls_open: int = 0,
):
    return _controls(
        _selected_log_path(controller_log),
        user,
        remote_host,
        _normalize_limit(limit),
        _normalize_page(page),
        _normalize_refresh(refresh),
        sort_by,
        sort_dir,
        bool(int(controls_open)),
    )


@rt("/dashboard-shell")  # noqa: F405
def dashboard_shell(
    controller_log: str = "",
    user: str = DEFAULT_USER,
    remote_host: str = DEFAULT_REMOTE_HOST,
    limit: int = DEFAULT_LIMIT,
    page: int = 1,
    refresh: int = DEFAULT_REFRESH,
    sort_by: str = "",
    sort_dir: str = "asc",
    state_filter: str = "",
    rule_filter: str = "",
):
    return _dashboard_shell(
        _selected_log_path(controller_log),
        user,
        remote_host,
        _normalize_limit(limit),
        _normalize_page(page),
        _normalize_refresh(refresh),
        sort_by,
        sort_dir,
        _normalize_state_filter(state_filter),
        _normalize_rule_filter(rule_filter),
    )


@rt("/job-modal")  # noqa: F405
def job_modal(
    smk_jobid: str,
    controller_log: str = "",
    user: str = DEFAULT_USER,
    remote_host: str = DEFAULT_REMOTE_HOST,
):
    run, _, _ = _load_run(
        controller_log,
        user,
        remote_host,
        smk_jobid=smk_jobid,
    )
    if run is None:
        return Div()  # noqa: F405

    row_lookup = {row["smk_jobid"]: row for row in job_rows(run, limit=0)}
    row = row_lookup.get(smk_jobid)
    job = run.jobs.get(smk_jobid)
    if row is None or job is None:
        return Div()  # noqa: F405

    header = Div(
        H3(row["rule"], cls="text-lg font-semibold"),  # noqa: F405
        Div(
            CodeSpan(row["state"]),  # noqa: F405
            P(f"SMK {row['smk_jobid']} | SLURM {row['slurm_jobid']}", cls="text-sm opacity-70"),
            cls="space-y-1",
        ),
        P(row["wildcards"], cls="text-sm break-all"),
        cls="space-y-2",
    )

    return Modal(  # noqa: F405
        _log_browser(_build_log_sources(run, row, job), smk_jobid=smk_jobid),
        header=header,
        footer=Div(
            ModalCloseButton("Close", cls=(ButtonT.secondary, "static"), submit=False),  # noqa: F405
            cls="flex justify-end",
        ),
        id="job-log-modal",
        hx_init=True,
        hx_open=True,
        dialog_cls="uk-width-5-6@xl uk-width-11-12 max-h-[92vh]",
        body_cls="space-y-4 overflow-y-auto",
        style="padding-top: 2vh; padding-bottom: 2vh;",
    )


serve(host=DEFAULT_HOST, port=DEFAULT_PORT, reload=False)  # noqa: F405
