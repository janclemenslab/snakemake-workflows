"""FastHTML dashboard for Snakemake controller logs."""

from __future__ import annotations

from functools import lru_cache
from html import escape
import json
import os
from pathlib import Path
import re
import signal
import subprocess
import sys
import threading
from urllib.parse import urlencode
import webbrowser

if __package__ in {None, ""}:
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from fasthtml.common import *  # noqa: F403
from monsterui.all import *  # noqa: F403

from fab_commands.monitor import (
    apply_squeue_state,
    cancel_slurm_run,
    controller_job_log_context,
    controller_logs,
    fetch_squeue_output,
    job_rows,
    parse_controller_log,
    read_log_tail,
    resolve_log_paths,
    rule_graph_data,
)
from fab_commands.submission import (
    launch_dashboard_submission,
    normalize_submission_targets,
    submission_metadata,
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

.fab-monitor-dag-shell {
    border: 1px solid #cbd5e1;
    border-radius: 1rem;
    background: linear-gradient(180deg, #f8fafc 0%, #f1f5f9 100%);
    padding: 0.9rem;
    overflow: auto;
}

.fab-monitor-dag-legend {
    display: flex;
    flex-wrap: wrap;
    gap: 0.75rem;
}

.fab-monitor-dag-legend-item {
    display: inline-flex;
    align-items: center;
    gap: 0.4rem;
    font-size: 0.8rem;
    color: #475569;
}

.fab-monitor-dag-swatch {
    width: 0.75rem;
    height: 0.75rem;
    border-radius: 999px;
    border: 1px solid transparent;
    flex: none;
}

.fab-submit-section-grid {
    display: grid;
    gap: 0.75rem;
    grid-template-columns: minmax(15rem, 19rem) minmax(0, 1fr);
}

.fab-submit-column {
    display: flex;
    flex-direction: column;
    gap: 0.75rem;
}

.fab-submit-checklist {
    display: grid;
    gap: 0.5rem;
}

.fab-submit-rule-option,
.fab-submit-experiment-option {
    display: flex;
    align-items: flex-start;
    gap: 0.65rem;
    border: 1px solid #cbd5e1;
    border-radius: 0.85rem;
    padding: 0.7rem 0.8rem;
    background: white;
    cursor: pointer;
}

.fab-submit-rule-option > input[type="checkbox"],
.fab-submit-experiment-option > input[type="checkbox"] {
    width: 1rem;
    min-width: 1rem;
    height: 1rem;
    margin: 0.15rem 0 0;
    padding: 0;
    flex: none;
    display: block;
    appearance: auto;
    -webkit-appearance: checkbox;
    accent-color: #0f766e;
    background: transparent;
    border: 0;
    box-shadow: none;
}

.fab-submit-rule-option:focus-visible,
.fab-submit-experiment-option:focus-visible {
    outline: 2px solid #0f766e;
    outline-offset: 2px;
}

.fab-submit-rule-option.is-selected,
.fab-submit-experiment-option.is-selected {
    border-color: #0f766e;
    background: #ecfeff;
}

.fab-submit-rule-option {
    min-height: 100%;
}

.fab-submit-experiment-option {
    justify-content: space-between;
}

.fab-submit-checklist-shell {
    border: 1px solid #cbd5e1;
    border-radius: 0.85rem;
    background: #f8fafc;
    padding: 0.65rem;
}

.fab-submit-scroll {
    max-height: 34rem;
    overflow-y: auto;
    padding-right: 0.2rem;
}

.fab-submit-preview {
    border: 1px solid #cbd5e1;
    border-radius: 0.85rem;
    background: #f8fafc;
    padding: 0.8rem;
}

.fab-submit-preview-pre {
    max-height: 14rem;
    overflow: auto;
    white-space: pre-wrap;
    font-family: ui-monospace, SFMono-Regular, Menlo, monospace;
    font-size: 12px;
    line-height: 1.45;
}

.fab-submit-confirm-window {
    width: min(52rem, 96vw);
    border: 1px solid #cbd5e1;
    border-radius: 1rem;
    background: white;
    box-shadow: 0 24px 80px rgba(15, 23, 42, 0.22);
    padding: 1rem 1.25rem;
}

.fab-submit-confirm-grid {
    display: grid;
    gap: 0.75rem;
    grid-template-columns: repeat(auto-fit, minmax(15rem, 1fr));
}

.fab-submit-confirm-shell {
    border: 1px solid #cbd5e1;
    border-radius: 0.85rem;
    background: #f8fafc;
    padding: 0.8rem;
}

.fab-submit-confirm-list {
    max-height: 14rem;
    overflow: auto;
    display: grid;
    gap: 0.35rem;
}

.fab-submit-confirm-item {
    font-size: 0.82rem;
    line-height: 1.45;
    color: #334155;
    word-break: break-word;
}

.fab-submit-confirm-item.is-mono {
    font-family: ui-monospace, SFMono-Regular, Menlo, monospace;
    font-size: 0.75rem;
}

.fab-submit-status-badge {
    display: inline-flex;
    align-items: center;
    justify-content: center;
    min-width: 4.5rem;
    padding: 0.2rem 0.55rem;
    border-radius: 999px;
    border: 1px solid #cbd5e1;
    background: #f8fafc;
    font-size: 0.72rem;
    font-weight: 600;
    text-transform: uppercase;
    letter-spacing: 0.02em;
}

.fab-submit-status-badge.is-finished {
    border-color: #86efac;
    background: #ecfdf5;
    color: #047857;
}

.fab-submit-status-badge.is-new {
    border-color: #cbd5e1;
    background: #f8fafc;
    color: #475569;
}

.fab-page-tabs {
    display: inline-flex;
    gap: 0.5rem;
}

.fab-page-tab {
    display: inline-flex;
    align-items: center;
    justify-content: center;
    min-width: 6.5rem;
    padding: 0.45rem 0.9rem;
    border: 1px solid #cbd5e1;
    border-radius: 999px;
    background: white;
    color: #334155;
    text-decoration: none;
    font-size: 0.9rem;
    font-weight: 600;
}

.fab-page-tab.is-active {
    border-color: #0f766e;
    background: #ecfeff;
    color: #115e59;
}

@media (max-width: 960px) {
    .fab-submit-section-grid {
        grid-template-columns: 1fr;
    }
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
  showLoading(templateId = 'job-modal-loading-template') {
    const container = document.getElementById('modal-container');
    const template = document.getElementById(templateId);
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

    let limit = parseInt(limitInput ? (limitInput.value || '100') : '100', 10);
    if (!Number.isFinite(limit) || limit < 1) {
      limit = 100;
    }

    let page = parseInt(pageInput ? (pageInput.value || '1') : '1', 10);
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
  controlsOpen() {
    const input = document.getElementById('controls_open');
    return Boolean(input && input.value === '1');
  },
  setControlsOpen(open) {
    const input = document.getElementById('controls_open');
    const panel = document.getElementById('controls-options-panel');
    const button = document.getElementById('controls-toggle-button');
    if (input) {
      input.value = open ? '1' : '0';
    }
    if (panel) {
      panel.style.display = open ? 'flex' : 'none';
      panel.setAttribute('aria-hidden', open ? 'false' : 'true');
    }
    if (button) {
      button.setAttribute('aria-expanded', open ? 'true' : 'false');
    }
  },
  toggleControlsOpen(event) {
    if (event) {
      event.preventDefault();
      event.stopPropagation();
    }
    this.pausePolling(2500);
    this.setControlsOpen(!this.controlsOpen());
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
  },
  submissionMetadata() {
    const metadataField = document.getElementById('submission-metadata');
    if (!metadataField) return null;
    const rawValue = metadataField.value || '';
    if (this._submissionMetadataRaw !== rawValue) {
      this._submissionMetadataRaw = rawValue;
      try {
        this._submissionMetadata = JSON.parse(rawValue);
      } catch (error) {
        this._submissionMetadata = null;
      }
    }
    return this._submissionMetadata;
  },
  escapeHtml(value) {
    return String(value == null ? '' : value).replace(/[&<>"']/g, (char) => (
      {
        '&': '&amp;',
        '<': '&lt;',
        '>': '&gt;',
        '"': '&quot;',
        "'": '&#39;'
      }[char] || char
    ));
  },
  submissionControlsSnapshot() {
    const values = {};
    const form = document.getElementById('submission-controls-form');
    if (!form) return values;
    new FormData(form).forEach((value, key) => {
      values[key] = String(value);
    });
    return values;
  },
  submissionSelectionSnapshot() {
    const form = document.getElementById('submission-form');
    if (!form) return null;

    const metadata = this.submissionMetadata() || {};
    const rulesByName = {};
    const experimentLookup = {};
    (metadata.rules || []).forEach((rule) => {
      rulesByName[rule.name] = rule;
    });
    (metadata.experiments || []).forEach((experiment) => {
      experimentLookup[experiment.id] = experiment;
    });

    const selectedRules = Array.from(
      form.querySelectorAll('[data-submit-rule="true"]:checked')
    ).map((input) => input.value);
    const selectedExperiments = Array.from(
      form.querySelectorAll('[data-submit-experiment="true"]:checked')
    ).map((input) => input.value);

    const selectedRuleItems = selectedRules.map((ruleName) => (
      rulesByName[ruleName] || { name: ruleName, label: ruleName }
    ));
    const selectedExperimentItems = selectedExperiments
      .map((experimentId) => experimentLookup[experimentId])
      .filter(Boolean);

    const targets = [];
    const seenTargets = new Set();
    selectedExperimentItems.forEach((experiment) => {
      selectedRules.forEach((ruleName) => {
        const ruleTargets = (experiment.targets && experiment.targets[ruleName]) || [];
        ruleTargets.forEach((target) => {
          if (!target || seenTargets.has(target)) return;
          seenTargets.add(target);
          targets.push(target);
        });
      });
    });

    const summary = document.getElementById('selection_summary');
    const selectionSummary = summary ? (summary.value || '') : '';

    return {
      selectedRules,
      selectedExperiments,
      selectedRuleItems,
      selectedExperimentItems,
      targets,
      selectionSummary
    };
  },
  experimentStatus(ruleStatus, selectedRules) {
    const effectiveRules = selectedRules.length
      ? selectedRules
      : Object.keys(ruleStatus || {});
    if (!effectiveRules.length) {
      return { state: 'new', done: 0, total: 0 };
    }
    const done = effectiveRules.filter((ruleName) => Boolean(ruleStatus && ruleStatus[ruleName])).length;
    return {
      state: done === effectiveRules.length ? 'finished' : 'new',
      done,
      total: effectiveRules.length
    };
  },
  syncCheckboxCard(checkbox) {
    if (!checkbox) return;
    const card = checkbox.closest('[data-checkbox-card="true"]');
    if (!card) return;
    card.classList.toggle('is-selected', checkbox.checked);
  },
  toggleCheckboxCard(inputId, event) {
    const target = event ? event.target : null;
    if (
      target &&
      typeof target.closest === 'function' &&
      target.closest('input, button, a, select, textarea')
    ) {
      return;
    }
    const checkbox = document.getElementById(inputId);
    if (!checkbox || checkbox.disabled) return;
    checkbox.checked = !checkbox.checked;
    checkbox.dispatchEvent(new Event('change', { bubbles: true }));
  },
  setSubmissionGroupChecked(selector, checked) {
    const form = document.getElementById('submission-form');
    if (!form) return;
    Array.from(form.querySelectorAll(selector)).forEach((checkbox) => {
      if (checkbox.disabled) return;
      checkbox.checked = checked;
      this.syncCheckboxCard(checkbox);
    });
    this.updateSubmissionPreview();
  },
  selectNewExperiments() {
    const form = document.getElementById('submission-form');
    if (!form) return;

    const metadata = this.submissionMetadata() || {};
    const experimentLookup = {};
    (metadata.experiments || []).forEach((experiment) => {
      experimentLookup[experiment.id] = experiment;
    });

    const selectedRules = Array.from(
      form.querySelectorAll('[data-submit-rule="true"]:checked')
    ).map((input) => input.value);

    Array.from(form.querySelectorAll('[data-submit-experiment="true"]')).forEach((checkbox) => {
      if (checkbox.disabled) return;
      const experiment = experimentLookup[checkbox.value];
      const status = this.experimentStatus(
        experiment ? (experiment.rule_status || {}) : {},
        selectedRules
      );
      checkbox.checked = status.state === 'new';
      this.syncCheckboxCard(checkbox);
    });

    this.updateSubmissionPreview();
  },
  renderSubmissionConfirmList(title, items, { mono = false, emptyText = 'None selected', limit = 12 } = {}) {
    const visibleItems = items.slice(0, limit);
    const extraCount = Math.max(0, items.length - visibleItems.length);
    const itemCls = mono ? 'fab-submit-confirm-item is-mono' : 'fab-submit-confirm-item';
    const listItems = visibleItems.length
      ? visibleItems.map((item) => (
          `<li class="${itemCls}">${this.escapeHtml(item)}</li>`
        )).join('')
      : `<li class="${itemCls} opacity-70">${this.escapeHtml(emptyText)}</li>`;
    const moreLine = extraCount
      ? `<li class="${itemCls} opacity-70">... (+${extraCount} more)</li>`
      : '';
    return `
      <div class="fab-submit-confirm-shell space-y-2">
        <p class="text-sm font-medium">${this.escapeHtml(title)}</p>
        <ul class="fab-submit-confirm-list">${listItems}${moreLine}</ul>
      </div>
    `;
  },
  submissionRequestFormData(mode) {
    const formData = new FormData();
    const controlsForm = document.getElementById('submission-controls-form');
    if (controlsForm) {
      new FormData(controlsForm).forEach((value, key) => {
        formData.append(key, value);
      });
    }
    if (mode === 'all') {
      formData.set('submit_mode', 'all');
      formData.set('targets_text', '');
      formData.set('selection_summary', 'Everything pending');
      return formData;
    }

    const submissionForm = document.getElementById('submission-form');
    if (submissionForm) {
      new FormData(submissionForm).forEach((value, key) => {
        formData.append(key, value);
      });
    }
    return formData;
  },
  setSubmissionStatusHtml(html) {
    const status = document.getElementById('submission-status');
    if (!status) return;
    status.innerHTML = html;
    if (typeof status.scrollIntoView === 'function') {
      status.scrollIntoView({ behavior: 'smooth', block: 'nearest' });
    }
  },
  submissionPendingStatus(message) {
    return `
      <div class="rounded-xl border border-slate-300 bg-slate-50 px-4 py-3 text-sm text-slate-700">
        ${this.escapeHtml(message)}
      </div>
    `;
  },
  submissionErrorStatus(message) {
    return `
      <div class="rounded-xl border border-red-300 bg-red-50 px-4 py-3 text-sm text-red-700">
        ${this.escapeHtml(message)}
      </div>
    `;
  },
  renderSubmissionModal({ title, subtitle = '', bodyHtml = '', footerHtml = '', allowClose = true }) {
    const container = document.getElementById('modal-container');
    if (!container) return;

    const overlayAttrs = allowClose
      ? 'onclick="window.fabMonitorDashboard.closeTransientModal();"'
      : '';
    const closeButton = allowClose
      ? `
        <button
          type="button"
          class="rounded-lg border border-slate-300 px-3 py-2 text-sm text-slate-700"
          onclick="window.fabMonitorDashboard.closeTransientModal();"
        >
          Close
        </button>
      `
      : '';

    container.innerHTML = `
      <div class="fab-monitor-loading-overlay" ${overlayAttrs}>
        <div class="fab-submit-confirm-window space-y-4" onclick="event.stopPropagation();">
          <div class="flex items-start justify-between gap-3">
            <div class="space-y-1">
              <p class="text-base font-semibold">${this.escapeHtml(title)}</p>
              <p class="text-sm opacity-70">${this.escapeHtml(subtitle)}</p>
            </div>
            ${closeButton}
          </div>
          ${bodyHtml}
          ${footerHtml}
        </div>
      </div>
    `;
  },
  showSubmissionConfirm(mode) {
    this.updateSubmissionPreview();

    const controls = this.submissionControlsSnapshot();
    const user = controls.user || '(not set)';
    const remoteHost = controls.remote_host || '(not set)';

    let title = 'Confirm Submission';
    let subtitle = '';
    let confirmLabel = 'Submit';
    let bodyHtml = '';

    if (mode === 'selection') {
      const snapshot = this.submissionSelectionSnapshot();
      if (!snapshot || !snapshot.targets.length) {
        this.updateSubmissionPreview();
        return;
      }

      const ruleItems = snapshot.selectedRuleItems.map((rule) => rule.label || rule.name || '');
      const experimentItems = snapshot.selectedExperimentItems.map((experiment) => (
        `${experiment.group} | ${experiment.label}`
      ));
      title = 'Confirm Selected Submission';
      subtitle = snapshot.selectionSummary || `${ruleItems.length} stages | ${experimentItems.length} experiments`;
      confirmLabel = `Submit ${snapshot.targets.length} target${snapshot.targets.length === 1 ? '' : 's'}`;
      bodyHtml = `
        <div class="fab-submit-confirm-grid">
          ${this.renderSubmissionConfirmList('Connection', [`User: ${user}`, `Host: ${remoteHost}`], { limit: 2 })}
          ${this.renderSubmissionConfirmList('Rules', ruleItems, { limit: 8 })}
        </div>
        <div class="fab-submit-confirm-grid">
          ${this.renderSubmissionConfirmList('Experiments', experimentItems, { limit: 12 })}
          ${this.renderSubmissionConfirmList('Targets', snapshot.targets, { mono: true, limit: 24 })}
        </div>
      `;
    } else {
      title = 'Confirm Full Submission';
      subtitle = 'Launch the workflow default target and submit everything pending.';
      confirmLabel = 'Submit Everything Pending';
      bodyHtml = `
        <div class="fab-submit-confirm-grid">
          ${this.renderSubmissionConfirmList('Connection', [`User: ${user}`, `Host: ${remoteHost}`], { limit: 2 })}
          ${this.renderSubmissionConfirmList(
            'Scope',
            ['This uses the workflow default target, equivalent to the current fab submit workflow.'],
            { limit: 1 }
          )}
        </div>
        <div class="fab-submit-confirm-shell space-y-2">
          <p class="text-sm font-medium">Note</p>
          <p class="fab-submit-confirm-item">
            Exact pending outputs are not enumerated here. Use Selective Submit when you want a file-level target preview before launch.
          </p>
        </div>
      `;
    }

    this.renderSubmissionModal({
      title,
      subtitle,
      bodyHtml,
      footerHtml: `
        <div class="flex justify-end gap-2">
          <button
            type="button"
            class="rounded-lg border border-slate-300 px-3 py-2 text-sm text-slate-700"
            onclick="window.fabMonitorDashboard.closeTransientModal();"
          >
            Cancel
          </button>
          <button
            type="button"
            class="rounded-lg bg-teal-700 px-3 py-2 text-sm font-semibold text-white"
            onclick="window.fabMonitorDashboard.confirmSubmission('${mode}');"
          >
            ${this.escapeHtml(confirmLabel)}
          </button>
        </div>
      `
    });
  },
  async confirmSubmission(mode) {
    this.pausePolling(2500);
    this.renderSubmissionModal({
      title: mode === 'all' ? 'Submitting Everything Pending' : 'Submitting Selected Targets',
      subtitle: 'Waiting for the controller to start and report its log path.',
      bodyHtml: `
        <div class="fab-submit-confirm-shell">
          <div class="flex items-center gap-3" role="status" aria-live="polite">
            <span class="fab-monitor-spinner" aria-hidden="true"></span>
            <div class="space-y-1">
              <p class="text-sm font-medium">Launching remote controller</p>
              <p class="fab-submit-confirm-item">
                The modal will update when the submission response includes the controller log name.
              </p>
            </div>
          </div>
        </div>
      `,
      allowClose: false
    });

    try {
      const response = await fetch('/submit-launch', {
        method: 'POST',
        body: this.submissionRequestFormData(mode),
        credentials: 'same-origin'
      });
      const html = await response.text();
      this.setSubmissionStatusHtml(html);
      this.renderSubmissionModal({
        title: mode === 'all' ? 'Full Submission Result' : 'Selected Submission Result',
        subtitle: 'The controller launch response is ready.',
        bodyHtml: `<div class="space-y-3">${html}</div>`,
        footerHtml: `
          <div class="flex justify-end">
            <button
              type="button"
              class="rounded-lg border border-slate-300 px-3 py-2 text-sm text-slate-700"
              onclick="window.fabMonitorDashboard.closeTransientModal();"
            >
              Close
            </button>
          </div>
        `
      });
    } catch (error) {
      const errorHtml = this.submissionErrorStatus(
        error && error.message
          ? error.message
          : 'The submission request failed before the server responded.'
      );
      this.setSubmissionStatusHtml(errorHtml);
      this.renderSubmissionModal({
        title: 'Submission Failed',
        subtitle: 'The server did not return a successful launch response.',
        bodyHtml: errorHtml,
        footerHtml: `
          <div class="flex justify-end">
            <button
              type="button"
              class="rounded-lg border border-slate-300 px-3 py-2 text-sm text-slate-700"
              onclick="window.fabMonitorDashboard.closeTransientModal();"
            >
              Close
            </button>
          </div>
        `
      });
    }
  },
  updateSubmissionPreview() {
    const form = document.getElementById('submission-form');
    if (!form) return;

    const metadata = this.submissionMetadata() || {};
    const experimentLookup = {};
    (metadata.experiments || []).forEach((experiment) => {
      experimentLookup[experiment.id] = experiment;
    });

    const selectedRules = Array.from(
      form.querySelectorAll('[data-submit-rule="true"]:checked')
    ).map((input) => input.value);
    const selectedExperiments = Array.from(
      form.querySelectorAll('[data-submit-experiment="true"]:checked')
    ).map((input) => input.value);

    const targets = [];
    const seenTargets = new Set();
    selectedExperiments.forEach((experimentId) => {
      const experiment = experimentLookup[experimentId];
      if (!experiment) return;
      selectedRules.forEach((ruleName) => {
        const ruleTargets = (experiment.targets && experiment.targets[ruleName]) || [];
        ruleTargets.forEach((target) => {
          if (!target || seenTargets.has(target)) return;
          seenTargets.add(target);
          targets.push(target);
        });
      });
    });

    const targetsField = document.getElementById('targets_text');
    if (targetsField) {
      targetsField.value = targets.join('\\n');
    }

    const summaryField = document.getElementById('selection_summary');
    if (summaryField) {
      const parts = [];
      if (selectedRules.length) {
        parts.push(`${selectedRules.length} stage${selectedRules.length === 1 ? '' : 's'}`);
      }
      if (selectedExperiments.length) {
        parts.push(`${selectedExperiments.length} experiment${selectedExperiments.length === 1 ? '' : 's'}`);
      }
      summaryField.value = parts.join(' | ');
    }

    const previewSummary = document.getElementById('submission-preview-summary');
    if (previewSummary) {
      if (!targets.length) {
        previewSummary.textContent = 'Select at least one rule and one experiment.';
      } else {
        previewSummary.textContent =
          `${targets.length} target${targets.length === 1 ? '' : 's'} ready for submission`;
      }
    }

    const previewTargets = document.getElementById('submission-preview-targets');
    if (previewTargets) {
      if (!targets.length) {
        previewTargets.textContent = '';
      } else {
        const visibleTargets = targets.slice(0, 24);
        const suffix = targets.length > visibleTargets.length
          ? `\\n... (+${targets.length - visibleTargets.length} more)`
          : '';
        previewTargets.textContent = visibleTargets.join('\\n') + suffix;
      }
    }

    const submitButton = document.getElementById('submit-selection-button');
    if (submitButton) {
      submitButton.disabled = targets.length === 0;
    }

    Array.from(form.querySelectorAll('[data-experiment-id]')).forEach((node) => {
      const experimentId = node.dataset.experimentId || '';
      const experiment = experimentLookup[experimentId];
      if (!experiment) return;
      const status = this.experimentStatus(experiment.rule_status || {}, selectedRules);
      const badge = node.querySelector('[data-experiment-status-badge]');
      if (badge) {
        badge.textContent = status.state;
        badge.classList.toggle('is-finished', status.state === 'finished');
        badge.classList.toggle('is-new', status.state !== 'finished');
      }
      const detail = node.querySelector('[data-experiment-status-detail]');
      if (detail) {
        if (!status.total) {
          detail.textContent = 'No rules selected';
        } else {
          detail.textContent = `${status.done}/${status.total} selected rule${status.total === 1 ? '' : 's'} done`;
        }
      }
    });
  },
  initSubmissionForm() {
    const form = document.getElementById('submission-form');
    if (!form) return;
    if (form.dataset.bound === '1') {
      this.updateSubmissionPreview();
      return;
    }

    form.dataset.bound = '1';

    Array.from(form.querySelectorAll('[data-submit-experiment="true"]')).forEach((checkbox) => {
      checkbox.addEventListener('change', () => {
        this.syncCheckboxCard(checkbox);
        this.updateSubmissionPreview();
      });
      this.syncCheckboxCard(checkbox);
    });

    Array.from(form.querySelectorAll('[data-submit-rule="true"]')).forEach((checkbox) => {
      checkbox.addEventListener('change', () => {
        this.syncCheckboxCard(checkbox);
        this.updateSubmissionPreview();
      });
      this.syncCheckboxCard(checkbox);
    });

    this.updateSubmissionPreview();
  }
};

if (!window.fabMonitorDashboardBound) {
  window.fabMonitorDashboardBound = true;
  document.addEventListener('DOMContentLoaded', () => {
    window.fabMonitorDashboard.renderJobsTable();
    window.fabMonitorDashboard.initSubmissionForm();
  });
  document.body.addEventListener('htmx:afterSwap', (event) => {
    if (event.target && event.target.id === 'dashboard-shell') {
      window.fabMonitorDashboard.renderJobsTable();
    }
    if (event.target && event.target.id === 'submission-page-shell') {
      window.fabMonitorDashboard.initSubmissionForm();
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


def _env_bool(name: str, default: bool = False) -> bool:
    value = os.environ.get(name)
    if value is None:
        return default
    return value.strip().lower() in {"1", "true", "yes", "on"}


DEFAULT_LIMIT = _env_int("FAB_MONITOR_LIMIT", 100, minimum=10)
DEFAULT_REFRESH = _env_int("FAB_MONITOR_REFRESH", 10, minimum=0)
DEFAULT_HOST = os.environ.get("FAB_MONITOR_BIND_HOST", "127.0.0.1").strip() or "127.0.0.1"
DEFAULT_PORT = _env_int("FAB_MONITOR_PORT", 2718, minimum=1)
DEFAULT_OPEN_BROWSER = _env_bool("FAB_MONITOR_OPEN_BROWSER", False)
DEFAULT_BROWSER_URL = os.environ.get("FAB_MONITOR_BROWSER_URL", "").strip()
DEFAULT_VERBOSE = _env_bool("FAB_MONITOR_VERBOSE", False)
PAGE_SIZE_OPTIONS = [50, 100, 200]
DASHBOARD_FORM_INCLUDE = "#controls-form, #jobs-filter-form"
DASHBOARD_REQUEST_SYNC = "#dashboard-shell:replace"


app, rt = fast_app(  # noqa: F405
    pico=False,
    live=False,
    title="Snakemake Monitor",
    hdrs=Theme.slate.headers(),  # noqa: F405
)


def _schedule_browser_open():
    if not DEFAULT_OPEN_BROWSER or not DEFAULT_BROWSER_URL:
        return

    def open_browser():
        try:
            webbrowser.open(DEFAULT_BROWSER_URL, new=2)
        except Exception:
            pass

    timer = threading.Timer(1.0, open_browser)
    timer.daemon = True
    timer.start()


def _schedule_server_shutdown(delay_seconds: float = 0.6):
    def shutdown():
        try:
            os.kill(os.getpid(), signal.SIGTERM)
        except Exception:
            os._exit(0)

    timer = threading.Timer(delay_seconds, shutdown)
    timer.daemon = True
    timer.start()


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


_DAG_STATE_STYLES = {
    "failed": {
        "fill": "#fef2f2",
        "stroke": "#fca5a5",
        "text": "#b91c1c",
    },
    "running": {
        "fill": "#eff6ff",
        "stroke": "#93c5fd",
        "text": "#1d4ed8",
    },
    "pending": {
        "fill": "#f8fafc",
        "stroke": "#cbd5e1",
        "text": "#475569",
    },
    "finished": {
        "fill": "#ecfdf5",
        "stroke": "#86efac",
        "text": "#047857",
    },
    "unknown": {
        "fill": "#f8fafc",
        "stroke": "#cbd5e1",
        "text": "#475569",
    },
}
_DOT_NODE_RE = re.compile(r'^\s*([\w]+)\[label\s*=\s*"([^"]+)"')
_DOT_EDGE_RE = re.compile(r"^\s*([\w]+)\s*->\s*([\w]+)")


def _dag_state_style(state: str) -> dict[str, str]:
    return _DAG_STATE_STYLES.get(state, _DAG_STATE_STYLES["unknown"])


def _dag_summary_text(node: dict[str, object]) -> str:
    status_counts = node.get("status_counts", {})
    if not isinstance(status_counts, dict):
        return f"{node.get('job_count', 0)} jobs"
    parts = [
        f"{count} {status}"
        for status, count in status_counts.items()
        if int(count)
    ]
    if not parts:
        return f"{node.get('job_count', 0)} jobs"
    return f"{node.get('job_count', 0)} jobs | " + ", ".join(parts)


def _dag_display_label(label: str, *, limit: int = 24) -> str:
    if len(label) <= limit:
        return label
    return f"{label[: limit - 1]}..."


def _workflow_snakefile(project_dir: Path) -> Path | None:
    candidates = [
        project_dir / "workflow" / "snakefile",
        project_dir / "Snakefile",
        project_dir / "snakefile",
    ]
    for candidate in candidates:
        if candidate.exists():
            return candidate
    return None


@lru_cache(maxsize=8)
def _workflow_rule_graph_cached(
    project_dir_str: str,
    snakefile_str: str,
    snakefile_mtime_ns: int,
) -> dict[str, list[dict[str, object]]]:
    del snakefile_mtime_ns
    project_dir = Path(project_dir_str)
    snakefile = Path(snakefile_str)
    snakefile_arg = (
        str(snakefile.relative_to(project_dir))
        if snakefile.is_relative_to(project_dir)
        else str(snakefile)
    )

    try:
        result = subprocess.run(
            ["snakemake", "-s", snakefile_arg, "--rulegraph"],
            cwd=project_dir,
            capture_output=True,
            text=True,
            timeout=20,
            check=False,
        )
    except (OSError, subprocess.SubprocessError):
        return {"nodes": [], "edges": []}

    if result.returncode != 0 or not result.stdout.strip():
        return {"nodes": [], "edges": []}

    labels: dict[str, str] = {}
    edge_ids: set[tuple[str, str]] = set()
    for line in result.stdout.splitlines():
        if match := _DOT_NODE_RE.match(line):
            node_id, label = match.groups()
            labels[node_id] = label
            continue
        if match := _DOT_EDGE_RE.match(line):
            edge_ids.add(match.groups())

    nodes = [
        {"id": label, "label": label}
        for label in sorted(labels.values(), key=str.casefold)
    ]

    edges = []
    seen_edges: set[tuple[str, str]] = set()
    for source_id, target_id in sorted(edge_ids):
        source = labels.get(source_id)
        target = labels.get(target_id)
        if not source or not target or (source, target) in seen_edges:
            continue
        seen_edges.add((source, target))
        edges.append({"source": source, "target": target})

    return {"nodes": nodes, "edges": edges}


def _workflow_rule_graph(project_dir: Path) -> dict[str, list[dict[str, object]]]:
    snakefile = _workflow_snakefile(project_dir)
    if snakefile is None:
        return {"nodes": [], "edges": []}
    return _workflow_rule_graph_cached(
        str(project_dir),
        str(snakefile),
        snakefile.stat().st_mtime_ns,
    )


def _dag_graph_for_run(run) -> tuple[dict[str, list[dict[str, object]]], bool]:
    observed_graph = rule_graph_data(run)
    observed_nodes = {
        str(node["id"]): node for node in observed_graph.get("nodes", [])
    }
    workflow_graph = _workflow_rule_graph(PROJECT_DIR)
    if not workflow_graph.get("nodes"):
        return observed_graph, False

    merged_nodes: list[dict[str, object]] = []
    seen: set[str] = set()
    for workflow_node in workflow_graph["nodes"]:
        node_id = str(workflow_node["id"])
        seen.add(node_id)
        observed = observed_nodes.get(node_id)
        merged_nodes.append(
            {
                "id": node_id,
                "label": node_id,
                "job_count": int(observed.get("job_count", 0)) if observed else 0,
                "state": str(observed.get("state", "unknown")) if observed else "unknown",
                "status_counts": dict(observed.get("status_counts", {})) if observed else {},
            }
        )

    for node_id in sorted(
        [node_id for node_id in observed_nodes if node_id not in seen],
        key=str.casefold,
    ):
        observed = observed_nodes[node_id]
        merged_nodes.append(
            {
                "id": node_id,
                "label": node_id,
                "job_count": int(observed.get("job_count", 0)),
                "state": str(observed.get("state", "unknown")),
                "status_counts": dict(observed.get("status_counts", {})),
            }
        )

    return {"nodes": merged_nodes, "edges": workflow_graph["edges"]}, True


def _rule_graph_svg(graph: dict[str, list[dict[str, object]]]) -> str:
    nodes = graph.get("nodes", [])
    edges = graph.get("edges", [])
    if not nodes:
        return ""

    node_lookup = {str(node["id"]): node for node in nodes}
    adjacency = {node_id: set() for node_id in node_lookup}
    parents = {node_id: set() for node_id in node_lookup}
    indegree = {node_id: 0 for node_id in node_lookup}

    for edge in edges:
        source = str(edge.get("source", ""))
        target = str(edge.get("target", ""))
        if source not in node_lookup or target not in node_lookup or source == target:
            continue
        if target in adjacency[source]:
            continue
        adjacency[source].add(target)
        parents[target].add(source)
        indegree[target] += 1

    queue = sorted(
        [node_id for node_id, degree in indegree.items() if degree == 0],
        key=str.casefold,
    )
    topo: list[str] = []
    while queue:
        current = queue.pop(0)
        topo.append(current)
        for target in sorted(adjacency[current], key=str.casefold):
            indegree[target] -= 1
            if indegree[target] == 0:
                queue.append(target)
        queue.sort(key=str.casefold)

    if len(topo) != len(node_lookup):
        remaining = sorted(
            [node_id for node_id in node_lookup if node_id not in topo],
            key=str.casefold,
        )
        topo.extend(remaining)

    levels: dict[str, int] = {}
    for node_id in topo:
        parent_levels = [levels[parent] + 1 for parent in parents[node_id] if parent in levels]
        levels[node_id] = max(parent_levels, default=0)

    column_count = max(levels.values(), default=0) + 1
    state_order = {"failed": 0, "running": 1, "pending": 2, "finished": 3, "unknown": 4}
    columns: dict[int, list[str]] = {level: [] for level in range(column_count)}
    for node_id, level in levels.items():
        columns.setdefault(level, []).append(node_id)
    for node_ids in columns.values():
        node_ids.sort(
            key=lambda node_id: (
                state_order.get(str(node_lookup[node_id].get("state", "unknown")), 99),
                node_id.casefold(),
            )
        )

    node_w = 220
    node_h = 72
    col_gap = 96
    row_gap = 34
    margin_x = 40
    margin_y = 28
    inner_height = max(
        (
            (len(node_ids) * node_h) + (max(len(node_ids) - 1, 0) * row_gap)
            for node_ids in columns.values()
        ),
        default=node_h,
    )
    width = (column_count * node_w) + (max(column_count - 1, 0) * col_gap) + (margin_x * 2)
    height = max(220, inner_height + (margin_y * 2))

    positions: dict[str, tuple[float, float]] = {}
    for level in range(column_count):
        node_ids = columns.get(level, [])
        column_height = (len(node_ids) * node_h) + (max(len(node_ids) - 1, 0) * row_gap)
        y = margin_y + ((inner_height - column_height) / 2 if inner_height > column_height else 0)
        x = margin_x + (level * (node_w + col_gap))
        for node_id in node_ids:
            positions[node_id] = (x, y)
            y += node_h + row_gap

    edge_paths = []
    for edge in edges:
        source = str(edge.get("source", ""))
        target = str(edge.get("target", ""))
        if source not in positions or target not in positions:
            continue
        source_x, source_y = positions[source]
        target_x, target_y = positions[target]
        x1 = source_x + node_w
        y1 = source_y + (node_h / 2)
        x2 = target_x
        y2 = target_y + (node_h / 2)
        dx = max(32.0, (x2 - x1) * 0.42)
        edge_paths.append(
            (
                f'<path d="M {x1:.1f} {y1:.1f} '
                f'C {x1 + dx:.1f} {y1:.1f}, {x2 - dx:.1f} {y2:.1f}, {x2:.1f} {y2:.1f}" '
                'fill="none" stroke="#94a3b8" stroke-width="2" '
                'marker-end="url(#dag-arrow)" />'
            )
        )

    node_blocks = []
    for node_id in topo:
        x, y = positions[node_id]
        node = node_lookup[node_id]
        style = _dag_state_style(str(node.get("state", "unknown")))
        label = str(node.get("label", node_id))
        meta = _dag_summary_text(node)
        tooltip = escape(f"{label}\n{meta}")
        label_text = escape(_dag_display_label(label))
        meta_text = escape(meta)
        node_blocks.append(
            f"""
<g>
  <title>{tooltip}</title>
  <rect x="{x:.1f}" y="{y:.1f}" width="{node_w}" height="{node_h}" rx="14"
        fill="{style['fill']}" stroke="{style['stroke']}" stroke-width="2" />
  <text x="{x + 14:.1f}" y="{y + 28:.1f}" fill="{style['text']}"
        font-family="ui-sans-serif, system-ui, sans-serif" font-size="14" font-weight="600">
    <tspan x="{x + 14:.1f}" dy="0">{label_text}</tspan>
    <tspan x="{x + 14:.1f}" dy="18" font-size="12" font-weight="500" fill="#475569">{meta_text}</tspan>
  </text>
</g>"""
        )

    return (
        f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {width} {height}" '
        f'width="{width}" height="{height}" role="img" aria-label="Workflow rule graph">'
        '<defs><marker id="dag-arrow" viewBox="0 0 10 10" refX="8" refY="5" '
        'markerWidth="8" markerHeight="8" orient="auto-start-reverse">'
        '<path d="M 0 0 L 10 5 L 0 10 z" fill="#94a3b8" /></marker></defs>'
        + "".join(edge_paths)
        + "".join(node_blocks)
        + "</svg>"
    )


def _dag_legend():
    items = [
        ("failed", "Failed"),
        ("running", "Running"),
        ("pending", "Pending"),
        ("finished", "Finished"),
        ("unknown", "Not materialized"),
    ]
    return Div(  # noqa: F405
        *[
            Span(  # noqa: F405
                Span(  # noqa: F405
                    cls="fab-monitor-dag-swatch",
                    style=(
                        f"background: {_dag_state_style(state)['fill']}; "
                        f"border-color: {_dag_state_style(state)['stroke']};"
                    ),
                ),
                label,
                cls="fab-monitor-dag-legend-item",
            )
            for state, label in items
        ],
        cls="fab-monitor-dag-legend",
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

    run_select = native_select(
        select_id="controller_log",
        name="controller_log",
        options=[
            (str(path), path.name, str(path) == controller_log) for path in logs
        ],
        title="Select controller log",
        class_name="uk-select uk-form-small",
        style="width: min(28rem, 100%);",
        disabled=not logs,
    )

    run_header = Div(
        H2(PROJECT_DIR.name, cls="text-lg font-semibold tracking-tight"),  # noqa: F405
        run_select,
        Button(  # noqa: F405
            UkIcon("refresh-cw", width=14, height=14),  # noqa: F405
            submit=False,
            cls=(ButtonT.ghost, ButtonT.sm),  # noqa: F405
            title="Refresh controller log list",
            hx_get="/controls",
            hx_target="#controls-card",
            hx_swap="outerHTML",
            hx_include="#controls-form",
        ),
        Button(  # noqa: F405
            "Options",
            submit=False,
            id="controls-toggle-button",
            cls=(ButtonT.ghost, ButtonT.sm, "ml-auto"),  # noqa: F405
            aria_expanded="true" if controls_open else "false",
            onclick="window.fabMonitorDashboard.toggleControlsOpen(event);",
        ),
        cls="flex flex-wrap items-center gap-2 min-w-0",
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

    return Card(  # noqa: F405
        Form(  # noqa: F405
            Input(type="hidden", id="page", name="page", value=str(page)),  # noqa: F405
            Input(type="hidden", id="sort_by", name="sort_by", value=sort_by),  # noqa: F405
            Input(type="hidden", id="sort_dir", name="sort_dir", value=sort_dir),  # noqa: F405
            Input(type="hidden", id="limit", name="limit", value=str(limit)),  # noqa: F405
            Input(type="hidden", id="controls_open", name="controls_open", value="1" if controls_open else "0"),  # noqa: F405
            run_header,
            Div(
                user_input,
                host_input,
                refresh_select,
                id="controls-options-panel",
                aria_hidden="false" if controls_open else "true",
                style="display: flex;" if controls_open else "display: none;",
                cls="flex flex-wrap xl:flex-nowrap gap-2 items-end pt-3",
            ),
            id="controls-form",
            cls="space-y-2",
        ),
        id="controls-card",
        body_cls="space-y-0 py-3",
    )


def _query_url(path: str, **params: str):
    filtered = {key: value for key, value in params.items() if value}
    if not filtered:
        return path
    return f"{path}?{urlencode(filtered)}"


def _page_tabs(
    active: str,
    *,
    controller_log: str = "",
    user: str = "",
    remote_host: str = "",
):
    logs_href = _query_url(
        "/",
        controller_log=controller_log,
        user=user,
        remote_host=remote_host,
    )
    submit_href = _query_url(
        "/submit",
        user=user,
        remote_host=remote_host,
    )
    return Div(  # noqa: F405
        A(  # noqa: F405
            "Logs",
            href=logs_href,
            cls=f"fab-page-tab{' is-active' if active == 'logs' else ''}",
        ),
        A(  # noqa: F405
            "Submit",
            href=submit_href,
            cls=f"fab-page-tab{' is-active' if active == 'submit' else ''}",
        ),
        cls="fab-page-tabs",
    )


def _page_header(
    active: str,
    *,
    controller_log: str = "",
    user: str = "",
    remote_host: str = "",
):
    shutdown_button = Button(  # noqa: F405
        "Shutdown",
        submit=False,
        cls=(ButtonT.secondary, ButtonT.sm),  # noqa: F405
        hx_post="/shutdown-server",
        hx_target="#server-action-status",
        hx_swap="innerHTML",
        hx_confirm="Shut down the dashboard server?",
    )
    return Div(  # noqa: F405
        DivFullySpaced(  # noqa: F405
            _page_tabs(
                active,
                controller_log=controller_log,
                user=user,
                remote_host=remote_host,
            ),
            shutdown_button,
        ),
        Div(id="server-action-status"),
        cls="space-y-2",
    )


def _submission_controls(user: str, remote_host: str):
    return Card(  # noqa: F405
        Form(  # noqa: F405
            Div(
                Div(
                    P("User", cls="text-sm font-medium"),
                    Input(  # noqa: F405
                        id="submit-user",
                        name="user",
                        value=user,
                        placeholder="bufi0980",
                        autocomplete="off",
                        cls="uk-input uk-form-small w-full",
                    ),
                    cls="space-y-2 min-w-0",
                    style="flex: 0.7 1 10rem;",
                ),
                Div(
                    P("Host", cls="text-sm font-medium"),
                    Input(  # noqa: F405
                        id="submit-remote-host",
                        name="remote_host",
                        value=remote_host,
                        placeholder="rosa.hpc.uni-oldenburg.de",
                        autocomplete="off",
                        cls="uk-input uk-form-small w-full",
                    ),
                    cls="space-y-2 min-w-0",
                    style="flex: 1.1 1 18rem;",
                ),
                cls="flex flex-wrap gap-2 items-end",
            ),
            id="submission-controls-form",
            cls="space-y-0",
        ),
        header=Div(  # noqa: F405
            H3("Connection", cls="text-base font-semibold"),  # noqa: F405
            P(
                "Submission uses the same remote account and host as `fab submit`.",
                cls="text-sm opacity-70",
            ),
            cls="space-y-1",
        ),
        body_cls="space-y-0 py-3",
    )


def _safe_dom_id(value: str) -> str:
    safe_value = re.sub(r"[^A-Za-z0-9_-]+", "-", value).strip("-")
    return safe_value or "item"


def _submission_rule_option(rule: dict[str, str]):
    checkbox_id = f"submission-rule-{_safe_dom_id(rule['name'])}"
    label_id = f"{checkbox_id}-label"
    description_id = f"{checkbox_id}-description"
    return Div(  # noqa: F405
        Input(  # noqa: F405
            id=checkbox_id,
            type="checkbox",
            value=rule["name"],
            checked=True,
            data_submit_rule="true",
            aria_labelledby=label_id,
            aria_describedby=description_id,
        ),
        Div(
            P(rule["label"], id=label_id, cls="text-sm font-medium"),
            P(rule["description"], id=description_id, cls="text-xs opacity-70"),
            cls="space-y-1 min-w-0",
        ),
        cls="fab-submit-rule-option",
        data_checkbox_card="true",
        onclick=f"window.fabMonitorDashboard.toggleCheckboxCard('{checkbox_id}', event);",
    )


def _submission_group_actions(selector: str, *, include_new: bool = False):
    buttons = [
        Button(  # noqa: F405
            "All",
            submit=False,
            cls=(ButtonT.ghost, ButtonT.sm),  # noqa: F405
            onclick=f"window.fabMonitorDashboard.setSubmissionGroupChecked('{selector}', true);",
        ),
        Button(  # noqa: F405
            "None",
            submit=False,
            cls=(ButtonT.ghost, ButtonT.sm),  # noqa: F405
            onclick=f"window.fabMonitorDashboard.setSubmissionGroupChecked('{selector}', false);",
        ),
    ]
    if include_new:
        buttons.insert(
            1,
            Button(  # noqa: F405
                "New",
                submit=False,
                cls=(ButtonT.ghost, ButtonT.sm),  # noqa: F405
                onclick="window.fabMonitorDashboard.selectNewExperiments();",
            ),
        )
    return Div(*buttons, cls="flex items-center gap-2")  # noqa: F405


def _experiment_status_summary(experiment: dict[str, object]):
    rule_status = experiment.get("rule_status", {})
    done = sum(1 for status in rule_status.values() if status)
    total = len(rule_status)
    state = "finished" if total and done == total else "new"
    detail = f"{done}/{total} rules done" if total else "No curated rules"
    return state, detail


def _submission_experiment_option(experiment: dict[str, object]):
    state, detail = _experiment_status_summary(experiment)
    badge_cls = f"fab-submit-status-badge {'is-finished' if state == 'finished' else 'is-new'}"
    checkbox_id = f"submission-experiment-{_safe_dom_id(str(experiment['id']))}"
    label_id = f"{checkbox_id}-label"
    detail_id = f"{checkbox_id}-detail"
    return Div(  # noqa: F405
        Input(  # noqa: F405
            id=checkbox_id,
            type="checkbox",
            value=str(experiment["id"]),
            data_submit_experiment="true",
            aria_labelledby=label_id,
            aria_describedby=detail_id,
        ),
        Div(
            Div(
                P(str(experiment["label"]), id=label_id, cls="text-sm font-medium break-all"),
                P(
                    f"{experiment['group']} | {experiment['path']}",
                    cls="text-xs opacity-70 break-all",
                ),
                P(
                    detail,
                    id=detail_id,
                    cls="text-xs opacity-70",
                    data_experiment_status_detail="true",
                ),
                cls="space-y-1 min-w-0",
            ),
            Span(state, cls=badge_cls, data_experiment_status_badge="true"),  # noqa: F405
            cls="flex items-start justify-between gap-3 w-full",
        ),
        cls="fab-submit-experiment-option",
        data_checkbox_card="true",
        data_experiment_id=str(experiment["id"]),
        onclick=f"window.fabMonitorDashboard.toggleCheckboxCard('{checkbox_id}', event);",
    )


def _submission_panel():
    metadata = submission_metadata(PROJECT_DIR)
    rules = metadata.get("rules", [])
    experiments = metadata.get("experiments", [])

    quick_submit = Card(  # noqa: F405
        P(str(metadata.get("quick_submit_help", "")), cls="text-sm opacity-75"),
        Form(  # noqa: F405
            Input(type="hidden", name="submit_mode", value="all"),  # noqa: F405
            Input(type="hidden", name="targets_text", value=""),  # noqa: F405
            Input(type="hidden", name="selection_summary", value="Everything pending"),  # noqa: F405
            Button(  # noqa: F405
                "Submit Everything Pending",
                submit=False,
                cls=(ButtonT.primary, ButtonT.sm),  # noqa: F405
                onclick="window.fabMonitorDashboard.showSubmissionConfirm('all');",
            ),
            hx_post="/submit-launch",
            hx_target="#submission-status",
            hx_swap="innerHTML",
            hx_include="#submission-controls-form",
        ),
        header=H3("Submit All", cls="text-base font-semibold"),  # noqa: F405
        body_cls="space-y-3",
    )

    selective_submit = Card(  # noqa: F405
        Form(  # noqa: F405
            Input(type="hidden", name="submit_mode", value="selection"),  # noqa: F405
            Input(type="hidden", id="selection_summary", name="selection_summary", value=""),  # noqa: F405
            Textarea(json.dumps(metadata), id="submission-metadata", hidden=True, style="display:none;"),  # noqa: F405
            Textarea("", id="targets_text", name="targets_text", hidden=True, style="display:none;"),  # noqa: F405
            P(str(metadata.get("selective_submit_help", "")), cls="text-sm opacity-75"),
            Div(
                Div(
                    DivFullySpaced(  # noqa: F405
                        P("Rules", cls="text-sm font-medium"),
                        _submission_group_actions('[data-submit-rule="true"]'),
                    ),
                    Div(
                        *[_submission_rule_option(rule) for rule in rules],
                        cls="fab-submit-checklist",
                    )
                    if rules
                    else P("No curated rules are available for this workflow.", cls="text-sm opacity-70"),  # noqa: F405
                    cls="fab-submit-column",
                ),
                Div(
                    Div(
                        DivFullySpaced(  # noqa: F405
                            P("Experiments", cls="text-sm font-medium"),
                            _submission_group_actions(
                                '[data-submit-experiment="true"]',
                                include_new=True,
                            ),
                        ),
                        P(
                            "Finished experiments stay selectable for reruns.",
                            cls="text-xs opacity-70",
                        ),
                        cls="space-y-1",
                    ),
                    Div(
                        Div(
                            *[_submission_experiment_option(experiment) for experiment in experiments],
                            cls="fab-submit-checklist fab-submit-scroll",
                        )
                        if experiments
                        else P("No experiments were discovered under `dat/`.", cls="text-sm opacity-70"),  # noqa: F405
                        cls="fab-submit-checklist-shell",
                    ),
                    cls="fab-submit-column",
                ),
                cls="fab-submit-section-grid",
            ),
            Div(
                P("Preview", cls="text-sm font-medium"),
                Div(
                    P(
                        "Select at least one rule and one experiment.",
                        id="submission-preview-summary",
                        cls="text-sm opacity-75",
                    ),
                    Pre(id="submission-preview-targets", cls="fab-submit-preview-pre"),  # noqa: F405
                    cls="fab-submit-preview space-y-2",
                ),
                cls="space-y-2",
            ),
            Div(
                Button(  # noqa: F405
                    "Submit Selection",
                    id="submit-selection-button",
                    submit=False,
                    cls=(ButtonT.primary, ButtonT.sm),  # noqa: F405
                    disabled=True,
                    onclick="window.fabMonitorDashboard.showSubmissionConfirm('selection');",
                ),
                cls="flex justify-end",
            ),
            id="submission-form",
            cls="space-y-4",
            hx_post="/submit-launch",
            hx_target="#submission-status",
            hx_swap="innerHTML",
            hx_include="#submission-controls-form, #submission-form",
        ),
        header=H3("Selective Submit", cls="text-base font-semibold"),  # noqa: F405
        body_cls="space-y-0",
    )

    return Div(  # noqa: F405
        quick_submit,
        selective_submit,
        Div(id="submission-status"),
        id="submission-page-shell",
        cls="space-y-3",
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

    dag_button = Button(  # noqa: F405
        "Workflow DAG",
        submit=False,
        cls=(ButtonT.secondary, ButtonT.sm),  # noqa: F405
        onclick="window.fabMonitorDashboard.showLoading('dag-modal-loading-template');",
        hx_get="/dag-modal",
        hx_target="#modal-container",
        hx_swap="innerHTML",
        hx_include="#controls-form",
        disabled=not run.jobs,
    )
    cancel_disabled = not run.slurm_run_id or run.is_terminal()
    cancel_title = (
        "The selected run is already complete."
        if run.is_terminal()
        else (
            f"Cancel SLURM run {run.slurm_run_id}"
            if run.slurm_run_id
            else "No SLURM run ID was found for the selected controller log."
        )
    )
    cancel_button = Button(  # noqa: F405
        "Cancel Run",
        submit=False,
        cls=(ButtonT.secondary, ButtonT.sm, "border-red-300 text-red-700 hover:bg-red-50"),  # noqa: F405
        title=cancel_title,
        hx_post="/cancel-run",
        hx_target="#run-action-status",
        hx_swap="innerHTML",
        hx_include="#controls-form",
        hx_confirm=(f"Cancel SLURM run {run.slurm_run_id}?" if run.slurm_run_id else None),
        disabled=cancel_disabled,
    )

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
            Div(
                Span("running ", Strong(str(counts["running"]))),  # noqa: F405
                Span("pending ", Strong(str(counts["submitted"] - counts["running"]))),  # noqa: F405
                Span("failed ", Strong(str(counts["failed"]))),  # noqa: F405
                cls="flex flex-wrap gap-4 text-sm whitespace-nowrap justify-end",
            ),
            Div(
                dag_button,
                cancel_button,
                cls="flex items-center justify-end gap-2 flex-wrap",
            ),
            cls="flex flex-col items-end gap-2",
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
    blocks.append(Div(id="run-action-status"))
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


def _dag_modal_loading_template():
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
                    P("Loading workflow DAG", cls="text-base font-semibold"),
                    P(
                        "Deriving the rule graph from the selected controller log.",
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
        onclick="window.fabMonitorDashboard.showLoading('job-modal-loading-template');",
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


def _dag_view(run):
    graph, used_workflow_graph = _dag_graph_for_run(run)
    node_count = len(graph.get("nodes", []))
    edge_count = len(graph.get("edges", []))
    if not node_count:
        return Card(  # noqa: F405
            P("No rule graph could be derived from the selected controller log.", cls="text-sm"),
            cls="border border-amber-300 bg-amber-50",
        )

    svg = _rule_graph_svg(graph)
    return Div(  # noqa: F405
        Div(
            P("Workflow DAG", cls="text-sm font-medium"),
            P(
                (
                    f"{node_count} rules | {edge_count} edges. "
                    + (
                        "Grey nodes have not been materialized in this run yet."
                        if used_workflow_graph
                        else "The graph is derived from jobs materialized in this run."
                    )
                ),
                cls="text-sm opacity-75",
            ),
            cls="space-y-1",
        ),
        _dag_legend(),
        Div(NotStr(svg), cls="fab-monitor-dag-shell"),
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


def _submission_status_card(
    *,
    ok: bool,
    title: str,
    detail: str = "",
    selection_summary: str = "",
    controller_log: str = "",
    stderr: str = "",
    user: str = "",
    remote_host: str = "",
):
    tone_cls = (
        "border border-emerald-300 bg-emerald-50"
        if ok
        else "border border-red-300 bg-red-50"
    )
    blocks = [P(title, cls="text-sm font-medium")]  # noqa: F405
    if selection_summary:
        blocks.append(P(selection_summary, cls="text-xs opacity-70"))
    if detail:
        blocks.append(P(detail, cls="text-xs opacity-70 break-all"))
    if controller_log:
        blocks.append(P(controller_log, cls="text-xs font-mono break-all"))
    if stderr:
        blocks.append(Pre(stderr, cls="fab-submit-preview-pre"))  # noqa: F405

    if ok:
        blocks.append(
            A(  # noqa: F405
                "Open Logs",
                href=_query_url(
                    "/",
                    controller_log=controller_log,
                    user=user,
                    remote_host=remote_host,
                ),
                cls="fab-page-tab",
            )
        )
    return Div(Card(*blocks, cls=tone_cls, body_cls="space-y-2"))  # noqa: F405


def _server_action_status(message: str, *, tone: str = "neutral"):
    tone_cls = {
        "success": "border border-emerald-300 bg-emerald-50",
        "danger": "border border-red-300 bg-red-50",
    }.get(tone, "border border-slate-300 bg-slate-50")
    return Card(  # noqa: F405
        P(message, cls="text-sm"),  # noqa: F405
        cls=tone_cls,
        body_cls="space-y-0 py-2",
    )


def _run_action_response(message: str, *, tone: str = "neutral", refresh_dashboard: bool = False):
    blocks = [_server_action_status(message, tone=tone)]
    if refresh_dashboard:
        blocks.append(
            Script(  # noqa: F405
                """
(() => {
  setTimeout(() => {
    if (
      window.fabMonitorDashboard &&
      typeof window.fabMonitorDashboard.refreshDashboardShell === 'function'
    ) {
      window.fabMonitorDashboard.refreshDashboardShell();
    }
  }, 700);
})();
                """.strip()
            )
        )
    return Div(*blocks)  # noqa: F405


def _server_shutdown_response():
    return Div(  # noqa: F405
        _server_action_status(
            "Dashboard server is shutting down. This tab should close automatically.",
            tone="danger",
        ),
        Script(  # noqa: F405
            """
(() => {
  const tryClose = () => {
    try { window.open('', '_self'); } catch (error) {}
    try { window.close(); } catch (error) {}
  };
  setTimeout(tryClose, 40);
  setTimeout(tryClose, 220);
  setTimeout(() => {
    tryClose();
    try { window.location.replace('about:blank'); } catch (error) {}
    tryClose();
  }, 500);
})();
            """.strip()
        ),
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
        _page_header("logs", controller_log=controller_log, user=user, remote_host=remote_host),
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
        Div(_dag_modal_loading_template(), id="dag-modal-loading-template", hidden=True),  # noqa: F405
        cls="space-y-3 py-3",
    )


@rt("/submit")  # noqa: F405
def submit_page(
    user: str = DEFAULT_USER,
    remote_host: str = DEFAULT_REMOTE_HOST,
):
    return Title(f"Snakemake Submit | {PROJECT_DIR.name}"), Container(  # noqa: F405
        Style(LOG_BROWSER_CSS),  # noqa: F405
        Script(LOG_BROWSER_SCRIPT),  # noqa: F405
        _page_header("submit", user=user, remote_host=remote_host),
        _submission_controls(user, remote_host),
        _submission_panel(),
        Div(id="modal-container"),  # noqa: F405
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


@app.post("/shutdown-server")  # noqa: F405
def shutdown_server():
    _schedule_server_shutdown()
    return _server_shutdown_response()


@app.post("/cancel-run")  # noqa: F405
def cancel_run(
    controller_log: str = "",
    user: str = DEFAULT_USER,
    remote_host: str = DEFAULT_REMOTE_HOST,
):
    run, _, _ = _load_run(controller_log, user, remote_host)
    if run is None:
        return _run_action_response("No controller run is selected.", tone="danger")

    if not run.slurm_run_id:
        return _run_action_response(
            "The selected controller log does not include a SLURM run ID.",
            tone="danger",
        )

    if run.is_terminal():
        return _run_action_response(
            "The selected run is already complete.",
            tone="neutral",
        )

    queue_user = (user or "").strip()
    queue_host = (remote_host or "").strip()
    if not queue_user:
        return _run_action_response(
            "Set an SSH user before canceling the run.",
            tone="danger",
        )
    if not queue_host:
        return _run_action_response(
            "Set an SSH host before canceling the run.",
            tone="danger",
        )

    ok, detail = cancel_slurm_run(
        user=queue_user,
        host=queue_host,
        run_id=run.slurm_run_id,
    )
    message = detail or (
        f"Cancellation requested for SLURM run {run.slurm_run_id}."
        if ok
        else f"Failed to cancel SLURM run {run.slurm_run_id}."
    )
    return _run_action_response(
        message,
        tone="success" if ok else "danger",
        refresh_dashboard=ok,
    )


@app.post("/submit-launch")  # noqa: F405
def submit_launch(
    user: str = DEFAULT_USER,
    remote_host: str = DEFAULT_REMOTE_HOST,
    submit_mode: str = "selection",
    targets_text: str = "",
    selection_summary: str = "",
):
    targets = []
    if submit_mode != "all":
        targets = normalize_submission_targets(targets_text)
        if not targets:
            return _submission_status_card(
                ok=False,
                title="No submission targets were selected.",
                detail="Choose at least one rule and one experiment.",
                user=user,
                remote_host=remote_host,
            )

    result = launch_dashboard_submission(
        project_dir=PROJECT_DIR,
        user=user,
        host=remote_host,
        targets=targets,
    )
    detail = result.detail
    return _submission_status_card(
        ok=result.ok,
        title=result.message,
        detail=detail,
        selection_summary=selection_summary,
        controller_log=result.controller_log,
        stderr="" if result.ok or result.stderr == detail else result.stderr,
        user=user,
        remote_host=remote_host,
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


@rt("/dag-modal")  # noqa: F405
def dag_modal(controller_log: str = ""):
    selected_log = _selected_log_path(controller_log)
    if not selected_log:
        return Div()  # noqa: F405

    run = parse_controller_log(Path(selected_log))
    graph, _ = _dag_graph_for_run(run)
    header = Div(
        H3("Workflow DAG", cls="text-lg font-semibold"),  # noqa: F405
        P(Path(selected_log).name, cls="text-sm opacity-70"),
        cls="space-y-1",
    )

    footer = Div(
        P(
            (
                f"{len(graph.get('nodes', []))} rules | "
                f"{len(graph.get('edges', []))} edges"
            ),
            cls="text-xs opacity-70",
        ),
        ModalCloseButton("Close", cls=(ButtonT.secondary, "static"), submit=False),  # noqa: F405
        cls="flex items-center justify-between gap-3",
    )

    return Modal(  # noqa: F405
        _dag_view(run),
        header=header,
        footer=footer,
        id="workflow-dag-modal",
        hx_init=True,
        hx_open=True,
        dialog_cls="uk-width-11-12 max-h-[94vh]",
        body_cls="space-y-4 overflow-y-auto",
        style="padding-top: 2vh; padding-bottom: 2vh;",
    )


_schedule_browser_open()

serve(  # noqa: F405
    host=DEFAULT_HOST,
    port=DEFAULT_PORT,
    reload=False,
    access_log=DEFAULT_VERBOSE,
)
