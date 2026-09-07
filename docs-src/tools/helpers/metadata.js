// Samplesheet builder for the metadata generator.
//
// Every dropdown is built from OSCAR_SCHEMA (helpers/schema.js), which is
// generated from lib/samplesheet.nf and lib/chemistry.nf. Hard-coding these
// lists is what let the tool drift and offer chemistries the pipeline rejects.

const FIELD_HINTS = {
  experiment_id: { placeholder: 'e.g. ILC_niche', title: 'Short name of the experiment' },
  historical_number: { placeholder: 'e.g. 1', title: 'Pool number: one per pool of cells loaded onto the chip' },
  replicate: { placeholder: 'e.g. A', title: 'Well within a pool: A, B, C...' },
  index: { placeholder: 'e.g. SI-TT-A1 or CATAGCCG', title: '10x kit code, or a raw index sequence' },
  n_donors: { placeholder: 'e.g. 2 or NA', title: 'Donors in this well. Above 1 runs vireo. Human only.' },
  adt_file: { placeholder: 'e.g. my_panel or NA', title: 'ADT CSV basename, without .csv' },
};

function optionsHTML(items, { placeholder = 'unselected' } = {}) {
  const opts = items.map((item) => {
    const value = typeof item === 'string' ? item : item.value;
    const label = typeof item === 'string' ? item : item.label;
    return `<option value="${value}">${label}</option>`;
  });
  return `<option value="" selected>${placeholder}</option>${opts.join('')}`;
}

// One <datalist> shared by every row, rather than 480 options duplicated per
// row. Lets people type a raw index sequence, which a <select> cannot.
function ensureIndexDatalist() {
  if (document.getElementById('index-codes')) return;
  const dl = document.createElement('datalist');
  dl.id = 'index-codes';
  dl.innerHTML = OSCAR_SCHEMA.indexCodes
    .map((code) => `<option value="${code}"></option>`)
    .join('');
  document.body.appendChild(dl);
}

function cell(inner) {
  return `<td>${inner}</td>`;
}

function buildRow(rowIndex) {
  const id = (field) => `${field}-${rowIndex}`;
  const hint = (field) => {
    const h = FIELD_HINTS[field] || {};
    return `placeholder="${h.placeholder || ''}" title="${h.title || ''}"`;
  };
  // Each control carries its column name, so the CSV writer reads by name
  // rather than by position. Reordering a column cannot silently corrupt output.
  const named = (field) => `data-column="${field}" name="${field}" id="${id(field)}"`;

  return [
    cell(`<select ${named('assay')} class="assay" required aria-label="Assay type">
            ${optionsHTML(OSCAR_SCHEMA.assays)}
          </select>`),
    cell(`<input type="text" ${named('experiment_id')} class="experiment_id" required
                 aria-label="Experiment ID" ${hint('experiment_id')}>`),
    cell(`<input type="number" min="1" step="1" ${named('historical_number')}
                 class="historical_number" required aria-label="Historical number"
                 ${hint('historical_number')}>`),
    cell(`<input type="text" ${named('replicate')} class="replicate" required
                 pattern="[A-Za-z0-9]+" aria-label="Replicate" ${hint('replicate')}>`),
    cell(`<select ${named('modality')} class="modality" required aria-label="Modality">
            ${optionsHTML(OSCAR_SCHEMA.modalities)}
          </select>`),
    cell(`<select ${named('chemistry')} class="chemistry" required aria-label="Chemistry">
            ${optionsHTML(OSCAR_SCHEMA.chemistries)}
          </select>`),
    cell(`<select ${named('index_type')} class="index_type" required aria-label="Index type">
            ${optionsHTML(OSCAR_SCHEMA.indexTypes)}
          </select>`),
    cell(`<input type="text" list="index-codes" ${named('index')} class="index" required
                 aria-label="Index" ${hint('index')}>`),
    cell(`<select ${named('species')} class="species" required aria-label="Species">
            ${optionsHTML(['human', 'mouse'])}
          </select>`),
    cell(`<input type="text" ${named('n_donors')} class="n_donors" required
                 aria-label="Number of donors" ${hint('n_donors')}>`),
    cell(`<input type="text" ${named('adt_file')} class="adt_file" required
                 aria-label="ADT file name" ${hint('adt_file')}>`),
    cell(`<div class="row-actions">
            <button type="button" class="row-duplicate" aria-label="Duplicate this row"
                    title="Duplicate row">&plus;</button>
            <button type="button" class="row-delete" aria-label="Delete this row"
                    title="Delete row">&times;</button>
          </div>`),
  ].join('');
}

function addRow(afterRow = null) {
  const container = document.getElementById('rowsContainer');
  const tr = document.createElement('tr');
  tr.innerHTML = buildRow(container.rows.length);
  if (afterRow) {
    afterRow.after(tr);
  } else {
    container.appendChild(tr);
  }
  return tr;
}

// Copy a row's values into a new row below it. Most samplesheets repeat a
// library across modalities, changing one or two fields.
function duplicateRow(sourceRow) {
  const clone = addRow(sourceRow);
  const from = sourceRow.querySelectorAll('[data-column]');
  const to = clone.querySelectorAll('[data-column]');
  from.forEach((src, i) => { to[i].value = src.value; });
  return clone;
}

function readRow(row) {
  const values = {};
  let valid = true;
  row.querySelectorAll('[data-column]').forEach((el) => {
    const value = el.value.trim();
    if (value === '') {
      el.classList.add('error');
      el.setAttribute('aria-invalid', 'true');
      valid = false;
    } else {
      el.classList.remove('error');
      el.removeAttribute('aria-invalid');
    }
    values[el.dataset.column] = value;
  });
  return { valid, values };
}

function escapeCSV(value) {
  return /[",\n]/.test(value) ? `"${value.replace(/"/g, '""')}"` : value;
}

function showMessage(text, kind = 'error') {
  const box = document.getElementById('messageContainer');
  box.textContent = text;
  box.className = `message message-${kind}`;
  box.style.display = text ? 'block' : 'none';
}

function generateCSV() {
  const rows = Array.from(document.querySelectorAll('#rowsContainer tr'));
  if (rows.length === 0) {
    showMessage('Add at least one row before downloading.');
    return;
  }

  const columns = OSCAR_SCHEMA.columns;
  const lines = [columns.join(',')];
  let firstBad = null;

  rows.forEach((row) => {
    const { valid, values } = readRow(row);
    if (!valid) {
      if (!firstBad) firstBad = row;
      return;
    }
    lines.push(columns.map((c) => escapeCSV(values[c] ?? 'NA')).join(','));
  });

  if (firstBad) {
    showMessage('Fill in every highlighted field before downloading.');
    firstBad.querySelector('.error')?.focus();
    return;
  }

  showMessage('');
  // BOM keeps Excel from mangling UTF-8.
  const blob = new Blob(['﻿' + lines.join('\n') + '\n'], {
    type: 'text/csv;charset=utf-8;',
  });
  const url = URL.createObjectURL(blob);
  const a = document.createElement('a');
  a.href = url;
  a.download = 'metadata.csv';
  document.body.appendChild(a);
  a.click();
  a.remove();
  URL.revokeObjectURL(url);
}

document.addEventListener('DOMContentLoaded', () => {
  ensureIndexDatalist();
  addRow();

  document.querySelector('.add-row-button').addEventListener('click', () => addRow());
  document.querySelectorAll('.generate-csv-button').forEach((button) => {
    button.addEventListener('click', generateCSV);
  });

  // Delegated so it covers rows added later.
  document.getElementById('rowsContainer').addEventListener('click', (event) => {
    const row = event.target.closest('tr');
    if (!row) return;
    if (event.target.matches('.row-delete')) {
      row.remove();
      showMessage('');
    } else if (event.target.matches('.row-duplicate')) {
      duplicateRow(row);
    }
  });
});
