// functions.js

// Utility function to escape CSV values
function escapeCSV(value) {
    const stringValue = String(value || '');
    if (stringValue.includes(',') || stringValue.includes('"') || stringValue.includes('\n')) {
        return `"${stringValue.replace(/"/g, '""')}"`;
    }
    return stringValue;
}

// Navigation bar.
//
// Built in place rather than fetched. The markup is six static links, and a
// failed fetch left the page with no navigation at all.
const NAV_LINKS = [
    { href: '../',                      text: 'Home' },
    { href: 'metadata_generator.html',  text: 'Metadata file generator' },
    { href: 'adt_generator.html',       text: 'Feature barcode generator' },
    { href: '../guide/samplesheet/',    text: 'Samplesheet guide' },
    { href: '../reference/parameters/', text: 'Parameters' },
];

function includeTopBar() {
    const target = document.getElementById('topbar-placeholder');
    if (!target) return;
    const nav = document.createElement('nav');
    nav.className = 'top-bar';
    nav.setAttribute('aria-label', 'Tools navigation');
    const here = window.location.pathname.split('/').pop();
    NAV_LINKS.forEach(({ href, text }) => {
        const a = document.createElement('a');
        a.href = href;
        a.textContent = text;
        if (href === here) a.setAttribute('aria-current', 'page');
        nav.appendChild(a);
    });
    target.replaceChildren(nav);
}

// Retained for the old entry points.
const indexIncludeTopBar = includeTopBar;
const pagesIncludeTopBar = includeTopBar;

function fetchLastCommitDate(owner, repo) {
    const url = `https://api.github.com/repos/${owner}/${repo}/commits`;
    return fetch(url)
        .then(response => response.json())
        .then(data => {
            if (data && data.length > 0) {
                return data[0].commit.author.date;
            } else {
                throw new Error('No commits found');
            }
        });
}

// Metadata generator logic lives in helpers/metadata.js.

// Functions for adt_generator.html

// Markers loaded from the selected TotalSeq CSV. Declared explicitly: it was
// an implicit global created by assignment, which breaks under strict mode.
let markersData = [];
const selectedPanelValues = { species: '', format: '', 'output-format': '' };

async function adtFetchCSV(file) {
    try {
        const response = await fetch(file);
        if (!response.ok) {
            throw new Error(`HTTP error! status: ${response.status}`);
        }
        const data = await response.text();
        return data;
    } catch (error) {
        console.error('Error fetching CSV:', error);
        adtShowMessage('Could not load marker data. Refresh and try again.');
        return '';
    }
}

function adtShowMessage(text, kind = 'error') {
    const box = document.getElementById('messageContainer');
    if (!box) return;
    box.textContent = text;
    box.className = `message message-${kind}`;
    box.hidden = !text;
}

function adtParseCSV(data) {
    const rows = data.split('\n').slice(1);
    return rows.map(row => {
        const [catalogue_number, totalseq_id, marker, clone, reactivity, barcode_sequence] = row.split(',');
        return { catalogue_number, totalseq_id, marker, clone, reactivity, barcode_sequence };
    });
}

async function adtFilterMarkers() {
    const format = document.getElementById('format').value;
    const species = document.getElementById('species').value;
    const searchInput = document.getElementById('search');

    if (!format || !species) return;
    
    // Show loading state
    searchInput.disabled = true;
    searchInput.placeholder = 'Loading markers...';
    
    try {
        const data = await adtFetchCSV(format);
        if (!data) {
            throw new Error('No data returned');
        }
        
        markersData = adtParseCSV(data);

        if (species === 'Human') {
            markersData = markersData.filter(marker => ['Human', 'Human and mouse', 'Isotype control', 'Epitope'].includes(marker.reactivity));
        } else if (species === 'Mouse') {
            markersData = markersData.filter(marker => ['Mouse', 'Human and mouse', 'Isotype control', 'Epitope'].includes(marker.reactivity));
        }

        // Call adtShowDropdown after data is filtered
        adtShowDropdown(format);
    } catch (error) {
        console.error('Error filtering markers:', error);
        markersData = [];
    } finally {
        // Restore input state
        searchInput.disabled = false;
        searchInput.placeholder = 'Marker, clone, or TotalSeq ID';
    }
}

function adtRemoveRow(element) {
    const row = element.parentElement.parentElement;
    row.remove();
}

function adtShowDropdown(format) {
    // Derive the TotalSeq letter once. Falls back to the select's own value so
    // a bare call still works, and to '' rather than throwing on no match.
    const source = format || document.getElementById('format').value || '';
    const formatLetter = (source.match(/totalseq_(\w)/i) || [, ''])[1].toUpperCase();

    const search = document.getElementById('search').value.toLowerCase();
    const dropdown = document.getElementById('dropdown-content');
    dropdown.innerHTML = '';

    if (search && markersData.length > 0) {
        const filteredMarkers = markersData.filter(marker => 
            marker.catalogue_number.toLowerCase().includes(search) ||
            marker.marker.toLowerCase().includes(search) ||
            marker.clone.toLowerCase().includes(search) ||
            marker.totalseq_id.toLowerCase().includes(search)
        );
        filteredMarkers.forEach(marker => {
            const button = document.createElement('button');
            button.type = 'button';
            button.className = 'dropdown-item';
            const text = `TotalSeq-${formatLetter}${marker.totalseq_id}, ${marker.marker}, ${marker.clone}`;
            button.textContent = text;
            button.addEventListener('click', () => {
                adtAddRow(marker.marker, marker.totalseq_id, marker.catalogue_number, marker.clone, marker.reactivity, marker.barcode_sequence);
                document.getElementById('search').value = `${formatLetter}${marker.totalseq_id}`;
                dropdown.innerHTML = '';
                dropdown.style.display = 'none';
            });
            dropdown.appendChild(button);
        });
        dropdown.style.display = filteredMarkers.length ? 'block' : 'none';
    } else {
        dropdown.style.display = 'none';
    }
}

function adtConfirmReset(type) {
    const speciesElement = document.getElementById('species');
    const formatElement = document.getElementById('format');
    const outputFormatElement = document.getElementById('output-format');
    const elements = { species: speciesElement, format: formatElement, 'output-format': outputFormatElement };
    const element = elements[type];
    const nextValue = element.value;
    const hasEntries = document.getElementById('rowsContainer').children.length > 0;

    if (hasEntries && nextValue !== selectedPanelValues[type] && !confirm('Changing this clears selected markers. Continue?')) {
        element.value = selectedPanelValues[type];
        return;
    }
    if (hasEntries && nextValue !== selectedPanelValues[type]) {
        document.getElementById('rowsContainer').replaceChildren();
        markersData = [];
        document.getElementById('search').value = '';
    }
    selectedPanelValues[type] = nextValue;
    if (type === 'species' || type === 'format') adtFilterMarkers();
}

function adtSubstituteCharacters(text) {
    // Comprehensive Greek letter to English letter mapping
    const greekMap = {
        'α': 'a', 'Α': 'A',
        'β': 'b', 'Β': 'B',
        'γ': 'g', 'Γ': 'G',
        'δ': 'd', 'Δ': 'D',
        'ε': 'e', 'Ε': 'E',
        'ζ': 'z', 'Ζ': 'Z',
        'η': 'h', 'Η': 'H',
        'θ': 'th', 'Θ': 'Th',
        'ι': 'i', 'Ι': 'I',
        'κ': 'k', 'Κ': 'K',
        'λ': 'l', 'Λ': 'L',
        'μ': 'm', 'Μ': 'M',
        'ν': 'n', 'Ν': 'N',
        'ξ': 'x', 'Ξ': 'X',
        'ο': 'o', 'Ο': 'O',
        'π': 'p', 'Π': 'P',
        'ρ': 'r', 'Ρ': 'R',
        'σ': 's', 'ς': 's', 'Σ': 'S',
        'τ': 't', 'Τ': 'T',
        'υ': 'u', 'Υ': 'U',
        'φ': 'ph', 'Φ': 'Ph',
        'χ': 'ch', 'Χ': 'Ch',
        'ψ': 'ps', 'Ψ': 'Ps',
        'ω': 'o', 'Ω': 'O'
    };
    
    return text
        .replace(/alpha/gi, 'a')
        .replace(/beta/gi, 'b')
        .replace(/gamma/gi, 'g')
        .replace(/delta/gi, 'd')
        .replace(/epsilon/gi, 'e')
        .replace(/zeta/gi, 'z')
        .replace(/eta/gi, 'h')
        .replace(/theta/gi, 'th')
        .replace(/iota/gi, 'i')
        .replace(/kappa/gi, 'k')
        .replace(/lambda/gi, 'l')
        .replace(/mu/gi, 'm')
        .replace(/nu/gi, 'n')
        .replace(/xi/gi, 'x')
        .replace(/omicron/gi, 'o')
        .replace(/pi/gi, 'p')
        .replace(/rho/gi, 'r')
        .replace(/sigma/gi, 's')
        .replace(/tau/gi, 't')
        .replace(/upsilon/gi, 'u')
        .replace(/phi/gi, 'ph')
        .replace(/chi/gi, 'ch')
        .replace(/psi/gi, 'ps')
        .replace(/omega/gi, 'o')
        .replace(/[αβγδεζηθικλμνξοπρςστυφχψωΑΒΓΔΕΖΗΘΙΚΛΜΝΞΟΠΡΣΤΥΦΧΨΩ]/g, match => greekMap[match] || match)
        .replace(/\s+/g, '')
        .replace(/\//g, '-')
        .replace(/\./g, '_')
        .replace(/,/g, '_');
}

function adtEditRow(event, button) {
    event.preventDefault();
    const row = button.closest('td');
    const span = row.querySelector('.corrected-name');
    const input = row.querySelector('.edit-input');
    const saveButton = row.querySelector('.save-button');

    input.value = span.textContent; // Set input value to current span text
    span.style.display = 'none';
    input.style.display = 'inline';
    button.style.display = 'none';
    saveButton.style.display = 'inline';
}

function adtSaveRow(event, button) {
    event.preventDefault();
    const row = button.closest('td');
    const span = row.querySelector('.corrected-name');
    const input = row.querySelector('.edit-input');
    const editButton = row.querySelector('.edit-button');

    span.textContent = input.value;
    span.style.display = 'inline';
    input.style.display = 'none';
    button.style.display = 'none';
    editButton.style.display = 'inline';
}

function adtAddRow(marker, totalseq_id, catalogue_number, clone, reactivity, barcode) {
    const tableBody = document.getElementById('rowsContainer');
    const row = document.createElement('tr');

    row.innerHTML = `
        <td>${totalseq_id}</td>
        <td>${marker}</td>
        <td>
            <div class="edit-container">
                <span class="corrected-name">${adtSubstituteCharacters(marker)}</span>
                <input type="text" class="edit-input" style="display:none;">
                <button type="button" class="edit-button" onclick="adtEditRow(event, this)">Edit</button>
                <button type="button" class="save-button" style="display:none;" onclick="adtSaveRow(event, this)">Save</button>
            </div>
        </td>
        <td>${catalogue_number}</td>
        <td>${clone}</td>
        <td>${reactivity}</td>
        <td>${barcode}</td>
        <td><button type="button" class="remove-button" onclick="adtRemoveRow(this)">Remove</button></td>
    `;

    tableBody.appendChild(row);
}

function adtGenerateCSV() {
    const outputFormat = document.getElementById('output-format').value;
    const format = document.getElementById('format').value;
    const csvNameInput = document.getElementById('csv-name-input').value.trim();
    
    // Validate all required fields
    if (!outputFormat) {
        adtShowMessage('Select a CSV output format.');
        return;
    }
    if (!format) {
        adtShowMessage('Select a TotalSeq panel.');
        return;
    }
    if (!csvNameInput) {
        adtShowMessage('Enter a file name.');
        return;
    }
    if (outputFormat === 'tapestri' && !format.includes('totalseq_d')) {
        adtShowMessage('Tapestri output requires TotalSeq-D.');
        return;
    }
    
    // Sanitise filename
    const csvName = csvNameInput.replace(/[^a-zA-Z0-9_-]/g, '_');

    const rows = Array.from(document.querySelectorAll('#rowsContainer tr'));
    if (rows.length === 0) {
        adtShowMessage('Add at least one marker before downloading.');
        return;
    }
    let csvContent = '';

    const hashtagRows = rows.filter(row => row.querySelector('td').textContent.includes('Hashtag'));
    const nonHashtagRows = rows.filter(row => !row.querySelector('td').textContent.includes('Hashtag'));

    nonHashtagRows.sort((a, b) => {
        const markerA = a.querySelector('td').textContent.toLowerCase();
        const markerB = b.querySelector('td').textContent.toLowerCase();
        return markerA.localeCompare(markerB);
    });

    const sortedRows = nonHashtagRows.concat(hashtagRows);

    if (outputFormat === 'cellranger') {
        csvContent = 'id,name,read,pattern,sequence,feature_type\n';
        sortedRows.forEach(row => {
            const cells = row.querySelectorAll('td');
            const marker = adtSubstituteCharacters(cells[1].textContent);
            const correctedName = cells[2].querySelector('.corrected-name').textContent;
            const barcode = cells[6].textContent.trim();  // Remove leading/trailing spaces
            let pattern = '';

            if (format.includes('totalseq_a')) {
                pattern = '5P(BC)';
            } else if (format.includes('totalseq_b') || format.includes('totalseq_c')) {
                pattern = '5PNNNNNNNNNN(BC)';
            }

            csvContent += `${correctedName},${correctedName},R2,${pattern},${barcode},Antibody Capture\n`;
        });
    } else if (outputFormat === 'kallisto') {
        csvContent = 'Feature Barcode name,Feature Barcode sequence\n';
        sortedRows.forEach(row => {
            const cells = row.querySelectorAll('td');
            const correctedName = cells[2].querySelector('.corrected-name').textContent;
            const barcode = cells[6].textContent.trim();  // Remove leading/trailing spaces
            csvContent += `${correctedName},${barcode}\n`;
        });
    } else if (outputFormat === 'tapestri' && format.includes('totalseq_d')) {
        csvContent = 'ID,Name,Sequence\n';
        sortedRows.forEach(row => {
            const cells = row.querySelectorAll('td');
            const correctedName = cells[2].querySelector('.corrected-name').textContent.trim(); // Ensure correct cell and trim
            const totalseq_id = cells[0].textContent.trim();  // Remove leading/trailing spaces
            const barcode = cells[6].textContent.trim();  // Remove leading/trailing spaces
    
            const new_id = `D${totalseq_id}`;
            csvContent += `${new_id},${correctedName},${barcode}\n`;
        });
    }

    // Add UTF-8 BOM for Excel compatibility
    const BOM = '\uFEFF';
    const blob = new Blob([BOM + csvContent], { type: 'text/csv;charset=utf-8;' });
    const link = document.createElement('a');
    const url = URL.createObjectURL(blob);
    link.setAttribute('href', url);
    link.setAttribute('download', `${csvName}.csv`);
    link.style.visibility = 'hidden';
    document.body.appendChild(link);
    link.click();
    document.body.removeChild(link);
    URL.revokeObjectURL(url);
    adtShowMessage(`${csvName}.csv downloaded.`, 'success');
}
