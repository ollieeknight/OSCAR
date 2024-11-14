// functions.js

// Top bar for index page (index.html)
function indexIncludeTopBar() {
    fetch('./helpers/topbar.html')
        .then(response => response.text())
        .then(data => {
            document.getElementById('topbar-placeholder').innerHTML = data;
            document.getElementById('home-link').href = 'index.html';
            document.getElementById('metadata-link').href = 'pages/metadata_generator.html';
            document.getElementById('adt-link').href = 'pages/adt_generator.html';
            document.getElementById('functions-link').href = 'pages/oscar_functions.html';
            document.getElementById('references-link').href = 'pages/references.html';
        })
        .catch(error => {
            console.error('Error fetching topbar:', error);
        });
}

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

// Top bar for pages within pages folder 
function pagesIncludeTopBar() {
    fetch('../helpers/topbar.html')
        .then(response => response.text())
        .then(data => {
            document.getElementById('topbar-placeholder').innerHTML = data;
            document.getElementById('home-link').href = '../index.html';
            document.getElementById('metadata-link').href = 'metadata_generator.html';
            document.getElementById('adt-link').href = 'adt_generator.html';
            document.getElementById('functions-link').href = 'oscar_functions.html';
            document.getElementById('references-link').href = 'references.html';
        })
        .catch(error => {
            console.error('Error fetching topbar:', error);
        });
}

// Functions for metadata_generator.html

// Function to fetch the row template from an external HTML file
async function metadataFetchRowTemplate() {
    try {
        const response = await fetch('../helpers/row_template.html');
        if (!response.ok) throw new Error('Failed to load row template');
        const template = await response.text();
        return template;
    } catch (error) {
        console.error('Error fetching row template:', error);
        return '';
    }
}

// Function to add a new row to the table
async function metadataAddRow() {
    const container = document.getElementById('rowsContainer');
    const template = await metadataFetchRowTemplate();
    if (template) {
        const newRow = document.createElement('tr');
        newRow.innerHTML = template;
        container.appendChild(newRow);
    } else {
        alert('Could not add row: template not loaded.');
    }
}

// Function to validate a row and collect its values
function metadataValidateRow(row) {
    const inputs = row.querySelectorAll('input, select');
    const values = Array.from(inputs).map(input => input.value || 'NA');
    let rowValid = true;

    inputs.forEach(input => {
        if (input.value === '') {
            input.classList.add('error'); // Add error class if input is empty
            rowValid = false;
        } else {
            input.classList.remove('error'); // Remove error class if input is not empty
        }
    });

    return { rowValid, values }; // Return validation status and values
}

// Function to download the generated CSV content
function metadataDownloadCSV(csvContent) {
    const blob = new Blob([csvContent], { type: 'text/csv;charset=utf-8;' });
    const url = URL.createObjectURL(blob);
    const a = document.createElement('a');
    a.href = url;
    a.download = 'metadata.csv'; // Output file name
    document.body.appendChild(a);
    a.click();
    document.body.removeChild(a);
    URL.revokeObjectURL(url); // Clean up the URL object
}

// Function to generate the CSV content from the table rows
function metadataGenerateCSV() {
    const rows = document.querySelectorAll('#rowsContainer tr');
    let csvContent = 'assay,experiment_id,historical_number,replicate,modality,chemistry,index_type,index,species,n_donors,adt_file\n';
    let allValid = true;

    rows.forEach(row => {
        const { rowValid, values } = metadataValidateRow(row);
        if (!rowValid) {
            allValid = false;
        } else {
            csvContent += values.join(',') + '\n';
        }
    });

    if (allValid) {
        metadataDownloadCSV(csvContent); // Download the CSV if all rows are valid
    } else {
        alert('Please complete all fields before generating the CSV.'); // Alert if any row is invalid
    }
}

// Functions for adt_generator.html

async function adtFetchCSV(file) {
    const response = await fetch(file);
    const data = await response.text();
    return data;
}

function adtParseCSV(data) {
    const rows = data.split('\n').slice(1);
    return rows.map(row => {
        const [catalogue_number, totalseq_id, marker, clone, reactivity, barcode] = row.split(',');
        return { catalogue_number, totalseq_id, marker, clone, reactivity, barcode };
    });
}

async function adtFilterMarkers() {
    const format = document.getElementById('format').value;
    const species = document.getElementById('species').value;
    const data = await adtFetchCSV(format);
    markersData = adtParseCSV(data);

    if (species === 'Human') {
        markersData = markersData.filter(marker => ['Human', 'Human and mouse', 'Isotype control', 'Epitope'].includes(marker.reactivity));
    } else if (species === 'Mouse') {
        markersData = markersData.filter(marker => ['Mouse', 'Human and mouse', 'Isotype control', 'Epitope'].includes(marker.reactivity));
    }

    // Call adtShowDropdown after data is filtered
    adtShowDropdown(); 
}

function adtAddRow(name, totalseq_id, catalogueNumber, clone, reactivity, barcode) {
    const container = document.getElementById('rowsContainer');
    const correctedName = adtSubstituteCharacters(name);
    const newRow = document.createElement('tr');
    newRow.innerHTML = `
        <td>${name}</td>
        <td>
            <div class="corrected-name-container">
                <span class="corrected-name">${correctedName}</span>
                <input type="text" class="edit-input" value="${correctedName}" style="display:none;">
                <button class="edit-button" onclick="adtEditRow(event, this)">Edit</button>
                <button class="save-btn" onclick="adtSaveRow(event, this)" style="display:none;">Save</button>
            </div>
        </td>
        <td>${totalseq_id}</td>
        <td>${catalogueNumber}</td>
        <td>${clone}</td>
        <td>${reactivity}</td>
        <td>${barcode}</td>
        <td><span class="remove-button" onclick="adtRemoveRow(this)">Remove</span></td>`;
    container.appendChild(newRow);
}

function adtRemoveRow(element) {
    const row = element.parentElement.parentElement;
    row.remove();
}

function adtShowDropdown() {
    const search = document.getElementById('search').value.toLowerCase();
    const dropdown = document.getElementById('dropdown-content');
    dropdown.innerHTML = '';

    if (search && markersData.length > 0) {
        const filteredMarkers = markersData.filter(marker => 
            marker.catalogue_number.toLowerCase().includes(search) ||
            marker.marker.toLowerCase().includes(search) ||
            marker.clone.toLowerCase().includes(search)
        );
        let maxWidth = 0;
        filteredMarkers.forEach(marker => {
            const div = document.createElement('div');
            const text = `${marker.marker}, clone ${marker.clone}, catalogue number ${marker.catalogue_number}`;
            const highlightedText = text.replace(new RegExp(search, 'gi'), match => `<strong>${match}</strong>`);
            div.innerHTML = highlightedText;
            div.onclick = () => {
                adtAddRow(marker.marker, marker.totalseq_id, marker.catalogue_number, marker.clone, marker.reactivity, marker.barcode);
                document.getElementById('search').value = '';
                dropdown.innerHTML = '';
            };
            dropdown.appendChild(div);
            const tempSpan = document.createElement('span');
            tempSpan.style.visibility = 'hidden';
            tempSpan.style.position = 'absolute';
            tempSpan.innerHTML = highlightedText;
            document.body.appendChild(tempSpan);
            const width = tempSpan.offsetWidth;
            document.body.removeChild(tempSpan);
            if (width > maxWidth) {
                maxWidth = width;
            }
        });
        dropdown.style.display = 'block';
        dropdown.style.width = `${maxWidth}px`;
    } else {
        dropdown.style.display = 'none';
    }
}

function adtConfirmReset(type) {
    const speciesElement = document.getElementById('species');
    const formatElement = document.getElementById('format');
    const outputFormatElement = document.getElementById('output-format');
    const table = document.getElementById('csvTable'); // Updated ID

    let currentValue;
    if (type === 'species') {
        currentValue = speciesElement.value;
        if (!speciesSelected) {
            speciesSelected = true;
            adtFilterMarkers();
            return;
        }
    } else if (type === 'format') {
        currentValue = formatElement.value;
        if (!formatSelected) {
            formatSelected = true;
            adtFilterMarkers();
            return;
        }
    } else if (type === 'output-format') {
        currentValue = outputFormatElement.value;
        if (!outputFormatSelected) {
            outputFormatSelected = true;
            return;
        }
    }

    const hasEntries = table.getElementsByTagName('tr').length > 1;

    if (currentValue !== '' && hasEntries) {
        const confirmed = confirm('Changing this will reset all choices. Is this alright?');
        if (!confirmed) {
            // Revert the selection to the previous value
            if (type === 'species') {
                speciesElement.value = '';
                speciesSelected = false;
            } else if (type === 'format') {
                formatElement.value = '';
                formatSelected = false;
            } else if (type === 'output-format') {
                outputFormatElement.value = '';
                outputFormatSelected = false;
            }
            return;
        } else {
            // Reset the form or perform necessary actions
            resetForm();
        }
    }
}

function resetForm() {
    // Clear all selections
    document.getElementById('species').value = '';
    document.getElementById('format').value = '';
    document.getElementById('output-format').value = '';
    
    // Clear the table rows
    const rowsContainer = document.getElementById('rowsContainer');
    rowsContainer.innerHTML = '';
    
    // Reset flags
    speciesSelected = false;
    formatSelected = false;
    outputFormatSelected = false;
}

function adtSubstituteCharacters(text) {
    return text
        .replace(/alpha/gi, 'a')
        .replace(/beta/gi, 'b')
        .replace(/gamma/gi, 'g')
        .replace(/delta/gi, 'd')
        .replace(/γ/gi, 'g')
        .replace(/δ/gi, 'd')
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
        <td>${marker}</td>
        <td>
            <div class="edit-container">
                <span class="corrected-name">${marker}</span>
                <input type="text" class="edit-input" style="display:none;">
                <button class="edit-button" onclick="adtEditRow(event, this)">Edit</button>
                <button class="save-button" style="display:none;" onclick="adtSaveRow(event, this)">Save</button>
            </div>
        </td>
        <td>${totalseq_id}</td>
        <td>${catalogue_number}</td>
        <td>${clone}</td>
        <td>${reactivity}</td>
        <td>${barcode}</td>
        <td><button class="remove-button" onclick="adtRemoveRow(this)">Remove</button></td>
    `;

    tableBody.appendChild(row);
}

function adtGenerateCSV() {
    const outputFormat = document.getElementById('output-format').value;
    if (!outputFormat) {
        alert('Please select an output format.');
        return;
    }

    const rows = Array.from(document.querySelectorAll('#rowsContainer tr'));
    const csvName = document.querySelector('.csv-name-input').value || 'adt_list';
    let csvContent = '';

    const hashtagRows = rows.filter(row => row.querySelector('td').textContent.includes('Hashtag'));
    const nonHashtagRows = rows.filter(row => !row.querySelector('td').textContent.includes('Hashtag'));

    nonHashtagRows.sort((a, b) => {
        const markerA = a.querySelector('td').textContent.toLowerCase();
        const markerB = b.querySelector('td').textContent.toLowerCase();
        return markerA.localeCompare(markerB);
    });

    const sortedRows = nonHashtagRows.concat(hashtagRows);

    const format = document.getElementById('format').value;

    if (outputFormat === 'cellranger') {
        csvContent = 'id,name,read,pattern,sequence,feature_type\n';
        sortedRows.forEach(row => {
            const cells = row.querySelectorAll('td');
            const marker = adtSubstituteCharacters(cells[0].textContent);
            const correctedName = cells[1].querySelector('.corrected-name').textContent;
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
            const correctedName = cells[1].querySelector('.corrected-name').textContent;
            const barcode = cells[6].textContent.trim();  // Remove leading/trailing spaces
            csvContent += `${correctedName},${barcode}\n`;
        });
    } else if (outputFormat === 'tapestri' && format.includes('totalseq_d')) {
        csvContent = 'ID,Name,Sequence\n';
        sortedRows.forEach(row => {
            const cells = row.querySelectorAll('td');
            const correctedName = cells[1].querySelector('.corrected-name').textContent.trim(); // Ensure correct cell and trim
            const totalseq_id = cells[2].textContent.trim();  // Remove leading/trailing spaces
            const barcode = cells[6].textContent.trim();  // Remove leading/trailing spaces
    
            const new_id = `D${totalseq_id}`;
            csvContent += `${new_id},${correctedName},${barcode}\n`;
        });
    }

    const blob = new Blob([csvContent], { type: 'text/csv;charset=utf-8;' });
    const link = document.createElement('a');
    const url = URL.createObjectURL(blob);
    link.setAttribute('href', url);
    link.setAttribute('download', `${csvName}.csv`);
    link.style.visibility = 'hidden';
    document.body.appendChild(link);
    link.click();
    document.body.removeChild(link);
}