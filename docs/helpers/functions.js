// functions.js

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

function pagesIncludeTopBar() {
    fetch('../helpers/topbar.html')
        .then(response => response.text())
        .then(data => {
            document.getElementById('topbar-placeholder').innerHTML = data;
            document.getElementById('home-link').href = 'index.html';
            document.getElementById('metadata-link').href = 'metadata_generator.html';
            document.getElementById('adt-link').href = 'adt_generator.html';
            document.getElementById('functions-link').href = 'oscar_functions.html';
            document.getElementById('references-link').href = 'references.html';
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
