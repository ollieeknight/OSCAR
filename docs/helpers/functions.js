// functions.js
function includeTopBar() {
    fetch('/docs/helpers/topbar.html')
        .then(response => response.text())
        .then(data => {
            document.getElementById('topbar-placeholder').innerHTML = data;
            // Set the href attributes
            document.getElementById('home-link').href = '/docs/index.html';
            document.getElementById('metadata-link').href = '/docs/pages/metadata_generator.html';
            document.getElementById('adt-link').href = '/docs/pages/adt_generator.html';
            document.getElementById('references-link').href = '/docs/pages/references.html';
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