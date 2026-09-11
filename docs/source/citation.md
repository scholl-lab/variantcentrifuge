# Citation and release archiving

VariantCentrifuge's `CITATION.cff` supplies the software citation shown by GitHub.
The `.zenodo.json` file supplies matching metadata for a future Zenodo archive.
No Zenodo DOI is currently available.

## Cite the version used

Record the release tag and, when available, its version-specific archive DOI.
Before archival, use the GitHub release URL. For development code, record the
full commit hash and its GitHub commit URL. An archive's concept DOI identifies
all versions; use the version-specific DOI when reproducibility requires a
particular release.

The current source version is **0.17.13**, which is not yet released. The latest
published release is [v0.17.12 (2026-07-19)](https://github.com/scholl-lab/variantcentrifuge/releases/tag/v0.17.12).
Citation metadata follows the current source version and omits release dates
until that version is released.

## Maintainer steps for the next release

1. Confirm the release contents and creator metadata. Set the same version in
   `variantcentrifuge/version.py`, `CITATION.cff`, and `.zenodo.json`. Add the
   actual release date as `date-released` in the CFF and `publication_date` in
   the JSON. Keep the title, description, authors, ORCIDs, affiliations,
   keywords, and MIT licence consistent across both files.
2. Validate the files before tagging:

   ```bash
   python -m pip install cffconvert
   cffconvert --validate
   python -m json.tool .zenodo.json > /dev/null
   ```

   Preserve UTF-8 text when editing names and affiliations. Zenodo gives
   `.zenodo.json` precedence when both metadata files exist, so update both.
3. A repository maintainer must connect GitHub to Zenodo and enable this
   repository, then publish the intended release following
   [Zenodo's GitHub archival guide](https://help.zenodo.org/docs/github/archive-software/github-upload/).
   Adding metadata alone does not create an archive. Existing tags and releases
   remain unchanged; to archive an older release, upload that exact release's
   source archive manually with metadata matching that version.
4. Inspect the published Zenodo record: verify the archived version, creators,
   licence, date, and repository link. Copy the real DOI and DOI badge supplied
   by Zenodo into the README, replacing the pending status. Record the
   version-specific DOI in the citation metadata for that release; do not carry
   it forward as the DOI of a different release.
5. For later releases, repeat the metadata checks and verify the new archive
   record. Keep the README citation guidance current.

These are manual maintainer actions. This guide does not activate an integration
or publish a release.
