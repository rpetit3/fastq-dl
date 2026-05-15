---
name: update-catalog
description: Regenerate catalog.json and llms.txt from the fastq-dl source code using the AST-based catalog generation script.
---

# Update Catalog

Regenerates `catalog.json` and `llms.txt` by running `.claude/skills/update-catalog/scripts/update_catalog.py`,
which uses AST parsing to extract metadata from the fastq-dl source code.

## Steps

1. **Check for unsaved changes** to `catalog.json` and `llms.txt`:
   ```bash
   git diff --name-only catalog.json llms.txt
   git diff --cached --name-only catalog.json llms.txt
   ```
   If either file has uncommitted changes, warn the user before proceeding. The script will
   overwrite both files.

2. **Run the catalog generation script**:
   ```bash
   python .claude/skills/update-catalog/scripts/update_catalog.py
   ```
   This parses all Python source files under `fastq_dl/`, extracts module metadata, function
   signatures, CLI options, dependencies, constants, and exception hierarchy, then writes:
   - `catalog.json` — machine-readable project metadata
   - `llms.txt` — AI-discovery document (markdown)

3. **Validate the output**:
   ```bash
   python -m json.tool catalog.json > /dev/null
   ```

4. **Show a summary of changes**:
   ```bash
   git diff --stat catalog.json llms.txt
   ```

5. **Optionally update CLAUDE.md**: If the script output indicates structural changes (new
   modules, new functions, changed exception hierarchy), review `CLAUDE.md` and update the
   relevant sections to match.

## When to Run

- After adding, removing, or renaming functions or modules
- After modifying CLI options (adding/removing Click decorators)
- After changing the exception hierarchy
- After updating dependencies in `pyproject.toml`
- Before a release to ensure documentation is current
