#!/usr/bin/env python3
"""Generate catalog.json and llms.txt from the fastq-dl source code.

Uses AST parsing to extract module metadata, function signatures, CLI options,
dependencies, constants, and exception hierarchy. Outputs:
  - catalog.json: machine-readable project metadata
  - llms.txt: AI-discovery document (markdown)

Usage:
    python .claude/skills/update-catalog/scripts/update_catalog.py
"""
import ast
import json
from datetime import datetime, timezone
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent.parent.parent
PACKAGE_DIR = PROJECT_ROOT / "fastq_dl"
CATALOG_PATH = PROJECT_ROOT / "catalog.json"
LLMS_PATH = PROJECT_ROOT / "llms.txt"


def parse_pyproject() -> dict:
    """Extract project metadata and dependencies from pyproject.toml."""
    pyproject_path = PROJECT_ROOT / "pyproject.toml"
    text = pyproject_path.read_text()

    # Simple TOML parsing for the fields we need (avoids tomllib dep on 3.10)
    meta = {}

    for line in text.splitlines():
        line = line.strip()
        if line.startswith("name ="):
            meta["name"] = line.split("=", 1)[1].strip().strip('"')
        elif line.startswith("version ="):
            meta["version"] = line.split("=", 1)[1].strip().strip('"')
        elif line.startswith("description ="):
            meta["description"] = line.split("=", 1)[1].strip().strip('"')
        elif line.startswith("license ="):
            meta["license"] = line.split("=", 1)[1].strip().strip('"')
        elif line.startswith("homepage ="):
            meta["homepage"] = line.split("=", 1)[1].strip().strip('"')
        elif line.startswith("repository ="):
            meta["repository"] = line.split("=", 1)[1].strip().strip('"')

    # Parse dependencies
    runtime_deps = {}
    dev_deps = {}
    in_runtime = False
    in_dev = False

    for line in text.splitlines():
        stripped = line.strip()
        if stripped == "[tool.poetry.dependencies]":
            in_runtime = True
            in_dev = False
            continue
        elif stripped == "[tool.poetry.group.dev.dependencies]":
            in_dev = True
            in_runtime = False
            continue
        elif stripped.startswith("["):
            in_runtime = False
            in_dev = False
            continue

        if (in_runtime or in_dev) and "=" in stripped and not stripped.startswith("#"):
            key, val = stripped.split("=", 1)
            key = key.strip()
            val = val.strip().strip('"')
            if key == "python":
                meta["python_requires"] = val
            elif in_runtime:
                runtime_deps[key] = val
            elif in_dev:
                dev_deps[key] = val

    meta["dependencies"] = {"runtime": runtime_deps, "dev": dev_deps}
    return meta


def get_docstring(node) -> str:
    """Extract docstring from an AST node."""
    ds = ast.get_docstring(node)
    return ds.strip() if ds else ""


def parse_docstring_args(docstring: str) -> dict[str, str]:
    """Parse the Args/Parameters section of a docstring to extract per-arg descriptions.

    Handles both NumPy-style ("Args:" or "Parameters:") docstrings with entries like:
        arg_name (type): Description text
        arg_name (type, optional): Description text
        arg_name: Description text
    Continuation lines (indented further) are appended to the previous arg.
    """
    if not docstring:
        return {}

    arg_descriptions = {}
    in_args_section = False
    current_arg = None
    args_indent = None

    for line in docstring.splitlines():
        stripped = line.strip()

        # Detect start of Args/Parameters section
        if stripped in ("Args:", "Parameters:", "Arguments:"):
            in_args_section = True
            args_indent = None
            continue

        # Detect end of section (new section header like "Returns:", "Raises:", etc.)
        if in_args_section and stripped and stripped.endswith(":") and not stripped.startswith("-"):
            word = stripped.rstrip(":")
            if word and word[0].isupper() and " " not in word:
                in_args_section = False
                current_arg = None
                continue

        if not in_args_section:
            continue

        if not stripped:
            continue

        # Determine indentation level
        leading = len(line) - len(line.lstrip())

        # First arg line sets the baseline indent
        if args_indent is None:
            args_indent = leading

        if leading == args_indent:
            # New argument line: "arg_name (type): description" or "arg_name: description"
            if ":" in stripped:
                # Split on first colon that's outside parentheses
                paren_depth = 0
                colon_pos = None
                for i, ch in enumerate(stripped):
                    if ch == "(":
                        paren_depth += 1
                    elif ch == ")":
                        paren_depth -= 1
                    elif ch == ":" and paren_depth == 0:
                        colon_pos = i
                        break

                if colon_pos is not None:
                    name_part = stripped[:colon_pos].strip()
                    desc_part = stripped[colon_pos + 1:].strip()
                    # Extract just the arg name (strip type annotation in parens)
                    arg_name = name_part.split("(")[0].strip().split(",")[0].strip()
                    current_arg = arg_name
                    arg_descriptions[current_arg] = desc_part
        elif leading > args_indent and current_arg:
            # Continuation line for current arg
            arg_descriptions[current_arg] += " " + stripped

    return arg_descriptions


def extract_function_info(func_node: ast.FunctionDef) -> dict:
    """Extract function metadata from an AST FunctionDef node."""
    docstring = get_docstring(func_node)
    doc_args = parse_docstring_args(docstring)

    args = []
    for arg in func_node.args.args:
        arg_info = {"name": arg.arg}
        if arg.annotation:
            arg_info["annotation"] = ast.unparse(arg.annotation)
        if arg.arg in doc_args:
            arg_info["description"] = doc_args[arg.arg]
        args.append(arg_info)

    # Defaults (aligned to end of args list)
    defaults = func_node.args.defaults
    if defaults:
        offset = len(args) - len(defaults)
        for i, default in enumerate(defaults):
            try:
                args[offset + i]["default"] = ast.unparse(default)
            except Exception:
                args[offset + i]["default"] = "..."

    returns = None
    if func_node.returns:
        returns = ast.unparse(func_node.returns)

    return {
        "name": func_node.name,
        "description": docstring,
        "args": args,
        "returns": returns,
        "lineno": func_node.lineno,
    }


def extract_class_info(class_node: ast.ClassDef) -> dict:
    """Extract class metadata from an AST ClassDef node."""
    bases = [ast.unparse(b) for b in class_node.bases]

    # Extract __init__ args for exception classes
    init_args = []
    for item in class_node.body:
        if isinstance(item, ast.FunctionDef) and item.name == "__init__":
            for arg in item.args.args:
                if arg.arg != "self":
                    arg_info = {"name": arg.arg}
                    if arg.annotation:
                        arg_info["annotation"] = ast.unparse(arg.annotation)
                    init_args.append(arg_info)
            # Defaults
            defaults = item.args.defaults
            if defaults:
                offset = len(init_args) - len(defaults)
                for i, default in enumerate(defaults):
                    try:
                        init_args[offset + i]["default"] = ast.unparse(default)
                    except Exception:
                        init_args[offset + i]["default"] = "..."
            break

    return {
        "name": class_node.name,
        "description": get_docstring(class_node),
        "bases": bases,
        "init_args": init_args,
        "lineno": class_node.lineno,
    }


def extract_constants(tree: ast.Module) -> dict:
    """Extract module-level constant assignments."""
    constants = {}
    for node in ast.iter_child_nodes(tree):
        if isinstance(node, ast.Assign):
            for target in node.targets:
                if isinstance(target, ast.Name) and target.id.isupper():
                    try:
                        constants[target.id] = ast.literal_eval(node.value)
                    except (ValueError, TypeError):
                        constants[target.id] = ast.unparse(node.value)
    return constants


def parse_module(filepath: Path) -> dict:
    """Parse a Python module and extract its metadata."""
    source = filepath.read_text()
    tree = ast.parse(source, filename=str(filepath))

    rel_path = filepath.relative_to(PROJECT_ROOT)
    module_name = (
        str(rel_path).replace("/", ".").replace("\\", ".").removesuffix(".py")
    )

    module_info = {
        "path": str(rel_path),
        "module": module_name,
        "description": get_docstring(tree),
        "functions": {},
        "classes": {},
        "constants": {},
    }

    for node in ast.iter_child_nodes(tree):
        if isinstance(node, ast.FunctionDef):
            func_info = extract_function_info(node)
            module_info["functions"][node.name] = func_info
        elif isinstance(node, ast.ClassDef):
            class_info = extract_class_info(node)
            module_info["classes"][node.name] = class_info

    # Extract constants from constants.py
    if filepath.name == "constants.py":
        module_info["constants"] = extract_constants(tree)

    return module_info


def extract_click_options(filepath: Path) -> list:
    """Extract Click option definitions from the CLI module."""
    source = filepath.read_text()
    tree = ast.parse(source, filename=str(filepath))

    options = []
    for node in ast.walk(tree):
        if not isinstance(node, ast.Call):
            continue
        if not isinstance(node.func, ast.Attribute):
            continue
        if node.func.attr != "option":
            continue

        option = {}
        # First positional args are the option names
        names = []
        for arg in node.args:
            if isinstance(arg, ast.Constant) and isinstance(arg.value, str):
                names.append(arg.value)

        if not names:
            continue

        # Find the long name and short name
        long_names = [n for n in names if n.startswith("--")]
        short_names = [n for n in names if n.startswith("-") and not n.startswith("--")]
        option["name"] = long_names[0] if long_names else names[0]
        if short_names:
            option["short"] = short_names[0]

        for kw in node.keywords:
            if kw.arg == "default":
                try:
                    option["default"] = ast.literal_eval(kw.value)
                except (ValueError, TypeError):
                    option["default"] = ast.unparse(kw.value)
            elif kw.arg == "help":
                if isinstance(kw.value, ast.Constant):
                    option["help"] = kw.value.value
            elif kw.arg == "required":
                if isinstance(kw.value, ast.Constant):
                    option["required"] = kw.value.value
            elif kw.arg == "is_flag":
                if isinstance(kw.value, ast.Constant):
                    option["is_flag"] = kw.value.value
            elif kw.arg == "type":
                try:
                    option["type"] = ast.unparse(kw.value)
                except Exception:
                    pass

        options.append(option)

    return options


def build_catalog() -> dict:
    """Build the complete catalog from source code."""
    meta = parse_pyproject()

    # Parse all Python modules
    modules = {}
    for py_file in sorted(PACKAGE_DIR.rglob("*.py")):
        if py_file.name == "__init__.py" and py_file.parent != PACKAGE_DIR:
            continue
        module_info = parse_module(py_file)
        modules[module_info["module"]] = module_info

    # Extract CLI options
    cli_path = PACKAGE_DIR / "cli" / "download.py"
    cli_options = extract_click_options(cli_path) if cli_path.exists() else []

    # Build a lookup from CLI option parameter names to their help/default
    # Click option names like "--group-by-experiment" map to Python arg "group_by_experiment"
    cli_param_lookup = {}
    for opt in cli_options:
        param_name = opt["name"].lstrip("-").replace("-", "_")
        # Handle renamed params (e.g. --ignore -> ignore_md5 via Click's parameter name)
        cli_param_lookup[param_name] = {
            "description": opt.get("help", ""),
            "default": opt.get("default"),
            "is_flag": opt.get("is_flag", False),
            "required": opt.get("required", False),
        }
    # Manual mappings for params whose Click name differs from the Python arg name
    if "ignore" in cli_param_lookup:
        cli_param_lookup["ignore_md5"] = cli_param_lookup["ignore"]

    # Enrich functions that receive CLI params but lack docstring-level arg descriptions
    cli_mod = modules.get("fastq_dl.cli.download", {})
    for func_name in ("_run_download", "_download_with_fallback", "fastqdl"):
        func_info = cli_mod.get("functions", {}).get(func_name)
        if not func_info:
            continue
        for arg in func_info["args"]:
            if "description" not in arg and arg["name"] in cli_param_lookup:
                info = cli_param_lookup[arg["name"]]
                if info["description"]:
                    arg["description"] = info["description"]
                if "default" not in arg and info["default"] is not None:
                    arg["default"] = str(info["default"])
                if "default" not in arg and info["is_flag"]:
                    arg["default"] = "False"

    # Build exception hierarchy from the exceptions module
    exceptions = {}
    exc_module = modules.get("fastq_dl.exceptions", {})
    for cls_name, cls_info in exc_module.get("classes", {}).items():
        exceptions[cls_name] = {
            "description": cls_info["description"],
            "bases": cls_info["bases"],
            "init_args": cls_info["init_args"],
        }

    # Extract constants
    constants = {}
    const_module = modules.get("fastq_dl.constants", {})
    constants = const_module.get("constants", {})

    catalog = {
        "name": meta.get("name", "fastq-dl"),
        "version": meta.get("version", "unknown"),
        "description": meta.get("description", ""),
        "repository": meta.get("repository", ""),
        "homepage": meta.get("homepage", ""),
        "license": meta.get("license", ""),
        "python_requires": meta.get("python_requires", ""),
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "modules": modules,
        "cli_options": cli_options,
        "dependencies": meta.get("dependencies", {}),
        "external_tools": [
            {"name": "wget", "purpose": "ENA FTP/HTTPS downloads"},
            {
                "name": "prefetch",
                "purpose": "SRA file download",
                "package": "sra-tools",
            },
            {
                "name": "fasterq-dump",
                "purpose": "SRA to FASTQ conversion",
                "package": "sra-tools",
            },
            {
                "name": "vdb-config",
                "purpose": "SRA configuration",
                "package": "sra-tools",
            },
            {"name": "pigz", "purpose": "Parallel gzip compression"},
        ],
        "exceptions": exceptions,
        "constants": constants,
    }

    return catalog


def generate_llms_txt(catalog: dict) -> str:
    """Generate llms.txt content from the catalog."""
    version = catalog["version"]
    description = catalog["description"]
    repo = catalog["repository"]

    # Count functions across all modules
    total_functions = sum(
        len(m.get("functions", {})) for m in catalog["modules"].values()
    )
    total_classes = sum(len(m.get("classes", {})) for m in catalog["modules"].values())

    # Build module index
    module_lines = []
    for mod_name, mod_info in sorted(catalog["modules"].items()):
        path = mod_info["path"]
        desc = mod_info.get("description", "").split("\n")[0] if mod_info.get("description") else ""
        funcs = list(mod_info.get("functions", {}).keys())
        classes = list(mod_info.get("classes", {}).keys())

        module_lines.append(f"### {mod_name}")
        module_lines.append(f"- **Path**: `{path}`")
        if desc:
            module_lines.append(f"- **Description**: {desc}")
        if funcs:
            module_lines.append(f"- **Functions**: {', '.join(f'`{f}`' for f in funcs)}")
        if classes:
            module_lines.append(f"- **Classes**: {', '.join(f'`{c}`' for c in classes)}")
        module_lines.append("")

    # Build CLI options table
    cli_lines = []
    for opt in catalog.get("cli_options", []):
        name = opt["name"]
        short = opt.get("short", "")
        help_text = opt.get("help", "")
        default = opt.get("default", "")
        required = opt.get("required", False)
        is_flag = opt.get("is_flag", False)

        if required:
            default_str = "required"
        elif is_flag:
            default_str = "off"
        elif default is not None and default != "":
            default_str = str(default)
        else:
            default_str = ""

        short_str = f" (`{short}`)" if short else ""
        cli_lines.append(f"| `{name}`{short_str} | {default_str} | {help_text} |")

    # Build exception list
    exc_lines = []
    for exc_name, exc_info in catalog.get("exceptions", {}).items():
        desc = exc_info.get("description", "")
        args = exc_info.get("init_args", [])
        arg_str = ""
        if args:
            arg_names = [a["name"] for a in args]
            arg_str = f" (attrs: {', '.join(arg_names)})"
        exc_lines.append(f"- **{exc_name}**: {desc}{arg_str}")

    # Build external tools list
    tool_lines = []
    for tool in catalog.get("external_tools", []):
        pkg = f" (from {tool['package']})" if "package" in tool else ""
        tool_lines.append(f"- `{tool['name']}` — {tool['purpose']}{pkg}")

    sections = [
        f"# fastq-dl",
        f"",
        f"> {description}",
        f"",
        f"- **Version**: {version}",
        f"- **Repository**: {repo}",
        f"- **License**: {catalog['license']}",
        f"- **Python**: {catalog['python_requires']}",
        f"- **Modules**: {len(catalog['modules'])} | **Functions**: {total_functions} | **Classes**: {total_classes}",
        f"",
        f"## Architecture",
        f"",
        f"fastq-dl follows a 3-layer architecture:",
        f"",
        f"1. **CLI Layer** (`fastq_dl/cli/download.py`): Click-based command-line interface that parses",
        f"   arguments, sets up logging, and orchestrates the download workflow.",
        f"",
        f"2. **Provider Layer** (`fastq_dl/providers/`): Handles metadata queries and file downloads from",
        f"   ENA and SRA. Includes automatic fallback between providers with configurable retry logic.",
        f"   - `generic.py` — Provider coordination, fallback orchestration",
        f"   - `ena.py` — ENA Data Warehouse API queries, wget-based FTP/HTTPS downloads",
        f"   - `sra.py` — SRA queries via pysradb, prefetch + fasterq-dump downloads",
        f"",
        f"3. **Utilities Layer** (`fastq_dl/utils.py`): Subprocess execution with retries, MD5 checksums,",
        f"   FASTQ merging, accession validation, and TSV output.",
        f"",
        f"### Data Flow",
        f"",
        f"```",
        f"CLI (accession input)",
        f" → validate_query() → accession type regex check",
        f" → get_run_info()   → query primary provider, fallback to secondary",
        f" → _download_with_fallback() → download FASTQs per run",
        f" → merge_runs() (optional) → concatenate grouped FASTQs",
        f" → write_tsv() → output metadata files",
        f"```",
        f"",
        f"## Module Index",
        f"",
    ]
    sections.extend(module_lines)
    sections.extend([
        f"## CLI Options",
        f"",
        f"| Option | Default | Description |",
        f"|--------|---------|-------------|",
    ])
    sections.extend(cli_lines)
    sections.extend([
        f"",
        f"## Exception Hierarchy",
        f"",
        f"All exceptions inherit from `FastqDLError`:",
        f"",
    ])
    sections.extend(exc_lines)
    sections.extend([
        f"",
        f"Exit codes: `1` = validation/provider/download error, `2` = empty/not found/missing, `3` = partial download",
        f"",
        f"## External Tool Dependencies",
        f"",
    ])
    sections.extend(tool_lines)
    sections.extend([
        f"",
        f"## Development",
        f"",
        f"```bash",
        f"poetry install          # Install dependencies",
        f"just test-cov           # Unit tests with 70% coverage minimum",
        f"just test-integration   # Integration tests (real API calls)",
        f"just fmt                # Format with ruff",
        f"just check              # Format + lint check",
        f"just build              # Build package",
        f"```",
        f"",
        f"## Machine-Readable Metadata",
        f"",
        f"See `catalog.json` for structured project metadata including detailed function signatures,",
        f"CLI option definitions, dependency versions, and the complete exception hierarchy.",
        f"",
    ])

    content = "\n".join(sections)

    return content


def main():
    print("Parsing source code...")
    catalog = build_catalog()

    print(f"Writing {CATALOG_PATH}...")
    with open(CATALOG_PATH, "w") as f:
        json.dump(catalog, f, indent=2, default=str)
        f.write("\n")

    print(f"Writing {LLMS_PATH}...")
    llms_content = generate_llms_txt(catalog)
    LLMS_PATH.write_text(llms_content)

    # Summary
    n_modules = len(catalog["modules"])
    n_functions = sum(
        len(m.get("functions", {})) for m in catalog["modules"].values()
    )
    n_options = len(catalog["cli_options"])
    n_exceptions = len(catalog["exceptions"])
    print(
        f"Done: {n_modules} modules, {n_functions} functions, "
        f"{n_options} CLI options, {n_exceptions} exceptions"
    )


if __name__ == "__main__":
    main()
