PROJECT := "fastq-dl"
OPEN := if os() == "macos" { "open" } else { "xdg-open" }
VERSION := `poetry version -s`

# format code with ruff
fmt:
    poetry run ruff format .
    poetry run ruff check --fix .

# check format of code with ruff
check-fmt:
    poetry run ruff format --check .

# lint code with ruff
lint:
    poetry run ruff check .

# install latest version with poetry
install:
    poetry install --no-interaction

# check formatting, linting, and tests
check: check-fmt lint

# prints out the commands to run to tag the release and push it
tag:
    @echo "Run \`git tag -a {{ VERSION }} -m <message>\` to tag the release"
    @echo "Then run \`git push origin {{ VERSION }}\` to push the tag"

# build a python release
build:
    poetry build --no-interaction

# run the unit tests
test args="":
    poetry run pytest {{ args }}

# run unit tests with coverage
test-cov args="":
    poetry run pytest tests/ \
        --cov=fastq_dl \
        --cov-report=term-missing \
        --cov-report=html \
        --cov-fail-under=70 \
        -m "not integration" \
        {{ args }}

# run only integration tests (makes real API calls)
test-integration args="":
    poetry run pytest tests/ \
        -m "integration" \
        {{ args }}

# open coverage report in browser
coverage-report:
    {{ OPEN }} htmlcov/index.html