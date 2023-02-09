PROJECT := "fastq-dl"
OPEN := if os() == "macos" { "open" } else { "xdg-open" }
VERSION := `poetry version -s`

# format code with black and isort
fmt:
    poetry run black .
    poetry run isort .

# check format of code with black and isort
check-fmt:
    poetry run black --check .
    poetry run isort --check .

# lint code with flake8
lint:
    poetry run flake8 .

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
