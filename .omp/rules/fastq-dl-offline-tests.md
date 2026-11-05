---
description: Unit tests stay offline; `responses` mocks the HTTP layer
condition: 'def test_'
globs: ["tests/**/*.py"]
---
CTX-FASTQDL-2 MUST mock every HTTP call in unit tests; network tests carry the `integration` marker.
Why: the unit suite has to run without network.
