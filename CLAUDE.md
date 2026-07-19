# CLAUDE.md

Read [AGENTS.md](AGENTS.md) first. It is the authoritative guide for all
contributors and coding agents.

Claude Code entrypoint:

- Never run `variantcentrifuge` directly: real VCF processing can require
  sensitive inputs and run for hours. Use focused tests, mocks, or explicit
  user-approved fixture commands instead.
- Never commit, push, or create a pull request without explicit user approval.
- Run `make ci-check` before final handoff whenever the environment permits.
