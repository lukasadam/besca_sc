# Git Commit Cheat Sheet

## Conventional Commits

Use the [Conventional Commits](https://www.conventionalcommits.org/) convention to keep Git history readable, consistent, and suitable for automated changelog generation.

### Commit message format

```text
<type>(<scope>): <short description>
```

Example:

```bash
git commit -m "feat(qc): add MAD-based outlier filtering"
```

- **Type:** Category of change (`feat`, `fix`, `docs`, etc.).
- **Scope:** Optional component affected (`qc`, `plotting`, `preprocessing`, etc.).
- **Description:** Short, imperative description of the change.

## Commit types

| Type       | Description                                | Example                                    |
| ---------- | ------------------------------------------ | ------------------------------------------ |
| `feat`     | New functionality                          | `feat(qc): add MAD-based filtering`        |
| `fix`      | Bug fix                                    | `fix(qc): handle zero MAD values`          |
| `refactor` | Restructure code without changing behavior | `refactor(qc): extract outlier detection`  |
| `perf`     | Performance improvement                    | `perf(qc): vectorize filtering operations` |
| `docs`     | Documentation changes                      | `docs: add Pixi cheatsheet`                |
| `test`     | Add or update tests                        | `test(qc): add edge-case tests`            |
| `build`    | Dependencies, packaging, build system      | `build: configure Pixi environments`       |
| `ci`       | Continuous integration changes             | `ci: add GitHub Actions tests`             |
| `chore`    | Repository maintenance                     | `chore: update gitignore`                  |
| `style`    | Formatting without changing behavior       | `style: apply Ruff formatting`             |
| `revert`   | Revert a previous commit                   | `revert: remove faulty normalization`      |

Note: Only `feat` and `fix` are specifically defined as standard types in the Conventional Commits specification. The others are widely used conventions.

## Breaking changes

Use `!` after the type or scope to indicate an incompatible change.

```bash
git commit -m "feat(api)!: rename filtering parameters"
```

For additional details, include a `BREAKING CHANGE:` footer in the commit body.

## Best practices

1. **Use imperative mood:** `add`, `fix`, `remove` rather than `added`, `fixed`, `removed`.
2. **Keep the subject concise:** Ideally under 50–72 characters.
3. **Commit one logical change at a time:** Avoid mixing unrelated features, bug fixes, and documentation changes.
4. **Use lowercase types:** `feat`, not `Feat`.
5. **Explain why when necessary:** Use the commit body for non-obvious design decisions.

## Useful Git commands

```bash
git status                        # Show changed files
git diff                          # Inspect unstaged changes
git add -p                        # Stage individual changes
git diff --cached                 # Review staged changes
git commit -m "fix(qc): handle NaN values"
git log --oneline                 # View commit history
git show HEAD                     # Inspect latest commit
git commit --amend                # Amend latest commit
```

## Quick decision guide

| Question                                         | Commit type |
| ------------------------------------------------ | ----------- |
| Did I add new functionality?                     | `feat`      |
| Did I fix incorrect behavior?                    | `fix`       |
| Did I reorganize code without changing behavior? | `refactor`  |
| Did I make existing code faster?                 | `perf`      |
| Did I only change documentation?                 | `docs`      |
| Did I add or modify tests?                       | `test`      |
| Did I change dependencies or packaging?          | `build`     |
| Did I change CI/CD configuration?                | `ci`        |
| Did I only change formatting?                    | `style`     |
| Did I perform general maintenance?               | `chore`     |
| Did I undo an earlier commit?                    | `revert`    |