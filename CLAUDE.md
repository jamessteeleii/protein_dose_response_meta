@AGENTS.md

# Claude Code-specific instructions

- Treat `AGENTS.md` as the canonical shared project policy. Apply all imported rules before these Claude-specific additions.
- Do not duplicate shared rules here. When the user asks to change project policy, edit `AGENTS.md`; reserve this file for Claude Code behavior only.
- Keep durable project facts and decisions in tracked project documentation, not solely in auto memory.
- Do not add `CLAUDE.local.md`, personal memory, or machine-specific settings to Git.
- If instructions appear missing or stale, ask the user to confirm this file under `/context` after starting a new Claude Code session.
- On Windows, preserve the `@AGENTS.md` import instead of replacing it with a symbolic link.
