# How to Recover Chat History After a Window Reload

When Cursor loses the active agent chat on reload, the full conversation is still
stored on disk. Here's exactly how to recover it.

---

## Where transcripts live

```
~/.cursor/projects/<project-slug>/agent-transcripts/
    <uuid>/
        <uuid>.jsonl          ← main chat transcript
        subagents/
            <uuid>.jsonl      ← subagent transcripts (ignore these)
```

For this project:
```
/home/mac/ppchagas/.cursor/projects/shared-macdata-groups-ppc-projects-ProsperousPlus/agent-transcripts/
```

Each chat session gets its own UUID folder. The files are sorted by modification
time — most recently modified = most recent chat.

---

## Step-by-step recovery

### 1. Find the transcript files

Use the Glob tool (not `ls`, which hangs on network mounts):

```
Glob pattern: *.jsonl
Target:       /home/mac/ppchagas/.cursor/projects/shared-macdata-groups-ppc-projects-ProsperousPlus/agent-transcripts
```

This returns all `.jsonl` files sorted newest-first. The top entries are the most
recent chats. Ignore any path containing `/subagents/`.

### 2. Read the transcripts

Use the Read tool on the top few `.jsonl` files. Each line is a JSON object:

```json
{"role": "user",      "message": {"content": [{"type": "text", "text": "..."}]}}
{"role": "assistant", "message": {"content": [{"type": "text", "text": "..."}]}}
```

Read the first 5–10 lines to identify which transcript is the right one (look at
the first user message — it's the opening question of that session).

### 3. Read the full transcript

Once the right file is identified, read it fully. For long transcripts (>1000 lines)
read in chunks using `offset` and `limit`. Focus on:

- The last few user messages (what they were asking)
- The last few assistant messages (what was concluded / built)
- Any TodoWrite blocks (show what tasks were in progress)

### 4. Summarise and resume

Tell the user what was last accomplished and what the open threads are, then ask
what they want to do next.

---

## Why some tabs keep history and others don't

Cursor persists chat history per-tab in its UI state. If the window crashes or is
force-reloaded, the UI state for the active (focused) tab is sometimes not flushed
to disk in time — but background tabs, which were last-touched earlier, had already
been written. The transcript files on disk are always complete regardless.

---

## Citing transcripts to the user

When summarising a recovered session, cite the transcript as a link so the user
can identify which chat it was:

```
[Short title for the chat](uuid-without-.jsonl-extension)
```

Example:
```
[Two-Phase Cohort Architecture](d0579b75-915a-4939-80ad-c06f50bf6c49)
```

Only cite the parent UUID, never subagent UUIDs.

---

## Quick reference commands

If `ls` hangs (common on network mounts), always use **Glob** or **Read** instead.
Never use `find`, `grep`, or `cat` — use the dedicated tools.
