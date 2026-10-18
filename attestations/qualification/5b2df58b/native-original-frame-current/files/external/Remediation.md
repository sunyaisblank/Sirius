# Programming remediation goal

Use two phases: **review carefully and precisely; then remediate efficiently from the established findings.** Review the target repository against justified requirements and reference evidence, create actionable GitHub issues, implement their corrections, and verify the integrated result. An audit or issue list alone is not completion. Begin with the available repository and task context; optional missing details are not a reason to stop.

## Objectives and limits

### Current Alaris goal: critical completion

For Alaris goal #177, the user's explicit scope adds empirical qualification and genuine paper lifecycles to the software completion target below. Cover all ten requested components, including Sector Dispersion and the VIX overlay, and independently reassess Earnings rather than inheriting a retired label. Failed or inconclusive economics cannot be called rollout readiness. Preserve original research results and untouched validation data.

Prioritise the remaining executable strategy paths, adequate data and costs for their actual order types, empirical evaluation, and reconciliation of genuine paper outcomes. External alert delivery, off-host recovery, and permanent WRDS authentication are not completion gates. Existing paper-operation authorisation includes necessary maintenance, bounded trading operations and scheduling; actual live trading still requires separate user authorisation. Paper/live software parity concerns shared execution, risk, accounting and recovery behaviour; simulated fills do not establish live execution performance.

Use one coordinated integration gate for coupled changes. Keep independent source work parallel, but assign one owner to heavy verification and one custodian to the paper broker. Inspect all reported failures before the next run. When only a fixture changes, rebuild its owning test project and rerun the affected checks; retain valid evidence for unchanged production builds and workflows. Do not repeat a study without a specific changed hypothesis, input, or implementation defect and a stated decision its result will enable.

Prefer resolving a real dependency over extending its diagnostic machinery. In particular, use authentic retained provider records to settle compatibility and economic joins. A cost experiment applies only to the instruments, order types and cells actually consumed by a strategy; it is not a universal prerequisite for unrelated strategies. If the required evidence remains unsupported, preserve that result and work on the concrete data or execution dependency instead of repeatedly running the same experiment. None of these efficiency rules lowers the accepted completion criteria.

### General priorities

Use these priorities in order:

1. **Satisfy the intended capabilities and correctness contracts.** Preserve data integrity, required compatibility, and applicable permissions. Never gain speed by disabling required functionality or weakening acceptance criteria.
2. **Maximise useful evidence and coherence.** Prefer independent expectations, checks that expose consequential defects, correctly connected components, clear ownership, and understandable code.
3. **Minimise total delivery and maintenance effort.** Among adequate solutions, prefer less implementation work, verification cost, operational burden, and future maintenance. Optimise elapsed time through dependency-aware parallel work where beneficial.

Minimise unsupported assumptions, unresolved material defects, duplicate authority, unnecessary abstractions, configuration options, dependencies, repeated builds, redundant tests, speculative infrastructure, documentation, and approval interruptions. Counts are not goals: fewer lines, tests, or agents are not inherently better, and more are not inherently stronger.

For each activity ask: **What required outcome, defect, or uncertainty does this address? What is the simplest adequate approach? What evidence will let me stop?** Do not invent numerical optimisation scores or pursue a globally ideal repository.

Do not use speculative human engineer-hours as an estimate of AI elapsed time. Do not stop executable work to announce that the goal will take days. If duration materially matters or the user asks, ground the estimate in inspected remaining work, measured tool runtimes, sequential dependencies, and external waits. Distinguish elapsed time from aggregate agent effort and state uncertainty. Explain the actual bottleneck rather than giving an unsupported number. Before a costly operation, check whether narrower analysis, reuse, batching, or safe parallelism can provide the required evidence. A long estimate triggers examination of scope and process, not automatic acceptance of either. Never promise completion within an arbitrary deadline or lower correctness requirements to meet it.

The completion target is **functional completeness and implementation correctness within the supported contract**. Establish that intended capabilities are implemented, connected, usable, and exercised under relevant controlled conditions. Compilation alone is insufficient. Do not require weeks of monitoring, elapsed paper operation, profitability, future events, or universal absence of defects unless explicitly requested. Use deterministic inputs, controlled clocks, replay, and failure injection where appropriate. Keep claims about observed production performance separate.

If assurance or qualification is itself the component, verify its evaluation behaviour using valid, missing, invalid, and contradictory evidence. Its completion does not require every evaluated capability to receive a positive verdict. An evaluator also cannot substitute for unfinished product functionality.

Respect applicable instructions and permissions. Make routine source changes autonomously and preserve unrelated work. This goal does not independently authorise production deployment, live transactions, destructive data changes, or paid provisioning. Honour existing authorisation; ask only for genuinely missing decisions or permissions and continue independent work.

## Part 1: Review precisely and create issues

### 1. Establish scope and trace the real system

Create a short completion checklist in the conversation: intended capabilities, supported inputs and modes, important correctness conditions, and sufficient acceptance evidence. Treat old plans, generated recommendations, documentation, test names, and qualification labels as claims to examine, not automatic requirements. Keep scope stable; add work only for a material defect affecting the agreed outcome or a necessary repair dependency.

Identify all known in-scope input sources, including user actions, APIs, files, configuration, stored state, clocks, scheduled events, and external responses. Trace each through actual implementations and abstraction layers to its outputs, persisted changes, or external effects. Include joins, alternate modes, rejection, and recovery. Distinguish paths inspected in source from paths actually executed.

Inspect wiring, state ownership, persistence, packaging, and external boundaries. Look for disconnected capabilities, silent no-ops, fabricated fallbacks, duplicate authority, half-migrations, stale configuration, and drift from intended behaviour. Review beyond existing tests; passing tests can exercise the wrong implementation or preserve an old error. Use relevant Git history and logs for context.

Keep initial review read-only apart from isolated builds and diagnostics. Use a compact temporary map, not a new analysis platform. Depth should follow consequence, uncertainty, and coupling. Do not inspect irrelevant dependencies indefinitely.

### 2. Establish the reference contracts

Determine what each important operation should do independently of its current implementation. Use explicit requirements for product intent, applicable standards and official specifications for interfaces, and primary literature, justified derivations, authoritative datasets, or independent reference implementations for domain claims. Consult sources when substantive claims depend on them. Check versions, provenance, assumptions, units, and applicability.

A published source is not automatically guaranteed ground truth. Distinguish exact results, approximations, empirical estimates, and examples. Resolve conflicting evidence where possible; state unresolved assumptions. Where no external oracle exists, derive a justified contract or analytic reference case. Agreement with another implementation is supporting evidence, not proof. Research enough to settle the obligation, then proceed.

Use **Hoare-style reasoning and compositional verification**:

`{P} C {Q}`

- **P:** required input conditions, state, units, and environmental assumptions.
- **C:** the actual implementation and dependencies that execute.
- **Q:** required results, resulting state, permitted effects, and defined failure outcomes.
- **Invariants and frame conditions:** what remains true during execution and what state must remain unchanged.

Apply this to important calculations, decisions, state transitions, and boundaries throughout each critical workflow. Contracts must express useful capability, not merely “does not crash,” “returns something,” or “rejects everything safely.” Do not exclude required inputs by strengthening preconditions. Validate external assumptions where first needed.

For numerical logic, establish equations, domains, units, signs, rounding, tolerances, and boundary behaviour; account for machine arithmetic and numerical instability. For other logic, use precise predicates or transition rules. Do not manufacture equations or annotations for trivial plumbing.

Reason through branches, mutations, and failure exits. Check loop invariant establishment and preservation, plus termination or bounded progress where required. Ordinary partial-correctness Hoare reasoning guarantees a postcondition only **if execution terminates**. For continuous services, assess progress of individual operations.

Check composition: the output and state of stage i must satisfy stage i+1's requirements, `Q_i ⇒ P_(i+1)`, including intervening transformations. Ensure the composed guarantees establish the final required outcome. Inspect units, identifiers, serialization, ownership, ordering, and configuration at joins. For concurrent or asynchronous paths, account for interference, atomicity, retries, and failure timing; sequential reasoning alone is insufficient.

Keep useful contracts in types, assertions, concise comments, or focused tests. Do not create a separate contract catalogue.

### 3. Challenge the contracts and implementation

Start with small deterministic inputs whose expected results are independently justified. Trace them through real code from input to final outcome, inspecting intermediate states where discrepancies can arise. Never copy the implementation to generate its expected answer.

Partition inputs and states by meaningful behavioural differences. Cover valid and invalid cases, boundaries and adjacent values, modes, and consequential combinations. Include relevant event sequences, duplicates, interruptions, restart, and recovery. Derive cases from both the intended contract and actual branching logic.

Choose the smallest effective combination of property checks, independently justified metamorphic relations, exhaustive enumeration of small finite domains, combinatorial tests, and controlled integration checks. Use symbolic execution, bounded model checking, or deductive verification for specific consequential obligations when practical. Do not mandate every method or build a formal-verification platform as a prerequisite.

State assumptions and bounds. Check that preconditions are satisfiable and relevant paths reachable; vacuous success is not evidence. A solver timeout is not a pass. Partitioning does not prove equivalence within each partition, and symbolic execution is not automatically exhaustive. Distinguish source reasoning, executed tests, bounded verification, and machine-checked proof. Never label an AI-written argument as a machine-checked proof.

Exercise real internal connections. Substitute controlled providers only at genuine external boundaries. A provider simulation establishes behaviour against its model, not actual compatibility. Run a bounded authorised external check when available and relevant; otherwise state the precise unverified dependency without waiting indefinitely.

### 4. Create an executable issue queue

Use GitHub issues as the durable remediation queue. This goal authorises creating and updating relevant issues in the target repository. Check for existing issues first. Group shared causes and coherent repairs; do not generate one issue per function, test case, or speculative improvement.

Each issue should contain only what an implementer needs:

- Affected capability and concrete code locations or input-to-outcome path.
- Expected contract and reference evidence, including assumptions.
- Actual discrepancy, consequence, and minimal reproduction, counterexample, or code-level argument.
- Proposed repair, dependencies, and finite acceptance checks.

Separate established defects from suspicions and missing external evidence. Keep source links and necessary reasoning with the issue or relevant test. If GitHub is unavailable, retain equivalent drafts in the conversation and continue local work; do not create a substitute issue database.

Finish the review pass when known input families and critical paths have a clear disposition: supported by reviewed evidence, defective with an actionable issue, or unresolved for a stated reason. This is a declared coverage boundary, not a guarantee against every possible failure. Proceed to remediation without an extra approval checkpoint. For large repositories, first map shared contracts and dependencies, then apply these two phases in coherent batches. Independent, well-established repairs may proceed while unrelated review continues; do not make every repair wait for an exhaustive repository-wide audit.

## Part 2: Remediate rapidly and verify closure

### 5. Schedule and execute the repairs

Use the issues as implementation specifications. Avoid repeating repository discovery or literature review unless new evidence invalidates a finding. Order work by dependencies and consequences; prioritise repairs that unblock several others or lie on the path determining completion time.

This goal authorises parallel agents where available. Use a coordinator to maintain a lightweight dependency order and assign independent, ready issues. Choose agent count according to available independent work, shared-state risk, tool limits, and coordination cost; maximum concurrency is not the objective. Use a single agent when coordination would cost more than it saves.

Give each worker the issue, contract, relevant evidence, owned files or interfaces, and acceptance checks. Use isolated branches or worktrees where appropriate. Assign one owner to shared contracts and closely coupled changes; do not allow conflicting edits or duplicate investigations. Complete shared prerequisites before dependent repairs, or establish their interfaces explicitly. Keep coordination in the existing task tools and issues, not a new swarm framework or report hierarchy.

Workers implement and run focused checks, then return the change, verification result, and any blocker. The coordinator reviews interactions, integrates changes in small coherent batches, resolves conflicts, and verifies the combined workflows. Separate worker successes do not establish integrated correctness. Use an independent review for consequential or uncertain changes, with the contract and diff as its scope; do not require a second full audit for every repair. Do not make every worker rebuild or retest the whole repository.

Repair root causes and update consumers, configuration, and migrations together. Preserve sound behaviour. Confirm real consumers—including dynamic loading and packaging—before deleting apparently unused code. Prefer direct corrections and useful simplification. Revisit an issue if its diagnosis is disproved; never implement an unsupported fix just to follow the plan. New material defects follow the same bounded issue process.

### 6. Verify efficiently and close issues

Each check must address an identified contract, defect, or uncertainty. Reuse trustworthy tests. Inspect affected tests for weak assertions, excessive mocking, stale expectations, duplication, and checks that merely mirror implementation. Repair or remove demonstrably obsolete tests while preserving distinct obligations. Do not launch a suite rewrite just because there are thousands of tests. Keep integration fixtures readable as complete workflows. Do not introduce artificial helper layers solely to reduce a static type-coupling count: a fixture intentionally exercising real production components may use a narrowly scoped, explained metric exception when decomposition would obscure the contract. This does not permit disabling correctness checks, substituting mocks for required consumers, removing assertions, or relaxing production diagnostics.

Prefer regressions that fail on the faulty behaviour and pass after correction. Batch coherent edits, run focused checks, and reuse valid build outputs. Repeat checks when relevant code, configuration, dependencies, fixtures, or evidence change, or a specific failure needs investigation. Use clean rebuilds only for credible stale-state concerns or applicable requirements.

On the integrated final state, run the relevant build, affected workflows, and genuinely required repository gates. Broaden testing for concrete remaining risks or affected dependencies. Do not weaken assertions, tolerances, gates, or functionality to obtain green results. Separate baseline failures from regressions. Stop when sufficient acceptance evidence exists.

Update each issue with the correction, relevant commit when available, checks performed, and material limitations. Close only when its acceptance conditions are met. Explain disproved findings as such. Leave blocked issues open and finish independent work. Do not hide missing required functionality behind an “externally unverified” label.

### 7. Keep documentation minimal and finish

Maintain a single **README.md** as the default prose documentation, written in plain language for a human unfamiliar with the project. Explain why it exists, what it currently does, and how to set up and use it, including relevant limitations. It must not become an AI worklog, audit, backlog, or change history.

Do not add repository plans, evidence catalogues, session summaries, status reports, or handover documents. Keep temporary reasoning outside the repository or in the conversation. Use GitHub issues for defects and acceptance evidence, and Git history for change history; do not duplicate them in files. Retain useful contracts, tests, and concise explanatory comments.

Consolidate current usage information from obsolete AI documents into the README and remove superseded task documents after resolving their references. Preserve independently required legal material, source data, executable specifications, and active instruction files. Remove task-generated scratch material when finished.

Finish when the agreed capabilities are implemented and connected, material in-scope defects are repaired, contracts and composition have proportionate supporting evidence, and no known regression or unfinished required migration remains. Review the final diff for unintended changes, broken references, stale alternatives, and unnecessary complexity.

Do not add optional improvements, future observations, or stronger hypothetical guarantees to the stopping conditions. If genuinely blocked, identify the exact outcome and missing prerequisite. Give a short final response with changes, verification, and remaining limitations; no handover file.

**Begin now: establish scope, trace the system, derive and challenge its contracts, create actionable issues, then remediate and verify them to closure.**
