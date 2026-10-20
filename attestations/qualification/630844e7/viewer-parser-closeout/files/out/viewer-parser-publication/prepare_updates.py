"""Prepare exact issue/PR body corrections; acceptance text stays authoritative."""
from pathlib import Path
import difflib
import json

ROOT = Path.cwd()
HERE = Path(__file__).resolve().parent
BEFORE = ROOT / 'out/viewer-constant-parser-closeout'
SOURCE = '630844e775a85e74eb5d9ed328a2b8019abb640d'
link = '[630844e7](https://github.com/sunyaisblank/Sirius/commit/' + SOURCE + ')'
evidence = 'https://github.com/sunyaisblank/Sirius/tree/079d204b968d7f1c9967b589d28de70017999b96/attestations/qualification/630844e7/viewer-parser-source-build-and-display'
history = 'https://github.com/sunyaisblank/Sirius/tree/e9b8928b6f7280c79f552b3bd2e5f3e46682ed2e/attestations/qualification/2c8ef086/current-build-domain-admission-and-sanitizer-stop'
updates = {}
for number in (44, 48, 53):
    original = json.loads((BEFORE / f'issue-{number}-before.json').read_text())
    body = original['body']
    if number == 44:
        start = body.index('Current development source is ')
        end = body.index('\n\n', start)
        replacement = ('Current development source is clean, committed and pushed at ' + link +
                       ' on `development/sanitizer-shader-constant-parser-2026-10-20`. The bounded native FP32 constructor remains integrated from 2c8ef086 after all eight original controls and the sole original 22/22 benefit gates passed. The subsequent change repairs only the viewer test helper that blocked the strict sanitizer build; production C++ and shaders are unchanged. Main remains 5e706073a3fd2f9afcf8efbcde2f094122050f6a; [draft PR #62](https://github.com/sunyaisblank/Sirius/pull/62) remains unmerged. Source/build milestones and finite controls do not establish full scientific, rendering or release qualification.\n\n'
                       'Current GCC/Clang and Debug/development ASan+UBSan builds pass all seven source targets; the two original display contracts pass on each, with unchanged 4,097 transfer samples, 2e-15 tolerance and film-range/alpha assertions. [Preserved builds and six display-case executions](' + evidence + '). No failed frame or rejected benefit trial was repeated. The two native-build domains were admitted into GCC/Clang at exact 2c8ef086 ([preserved admission](' + history + ')); current 630844e7 admits 0/8. The trailing 5b2df58b configuration note remains a historical snapshot. Full Mandatory, scientific, provider, physical and installed-release acceptance stays open.')
        body = body[:start] + replacement + body[end:]
        body = body.replace('The integrated source is the exact reviewed constructor candidate:',
                            'The retained-shader integration at 2c8ef086 is the exact reviewed constructor candidate:', 1)
        assert body[body.index('Acceptance remains:'):] == original['body'][original['body'].index('Acceptance remains:'):]
    elif number == 48:
        tail = body[body.index('Affected source includes'):]
        domains = body[body.index('- macos-native-build'):body.index('Affected source includes')]
        body = ('Complete retained Vulkan runtime, scientific/provider qualification and release admission remain open. Source/build and finite controls do not establish the supported runtime estate.\n\n'
                'Current development source is clean, committed and pushed at ' + link + ' on `development/sanitizer-shader-constant-parser-2026-10-20`. Main remains 5e706073a3fd2f9afcf8efbcde2f094122050f6a; [draft PR #62](https://github.com/sunyaisblank/Sirius/pull/62) remains unmerged.\n\n'
                'All seven executables build with the existing strict GCC/Clang presets and genuine Debug/development ASan+UBSan preset. GCC14 had stopped the prior sanitizer build inside the viewer shader-constant regex helper; a small decimal/fraction parser fixes that compile blocker without changing production C++, shaders, original display assertions, contracts or warnings-as-errors. The two original display contracts pass on each build: six executions, zero failures/errors/skips/disabled cases. All 35 retained SPVs, three generated headers and 19 physics outputs match across the three presets. Live registration is 1,152 per preset; this is registration only. [Preserved source/build/display evidence](' + evidence + '). The failed harness/readers and their exact corrections are retained; no full Mandatory, scientific/provider or render workload ran.\n\n'
                'The exact 2c8ef086 macOS/Windows build-only records from integration CI 39116392750 were verified, restored from their original archives and actually admitted by both GCC/Clang consumers: 2/8 domains, six pending, with both original authority cases passing on each compiler. [Preserved admission and original sanitizer compile stop](' + history + '). This is historical exact-revision evidence; it is not relabelled for 630844e7. The earlier five host contract cases remain bounded to their recorded 5b2df58b execution, and are not counted as fresh current-source sanitizer science.\n\n'
                'Current 630844e7 configurations admit **0/8 external domains**; all eight remain pending. Overall readiness and release readiness remain false. Qualification-mode GCC/Clang evidence-generation readiness is false because complete Mandatory evidence is absent. Debug/development sanitizer mode intentionally allows evidence generation and does not require a Mandatory receipt; that is not release or full-estate admission. The existing sanitizer test-preset ASan/LSan/UBSan environment, including its existing Vulkan-initialisation suppression, was retained unchanged. The domains are:\n\n' + domains + tail)
        assert body[body.index('Required acceptance:'):] == original['body'][original['body'].index('Required acceptance:'):]
    else:
        start = body.index('Current product/scientific source is ')
        end = body.index('\n\n', start)
        body = body[:start] + ('Current development source is committed and pushed at ' + link +
               ' on `development/sanitizer-shader-constant-parser-2026-10-20`. [Draft PR #62](https://github.com/sunyaisblank/Sirius/pull/62) remains at 4bcd987 on `development/retained-dense-step-control-2026-10-06`; later development is not integrated into that draft. Main remains `5e706073a3fd2f9afcf8efbcde2f094122050f6a`. No further CPU endpoint-ownership defect has been established by the recorded reviews. The current test-helper compile repair and six display checks ([evidence](' + evidence + ')) are not a new execution of the scientific results below.') + body[end:]
        for old, new in (
            ('At current 4bcd987,', 'At 4bcd987,'),
            ('same current primary GCC producers', 'primary GCC producers recorded for that revision'),
            ('three prior current executions', 'three prior executions at that revision'),
            ('finite current-GCC evidence', 'finite GCC evidence at the revisions identified in the linked result'),
            ('the current oracle producer', 'the recorded oracle producer'),
            ('[Current scientific result and preserved evidence]', '[Recorded scientific result and preserved evidence]'),
            ('[Current Page–Thorne result and preserved evidence]', '[Recorded Page–Thorne result and preserved evidence]'),
            ('Current qualification is 2/8 build-only domains, six pending and readiness/release false.',
             'Historical exact 2c8ef086 build-domain admission reached 2/8; current 630844e7 configurations admit 0/8. Overall readiness and release readiness remain false.')):
            body = body.replace(old, new)
    assert body != original['body'] and original['state'] == 'OPEN'
    (HERE / f'issue-{number}-before.json').write_text(json.dumps(original, indent=2) + '\n')
    (HERE / f'issue-{number}-body.txt').write_text(body)
    updates[str(number)] = {'before_updated_at': original['updatedAt'], 'characters_before': len(original['body']), 'characters_after': len(body), 'acceptance_preserved': True}
    (HERE / f'issue-{number}.diff').write_text(''.join(difflib.unified_diff(original['body'].splitlines(True), body.splitlines(True))))

original = json.loads((BEFORE / 'pr-62-before.json').read_text())
assert original['state'] == 'OPEN' and original['isDraft'] and original['headRefOid'] == '4bcd9875cd9308deabed06a5fa2f17d2321b0cbc'
body = original['body'].replace('[2c8ef086](https://github.com/sunyaisblank/Sirius/commit/2c8ef08693d46f57d5dea14efdbddd4874dac684)', link)
body = body.replace('`development/native-constant-constructor-2026-10-18`', '`development/sanitizer-shader-constant-parser-2026-10-20`')
(HERE / 'pr-62-before.json').write_text(json.dumps(original, indent=2) + '\n')
(HERE / 'pr-62-body.txt').write_text(body)
(HERE / 'pr-62.diff').write_text(''.join(difflib.unified_diff(original['body'].splitlines(True), body.splitlines(True))))
(HERE / 'prepared.json').write_text(json.dumps({'source_revision': SOURCE, 'evidence_commit': '079d204b968d7f1c9967b589d28de70017999b96', 'issue_updates': updates, 'PR_head_base_draft_and_validation_unchanged': True}, indent=2) + '\n')
print(json.dumps(updates))
