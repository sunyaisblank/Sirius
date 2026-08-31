#!/usr/bin/env python3
"""Shared parsing and integrity checks for first-party GoogleTest sources."""

from dataclasses import dataclass
import re
from pathlib import Path


TEST_PATTERN = re.compile(
    r"\bTEST(?:_F)?\(\s*(?P<suite>\w+)\s*,\s*(?P<case>\w+)\s*\)\s*\{"
)
UNSUPPORTED_PATTERN = re.compile(
    r"\b(?P<macro>TEST_P|TYPED_TEST|TYPED_TEST_P)\s*\("
)
INACTIVE_PREPROCESSOR_PATTERN = re.compile(r"(?m)^\s*#\s*if\s+0(?:\D|$)")
DIRECT_POSTCONDITION_PATTERN = re.compile(r"\b(?:EXPECT|ASSERT)_[A-Z0-9_]+\s*\(")
SKIP_PATTERN = re.compile(r"\bGTEST_SKIP\s*\(")

# These cases can dispatch CPU/Vulkan work or publish a completed frame. Their
# exact identity is shared with the CTest label generator, and source governance
# requires them to live under tests/render so the no-render application target
# cannot acquire them through a source-list edit.
APP_RENDERING_TESTS = frozenset({
    "RenderCommandParse.ExplicitGpuRequestRunsVulkanWhenDevicePresent",
    "RenderCommandParse.BackendVulkanDeclinesMetricOffTheRenderPath",
    "RenderCommandParse.ReusedCommandDoesNotRetainAnEarlierGpuRequest",
    "RenderCommandParse.CliCpuOverridesLowerLayerVulkanBackend",
    "ViewCommandOperational.HeadlessRefinementProducesASynchronisedFrame",
    "ViewCommandOperational.VulkanRefinementPublishesProgressiveFrames",
})

# These public Execute calls deliberately exercise a parse or validation
# rejection. Any new unguarded command execution in tests/app requires an
# explicit policy decision instead of silently becoming a render-capable test.
APP_NON_RENDER_RENDER_COMMAND_EXECUTIONS = frozenset({
    "RenderCommandParse.UnknownMetricFailsValidation",
    "RenderCommandParse.UnknownOptionRejected",
    "RenderCommandParse.TrailingNumericGarbageRejected",
    "RenderCommandParse.NonFiniteNumericValueRejected",
    "RenderCommandParse.UnexpectedPositionalArgumentRejected",
    "RenderCommandParse.MetricOverrideDefaultsMassWithoutDiscardingExplicitInput",
    "RenderCommandParse.RetiredBackendNamesAreNotRemapped",
    "RenderCommandParse.MalformedCameraBetaRejected",
})
APP_NON_RENDER_VIEW_COMMAND_EXECUTIONS = frozenset({
    "ViewCommandOperational.StrictParsingAndSessionProjection",
    "ViewCommandOperational.RelativisticJetsDeclineBeforeViewerInitialisation",
})

_OBVIOUS_NO_OP_PATTERNS = (
    re.compile(
        r"\b(?:EXPECT|ASSERT)_TRUE\s*\(\s*(?:true|1)\s*\)"
        r"|\b(?:EXPECT|ASSERT)_FALSE\s*\(\s*(?:false|0)\s*\)"
    ),
    re.compile(r"\bSUCCEED\s*\(\s*\)"),
    re.compile(
        r"\b(?:EXPECT|ASSERT)_(?:EQ|NE|LT|LE|GT|GE|DOUBLE_EQ|FLOAT_EQ)\s*"
        r"\(\s*(?P<operand>[A-Za-z_]\w*(?:\.[A-Za-z_]\w*)*|[-+]?\d+(?:\.\d+)?)"
        r"\s*,\s*(?P=operand)\s*\)"
    ),
    re.compile(
        r"\b(?:EXPECT|ASSERT)_NEAR\s*"
        r"\(\s*(?P<operand>[A-Za-z_]\w*(?:\.[A-Za-z_]\w*)*|[-+]?\d+(?:\.\d+)?)"
        r"\s*,\s*(?P=operand)\s*,"
    ),
)


@dataclass(frozen=True)
class SourceTest:
    suite: str
    case: str
    path: Path
    body_code: str

    @property
    def name(self):
        return f"{self.suite}.{self.case}"

    @property
    def has_direct_postcondition(self):
        return DIRECT_POSTCONDITION_PATTERN.search(self.body_code) is not None

    @property
    def uses_skip(self):
        return SKIP_PATTERN.search(self.body_code) is not None


def _mask_cpp_non_code(content):
    """Replace comments and literals with spaces while preserving offsets."""
    masked = list(content)
    length = len(content)
    index = 0

    def blank(start, end):
        for offset in range(start, end):
            if content[offset] not in "\r\n":
                masked[offset] = " "

    while index < length:
        if content.startswith("//", index):
            end = content.find("\n", index + 2)
            if end < 0:
                end = length
            blank(index, end)
            index = end
            continue
        if content.startswith("/*", index):
            end = content.find("*/", index + 2)
            if end < 0:
                end = length
            else:
                end += 2
            blank(index, end)
            index = end
            continue
        if content.startswith('R"', index):
            delimiter_end = content.find("(", index + 2, min(length, index + 19))
            if delimiter_end >= 0:
                delimiter = content[index + 2:delimiter_end]
                terminator = ")" + delimiter + '"'
                end = content.find(terminator, delimiter_end + 1)
                end = length if end < 0 else end + len(terminator)
                blank(index, end)
                index = end
                continue
        if content[index] in {'"', "'"}:
            quote = content[index]
            end = index + 1
            while end < length:
                if content[end] == "\\":
                    end += 2
                    continue
                end += 1
                if content[end - 1] == quote:
                    break
            blank(index, min(end, length))
            index = end
            continue
        index += 1
    return "".join(masked)


def _matching_brace(masked, opening):
    depth = 0
    for index in range(opening, len(masked)):
        if masked[index] == "{":
            depth += 1
        elif masked[index] == "}":
            depth -= 1
            if depth == 0:
                return index
    return None


def _display_path(path, root):
    try:
        return path.relative_to(root)
    except ValueError:
        return path


def _enforce_execution_boundary(record, display_path):
    """Keep render-capable execution out of the application test target."""
    parts = display_path.parts
    under_app = len(parts) >= 2 and parts[0:2] == ("tests", "app")
    under_render = len(parts) >= 2 and parts[0:2] == ("tests", "render")

    if record.name in APP_RENDERING_TESTS and not under_render:
        raise ValueError(
            f"{display_path} places rendering test {record.name} outside tests/render"
        )
    if not under_app:
        return

    body = record.body_code
    if (re.search(r"\bInteractiveViewer\b", body) and
            re.search(r"\.\s*Start\s*\(", body)):
        raise ValueError(
            f"{display_path} lets application test {record.name} start a rendering viewer"
        )
    if (re.search(r"\bRenderSession\b", body) and
            re.search(r"\.\s*Execute\s*\(", body)):
        raise ValueError(
            f"{display_path} lets application test {record.name} execute a render session"
        )

    invokes_execute = re.search(r"\.\s*Execute\s*\(", body) is not None
    if re.search(r"\bRenderCommand\b", body) and invokes_execute:
        parse_only = "kStopSentinel" in body
        if not parse_only and record.name not in APP_NON_RENDER_RENDER_COMMAND_EXECUTIONS:
            raise ValueError(
                f"{display_path} has unclassified RenderCommand execution in {record.name}"
            )
    if (re.search(r"\bViewCommand\b", body) and invokes_execute and
            record.name not in APP_NON_RENDER_VIEW_COMMAND_EXECUTIONS):
        raise ValueError(
            f"{display_path} has unclassified ViewCommand execution in {record.name}"
        )


def _enforce_no_app_execution_helper(masked, content, display_path, test_spans):
    """Reject an app-side helper that hides Start/Execute outside a test body."""
    parts = display_path.parts
    if len(parts) < 2 or parts[0:2] != ("tests", "app"):
        return
    outside_tests = list(masked)
    for start, end in test_spans:
        for offset in range(start, end):
            if outside_tests[offset] not in "\r\n":
                outside_tests[offset] = " "
    hidden = re.search(r"\.\s*(?:Start|Execute)\s*\(", "".join(outside_tests))
    if hidden:
        line = content.count("\n", 0, hidden.start()) + 1
        raise ValueError(
            f"{display_path}:{line} hides app-side Start/Execute in a test helper"
        )


def collect_source_tests(test_dir, root):
    """Return exact test records and reject source shapes outside policy."""
    records = {}
    for path in sorted(test_dir.rglob("*_test.cpp")):
        content = path.read_text(encoding="utf-8")
        masked = _mask_cpp_non_code(content)
        display_path = _display_path(path, root)
        test_spans = []

        unsupported = UNSUPPORTED_PATTERN.search(masked)
        if unsupported:
            raise ValueError(
                f"{display_path} uses ungoverned {unsupported.group('macro')}"
            )
        inactive = INACTIVE_PREPROCESSOR_PATTERN.search(masked)
        if inactive:
            line = content.count("\n", 0, inactive.start()) + 1
            raise ValueError(
                f"{display_path}:{line} contains an inactive #if 0 source block"
            )

        for pattern in _OBVIOUS_NO_OP_PATTERNS:
            no_op = pattern.search(masked)
            if no_op:
                line = content.count("\n", 0, no_op.start()) + 1
                expression = content[no_op.start():no_op.end()].strip()
                raise ValueError(
                    f"{display_path}:{line} contains obvious no-op assertion {expression!r}"
                )

        for match in TEST_PATTERN.finditer(masked):
            opening = match.end() - 1
            closing = _matching_brace(masked, opening)
            if closing is None:
                line = content.count("\n", 0, opening) + 1
                raise ValueError(f"{display_path}:{line} has an unterminated test body")
            record = SourceTest(
                suite=match.group("suite"),
                case=match.group("case"),
                path=path,
                body_code=masked[opening + 1:closing],
            )
            if record.suite.startswith("DISABLED_") or record.case.startswith("DISABLED_"):
                line = content.count("\n", 0, match.start()) + 1
                raise ValueError(
                    f"{display_path}:{line} contains disabled test {record.name}"
                )
            if record.name in records:
                raise ValueError(
                    f"duplicate test identity {record.name}: "
                    f"{_display_path(records[record.name].path, root)} and {display_path}"
                )
            if not record.has_direct_postcondition:
                line = content.count("\n", 0, opening) + 1
                raise ValueError(
                    f"{display_path}:{line} test {record.name} has no direct "
                    "EXPECT_/ASSERT_ postcondition"
                )
            _enforce_execution_boundary(record, display_path)
            records[record.name] = record
            test_spans.append((opening, closing + 1))
        _enforce_no_app_execution_helper(masked, content, display_path, test_spans)
    return records
