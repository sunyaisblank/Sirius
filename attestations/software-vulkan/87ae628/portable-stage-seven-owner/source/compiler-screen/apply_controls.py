import re, sys
from pathlib import Path
def portable_normal_sum_controls(text):
    """Admit guarded binary32 FastTwoSum and positive upward-bound sums."""
    if (re.findall(r"OpCapability (\S+)", text) != ["Shader"] or
            set(re.findall(r"OpTypeInt (\d+) [01]", text)) != {"32"} or
            re.findall(r"OpTypeFloat (\d+)", text) != ["32"]):
        raise ValueError("portable normal sum introduced an optional type or capability")
    if re.search(r"\b(?:FPFastMathMode|FPFastMathDefault|RelaxedPrecision|OpFmaKHR)\b", text):
        raise ValueError("portable normal sum introduced relaxed or fused arithmetic")
    if re.search(r"OpExecutionMode\S* .* (?:Denorm|RoundingMode|SignedZeroInfNan)", text):
        raise ValueError("portable normal sum introduced an unqualified float mode")
    float_type = re.search(r"(%\S+) = OpTypeFloat 32", text)[1]
    operations = re.findall(r"(%\S+) = (OpF(?!unction)\S+) (%\S+) ", text)
    bodies = [body for body in re.findall(r"%\S+ = OpFunction .*?OpFunctionEnd", text, re.S)
              if re.search(r"= OpF(?!unction)", body)]
    expected = []
    if not bodies:
        raise ValueError("portable normal sum lost its guarded arithmetic")
    for body in bodies:
        if re.match(r"%RPUpAdd(?:_\d+)? = OpFunction ", body):
            if re.findall(r"= (OpF(?!unction)\S+) ", body) != ["OpFAdd"]:
                raise ValueError("normal upward bound changed its single add body")
            expected.append("OpFAdd")
            continue
        if (not re.match(r"%RetainedSum32(?:_\d+)? = OpFunction ", body) or
                re.findall(r"= (OpF(?!unction)\S+) ", body) != ["OpFAdd", "OpFSub", "OpFSub"]):
            raise ValueError("portable normal sum changed its ordered three-operation body")
        expected.extend(["OpFAdd", "OpFSub", "OpFSub"])
        links = re.findall(r"(%\S+) = OpF(?:Add|Sub) %\S+ (%\S+) (%\S+)", body)
        copies = dict(re.findall(r"(%\S+) = OpCopyObject " + re.escape(float_type) + r" (%\S+)", body))
        def origin(value):
            seen = set()
            while value in copies:
                if value in seen:
                    raise ValueError("portable normal sum has a cyclic copy chain")
                seen.add(value)
                value = copies[value]
            return value
        total, a, b = links[0]
        displaced, left, right = links[1]
        _, residual_left, residual_right = links[2]
        if (origin(left) != total or origin(right) != origin(a) or
                origin(residual_left) != origin(b) or origin(residual_right) != displaced):
            raise ValueError("portable normal sum changed its FastTwoSum dependency chain")
    decorated = set(re.findall(r"OpDecorate (%\S+) NoContraction", text))
    if ([op for _, op, _ in operations] != expected or
            any(value not in decorated or kind != float_type for value, _, kind in operations)):
        raise ValueError("portable normal sum lost binary32 operation ordering or NoContraction")
    if (re.search(r"\bOpConvert(?:FTo[SU]|[SU]ToF)\b|\bOpDot\b|\bOpQuantizeToF16\b|"
                  r"\bOp(?:Vector|Matrix)Times\w+\b|\bOpOuterProduct\b", text) or
            re.search(r"= OpExtInst " + re.escape(float_type) + r" ", text)):
        raise ValueError("portable normal sum introduced additional floating arithmetic")
    entry = re.search(r"OpEntryPoint GLCompute (%\S+)", text)[1]
    text = text.replace("OpCapability Shader", "OpCapability Shader\n               OpCapability RoundingModeRTE", 1)
    local_size = re.search(r"^.*OpExecutionMode " + re.escape(entry) + r" LocalSize.*$", text, re.M)[0]
    return text.replace(local_size, local_size + "\n               OpExecutionMode " + entry + " RoundingModeRTE 32", 1)

source, destination = map(Path, sys.argv[1:])
destination.write_text(portable_normal_sum_controls(source.read_text()))
