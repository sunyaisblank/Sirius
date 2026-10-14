import json
import re
import struct

def interface_signature(text):
    """Normalize IDs only for the complete descriptor/shared type closure."""
    definitions = dict(re.findall(r"^\s*(%\S+) = (.+)$", text, re.M))
    names = dict(re.findall(r'OpName (%\S+) "([^"]*)"', text))
    decorations = {}
    for owner, suffix in re.findall(r"Op(?:Member)?Decorate (%\S+) (.+)$", text, re.M):
        decorations.setdefault(owner, []).append(suffix)

    def signature(identifier, active=()):
        assert identifier not in active, "cyclic interface type"
        line = definitions[identifier]
        assert line.startswith(("OpType", "OpConstant ")), line
        values = [signature(word, active + (identifier,)) if word.startswith("%") else word
                  for word in line.split()]
        return [values, sorted(decorations.get(identifier, []))]

    answer = []
    for identifier, line in definitions.items():
        words = line.split()
        if words[0] == "OpVariable" and words[2] in (
                "StorageBuffer", "Uniform", "Workgroup", "PushConstant"):
            answer.append([names.get(identifier, ""), words[2], signature(words[1]),
                           sorted(decorations.get(identifier, []))])
    assert answer
    return sorted(answer, key=lambda row: json.dumps(row))
