"""Private header and exact ELF payload joins, adapted from preserved b66 tools."""
from pathlib import Path
import hashlib,re,struct
NAMESPACE="retained_program_literal_"

def require(ok,message):
    if not ok: raise ValueError(message)

ARRAY = re.compile(
    r"inline constexpr std::array<std::uint32_t, (\d+)> (\w+)\{\{(.*?)\}\};", re.S
)


def identity(path):
    path = Path(path).resolve(strict=True)
    data = path.read_bytes()
    return {"path": str(path), "bytes": len(data), "sha256": hashlib.sha256(data).hexdigest()}


def replace_once(text, before, after):
    assert text.count(before) == 1, before[:160]
    return text.replace(before, after, 1)


def arrays(text):
    result = {}
    for match in ARRAY.finditer(text):
        count, name, body = match.groups()
        assert name not in result, name
        values = [int(x.strip().removesuffix("u"), 0) for x in body.split(",") if x.strip()]
        assert len(values) == int(count), name
        result[name] = struct.pack("<" + "I" * len(values), *values)
    assert len(result) == 40
    assert sum(name.endswith("Program") for name in result) == 7
    assert sum(name.endswith("Shader") for name in result) == 33
    return result


def replace_payload(text, name, payload):
    assert len(payload) % 4 == 0
    words = struct.unpack("<" + "I" * (len(payload) // 4), payload)
    matches = [m for m in ARRAY.finditer(text) if m.group(2) == name]
    assert len(matches) == 1, name
    match = matches[0]
    lines = [",".join(str(v) + "u" for v in words[i:i + 16]) + "," for i in range(0, len(words), 16)]
    replacement = (
        f"inline constexpr std::array<std::uint32_t, {len(words)}> {name}{{{{\n"
        + "\n".join(lines) + "\n}};"
    )
    return text[:match.start()] + replacement + text[match.end():]


def elf_arrays(exe, arrays, namespace, allowed_production=None):
    # Read each defined object at its actual ELF symbol address, not by a
    # payload substring search or a symbol-name-only check.
    raw = exe.read_bytes()
    require(raw[:6] == b"\x7fELF\x02\x01", "Expected little-endian ELF64")
    offset = struct.unpack_from("<Q", raw, 40)[0]
    entry, count = struct.unpack_from("<HH", raw, 58)
    require(entry == 64 and count > 0 and offset + entry * count <= len(raw),
              "ELF section table outside file")
    sections = [struct.unpack_from("<IIQQQQIIQQ", raw, offset + entry * i) for i in range(count)]
    symbols = {}
    for section in sections:
        if section[1] != 2:
            continue
        require(section[6] < count and section[9] == 24 and section[5] % 24 == 0,
                  "ELF symbol table malformed")
        strings_section = sections[section[6]]
        require(strings_section[1] == 3 and strings_section[4] + strings_section[5] <= len(raw),
                  "ELF symbol string table outside file")
        strings = raw[strings_section[4]:strings_section[4] + strings_section[5]]
        require(section[4] + section[5] <= len(raw), "ELF symbol table outside file")
        for at in range(section[4], section[4] + section[5], 24):
            name, info, _, index, address, size = struct.unpack_from("<IBBHQQ", raw, at)
            require(name < len(strings), "ELF symbol name outside string table")
            end = strings.find(b"\0", name)
            require(end >= 0, "ELF symbol name unterminated")
            symbol = strings[name:end].decode()
            if info & 15 == 1 and 0 < index < count:
                symbols.setdefault(symbol, []).append((index, address, size))
    production_prefix = "_ZN6sirius7backend16retained_program"
    expected_normal = arrays if namespace == "retained_program" else (allowed_production or {})
    expected_symbols = {production_prefix + str(len(name)) + name + "E" for name in expected_normal}
    require({name for name in symbols if name.startswith(production_prefix)} == expected_symbols,
            "Only the original test object's exact normal-reference shader objects may coexist")
    records = {}
    for name, payload in arrays.items():
        symbol = "_ZN6sirius7backend" + str(len(namespace)) + namespace + str(len(name)) + name + "E"
        matches = symbols.get(symbol, [])
        require(len(matches) == 1, "Normal defined array object absent/ambiguous: " + name)
        index, address, size = matches[0]
        section = sections[index]
        relative = address - section[3]
        require(section[1] != 8 and size == len(payload) and relative >= 0 and
                  relative + size <= section[5] and section[4] + relative + size <= len(raw),
                  "ELF array range differs: " + name)
        at = section[4] + relative
        require(raw[at:at + size] == payload, "ELF symbol-address payload differs: " + name)
        records[name] = {"symbol": symbol, "address": address, "file_offset": at,
                         "bytes": size, "sha256": hashlib.sha256(payload).hexdigest()}
    if namespace.startswith(NAMESPACE):
        other = "candidate" if namespace.endswith("baseline") else "baseline"
        require((NAMESPACE + other).encode() not in raw, "Other private factory linked")
    return records

